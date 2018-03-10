!-------------------------------------------------------------------------------
! MODULE: HamiltonianActions
!> @brief
!> Calculate effective field by applying the derivative of the Hamiltonian
!> @details The effective field, \f$\mathbf{B}_i\f$, on an atom \f$\textit{i}\f$, is calculated from
!> \f$ \mathbf{B}_i=-\frac{\partial \mathbf{H}}{\partial \mathbf{m}_i},\f$ where primarily the part of
!> the Hamiltonian, \f$\mathbf{H}\f$, which represents interatomic exchange interactions,
!> \f$\mathbf{H}_\mathrm{ex}\f$, are considered. For this we use the classical Heisenberg Hamiltonian,
!> \f$ \mathbf{H}_\mathrm{ex}=-\frac{1}{2}\sum_{i\neq j}J_{ij}\mathbf{m}_i\cdot\mathbf{m}_j,\f$ where
!> \f$i\f$ and \f$j\f$ are atomic indices and \f$J_{ij}\f$ is the strength of the exchange interaction,
!> which is calculated from first principles theory.
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module HamiltonianActions
   use Parameters
   use Profiling
   use HamiltonianData


contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: effective_field
   !> @brief
   !> Calculate effective field by applying the derivative of the Hamiltonian
   !> @author Anders Bergman, Lars Bergqvist, Johan Hellsvik
   !> @todo Check consistency of terms wrt the input parameters, especially the anisotropies
   !> @todo Replace moment unit vectors emom with full length vectors emomM
   !ham%dm_vect !< DM vector \f$H_{DM}=\sum{D_{ij}\dot(m_i \times m_j)}\f$
   !> @todo Check the sign of the dipolar field
   !-----------------------------------------------------------------------------
   subroutine effective_field(Natom, Mensemble, start_atom, stop_atom, &
         do_jtensor, exc_inter, do_dm, do_pd, do_biqdm, do_bq, &
         do_dip,emomM, mmom, external_field,time_external_field, beff, beff1, beff2, &
         OPT_flag, max_no_constellations, maxNoConstl, &
         unitCellType, constlNCoup, constellations, constellationsNeighType, mult_axis,&
         energy)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: start_atom !< Atom to start loop for
      integer, intent(in) :: stop_atom !< Atom to end loop for
      integer, intent(in) :: do_jtensor   !<  Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
      character(len=1),intent(in) :: exc_inter !< Interpolation of Jij (Y/N)
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_biqdm   !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom  !< Current magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: time_external_field !< External time-dependent magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff1 !< Internal effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff2 !< External field from application of Hamiltonian
      real(dblprec), intent(inout) :: energy !< Total energy

      character(len=1), intent(in) :: mult_axis

      !! +++ Optimization Routines Variables +++ !!
      !! ======================================= !!
      integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
      integer, dimension(:), intent(in) :: maxNoConstl
      integer, dimension(:,:), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup
      real(dblprec), dimension(:,:,:), intent(in) :: constellations
      logical, intent(in) :: OPT_flag
      integer, dimension(:,:,:), intent(in) :: constellationsNeighType
      real(dblprec), dimension(3,max_no_constellations,Mensemble) :: beff1_constellations
      !! +++ End Region +++ !!

      !.. Local scalars
      integer :: i,k
      real(dblprec), dimension(3) :: tfield, beff_s, beff_q


      !.. Executable statements
      if(OPT_flag) call pre_optimize()


      !reduction(+:energy)
      !$omp parallel do default(shared) schedule(static) private(i,k,beff_s,beff_q,tfield) collapse(2) reduction(+:energy)
      do k=1, Mensemble
         do i=start_atom, stop_atom

            beff_s=0.0d0
            beff_q=0.0d0

            if(do_jtensor/=1) then

               ! Heisenberg exchange term
               if(exc_inter=='N') then
                  call heisenberg_field(i,k,beff_s)
               else
                  call heisenberg_rescaling_field(i,k,beff_s)
               endif

               ! Dzyaloshinskii-Moriya term
               if(do_dm==1) call dzyaloshinskii_moriya_field(i, k, beff_s)
            else
               call tensor_field(i, k, beff_s)
            end if

            ! Pseudo-Dipolar term
            if(do_pd==1) call pseudo_dipolar_field(i, k, beff_s)

            ! BIQDM term
            if(do_biqdm==1) call dzyaloshinskii_moriya_bq_field(i, k, beff_s)

            ! Biquadratic exchange term
            if(do_bq==1) call biquadratic_field(i, k, beff_s)

            ! Dipolar term
            if(do_dip==1) call dipolar_field(i, k, beff_s)

            ! Anisotropy
            if (taniso(i)==1) then
               ! Uniaxial anisotropy
               call uniaxial_anisotropy_field(i, k, beff_q)
            elseif (taniso(i)==2) then
               ! Cubic anisotropy
               call cubic_anisotropy_field(i, k, beff_q)
            elseif (taniso(i)==7)then
               ! Uniaxial and cubic anisotropy
               call uniaxial_anisotropy_field(i, k, beff_q)
               tfield=0.0d0
               call cubic_anisotropy_field(i, k, tfield)
               beff_q=beff_q+tfield*sb(i)
            endif


            beff1(1:3,i,k)= beff_s
            beff2(1:3,i,k)= beff_q+external_field(1:3,i,k)+time_external_field(1:3,i,k)
            beff(1:3,i,k) = beff1(1:3,i,k)+beff2(1:3,i,k)

            tfield=0.5d0*(beff_s+1.0d0*beff_q+external_field(1:3,i,k))
            energy=energy - emomM(1,i,k)*tfield(1)-emomM(2,i,k)*tfield(2)-emomM(3,i,k)*tfield(3)
         end do
      end do
      !$omp end parallel do

   contains

      !---------------pre_optimize---------------!
      !> Preoptimization
      subroutine pre_optimize()
         implicit none
         integer :: k,i,j

         beff1_constellations(:,:,:) = 0.0d0
         do k=1,Mensemble
            ! Compute the Heissenberg exchange term for the (representative) unit cells in each constellation
            do i=1,maxNoConstl(k)
               ! Computation of the exchange term, going to max_no_neigh since the exchange factor,
               ! governed by constlNCoup, is zero for iterations beyond the neighbourhood of atom i
               do j=1,max_no_neigh
                  beff1_constellations(:,i,k) = beff1_constellations(:,i,k) + &
                     constlNCoup(j,i,k)*constellations(:,constellationsNeighType(j,i,k),k)
               end do
            end do
         end do

      end subroutine pre_optimize



      !---------------heisenberg_field---------------!
      !> Heisenberg
      subroutine heisenberg_field(i, k, field)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field

         integer :: j

         ! WARNING LSF is not working for ASD ncoup set to first configuration only ASK FAN
         !Exchange term, run optimization or use the standard scheme
         if(OPT_flag) then
            do j=1,nlistsize(i)
               if(unitCellType(i,k).ge.1) then ! If in a constellation
                  field = field+beff1_constellations(:,unitCellType(i,k),k)
               endif
            end do
         else
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
            !$omp simd reduction(+:field)
#endif
            do j=1,nlistsize(i)
               !           !DIR$ vector always aligned
               field = field + ncoup(j,i,1)*emomM(:,nlist(j,i),k)
            end do
         endif
      end subroutine heisenberg_field


      !---------------heisenberg_rescaling_field---------------!
      !> Heisenberg
      subroutine heisenberg_rescaling_field(i, k, field)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field

         integer :: j
         real(dblprec) :: excscale,maxfield
         real(dblprec),dimension(3) :: eff_field

         eff_field=0.d0 ; maxfield=0.d0

#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
         !      !$omp simd reduction(+:eff_field,maxfield)
#endif
         do j=1,nlistsize(i)
            eff_field = eff_field + emomM(:,nlist(j,i),k)
            maxfield = maxfield + mmom(nlist(j,i),k)
         enddo
         excscale=sqrt(eff_field(1)**2+eff_field(2)**2+eff_field(3)**2)/maxfield
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
         !$omp simd reduction(+:field)
#endif
         do j=1,nlistsize(i)
            field = field + ((excscale*ncoup(j,i,1)+(1.d0-excscale)*ncoupD(j,i,1)))*emomM(:,nlist(j,i),k)
         end do
      end subroutine heisenberg_rescaling_field


      !---------------tensor_heisenberg_field---------------!
      !> @brief Calculates heisenberg, DM and anisotropy through one tensor
      !>
      !> @Date 09/15/2014 - Thomas Nystrand
      !> - Added 0.5*m_i*I part
      subroutine tensor_field(i, k, field)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j ! Neighbourlist index
         integer :: x ! Exchange index


         !Exchange term
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
         !$omp simd reduction(+:field)
#endif
         do j=1,nlistsize(i)
            x = nlist(j,i);
            ! Matrix tensor multiplication: f = f+0.5*I*m_i+0.5*m_i*I
            ! Faster then: field = field + 0.5*MATMUL(j_tens(:,:,j,i),emomM(:,x,k)) + 0.5*MATMUL(emomM(:,x,k),j_tens(:,:,j,i))
            !field = field + ncoup(j,i)*emomM(:,x,k)
            field(1) = field(1) &
               + 0.5*(        &
               + j_tens(1,1,j,i)*emomM(1,x,k) + j_tens(1,2,j,i)*emomM(2,x,k) + j_tens(1,3,j,i)*emomM(3,x,k) &
               + emomM(1,x,k)*j_tens(1,1,j,i) + emomM(2,x,k)*j_tens(2,1,j,i) + emomM(3,x,k)*j_tens(3,1,j,i) &
               )
            field(2) = field(2) &
               + 0.5*(        &
               + j_tens(2,1,j,i)*emomM(1,x,k) + j_tens(2,2,j,i)*emomM(2,x,k) + j_tens(2,3,j,i)*emomM(3,x,k) &
               + emomM(1,x,k)*j_tens(1,2,j,i) + emomM(2,x,k)*j_tens(2,2,j,i) + emomM(3,x,k)*j_tens(3,2,j,i) &
               )
            field(3) = field(3) &
               + 0.5*(        &
               + j_tens(3,1,j,i)*emomM(1,x,k) + j_tens(3,2,j,i)*emomM(2,x,k) + j_tens(3,3,j,i)*emomM(3,x,k) &
               + emomM(1,x,k)*j_tens(1,3,j,i) + emomM(2,x,k)*j_tens(2,3,j,i) + emomM(3,x,k)*j_tens(3,3,j,i) &
               )
         end do
         !stop
      end subroutine tensor_field




      !---------------dzyaloshinskii_moriya_field---------------!
      !> DM-field
      subroutine dzyaloshinskii_moriya_field(i, k, field)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j

         ! Dzyaloshinskii_moriya term
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
         !$omp simd reduction(+:field)
#endif
         do j=1,dmlistsize(i)
            field(1) = field(1) - dm_vect(3,j,i)*emomM(2,dmlist(j,i),k) +&
               dm_vect(2,j,i)*emomM(3,dmlist(j,i),k)
            field(2) = field(2) - dm_vect(1,j,i)*emomM(3,dmlist(j,i),k) +&
               dm_vect(3,j,i)*emomM(1,dmlist(j,i),k)
            field(3) = field(3) - dm_vect(2,j,i)*emomM(1,dmlist(j,i),k) +&
               dm_vect(1,j,i)*emomM(2,dmlist(j,i),k)
         end do

      end subroutine dzyaloshinskii_moriya_field



      !---------------pseudo_dipolar_field---------------!
      !> PD-field
      subroutine pseudo_dipolar_field(i, k, field)
         !
         !.. Implicit declarations
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j

         ! Pseudo-Dipolar term
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
         !$omp simd reduction(+:field)
#endif
         do j=1,pdlistsize(i)
            field(1) = field(1) + pd_vect(1,j,i)*emomM(1,pdlist(j,i),k) +&
               pd_vect(4,j,i)*emomM(2,pdlist(j,i),k) +&
               pd_vect(5,j,i)*emomM(3,pdlist(j,i),k)
            field(2) = field(2) + pd_vect(4,j,i)*emomM(1,pdlist(j,i),k) +&
               pd_vect(2,j,i)*emomM(2,pdlist(j,i),k) +&
               pd_vect(6,j,i)*emomM(3,pdlist(j,i),k)
            field(3) = field(3) + pd_vect(5,j,i)*emomM(1,pdlist(j,i),k) +&
               pd_vect(6,j,i)*emomM(2,pdlist(j,i),k) +&
               pd_vect(3,j,i)*emomM(3,pdlist(j,i),k)
         end do

      end subroutine pseudo_dipolar_field




      !---------------dzyaloshinskii_moriya_bq_field---------------!
      !> DM-BQ field
      subroutine dzyaloshinskii_moriya_bq_field(i, k, field)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j
         real(dblprec), dimension(3) :: dot !< Work array

         ! BIQDM term
         do j=1,biqdmlistsize(i)
            dot(1) = emomM(1,i,k)*emomM(2,biqdmlist(j,i),k)-&
               emomM(2,i,k)*emomM(1,biqdmlist(j,i),k)
            dot(2) = emomM(2,i,k)*emomM(3,biqdmlist(j,i),k)-&
               emomM(3,i,k)*emomM(2,biqdmlist(j,i),k)
            dot(3) = emomM(3,i,k)*emomM(1,biqdmlist(j,i),k)-&
               emomM(1,i,k)*emomM(3,biqdmlist(j,i),k)
            field(1) = field(1) + 2.0d0*biqdm_vect(1,j,i)*(&
               dot(1)*emomM(3,biqdmlist(j,i),k)-&
               dot(2)*emomM(2,biqdmlist(j,i),k))
            field(2) = field(2) + 2.0d0*biqdm_vect(1,j,i)*(&
               dot(2)*emomM(1,biqdmlist(j,i),k)-&
               dot(3)*emomM(3,biqdmlist(j,i),k))
            field(3) = field(3) + 2.0d0*biqdm_vect(1,j,i)*(&
               dot(3)*emomM(2,biqdmlist(j,i),k)-&
               dot(1)*emomM(1,biqdmlist(j,i),k))
         end do
      end subroutine dzyaloshinskii_moriya_bq_field



      !---------------biquadratic_field---------------!
      !> BQ-field
      subroutine biquadratic_field(i, k, field)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j
         real(dblprec) :: dot

         ! Biquadratic exchange term
         do j=1,bqlistsize(i)
            dot=emomM(1,bqlist(j,i),k)*emomM(1,i,k)+&
               emomM(2,bqlist(j,i),k)*emomM(2,i,k)+&
               emomM(3,bqlist(j,i),k)*emomM(3,i,k)
            field = field + 2.0d0*j_bq(j,i)*dot*emomM(1:3,bqlist(j,i),k)
         end do
      end subroutine biquadratic_field




      !---------------dipolar_field---------------!
      !> DP-field
      subroutine dipolar_field(i, k, field)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         integer :: j

#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
         !$omp simd reduction(+:field)
#endif
         do j=1,Natom
            field(:) = field(:) + Qdip(1,:,j,i)*emomM(1,j,k) + Qdip(2,:,j,i)*emomM(2,j,k) + Qdip(3,:,j,i)*emomM(3,j,k)
         end do
      end subroutine dipolar_field

      !---------------uniaxial_anisotropy---------------!
      !> Uniaxial anisotropy
      subroutine uniaxial_anisotropy_field(i, k, field)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !
         real(dblprec) :: tt1,tt2,tt3
         real(dblprec) :: tt1_d,tt2_d,tt3_d

         ! Uniaxial anisotropy
         ! cos(theta)
         tt1=emomM(1,i,k)*eaniso(1,i)+emomM(2,i,k)*eaniso(2,i)+emomM(3,i,k)*eaniso(3,i)
         ! k1 + 2*k2*sin^2(theta) = k1 + 2*k2*(1-cos^2(theta))
         tt2=kaniso(1,i)+2.0d0*kaniso(2,i)*(1.0d0-tt1*tt1)
         ! 2 * cos(theta)* [k1 + 2*k2*sin^2(theta)]
         tt3= 2.0d0*tt1*tt2

         if (mult_axis=='Y') then
            ! Uniaxial anisotropy
            ! cos(theta)
            tt1_d=emomM(1,i,k)*eaniso_diff(1,i)+emomM(2,i,k)*eaniso_diff(2,i)+emomM(3,i,k)*eaniso_diff(3,i)
            ! k1 + 2*k2*sin^2(theta) = k1 + 2*k2*(1-cos^2(theta))
            tt2_d=kaniso_diff(1,i)+2.0d0*kaniso_diff(2,i)*(1.0d0-tt1_d*tt1_d)
            ! 2 * cos(theta)* [k1 + 2*k2*sin^2(theta)]
            tt3_d= 2.0d0*tt1_d*tt2_d

            field  = field - tt3*eaniso(1:3,i)-tt3_d*eaniso_diff(1:3,i)
         else

            field  = field - tt3*eaniso(1:3,i)
         endif

      end subroutine uniaxial_anisotropy_field




      !---------------cubic_anisotropy_field---------------!
      !> Cubic anisotropy
      subroutine cubic_anisotropy_field(i, k, field)
         implicit none

         integer, intent(in) :: i !< Atom to calculate effective field for
         integer, intent(in) :: k !< Current ensemble
         real(dblprec), dimension(3), intent(inout) :: field !< Effective field
         !

         field(1) = field(1)  &
            + 2.0d0*kaniso(1,i)*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
            + 2.0d0*kaniso(2,i)*emomM(1,i,k)*(emomM(2,i,k)**2*emomM(3,i,k)**2)
         field(2) = field(2)  &
            + 2.0d0*kaniso(1,i)*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
            + 2.0d0*kaniso(2,i)*emomM(2,i,k)*(emomM(3,i,k)**2*emomM(1,i,k)**2)
         field(3) = field(3)  &
            + 2.0d0*kaniso(1,i)*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
            + 2.0d0*kaniso(2,i)*emomM(3,i,k)*(emomM(1,i,k)**2*emomM(2,i,k)**2)

         if (mult_axis=='Y') then

            field(1) = field(1)  &
               + 2.0d0*kaniso_diff(1,i)*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
               + 2.0d0*kaniso_diff(2,i)*emomM(1,i,k)*emomM(2,i,k)**2*emomM(3,i,k)**2
            field(2) = field(2)  &
               + 2.0d0*kaniso_diff(1,i)*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
               + 2.0d0*kaniso_diff(2,i)*emomM(2,i,k)*emomM(3,i,k)**2*emomM(1,i,k)**2
            field(3) = field(3)  &
               + 2.0d0*kaniso_diff(1,i)*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
               + 2.0d0*kaniso_diff(2,i)*emomM(3,i,k)*emomM(1,i,k)**2*emomM(2,i,k)**2

         endif

      end subroutine cubic_anisotropy_field


   end subroutine effective_field

end module HamiltonianActions
