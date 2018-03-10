!> Calculate effective field by applying the derivative of the Hamiltonian
!> \details The effective field, \f$\mathbf{B}_i\f$, on an atom \f$\textit{i}\f$, is calculated from
!> \f$ \mathbf{B}_i=-\frac{\partial \mathbf{H}}{\partial \mathbf{m}_i},\f$ where primarily the part of
!> the Hamiltonian, \f$\mathbf{H}\f$, which represents interatomic exchange interactions,
!> \f$\mathbf{H}_\mathrm{ex}\f$, are considered. For this we use the classical Heisenberg Hamiltonian,
!> \f$ \mathbf{H}_\mathrm{ex}=-\frac{1}{2}\sum_{i\neq j}J_{ij}\mathbf{m}_i\cdot\mathbf{m}_j,\f$ where
!> \f$i\f$ and \f$j\f$ are atomic indices and \f$J_{ij}\f$ is the strength of the exchange interaction,
!> which is calculated from first principles theory.
!!
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module ApplyHamiltonian
   use Parameters
   use Profiling
   implicit none

   private


contains


   !> Calculate effective field by applying the derivative of the Hamiltonian
   !! @todo Check consistency of terms wrt the input parameters, especially the anisotropies
   !! @todo Replace moment unit vectors emom with full length vectors emomM
   !dm_vect !< DM vector \f$H_{DM}=\sum{D_{ij}\dot(m_i \times m_j)}\f$
   subroutine heisge(Natom, Mensemble, max_no_neigh, ncoup, nlist, nlistsize, &
         do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, &
         do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
         do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, & !for BIQDM
         do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
         taniso, eaniso, kaniso, sb, beff, beff1, beff2, &
         emomM, emom, external_field, do_dip, Qdip)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      real(dblprec), dimension(max_no_neigh,Natom), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
      integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
      integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
      integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
      integer, intent(in) :: do_biqdm   !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange vector
      integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
      integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
      real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
      integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
      integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between the anisotropies
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff1 !< Internal effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff2 !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction
      real(dblprec) :: xu1,xu2

      !.. Local scalars
      integer :: i,j,k,ik
      real(dblprec), dimension(Mensemble) :: tt,tt1,tt2,tt3,ttx,tty,ttz
      real(dblprec), dimension(Mensemble) :: bqmdot
      real(dblprec), dimension(Mensemble) :: sxy,syz,szx

      !.. Executable statements

      !$omp parallel do default(shared) private(ik,i,k,j,ttx,tty,ttz,tt,tt1,tt2,tt3,xu1,xu2) schedule(static)
      do ik=1,Natom*Mensemble
         i=mod(ik-1,Natom)+1
         k=int((ik-1)/Natom)+1

         ttx(k)=0.0d0
         tty(k)=0.0d0
         ttz(k)=0.0d0

         !Exchange term
         do j=1,nlistsize(i)
            !       do k=1,Mensemble
            ttx(k) = ttx(k) + ncoup(j,i)*emomM(1,nlist(j,i),k)
            tty(k) = tty(k) + ncoup(j,i)*emomM(2,nlist(j,i),k)
            ttz(k) = ttz(k) + ncoup(j,i)*emomM(3,nlist(j,i),k)
            !       end do
         end do

         ! Dzyaloshinskii-Moriya term
         if(do_dm==1) then
            do j=1,dmlistsize(i)
               ttx(k) = ttx(k) - dm_vect(3,j,i)*emomM(2,dmlist(j,i),k) +&
                  dm_vect(2,j,i)*emomM(3,dmlist(j,i),k)
               tty(k) = tty(k) - dm_vect(1,j,i)*emomM(3,dmlist(j,i),k) +&
                  dm_vect(3,j,i)*emomM(1,dmlist(j,i),k)
               ttz(k) = ttz(k) - dm_vect(2,j,i)*emomM(1,dmlist(j,i),k) +&
                  dm_vect(1,j,i)*emomM(2,dmlist(j,i),k)
            end do
         end if

         ! Pseudo-Dipolar term
         if(do_pd==1) then
            do j=1,pdlistsize(i)
               ttx(k) = ttx(k) + pd_vect(1,j,i)*emomM(1,pdlist(j,i),k) +&
                  pd_vect(4,j,i)*emomM(2,pdlist(j,i),k) +&
                  pd_vect(5,j,i)*emomM(3,pdlist(j,i),k)
               tty(k) = tty(k) + pd_vect(4,j,i)*emomM(1,pdlist(j,i),k) +&
                  pd_vect(2,j,i)*emomM(2,pdlist(j,i),k) +&
                  pd_vect(6,j,i)*emomM(3,pdlist(j,i),k)
               ttz(k) = ttz(k) + pd_vect(5,j,i)*emomM(1,pdlist(j,i),k) +&
                  pd_vect(6,j,i)*emomM(2,pdlist(j,i),k) +&
                  pd_vect(3,j,i)*emomM(3,pdlist(j,i),k)
            end do
         end if

         ! BIQDM term
         if(do_biqdm==1) then
            do j=1,biqdmlistsize(i)
               sxy(k) = emomM(1,i,k)*emomM(2,biqdmlist(j,i),k)-&
                  emomM(2,i,k)*emomM(1,biqdmlist(j,i),k)
               syz(k) = emomM(2,i,k)*emomM(3,biqdmlist(j,i),k)-&
                  emomM(3,i,k)*emomM(2,biqdmlist(j,i),k)
               szx(k) = emomM(3,i,k)*emomM(1,biqdmlist(j,i),k)-&
                  emomM(1,i,k)*emomM(3,biqdmlist(j,i),k)
               ttx(k) = ttx(k) + 2.0d0*biqdm_vect(1,j,i)*(&
                  szx(k)*emomM(3,biqdmlist(j,i),k)-&
                  sxy(k)*emomM(2,biqdmlist(j,i),k))
               tty(k) = tty(k) + 2.0d0*biqdm_vect(1,j,i)*(&
                  sxy(k)*emomM(1,biqdmlist(j,i),k)-&
                  syz(k)*emomM(3,biqdmlist(j,i),k))
               ttz(k) = ttz(k) + 2.0d0*biqdm_vect(1,j,i)*(&
                  syz(k)*emomM(2,biqdmlist(j,i),k)-&
                  szx(k)*emomM(1,biqdmlist(j,i),k))
            end do
         end if

         ! Biquadratic exchange term
         if(do_bq==1) then
            do j=1,bqlistsize(i)
               bqmdot(k)=emomM(1,bqlist(j,i),k)*emomM(1,i,k)+&
                  emomM(2,bqlist(j,i),k)*emomM(2,i,k)+&
                  emomM(3,bqlist(j,i),k)*emomM(3,i,k)
               ttx(k) = ttx(k)+ 2.0d0*j_bq(j,i)*bqmdot(k)*emomM(1,bqlist(j,i),k)
               tty(k) = tty(k)+ 2.0d0*j_bq(j,i)*bqmdot(k)*emomM(2,bqlist(j,i),k)
               ttz(k) = ttz(k)+ 2.0d0*j_bq(j,i)*bqmdot(k)*emomM(3,bqlist(j,i),k)
            end do
         end if

         ! Dipolar term
         if(present(Qdip)) then
            if(do_dip==1) then
               do j=1,Natom
                  ttx(k) = ttx(k) + Qdip(1,1,j,i)*emomM(1,j,k) + Qdip(2,1,j,i)*emomM(2,j,k) + Qdip(3,1,j,i)*emomM(3,j,k)
                  tty(k) = tty(k) + Qdip(1,2,j,i)*emomM(1,j,k) + Qdip(2,2,j,i)*emomM(2,j,k) + Qdip(3,2,j,i)*emomM(3,j,k)
                  ttz(k) = ttz(k) + Qdip(1,3,j,i)*emomM(1,j,k) + Qdip(2,3,j,i)*emomM(2,j,k) + Qdip(3,3,j,i)*emomM(3,j,k)
               end do
            end if
         end if

         ! Sum effective fields
         beff1(1,i,k) = ttx(k)
         beff1(2,i,k) = tty(k)
         beff1(3,i,k) = ttz(k)
         ttx(k) = 0d0
         tty(k) = 0d0
         ttz(k) = 0d0

         ! ! Add anisotropies and external fields
         ! Anisotropy
         if (taniso(i)==1) then
            ! Uniaxial anisotropy
            tt1(k)=emomM(1,i,k)*eaniso(1,i)+emomM(2,i,k)*eaniso(2,i)+emomM(3,i,k)*eaniso(3,i)

            tt2(k)=kaniso(1,i)+2.0d0*kaniso(2,i)*(1-tt1(k)*tt1(k))
            !
            tt3(k)= 2.0d0*tt1(k)*tt2(k)

            ttx(k)  = ttx(k) - tt3(k)*eaniso(1,i)
            tty(k)  = tty(k) - tt3(k)*eaniso(2,i)
            ttz(k)  = ttz(k) - tt3(k)*eaniso(3,i)

         elseif (taniso(i)==2) then
            ! Cubic anisotropy
            ttx(k) = ttx(k)  &
               + 2.0d0*kaniso(1,i)*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
               + 2.0d0*kaniso(2,i)*emomM(1,i,k)*emomM(2,i,k)**2*emomM(3,i,k)**2
            tty(k) = tty(k)  &
               + 2.0d0*kaniso(1,i)*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
               + 2.0d0*kaniso(2,i)*emomM(2,i,k)*emomM(3,i,k)**2*emomM(1,i,k)**2
            ttz(k) = ttz(k)  &
               + 2.0d0*kaniso(1,i)*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
               + 2.0d0*kaniso(2,i)*emomM(3,i,k)*emomM(1,i,k)**2*emomM(2,i,k)**2
         elseif (taniso(i)==3) then
            tt1(k)=emom(1,i,k)*eaniso(1,i)+emom(2,i,k)*eaniso(2,i)+emom(3,i,k)*eaniso(3,i)
            ! only for -45 < theta < 45 and 135<theta<225
            if(abs(tt1(k))>0.70710678118654d0) then
               tt2(k)=16.0d0*tt1(k)**3-8.0d0*tt1(k)
               ! *K
               tt3(k)=kaniso(1,i)*tt2(k)
               !
               ttx(k)  = ttx(k) - tt3(k)*eaniso(1,i)
               tty(k)  = tty(k) - tt3(k)*eaniso(2,i)
               ttz(k)  = ttz(k) - tt3(k)*eaniso(3,i)
            end if
         elseif (taniso(i)==7)then
            !    Uniaxial
            ! K1*(sin theta)^2+K2*(sin theta)^4
            ! cos(theta)
            tt1(k)=emomM(1,i,k)*eaniso(1,i)+emomM(2,i,k)*eaniso(2,i)+emomM(3,i,k)*eaniso(3,i)
            ! cos(theta)/m_i
            ! tt1(k)=mmomi(i,k)*tt(k)
            ! K1+2*K2*(m.ek)^2
            tt2(k)=kaniso(1,i)+2.0d0*kaniso(2,i)*tt1(k)**2
            !
            tt3(k)=2.0d0*tt1(k)*tt2(k)

            ttx(k)  = ttx(k) - tt3(k)*eaniso(1,i)
            tty(k)  = tty(k) - tt3(k)*eaniso(2,i)
            ttz(k)  = ttz(k) - tt3(k)*eaniso(3,i)

            !    +  Cubic
            !     The Cubic Anisotropy constant = Uniaxial constant x sb
            xu1=kaniso(1,i)*sb(i)
            xu2=kaniso(2,i)*sb(i)
            !
            ttx(k) = ttx(k)  &
               + 2.0d0*xu1*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
               + 2.0d0*xu2*emomM(1,i,k)*emomM(2,i,k)**2*emomM(3,i,k)**2

            tty(k) = tty(k)  &
               + 2.0d0*xu1*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
               + 2.0d0*xu2*emomM(2,i,k)*emomM(3,i,k)**2*emomM(1,i,k)**2

            ttz(k) = ttz(k)  &
               + 2.0d0*xu1*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
               + 2.0d0*xu2*emomM(3,i,k)*emomM(1,i,k)**2*emomM(2,i,k)**2
         endif

         ! Add static global field
         beff2(1,i,k) = ttx(k)+external_field(1,i,k)
         beff2(2,i,k) = tty(k)+external_field(2,i,k)
         beff2(3,i,k) = ttz(k)+external_field(3,i,k)
         beff(1,i,k) = beff1(1,i,k)+beff2(1,i,k)
         beff(2,i,k) = beff1(2,i,k)+beff2(2,i,k)
         beff(3,i,k) = beff1(3,i,k)+beff2(3,i,k)
      end do
      !$omp end parallel do

   end subroutine heisge


   !> Calculate effective field by applying the derivative of the Hamiltonian
   ! using tensorial (SKKR) exchange.
   !! @todo Check consistency of terms wrt the input parameters, especially the anisotropies
   !! @todo Replace moment unit vectors emom with full length vectors emomM
   !dm_vect !< DM vector \f$H_{DM}=\sum{D_{ij}\dot(m_i \times m_j)}\f$
   subroutine heisge_tens(Natom, Mensemble, max_no_neigh, j_tens, nlist, nlistsize, &
         do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, &
         do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
         taniso, eaniso, kaniso, sb, beff, beff1, beff2, &
         emomM, external_field)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      real(dblprec), dimension(3,3,max_no_neigh,Natom), intent(in) :: j_tens !< Tensorial exchange couplings (SKKR)
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, intent(in) :: do_biqdm   !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange vector
      integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
      integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
      real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
      integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
      integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between the anisotropies (Cubic and Uniaxial)
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff1 !< Internal effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff2 !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field
      real(dblprec) :: xu1,xu2

      !.. Local scalars
      integer :: i,j,k
      real(dblprec), dimension(Mensemble) :: tt,tt1,tt2,tt3,ttx,tty,ttz
      real(dblprec), dimension(Mensemble) :: bqmdot
      real(dblprec), dimension(Mensemble) :: sxy,syz,szx

      !.. Executable statements


      ! Add contribution from exchange
      !$omp parallel do default(shared) private(i,j,k,ttx,tty,ttz,tt,tt1,tt2,tt3,xu1,xu2) schedule(static)
      do i=1, Natom
         ttx=0.0d0
         tty=0.0d0
         ttz=0.0d0
         !Exchange term
         do j=1,nlistsize(i)
            do k=1,Mensemble
               ttx(k) = ttx(k) + j_tens(1,1,j,i)*emomM(1,nlist(j,i),k) + j_tens(1,2,j,i)*emomM(2,nlist(j,i),k) &
                  + j_tens(1,3,j,i)*emomM(3,nlist(j,i),k)
               tty(k) = tty(k) + j_tens(2,1,j,i)*emomM(1,nlist(j,i),k) + j_tens(2,2,j,i)*emomM(2,nlist(j,i),k) &
                  + j_tens(2,3,j,i)*emomM(3,nlist(j,i),k)
               ttz(k) = ttz(k) + j_tens(3,1,j,i)*emomM(1,nlist(j,i),k) + j_tens(3,2,j,i)*emomM(2,nlist(j,i),k) &
                  + j_tens(3,3,j,i)*emomM(3,nlist(j,i),k)
            end do
         end do

         ! BIQDM term
         if(do_biqdm==1) then
            do j=1,biqdmlistsize(i)
               do k=1,Mensemble
                  sxy(k) = emomM(1,i,k)*emomM(2,biqdmlist(j,i),k)-&
                     emomM(2,i,k)*emomM(1,biqdmlist(j,i),k)
                  syz(k) = emomM(2,i,k)*emomM(3,biqdmlist(j,i),k)-&
                     emomM(3,i,k)*emomM(2,biqdmlist(j,i),k)
                  szx(k) = emomM(3,i,k)*emomM(1,biqdmlist(j,i),k)-&
                     emomM(1,i,k)*emomM(3,biqdmlist(j,i),k)
                  ttx(k) = ttx(k) + 2.0d0*biqdm_vect(1,j,i)*(&
                     szx(k)*emomM(3,biqdmlist(j,i),k)-&
                     sxy(k)*emomM(2,biqdmlist(j,i),k))
                  tty(k) = tty(k) + 2.0d0*biqdm_vect(1,j,i)*(&
                     sxy(k)*emomM(1,biqdmlist(j,i),k)-&
                     syz(k)*emomM(3,biqdmlist(j,i),k))
                  ttz(k) = ttz(k) + 2.0d0*biqdm_vect(1,j,i)*(&
                     syz(k)*emomM(2,biqdmlist(j,i),k)-&
                     szx(k)*emomM(1,biqdmlist(j,i),k))
               end do
            end do
         end if

         ! Biquadratic exchange term
         if(do_bq==1) then
            do j=1,bqlistsize(i)
               do k=1,Mensemble
                  bqmdot(k)=emomM(1,bqlist(j,i),k)*emomM(1,i,k)+&
                     emomM(2,bqlist(j,i),k)*emomM(2,i,k)+&
                     emomM(3,bqlist(j,i),k)*emomM(3,i,k)
                  ttx(k) = ttx(k) + 2.0d0*j_bq(j,i)*bqmdot(k)*emomM(1,bqlist(j,i),k)
                  tty(k) = tty(k) + 2.0d0*j_bq(j,i)*bqmdot(k)*emomM(2,bqlist(j,i),k)
                  ttz(k) = ttz(k) + 2.0d0*j_bq(j,i)*bqmdot(k)*emomM(3,bqlist(j,i),k)
               end do
            end do
         end if

         ! Sum effective fields
         do k=1,Mensemble
            beff1(1,i,k) = ttx(k)
            beff1(2,i,k) = tty(k)
            beff1(3,i,k) = ttz(k)
            ttx(k) = 0d0
            tty(k) = 0d0
            ttz(k) = 0d0
         end do
      end do
      !$omp end parallel do

      ! Add anisotropies and external fields
      !$omp parallel do default(shared) private(i,j,k,ttx,tty,ttz,tt,tt1,tt2,tt3) schedule(static)
      do i=1, Natom
         ttx = 0d0
         tty = 0d0
         ttz = 0d0

         ! Anisotropy
         if (taniso(i)==1) then
            ! Uniaxial anisotropy
            do k=1,Mensemble
               ! K1*(m.ek)^2+K2*(m.ek)^4
               !
               ! cos(theta)
               tt1(k)=emomM(1,i,k)*eaniso(1,i)+emomM(2,i,k)*eaniso(2,i)+emomM(3,i,k)*eaniso(3,i)
               ! K1+2*K2*(1-cos2(theta))
               tt2(k)=kaniso(1,i)+2.0d0*kaniso(2,i)*(1-tt(k)*tt(k))
               !
               tt3(k)=2.0d0*tt(k)*tt2(k)

               ttx(k)  = ttx(k) - tt3(k)*eaniso(1,i)
               tty(k)  = tty(k) - tt3(k)*eaniso(2,i)
               ttz(k)  = ttz(k) - tt3(k)*eaniso(3,i)

            end do
         elseif (taniso(i)==2) then
            ! Cubic anisotropy
            do k=1,Mensemble
               ttx(k) = ttx(k)  &
                  + 2.0d0*kaniso(1,i)*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
                  + 2.0d0*kaniso(2,i)+emomM(1,i,k)*emomM(2,i,k)**2*emomM(3,i,k)**2

               tty(k) = tty(k)  &
                  + 2.0d0*kaniso(1,i)*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
                  + 2.0d0*kaniso(2,i)+emomM(2,i,k)*emomM(3,i,k)**2*emomM(1,i,k)**2

               ttz(k) = ttz(k)  &
                  + 2.0d0*kaniso(1,i)*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
                  + 2.0d0*kaniso(2,i)+emomM(3,i,k)*emomM(1,i,k)**2*emomM(2,i,k)**2

            end do
            !!     endif
         elseif (taniso(i)==7)then
            !    Uniaxial
            do k=1,Mensemble
               ! K1*(m.ek)^2+K2*(m.ek)^4
               !
               ! cos(theta)
               tt1(k)=emomM(1,i,k)*eaniso(1,i)+emomM(2,i,k)*eaniso(2,i)+emomM(3,i,k)*eaniso(3,i)
               ! K1+2*K2*(1-cos2(theta))
               tt2(k)=kaniso(1,i)+2.0d0*kaniso(2,i)*(1-tt(k)*tt(k))
               !
               tt3(k)=2.0d0*tt(k)*tt2(k)

               ttx(k)  = ttx(k) - tt3(k)*eaniso(1,i)
               tty(k)  = tty(k) - tt3(k)*eaniso(2,i)
               ttz(k)  = ttz(k) - tt3(k)*eaniso(3,i)

            end do
            !    +  Cubic
            !     The Cubic Anisotropy constant = Uniaxial constant x sb
            xu1=kaniso(1,i)*sb(i)
            xu2=kaniso(2,i)*sb(i)
            !
            do k=1,Mensemble
               ttx(k) = ttx(k)  &
                  + 2.0d0*xu1*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
                  + 2.0d0*xu2+emomM(1,i,k)*emomM(2,i,k)**2*emomM(3,i,k)**2

               tty(k) = tty(k)  &
                  + 2.0d0*xu1*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
                  + 2.0d0*xu2+emomM(2,i,k)*emomM(3,i,k)**2*emomM(1,i,k)**2

               ttz(k) = ttz(k)  &
                  + 2.0d0*xu1*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
                  + 2.0d0*xu2+emomM(3,i,k)*emomM(1,i,k)**2*emomM(2,i,k)**2

            end do
         endif
         ! Add static global field
         do k=1,Mensemble
            beff2(1,i,k) = ttx(k)+external_field(1,i,k)
            beff2(2,i,k) = tty(k)+external_field(2,i,k)
            beff2(3,i,k) = ttz(k)+external_field(3,i,k)
            beff(1,i,k) = beff1(1,i,k)+beff2(1,i,k)
            beff(2,i,k) = beff1(2,i,k)+beff2(2,i,k)
            beff(3,i,k) = beff1(3,i,k)+beff2(3,i,k)
         end do
      end do
      !$omp end parallel do

   end subroutine heisge_tens

   !> Calculate effective field by applying the derivative of the Hamiltonian
   !! @todo Check consistency of terms wrt the input parameters, especially the anisotropies
   !! @todo Replace moment unit vectors emom with full length vectors emomM
   !dm_vect !< DM vector \f$H_{DM}=\sum{D_{ij}\dot(m_i \times m_j)}\f$
   subroutine heisge_jij(Natom, Mensemble, max_no_neigh, ncoup, nlist, nlistsize, &
         beff, emomM, emom, external_field )
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      real(dblprec), dimension(max_no_neigh,Natom), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field

      !.. Local scalars
      integer :: i,j,k,ik
      real(dblprec) :: xu1,xu2
      real(dblprec), dimension(Mensemble) :: tt,tt1,tt2,tt3,ttx,tty,ttz
      real(dblprec), dimension(Mensemble) :: bqmdot
      real(dblprec), dimension(Mensemble) :: sxy,syz,szx

      !.. Executable statements

      !$omp parallel do default(shared) private(ik,i,k,j,ttx,tty,ttz) schedule(static)
      do ik=1,Natom*Mensemble
         i=mod(ik-1,Natom)+1
         k=int((ik-1)/Natom)+1

         ttx(k)=0.0d0
         tty(k)=0.0d0
         ttz(k)=0.0d0

         !Exchange term
         do j=1,nlistsize(i)
            !       do k=1,Mensemble
            ttx(k) = ttx(k) + ncoup(j,i)*emomM(1,nlist(j,i),k)
            tty(k) = tty(k) + ncoup(j,i)*emomM(2,nlist(j,i),k)
            ttz(k) = ttz(k) + ncoup(j,i)*emomM(3,nlist(j,i),k)
            !       end do
         end do

         ! Sum effective fields
         beff(1,i,k)  = ttx(k) + external_field(1,i,k)
         beff(2,i,k)  = tty(k) + external_field(2,i,k)
         beff(3,i,k)  = ttz(k) + external_field(3,i,k)


      end do
      !$omp end parallel do

   end subroutine heisge_jij

   !> Calculate effective field by applying the derivative of the Hamiltonian
   !! @todo Check consistency of terms wrt the input parameters, especially the anisotropies
   !! @todo Replace moment unit vectors emom with full length vectors emomM
   !dm_vect !< DM vector \f$H_{DM}=\sum{D_{ij}\dot(m_i \times m_j)}\f$
   subroutine heisge_jij_improve(Natom, Mensemble, max_no_neigh, ncoup, nlist, nlistsize, &
         beff, emomM, emom, external_field )
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      real(dblprec), dimension(max_no_neigh,Natom), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field

      !.. Local scalars
      integer :: i,j,k,ik
      real(dblprec) :: tx,ty,tz

      integer :: n,nls
      real(dblprec) :: nc

      !.. Executable statements

      !$omp parallel do default(shared) private(ik,i,k,j,tx,ty,tz,nls,n,nc) schedule(static)
      do ik=1,Natom*Mensemble
         i=mod(ik-1,Natom)+1
         k=int((ik-1)/Natom)+1

         tx = 0.0d0
         ty = 0.0d0
         tz = 0.0d0

         !Exchange term
         nls=nlistsize(i)
         do j=1,nls
            nc = ncoup(j,i)
            n  = nlist(j,i)
            tx = tx + nc*emomM(1,n,k)
            ty = ty + nc*emomM(2,n,k)
            tz = tz + nc*emomM(3,n,k)
         end do

         ! Sum effective fields
         beff(1,i,k)  = tx + external_field(1,i,k)
         beff(2,i,k)  = ty + external_field(2,i,k)
         beff(3,i,k)  = tz + external_field(3,i,k)

      end do
      !$omp end parallel do

   end subroutine heisge_jij_improve


   !> Calculate effective field for SMP by applying the derivative of the Hamiltonian
   !! @todo At first only for Heisenberg exchange, fix also for other terms
   subroutine heisge_sph(Natom, Mensemble, max_no_neigh, ncoup, nlist, nlistsize, &
         do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, &
         do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
         do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, & !for BIQDM
         do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
         taniso, eaniso, kaniso, sb, beff, beff1, beff2, &
         emomM, emom, external_field, do_dip, Qdip)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      real(dblprec), dimension(max_no_neigh,Natom), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
      integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
      integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
      integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
      integer, intent(in) :: do_biqdm   !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange vector
      integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
      integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
      real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
      integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
      integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between the anisotropies
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff1 !< Internal effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff2 !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction
      real(dblprec) :: xu1,xu2

      !.. Local scalars
      integer :: i,j,k,ik
      real(dblprec), dimension(Mensemble) :: tt,tt1,tt2,tt3,ttx,tty,ttz
      real(dblprec), dimension(Mensemble) :: bqmdot
      real(dblprec), dimension(Mensemble) :: sxy,syz,szx

      !.. Executable statements

      !$omp parallel do default(shared) private(ik,i,k,j,ttx,tty,ttz,tt,tt1,tt2,tt3,xu1,xu2) schedule(static)
      do ik=1,Natom*Mensemble
         i=mod(ik-1,Natom)+1
         k=int((ik-1)/Natom)+1

         ttx(k)=0.0d0
         tty(k)=0.0d0
         ttz(k)=0.0d0

         !Exchange term
         do j=1,nlistsize(i)
            ttx(k) = ttx(k) + ncoup(j,i)*emomM(1,nlist(j,i),k)
            tty(k) = tty(k) + ncoup(j,i)*emomM(2,nlist(j,i),k)
            ttz(k) = ttz(k) + ncoup(j,i)*emomM(3,nlist(j,i),k)
         end do

         ! Dzyaloshinskii-Moriya term
         if(do_dm==1) then
            do j=1,dmlistsize(i)
               ttx(k) = ttx(k) - dm_vect(3,j,i)*emomM(2,dmlist(j,i),k) +&
                  dm_vect(2,j,i)*emomM(3,dmlist(j,i),k)
               tty(k) = tty(k) - dm_vect(1,j,i)*emomM(3,dmlist(j,i),k) +&
                  dm_vect(3,j,i)*emomM(1,dmlist(j,i),k)
               ttz(k) = ttz(k) - dm_vect(2,j,i)*emomM(1,dmlist(j,i),k) +&
                  dm_vect(1,j,i)*emomM(2,dmlist(j,i),k)
            end do
         end if

         ! Pseudo-Dipolar term
         if(do_pd==1) then
            do j=1,pdlistsize(i)
               ttx(k) = ttx(k) + pd_vect(1,j,i)*emomM(1,pdlist(j,i),k) +&
                  pd_vect(4,j,i)*emomM(2,pdlist(j,i),k) +&
                  pd_vect(5,j,i)*emomM(3,pdlist(j,i),k)
               tty(k) = tty(k) + pd_vect(4,j,i)*emomM(1,pdlist(j,i),k) +&
                  pd_vect(2,j,i)*emomM(2,pdlist(j,i),k) +&
                  pd_vect(6,j,i)*emomM(3,pdlist(j,i),k)
               ttz(k) = ttz(k) + pd_vect(5,j,i)*emomM(1,pdlist(j,i),k) +&
                  pd_vect(6,j,i)*emomM(2,pdlist(j,i),k) +&
                  pd_vect(3,j,i)*emomM(3,pdlist(j,i),k)
            end do
         end if

         ! BIQDM term
         if(do_biqdm==1) then
            do j=1,biqdmlistsize(i)
               sxy(k) = emomM(1,i,k)*emomM(2,biqdmlist(j,i),k)-&
                  emomM(2,i,k)*emomM(1,biqdmlist(j,i),k)
               syz(k) = emomM(2,i,k)*emomM(3,biqdmlist(j,i),k)-&
                  emomM(3,i,k)*emomM(2,biqdmlist(j,i),k)
               szx(k) = emomM(3,i,k)*emomM(1,biqdmlist(j,i),k)-&
                  emomM(1,i,k)*emomM(3,biqdmlist(j,i),k)
               ttx(k) = ttx(k) + 2.0d0*biqdm_vect(1,j,i)*(&
                  szx(k)*emomM(3,biqdmlist(j,i),k)-&
                  sxy(k)*emomM(2,biqdmlist(j,i),k))
               tty(k) = tty(k) + 2.0d0*biqdm_vect(1,j,i)*(&
                  sxy(k)*emomM(1,biqdmlist(j,i),k)-&
                  syz(k)*emomM(3,biqdmlist(j,i),k))
               ttz(k) = ttz(k) + 2.0d0*biqdm_vect(1,j,i)*(&
                  syz(k)*emomM(2,biqdmlist(j,i),k)-&
                  szx(k)*emomM(1,biqdmlist(j,i),k))
            end do
         end if

         ! Biquadratic exchange term
         if(do_bq==1) then
            do j=1,bqlistsize(i)
               bqmdot(k)=emomM(1,bqlist(j,i),k)*emomM(1,i,k)+&
                  emomM(2,bqlist(j,i),k)*emomM(2,i,k)+&
                  emomM(3,bqlist(j,i),k)*emomM(3,i,k)
               ttx(k) = ttx(k)+ 2.0d0*j_bq(j,i)*bqmdot(k)*emomM(1,bqlist(j,i),k)
               tty(k) = tty(k)+ 2.0d0*j_bq(j,i)*bqmdot(k)*emomM(2,bqlist(j,i),k)
               ttz(k) = ttz(k)+ 2.0d0*j_bq(j,i)*bqmdot(k)*emomM(3,bqlist(j,i),k)
            end do
         end if

         ! Dipolar term
         if(present(Qdip)) then
            if(do_dip==1) then
               do j=1,Natom
                  ttx(k) = ttx(k) + Qdip(1,1,j,i)*emomM(1,j,k) + Qdip(2,1,j,i)*emomM(2,j,k) + Qdip(3,1,j,i)*emomM(3,j,k)
                  tty(k) = tty(k) + Qdip(1,2,j,i)*emomM(1,j,k) + Qdip(2,2,j,i)*emomM(2,j,k) + Qdip(3,2,j,i)*emomM(3,j,k)
                  ttz(k) = ttz(k) + Qdip(1,3,j,i)*emomM(1,j,k) + Qdip(2,3,j,i)*emomM(2,j,k) + Qdip(3,3,j,i)*emomM(3,j,k)
               end do
            end if
         end if

         ! Sum effective fields
         beff1(1,i,k) = ttx(k)
         beff1(2,i,k) = tty(k)
         beff1(3,i,k) = ttz(k)
         ttx(k) = 0d0
         tty(k) = 0d0
         ttz(k) = 0d0

         ! ! Add anisotropies and external fields
         ! Anisotropy
         if (taniso(i)==1) then
            ! Uniaxial anisotropy
            tt1(k)=emomM(1,i,k)*eaniso(1,i)+emomM(2,i,k)*eaniso(2,i)+emomM(3,i,k)*eaniso(3,i)

            tt2(k)=kaniso(1,i)+2.0d0*kaniso(2,i)*(1-tt1(k)*tt1(k))
            !
            tt3(k)= 2.0d0*tt1(k)*tt2(k)

            ttx(k)  = ttx(k) - tt3(k)*eaniso(1,i)
            tty(k)  = tty(k) - tt3(k)*eaniso(2,i)
            ttz(k)  = ttz(k) - tt3(k)*eaniso(3,i)

         elseif (taniso(i)==2) then
            ! Cubic anisotropy
            ttx(k) = ttx(k)  &
               + 2.0d0*kaniso(1,i)*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
               + 2.0d0*kaniso(2,i)*emomM(1,i,k)*emomM(2,i,k)**2*emomM(3,i,k)**2
            tty(k) = tty(k)  &
               + 2.0d0*kaniso(1,i)*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
               + 2.0d0*kaniso(2,i)*emomM(2,i,k)*emomM(3,i,k)**2*emomM(1,i,k)**2
            ttz(k) = ttz(k)  &
               + 2.0d0*kaniso(1,i)*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
               + 2.0d0*kaniso(2,i)*emomM(3,i,k)*emomM(1,i,k)**2*emomM(2,i,k)**2
         elseif (taniso(i)==3) then
            tt1(k)=emom(1,i,k)*eaniso(1,i)+emom(2,i,k)*eaniso(2,i)+emom(3,i,k)*eaniso(3,i)
            ! only for -45 < theta < 45 and 135<theta<225
            if(abs(tt1(k))>0.70710678118654d0) then
               tt2(k)=16.0d0*tt1(k)**3-8.0d0*tt1(k)
               ! *K
               tt3(k)=kaniso(1,i)*tt2(k)
               !
               ttx(k)  = ttx(k) - tt3(k)*eaniso(1,i)
               tty(k)  = tty(k) - tt3(k)*eaniso(2,i)
               ttz(k)  = ttz(k) - tt3(k)*eaniso(3,i)
            end if
         elseif (taniso(i)==7)then
            !    Uniaxial
            ! K1*(sin theta)^2+K2*(sin theta)^4
            ! cos(theta)
            tt1(k)=emomM(1,i,k)*eaniso(1,i)+emomM(2,i,k)*eaniso(2,i)+emomM(3,i,k)*eaniso(3,i)
            ! cos(theta)/m_i
            ! tt1(k)=mmomi(i,k)*tt(k)
            ! K1+2*K2*(m.ek)^2
            tt2(k)=kaniso(1,i)+2.0d0*kaniso(2,i)*tt1(k)**2
            !
            tt3(k)=2.0d0*tt1(k)*tt2(k)

            ttx(k)  = ttx(k) - tt3(k)*eaniso(1,i)
            tty(k)  = tty(k) - tt3(k)*eaniso(2,i)
            ttz(k)  = ttz(k) - tt3(k)*eaniso(3,i)

            !    +  Cubic
            !     The Cubic Anisotropy constant = Uniaxial constant x sb
            xu1=kaniso(1,i)*sb(i)
            xu2=kaniso(2,i)*sb(i)
            !
            ttx(k) = ttx(k)  &
               + 2.0d0*xu1*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
               + 2.0d0*xu2*emomM(1,i,k)*emomM(2,i,k)**2*emomM(3,i,k)**2

            tty(k) = tty(k)  &
               + 2.0d0*xu1*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
               + 2.0d0*xu2*emomM(2,i,k)*emomM(3,i,k)**2*emomM(1,i,k)**2

            ttz(k) = ttz(k)  &
               + 2.0d0*xu1*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
               + 2.0d0*xu2*emomM(3,i,k)*emomM(1,i,k)**2*emomM(2,i,k)**2
         endif

         ! Add static global field
         beff2(1,i,k) = ttx(k)+external_field(1,i,k)
         beff2(2,i,k) = tty(k)+external_field(2,i,k)
         beff2(3,i,k) = ttz(k)+external_field(3,i,k)
         beff(1,i,k) = beff1(1,i,k)+beff2(1,i,k)
         beff(2,i,k) = beff1(2,i,k)+beff2(2,i,k)
         beff(3,i,k) = beff1(3,i,k)+beff2(3,i,k)
      end do
      !$omp end parallel do

   end subroutine heisge_sph


end module ApplyHamiltonian
