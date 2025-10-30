!-------------------------------------------------------------------------------
! MODULE: GNEB
!> Routines for the geodesic nudged elastic band (GNEB) method
!> @todo One can change the im_ene arrays so that the ones from the new energy routine are used
!> @author Pavel Bessarab
!> @updated by Nikos Ntallis
!> @im_ene array has been removed. After the minimazation the ene structure holds the energy per image. 
!> @if plotenergy flag is used energies/atom are printed as in the main code.
!> 
!-------------------------------------------------------------------------------
module GNEB
   use Parameters
   use HamiltonianActions
   use VPO
   use Energy, only: ene, calc_energy,allocate_energies
   use Pathinit, only: save_path
   use Constants

   implicit none

   ! .. Global allocatable arrays
   real(dblprec), dimension(:), allocatable :: u !< Total energy of the sample per image
   real(dblprec), dimension(:), allocatable :: fpp
   real(dblprec), dimension(:), allocatable :: pathl
   real(dblprec), dimension(:,:), allocatable :: coo
   real(dblprec), dimension(:,:), allocatable :: ang
   real(dblprec), dimension(:,:), allocatable :: tau
   real(dblprec), dimension(:,:), allocatable :: gsp
   real(dblprec), dimension(:,:), allocatable :: tau_i
   real(dblprec), dimension(:,:,:), allocatable :: ax
   real(dblprec), dimension(:,:,:), allocatable :: vel
   real(dblprec), dimension(:,:,:), allocatable :: beff !< Total effective field from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: beff0
   real(dblprec), dimension(:,:,:), allocatable :: beff1 !< Internal effective field from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: beff2 !< External field from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: time_external_field !< External time-dependent magnetic field, zero here

   private

   public :: save_en, gneb_mep, gneb_ci_mep, find_sp, find_sp_conf

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: gneb_mep
   !> Evolves the path to the nearest MEP according to the GNEB+VPO algorithm
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine gneb_mep(nim,nHam,Natom,every,do_dm,do_pd,do_bq,do_ring,do_chir,do_dip,itrmax,&
      Nchmax,conf_num,do_biqdm,Num_macro,do_jtensor,plotenergy,do_anisotropy,       &
      max_no_constellations,mass,ftol,kappa,delta_t,simid,do_lsf,en_zero,fixed_if,  &
      mult_axis,exc_inter,lsf_field,lsf_interpolate,OPT_flag,cell_index,            &
      macro_nlistsize,mmom,emom,emomM_macro,external_field,maxNoConstl,unitCellType,&
      constellationsNeighType,constlNCoup,constellations,tenergy,emomM,rx,&
      dene,NA,N1,N2,N3,mode,do_mom_legacy)

      implicit none

      ! .. Input variables
      integer, intent(in) :: NA             !< Number of atoms in one cell
      integer, intent(in) :: N1             !< Number of cell repetitions in x direction
      integer, intent(in) :: N2             !< Number of cell repetitions in y direction
      integer, intent(in) :: N3             !< Number of cell repetitions in z direction
      integer, intent(in) :: nim            !< Number of images in GNEB calculations
      integer, intent(in) :: nHam           !< Number of atoms in Hamiltonian
      integer, intent(in) :: Natom          !< Number of atoms in system
      integer, intent(in) :: every          !< Save path every 'every' step
      integer, intent(in) :: do_dm          !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd          !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq          !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_ring        !< Add four-spin ring (4SR) term to Hamiltonian (0/1)      
      integer, intent(in) :: do_dip         !< Calculate dipole-dipole contribution (0/1)
      integer, intent(in) :: itrmax         !< Maximum number of iterations
      integer, intent(in) :: Nchmax         !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_chir        !< Add scalar chirality exchange (CHIR) term to Hamiltonian (0/1)
      integer, intent(in) :: conf_num       !< number of configurations for LSF
      integer, intent(in) :: do_biqdm       !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: Num_macro      !< Number of macrocells in the system
      integer, intent(in) :: do_jtensor     !< Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
      integer, intent(in) :: plotenergy     !< Calculate and plot energy (0/1)
      integer, intent(in) :: do_anisotropy  !< Read anisotropy data (1/0)
      integer, intent(in) :: max_no_constellations !< The maximum (global) length of the constellation matrix
      real(dblprec), intent(in) :: mass      !< mass of the point
      real(dblprec), intent(in) :: ftol      !< Tolerance
      real(dblprec), intent(in) :: kappa     !< spring constant
      real(dblprec), intent(in) :: delta_t   !< timestep
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid     !< Name of simulation
      character(len=1), intent(in) :: do_lsf    !< Including LSF energy
      character(len=1), intent(in) :: en_zero   !< Level of zero energy
      character(len=1), intent(in) :: fixed_if  !< Fix endpoints of the path (Y) or allow them to move along energy isocontours (N)
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy
      character(len=1), intent(in) :: exc_inter !< Interpolation of Jij (Y/N)
      character(len=1), intent(in) :: lsf_field       !< LSF field contribution (Local/Total)
      character(len=1), intent(in) :: do_mom_legacy   !< Flag to print/read moments in legacy output
      character(len=1), intent(in) :: lsf_interpolate !< Interpolate LSF or not
      logical, intent(in) :: OPT_flag  !< Optimization flag
      integer, dimension(Natom), intent(in) :: cell_index    !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      real(dblprec), dimension(Natom,nim), intent(in) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Num_macro,nim), intent(in) :: emomM_macro    !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,Natom,nim), intent(in)     :: external_field !< External magnetic field
      integer, dimension(:), intent(in) :: maxNoConstl   !< Number of existing entries in for each ensemble in the constellation matrix
      integer, dimension(:,:), intent(in) :: unitCellType !< Array of constellation id and classification (core, boundary, or noise) per atom
      integer, dimension(:,:,:), intent(in) :: constellationsNeighType !< Every constellation atom´s neighbour type atoms.  This will tell which values in the constellation to look for during the applyhamiltionian step
      real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup  !< Couplings for saved constellations
      real(dblprec), dimension(:,:,:), intent(in) :: constellations !< Saved fixed unit cell configurations, these represent a configuration present in the domain with same configuration of unit cells in the neighborhood
      ! .. In/out variables
      real(dblprec), intent(inout) :: tenergy !< Total energy
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: emomM  !< Current magnetic moment vector
      ! .. Output variables
      real(dblprec), dimension(nim), intent(out) :: rx   !< Reaction coordinate
      real(dblprec), dimension(nim), intent(out) :: dene !< Derivative of the energy with respect to reaction coordinate
      !real(dblprec), dimension(nim), intent(out) :: im_ene  !< Energy of the images

      !.. Local variables
      integer :: ii,jj,itr,imax,kk,ie,ib
      real(dblprec) :: u0,fp,fp_i,ftol_mRy,fchk,u1,unim,ff,gf,fc,fcinv
      character(len=35) :: filn

      fc=mry/mub
      fcinv=mub/mry

      ! Call to allocate and initialize the arrays needed for the GNEB calculation
      call allocate_GNEB_arrays(1,nim,Natom)


      ftol_mRy = ftol*fc

      fchk=1+ftol_mRy
      itr=1

      call the_path(nim,Natom,emomM,pathl)

      !------------------------------------------------------------------------------
      ! Pre calculation of the field and energy
      !------------------------------------------------------------------------------
      print*,'! Calculation of the effective field'
      call timing(0,'Hamiltonian   ','ON')
      call effective_field(Natom,nim,1,Natom,emomM,mmom,external_field,&
         time_external_field,beff, beff1, beff2,OPT_flag, max_no_constellations,     &
         maxNoConstl, unitCellType, constlNCoup, constellations,                     &
         constellationsNeighType,tenergy,Num_macro,cell_index,            &
         emomM_macro,macro_nlistsize,NA,N1,N2,N3)
      call timing(0,'Hamiltonian   ','OF')


      call calc_energy(nHam,1,Natom,Nchmax, &
         conf_num,nim,Natom,Num_macro,1,plotenergy,              &
         0.0_dblprec,delta_t,do_lsf,lsf_field,     &
         lsf_interpolate,'N',simid,cell_index,macro_nlistsize,mmom,emom,emomM,       &
         emomM_macro,external_field,time_external_field,max_no_constellations,       &
         maxNoConstl,unitCellType,constlNCoup,constellations,OPT_flag,               &
         constellationsNeighType,tenergy,NA,N1,N2,N3)

       u(:) = ene%energy(:)*Natom

      call convert_force(Natom,nim,mmom,emom,beff)

      if (en_zero=='I') then
         u0 = u(1)
      elseif (en_zero=='F') then
         u0 = u(nim)
      else
         u0 = 0.0_dblprec
      end if
      u1 = u(1)
      unim = u(nim)

      call tang(nim,Natom,1,mmom,emomM,tau_i)
      fp_i = 0.0_dblprec
      !$omp parallel do default(shared), private(jj,kk), reduction(+:fp_i)
      !!, collapse(2)
      do jj=1,Natom
         do kk=1,3
            fp_i = fp_i + beff(kk,jj,1)*tau_i(kk,jj)
         end do
      end do
      !$omp end parallel do
      fpp(1) = fp_i

      if (fixed_if=='N') then
         ff = 0.0_dblprec
         gf = 0.0_dblprec
         !$omp parallel do default(shared), private(jj,kk), reduction(+:ff,gf)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               ff = ff + beff(kk,jj,1)*beff(kk,jj,1)
               gsp(kk,jj) = kappa*tau_i(kk,jj)*(pathl(2)-pathl(1))
               gf = gf + gsp(kk,jj)*beff(kk,jj,1)
            end do
         end do
         !$omp end parallel do
         ff = sqrt(ff)
         !$omp parallel do default(shared), private(kk,jj)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               beff(kk,jj,1) = gsp(kk,jj) - (gf/ff-u(1)+u1)*beff(kk,jj,1)/ff
            end do
         end do
         !$omp end parallel do
      end if

      do ii=2,nim-1

         call tang_spec(nim,Natom,ii,mmom,emomM,u,tau)
         call tang(nim,Natom,ii,mmom,emomM,tau_i)

         fp = 0.0_dblprec
         fp_i = 0.0_dblprec

         !$omp parallel do default(shared), private(jj,kk), reduction(+:fp,fp_i)
         !!,collapse(2)
         do jj=1,Natom
            do kk=1,3
               fp = fp + beff(kk,jj,ii)*tau(kk,jj)
               fp_i = fp_i + beff(kk,jj,ii)*tau_i(kk,jj)
            end do
         end do
         !$omp end parallel do
         fpp(ii) = fp_i

         !$omp parallel do default(shared), private(kk,jj)
         !!,collapse(2)
         do jj=1,Natom
            do kk=1,3
!             print*, ii
               beff(kk,jj,ii) = beff(kk,jj,ii) - tau(kk,jj)*fp +                    &
                  kappa*tau(kk,jj)*(pathl(ii+1)+pathl(ii-1)-2.0_dblprec*pathl(ii))
            end do
         end do
         !$omp end parallel do
      end do

      call tang(nim,Natom,nim,mmom,emomM,tau_i)
      fp_i = 0.0_dblprec
      !$omp parallel do default(shared), private(jj,kk),reduction(+:fp_i)
      !!, collapse(2)
      do jj=1,Natom
         do kk=1,3
            fp_i = fp_i + beff(kk,jj,nim)*tau_i(kk,jj)
         end do
      end do
      !$omp end parallel do
      fpp(nim) = fp_i

      if (fixed_if=='N') then
         ff = 0.0_dblprec
         gf = 0.0_dblprec
         !$omp parallel do default(shared), private(jj,kk),reduction(+:ff,gf)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               ff = ff + beff(kk,jj,nim)*beff(kk,jj,nim)
               gsp(kk,jj) = kappa*tau_i(kk,jj)*(pathl(nim)-pathl(nim-1))
               gf = gf + gsp(kk,jj)*beff(kk,jj,nim)
            end do
         end do
         !$omp end parallel do
         ff = sqrt(ff)
         !$omp parallel do default(shared), private(kk,jj)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               beff(kk,jj,nim) = gsp(kk,jj) - (gf/ff-u(nim)+unim)*beff(kk,jj,nim)/ff
            end do
         end do
         !$omp end parallel do
      end if

      fchk=0.0_dblprec
      imax = 1
      do ii=1,nim
         do jj=1,Natom
            do kk=1,3
               if (abs(beff(kk,jj,ii))>fchk) then
                  fchk = abs(beff(kk,jj,ii))
                  imax = ii
               end if
            end do
         end do
      end do
      itr=1
      filn = 'force_mep.'//trim(adjustl(simid))//'.out'
      open(ofileno,file = filn,access = 'sequential',action='write',status='replace')
      write(ofileno,'(a12,2x,a16,2x,a10)')"# Iter.","Force","Max Image"
      close(ofileno)

      open(ofileno,file = filn, access = 'sequential', action = 'write',            &
         status = 'old',position = 'append')
      write(ofileno,'(i12,2x,es16.8E3,2x,i10)',advance = 'no') itr,fchk*fcinv,imax
      close(ofileno)
      call save_path(Natom,nim,simid,3,emom,mmom,mode,do_mom_legacy)
      if (fixed_if == 'N') then
            ib = 1
            ie = nim-1
         else
            ib = 2
            ie = nim-1
      end if
      !------------------------------------------------------------------------------
      ! MAIN LOOP
      !------------------------------------------------------------------------------
      do while ((fchk>ftol_mRy).and.(itr<=itrmax))
         do ii=ib,ie
            call calc_axis_all(Natom,emomM(:,:,ii),beff(:,:,ii),ax(:,:,ii))
            call quick_min_coo(Natom,emomM(:,:,ii),vel(:,:,ii),beff(:,:,ii),coo,    &
               mass,delta_t)
            call calc_ang_all(Natom,emomM(:,:,ii),coo,ang(:,ii))
            do jj=1,Natom
               beff0(:,jj,ii) = beff(:,jj,ii)
               emomM(:,jj,ii) = coo(:,jj)
               emom(:,jj,ii)  = emomM(:,jj,ii)/mmom(jj,ii)
            end do
         end do

         call timing(0,'Hamiltonian   ','ON')
         ! Calculation of the effective field
         call effective_field(Natom,nim,1,Natom,emomM,mmom,external_field, &
            time_external_field,beff, beff1, beff2,OPT_flag, max_no_constellations,      &
            maxNoConstl, unitCellType, constlNCoup, constellations,                      &
            constellationsNeighType, tenergy,Num_macro,cell_index,             &
            emomM_macro,macro_nlistsize,NA,N1,N2,N3)
         call timing(0,'Hamiltonian   ','OF')

         ! Calculation of the total energy
        
         call calc_energy(nHam,itr,Natom,Nchmax,&
            conf_num,nim,Natom,Num_macro,1,plotenergy,               &
            0.0_dblprec,delta_t,do_lsf,lsf_field,      &
            lsf_interpolate,'N',simid,cell_index,macro_nlistsize,mmom,emom,emomM,        &
            emomM_macro,external_field,time_external_field,max_no_constellations,        &
            maxNoConstl,unitCellType,constlNCoup,constellations,OPT_flag,                &
            constellationsNeighType,tenergy,NA,N1,N2,N3)

         u(:) = ene%energy(:)*Natom

         call convert_force(Natom,nim,mmom,emom,beff)

         call the_path(nim,Natom,emomM,pathl)

         call tang(nim,Natom,1,mmom,emomM,tau_i)
         fp_i = 0.0_dblprec
         !$omp parallel do default(shared), private(jj,kk), reduction(+:fp_i)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               fp_i = fp_i + beff(kk,jj,1)*tau_i(kk,jj)
            end do
         end do
         !$omp end parallel do
         fpp(1) = fp_i

         call tang(nim,Natom,nim,mmom,emomM,tau_i)
         fp_i = 0.0_dblprec
         !$omp parallel do default(shared), private(jj,kk), reduction(+:fp_i)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               fp_i = fp_i + beff(kk,jj,nim)*tau_i(kk,jj)
            end do
         end do
         !$omp end parallel do
         fpp(nim) = fp_i

         if (fixed_if=='N') then
            call tang(nim,Natom,1,mmom,emomM,tau_i)
            ff = 0.0_dblprec
            gf = 0.0_dblprec
            !$omp parallel do default(shared), private(jj,kk), reduction(+:ff,gf)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  ff = ff + beff(kk,jj,1)*beff(kk,jj,1)
                  gsp(kk,jj) = kappa*tau_i(kk,jj)*(pathl(2)-pathl(1))
                  gf = gf + gsp(kk,jj)*beff(kk,jj,1)
               end do
            end do
            !$omp end parallel do
            ff = sqrt(ff)
            !$omp parallel do default(shared), private(jj,kk)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  beff(kk,jj,1) = gsp(kk,jj) - (gf/ff-u(1)+u1)*beff(kk,jj,1)/ff
               end do
            end do
            !$omp end parallel do
         end if

         do ii=2,nim-1
            call tang_spec(nim,Natom,ii,mmom,emomM,u,tau)
            call tang(nim,Natom,ii,mmom,emomM,tau_i)

            fp = 0.0_dblprec
            fp_i = 0.0_dblprec

            !$omp parallel do default(shared), private(jj,kk),reduction(+:fp,fp_i)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  fp = fp + beff(kk,jj,ii)*tau(kk,jj)
                  fp_i = fp_i + beff(kk,jj,ii)*tau_i(kk,jj)
               end do
            end do
            !$omp end parallel do
            fpp(ii) = fp_i

            !$omp parallel do default(shared), private(jj,kk)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  beff(kk,jj,ii) = beff(kk,jj,ii) - tau(kk,jj)*fp +                 &
                     kappa*tau(kk,jj)*(pathl(ii+1)+pathl(ii-1)-2.0_dblprec*pathl(ii))
               end do
            end do
            !$omp end parallel do
         end do

         if (fixed_if=='N') then
            call tang(nim,Natom,nim,mmom,emomM,tau_i)
            ff = 0.0_dblprec
            gf = 0.0_dblprec
            !$omp parallel do default(shared), private(jj,kk),reduction(+:ff,gf)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  ff = ff + beff(kk,jj,nim)*beff(kk,jj,nim)
                  gsp(kk,jj) = kappa*tau_i(kk,jj)*(pathl(nim)-pathl(nim-1))
                  gf = gf + gsp(kk,jj)*beff(kk,jj,nim)
               end do
            end do
            !$omp end parallel do
            ff = sqrt(ff)
            !$omp parallel do default(shared), private(jj,kk)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  beff(kk,jj,nim) = gsp(kk,jj) - (gf/ff-u(nim)+unim)*beff(kk,jj,nim)/ff
               end do
            end do
            !$omp end parallel do
         end if

         do ii=ib,ie
            call quick_min_vel(Natom,vel(:,:,ii),beff0(:,:,ii),beff(:,:,ii),        &
               ax(:,:,ii),ang(:,ii),mass,delta_t)
         end do

         call project_vel_gneb(Natom,nim,vel,beff)

         imax=1
         fchk=0.0_dblprec

         do ii=ib,ie
            do jj=1,Natom
               do kk=1,3
                  if (abs(beff(kk,jj,ii))>fchk) then
                     fchk = abs(beff(kk,jj,ii))
                     imax = ii
                  end if
               end do
            end do
         end do

         itr=itr+1

         if (mod(itr,every).eq.0) then
            call prn_gneb_progress(itr, itrmax, fchk*fcinv, imax,'N',0)

            open(ofileno,file = filn, access = 'sequential', action = 'write',      &
               status = 'old',position = 'append')
            write(ofileno,'(i12,2x,es16.8E3,2x,i3)',advance = 'no') itr,fchk*fcinv, &
               imax

            close(ofileno)

            call save_path(Natom,nim,simid,3,emom,mmom,mode,do_mom_legacy)
         end if
      end do

      if (itr>itrmax) then
         write(*,'(a)') 'WARNING: exceeded maximum number of iterations in GNEB'
      end if

      do ii=1,nim
         ene%energy(ii) = (u(ii)-u0)*Natom  !!!
         dene(ii) = -fpp(ii)*fcinv
         rx(ii) = pathl(ii)
      end do
      ! Call to deallocate the arrays needed for the GNEB calculation
      call allocate_GNEB_arrays(-1,nim,Natom)

   end subroutine gneb_mep

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: gneb_ci_mep
   !> @brief Evolves the path to the nearest MEP according to the CI-GNEB+VPO algorithm
   !> @details This tries to find the MEP by using a climbing image approach
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine gneb_ci_mep(nim,nHam,Natom,every,do_dm,do_pd,do_bq,do_ring,do_chir,do_dip,&
      itrmax,Nchmax,conf_num,do_biqdm,Num_macro,do_jtensor,plotenergy,do_anisotropy,&
      max_no_constellations,mass,ftol,kappa,delta_t,simid,do_lsf,en_zero,fixed_if,  &
      mult_axis,exc_inter,lsf_field,lsf_interpolate,OPT_flag, cell_index,           &
      macro_nlistsize,mmom,emom,emomM_macro,external_field,maxNoConstl,unitCellType,&
      constellationsNeighType,constlNCoup,constellations,tenergy,emomM,ci,rx,&
      dene,NA,N1,N2,N3,mode,do_mom_legacy)

      implicit none

      ! .. Input variables
      integer, intent(in) :: NA              !< Number of atoms in one cell
      integer, intent(in) :: N1              !< Number of cell repetitions in x direction
      integer, intent(in) :: N2              !< Number of cell repetitions in y direction
      integer, intent(in) :: N3              !< Number of cell repetitions in z direction
      integer, intent(in) :: nim             !< Number of images in GNEB calculations
      integer, intent(in) :: nHam            !< Number of atoms in Hamiltonian
      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: every           !< Save path every 'every' step
      integer, intent(in) :: do_dm           !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd           !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq           !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_ring         !< Add four-spin ring (4SR) term to Hamiltonian (0/1)      
      integer, intent(in) :: do_dip          !< Calculate dipole-dipole contribution (0/1)
      integer, intent(in) :: itrmax          !< Maximum number of iterations
      integer, intent(in) :: Nchmax          !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_chir         !< Add scalar chirality exchange (CHIR) term to Hamiltonian (0/1)
      integer, intent(in) :: conf_num        !< number of configurations for LSF
      integer, intent(in) :: do_biqdm        !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: Num_macro       !< Number of macrocells in the system
      integer, intent(in) :: do_jtensor      !< Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
      integer, intent(in) :: plotenergy      !< Calculate and plot energy (0/1)
      integer, intent(in) :: do_anisotropy   !< Read anisotropy data (1/0)
      integer, intent(in) :: max_no_constellations !< The maximum (global) length of the constellation matrix
      real(dblprec), intent(in) :: mass      !< mass of the point
      real(dblprec), intent(in) :: ftol      !< Tolerance
      real(dblprec), intent(in) :: kappa     !< spring constant
      real(dblprec), intent(in) :: delta_t   !< timestep
      character(len=1), intent(in) :: mode      !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid     !< Name of simulation
      character(len=1), intent(in) :: do_lsf    !< Including LSF energy
      character(len=1), intent(in) :: en_zero   !< Level of zero energy
      character(len=1), intent(in) :: fixed_if  !< Fix endpoints of the path (Y) or allow them to move along energy isocontours (N)
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy 
      character(len=1), intent(in) :: exc_inter !< Interpolation of Jij (Y/N)
      character(len=1), intent(in) :: lsf_field     !< LSF field contribution (Local/Total)
      character(len=1), intent(in) :: do_mom_legacy  !< Flag to print/read moments in legacy output
      character(len=1), intent(in) :: lsf_interpolate    !< Interpolate LSF or not
      logical, intent(in) :: OPT_flag !< Optimization flag
      integer, dimension(Natom), intent(in) :: cell_index    !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      real(dblprec), dimension(Natom,nim), intent(in) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Num_macro,nim), intent(in) :: emomM_macro    !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,Natom,nim), intent(in)     :: external_field !< External magnetic field
      integer, dimension(:), intent(in) :: maxNoConstl   !< Number of existing entries in for each ensemble in the constellation matrix
      integer, dimension(:,:), intent(in) :: unitCellType !< Array of constellation id and classification (core, boundary, or noise) per atom
      integer, dimension(:,:,:), intent(in) :: constellationsNeighType !< Every constellation atom´s neighbour type atoms.  This will tell which values in the constellation to look for during the applyhamiltionian step
      real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup  !< Couplings for saved constellations
      real(dblprec), dimension(:,:,:), intent(in) :: constellations  !< Saved fixed unit cell configurations, these represent a configuration present in the domain with same configuration of unit cells in the neighborhood
      ! .. In/out variables
      real(dblprec), intent(inout) :: tenergy !< Total energy
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: emom     !< Current unit moment vector
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: emomM  !< Current magnetic moment vector
      ! .. Output variables
      integer, intent(out) :: ci !< Index of the climbing image
      real(dblprec), dimension(nim), intent(out) :: rx   !< Reaction coordinate
      real(dblprec), dimension(nim), intent(out) :: dene !< Derivative of the energy with respect to reaction coordinate
      !real(dblprec), dimension(nim), intent(out) :: im_ene  !< Energy of the images

      ! .. Local variables
      integer :: ii,jj,itr,imax,kk,ib,ie
      real(dblprec) :: u0, fp,fp_i,ftol_mRy,fchk,u1,unim,ff,gf,fc,fcinv
      character(len=35) :: filn

      fc=mry/mub
      fcinv=mub/mry

      ! Call to allocate and initialize the arrays needed for the GNEB calculation
      call allocate_GNEB_arrays(1,nim,Natom)

      ftol_mRy = ftol*fc

      fchk=1+ftol_mRy
      itr=1

      call the_path(nim,Natom,emomM,pathl)

      !------------------------------------------------------------------------------
      ! Pre calculation of the field and energy
      !------------------------------------------------------------------------------
      ! Calculation of the effective field
      call timing(0,'Hamiltonian   ','ON')
      call effective_field(Natom,nim,1,Natom,emomM,mmom,external_field,&
         time_external_field,beff,beff1,beff2,OPT_flag,max_no_constellations,        &
         maxNoConstl,unitCellType,constlNCoup,constellations,                        &
         constellationsNeighType,tenergy,Num_macro,cell_index,emomM_macro, &
         macro_nlistsize,NA,N1,N2,N3)
      call timing(0,'Hamiltonian   ','OF')

      ! Calculate the total energy of the system
      call calc_energy(nHam,1,Natom,Nchmax, &
         conf_num,nim,Natom,Num_macro,1,plotenergy,              &
         0.0_dblprec,delta_t,do_lsf,lsf_field,     &
         lsf_interpolate,'N',simid,cell_index,macro_nlistsize,mmom,emom,emomM,       &
         emomM_macro,external_field,time_external_field,max_no_constellations,       &
         maxNoConstl,unitCellType,constlNCoup,constellations,OPT_flag,               &
         constellationsNeighType,tenergy,NA,N1,N2,N3)
!!!!!!!!!!!!!!!!!!!
      u(:) = ene%energy(:)*Natom

      call convert_force(Natom,nim,mmom,emom,beff)

      ci =1

      do ii=1,nim
         if (u(ii)>u(ci)) then
            ci = ii
         end if
      end do

      if (en_zero=='I') then
         u0 = u(1)
      elseif (en_zero=='F') then
         u0 = u(nim)
      else
         u0 = 0.0_dblprec
      end if

      u1 = u(1)
      unim = u(nim)

      call tang(nim,Natom,1,mmom,emomM,tau_i)
      fp_i = 0.0_dblprec
      !$omp parallel do default(shared), private(jj,kk),reduction(+:fp_i)
      !!, collapse(2)
      do jj=1,Natom
         do kk=1,3
            fp_i = fp_i + beff(kk,jj,1)*tau_i(kk,jj)
         end do
      end do
      !$omp end parallel do
      fpp(1) = fp_i

      if (fixed_if=='N') then
         ff = 0.0_dblprec
         gf = 0.0_dblprec
         !$omp parallel do default(shared), private(jj,kk),reduction(+:ff,gf)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               ff = ff + beff(kk,jj,1)*beff(kk,jj,1)
               gsp(kk,jj) = kappa*tau_i(kk,jj)*(pathl(2)-pathl(1))
               gf = gf + gsp(kk,jj)*beff(kk,jj,1)
            end do
         end do
         !$omp end parallel do
         ff = sqrt(ff)
         !$omp parallel do default(shared), private(jj,kk)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               beff(kk,jj,1) = gsp(kk,jj) - (gf/ff-u(1)+u1)*beff(kk,jj,1)/ff
            end do
         end do
         !$omp end parallel do
      end if

      do ii=2,nim-1
         call tang_spec(nim,Natom,ii,mmom,emomM,u,tau)
         call tang(nim,Natom,ii,mmom,emomM,tau_i)

         fp = 0.0_dblprec
         fp_i = 0.0_dblprec

         !$omp parallel do default(shared), private(jj,kk), reduction(+:fp,fp_i)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               fp = fp + beff(kk,jj,ii)*tau(kk,jj)
               fp_i = fp_i + beff(kk,jj,ii)*tau_i(kk,jj)
            end do
         end do
         !$omp end parallel do
         fpp(ii) = fp_i

         if (ii==ci) then
            !$omp parallel do default(shared), private(jj,kk)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  beff(kk,jj,ii) = beff(kk,jj,ii) - 2.0_dblprec*tau(kk,jj)*fp
               end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do default(shared), private(jj,kk)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  beff(kk,jj,ii) = beff(kk,jj,ii) - tau(kk,jj)*fp +                 &
                     kappa*tau(kk,jj)*(pathl(ii+1)+pathl(ii-1)-2.0_dblprec*pathl(ii))
               end do
            end do
            !$omp end parallel do
         end if
      end do

      call tang(nim,Natom,nim,mmom,emomM,tau_i)
      fp_i = 0.0_dblprec
      !$omp parallel do default(shared), private(jj,kk),reduction(+:fp_i)
      !!, collapse(2)
      do jj=1,Natom
         do kk=1,3
            fp_i = fp_i + beff(kk,jj,nim)*tau_i(kk,jj)
         end do
      end do
      !$omp end parallel do
      fpp(nim) = fp_i

      if (fixed_if=='N') then
         ff = 0.0_dblprec
         gf = 0.0_dblprec
         !$omp parallel do default(shared), private(jj,kk), reduction(+:ff,gf)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               ff = ff + beff(kk,jj,nim)*beff(kk,jj,nim)
               gsp(kk,jj) = kappa*tau_i(kk,jj)*(pathl(nim)-pathl(nim-1))
               gf = gf + gsp(kk,jj)*beff(kk,jj,nim)
            end do
         end do
         !$omp end parallel do
         ff = sqrt(ff)
         !$omp parallel do default(shared), private(jj,kk)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               beff(kk,jj,nim) = gsp(kk,jj) - (gf/ff-u(nim)+unim)*beff(kk,jj,nim)/ff
            end do
         end do
         !$omp end parallel do
      end if

      fchk=0.0_dblprec
      imax = 1
      do ii=1,nim
         do jj=1,Natom
            do kk=1,3
               if (abs(beff(kk,jj,ii))>fchk) then
                  fchk = abs(beff(kk,jj,ii))
                  imax = ii
               end if
            end do
         end do
      end do
      itr=1

      filn = 'force_mep.'//trim(adjustl(simid))//'.out'
      open(ofileno,file = filn,access = 'sequential',action='write',status='replace')
      write(ofileno,'(a12,2x,a16,2x,a10,2x,a10)') "# Iter.","Force","Max Image","CI"
      close(ofileno)

      open(ofileno,file = filn, access = 'sequential', action = 'write',            &
         status = 'old',position = 'append')
      write(ofileno,'(i12,2x,es16.8E3,2x,i10,2x,i10)',advance = 'no') itr,          &
         fchk*fcinv,imax,ci

      close(ofileno)

      call save_path(Natom,nim,simid,3,emom,mmom,mode,do_mom_legacy)

      if (fixed_if == 'N') then
         ib = 1
         ie = nim-1
      else
         ib = 2
         ie = nim-1
      end if

      !------------------------------------------------------------------------------
      ! MAIN LOOP
      !------------------------------------------------------------------------------
      do while ((fchk>ftol_mRy).and.(itr<=itrmax))
         ! Set the cimbling image to be the first one
         ci =1
         do ii=ib,ie
            call calc_axis_all(Natom,emomM(:,:,ii),beff(:,:,ii),ax(:,:,ii))
            call quick_min_coo(Natom,emomM(:,:,ii),vel(:,:,ii),beff(:,:,ii),coo,    &
               mass,delta_t)
            call calc_ang_all(Natom,emomM(:,:,ii),coo,ang(:,ii))
            do jj=1,Natom
               beff0(:,jj,ii) = beff(:,jj,ii)
               emomM(:,jj,ii) = coo(:,jj)
               emom(:,jj,ii)  = emomM(:,jj,ii)/mmom(jj,ii)
            end do
         end do

         ! Calculation of the effective field
         call timing(0,'Hamiltonian   ','ON')
         call effective_field(Natom,nim,1,Natom,emomM,mmom,external_field,&
            time_external_field,beff, beff1, beff2, OPT_flag, max_no_constellations,    &
            maxNoConstl,unitCellType, constlNCoup, constellations,                      &
            constellationsNeighType, tenergy,Num_macro,cell_index,            &
            emomM_macro,macro_nlistsize,NA,N1,N2,N3)
         call timing(0,'Hamiltonian   ','OF')

         ! Calculate the total energy per spin of the system
         call calc_energy(nHam,itr,Natom,Nchmax,  &
            conf_num,nim,Natom,Num_macro,1,plotenergy,   &
            0.0_dblprec,delta_t,do_lsf,lsf_field, &
            lsf_interpolate,'N',simid,cell_index,macro_nlistsize,mmom,emom,emomM,   &
            emomM_macro,external_field,time_external_field,max_no_constellations,   &
            maxNoConstl,unitCellType,constlNCoup,constellations,OPT_flag,           &
            constellationsNeighType,tenergy,NA,N1,N2,N3)

         ! Calculate the total energy of the sample (Not per spin)
         u(:) = ene%energy(:)*Natom

         call convert_force(Natom,nim,mmom,emom,beff)
         ! If the energy of a given image is lower than the one of the climbing image
         ! the new climbing image is that one with lower energy
         do ii=2,nim-1
            if (u(ii)>u(ci)) then
               ci = ii
            end if
         end do

         call the_path(nim,Natom,emomM,pathl)

         call tang(nim,Natom,1,mmom,emomM,tau_i)
         fp_i = 0.0_dblprec
         !$omp parallel do default(shared), private(jj,kk),reduction(+:fp_i)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               fp_i = fp_i + beff(kk,jj,1)*tau_i(kk,jj)
            end do
         end do
         !$omp end parallel do
         fpp(1) = fp_i

         call tang(nim,Natom,nim,mmom,emomM,tau_i)
         fp_i = 0.0_dblprec
         !$omp parallel do default(shared), private(jj,kk),reduction(+:fp_i)
         !!, collapse(2)
         do jj=1,Natom
            do kk=1,3
               fp_i = fp_i + beff(kk,jj,nim)*tau_i(kk,jj)
            end do
         end do
         !$omp end parallel do
         fpp(nim) = fp_i

         if (fixed_if=='N') then
            call tang(nim,Natom,1,mmom,emomM,tau_i)
            ff = 0.0_dblprec
            gf = 0.0_dblprec
            !$omp parallel do default(shared), private(jj,kk), reduction(+:ff,gf)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  ff = ff + beff(kk,jj,1)*beff(kk,jj,1)
                  gsp(kk,jj) = kappa*tau_i(kk,jj)*(pathl(2)-pathl(1))
                  gf = gf + gsp(kk,jj)*beff(kk,jj,1)
               end do
            end do
            !$omp end parallel do
            ff = sqrt(ff)
            !$omp parallel do default(shared), private(jj,kk)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  beff(kk,jj,1) = gsp(kk,jj) - (gf/ff-u(1)+u1)*beff(kk,jj,1)/ff
               end do
            end do
            !$omp end parallel do
         end if

         do ii=2,nim-1
            call tang_spec(nim,Natom,ii,mmom,emomM,u,tau)
            call tang(nim,Natom,ii,mmom,emomM,tau_i)
            fp = 0.0_dblprec
            fp_i = 0.0_dblprec
            !$omp parallel do default(shared), private(jj,kk), reduction(+:fp,fp_i)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  fp = fp + beff(kk,jj,ii)*tau(kk,jj)
                  fp_i = fp_i + beff(kk,jj,ii)*tau_i(kk,jj)
               end do
            end do
            !$omp end parallel do
            fpp(ii) = fp_i

            if (ii==ci) then
               !$omp parallel do default(shared), private(jj,kk)
               !!, collapse(2)
               do jj=1,Natom
                  do kk=1,3
                     beff(kk,jj,ii) = beff(kk,jj,ii) - 2.0_dblprec*tau(kk,jj)*fp
                  end do
               end do
               !$omp end parallel do
            else
               !$omp parallel do default(shared), private(jj,kk)
               !!, collapse(2)
               do jj=1,Natom
                  do kk=1,3
                     beff(kk,jj,ii) = beff(kk,jj,ii) - tau(kk,jj)*fp +              &
                     kappa*tau(kk,jj)*(pathl(ii+1)+pathl(ii-1)-2.0_dblprec*pathl(ii))
                  end do
               end do
               !$omp end parallel do
            end if
         end do

         if (fixed_if=='N') then
            call tang(nim,Natom,nim,mmom,emomM,tau_i)
            ff = 0.0_dblprec
            gf = 0.0_dblprec
            !$omp parallel do default(shared), private(jj,kk), reduction(+:ff,gf)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  ff = ff + beff(kk,jj,nim)*beff(kk,jj,nim)
                  gsp(kk,jj) = kappa*tau_i(kk,jj)*(pathl(nim)-pathl(nim-1))
                  gf = gf + gsp(kk,jj)*beff(kk,jj,nim)
               end do
            end do
            !$omp end parallel do
            ff = sqrt(ff)
            !$omp parallel do default(shared), private(jj,kk)
            !!, collapse(2)
            do jj=1,Natom
               do kk=1,3
                  beff(kk,jj,nim) = gsp(kk,jj) - (gf/ff-u(nim)+unim)*beff(kk,jj,nim)/ff
               end do
            end do
            !$omp end parallel do
         end if

         do ii=ib,ie
            call quick_min_vel(Natom,vel(:,:,ii),beff0(:,:,ii),beff(:,:,ii),        &
               ax(:,:,ii),ang(:,ii),mass,delta_t)
         end do

         call project_vel_gneb(Natom,nim,vel,beff)

         imax=1
         fchk=0.0_dblprec

         do ii=ib,ie
            do jj=1,Natom
               do kk=1,3
                  if (abs(beff(kk,jj,ii))>fchk) then
                     fchk = abs(beff(kk,jj,ii))
                     imax = ii
                  end if
               end do
            end do
         end do

         itr=itr+1

         if (mod(itr,every).eq.0) then
            call prn_gneb_progress(itr, itrmax, fchk*fcinv, imax,'Y',ci)
            call save_en(nim,pathl,ene%energy(:),-fpp*fcinv,pathl(nim),             &
               'en_path.restart.out','Y')
            open(ofileno,file = filn, access = 'sequential', action = 'write',      &
               status = 'old',position = 'append')
            write(ofileno,'(i12,2x,es16.8E3,2x,i10,2x,i10)',advance = 'no') itr,    &
               fchk*fcinv,imax,ci
            close(ofileno)
            call save_path(Natom,nim,simid,3,emom,mmom,mode,do_mom_legacy)
         end if
      end do

      if (itr>itrmax) then
         write(*,'(a)') 'WARNING: exceeded maximum number of iterations in CI-GNEB'
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ii=1,nim
         ene%energy(ii) = (u(ii)-u0)
         dene(ii) = -fpp(ii)*fcinv
         rx(ii) = pathl(ii)
      end do

      ! Call to deallocate the arrays needed for the GNEB calculation
      call allocate_GNEB_arrays(-1,nim,Natom)

   end subroutine gneb_ci_mep

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: the_path
   !> @brief Calculate poly-geodesic length of the path
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine the_path(nim,PP,path,pathlen)
      
      implicit none
      
      integer, intent(in) :: PP
      integer, intent(in) :: nim !! Number of images
      real(dblprec), dimension(3,PP,nim) ,intent(in) :: path
      ! .. Output variables
      real(dblprec), dimension(nim) ,intent(out) :: pathlen
      ! .. Local variables
      real(dblprec) :: tmp,l1
      integer :: ii,jj

      pathlen(1) = 0.0_dblprec
      do ii=2,nim
         tmp = 0.0_dblprec
         do jj=1,PP
            l1 = calc_ang(path(:,jj,ii-1),path(:,jj,ii))
            tmp = tmp+l1*l1
         end do
         pathlen(ii) = pathlen(ii-1) + sqrt(tmp)
      end do

   end subroutine the_path

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: tang_spec
   !> @brief Estimate tangent to the path according to the special definition: http://arxiv.org/abs/1502.05065
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine tang_spec(nim,PP,im,mmom,emomM,u,tau)
      implicit none

      integer, intent(in) :: im
      integer, intent(in) :: PP
      integer, intent(in) :: nim !< Number of images
      real(dblprec), dimension(nim), intent(in) :: u
      real(dblprec), dimension(PP,nim), intent(in) :: mmom !< Magnitude of the magnetic moments
      real(dblprec), dimension(3,PP,nim), intent(in) :: emomM !< Current magnetic moment vector
      ! .. Output variables
      real(dblprec), dimension(3,PP), intent(out) :: tau

      ! .. Local variables
      integer :: ii,kk
      real(dblprec) :: u1, u2, u3,dumin,dumax,tmp,mmomx,mtau,mtaum,mtaup
      real(dblprec), dimension(3,PP) :: taup,taum

      u1=u(im-1)
      u2=u(im)
      u3=u(im+1)
      if (u3>u2.and.u2>u1) then
         do ii=1,PP
            tau(:,ii) = emomM(:,ii,im+1)-emomM(:,ii,im)
         end do
      elseif (u1>u2.and.u2>u3) then
         do ii=1,PP
            tau(:,ii) = emomM(:,ii,im)-emomM(:,ii,im-1)
         end do
      else
         mtaup = 0.0_dblprec
         mtaum = 0.0_dblprec
         do ii=1,PP
            taup(:,ii) = emomM(:,ii,im+1)-emomM(:,ii,im)
            taum(:,ii) = emomM(:,ii,im)-emomM(:,ii,im-1)
            do kk=1,3
               mtaum = mtaum + taum(kk,ii)*taum(kk,ii)
               mtaup = mtaup + taup(kk,ii)*taup(kk,ii)
            end do
         end do

         mtaum = sqrt(mtaum)
         mtaup = sqrt(mtaup)
         do ii=1,PP
            do kk=1,3
               taum(kk,ii) = taum(kk,ii)/mtaum
               taup(kk,ii) = taup(kk,ii)/mtaup
            end do
         end do

         dumax = abs(u3-u2)
         dumin = abs(u1-u2)
         if (dumax<dumin) then
            tmp=dumax
            dumax=dumin
            dumin=tmp
         end if
         if (u3>u1) then
            do ii=1,PP
               tau(:,ii)=dumax*taup(:,ii)+dumin*taum(:,ii)
            end do
         else
            do ii=1,PP
               tau(:,ii)=dumin*taup(:,ii)+dumax*taum(:,ii)
            end do
         end if
      end if

      mtau = 0.0_dblprec
      do ii=1,PP
         mmomx = mmom(ii,im)**2
         tmp = 0.0_dblprec
         do kk=1,3
            tmp = tmp+tau(kk,ii)*emomM(kk,ii,im)/mmomx
         end do
         tau(:,ii) = tau(:,ii) - tmp*emomM(:,ii,im)
         do kk=1,3
            mtau = mtau + tau(kk,ii)*tau(kk,ii)
         end do
      end do
      mtau = sqrt(mtau)
      do ii=1,PP
         do kk=1,3
            tau(kk,ii) = tau(kk,ii)/mtau
         end do
      end do

   end subroutine tang_spec

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: tang_spec
   !> @brief Estimate tangent to the path
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine tang(nim,PP,im,mmom,emomM,tau)

      implicit none

      integer, intent(in) :: im
      integer, intent(in) :: PP
      integer, intent(in) :: nim !< Number of images
      real(dblprec), dimension(PP,nim), intent(in) :: mmom !< Magnitude of the magnetic moments
      real(dblprec), dimension(3,PP,nim), intent(in) :: emomM !< Current magnetic moment vector
      ! .. Output variables
      real(dblprec), dimension(3,PP), intent(out) :: tau
      !.. Local variables
      integer :: ii,kk
      real(dblprec) :: tmp,mmomx,mtau

      if (im==1) then
         do ii=1,PP
            tau(:,ii) = emomM(:,ii,im+1)-emomM(:,ii,im)
         end do
      elseif (im==nim) then
         do ii=1,PP
            tau(:,ii) = emomM(:,ii,im)-emomM(:,ii,im-1)
         end do
      else
         do ii=1,PP
            tau(:,ii) = emomM(:,ii,im+1)-emomM(:,ii,im-1)
         end do
      end if

      mtau = 0.0_dblprec
      do ii=1,PP
         mmomx = mmom(ii,im)**2
         tmp = 0.0_dblprec
         do kk=1,3
            tmp = tmp+tau(kk,ii)*emomM(kk,ii,im)/mmomx
         end do
         tau(:,ii) = tau(:,ii) - tmp*emomM(:,ii,im)
         do kk=1,3
            mtau = mtau + tau(kk,ii)*tau(kk,ii)
         end do
      end do
      mtau = sqrt(mtau)
      do ii=1,PP
         do kk=1,3
            tau(kk,ii) = tau(kk,ii)/mtau
         end do
      end do

   end subroutine tang

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: find_SP
   !> @brief Find the saddle point
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine find_SP(nim,rx,coef,rx0)

      use ErrorHandling,        only : ErrorHandling_ERROR

      implicit none

      integer, intent(in) :: nim !< Number of images
      real(dblprec), dimension(nim), intent(in) :: rx       !< Reaction coordinate
      real(dblprec), dimension(4,nim), intent(in) :: coef   !< Coefficients of the piecewise Hermite polynomials
      ! .. Output variables
      real(dblprec), intent(out) :: rx0
      ! .. Local variables
      integer :: ii, i1,ci
      real(dblprec) :: eps=epsilon(rx0),d,dl,x1,x2,q

      ci = 1
      do ii=1,nim
         if (coef(1,ii)>coef(1,ci)) then
            ci = ii
         end if
      end do

      i1=1

      if ((ci==1).or.(ci==nim)) then
         rx0 = rx(ci)
      else
         if (abs(coef(2,ci))<eps) then
            rx0 = rx(ci)
         else
            if (coef(2,ci)>0.0_dblprec) then
               i1 = ci
            elseif (coef(2,ci)<0.0_dblprec) then
               i1 = ci-1
            end if
            dl = rx(i1+1)-rx(i1)
            d = coef(3,i1)*coef(3,i1)-3.0_dblprec*coef(4,i1)*coef(2,i1)

            if (d<0.0_dblprec) then
               call ErrorHandling_ERROR("Energy maximum has not been found!")
            else
               q = -(coef(3,i1)+sign(1.0_dblprec,coef(3,i1))*sqrt(d))
               x1 = q/coef(4,i1)
               x2 = coef(2,i1)/q
               if ((x1>0.0_dblprec).and.(x1<dl)) then
                  rx0 = rx(i1)+x1
               elseif ((x2>0.0_dblprec).and.(x2<dl)) then
                  rx0 = rx(i1)+x2
               else
                  call ErrorHandling_ERROR("Energy maximum has not been found!")
               end if
            end if
         end if
      end if

   end subroutine find_SP

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: find_SP_conf_one
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine find_SP_conf_one(ni,nf,l1,l2,lsp,nsp)

      use RandomNumbers, only : rng_uniform
      use constants
      use math_functions, only: f_normalize_vec

      implicit none

      real(dblprec), intent(in) :: l1,l2,lsp
      real(dblprec), dimension(3), intent(in) :: ni,nf
      ! .. Output variables
      real(dblprec), dimension(3), intent(out) :: nsp
      ! .. Local variables
      integer :: ii,jj
      real(dblprec) :: dl, theta, angle, tmp, pr, eps=epsilon(angle)
      real(dblprec), dimension(3) :: ax,vec
      real(dblprec), dimension(2) :: rn

      dl = lsp-l1

      angle = calc_ang(ni,nf)
      if (angle<eps) then
         tmp=0.0_dblprec
         do ii = 1,3
            vec(ii) = (nf(ii) - ni(ii))*dl/(l2-l1)
            nsp(ii) = ni(ii) + vec(ii)
            tmp = tmp + nsp(ii)*nsp(ii)
         end do
         tmp = sqrt(tmp)
         nsp(:) = nsp(:)/tmp

      elseif (abs(angle-pi)<eps) then
         tmp = 0.0_dblprec
         do while (tmp<eps)
            pr = 0.0_dblprec
            do jj=1,3
               call rng_uniform(rn,2)
               ax(jj) = sign(1.0_dblprec,2.0_dblprec*rn(1)-1.0_dblprec)*(rn(2)+1.0_dblprec)
               pr = pr + ax(jj)*ni(jj)
            end do
            ax(:) = ax(:) - pr*ni(:)
            tmp = norm2(ax)
         end do
         tmp = norm2(ax)
         ax(:) = ax(:)/tmp
         theta = pi*dl/(l2-l1)

         nsp(1) = ni(1)*cos(theta) + sin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
         nsp(2) = ni(2)*cos(theta) - sin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
         nsp(3) = ni(3)*cos(theta) + sin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
         nsp = f_normalize_vec(nsp,3)

      else
         ax = calc_axis(ni,nf)
         theta = angle*dl/(l2-l1)
         nsp(1) = ni(1)*cos(theta) + sin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
         nsp(2) = ni(2)*cos(theta) - sin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
         nsp(3) = ni(3)*cos(theta) + sin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
         nsp = f_normalize_vec(nsp,3)

      end if

   end subroutine find_SP_conf_one

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: find_SP_conf
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine find_SP_conf(Natom,ni,nf,l1,l2,lsp,nsp)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), intent(in) :: l1,l2,lsp
      real(dblprec), dimension(3,Natom), intent(in) :: ni,nf
      ! .. Output variables
      real(dblprec), dimension(3,Natom), intent(out) :: nsp
      ! .. Local variables
      integer :: ii

      do ii=1,Natom
         call find_SP_conf_one(ni(:,ii),nf(:,ii),l1,l2,lsp,nsp(:,ii))
      end do

   end subroutine find_SP_conf

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: save_en
   !> @brief Print the path to file. Ensembles correspond to images in GNEB method
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine save_en(nim,x,y,dy,x0,filn,do_norm_rx)

      use ErrorHandling,        only : ErrorHandling_ERROR
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: nim       !< Number of samples
      real(dblprec), intent(in) :: x0  !< Normalization for the reaction coordinate
      character(*), intent(in) :: filn                !< filename
      character(len=1), intent(in) :: do_norm_rx      !< normalize reaction coordinate (Y/N)
      real(dblprec), dimension(nim), intent(in) :: x  !< Reaction coordinate
      real(dblprec), dimension(nim), intent(in) :: y  !< Energy
      real(dblprec), dimension(nim), intent(in) :: dy !< Derivative of the energy wrt x

      integer :: ii
      real(dblprec) :: norm

      norm=1.0_dblprec
      if (do_norm_rx=='Y') then
         norm = x0
      elseif (do_norm_rx=='N') then
         norm = 1.0_dblprec
      else
         call ErrorHandling_ERROR("Invalid value for do_norm_rx!")
      end if

      open(ofileno, file=filn, access = 'sequential',action = 'write', status = 'replace')
      write(ofileno,'(a6,a16,2x,a16,2x,a16)')"#Itr.","React_coord","Ene", "dEne/dReact_coord"
      do ii=1,nim
         write (ofileno,10002) ii,x(ii)/norm, y(ii), dy(ii)
      end do
      close(ofileno)
      return
      call ErrorHandling_ERROR("Error writing energy along the path to file!")
      10002 format (i6,es16.8E3,2x,es16.8E3,2x,es16.8E3)

   end subroutine save_en

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_gneb_progress
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine prn_gneb_progress(itr, itrmax, fchk, imax,do_ci,ci)

      implicit none

      integer, intent(in) :: ci !< Climbing image
      integer, intent(in) :: itr !< Current iteration
      integer, intent(in) :: imax
      integer, intent(in) :: itrmax !< Maxinum number of iterations
      real(dblprec), intent(in) :: fchk   !< Current force
      character(len=1), intent(in) :: do_ci
      character(35) :: num,num1,num2,num3
      write(num,'(i8)') idnint(real(itr,dblprec)*100/real(itrmax,dblprec))
      write(num1,'(es16.8E3)') fchk
      write(num2,'(i8)') imax
      write(num3,'(i8)') ci

      if (do_ci.eq.'Y') then
         write (*,'(2x,8a)') 'MP  ',trim(adjustl(num)),'% of itrmax.   fchk: ',trim(adjustl(num1)),'   imax: ',trim(adjustl(num2)),'   ci: ',trim(adjustl(num3))
      else
         write (*,'(2x,6a)') 'MP  ',trim(adjustl(num)),'% of itrmax.   fchk: ',trim(adjustl(num1)),'   imax: ',trim(adjustl(num2))
      end if

   end subroutine prn_gneb_progress

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: allocate_GNEB_arrays
   !> @brief Allocation/deallocation of the arrays needed for the GNEB calculation
   !> @author Jonathan Chico
   !> @date 27.03.2018
   !---------------------------------------------------------------------------------
   subroutine allocate_GNEB_arrays(flag,nim,Natom)

      implicit none

      integer, intent(in) :: flag   !< Allocate deallocate flag
      integer, intent(in) :: nim    !< Number og GNEB images, corresponding to Mensembles
      integer, intent(in) :: Natom  !< Number of atoms in the system

      integer :: i_stat,i_all

      ! Allocate the GNEB arrays
      if (flag>0) then
         allocate(fpp(nim),stat=i_stat)
         call memocc(i_stat,product(shape(fpp))*kind(fpp),'fpp','allocate_GNEB_arrays')
         fpp=0.0_dblprec
         allocate(pathl(nim),stat=i_stat)
         call memocc(i_stat,product(shape(pathl))*kind(pathl),'pathl','allocate_GNEB_arrays')
         pathl=0.0_dblprec
         allocate(vel(3,Natom,nim),stat=i_stat)
         call memocc(i_stat,product(shape(vel))*kind(vel),'vel','allocate_GNEB_arrays')
         vel=0.0_dblprec
         allocate(coo(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(coo))*kind(coo),'coo','allocate_GNEB_arrays')
         coo=0.0_dblprec
         allocate(ax(3,Natom,nim),stat=i_stat)
         call memocc(i_stat,product(shape(ax))*kind(ax),'ax','allocate_GNEB_arrays')
         ax=0.0_dblprec
         allocate(ang(Natom,nim),stat=i_stat)
         call memocc(i_stat,product(shape(ang))*kind(ang),'ang','allocate_GNEB_arrays')
         ang=0.0_dblprec
         allocate(beff(3,Natom,nim),stat=i_stat)
         call memocc(i_stat,product(shape(beff))*kind(beff),'beff','allocate_GNEB_arrays')
         beff=0.0_dblprec
         allocate(beff1(3,Natom,nim),stat=i_stat)
         call memocc(i_stat,product(shape(beff1))*kind(beff1),'beff1','allocate_GNEB_arrays')
         beff1=0.0_dblprec
         allocate(beff2(3,Natom,nim),stat=i_stat)
         call memocc(i_stat,product(shape(beff2))*kind(beff2),'beff2','allocate_GNEB_arrays')
         beff2=0.0_dblprec
         allocate(beff0(3,Natom,nim),stat=i_stat)
         call memocc(i_stat,product(shape(beff0))*kind(beff0),'beff0','allocate_GNEB_arrays')
         beff0=0.0_dblprec
         allocate(u(nim),stat=i_stat)
         call memocc(i_stat,product(shape(u))*kind(u),'u','allocate_GNEB_arrays')
         u=0.0_dblprec
         allocate(time_external_field(3,Natom,nim),stat=i_stat)
         call memocc(i_stat,product(shape(time_external_field))*kind(time_external_field),'time_external_field','allocate_GNEB_arrays')
         time_external_field=0.0_dblprec
         allocate(tau(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(tau))*kind(tau),'tau','allocate_GNEB_arrays')
         tau=0.0_dblprec
         allocate(tau_i(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(tau_i))*kind(tau_i),'tau_i','allocate_GNEB_arrays')
         tau_i=0.0_dblprec
         allocate(gsp(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(gsp))*kind(gsp),'gsp','allocate_GNEB_arrays')
         gsp=0.0_dblprec
      else
         if (allocated(u)) then 
            i_all=-product(shape(u))*kind(u)
            deallocate(u,stat=i_stat)
            call memocc(i_stat,i_all,'u','allocate_GNEB_arrays')
         endif

         if (allocated(ax)) then 
            i_all=-product(shape(ax))*kind(ax)
            deallocate(ax,stat=i_stat)
            call memocc(i_stat,i_all,'ax','allocate_GNEB_arrays')
         endif
         
         if (allocated(gsp)) then 
            i_all=-product(shape(gsp))*kind(gsp)
            deallocate(gsp,stat=i_stat)
            call memocc(i_stat,i_all,'gsp','allocate_GNEB_arrays')
         endif

         if (allocated(vel)) then 
            i_all=-product(shape(vel))*kind(vel)
            deallocate(vel,stat=i_stat)
            call memocc(i_stat,i_all,'vel','allocate_GNEB_arrays')
         endif

         if (allocated(coo)) then
            i_all=-product(shape(coo))*kind(coo)
            deallocate(coo,stat=i_stat)
            call memocc(i_stat,i_all,'coo','allocate_GNEB_arrays')
         endif

         if (allocated(ang)) then
            i_all=-product(shape(ang))*kind(ang)
            deallocate(ang,stat=i_stat)
            call memocc(i_stat,i_all,'ang','allocate_GNEB_arrays')
         endif
            
         if (allocated(fpp)) then
            i_all=-product(shape(fpp))*kind(fpp)
            deallocate(fpp,stat=i_stat)
            call memocc(i_stat,i_all,'fpp','allocate_GNEB_arrays')
         endif

         if (allocated(tau)) then
            i_all=-product(shape(tau))*kind(tau)
            deallocate(tau,stat=i_stat)
            call memocc(i_stat,i_all,'tau','allocate_GNEB_arrays')
         endif

         if (allocated(beff)) then
            i_all=-product(shape(beff))*kind(beff)
            deallocate(beff,stat=i_stat)
            call memocc(i_stat,i_all,'beff','allocate_GNEB_arrays')
         endif

         if (allocated(beff0)) then
            i_all=-product(shape(beff0))*kind(beff0)
            deallocate(beff0,stat=i_stat)
            call memocc(i_stat,i_all,'beff0','allocate_GNEB_arrays')
         endif

         if (allocated(beff1)) then
            i_all=-product(shape(beff1))*kind(beff1)
            deallocate(beff1,stat=i_stat)
            call memocc(i_stat,i_all,'beff1','allocate_GNEB_arrays')
         endif

         if (allocated(beff2)) then
            i_all=-product(shape(beff2))*kind(beff2)
            deallocate(beff2,stat=i_stat)
            call memocc(i_stat,i_all,'beff2','allocate_GNEB_arrays')
         endif

         if (allocated(pathl)) then
            i_all=-product(shape(pathl))*kind(pathl)
            deallocate(pathl,stat=i_stat)
            call memocc(i_stat,i_all,'pathl','allocate_GNEB_arrays')
         endif

         if (allocated(tau_i)) then
            i_all=-product(shape(tau_i))*kind(tau_i)
            deallocate(tau_i,stat=i_stat)
            call memocc(i_stat,i_all,'tau_i','allocate_GNEB_arrays')
         endif

         if (allocated(time_external_field)) then
            i_all=-product(shape(time_external_field))*kind(time_external_field)
            deallocate(time_external_field,stat=i_stat)
            call memocc(i_stat,i_all,'time_external_field','allocate_GNEB_arrays')
         endif

      endif

   end subroutine allocate_GNEB_arrays

end module GNEB
