!-------------------------------------------------------------------------------
! MODULE: energyminima
!> @brief
!> Routines for minimizing the effective fields
!> @author
!> Pavel Bessarab
!-------------------------------------------------------------------------------
module energyminima
   use Parameters
   use HamiltonianActions
   use VPO
   use OSO
   use Constants
   use Energy
   use Pathinit, only: save_path

   implicit none

   ! .. Global allocatable variables
   real(dblprec), dimension(:), allocatable :: fchk
   real(dblprec), dimension(:,:), allocatable :: coo
   real(dblprec), dimension(:,:), allocatable :: ang
   real(dblprec), dimension(:,:,:), allocatable :: ax
   real(dblprec), dimension(:,:,:), allocatable :: vel
   real(dblprec), dimension(:,:,:), allocatable :: beff !< Total effective field from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: beff1 !< Internal effective field from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: beff2 !< External field from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: beff0
   real(dblprec), dimension(:,:,:), allocatable :: time_external_field !< External time-dependent magnetic field, zero here

   private

   public :: vpo_min, save_ifmom, vpo_min_single
contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: vpo_min
   !> @brief
   !> Evolves the system to the nearest energy minimum according to the VPO algorithm
   !> @author
   !> Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine vpo_min(nHam,every,Natom,Nchmax, &
      itrmax,OPT_flag,conf_num,Mensemble,Num_macro,plotenergy,  &
      max_no_constellations,ftol,mass,delta_t,simid,do_lsf,lsf_field, &
      lsf_interpolate,cell_index,macro_nlistsize,mmom,emom,     &
      emomM_macro,external_field,maxNoConstl,unitCellType,constlNCoup,              &
      constellations,constellationsNeighType,energy,emomM,NA,N1,N2,N3,mode,         &
      do_mom_legacy)

      implicit none

      integer, intent(in) :: NA           !< Number of atoms in one cell
      integer, intent(in) :: N1           !< Number of cell repetitions in x direction
      integer, intent(in) :: N2           !< Number of cell repetitions in y direction
      integer, intent(in) :: N3           !< Number of cell repetitions in z direction
      integer, intent(in) :: nHam         !< Number of atoms in Hamiltonian
      integer, intent(in) :: every        !< Save path every 'every' step
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Nchmax       !< Max number of chemical components on each site in cell
      integer, intent(in) :: itrmax       !< Maximum number of iterations
      logical, intent(in) :: OPT_flag
      integer, intent(in) :: conf_num     !< number of configurations for LSF
      integer, intent(in) :: Mensemble    !< Number of images in GNEB calculations
      integer, intent(in) :: Num_macro    !< Number of macrocells in the system
      integer, intent(in) :: plotenergy   !< Calculate and plot energy (0/1)
      integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
      real(dblprec), intent(in) :: ftol               !< Tolerance
      real(dblprec), intent(in) :: mass               !< mass of the point
      real(dblprec), intent(in) :: delta_t            !< timestep
      character(len=1), intent(in) :: mode            !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid           !< Name of simulation
      character(len=1), intent(in) :: do_lsf          !< Including LSF energy
      character(len=1), intent(in) :: lsf_field       !< LSF field contribution (Local/Total)
      character(len=1), intent(in) :: do_mom_legacy    !< Flag to print/read moments in legacy output
      character(len=1), intent(in) :: lsf_interpolate !< Interpolate LSF or not
      integer, dimension(Natom), intent(in) :: cell_index    !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field
      integer, dimension(:), intent(in) :: maxNoConstl
      integer, dimension(:,:), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup
      real(dblprec), dimension(:,:,:), intent(in) :: constellations
      integer, dimension(:,:,:), intent(in) :: constellationsNeighType

      ! .. In/out variables
      real(dblprec), intent(inout) :: energy
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom     !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector

      ! .. Local variables
      integer :: ii,jj,itr,imax,kk
      real(dblprec) :: fchkmax
      real(dblprec) :: init_ene !! Energy of the initial state
      real(dblprec) :: ftol_mRy,fc,fcinv
      character(3) :: num,num3
      character(35) :: num2,filn

      fc=mry/mub
      fcinv=mub/mry

      ftol_mRy = ftol*fc

      ! Call to allocate and initialize the arrays needed for the VPO calculation
      call allocate_VPO_arrays(1,Natom,Mensemble)

      ! Initialization of the variables
      fchkmax=1+ftol_mRy
      itr=1
      !-------------------------------------------------------------------------
      ! Pre calculation of the field and energy
      !-------------------------------------------------------------------------
      ! Calculation of the effective field
      call effective_field(Natom,Mensemble,1,nHam,emomM,mmom,    &
         external_field,time_external_field,beff,beff1,beff2,OPT_flag,              &
         max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations, &
         constellationsNeighType,energy,Num_macro,cell_index,emomM_macro, &
         macro_nlistsize,NA,N1,N2,N3)

      ! Calculation of the total energy
      call calc_energy(nHam,1,Natom,Nchmax,&
         conf_num,Mensemble,Natom,Num_macro,1,plotenergy,       &
         0.0_dblprec,delta_t,do_lsf,lsf_field,    &
         lsf_interpolate,'N',simid,cell_index,macro_nlistsize,mmom,emom,emomM,      &
         emomM_macro,external_field,time_external_field,max_no_constellations,      &
         maxNoConstl,unitCellType,constlNCoup,constellations,OPT_flag,              &
         constellationsNeighType,energy,NA,N1,N2,N3)

      ! Storing the energies for the initial and final images
      ! In VPO only the initial and final states are relaxed

      call convert_force(Natom,Mensemble,mmom,emom,beff)

      ! Storing the energies for the initial image
      init_ene = ene%energy(1)
      filn = 'force_min.'//trim(adjustl(simid))//'.out'
      open(ofileno,file = filn,access = 'sequential',action='write',status='replace')
      write(ofileno,'(a12)',advance='no')"# Iter."
      do ii=1,Mensemble
         write(ofileno,'(2x,a16,2x,a16)',advance='no')"Force","E_final-E_ini"
      enddo
      close(ofileno)

      fchkmax = 0.0_dblprec
      imax = 1
      do ii=1,Mensemble
         fchk(ii) = 0.0_dblprec
         do jj=1,Natom
            do kk=1,3
               if (abs(beff(kk,jj,ii))>fchk(ii)) then
                  fchk(ii) = abs(beff(kk,jj,ii))
               end if
            end do
         end do
         if (fchk(ii)>fchkmax) then
            fchkmax = fchk(ii)
            imax = ii
         end if
      end do

      open(ofileno,file = filn, access = 'sequential', action = 'write',            &
         status = 'old',position = 'append')
      write(ofileno,'(i12)',advance = 'no') itr
      do ii=1,Mensemble
         write(ofileno,'(2x,es16.8E3,2x,es16.8E3)',advance = 'no')fchk(ii)*fcinv,   &
            (ene%energy(ii)-init_ene)*fcinv
      end do
      close(ofileno)

      call save_path(Natom,Mensemble,simid,2,emom,mmom,mode,do_mom_legacy)
      !------------------------------------------------------------------------------
      ! MAIN LOOP
      !------------------------------------------------------------------------------
      do while ((fchkmax>ftol_mRy).and.(itr<=itrmax))
         do ii=1,Mensemble
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
         call effective_field(Natom,Mensemble,1,nHam,emomM,mmom, &
            external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
            max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
            constellations,constellationsNeighType,energy,Num_macro,      &
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

         call convert_force(Natom,Mensemble,mmom,emom,beff)
         do ii=1,Mensemble
            call quick_min_vel(Natom,vel(:,:,ii),beff0(:,:,ii),beff(:,:,ii),        &
               ax(:,:,ii),ang(:,ii),mass,delta_t)
         end do

         call project_vel_gneb(Natom,Mensemble,vel,beff)

         fchkmax = 0.0_dblprec
         imax = 1
         do ii=1,Mensemble
            fchk(ii) = 0.0_dblprec
            do jj=1,Natom
               do kk=1,3
                  if (abs(beff(kk,jj,ii))>fchk(ii)) then
                     fchk(ii) = abs(beff(kk,jj,ii))
                  end if
               end do
            end do
            if (fchk(ii)>fchkmax) then
               fchkmax = fchk(ii)
               imax = ii
            end if
         end do

         itr=itr+1

         if (mod(itr,every).eq.0) then
            ! Calculation of the total energy
            call calc_energy(nHam,itr,Natom,Nchmax,&
               conf_num,Mensemble,Natom,Num_macro,1,plotenergy, &
               0.0_dblprec,delta_t,do_lsf,        &
               lsf_field,lsf_interpolate,'N',simid,cell_index,macro_nlistsize,mmom, &
               emom,emomM,emomM_macro,external_field,time_external_field,           &
               max_no_constellations,maxNoConstl,unitCellType,constlNCoup,          &
               constellations,OPT_flag,constellationsNeighType,energy,NA,N1,N2,N3)

            write(num,'(i3)') idnint(real(itr,dblprec)/real(itrmax,dblprec)*100.0_dblprec)
            write(num2,'(es16.8E3)') fchkmax*mub/mry
            write(num3,'(i3)') imax
            write(*,'(2x,7a)')'IP   ',trim(adjustl(num)),'% of itrmax.   ','fchk: ',&
               trim(adjustl(num2)),'   imax: ',trim(adjustl(num3))

            open(ofileno,file = filn, access = 'sequential', action = 'write',      &
               status = 'old',position = 'append')
            write(ofileno,'(i12)',advance = 'no') itr
            do ii=1,Mensemble
               write(ofileno,'(2x,es16.8E3,2x,es16.8E3)',advance = 'no')            &
                  fchk(ii)*fcinv,(ene%energy(ii)-init_ene)*fcinv
            end do
            close(ofileno)

            call save_path(Natom, Mensemble, simid,2,emom,mmom,mode,do_mom_legacy)
         end if
      end do

      if (itr>itrmax) then
         write(*,'(a)') 'WARNING: exceeded maximum iterations in GNEB VPO'
      end if
      ! Call to deallocate the arrays needed for the VPO calculation
      call allocate_VPO_arrays(-1,Natom,Mensemble)

   end subroutine vpo_min

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: save_ifmom
   !> @brief
   !> Print initial and final configuration to file
   !> @author
   !> Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine save_ifmom(Natom,Mensemble,simid,prn_mode,mmom,emom,mode,do_mom_legacy)
      !
      use Restart, only: prn_mag_conf
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: prn_mode  !< 0 - states before relaxation, 1 - states after relaxation
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=1), intent(in) :: mode            !< Type of simulation
      character(len=8), intent(in) :: simid           !< Name of simulation
      character(len=1), intent(in)  :: do_mom_legacy  !< Flag to print/read moments in legacy output
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector

      integer :: i, j
      character(len=35) :: filn

      !------------------------------------------------------------------------------
      ! Print the initial and final configurations for GNEB in the new format
      !------------------------------------------------------------------------------
      if (do_mom_legacy.ne.'Y') then
         if (prn_mode==0) then
            call prn_mag_conf(Natom,0,2,'M',simid,mmom(:,1:Mensemble:(Mensemble-1)),&
               emom(:,:,1:Mensemble:(Mensemble-1)),'_if_in',mode)
         elseif (prn_mode==1) then
            call prn_mag_conf(Natom,0,2,'M',simid,mmom(:,1:Mensemble:(Mensemble-1)),&
               emom(:,:,1:Mensemble:(Mensemble-1)),'_if',mode)
         else
            write (*,*) "Error writing initial and final states to file"
            stop
         end if
      else
      !------------------------------------------------------------------------------
      ! Print the initial and final configurations for GNEB in the legacy format
      !------------------------------------------------------------------------------
         if (prn_mode==0) then
            filn = 'moment_if_in.'//trim(adjustl(simid))//'.out'
         elseif (prn_mode==1) then
            filn = 'moment_if.'//trim(adjustl(simid))//'.out'
         else
            write (*,*) "Error writing initial and final states to file"
            stop
         end if
         open(ofileno, file=filn, access = 'sequential',action = 'write', status = 'replace')
         do i=1, Mensemble,(Mensemble-1)
            do j=1, Natom
               write (ofileno,10002) i, j, emom(1,j,i), emom(2,j,i), emom(3,j,i), mmom(j,i)
            end do
         end do
         close(ofileno)
         return
         write (*,*) "Error writing initial and final states to file"
      endif
10002 format (i8,2x,i8,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)

   end subroutine save_ifmom

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_VPO_arrays
   !> @brief Allocation/deallocation of the arrays needed for the VPO calculation
   !> @author Jonathan Chico
   !> @date 27.03.2018
   !----------------------------------------------------------------------------
   subroutine allocate_VPO_arrays(flag,Natom,Mensemble)

      implicit none

      integer, intent(in) :: flag      !< Flag to allocate/deallocate arrays
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles

      integer :: i_stat,i_all

      if (flag>0) then
         allocate(fchk(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(fchk))*kind(fchk),'fchk','allocate_VPO_arrays')
         fchk=0.0_dblprec
         allocate(vel(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(vel))*kind(vel),'vel','allocate_VPO_arrays')
         vel=0.0_dblprec
         allocate(coo(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(coo))*kind(coo),'coo','allocate_VPO_arrays')
         coo=0.0_dblprec
         allocate(ax(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ax))*kind(ax),'ax','allocate_VPO_arrays')
         ax=0.0_dblprec
         allocate(ang(Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ang))*kind(ang),'ang','allocate_VPO_arrays')
         ang=0.0_dblprec
         allocate(beff(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(beff))*kind(beff),'beff','allocate_VPO_arrays')
         beff=0.0_dblprec
         allocate(beff1(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(beff1))*kind(beff1),'beff1','allocate_VPO_arrays')
         beff1=0.0_dblprec
         allocate(beff2(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(beff2))*kind(beff2),'beff2','allocate_VPO_arrays')
         beff2=0.0_dblprec
         allocate(beff0(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(beff0))*kind(beff0),'beff0','allocate_VPO_arrays')
         beff0=0.0_dblprec
         allocate(time_external_field(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(time_external_field))*kind(time_external_field),'time_external_field','allocate_VPO_arrays')
         time_external_field=0.0_dblprec
      else

         i_all=-product(shape(ax))*kind(ax)
         deallocate(ax,stat=i_stat)
         call memocc(i_stat,i_all,'ax','allocate_VPO_arrays')

         i_all=-product(shape(vel))*kind(vel)
         deallocate(vel,stat=i_stat)
         call memocc(i_stat,i_all,'vel','allocate_VPO_arrays')

         i_all=-product(shape(coo))*kind(coo)
         deallocate(coo,stat=i_stat)
         call memocc(i_stat,i_all,'coo','allocate_VPO_arrays')

         i_all=-product(shape(ang))*kind(ang)
         deallocate(ang,stat=i_stat)
         call memocc(i_stat,i_all,'ang','allocate_VPO_arrays')

         i_all=-product(shape(fchk))*kind(fchk)
         deallocate(fchk,stat=i_stat)
         call memocc(i_stat,i_all,'fchk','allocate_VPO_arrays')

         i_all=-product(shape(beff))*kind(beff)
         deallocate(beff,stat=i_stat)
         call memocc(i_stat,i_all,'beff','allocate_VPO_arrays')

         i_all=-product(shape(beff0))*kind(beff0)
         deallocate(beff0,stat=i_stat)
         call memocc(i_stat,i_all,'beff0','allocate_VPO_arrays')

         i_all=-product(shape(beff1))*kind(beff1)
         deallocate(beff1,stat=i_stat)
         call memocc(i_stat,i_all,'beff1','allocate_VPO_arrays')

         i_all=-product(shape(beff2))*kind(beff2)
         deallocate(beff2,stat=i_stat)
         call memocc(i_stat,i_all,'beff2','allocate_VPO_arrays')

         i_all=-product(shape(time_external_field))*kind(time_external_field)
         deallocate(time_external_field,stat=i_stat)
         call memocc(i_stat,i_all,'time_external_field','allocate_VPO_arrays')

      endif


   end subroutine allocate_VPO_arrays

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: vpo_min_single
   !> @brief
   !> Evolves the system to the nearest energy minimum according to the VPO algorithm
   !> Modified from original vpo_min to work outside the GNEB image concept (AB)
   !> @author
   !> Pavel Bessarab, Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine vpo_min_single(nHam,every,Natom,Nchmax, &
      itrmax,OPT_flag,conf_num,Mensemble,Num_macro,plotenergy,  &
      max_no_constellations,ftol,mass,delta_t,simid,do_lsf,lsf_field, &
      lsf_interpolate,cell_index,macro_nlistsize,mmom,emom,     &
      emomM_macro,external_field,maxNoConstl,unitCellType,constlNCoup,              &
      constellations,constellationsNeighType,energy,emomM,NA,N1,N2,N3,mode,         &
      do_mom_legacy)

      implicit none

      integer, intent(in) :: NA           !< Number of atoms in one cell
      integer, intent(in) :: N1           !< Number of cell repetitions in x direction
      integer, intent(in) :: N2           !< Number of cell repetitions in y direction
      integer, intent(in) :: N3           !< Number of cell repetitions in z direction
      integer, intent(in) :: nHam         !< Number of atoms in Hamiltonian
      integer, intent(in) :: every        !< Save path every 'every' step
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Nchmax       !< Max number of chemical components on each site in cell
      integer, intent(in) :: itrmax       !< Maximum number of iterations
      logical, intent(in) :: OPT_flag
      integer, intent(in) :: conf_num     !< number of configurations for LSF
      integer, intent(in) :: Mensemble    !< Number of images in GNEB calculations
      integer, intent(in) :: Num_macro    !< Number of macrocells in the system
      integer, intent(in) :: plotenergy   !< Calculate and plot energy (0/1)
      integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
      real(dblprec), intent(in) :: ftol               !< Tolerance
      real(dblprec), intent(in) :: mass               !< mass of the point
      real(dblprec), intent(in) :: delta_t            !< timestep
      character(len=1), intent(in) :: mode            !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid           !< Name of simulation
      character(len=1), intent(in) :: do_lsf          !< Including LSF energy
      character(len=1), intent(in) :: lsf_field       !< LSF field contribution (Local/Total)
      character(len=1), intent(in) :: do_mom_legacy    !< Flag to print/read moments in legacy output
      character(len=1), intent(in) :: lsf_interpolate !< Interpolate LSF or not
      integer, dimension(Natom), intent(in) :: cell_index    !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field
      integer, dimension(:), intent(in) :: maxNoConstl
      integer, dimension(:,:), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup
      real(dblprec), dimension(:,:,:), intent(in) :: constellations
      integer, dimension(:,:,:), intent(in) :: constellationsNeighType

      ! .. In/out variables
      real(dblprec), intent(inout) :: energy
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom     !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector

      ! .. Local variables
      integer :: ii,jj,itr,imax,kk
      real(dblprec) :: fchkmax
      real(dblprec) :: init_ene !! Energy of the initial state
      real(dblprec) :: ftol_mRy,fc,fcinv
      character(3) :: num,num3
      character(35) :: num2,filn

      fc=mry/mub
      fcinv=mub/mry

      ftol_mRy = ftol*fc

      ! Call to allocate and initialize the arrays needed for the VPO calculation
      call allocate_VPO_arrays(1,Natom,Mensemble)

      ! Initialization of the variables
      fchkmax=1+ftol_mRy
      itr=1
      !-------------------------------------------------------------------------
      ! Pre calculation of the field and energy
      !-------------------------------------------------------------------------
      ! Calculation of the effective field
      call effective_field(Natom,Mensemble,1,nHam,emomM,mmom,    &
         external_field,time_external_field,beff,beff1,beff2,OPT_flag,              &
         max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations, &
         constellationsNeighType,energy,Num_macro,cell_index,emomM_macro, &
         macro_nlistsize,NA,N1,N2,N3)

      ! Convert effective fields to VPO force equivalents
      call convert_force(Natom,Mensemble,mmom,emom,beff)

      ! Storing the energies for the initial image
      init_ene = energy/Natom
      filn = 'force_min.'//trim(adjustl(simid))//'.out'
      open(ofileno,file = filn,access = 'sequential',action='write',status='replace')
      write(ofileno,'(a12)',advance='no')"# Iter."
      do ii=1,Mensemble
         write(ofileno,'(2x,a16,2x,a16)',advance='no')"Force","E_final-E_ini"
      enddo
      close(ofileno)

      fchkmax = 0.0_dblprec
      imax = 1
      do ii=1,Mensemble
         fchk(ii) = 0.0_dblprec
         do jj=1,Natom
            do kk=1,3
               if (abs(beff(kk,jj,ii))>fchk(ii)) then
                  fchk(ii) = abs(beff(kk,jj,ii))
               end if
            end do
         end do
         if (fchk(ii)>fchkmax) then
            fchkmax = fchk(ii)
            imax = ii
         end if
      end do

      open(ofileno,file = filn, access = 'sequential', action = 'write',            &
         status = 'old',position = 'append')
      write(ofileno,'(i12)',advance = 'no') itr
      do ii=1,Mensemble
         write(ofileno,'(2x,es16.8E3,2x,es16.8E3)',advance = 'no')fchk(ii)*fcinv,   &
            energy/Natom - init_ene
      end do
      close(ofileno)

      !------------------------------------------------------------------------------
      ! MAIN LOOP
      !------------------------------------------------------------------------------
      do while ((fchkmax>ftol_mRy).and.(itr<=itrmax))
         do ii=1,Mensemble
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
         call effective_field(Natom,Mensemble,1,nHam,emomM,mmom, &
            external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
            max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
            constellations,constellationsNeighType,energy,Num_macro,      &
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

         call convert_force(Natom,Mensemble,mmom,emom,beff)
         do ii=1,Mensemble
            call quick_min_vel(Natom,vel(:,:,ii),beff0(:,:,ii),beff(:,:,ii),        &
               ax(:,:,ii),ang(:,ii),mass,delta_t)
         end do

         call project_vel_gneb(Natom,Mensemble,vel,beff)

         fchkmax = 0.0_dblprec
         imax = 1
         do ii=1,Mensemble
            fchk(ii) = 0.0_dblprec
            do jj=1,Natom
               do kk=1,3
                  if (abs(beff(kk,jj,ii))>fchk(ii)) then
                     fchk(ii) = abs(beff(kk,jj,ii))
                  end if
               end do
            end do
            if (fchk(ii)>fchkmax) then
               fchkmax = fchk(ii)
               imax = ii
            end if
         end do

         itr=itr+1

         write(200,'(A,i7, A,F12.8,A,F10.4)') 'Iteration', itr, ' : energy =', energy/Natom, ' , |g|_inf =', fchkmax*mub/mry
         if (mod(itr,every).eq.0) then
            ! Calculation of the total energy
            call effective_field(Natom,Mensemble,1,nHam,emomM,mmom, &
               external_field,time_external_field,beff,beff1,beff2,OPT_flag,           &
               max_no_constellations,maxNoConstl,unitCellType,constlNCoup,             &
               constellations,constellationsNeighType,energy,Num_macro,      &
               cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

            write(num,'(i3)') idnint(real(itr,dblprec)/real(itrmax,dblprec)*100.0_dblprec)
            write(num2,'(es16.8E3)') fchkmax*mub/mry
            write(num3,'(i3)') imax
            write(*,'(2x,7a)')'IP   ',trim(adjustl(num)),'% of itrmax.   ','fchk: ',&
               trim(adjustl(num2)),'   imax: ',trim(adjustl(num3))

            open(ofileno,file = filn, access = 'sequential', action = 'write',      &
               status = 'old',position = 'append')
            write(ofileno,'(i12)',advance = 'no') itr
            do ii=1,Mensemble
               write(ofileno,'(2x,es16.8E3,2x,es16.8E3)',advance = 'no')fchk(ii)*fcinv,   &
                energy/Natom - init_ene
            end do
            close(ofileno)

         end if
      end do
      print *, 'Minimization done, energy ', energy/Natom, 'mRy/atom'

      if (itr>itrmax) then
         write(*,'(a)') 'WARNING: exceeded maximum iterations in GNEB VPO'
      end if
      ! Call to deallocate the arrays needed for the VPO calculation
      call allocate_VPO_arrays(-1,Natom,Mensemble)

   end subroutine vpo_min_single
end module energyminima
