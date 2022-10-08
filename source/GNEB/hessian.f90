!------------------------------------------------------------------------------------
! MODULE: Hessian
!> @brief Routines used for the calculation of the Hessian
!> @details The Hessian defined as \f$ H=\frac{\partial^2 \mathcal{H}}{\partial m_i^\mu \partial m_j^\nu}\f$
!> @notes Re-structured by Jonathan Chico
!> @author Pavel Bessarab
!------------------------------------------------------------------------------------
module Hessian
   use Parameters
   use Profiling
   use Ewaldmom
   use Constants

   implicit none

   real(dblprec), dimension(:,:), allocatable :: hess       !< Matrix of the second order derivatives of the Hamiltonian with respect to the magnetic moments
   real(dblprec), dimension(:,:), allocatable :: hess_proj
   real(dblprec), dimension(:,:), allocatable :: vec_sp
   real(dblprec), dimension(:,:), allocatable :: vec_curr

   public

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: Hessian wraper
   !> @brief Wrapper routine for the calculation of the Hessian for the calculation
   !> of the prefactor in Harmonic Transition State Theory (HTST)
   !> @details This wrapper calculated the Hessian matrix for the intial, final and
   !> saddle point states. The eigenvalues of the Hessian are then used to calculate
   !> the prefactor of the Arrhenius process, i.e. in KMC terms, the attempt frequency.
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine hessian_wrapper(N1,N2,N3,nim,nHam,Natom,do_dm,do_pd,do_bq,do_ring,do_chir, &
      do_dip,do_biqdm,Num_macro,do_jtensor,plotenergy,max_no_neigh,do_anisotropy,   &
      max_no_dmneigh,max_no_constellations,eig_0,BC1,BC2,BC3,simid,do_lsf,is_afm,   &
      mult_axis,exc_inter,lsf_field,do_hess_sp,do_hess_ini,do_hess_fin,             &
      lsf_interpolate,OPT_flag,aHam,taniso,nlistsize,dmlistsize,cell_index,         &
      macro_nlistsize,nlist,dmlist,C1,C2,C3,ene0,im_ene,coord,eaniso,kaniso,emomsp, &
      ncoup,mmom,emomM_macro,dm_vect,external_field,maxNoConstl,unitCellType,       &
      constellationsNeighType,constlNCoup,constellations,emomM,beff,beff1,beff2,NA, &
      emom)

      use HamiltonianActions, only : effective_field

      implicit none

      integer, intent(in) :: NA              !< Number of atoms in one cell
      integer, intent(in) :: N1              !< Number of cell repetitions in x direction
      integer, intent(in) :: N2              !< Number of cell repetitions in y direction
      integer, intent(in) :: N3              !< Number of cell repetitions in z direction
      integer, intent(in) :: nim             !< Number of images
      integer, intent(in) :: nHam            !< Number of atoms in Hamiltonian
      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: do_dm           !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd           !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq           !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_ring         !< Add four-spin ring (4SR) term to Hamiltonian (0/1)
      integer, intent(in) :: do_chir         !< Add scalar chirality exchange (CHIR) term to Hamiltonian (0/1)
      integer, intent(in) :: do_dip          !< Calculate dipole-dipole contribution (0/1)
      integer, intent(in) :: do_biqdm        !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: Num_macro       !< Number of macrocells in the system
      integer, intent(in) :: do_jtensor      !< Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
      integer, intent(in) :: plotenergy      !< Calculate and plot energy (0/1)
      integer, intent(in) :: max_no_neigh    !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: do_anisotropy   !< Read anisotropy data (1/0)
      integer, intent(in) :: max_no_dmneigh  !< Calculated number of neighbours with DM interactions
      integer, intent(in) :: max_no_constellations !< The maximum (global) length of the constellation matrix
      real(dblprec), intent(in) :: eig_0  !< Threshold value for "zero" eigenvalue
      character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
      character(len=8), intent(in) :: simid     !< Name of simulation
      character(len=1), intent(in) :: do_lsf    !< Including LSF energy
      character(len=1), intent(in) :: is_afm    !< Flag to indicate if the state is Antiferromagnetic
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      character(len=1), intent(in) :: exc_inter !< Interpolation of Jij (Y/N)
      character(len=1), intent(in) :: lsf_field    !< LSF field contribution (Local/Total)
      character(len=1), intent(in) :: do_hess_sp   !< Calculate the Hessian at the saddle point
      character(len=1), intent(in) :: do_hess_ini  !< Calculate the Hessian at the initial configuration
      character(len=1), intent(in) :: do_hess_fin  !< Calculate the Hessian at the fianl configuration
      character(len=1), intent(in) :: lsf_interpolate    !< Interpolate LSF or not
      logical, intent(in) :: OPT_flag !< Optimization flag
      integer, dimension(Natom), intent(in) :: aham      !< Hamiltonian look-up table
      integer, dimension(Natom), intent(in) :: taniso    !< Type of anisotropy (0-2)
      integer, dimension(nHam), intent(in) :: nlistsize  !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, dimension(Natom), intent(in) :: cell_index    !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist    !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist !< List of neighbours for DM
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(1), intent(in) :: ene0 !< Zero value of the energy
      real(dblprec), dimension(nim), intent(in) :: im_ene   !< Energies of the images
      real(dblprec), dimension(Natom,nim), intent(in)          :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(max_no_neigh,nHam), intent(in)  :: ncoup    !< Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom), intent(in)            :: coord    !< Coordinates of atoms
      real(dblprec), dimension(3,Natom), intent(in)            :: eaniso   !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in)            :: kaniso   !< Anisotropy constant
      real(dblprec), dimension(3,Natom), intent(in)            :: emomsp   !< Spin configuation at the saddle point
      real(dblprec), dimension(3,Num_macro,nim), intent(in)    :: emomM_macro    !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,max_no_dmneigh,nHam), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      real(dblprec), dimension(3,Natom,nim), intent(in)           :: external_field !< External magnetic field
      integer, dimension(:), intent(in) :: maxNoConstl !< Number of existing entries in for each ensemble in the constellation matrix
      integer, dimension(:,:), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      integer, dimension(:,:,:), intent(in) :: constellationsNeighType !< Every constellation atom's neighbour type atoms.  This will tell which values in the constellation to look for during the applyhamiltionian step
      real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup !< Couplings for saved constellations
      real(dblprec), dimension(:,:,:), intent(in) :: constellations !< Saved fixed unit cell configurations, these represent a configuration present in the domain with same configuration of unit cells in the neighborhood
      ! .. In/out variables
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: emomM    !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: emom    !< Current unit moment vector
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: beff    !< Current magnetic fields
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: beff1   !< Current magnetic fields
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: beff2   !< Current magnetic fields

      ! Local variables
      integer :: info      !< General flag of the status of the calculation
      integer :: info_sp   !< Flag of the status of the eigenvalues for the Saddle state
      integer :: info_ini  !< Flag of the status of the eigenvalues for the Initial state
      integer :: info_fin  !< Flag of the status of the eigenvalues for the Final state
      integer :: nvec_sp   !< Number of eignevectors for the Saddle point state
      integer :: nvec_ini  !< Number of eignevectors for the Initial point state
      integer :: nvec_fin  !< Number of eignevectors for the Final point state
      integer :: ii, jj,i_stat,i_all
      real(dblprec) :: denergy   !< Mean value of the Total energy
      real(dblprec) :: pref_if   !< Prefactor of the transition Intial to Final state
      real(dblprec) :: pref_fi   !< Prefactor of the transition Final to Initial state
      real(dblprec) :: expo_if   !< Exponential of the transition Intial to Final state
      real(dblprec) :: expo_fi   !< Exponential of the transition Final to Initial state
      real(dblprec) :: vol_sp    !< Volume of the Saddle point state
      real(dblprec) :: vol_ini   !< Volume of the Initial state
      real(dblprec) :: vol_fin   !< Volume of the Final state
      real(dblprec) :: fcinv
      character(len=35) :: beff_filn
      character(len=10) :: state
      real(dblprec), dimension(3,Natom,nim) :: tef !< External time-dependent magnetic field, zero here

      ! .. Local allocatable arrays
      integer, dimension(:), allocatable :: NNN !< Number of repetitions of the unit cell
      real(dblprec), dimension(:), allocatable :: eig_sp    !< Eigenvalues of the Saddle point state
      real(dblprec), dimension(:), allocatable :: eig_ini   !< Eigenvalues of the Initial state
      real(dblprec), dimension(:), allocatable :: eig_fin   !< Eigenvalues of the Final state
      real(dblprec), dimension(:), allocatable :: vel_par   !< Velocity parameters
      real(dblprec), dimension(:,:), allocatable :: bass !< Matrix of the lattice vectors
      real(dblprec), dimension(:,:), allocatable :: basis
      character(len=1), dimension(:), allocatable :: bou_con !< Boundary conditions

      info_ini = 0
      info_fin = 0
      info_sp  = 0
      nvec_ini = 0
      nvec_fin = 0
      nvec_sp  = 0

      vol_ini = 1.0_dblprec
      vol_fin = 1.0_dblprec
      vol_sp  = 1.0_dblprec

      ! Set the time dependent fields to zero
      tef=0.0_dblprec
      ! Allocation of local arrays
      allocate(bass(3,3),stat=i_stat)
      call memocc(i_stat,product(shape(bass))*kind(bass),'bass','hessian_wrapper')
      bass=0.0_dblprec
      allocate(basis(3*Natom,2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(basis))*kind(basis),'basis','hessian_wrapper')
      basis=0.0_dblprec
      allocate(NNN(3),stat=i_stat)
      call memocc(i_stat,product(shape(NNN))*kind(NNN),'NNN','hessian_wrapper')
      NNN=0
      allocate(bou_con(3),stat=i_stat)
      call memocc(i_stat,product(shape(bou_con))*kind(bou_con),'bou_con','hessian_wrapper')

      fcinv=mub/mry

      do ii=1,3
         bass(ii,1) = C1(ii)
         bass(ii,2) = C2(ii)
         bass(ii,3) = C3(ii)
      end do
      bou_con(1) = BC1
      bou_con(2) = BC2
      bou_con(3) = BC3
      NNN(1) = N1
      NNN(2) = N2
      NNN(3) = N3

      ! Set the time dependent fields to zero
      tef=0.0_dblprec
      !------------------------------------------------------------------------------
      !Calculate Hessian matrix
      !------------------------------------------------------------------------------
      write (*,'(1x,a)') "     Hessian calculations in progress     "
      write (*,'(1x,a)') "------------------------------------------"
      !------------------------------------------------------------------------------
      ! Initial state Hessian calculation
      !------------------------------------------------------------------------------
      if (do_hess_ini=='Y') then
         allocate(eig_ini(2*Natom),stat=i_stat)
         call memocc(i_stat,product(shape(eig_ini))*kind(eig_ini),'eig_ini','hessian_wrapper')
         eig_ini=0.0_dblprec
         state='Initial'
         call timing(0,'Hamiltonian   ','ON')
         ! Calculate effective fields at the initial point
         beff(:,:,1)=0.0_dblprec
         call effective_field(Natom,1,1,Natom,emomM(:,:,1),mmom(:,1), &
            external_field(:,:,1),tef(:,:,1),beff(:,:,1),beff1(:,:,1),beff2(:,:,1),   &
            OPT_flag,max_no_constellations,maxNoConstl,unitCellType,constlNCoup,      &
            constellations,constellationsNeighType,denergy,Num_macro,       &
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
         call timing(0,'Hamiltonian   ','OF')
         call hessian_eigenvalues(nHam,Natom,do_dm,max_no_neigh,max_no_dmneigh,     &
            eig_0,is_afm,simid,state,NNN,aham,taniso,nlistsize,dmlistsize,nlist,    &
            dmlist,bass,basis,coord,eaniso,kaniso,ncoup,emomM(:,:,1),dm_vect,       &
            bou_con,nvec_ini,info_ini,vol_ini,eig_ini,beff(:,:,1),emom(:,:,1))
      end if
      !------------------------------------------------------------------------------
      ! Final state Hessian calculation
      !------------------------------------------------------------------------------
      if (do_hess_fin=='Y') then
         allocate(eig_fin(2*Natom),stat=i_stat)
         call memocc(i_stat,product(shape(eig_fin))*kind(eig_fin),'eig_fin','hessian_wrapper')
         eig_fin=0.0_dblprec
         state='Final'
         call timing(0,'Hamiltonian   ','ON')
         ! Calculate effective fields at the final point
         beff(:,:,nim)=0.0_dblprec
         call effective_field(Natom,1,1,Natom,emomM(:,:,nim),mmom(:,nim),     &
            external_field(:,:,nim),tef(:,:,nim),beff(:,:,nim),beff1(:,:,nim),beff2(:,:,nim), &
            OPT_flag,max_no_constellations,maxNoConstl,unitCellType,constlNCoup,    &
            constellations,constellationsNeighType,denergy,Num_macro,     &
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
         call timing(0,'Hamiltonian   ','OF')
         call hessian_eigenvalues(nHam,Natom,do_dm,max_no_neigh,max_no_dmneigh,     &
            eig_0,is_afm,simid,state,NNN,aham,taniso,nlistsize,dmlistsize,nlist,    &
            dmlist,bass,basis,coord,eaniso,kaniso,ncoup,emomM(:,:,nim),dm_vect,     &
            bou_con,nvec_fin,info_fin,vol_fin,eig_fin,beff(:,:,nim),emom(:,:,nim))
      end if
      !------------------------------------------------------------------------------
      ! Saddle point Hessian calculation
      !------------------------------------------------------------------------------
      if (do_hess_sp=='Y') then
         ! Find the moments at the saddle point
         do ii=1,Natom
            emomM(:,ii,1) = emomsp(:,ii)*mmom(ii,1)
         end do

         call timing(0,'Hamiltonian   ','ON')
         ! Calculate effective fields at the saddle point
         call effective_field(Natom,1,1,Natom,emomM(:,:,1),mmom(:,1),&
            external_field(:,:,1),tef(:,:,1),beff(:,:,1),beff1(:,:,1),beff2(:,:,1),  &
            OPT_flag,max_no_constellations,maxNoConstl,unitCellType,constlNCoup,     &
            constellations,constellationsNeighType,denergy,Num_macro,      &
            cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
         call timing(0,'Hamiltonian   ','OF')
         ! Save effective fields at the SP
         beff_filn = 'beff_sp.'//trim(adjustl(simid))//'.out'
         open(ofileno, file=beff_filn, access = 'sequential',action = 'write', status = 'replace')
         write(ofileno,'(a8,2x,a16,2x,a16,2x,a16)') "# Site","B_x","B_y","B_z"
         do jj=1,Natom
            write(ofileno,'(i8,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)',advance = 'yes') &
            jj,beff(1,jj,1)*fcinv,beff(2,jj,1)*fcinv,beff(3,jj,1)*fcinv
         end do
         close(ofileno)

         allocate(eig_sp(2*Natom),stat=i_stat)
         call memocc(i_stat,product(shape(eig_sp))*kind(eig_sp),'eig_sp','hessian_wrapper')
         eig_sp=0.0_dblprec
         state='Saddle'
         ! Calculate the eigenvalues and eigenvectors of the Hessian at the saddle point
         call hessian_SP_eigenvalues(nHam,Natom,do_dm,max_no_neigh,max_no_dmneigh,  &
            eig_0,is_afm,simid,state,NNN,aham,taniso,nlistsize,dmlistsize,nlist,    &
            dmlist,bass,basis,coord,eaniso,kaniso,ncoup,emomM(:,:,1),dm_vect,       &
            bou_con,nvec_sp,info_sp,vol_sp,eig_sp,vel_par,beff(:,:,1),emomsp(:,:),  &
            mmom(:,1))

      endif
      !------------------------------------------------------------------------------
      ! If the calculation of the eigenvalues has been successful for the initial, and saddle point
      ! states then calculate the rate of the process of going from one state to the other
      !------------------------------------------------------------------------------
      if ((do_hess_ini=='Y').and.(do_hess_sp=='Y').and.(info_sp==0).and.(info_ini==0)) then
         write (*,'(2x,a)',advance='no') "Calculation of the prefactor for intial ---> final state transition "
         call calc_pref(Natom,nvec_ini,nvec_sp-1,vol_sp,vol_ini,eig_ini,eig_sp,vel_par,pref_if,expo_if,info)
         write (*,'(x,a)') "done"
         if (info.ne.0) then
            write(*,'(2x, a)') 'WARNING: Failed to calculate the prefactor for intial ---> final state transition!'
         else
            open(ofileno, file='rate_if.'//trim(adjustl(simid))//'.out', access = 'sequential',action = 'write', status = 'replace')
            write(ofileno,'(a16,2x,a16,2x,a16,2x,a16,2x,a16)') "# Barrier","Prefactor","Exp Trans","Vol_ini","Vol_sp"
            write(ofileno,'(es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)') ene0(1)-im_ene(1),pref_if,expo_if,vol_ini,vol_sp
            close(ofileno)
         end if
         write (*,'(2x,a)') "-----------------------------------------"
      else
         open(ofileno, file='rate_if.'//trim(adjustl(simid))//'.out', access = 'sequential',action = 'write', status = 'replace')
         write(ofileno,'(a16)') "# Barrier"
         write(ofileno,'(es16.8E3)') ene0(1)-im_ene(1)
         close(ofileno)
      end if

      !------------------------------------------------------------------------------
      ! If the calculation of the eigenvalues has been successful for the final, and saddle point
      ! states then calculate the rate of the process of going from one state to the other
      !------------------------------------------------------------------------------
      if ((do_hess_fin=='Y').and.(do_hess_sp=='Y').and.(info_sp==0).and.(info_fin==0)) then
         write (*,'(2x,a)',advance='no') "Calculation of the prefactor for final ---> initial state transition"
         call calc_pref(Natom,nvec_fin,nvec_sp-1,vol_sp,vol_fin,eig_fin,eig_sp,vel_par,pref_fi,expo_fi,info)
         write (*,'(x,a)') "done"
         if (info.ne.0) then
            write(*,'(2x,a)') 'WARNING: Failed to calculate the prefactor for final ---> initial state transition!'
         else
            open(ofileno, file='rate_fi.'//trim(adjustl(simid))//'.out', access = 'sequential',action = 'write', status = 'replace')
            write(ofileno,'(a16,2x,a16,2x,a16,2x,a16,2x,a16)') "# Barrier","Prefactor","Exp Trans","Vol_fin","Vol_sp"
            write(ofileno,'(es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)') ene0(1)-im_ene(nim),pref_fi,expo_fi,vol_fin,vol_sp
            close(ofileno)
         end if
         write (*,'(2x,a)') "-----------------------------------------"
      else
         open(ofileno, file='rate_fi.'//trim(adjustl(simid))//'.out', access = 'sequential',action = 'write', status = 'replace')
         write(ofileno,'(a16)') "# Barrier"
         write(ofileno,'(es16.8E3)') ene0(1)-im_ene(nim)
         close(ofileno)

      end if

      write(*,'(2x,a,x,es16.8)') "Initial State Volume :", vol_ini
      write(*,'(2x,a,x,es16.8)') "Final State Volume   :", vol_fin
      write(*,'(2x,a,x,es16.8)') "Saddle State Volume  :", vol_sp

      ! Deallocation of the needed arrays for the Hessian wrapper
      if (allocated(eig_ini)) then
         i_all=-product(shape(eig_ini))*kind(eig_ini)
         deallocate(eig_ini,stat=i_stat)
         call memocc(i_stat,i_all,'eig_ini','hessian_wrapper')
      endif
      if (allocated(eig_fin)) then
         i_all=-product(shape(eig_fin))*kind(eig_fin)
         deallocate(eig_fin,stat=i_stat)
         call memocc(i_stat,i_all,'eig_fin','hessian_wrapper')
      endif
      if (allocated(eig_sp)) then
         i_all=-product(shape(eig_sp))*kind(eig_sp)
         deallocate(eig_sp,stat=i_stat)
         call memocc(i_stat,i_all,'eig_sp','hessian_wrapper')
      endif
      if (allocated(vel_par)) then
         i_all=-product(shape(vel_par))*kind(vel_par)
         deallocate(vel_par,stat=i_stat)
         call memocc(i_stat,i_all,'vel_par','hessian_wrapper')
      end if
      i_all=-product(shape(bass))*kind(bass)
      deallocate(bass,stat=i_stat)
      call memocc(i_stat,i_all,'bass','hessian_wrapper')
      i_all=-product(shape(basis))*kind(basis)
      deallocate(basis,stat=i_stat)
      call memocc(i_stat,i_all,'basis','hessian_wrapper')
      i_all=-product(shape(bou_con))*kind(bou_con)
      deallocate(bou_con,stat=i_stat)
      call memocc(i_stat,i_all,'bou_con','hessian_wrapper')
      i_all=-product(shape(NNN))*kind(NNN)
      deallocate(NNN,stat=i_stat)
      call memocc(i_stat,i_all,'NNN','hessian_wrapper')

   end subroutine hessian_wrapper

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: hessian_eigenvalues
   !> @brief Wrapper routine for the determination of the Hessian and its eigenvalues
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine hessian_eigenvalues(nHam,Natom,do_dm,max_no_neigh,max_no_dmneigh,     &
      eig_0,is_afm,simid,state,NNN,aham,taniso,nlistsize,dmlistsize,nlist,dmlist,   &
      bass,basis,coord,eaniso,kaniso,ncoup,emomM,dm_vect,bou_con,nvec_curr,         &
      info_curr,vol_curr,eig_curr,beff,emom)

      implicit none

      integer, intent(in) :: nham            !< number of atoms in hamiltonian
      integer, intent(in) :: Natom           !< number of atoms in system
      integer, intent(in) :: do_dm           !< add dzyaloshinskii-moriya (dm) term to hamiltonian (0/1)
      integer, intent(in) :: max_no_neigh    !< calculated maximum of neighbours for exchange
      integer, intent(in) :: max_no_dmneigh  !< calculated number of neighbours with dm interactions
      real(dblprec), intent(in) :: eig_0     !< Threshold value for "zero" eigenvalue
      character(len=8), intent(in)  :: simid    !< name of simulation
      character(len=10), intent(in) :: state    !< which is the state that is being calculated (initial/final)
      character(len=1), intent(in)  :: is_afm   !< Flag to indicate if the state is Antiferromagnetic
      integer, dimension(3), intent(in) :: NNN  !< Number of repetitions of the unit cell
      integer, dimension(Natom), intent(in) :: aham      !< hamiltonian look-up table
      integer, dimension(Natom), intent(in) :: taniso    !< type of anisotropy (0-2)
      integer, dimension(nham), intent(in) :: nlistsize  !< size of neighbour list for heisenberg exchange couplings
      integer, dimension(nham), intent(in) :: dmlistsize !< size of neighbour list for dm
      integer, dimension(max_no_neigh,natom), intent(in) :: nlist    !< neighbour list for heisenberg exchange couplings
      integer, dimension(max_no_dmneigh,natom), intent(in) :: dmlist !< list of neighbours for dm
      real(dblprec), dimension(3,3), intent(in) :: bass !< Matrix of the lattice vectors
      real(dblprec), dimension(3,Natom), intent(in) :: beff    !< current magnetic moment vector
      real(dblprec), dimension(3,Natom), intent(in) :: emomm   !< current magnetic moment vector
      real(dblprec), dimension(3,Natom), intent(in) :: emom   !< current magnetic moment vector
      real(dblprec), dimension(3,Natom), intent(in) :: coord   !< Coordinates of atoms
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso  !< unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso  !< anisotropy constant
      real(dblprec), dimension(max_no_neigh,nham), intent(in) :: ncoup  !< heisenberg exchange couplings
      real(dblprec), dimension(3,max_no_dmneigh,nham), intent(in) :: dm_vect !< dzyaloshinskii-moriya exchange vector
      character(len=1), dimension(3), intent(in) :: bou_con !< Boundary conditions
      ! .. Output variables
      integer, intent(out) :: nvec_curr   !< Number of eigenvectors of the current image
      integer, intent(out) :: info_curr   !< Information flag of the current image
      !.. In/out variables
      real(dblprec), intent(inout) :: vol_curr  !< Volume of the zero modes in the current image
      real(dblprec), dimension(2*Natom), intent(inout) :: eig_curr   !< Eigenvalues of the current image
      real(dblprec), dimension(3*Natom,2*Natom), intent(inout) :: basis
      ! .. Local variables
      integer :: i_stat, i_all
      integer :: ii, jj, kk
      integer :: info
      real(dblprec) :: fc, fcinv, tmp
      character(len=256) :: msg
      character(len=35) :: num
      logical :: do_hess_calc

      do_hess_calc=.true.

      fc=mry/mub
      fcinv=mub/mry

      allocate(hess(3*Natom,3*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(hess))*kind(hess),'hess','hessian_eigenvalues')
      hess=0.0_dblprec

      allocate(hess_proj(2*Natom,2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(hess_proj))*kind(hess_proj),'hess_proj','hessian_eigenvalues')
      hess_proj=0.0_dblprec

      write(*,'(2x,a,x,a)') state, 'state'
      write(*,'(2x,a)') 'Evaluation of the Hessian   '
      call calc_hess_lambda(Natom, nHam,emomM,beff, max_no_neigh, nlistsize,nlist,  &
         ncoup,aham,do_dm, max_no_dmneigh,dmlistsize, dmlist, dm_vect,taniso,eaniso,&
         kaniso, hess)
      write(*,'(2x,a)') 'Done'

      write(*,'(2x,a)',advance='no') 'Projection of the Hessian on the tangent space '
      call calc_basis(Natom,emom,basis)

      call project_hess(Natom,basis,hess,hess_proj)

      write (*,'(x,a)') "done"
      write (*,'(2x,a)',advance='no') "Evaluation of the eigenvalues "
      call calc_hess_eig(Natom,hess_proj,eig_curr,info)
      write (*,'(2x,a)') "Done"
      if (info.ne.0) then
         write(*,'(2x,a)') 'WARNING: Failed to calculate eigenvalues!'
         info_curr = info
      end if

      open(ofileno, file='eig_'//trim(adjustl(state))//'.'//trim(adjustl(simid))//'.out',&
      access = 'sequential',action = 'write', status = 'replace')
      write(ofileno,'(a8,2x,a16)') "# Ent.","Eigenvalue"
      do ii=1,2*Natom
         write(ofileno,'(i8,2x,es16.8E3)',advance = 'yes') ii,eig_curr(ii)*fcinv
      end do
      close(ofileno)

      ! Check if these conditions are fulfilled for calculation of the extremum points
      ii=1
      do while (ii<=2*Natom .and. do_hess_calc)
         if ((eig_curr(ii)<0.0_dblprec).and.(abs(eig_curr(ii)*fcinv)>eig_0)) then
            write(num,'(es16.8E3)') eig_curr(ii)*fcinv
            msg = 'WARNING: '//trim(adjustl(state))//' state is not a minimum, because '//trim(adjustl(num))//'<0!'
            write(*,'(2x,a)') trim(adjustl(msg))
            info_curr = 1
            do_hess_calc=.false.
         end if
         ii=ii+1
      end do

      ! If the tests are passed continue with the calculation
      if (do_hess_calc) then
         write (*,'(2x,a,x,a,x,a)') 'Checking for zero modes at the',state,'point'
         nvec_curr = 0
         do ii=1,2*Natom
            if (abs(eig_curr(ii)*fcinv).le.eig_0) then
               nvec_curr = nvec_curr+1
            end if
         end do

         if (nvec_curr==0) then
            write (*,'(2x,a)') "No zero modes have been found."
         else
            allocate(vec_curr(3*Natom,nvec_curr),stat=i_stat)
            call memocc(i_stat,product(shape(vec_curr))*kind(vec_curr),'vec_curr','hessian_eigenvalues')
            vec_curr=0.0_dblprec

            write(num,'(i8)') nvec_curr
            msg = trim(adjustl(num))//' zero modes have been found.'
            write (*,'(2x,a)') trim(adjustl(msg))
            write (*,'(2x,a)',advance='no') "Calculation of zero-mode eigenvectors "
            call calc_hess_V(Natom,nvec_curr,hess_proj,eig_curr,basis,vec_curr,info)
            write (*,'(x,a)') "done"
            if (info.ne.0) then
               write(*,'(2x,a)') 'WARNING: Failed to calculate zero-mode eigenvectors!'
               info_curr = info
            else
               do kk=1,nvec_curr
                  tmp = 0.0_dblprec
                  !$omp parallel do default(shared), private(ii,jj), reduction(+:tmp), collapse(2)
                  do ii=1,3*Natom
                     do jj=1,3*Natom
                        tmp = tmp + vec_curr(ii,kk)*hess(ii,jj)*vec_curr(jj,kk)
                     end do
                  end do
                  !$omp end parallel do
                  write(*,'(3x,a,x,es16.8,x,es16.8)')'compare eig:',eig_curr(kk)*fcinv,tmp*fcinv
               end do
            end if
            if (is_afm=='Y') then
               call calc_zero_vol_afm(Natom,NNN,emomM,coord,bass,bou_con,vol_curr)
            else
               call calc_zero_vol(Natom,NNN,emomM,coord,bass,bou_con,vol_curr)
            end if

            i_all=-product(shape(vec_curr))*kind(vec_curr)
            deallocate(vec_curr,stat=i_stat)
            call memocc(i_stat,i_all,'vec_curr','hessian_eigenvalues')
         end if
      endif

      write (*,'(2x,a)') "-----------------------------------------"

      eig_curr(:) = eig_curr(:)*fcinv

      i_all=-product(shape(hess))*kind(hess)
      deallocate(hess,stat=i_stat)
      call memocc(i_stat,i_all,'hess','hessian_eigenvalues')

      i_all=-product(shape(hess_proj))*kind(hess_proj)
      deallocate(hess_proj,stat=i_stat)
      call memocc(i_stat,i_all,'hess_proj','hessian_eigenvalues')

   end subroutine hessian_eigenvalues

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: hessian_SP_eigenvalues
   !> @brief Wrapper routine for the determination of the Hessian and its eigenvalues
   !> for the saddle point
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine hessian_SP_eigenvalues(nHam,Natom,do_dm,max_no_neigh,max_no_dmneigh,  &
      eig_0,is_afm,simid,state,NNN,aham,taniso,nlistsize,dmlistsize, nlist,dmlist,  &
      bass,basis,coord,eaniso,kaniso,ncoup,emomM,dm_vect,bou_con,nvec_curr,         &
      info_curr,vol_curr,eig_curr,vel_par,beff,emom,mmom)

      implicit none

      integer, intent(in) :: nham            !< number of atoms in hamiltonian
      integer, intent(in) :: Natom           !< number of atoms in system
      integer, intent(in) :: do_dm           !< add dzyaloshinskii-moriya (dm) term to hamiltonian (0/1)
      integer, intent(in) :: max_no_neigh    !< calculated maximum of neighbours for exchange
      integer, intent(in) :: max_no_dmneigh  !< calculated number of neighbours with dm interactions
      real(dblprec), intent(in) :: eig_0     !< Threshold value for "zero" eigenvalue
      character(len=8), intent(in)  :: simid    !< name of simulation
      character(len=10), intent(in) :: state    !< which is the state that is being calculated (initial/final)
      character(len=1), intent(in)  :: is_afm   !< Flag to indicate if the state is Antiferromagnetic
      integer, dimension(3), intent(in) :: NNN
      integer, dimension(Natom), intent(in) :: aham      !< hamiltonian look-up table
      integer, dimension(Natom), intent(in) :: taniso    !< type of anisotropy (0-2)
      integer, dimension(nham), intent(in) :: nlistsize  !< size of neighbour list for heisenberg exchange couplings
      integer, dimension(nham), intent(in) :: dmlistsize !< size of neighbour list for dm
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist    !< neighbour list for heisenberg exchange couplings
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist !< list of neighbours for dm
      real(dblprec), dimension(Natom), intent(in) :: mmom   !< Magnitude of magnetic moments
      real(dblprec), dimension(3,3), intent(in) :: bass
      real(dblprec), dimension(3,Natom), intent(in) :: beff    !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom), intent(in) :: emom    !< current magnetic moment vector
      real(dblprec), dimension(3,Natom), intent(in) :: coord   !< Coordinates of atoms
      real(dblprec), dimension(3,Natom), intent(in) :: emomM   !< current magnetic moment vector
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso  !< unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso  !< anisotropy constant
      real(dblprec), dimension(max_no_neigh,nham), intent(in) :: ncoup        !< heisenberg exchange couplings
      real(dblprec), dimension(3,max_no_dmneigh,nham), intent(in) :: dm_vect  !< dzyaloshinskii-moriya exchange vector
      character(len=1), dimension(3), intent(in) :: bou_con !< Boundary conditions
      ! .. Output variables
      integer, intent(out) :: nvec_curr   !< Number of eigenvectors of the current image
      integer, intent(out) :: info_curr   !< Information flag of the current image
      !.. In/out variables
      real(dblprec), intent(inout) :: vol_curr  !< Volume of the zero modes in the current image
      real(dblprec), dimension(2*Natom), intent(inout) :: eig_curr   !< Eigenvalues of the current image
      real(dblprec), dimension(3*Natom,2*Natom), intent(inout) :: basis
      ! .. In/out allocatable variables
      real(dblprec), dimension(:), allocatable, intent(inout) :: vel_par !< velocity Parameters

      ! .. Local variables
      integer :: i_stat, i_all
      integer :: ii, jj, kk
      integer :: info
      real(dblprec) :: fc, fcinv, tmp
      character(len=256) :: msg
      character(len=35) :: num
      logical :: do_hess_calc

      do_hess_calc=.true.

      fc=mry/mub
      fcinv=mub/mry

      allocate(hess(3*Natom,3*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(hess))*kind(hess),'hess','hessian_SP_eigenvalues')
      hess=0.0_dblprec

      allocate(hess_proj(2*Natom,2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(hess_proj))*kind(hess_proj),'hess_proj','hessian_SP_eigenvalues')
      hess_proj=0.0_dblprec

      write(*,'(2x,a,x,a)') state, 'state'
      write(*,'(2x,a)') 'Evaluation of the Hessian   '
      call calc_hess_lambda(Natom, nHam,emomM,beff,max_no_neigh, nlistsize,nlist,   &
         ncoup,aham,do_dm, max_no_dmneigh, dmlistsize,dmlist, dm_vect,taniso,eaniso,&
         kaniso, hess)
      write(*,'(2x,a)') 'Done'

      write(*,'(2x,a)',advance='no') 'Projection of the Hessian on the tangent space '
      call calc_basis(Natom,emom,basis)

      call project_hess(Natom,basis,hess,hess_proj)

      write (*,'(x,a)') "done"
      write (*,'(2x,a)',advance='no') "Evaluation of the eigenvalues "
      call calc_hess_eig(Natom,hess_proj,eig_curr,info)
      write (*,'(2x,a)') "Done"
      if (info.ne.0) then
         write(*,'(2x,a)') 'WARNING: Failed to calculate eigenvalues!'
         info_curr = info
      end if

      open(ofileno, file='eig_'//trim(adjustl(state))//'.'//trim(adjustl(simid))//'.out',&
      access = 'sequential',action = 'write', status = 'replace')
      write(ofileno,'(a8,2x,a16)') "# Ent.","Eigenvalue"
      do ii=1,2*Natom
         write(ofileno,'(i8,2x,es16.8E3)',advance = 'yes') ii,eig_curr(ii)*fcinv
      end do
      close(ofileno)

      !------------------------------------------------------------------------------
      ! Check if the selected image is a saddle point, if it is then the minimum
      ! value is smaller than "zero" and the second smallest eigenvalue must be
      ! larger than "zero"
      !------------------------------------------------------------------------------
      ! Check if the minimum value of the eigenvalues is larger than the "zero" value
      if ((eig_curr(1)>0.0_dblprec).or.(abs(eig_curr(1)*fcinv).le.eig_0)) then
         write(num,'(es16.8E3)') eig_curr(1)*fcinv
         if (eig_curr(1)>0.0_dblprec) then
            msg = 'WARNING: Saddle point has not been found, because '//trim(adjustl(num))//'>0!'
         else
            msg = 'WARNING: Saddle point has not been found, because '//trim(adjustl(num))//'=0!'
         end if
         write(*,'(a)') trim(adjustl(msg))
         info_curr = 1
         nvec_curr = 0
         do_hess_calc=.false.
      end if

      ! Check if the second smallest value of the eigenvalues is smaller than the "zero" value
      if ((eig_curr(2)<0.0_dblprec).and.(abs(eig_curr(2)*fcinv)>eig_0)) then
         write(num,'(es16.8E3)') eig_curr(2)*fcinv
         msg = 'WARNING: Second eigenvalue is negative at the Saddle point, '//trim(adjustl(num))//'<0!'
         write(*,'(a)') trim(adjustl(msg))
         info_curr = 1
         nvec_curr = 0
         do_hess_calc=.false.
      end if

      if (do_hess_calc) then

         write (*,'(2x,a)') 'Checking for zero modes at the Saddle point '
         nvec_curr = 1
         do ii=2,2*Natom
            if (abs(eig_curr(ii)*fcinv).le.eig_0) then
               nvec_curr = nvec_curr+1
            end if
         end do

         if (nvec_curr==1) then
            write (*,'(2x,a)') "No zero modes have been found."
         else
            write(num,'(I8)') nvec_curr-1
            msg = trim(adjustl(num))//' zero modes have been found.'
            write (*,'(2x,a)') trim(adjustl(msg))
            if (is_afm=='Y') then
               call calc_zero_vol_afm(Natom,NNN,emomM,coord,bass,bou_con,vol_curr)
            else
               call calc_zero_vol(Natom,NNN,emomM,coord,bass,bou_con,vol_curr)
            end if
         end if

         write (*,'(2x,a)',advance='no') "Calculation of the eigenvectors "
         allocate(vec_sp(2*Natom,2*Natom),stat=i_stat)
         call memocc(i_stat,product(shape(vec_sp))*kind(vec_sp),'vec_sp','hessian_SP_eigenvalues')
         vec_sp=0.0_dblprec

         call calc_hess_eig_VV(Natom,hess_proj,vec_sp,info)
         write (*,'(x,a)') "Done"
         if (info.ne.0) then

            write(*,'(2x,a)') 'WARNING: Failed to calculate eigenvectors!'
            info_curr = info
         else
            do kk=1,3
               tmp = 0.0_dblprec
               !$omp parallel do default(shared), private(ii,jj), reduction(+:tmp), collapse(2)
               do ii=1,2*Natom
                  do jj=1,2*Natom
                     tmp = tmp + vec_sp(ii,kk)*hess_proj(ii,jj)*vec_sp(jj,kk)
                  end do
               end do
               !$omp end parallel do
               write(*,'(3x,a,x,es16.8,x,es16.8)')'compare eig:',eig_curr(kk)*fcinv,tmp*fcinv
            end do
            write (*,'(2x,a)',advance='no') "Calculation of the velocity parameters "

            allocate(vel_par(2*Natom),stat=i_stat)
            call memocc(i_stat,product(shape(vel_par))*kind(vel_par),'vel_par','hessian_SP_eigenvalues')
            vel_par=0.0_dblprec

            call calc_vel_par(Natom,emomM,hess,basis,vec_sp,vel_par,mmom)
            open(ofileno, file='vel.'//trim(adjustl(simid))//'.out', access = 'sequential',action = 'write', status = 'replace')
            do ii=1,2*Natom
               vel_par(ii) = vel_par(ii)*fcinv
               write(ofileno,'(es16.8E3)',advance = 'yes') vel_par(ii)
            end do
            close(ofileno)
            write (*,'(x,a)') "done"
         end if

         if (nvec_curr>0) then
            i_all=-product(shape(vec_sp))*kind(vec_sp)
            deallocate(vec_sp,stat=i_stat)
            call memocc(i_stat,i_all,'vec_sp','hessian_SP_eigenvalues')
         endif

      endif

      write (*,'(2x,a)') "-----------------------------------------"

      eig_curr(:) = eig_curr(:)*fcinv

      i_all=-product(shape(hess))*kind(hess)
      deallocate(hess,stat=i_stat)
      call memocc(i_stat,i_all,'hess','hessian_SP_eigenvalues')

      i_all=-product(shape(hess_proj))*kind(hess_proj)
      deallocate(hess_proj,stat=i_stat)
      call memocc(i_stat,i_all,'hess_proj','hessian_SP_eigenvalues')

   end subroutine hessian_SP_eigenvalues

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_hess
   !> @brief Calculates a matrix of second derivatives of the energy with respect to
   !> the Cartesian coordinates of normalized moments
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_hess(Natom,nHam,emomM,max_no_neigh,nlistsize, nlist,ncoup,aham,  &
      do_dm,max_no_dmneigh,dmlistsize,dmlist,dm_vect,taniso,eaniso,kaniso,hess)
      !
      use Constants
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: do_dm !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, intent(in) :: max_no_neigh   !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      integer, dimension(Natom),intent(in) :: aham !< Hamiltonian look-up table
      integer, dimension(Natom),intent(in) :: taniso      !< Type of anisotropy (0-2)
      integer, dimension(nHam), intent(in) :: nlistsize  !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist    !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist !< List of neighbours for DM
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(3,Natom), intent(in) :: emomM    !< Current magnetic moment vector
      real(dblprec), dimension(max_no_neigh,nHam), intent(in) :: ncoup  !< Heisenberg exchange couplings
      real(dblprec), dimension(3,max_no_dmneigh,nHam), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      !.. Output variables
      real(dblprec), dimension(3*Natom,3*Natom), intent(out) :: hess !< Hessian of the system

      !.. Local scalars
      integer :: ii, jj, ih
      real(dblprec) :: fcinv,tmpi,tmpj

      ! Initialize energy variables
      hess = 0.0_dblprec

      ! Factor for energy scale
      fcinv = mub/mry

      ! Loop over atoms
      !$omp parallel do default(shared), private(ii,jj,ih,tmpi,tmpj)
      do ii=1, Natom
         ih=aham(ii)
         tmpi = 0.0_dblprec
         tmpi = sum(emomM(:,ii)*emomM(:,ii))
         tmpi = sqrt(tmpi)
         !---------------------------------------------------------------------------
         ! Exchange interactions
         !---------------------------------------------------------------------------
         do jj=1,nlistsize(ih)
            tmpj = 0.0_dblprec
            tmpj = sum(emomM(:,nlist(jj,ii))*emomM(:,nlist(jj,ii)))
            tmpj = sqrt(tmpj)
            hess(3*ii-2,3*nlist(jj,ii)-2) = hess(3*ii-2,3*nlist(jj,ii)-2) - tmpi*tmpj*ncoup(jj,ih)
            hess(3*ii-1,3*nlist(jj,ii)-1) = hess(3*ii-1,3*nlist(jj,ii)-1) - tmpi*tmpj*ncoup(jj,ih)
            hess(3*ii,3*nlist(jj,ii)) = hess(3*ii,3*nlist(jj,ii)) - tmpi*tmpj*ncoup(jj,ih)
         end do
         if(do_dm==1) then
            !------------------------------------------------------------------------
            ! Dzyaloshinskii-Moriya interactions
            !------------------------------------------------------------------------
            do jj=1,dmlistsize(ih)
               tmpj = 0.0_dblprec
               tmpj = sum(emomM(:,dmlist(jj,ii))*emomM(:,dmlist(jj,ii)))
               tmpj = sqrt(tmpj)
               hess(3*dmlist(jj,ii),3*ii-1) = hess(3*dmlist(jj,ii),3*ii-1) + tmpi*tmpj*dm_vect(1,jj,ih)
               hess(3*dmlist(jj,ii)-1,3*ii) = hess(3*dmlist(jj,ii)-1,3*ii) - tmpi*tmpj*dm_vect(1,jj,ih)
               hess(3*dmlist(jj,ii)-2,3*ii) = hess(3*dmlist(jj,ii)-2,3*ii) + tmpi*tmpj*dm_vect(2,jj,ih)
               hess(3*dmlist(jj,ii),3*ii-2) = hess(3*dmlist(jj,ii),3*ii-2) - tmpi*tmpj*dm_vect(2,jj,ih)
               hess(3*dmlist(jj,ii)-1,3*ii-2) = hess(3*dmlist(jj,ii)-1,3*ii-2) + tmpi*tmpj*dm_vect(3,jj,ih)
               hess(3*dmlist(jj,ii)-2,3*ii-1) = hess(3*dmlist(jj,ii)-2,3*ii-1) - tmpi*tmpj*dm_vect(3,jj,ih)
            end do
         end if
         !---------------------------------------------------------------------------
         ! Uniaxial anisotropy
         !---------------------------------------------------------------------------
         if (taniso(ii)==1) then
            hess(3*ii-2,3*ii-2) = 2.0_dblprec*kaniso(1,ii)*eaniso(1,ii)*eaniso(1,ii)*tmpi*tmpi
            hess(3*ii-1,3*ii-1) = 2.0_dblprec*kaniso(1,ii)*eaniso(2,ii)*eaniso(2,ii)*tmpi*tmpi
            hess(3*ii,3*ii) = 2.0_dblprec*kaniso(1,ii)*eaniso(3,ii)*eaniso(3,ii)*tmpi*tmpi
         endif
      end do
      !$omp end parallel do

   end subroutine calc_hess

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_hess_lambda
   !> @brief Calculates a matrix of second derivatives of the Lagrange function
   !> with respect to Cartesian coordinates of normalized moments
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_hess_lambda(Natom,nHam, emomM,beff,max_no_neigh,nlistsize,nlist, &
      ncoup,aham,do_dm,max_no_dmneigh,dmlistsize,dmlist,dm_vect,taniso,eaniso,      &
      kaniso,hess_lambda)
      !
      use Constants
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: nHam    !< Number of atoms in Hamiltonian
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: Natom   !< Number of atoms in system
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: max_no_dmneigh    !< Calculated number of neighbours with DM interactions
      integer, dimension(Natom),intent(in) :: aham       !< Hamiltonian look-up table
      integer, dimension(Natom),intent(in) :: taniso     !< Type of anisotropy (0-2)
      integer, dimension(nHam), intent(in) :: nlistsize  !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist    !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist !< List of neighbours for DM
      real(dblprec), dimension(3,Natom), intent(in) :: beff   !< Effective magnetic field
      real(dblprec), dimension(3,Natom), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(max_no_neigh,nHam), intent(in) :: ncoup  !< Heisenberg exchange couplings
      real(dblprec), dimension(3,max_no_dmneigh,nHam), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      ! .. Output variables
      real(dblprec), dimension(3*Natom,3*Natom), intent(out) :: hess_lambda !< Hessian of the system

      !.. Local variables
      integer :: ii,jj
      integer :: i_stat,i_all
      real(dblprec) :: hes_scale_fac
      real(dblprec), dimension(:), allocatable :: lambda
      real(dblprec), dimension(:,:), allocatable :: hess

      allocate(lambda(Natom),stat=i_stat)
      call memocc(i_stat,product(shape(lambda))*kind(lambda),'lambda','calc_hess_lambda')
      lambda=0.0_dblprec
      allocate(hess(3*Natom,3*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(hess))*kind(hess),'hess','calc_hess_lambda')
      hess=0.0_dblprec

      call calc_hess(Natom,nHam,emomM,max_no_neigh,nlistsize,nlist,ncoup,aham,do_dm,&
         max_no_dmneigh,dmlistsize,dmlist,dm_vect,taniso, eaniso, kaniso, hess)

      call calc_lambda(Natom,emomM,beff,lambda)

      hes_scale_fac=ry_ev*mub/mry

      do ii=1,3
         write(*,'(3x,a,x,es16.8,es16.8,es16.8)') 'Hess:', hess(ii,1)*hes_scale_fac,&
            hess(ii,2)*hes_scale_fac,hess(ii,3)*hes_scale_fac
      end do

      call dcopy(9*Natom*Natom,hess,1,hess_lambda,1)

      !$omp parallel do default(shared), private(ii,jj), collapse(2)
      do ii=1,Natom
         do jj=1,3
            hess_lambda(3*ii-3+jj,3*ii-3+jj) = hess_lambda(3*ii-3+jj,3*ii-3+jj) + lambda(ii)
         end do
      end do
      !$omp end parallel do
      i_all=-product(shape(lambda))*kind(lambda)
      deallocate(lambda,stat=i_stat)
      call memocc(i_stat,i_all,'lambda','calc_hess_lambda')

      i_all=-product(shape(hess))*kind(hess)
      deallocate(hess,stat=i_stat)
      call memocc(i_stat,i_all,'hess','calc_hess_lambda')

   end subroutine calc_hess_lambda

   !---------------------------------------------------------------------------------
   ! subroutine calc_lambda
   !> @brief calculates lagrange multipliers
   !> @author pavel bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_lambda(Natom,emomM,beff,lambda)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3,Natom), intent(in) :: beff  !< Effective magnetic field
      real(dblprec), dimension(3,Natom), intent(in) :: emomM !< Magnetic moment vector
      real(dblprec), dimension(Natom), intent(out) :: lambda !< Lagrange multipliers
      !.. Local variables
      integer :: ii
      real(dblprec) :: tmp1,tmp2,tmp3

      !$omp parallel do default(shared), private(ii,tmp1,tmp2,tmp3)
      do ii=1,Natom
         tmp1 = 0.0_dblprec
         tmp2 = 0.0_dblprec
         tmp3 = 0.0_dblprec
         !tmp1 = sum(beff(:,ii))
         !tmp2 = sum(emomM(:,ii))
         !tmp3 = norm2(emomM(:,ii))
         !lambda(ii) = tmp1*tmp3*tmp3/tmp2
         lambda(ii) = emomM(1,ii)*beff(1,ii)+emomM(2,ii)*beff(2,ii)+emomM(3,ii)*beff(3,ii)
      end do
      !$omp end parallel do
   end subroutine calc_lambda

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: project_hess
   !> @brief Calculates the Projection of the Hessian
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine project_hess(Natom,basis,hess_lambda,hess_proj)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3*Natom,2*Natom), intent(in) :: basis    !< basis
      real(dblprec), dimension(3*Natom,3*Natom), intent(in) :: hess_lambda !< Hessian of the Lagrange function
      ! .. Output variables
      real(dblprec), dimension(2*Natom,2*Natom), intent(out) :: hess_proj !< projected Hessian
      ! .. Local variables
      integer :: ii,jj
      integer :: i_stat,i_all
      ! .. Local allocatable variables
      real(dblprec), dimension(:,:), allocatable :: aux
      real(dblprec), dimension(:,:), allocatable :: aux1
      real(dblprec), dimension(:,:), allocatable :: aux2

      allocate(aux1(2,3),stat=i_stat)
      call memocc(i_stat,product(shape(aux1))*kind(aux1),'aux1','project_hess')
      aux1=0.0_dblprec
      allocate(aux2(2,2),stat=i_stat)
      call memocc(i_stat,product(shape(aux2))*kind(aux2),'aux2','project_hess')
      aux2=0.0_dblprec
      allocate(aux(2*Natom,2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(aux))*kind(aux),'aux','project_hess')
      aux1=0.0_dblprec
      !$omp parallel do default(shared), private(ii,jj,aux1,aux2), collapse(2)
      do ii = 1,Natom
         do jj = 1,Natom
            call dgemm('T','N',2,3,3,1.0_dblprec,basis(3*ii-2:3*ii,2*ii-1:2*ii),3,  &
               hess_lambda(3*ii-2:3*ii,3*jj-2:3*jj),3,0.0_dblprec,aux1,2)
            call dgemm('N','N',2,2,3,1.0_dblprec,aux1,2,                            &
               basis(3*jj-2:3*jj,2*jj-1:2*jj),3,0.0_dblprec,aux2,2)

            aux(2*ii-1,2*jj-1)   = aux2(1,1)
            aux(2*ii-1,2*jj)     = aux2(1,2)
            aux(2*ii,2*jj-1)     = aux2(2,1)
            aux(2*ii,2*jj)       = aux2(2,2)
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared), private(ii,jj), collapse(2)
      do ii=1,2*Natom
         do jj=1,2*Natom
            hess_proj(ii,jj) = 0.50_dblprec*(aux(ii,jj)+aux(jj,ii))
         end do
      end do
      !$omp end parallel do
      i_all=-product(shape(aux))*kind(aux)
      deallocate(aux,stat=i_stat)
      call memocc(i_stat,i_all,'aux','project_hess')
      i_all=-product(shape(aux1))*kind(aux1)
      deallocate(aux1,stat=i_stat)
      call memocc(i_stat,i_all,'aux1','project_hess')
      i_all=-product(shape(aux2))*kind(aux2)
      deallocate(aux2,stat=i_stat)
      call memocc(i_stat,i_all,'aux2','project_hess')

   end subroutine project_hess

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_basis
   !> @brief Calculate basis for the tangent space
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_basis(Natom,emom,basis)

      use math_functions, only: f_normalize_vec, f_cross_product
      use RandomNumbers, only : rng_uniform,rng_init

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3,Natom), intent(in) :: emom    !< Magnetic moment vector
      ! .. Output variables
      real(dblprec), dimension(3*Natom,2*Natom), intent(out) :: basis
      ! .. Local variables
      integer :: ii,jj
      integer :: i_stat,i_all
      real(dblprec) :: tmp
      real(dblprec) :: eps
      ! .. Local allocatable variables
      real(dblprec), dimension(:), allocatable :: xx
      real(dblprec), dimension(:), allocatable :: yy
      real(dblprec), dimension(:), allocatable :: rand_numbers

      tmp=0.0_dblprec
      eps=epsilon(tmp)

      allocate(xx(3),stat=i_stat)
      call memocc(i_stat,product(shape(xx))*kind(xx),'xx','calc_basis')
      xx=0.0_dblprec
      allocate(yy(3),stat=i_stat)
      call memocc(i_stat,product(shape(yy))*kind(yy),'yy','calc_basis')
      yy=0.0_dblprec
      allocate(rand_numbers(3),stat=i_stat)
      call memocc(i_stat,product(shape(rand_numbers))*kind(rand_numbers),'rand_numbers','calc_basis')
      rand_numbers=0.0_dblprec

      basis = 0.0_dblprec

      do ii=1,Natom
         tmp = 0.0_dblprec
         call rng_uniform(rand_numbers,3)
         do jj=1,3
            xx(jj) = 2.0_dblprec*rand_numbers(jj)-1.0_dblprec
            tmp = tmp + xx(jj)*emom(jj,ii)
         end do
         do jj=1,3
            xx(jj) = xx(jj) - tmp*emom(jj,ii)
         end do
         tmp = norm2(xx)
         do while (tmp<eps)
            call rng_uniform(rand_numbers,3)
            tmp = 0.0_dblprec
            do jj=1,3
               xx(jj) = 2.0_dblprec*rand_numbers(jj)-1.0_dblprec
               tmp = tmp + xx(jj)*emom(jj,ii)
            end do
            do jj=1,3
               xx(jj) = xx(jj) - tmp*emom(jj,ii)
            end do
            tmp = norm2(xx)
         end do

         yy = f_cross_product(emom(:,ii),xx)

         do jj=1,3
            basis(3*ii-3+jj,2*ii-1) = xx(jj)
            basis(3*ii-3+jj,2*ii)   = yy(jj)
         end do
      end do

      !$omp parallel do default(shared), private(ii)
      do ii=1,2*Natom
         basis(:,ii)=f_normalize_vec(basis(:,ii),3*Natom)
      end do
      !$omp end parallel do

      i_all=-product(shape(xx))*kind(xx)
      deallocate(xx,stat=i_stat)
      call memocc(i_stat,i_all,'xx','calc_hess')
      i_all=-product(shape(yy))*kind(yy)
      deallocate(yy,stat=i_stat)
      call memocc(i_stat,i_all,'yy','calc_hess')
      i_all=-product(shape(rand_numbers))*kind(rand_numbers)
      deallocate(rand_numbers,stat=i_stat)
      call memocc(i_stat,i_all,'rand_numbers','calc_hess')

   end subroutine calc_basis

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_hess_eig
   !> @brief Calculates eigenvalues of the projected Hessian
   !> @author Pavel Bessarab
   !> @note Changed eigensolver to dsyevr by Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine calc_hess_eig(Natom,hess,eig,info)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(2*Natom,2*Natom), intent(in) :: hess !< projected Hessian
      ! .. Output variables
      integer, intent(out) :: info
      real(dblprec), dimension(2*Natom), intent(out) :: eig !< eigenvalues
      ! .. Local variables

      call eigen_driver(hess,2*Natom,eig,info)

   end subroutine calc_hess_eig

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_hess_eig_VV
   !> @brief Calculates first eigenvalues and eigenvectors of the projected Hessian
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_hess_eig_VV(Natom,hess,vec,info)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(2*Natom,2*Natom), intent(in) :: hess!< projected Hessian

      ! .. Output variables
      integer, intent(out) :: info
      real(dblprec), dimension(2*Natom,2*Natom), intent(out) :: vec !< eigenvalues

      ! .. Local variables
      integer :: i_stat,i_all

      ! .. Local allocatable variables
      real(dblprec), dimension(:), allocatable :: eig

      allocate(eig(2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(eig))*kind(eig),'eig','calc_hess_eig_VV')
      
      call eigen_driver(hess,2*Natom,eig,info,vec)
      
      i_all=-product(shape(eig))*kind(eig)
      deallocate(eig,stat=i_stat)
      call memocc(i_stat,i_all,'eig','calc_hess_eig_VV')

   end subroutine calc_hess_eig_VV

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_hess_V
   !> @brief Calculates first eigenvectors of the projected Hessian
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_hess_V(Natom,m,hess,eig,basis,vec,info,vec2)

      implicit none

      integer, intent(in) :: m
      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(2*Natom), intent(in) :: eig !< eigenvalues
      real(dblprec), dimension(2*Natom,2*Natom), intent(in) :: hess !< projected Hessian
      real(dblprec), dimension(3*Natom,2*Natom), intent(in) :: basis !< basis
      ! .. Output variables
      integer, intent(out) :: info
      real(dblprec), dimension(3*Natom,m), intent(out) :: vec !< eigenvectors in original basis
      real(dblprec), optional, dimension(2*Natom,m), intent(out) :: vec2
      ! .. Local variables
      integer :: lwork,info_i,ii,jj
      integer :: i_stat, i_all
      integer(8) :: lwmax
      real(dblprec), dimension(1) :: ww

      ! .. Local allocatable variables
      integer, dimension(:), allocatable :: iwork
      integer, dimension(:), allocatable :: iblock
      integer, dimension(:), allocatable :: isplit
      integer, dimension(:), allocatable :: ifailv
      real(dblprec), dimension(:), allocatable :: d
      real(dblprec), dimension(:), allocatable :: e
      real(dblprec), dimension(:), allocatable :: tau
      real(dblprec), dimension(:), allocatable :: work
      real(dblprec), dimension(:,:), allocatable :: hwork
      real(dblprec), dimension(:,:), allocatable :: vecaux


      lwmax = 5*Natom*Natom

      allocate(iwork(6*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(iwork))*kind(iwork),'iwork','calc_hess_V')
      iwork=0
      allocate(iblock(2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(iblock))*kind(iblock),'iblock','calc_hess_V')
      iblock=0
      allocate(isplit(2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(isplit))*kind(isplit),'isplit','calc_hess_V')
      isplit=0
      allocate(ifailv(m),stat=i_stat)
      call memocc(i_stat,product(shape(ifailv))*kind(ifailv),'ifailv','calc_hess_V')
      ifailv=0
      allocate(d(2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(d))*kind(d),'d','calc_hess_V')
      d=0.0_dblprec
      allocate(e(2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(e))*kind(e),'e','calc_hess_V')
      e=0.0_dblprec
      allocate(tau(2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(tau))*kind(tau),'tau','calc_hess_V')
      tau=0.0_dblprec
      allocate(hwork(2*Natom,2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(hwork))*kind(hwork),'hwork','calc_hess_V')
      hwork=0.0_dblprec
      allocate(vecaux(2*Natom,m),stat=i_stat)
      call memocc(i_stat,product(shape(vecaux))*kind(vecaux),'vecaux','calc_hess_V')
      vecaux=0.0_dblprec

      !!! Possible change to eigen_driver
      !!! e=eig
      !!! call eigen_driver(hess,2*Natom,e,info,vec,m)

      !$omp parallel do default(shared), private(ii,jj)
      do ii=1,2*Natom
         do jj=1,2*Natom
            hwork(ii,jj) = hess(ii,jj)
         end do
         iblock(ii) = 1
         isplit(ii) = 2*Natom
      end do
      !$omp end parallel do

      lwork = -1

      call dsytrd('U',2*Natom,hwork,2*Natom,d,e,tau,ww,lwork,info_i)

      lwork = int(ww(1))
      allocate(work(lwork),stat=i_stat)
      call memocc(i_stat,product(shape(work))*kind(work),'work','calc_hess_V')

      call dsytrd('U',2*Natom,hwork,2*Natom,d,e,tau,work,lwork,info_i)

      info = info_i

      call dstein(2*Natom,d,e,m,eig,iblock,isplit,vecaux,2*Natom,work,iwork,ifailv, &
         info_i)
      if (info_i.ne.0) then
         info = info_i
      end if

      !lwork = -1
      !call dormtr('L','U','N',2*Natom,m,hwork,2*Natom,tau,vecaux,2*Natom,work,lwork,&
      !   info_i)
      !lwork = lwmax !min(lwmax,int(work(1)))
      call dormtr('L','U','N',2*Natom,m,hwork,2*Natom,tau,vecaux,2*Natom,work,lwork,&
         info_i)
      if (info_i.ne.0) then
         info = info_i
      end if

      call dgemm('N','N',3*Natom,m,2*Natom,1.0_dblprec,basis,3*Natom,vecaux,2*Natom,&
         0.0_dblprec,vec,3*Natom)

      if (present(vec2)) then
         vec2 = vecaux
      end if

      ! Deallocation of arrays
      i_all=-product(shape(d))*kind(d)
      deallocate(d,stat=i_stat)
      call memocc(i_stat,i_all,'d','calc_hess_V')
      i_all=-product(shape(e))*kind(e)
      deallocate(e,stat=i_stat)
      call memocc(i_stat,i_all,'e','calc_hess_V')
      i_all=-product(shape(tau))*kind(tau)
      deallocate(tau,stat=i_stat)
      call memocc(i_stat,i_all,'tau','calc_hess_V')
      i_all=-product(shape(ifailv))*kind(ifailv)
      deallocate(ifailv,stat=i_stat)
      call memocc(i_stat,i_all,'ifailv','calc_hess_V')
      i_all=-product(shape(work))*kind(work)
      deallocate(work,stat=i_stat)
      call memocc(i_stat,i_all,'work','calc_hess_V')
      i_all=-product(shape(hwork))*kind(hwork)
      deallocate(hwork,stat=i_stat)
      call memocc(i_stat,i_all,'hwork','calc_hess_V')
       i_all=-product(shape(vecaux))*kind(vecaux)
       deallocate(vecaux,stat=i_stat)
       call memocc(i_stat,i_all,'vecaux','calc_hess_V')
      i_all=-product(shape(iwork))*kind(iwork)
      deallocate(iwork,stat=i_stat)
      call memocc(i_stat,i_all,'iwork','calc_hess_V')
      i_all=-product(shape(iblock))*kind(iblock)
      deallocate(iblock,stat=i_stat)
      call memocc(i_stat,i_all,'iblock','calc_hess_V')
      i_all=-product(shape(isplit))*kind(isplit)
      deallocate(isplit,stat=i_stat)
      call memocc(i_stat,i_all,'isplit','calc_hess_V')

   end subroutine calc_hess_V

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_vel_mat
   !> @brief Calculates velocity matrix
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_vel_mat1(Natom,emomM,hl,vel,mmom)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(Natom), intent(in) :: mmom !< MMagnitude of the magnetic moments
      real(dblprec), dimension(3,Natom), intent(in) :: emomM !< Magnetic moment vector
      real(dblprec), dimension(3*Natom,3*Natom), intent(in) :: hl !< Hessian
      real(dblprec), dimension(3*Natom,3*Natom), intent(out) :: vel !< velocity
      ! .. Local variables
      integer :: ii,jj
      real(dblprec) :: mmom2,beffx,beffy,beffz,mmomi,mmomj

      vel = 0.0_dblprec

      !$omp parallel do default(shared), private(ii,jj,mmom2,mmomi,mmomj,beffx,beffy,beffz)
      do ii=1,Natom
         mmomi = mmom(ii)
         mmom2 = mmomi**2
         beffx = 0.0_dblprec
         beffy = 0.0_dblprec
         beffz = 0.0_dblprec
         do jj=1,Natom
            mmomj = mmom(jj)
            vel(3*ii-2,3*jj-2) = 0.50_dblprec*(emomM(2,ii)*(hl(3*jj-2,3*ii)+hl(3*ii,3*jj-2))-emomM(3,ii)*(hl(3*jj-2,3*ii-1)+hl(3*ii-1,3*jj-2)))/mmom2
            vel(3*ii-2,3*jj-1) = 0.50_dblprec*(emomM(2,ii)*(hl(3*jj-1,3*ii)+hl(3*ii,3*jj-1))-emomM(3,ii)*(hl(3*jj-1,3*ii-1)+hl(3*ii-1,3*jj-1)))/mmom2
            vel(3*ii-2,3*jj) = 0.50_dblprec*(emomM(2,ii)*(hl(3*jj,3*ii)+hl(3*ii,3*jj))-emomM(3,ii)*(hl(3*jj,3*ii-1)+hl(3*ii-1,3*jj)))/mmom2

            vel(3*ii-1,3*jj-2) = 0.50_dblprec*(emomM(3,ii)*(hl(3*jj-2,3*ii-2)+hl(3*ii,3*jj-2))-emomM(1,ii)*(hl(3*jj-2,3*ii)+hl(3*ii,3*jj-2)))/mmom2
            vel(3*ii-1,3*jj-1) = 0.50_dblprec*(emomM(3,ii)*(hl(3*jj-1,3*ii-2)+hl(3*ii-2,3*jj-1))-emomM(1,ii)*(hl(3*jj-1,3*ii)+hl(3*ii,3*jj-1)))/mmom2
            vel(3*ii-1,3*jj) = 0.50_dblprec*(emomM(3,ii)*(hl(3*jj,3*ii-2)+hl(3*ii-2,3*jj))-emomM(1,ii)*(hl(3*jj,3*ii)+hl(3*ii,3*jj)))/mmom2

            vel(3*ii,3*jj-2) = 0.50_dblprec*(emomM(1,ii)*(hl(3*jj-2,3*ii-1)+hl(3*ii-1,3*jj-2))-emomM(2,ii)*(hl(3*jj-2,3*ii-2)+hl(3*ii-2,3*jj-2)))/mmom2
            vel(3*ii,3*jj-1) = 0.50_dblprec*(emomM(1,ii)*(hl(3*jj-1,3*ii-1)+hl(3*ii-1,3*jj-1))-emomM(2,ii)*(hl(3*jj-1,3*ii-2)+hl(3*ii-2,3*jj-1)))/mmom2
            vel(3*ii,3*jj) = 0.50_dblprec*(emomM(1,ii)*(hl(3*jj,3*ii-1)+hl(3*ii-1,3*jj))-emomM(2,ii)*(hl(3*jj,3*ii-2)+hl(3*ii-2,3*jj)))/mmom2

            beffx = beffx - (emomM(1,jj)*(hl(3*jj-2,3*ii-2)+hl(3*ii-2,3*jj-2))&
                  +emomM(2,jj)*(hl(3*jj-1,3*ii-2)+hl(3*ii-2,3*jj-1))&
                  +emomM(3,jj)*(hl(3*jj,3*ii-2)+hl(3*ii-2,3*jj)))/mmomj
            beffy = beffy - (emomM(1,jj)*(hl(3*jj-2,3*ii-1)+hl(3*ii-1,3*jj-2))&
                  +emomM(2,jj)*(hl(3*jj-1,3*ii-1)+hl(3*ii-1,3*jj-1))&
                  +emomM(3,jj)*(hl(3*jj,3*ii-1)+hl(3*ii-1,3*jj)))/mmomj
            beffz = beffz - (emomM(1,jj)*(hl(3*jj-2,3*ii)+hl(3*ii,3*jj-2))&
                  +emomM(2,jj)*(hl(3*jj-1,3*ii)+hl(3*ii,3*jj-1))&
                  +emomM(3,jj)*(hl(3*jj,3*ii)+hl(3*ii,3*jj)))/mmomj
         end do
         beffx = 0.50_dblprec*beffx/mmomi
         beffy = 0.50_dblprec*beffy/mmomi
         beffz = 0.50_dblprec*beffz/mmomi
         vel(3*ii-2,3*ii-1)   = vel(3*ii-2,3*ii-1)-beffz
         vel(3*ii-2,3*ii)     = vel(3*ii-2,3*ii)+beffy
         vel(3*ii-1,3*ii-2)   = vel(3*ii-1,3*ii-2)+beffz
         vel(3*ii-1,3*ii)     = vel(3*ii-1,3*ii)-beffx
         vel(3*ii,3*ii-2)     = vel(3*ii,3*ii-2)-beffy
         vel(3*ii,3*ii-1)     = vel(3*ii,3*ii-1)+beffx
      end do
      !$omp end parallel do

   end subroutine calc_vel_mat1

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: project_vel_mat
   !> @brief Projects velocity matrix on the tangent space
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine project_vel_mat(Natom,basis,vel,vel_proj)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3*Natom,3*Natom), intent(in) :: vel   !< velocity matrix
      real(dblprec), dimension(3*Natom,2*Natom), intent(in) :: basis !< basis
      ! .. Output variables
      real(dblprec), dimension(2*Natom,2*Natom), intent(out) :: vel_proj !< projected velocity
      ! .. Local variables
      integer :: ii,jj
      integer :: i_stat,i_all
      ! .. Local allocatable variables
      real(dblprec), dimension(:,:), allocatable :: aux1,aux2

      allocate(aux1(2,3),stat=i_stat)
      call memocc(i_stat,product(shape(aux1))*kind(aux1),'aux1','project_vel_mat')
      aux1=0.0_dblprec
      allocate(aux2(2,2),stat=i_stat)
      call memocc(i_stat,product(shape(aux2))*kind(aux2),'aux2','project_vel_mat')
      aux1=0.0_dblprec

      vel_proj = 0.0_dblprec

      !$omp parallel do default(shared), private(ii,jj,aux1,aux2), collapse(2)
      do ii = 1,Natom
         do jj = 1,Natom
            call dgemm('T','N',2,3,3,1.0_dblprec,basis(3*ii-2:3*ii,2*ii-1:2*ii),3,   &
               vel(3*ii-2:3*ii,3*jj-2:3*jj),3,0.0_dblprec,aux1,2)
            call dgemm('N','N',2,2,3,1.0_dblprec,aux1,2,                            &
               basis(3*jj-2:3*jj,2*jj-1:2*jj),3,0.0_dblprec,aux2,2)

            vel_proj(2*ii-1,2*jj-1) = aux2(1,1)
            vel_proj(2*ii-1,2*jj)   = aux2(1,2)
            vel_proj(2*ii,2*jj-1)   = aux2(2,1)
            vel_proj(2*ii,2*jj)     = aux2(2,2)
         end do
      end do
      !$omp end parallel do

      i_all=-product(shape(aux1))*kind(aux1)
      deallocate(aux1,stat=i_stat)
      call memocc(i_stat,i_all,'aux1','project_vel_mat')
      i_all=-product(shape(aux2))*kind(aux2)
      deallocate(aux2,stat=i_stat)
      call memocc(i_stat,i_all,'aux2','project_vel_mat')

   end subroutine project_vel_mat

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_vel
   !> @brief Calculates the velocity coefficients
   !> @details The velocity coefficients are the expansion coefficients of the linearlized
   !> LL equation, which are dubbed `a_i` in Bessarab et al. Sci. Rep. Vol. 8, 3433 (2018) 
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_vel_par(Natom,emomM,hl,basis,vec,vel_par,mmom)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(Natom), intent(in) :: mmom   !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom), intent(in) :: emomM         !< Magnetic moment vector
      real(dblprec), dimension(3*Natom,3*Natom), intent(in) :: hl    !< Hessian
      real(dblprec), dimension(3*Natom,2*Natom), intent(in) :: basis !< basis
      real(dblprec), dimension(2*Natom,2*Natom), intent(in) :: vec   !< eigenbasis
      ! .. Output variables
      real(dblprec), dimension(2*Natom), intent(inout) :: vel_par    !< Velocity coefficients
      ! .. Local allocatable variables
      real(dblprec), dimension(:,:), allocatable :: vel              !< velocity matrix
      real(dblprec), dimension(:,:), allocatable :: aux
      real(dblprec), dimension(:,:), allocatable :: vel_proj         !< projected velocity matrix
      ! .. Local variables
      integer :: ii,jj
      integer :: i_stat,i_all

      allocate(vel(3*Natom,3*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(vel))*kind(vel),'vel','calc_vel_par')
      vel=0.0_dblprec
      allocate(vel_proj(2*Natom,2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(vel_proj))*kind(vel_proj),'vel_proj','calc_vel_par')
      vel_proj=0.0_dblprec
      allocate(aux(2*Natom,2*Natom),stat=i_stat)
      call memocc(i_stat,product(shape(aux))*kind(aux),'aux','calc_vel_par')
      aux=0.0_dblprec

      call calc_vel_mat1(Natom,emomM,hl,vel,mmom)

      call project_vel_mat(Natom,basis,vel,vel_proj)

      call dgemm('N','N',2*Natom,2*Natom,2*Natom,1.0_dblprec,vel_proj,2*Natom,vec,  &
         2*Natom,0.0_dblprec,aux,2*Natom)

      vel_par=0.0_dblprec
      !$omp parallel do default(shared), private(ii,jj), collapse(2)
      do ii=1,2*Natom
         do jj=1,2*Natom
            vel_par(ii) = vel_par(ii) + vec(jj,1)*aux(jj,ii)
         end do
      end do
      !$omp end parallel do

      i_all=-product(shape(aux))*kind(aux)
      deallocate(aux,stat=i_stat)
      call memocc(i_stat,i_all,'aux','calc_vel_par')
      i_all=-product(shape(vel))*kind(vel)
      deallocate(vel,stat=i_stat)
      call memocc(i_stat,i_all,'vel','calc_vel_par')
      i_all=-product(shape(vel_proj))*kind(vel_proj)
      deallocate(vel_proj,stat=i_stat)
      call memocc(i_stat,i_all,'vel_proj','calc_vel_par')

   end subroutine calc_vel_par

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_pref
   !> @brief Calculates prefactor for the determination of the Arrhenius law 
   !> governing the transition between states
   !> @details The prefactor, \f$ \nu\f$, is calculated such that it can describe the transition
   !> between different states given by the Arrhenius law \f$ \tau = \nu \exp{\frac{-\Delta E}{K_BT}} \f$
   !> where \f$ \Delta E\f$ is the energy barrier between the states and \f$\tau\f$ is the transition rate.
   !> 
   !> The prefactor can be obtained by \f$\displaystyle \nu = \frac{\gamma}{2\pi}\frac{V_{sp}}{V_{min}}\left(2\pi K_BT\right)^\frac{P_{min}-P_{SP}}{2}\sqrt{\sum_i\frac{a_i^2}{\epsilon_i}}\sqrt{\frac{\det H_{min}}{\det H_{SP}}}\f$
   !>
   !> where \f$\det H_k\f$ is the determinant of the Hessian matrix, \f$ V_k\f$ is the volume of
   !> the associated Goldstone mode, \f$P_k\f$ is the number of Goldstone modes, \f$\epsilon_i\f$ are the
   !> eigenvalues of teh Hessian at the saddle point, \f$a_i\f$ are expansion coefficients in the linearized 
   !> equation for the unstable mode derived from the Landau-Lifshitz equations of motion for 2N degrees of 
   !> freedom of the system. The modes \f$k\f$ indicate wheter the present states is a minimum \f$\min\f$ or
   !> the saddle point \f$SP\f$.
   !>
   !> The prefactor is in units of \f$\displaystyle\frac{1}{sT}\left(mRy\right)^\frac{P_{min}-P_{SP}}{2} \f$
   !>
   !> For more details look at Bessarab et al. Sci. Rep. Vol. 8, 3433 (2018)
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_pref(Natom,nmin,nsp,vs,vm,eig_min,eig_sp,vel_par,pref,expo,info)

      implicit none

      integer, intent(in) :: nsp    !< Number of eigenvectors in the saddle point
      integer, intent(in) :: nmin   !< Number of eignevectors for a state
      integer, intent(in) :: Natom  !< Number of atoms in system
      real(dblprec), intent(in) :: vs  !< Volume of the saddle point
      real(dblprec), intent(in) :: vm  !< Volume of a state
      real(dblprec), dimension(2*Natom), intent(in) :: eig_sp  !< Eigenvalues of the saddle point
      real(dblprec), dimension(2*Natom), intent(in) :: eig_min !< Eigenvalues of either Initial/Final states
      real(dblprec), dimension(2*Natom), intent(in) :: vel_par  !< Velocity parameter
      ! .. Output variables
      integer, intent(out) :: info
      real(dblprec), intent(out) :: pref  !< Prefactor for the determination of the rate
      real(dblprec), intent(out) :: expo  !< Exponent \f$\frac{P_{min}-P_{SP}}{2}\f$ Needed to calculate the transition rates
      ! .. Local variables
      integer :: ii,jj
      real(dblprec) :: s,m,kb_mRy,hb_mRy,me

      kb_mRy = k_bolt_ev*1000.0_dblprec/ry_ev
      hb_mRy = hbar_mev/ry_ev

      expo = 0.50_dblprec*real(nmin-nsp,dblprec)
      info = 0

      me = 1.0_dblprec

      ! Calculation of the factor (2piK_B)^{Pmin-Psp}/2
      do ii=1,(nmin-nsp)
         me = me*2.0_dblprec*pi*kb_mRy
      end do

      me = sqrt(me)

      ! Calculation of the factor sqrt(\det H_min / \det H_SP)
      m = 1.0_dblprec
      if (nmin>(nsp+1)) then
         do ii = 1,(nmin-nsp-1)
            m = m/sqrt(eig_sp(ii+1+nsp))
         end do
         jj = nmin+1
      elseif (nmin<(nsp+1)) then
         do ii = 1,(nsp+1-nmin)
            m = m*sqrt(eig_min(ii+nmin))
         end do
         jj = nsp+2
      else
         m = 1.0_dblprec
         jj = nmin+1
      end if

      if ((jj>2*Natom).or.(nsp+2>2*Natom)) then
         info = jj
      end if

      do ii=jj,2*Natom
         m = m*sqrt(eig_min(ii)/eig_sp(ii))
      end do

      s = 0.0_dblprec
      ! This seems to be the calculation of the "a_i*a_i/epsilon_i" term in the prefactor
      ! So the velocity parameters are the expansion coefficients for the unstable solution
      ! of the LL equation
      do ii = (nsp+2),2*Natom
         s = s + vel_par(ii)*vel_par(ii)/eig_sp(ii)
      end do
      s = sqrt(s)
      pref = g_e_abs*m*s*me*vs/(hb_mRy*2.0_dblprec*pi*vm)

   end subroutine calc_pref

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_pref
   !> @brief Calculates volume of zero mode
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_zero_vol(Natom,NN,emomM,pos,basis,bc,vol)

      use VPO, only: calc_ang

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, dimension(3), intent(in) :: NN   !< Number of repetitions of the unit cell
      real(dblprec), dimension(3,Natom), intent(in)   :: pos !< Coordinates of atoms
      real(dblprec), dimension(3,3), intent(in)       :: basis !< Basis
      real(dblprec), dimension(3,Natom), intent(in)   :: emomM !< Magnetic moment vector
      character(len=1), dimension(3), intent(in) :: bc !< Boundary conditions
      ! .. Output variables
      real(dblprec), intent(out) :: vol !< volumen of the zero mode
      ! .. Local allocatable arrays
      integer, dimension(:), allocatable :: kk
      integer, dimension(:), allocatable :: ipiv
      real(dblprec), dimension(:), allocatable :: v
      real(dblprec), dimension(:), allocatable :: ll
      real(dblprec), dimension(:), allocatable :: auxv
      real(dblprec), dimension(:,:), allocatable :: auxM
      real(dblprec), dimension(:,:), allocatable :: auxb
      ! .. Local variables
      integer :: ii,jj,k,l,info,num,i_stat,i_all
      real(dblprec) :: ll1,ll2,ll3
      character(len=8) :: d

      allocate(v(3),stat=i_stat)
      call memocc(i_stat,product(shape(v))*kind(v),'v','calc_zero_vol')
      v=0.0_dblprec
      allocate(ll(3),stat=i_stat)
      call memocc(i_stat,product(shape(ll))*kind(ll),'ll','calc_zero_vol')
      ll=0.0_dblprec
      allocate(kk(3),stat=i_stat)
      call memocc(i_stat,product(shape(kk))*kind(kk),'kk','calc_zero_vol')
      kk=0
      allocate(auxv(3),stat=i_stat)
      call memocc(i_stat,product(shape(auxv))*kind(auxv),'auxv','calc_zero_vol')
      auxv=0.0_dblprec
      allocate(ipiv(3),stat=i_stat)
      call memocc(i_stat,product(shape(ipiv))*kind(ipiv),'ipiv','calc_zero_vol')
      ipiv=0
      allocate(auxb(3,3),stat=i_stat)
      call memocc(i_stat,product(shape(auxb))*kind(auxb),'auxb','calc_zero_vol')
      auxb=0.0_dblprec
      allocate(auxM(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(auxM))*kind(auxM),'auxM','calc_zero_vol')
      auxM=0.0_dblprec

      vol = 1.0_dblprec

      !!$omp parallel do default(shared), private(ii,jj), collapse(2)
      do ii=1,3
         do jj=1,3
            auxb(ii,jj) = basis(ii,jj)
         end do
      end do
      !!$omp end parallel do

      call dgetrf(3,3,auxb,3,ipiv,info)

      if (info.ne.0) then
         write(6,'(a)') 'WARNING: Failed to calculate inverse of the basis!'
         return
      end if

      do ii=1,3
         kk(ii) = 0
         ll(ii) = 0.0_dblprec
      end do
      k = 0
      do ii=1,3
         if (bc(ii)=='P') then
            k = k+1
            kk(k) = ii

            do jj=1,Natom
               do l=1,3
                  v(l) = pos(l,jj) + basis(l,ii)
               end do
               call dgetrs('N',3,1,auxb,3,ipiv,v,3,info)
               if (info.ne.0) then
                  write(6,'(a)') 'WARNING: Failed to calculate inverse of the basis!'
                  return
               end if
               do l=1,3
                  if (nint(v(l))>NN(l)-1) then
                     v(l) = 0.0_dblprec
                  elseif (nint(v(l))<0) then
                     v(l) = real(NN(1),dblprec)-1.0_dblprec
                  end if
                  auxv(l) = v(l)
               end do

               call dgemv('N',3,3,1.0_dblprec,basis,3,auxv,1,0.0_dblprec,v,1)

               num = find_pos(Natom,v,pos)

               if (num<0) then
                  write(d,'(I8)') ii
                  write(6,'(a)') 'WARNING: Failed to move system along zero mode'//trim(adjustl(d))//'!'
                  return
               end if

               do l=1,3
                  auxM(l,num) = emomM(l,jj)
               end do
            end do

            do jj=1,Natom
               ll1 = calc_ang(emomM(:,jj),auxM(:,jj))
               ll(k) = ll(k) + ll1*ll1
            end do
            ll(k) = sqrt(ll(k))
            write(*,'(2x,a,2x,es16.8)')'zero mode length:',ll(k)
         end if
      end do

      !!$omp parallel do default(shared), private(ii,jj), collapse(2)
      do ii=1,3
         do jj=1,3
            auxb(ii,jj) = basis(ii,jj)/norm2(basis(:,jj))
         end do
      end do
      !!$omp end parallel do

      if (k==1) then
         vol = ll(1)
      elseif (k==2) then
         ll1 = auxb(2,kk(1))*auxb(3,kk(2))-auxb(3,kk(1))*auxb(2,kk(2))
         ll2 = auxb(1,kk(1))*auxb(3,kk(2))-auxb(3,kk(1))*auxb(1,kk(2))
         ll3 = auxb(1,kk(1))*auxb(2,kk(2))-auxb(2,kk(1))*auxb(1,kk(2))
         vol = ll(1)*ll(2)*sqrt(ll1*ll1+ll2*ll2+ll3*ll3)
      elseif (k==3) then
         ll1 = auxb(2,kk(1))*auxb(3,kk(2))-auxb(3,kk(1))*auxb(2,kk(2))
         ll2 = auxb(1,kk(1))*auxb(3,kk(2))-auxb(3,kk(1))*auxb(1,kk(2))
         ll3 = auxb(1,kk(1))*auxb(2,kk(2))-auxb(2,kk(1))*auxb(1,kk(2))
         vol = ll(1)*ll(2)*ll(3)*abs(auxb(1,kk(3))*ll1+auxb(2,kk(3))*ll2+auxb(3,kk(3))*ll3)
      else
         vol = 1.0_dblprec
      end if

      i_all=-product(shape(v))*kind(v)
      deallocate(v,stat=i_stat)
      call memocc(i_stat,i_all,'v','calc_zero_vol')
      i_all=-product(shape(ll))*kind(ll)
      deallocate(ll,stat=i_stat)
      call memocc(i_stat,i_all,'ll','calc_zero_vol')
      i_all=-product(shape(kk))*kind(kk)
      deallocate(kk,stat=i_stat)
      call memocc(i_stat,i_all,'kk','calc_zero_vol')
      i_all=-product(shape(ipiv))*kind(ipiv)
      deallocate(ipiv,stat=i_stat)
      call memocc(i_stat,i_all,'ipiv','calc_zero_vol')
      i_all=-product(shape(auxv))*kind(auxv)
      deallocate(auxv,stat=i_stat)
      call memocc(i_stat,i_all,'auxv','calc_zero_vol')
      i_all=-product(shape(auxb))*kind(auxb)
      deallocate(auxb,stat=i_stat)
      call memocc(i_stat,i_all,'auxb','calc_zero_vol')
      i_all=-product(shape(auxM))*kind(auxM)
      deallocate(auxM,stat=i_stat)
      call memocc(i_stat,i_all,'auxv','calc_zero_vol')

   end subroutine calc_zero_vol

   !----------------------------------------------------------------------------
   ! SUBROUTINE: calc_pref
   !> @brief Calculates volume of zero mode in AFM
   !> @author Pavel Bessarab
   !----------------------------------------------------------------------------
   subroutine calc_zero_vol_afm(Natom,NN,emomM,pos,basis,bc,vol)

      use VPO, only: calc_ang

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, dimension(3), intent(in) :: NN !< Number of repetitions of the unit cell
      real(dblprec), dimension(3,Natom), intent(in)   :: pos !< Coordinates of atoms
      real(dblprec), dimension(3,3), intent(in)       :: basis !< Basis
      real(dblprec), dimension(3,Natom), intent(in)   :: emomM !< Magnetic moment vector
      character(len=1), dimension(3), intent(in) :: bc   !< Boundary conditions
      ! .. Output variables
      real(dblprec), intent(out) :: vol !< Volume of the zero mode in AFM configurration
      ! .. Local allocatable arrays
      integer, dimension(:), allocatable :: kk
      integer, dimension(:), allocatable :: ipiv
      real(dblprec), dimension(:), allocatable :: v
      real(dblprec), dimension(:), allocatable :: ll
      real(dblprec), dimension(:), allocatable :: auxv
      real(dblprec), dimension(:,:), allocatable :: auxM
      real(dblprec), dimension(:,:), allocatable :: auxb
      ! .. Local variables
      integer :: ii,jj,k,l,info,num,i_stat,i_all
      character(len=8) :: d
      real(dblprec) :: ll1,ll2,ll3

      allocate(v(3),stat=i_stat)
      call memocc(i_stat,product(shape(v))*kind(v),'v','calc_zero_vol_afm')
      v=0.0_dblprec
      allocate(ll(3),stat=i_stat)
      call memocc(i_stat,product(shape(ll))*kind(ll),'ll','calc_zero_vol_afm')
      ll=0.0_dblprec
      allocate(kk(3),stat=i_stat)
      call memocc(i_stat,product(shape(kk))*kind(kk),'kk','calc_zero_vol_afm')
      kk=0
      allocate(auxv(3),stat=i_stat)
      call memocc(i_stat,product(shape(auxv))*kind(auxv),'auxv','calc_zero_vol_afm')
      auxv=0.0_dblprec
      allocate(ipiv(3),stat=i_stat)
      call memocc(i_stat,product(shape(ipiv))*kind(ipiv),'ipiv','calc_zero_vol_afm')
      ipiv=0
      allocate(auxb(3,3),stat=i_stat)
      call memocc(i_stat,product(shape(auxb))*kind(auxb),'auxb','calc_zero_vol_afm')
      auxb=0.0_dblprec
      allocate(auxM(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(auxM))*kind(auxM),'auxM','calc_zero_vol_afm')
      auxM=0.0_dblprec

      vol = 1.0_dblprec

      !!$omp parallel do default(shared), private(ii,jj), collapse(2)
      do ii=1,3
         do jj=1,3
            auxb(ii,jj) = basis(ii,jj)
         end do
      end do
      !!$omp end parallel do

      call dgetrf(3,3,auxb,3,ipiv,info)

      if (info.ne.0) then
         write(6,'(a)') 'WARNING: Failed to calculate inverse of the basis!'
         return
      end if

      k = 0
      do ii=1,3
         if (bc(ii)=='P') then
            k = k+1
            kk(k) = ii

            do jj=1,Natom
               do l=1,3
                  v(l) = pos(l,jj) + basis(l,ii)
               end do
               call dgetrs('N',3,1,auxb,3,ipiv,v,3,info)
               if (info.ne.0) then
                  write(6,'(a)') 'WARNING: Failed to calculate inverse of the basis!'
                  return
               end if
               do l=1,3
                  if (nint(v(l))>NN(l)-1) then
                     v(l) = 0.0_dblprec
                  elseif (nint(v(l))<0) then
                     v(l) = real(NN(1),dblprec)-1.0_dblprec
                  end if
                  auxv(l) = v(l)
               end do

               call dgemv('N',3,3,1.0_dblprec,basis,3,auxv,1,0.0_dblprec,v,1)

               num = find_pos(Natom,v,pos)

               if (num<0) then
                  write(d,'(I8)') ii
                  write(6,'(a)') 'WARNING: Failed to move system along zero mode'//trim(adjustl(d))//'!'
                  return
               end if

               do l=1,3
                  auxM(l,num) = -emomM(l,jj)
               end do
            end do

            do jj=1,Natom
               ll1 = calc_ang(emomM(:,jj),auxM(:,jj))
               ll(k) = ll(k) + ll1*ll1
            end do
            ll(k) = sqrt(ll(k))
            write(*,'(2x,a,2x,es16.8)')'zero mode length:',ll(k)
         end if
      end do

      !!$omp parallel do default(shared), private(ii,jj), collapse(2)
      do jj=1,3
         do ii=1,3
            auxb(ii,jj) = basis(ii,jj)/norm2(basis(:,jj))
         end do
      end do
      !!$omp end parallel do

      if (k==1) then
         vol = ll(1)
      elseif (k==2) then
         ll1 = auxb(2,kk(1))*auxb(3,kk(2))-auxb(3,kk(1))*auxb(2,kk(2))
         ll2 = auxb(1,kk(1))*auxb(3,kk(2))-auxb(3,kk(1))*auxb(1,kk(2))
         ll3 = auxb(1,kk(1))*auxb(2,kk(2))-auxb(2,kk(1))*auxb(1,kk(2))
         vol = ll(1)*ll(2)*sqrt(ll1*ll1+ll2*ll2+ll3*ll3)
      elseif (k==3) then
         ll1 = auxb(2,kk(1))*auxb(3,kk(2))-auxb(3,kk(1))*auxb(2,kk(2))
         ll2 = auxb(1,kk(1))*auxb(3,kk(2))-auxb(3,kk(1))*auxb(1,kk(2))
         ll3 = auxb(1,kk(1))*auxb(2,kk(2))-auxb(2,kk(1))*auxb(1,kk(2))
         vol = ll(1)*ll(2)*ll(3)*abs(auxb(1,kk(3))*ll1+auxb(2,kk(3))*ll2+auxb(3,kk(3))*ll3)
      else
         vol = 1.0_dblprec
      end if

      i_all=-product(shape(v))*kind(v)
      deallocate(v,stat=i_stat)
      call memocc(i_stat,i_all,'v','calc_zero_vol_afm')
      i_all=-product(shape(ll))*kind(ll)
      deallocate(ll,stat=i_stat)
      call memocc(i_stat,i_all,'ll','calc_zero_vol_afm')
      i_all=-product(shape(kk))*kind(kk)
      deallocate(kk,stat=i_stat)
      call memocc(i_stat,i_all,'kk','calc_zero_vol_afm')
      i_all=-product(shape(ipiv))*kind(ipiv)
      deallocate(ipiv,stat=i_stat)
      call memocc(i_stat,i_all,'ipiv','calc_zero_vol_afm')
      i_all=-product(shape(auxv))*kind(auxv)
      deallocate(auxv,stat=i_stat)
      call memocc(i_stat,i_all,'auxv','calc_zero_vol_afm')
      i_all=-product(shape(auxb))*kind(auxb)
      deallocate(auxb,stat=i_stat)
      call memocc(i_stat,i_all,'auxb','calc_zero_vol_afm')
      i_all=-product(shape(auxM))*kind(auxM)
      deallocate(auxM,stat=i_stat)
      call memocc(i_stat,i_all,'auxv','calc_zero_vol_afm')

   end subroutine calc_zero_vol_afm

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: find_pos
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   function find_pos(Natom,v,pos)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3), intent(in) :: v
      real(dblprec), dimension(3,Natom), intent(in) :: pos
      ! .. Local variables
      integer :: ii
      integer :: find_pos
      real(dblprec) :: e

      e = 0.00000010_dblprec

      find_pos = -1

      do ii=1,Natom
         if ((abs(v(1)-pos(1,ii))<e).and.(abs(v(2)-pos(2,ii))<e).and.(abs(v(3)-pos(3,ii))<e)) then
            find_pos = ii
            return
         end if
      end do

   end function find_pos

    !---------------------------------------------------------------------------------
    ! SUBROUTINE: eigen_driver
    !> @brief 
    !> @author Anders Bergman
    !---------------------------------------------------------------------------------
    subroutine eigen_driver(inmat,lda,eigvals,info,eigvecs,m)
 
       implicit none
 
       integer, intent(in) :: lda !< dimension (1d) of matrix
       real(dblprec), dimension(lda,lda), intent(in) :: inmat !< input matrix
       integer, intent(in), optional :: m !< number of eigenvalues to calculate
       ! .. Output variables
       real(dblprec), dimension(lda), intent(out) :: eigvals !< eigenvalues
       integer, intent(out) :: info
       real(dblprec), dimension(lda,*), intent(out), optional :: eigvecs !< eigenvectors

       ! .. Local variables
       integer :: lwork, liwork
       integer :: i_stat,i_all
       integer :: n_eig
       real(dblprec), dimension(1) :: ww
       integer, dimension(1) :: iww
       real(dblprec) :: vl,vu,abstol
       integer :: il, iu
 
       ! .. Local allocatable variables
       real(dblprec), dimension(:), allocatable :: work
       integer, dimension(:), allocatable :: iwork, isuppz
       real(dblprec), dimension(:,:), allocatable :: hwork, zwork
 
       allocate(isuppz(2*lda),stat=i_stat)
       call memocc(i_stat,product(shape(isuppz))*kind(isuppz),'isuppz','eigen_driver')
       allocate(hwork(lda,lda),stat=i_stat)
       call memocc(i_stat,product(shape(hwork))*kind(hwork),'hwork','eigen_driver')
 
       call dcopy(lda*lda,inmat,1,hwork,1)
 
       vl=0.0_dblprec;vu=0.0_dblprec;abstol=epsilon(vl)
       il=1;iu=1;
       lwork=-1
       liwork=-1
       if(present(eigvecs)) then
          if(present(m)) then
             iu=m
             call dsyevr('V', 'I', 'U', lda, hwork, lda, vl, vu, il, iu, abstol, n_eig, eigvals, eigvecs, lda, isuppz, ww, lwork, iww, liwork, info)
          else
             call dsyevr('V', 'A', 'U', lda, hwork, lda, vl, vu, il, iu, abstol, n_eig, eigvals, eigvecs, lda, isuppz, ww, lwork, iww, liwork, info)
          end if
       else
          allocate(zwork(lda,lda),stat=i_stat)
          call memocc(i_stat,product(shape(zwork))*kind(zwork),'zwork','eigen_driver')
          call dsyevr('N', 'A', 'U', lda, hwork, lda, vl, vu, il, iu, abstol, n_eig, eigvals, zwork, lda, isuppz, ww, lwork, iww, liwork, info)
       end if

          lwork = int(ww(1))
          allocate(work(lwork),stat=i_stat)
          call memocc(i_stat,product(shape(work))*kind(work),'work','eigen_driver')
          liwork = int(iww(1))
          write(*,'(x,i8,x,i8,x,i4)') lwork,liwork,info
          allocate(iwork(liwork),stat=i_stat)
          call memocc(i_stat,product(shape(iwork))*kind(iwork),'iwork','eigen_driver')

       if(present(eigvecs)) then
          if(present(m)) then
             iu=m
             call dsyevr('V', 'I', 'U', lda, hwork, lda, vl, vu, il, iu, abstol, n_eig, eigvals, eigvecs, lda, isuppz, work, lwork, iwork, liwork, info)
          else
             call dsyevr('V', 'A', 'U', lda, hwork, lda, vl, vu, il, iu, abstol, n_eig, eigvals, eigvecs, lda, isuppz, work, lwork, iwork, liwork, info)
          end if
       else
          call dsyevr('N', 'A', 'U', lda, hwork, lda, vl, vu, il, iu, abstol, n_eig, eigvals, zwork, lda, isuppz, work, lwork, iwork, liwork, info)
          i_all=-product(shape(zwork))*kind(zwork)
          deallocate(zwork,stat=i_stat)
          call memocc(i_stat,i_all,'zwork','eigen_driver')
       end if

 
       i_all=-product(shape(work))*kind(work)
       deallocate(work,stat=i_stat)
       call memocc(i_stat,i_all,'work','eigen_driver')
       i_all=-product(shape(hwork))*kind(hwork)
       deallocate(hwork,stat=i_stat)
       call memocc(i_stat,i_all,'hwork','eigen_driver')
       i_all=-product(shape(iwork))*kind(iwork)
       deallocate(iwork,stat=i_stat)
       call memocc(i_stat,i_all,'iwork','eigen_driver')
       i_all=-product(shape(isuppz))*kind(isuppz)
       deallocate(isuppz,stat=i_stat)
       call memocc(i_stat,i_all,'isuppz','eigen_driver')
 
    end subroutine eigen_driver

end module Hessian
