!------------------------------------------------------------------------------------
!> @brief
!> Intended to use for printing all the relevant information relating
!> the ASD parameters to micromagnetic variables
!
!> @author
!> Jonathan Chico
!> @copyright
!> GNU Public License
!------------------------------------------------------------------------------------
module prn_micromagnetic
   use Parameters
   use Profiling
   use Stiffness

   implicit none

   private

   public :: stiffness_wrapper, init_stiffness

contains

   !---------------------------------------------------------------------------------
   !> @brief
   !> Wrapper for writing all the micromagnetic information
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine stiffness_wrapper(NT,NA,N1,N2,N3,hdim,mconf,Natom,nHam,Nchmax,        &
      conf_num,do_ralloy,Natom_full,max_no_neigh,Nch,anumb,atype,aham,nlistsize,    &
      atype_ch,asite_ch,achem_ch,nlist,alat,C1,C2,C3,coord,chconc,ammom_inp,ncoup,  &
      max_no_dmneigh,dmlistsize,dmlist,dm_vect,do_anisotropy,anisotropy,simid)

      implicit none

      !.. Input variables
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: hdim   !< Number of elements in Hamiltonian element (scalar or vector)
      integer, intent(in) :: mconf  !< LSF ground state configuration
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: nHam   !< Number of atoms in Hamiltonian
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: conf_num !< Number of LSF configurations
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: max_no_dmneigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: do_anisotropy  !< Read anisotropy data (1/0)
      integer, dimension(NA), intent(in) :: Nch !< Number of chemical components on each site in cell
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, dimension(Natom), intent(in) :: aham  !< Hamiltonian look-up table
      integer, dimension(nHam), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist  !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist !< Neighbour list for Heisenberg exchange couplings

      real(dblprec), intent(in) :: alat !< Lattice parameter
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(NA,Nchmax), intent(in) :: chconc !< Chemical concentration on sites
      real(dblprec), dimension(NA,6,Nchmax), intent(in) :: anisotropy !< Input data for anisotropies
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(hdim,max_no_neigh,nHam,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(3,max_no_dmneigh,nHam), intent(in) :: dm_vect !< Heisenberg exchange couplings

      character(len=8), intent(in) :: simid !< Name of simulation

      ! Allocate variables needed for the calculation of the stiffness
      call allocate_stiffness(NA,natom,1,do_stiffness,do_dm_stiffness)

      write(*,'(1x,a)',advance='no') 'Calculate the ferromagnetic stiffness'
      ! Calculate the excnage stiffness
      call ferro_stiffness(NT,NA,N1,N2,N3,1,mconf,Natom,nHam, Nchmax,eta_max,       &
         eta_min,conf_num,do_ralloy,Natom_full,max_no_neigh,Nch,anumb,atype,        &
         nlistsize,atype_ch,asite_ch,achem_ch,nlist,alat,C1,C2,C3,coord,chconc,     &
         ammom_inp,ncoup,aham,A_xc,M_sat,Dxc_fit,cell_vol,A_xc_lsq,D_err_fit,       &
         Dxc_fit_lsq,J0_matrix,D_xc_stiffness_matrix,A_xc_stiffness_matrix,         &
         D_xc_stiffness_matrix_lsq,A_xc_stiffness_matrix_lsq)
      write(*,'(1x,a)') 'done'

      ! Check if the DMI stiffness is also to be calculated
      if (do_dm_stiffness=='Y') then

         write(*,'(1x,a)',advance='no') 'Calculate the DMI spiralization'
         call DMI_stiffness(NT,NA,N1,N2,N3,Natom,nHam,Nchmax,eta_max,eta_min,       &
            max_no_dmneigh,anumb,dmlistsize,dmlist,alat,coord,ammom_inp,dm_vect,    &
            DM0_mat,DM0_mat_lsq,aham)
         write(*,'(1x,a)') 'done'
      endif

     if (do_ralloy==1) then
      write(*,'(1x,a)',advance='no') 'Calculate the ferromagnetic stiffness of random alloy'
      ! Calculate the excnage stiffness
      call ferro_random_stiffness(NT,NA,N1,N2,N3,1,mconf,Natom,nHam, Nchmax,eta_max,&
         eta_min,conf_num,do_ralloy,Natom_full,max_no_neigh,Nch,anumb,atype,        &
         nlistsize,atype_ch,asite_ch,achem_ch,nlist,alat,C1,C2,C3,coord,chconc,     &
         ammom_inp,ncoup,aham,Axc_fit_alloy,Dxc_fit_alloy,Tc_alloy,J0_matrix_alloy)
      write(*,'(1x,a)') 'done'
     endif

      ! Print the micromagnetic information
      call prn_micro_wrapper(NA,natom,Nchmax,eta_min,eta_max,do_anisotropy,A_xc,    &
         M_sat,Dxc_fit,cell_vol,A_xc_lsq,D_err_fit,Dxc_fit_lsq,                     &
         D_xc_stiffness_matrix,A_xc_stiffness_matrix,D_xc_stiffness_matrix_lsq,     &
         A_xc_stiffness_matrix_lsq,anisotropy,do_dm_stiffness,simid,J0_matrix,      &
         DM0_mat,DM0_mat_lsq,prn_J0_matrix,Axc_fit_alloy,Dxc_fit_alloy,Tc_alloy,    &
         J0_matrix_alloy,do_ralloy)


      ! Deallocate the arrays for the stiffness
      call allocate_stiffness(NA,natom,-1,do_stiffness,do_dm_stiffness)

   end subroutine stiffness_wrapper

   !---------------------------------------------------------------------------------
   !> @brief
   !> Initializtion of the stiffnes default variables
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine init_stiffness()

      implicit none

      eta_max=0
      eta_min=0

      do_stiffness='N'
      do_dm_stiffness='N'
      prn_J0_matrix='N'

      A_xc        = 0.0_dblprec
      M_sat       = 0.0_dblprec
      Dxc_fit     = 0.0_dblprec
      cell_vol    = 0.0_dblprec
      A_xc_lsq    = 0.0_dblprec
      D_err_fit   = 0.0_dblprec
      Dxc_fit_lsq = 0.0_dblprec

   end subroutine init_stiffness

   !--------------------------------------------------------------------------------
   !> @brief
   !>  Allocate the necessary arrays to calculate the stiffness
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine allocate_stiffness(NA,natom,flag,do_stiffness,do_dm_stiffness)

      implicit none

      ! .. Input variables
      integer, intent(in) :: NA !< Number of atoms in one cell
      integer, intent(in) :: natom !< Number of atoms
      integer, intent(in) :: flag !< Flag for allocation for arrays
      character(len=1), intent(in) :: do_stiffness !< Calculate the spin wave stiffness of a ferromagnet
      character(len=1), intent(in) :: do_dm_stiffness !< Calculate the DMI spiralization

      ! .. Local variables
      integer :: i_stat, i_all

      if (flag.gt.0) then

         ! Allocation of arrays for the exchange stiffness
         if ( do_stiffness.eq.'Y') then
            allocate(J0_matrix(NA,NA),stat=i_stat)
            call memocc(i_stat,product(shape(J0_matrix))*kind(J0_matrix),'J0_matrix','allocate_stiffness')
            J0_matrix=0.0_dblprec
            allocate(D_xc_stiffness_matrix(3,3),stat=i_stat)
            call memocc(i_stat,product(shape(D_xc_stiffness_matrix))*kind(D_xc_stiffness_matrix),'D_xc_stiffness_matrix','allocate_stiffness')
            D_xc_stiffness_matrix=0.0_dblprec
            allocate(A_xc_stiffness_matrix(3,3),stat=i_stat)
            call memocc(i_stat,product(shape(A_xc_stiffness_matrix))*kind(A_xc_stiffness_matrix),'A_xc_stiffness_matrix','allocate_stiffness')
            A_xc_stiffness_matrix=0.0_dblprec
            allocate(D_xc_stiffness_matrix_lsq(3,3),stat=i_stat)
            call memocc(i_stat,product(shape(D_xc_stiffness_matrix_lsq))*kind(D_xc_stiffness_matrix_lsq),'D_xc_stiffness_matrix_lsq','allocate_stiffness')
            D_xc_stiffness_matrix_lsq=0.0_dblprec
            allocate(A_xc_stiffness_matrix_lsq(3,3),stat=i_stat)
            call memocc(i_stat,product(shape(A_xc_stiffness_matrix_lsq))*kind(A_xc_stiffness_matrix_lsq),'A_xc_stiffness_matrix_lsq','allocate_stiffness')
            A_xc_stiffness_matrix_lsq=0.0_dblprec
            ! Allocation of the DM spiralization matrix
            if (do_dm_stiffness.eq.'Y') then
               allocate(DM0_mat(3,3),stat=i_stat)
               call memocc(i_stat,product(shape(DM0_mat))*kind(DM0_mat),'DM0_mat','allocate_stiffness')
               DM0_mat=0.0_dblprec
               allocate(DM0_mat_lsq(3,3),stat=i_stat)
               call memocc(i_stat,product(shape(DM0_mat_lsq))*kind(DM0_mat_lsq),'DM0_mat_lsq','allocate_stiffness')
               DM0_mat=0.0_dblprec
            endif
!            if (do_ralloy==1) then
              allocate(J0_matrix_alloy(NA,NA,natom),stat=i_stat)
              call memocc(i_stat,product(shape(J0_matrix_alloy))*kind(J0_matrix_alloy),'J0_matrix_alloy','allocate_stiffness')
              J0_matrix_alloy=0.0_dblprec
              allocate(Axc_fit_alloy(natom),stat=i_stat)
              call memocc(i_stat,product(shape(Axc_fit_alloy))*kind(Axc_fit_alloy),'Axc_fit_alloy','allocate_stiffness')
              Axc_fit_alloy=0.0_dblprec
              allocate(Dxc_fit_alloy(natom),stat=i_stat)
              call memocc(i_stat,product(shape(Dxc_fit_alloy))*kind(Dxc_fit_alloy),'Dxc_fit_alloy','allocate_stiffness')
              Dxc_fit_alloy=0.0_dblprec
              allocate(Tc_alloy(natom),stat=i_stat)
              call memocc(i_stat,product(shape(Tc_alloy))*kind(Tc_alloy),'Tc_alloy','allocate_stiffness')
              Tc_alloy=0.0_dblprec
!            endif
         endif
      else
         ! Deallocate arrays for the exchange stiffness
         if ( do_stiffness.eq.'Y') then
            i_all=-product(shape(J0_matrix))*kind(J0_matrix)
            deallocate(J0_matrix,stat=i_stat)
            call memocc(i_stat,i_all,'J0_matrix','allocate_stiffness')
            i_all=-product(shape(D_xc_stiffness_matrix))*kind(D_xc_stiffness_matrix)
            deallocate(D_xc_stiffness_matrix,stat=i_stat)
            call memocc(i_stat,i_all,'D_xc_stiffness_matrix','allocate_stiffness')
            i_all=-product(shape(A_xc_stiffness_matrix))*kind(A_xc_stiffness_matrix)
            deallocate(A_xc_stiffness_matrix,stat=i_stat)
            call memocc(i_stat,i_all,'A_xc_stiffness_matrix','allocate_stiffness')
            i_all=-product(shape(D_xc_stiffness_matrix_lsq))*kind(D_xc_stiffness_matrix_lsq)
            deallocate(D_xc_stiffness_matrix_lsq,stat=i_stat)
            call memocc(i_stat,i_all,'D_xc_stiffness_matrix_lsq','allocate_stiffness')
            i_all=-product(shape(A_xc_stiffness_matrix_lsq))*kind(A_xc_stiffness_matrix_lsq)
            deallocate(A_xc_stiffness_matrix_lsq,stat=i_stat)
            call memocc(i_stat,i_all,'A_xc_stiffness_matrix_lsq','allocate_stiffness')
            ! Deallocate DMI spiralization arrays
            if (do_dm_stiffness.eq.'Y') then
               i_all=-product(shape(DM0_mat))*kind(DM0_mat)
               deallocate(DM0_mat,stat=i_stat)
               call memocc(i_stat,i_all,'DM0_mat','allocate_stiffness')
               i_all=-product(shape(DM0_mat_lsq))*kind(DM0_mat_lsq)
               deallocate(DM0_mat_lsq,stat=i_stat)
               call memocc(i_stat,i_all,'DM0_mat_lsq','allocate_stiffness')
            endif
!            if(do_ralloy==1) then
               i_all=-product(shape(J0_matrix_alloy))*kind(J0_matrix_alloy)
               deallocate(J0_matrix_alloy,stat=i_stat)
               call memocc(i_stat,i_all,'J0_matrix_alloy','allocate_stiffness')
               i_all=-product(shape(Axc_fit_alloy))*kind(Axc_fit_alloy)
               deallocate(Axc_fit_alloy,stat=i_stat)
               call memocc(i_stat,i_all,'Axc_fit_alloy','allocate_stiffness')
               i_all=-product(shape(Dxc_fit_alloy))*kind(Dxc_fit_alloy)
               deallocate(Dxc_fit_alloy,stat=i_stat)
               call memocc(i_stat,i_all,'Dxc_fit_alloy','allocate_stiffness')
               i_all=-product(shape(Tc_alloy))*kind(Tc_alloy)
               deallocate(Tc_alloy,stat=i_stat)
               call memocc(i_stat,i_all,'Tc_alloy','allocate_stiffness')
!            endif
         endif
      endif


   end subroutine allocate_stiffness

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing wrapper for micromagnetic variables
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine prn_micro_wrapper(NA,natom,Nchmax,eta_min,eta_max,do_anisotropy,A_xc, &
         M_sat,Dxc_fit,cell_vol,A_xc_lsq,D_err_fit,Dxc_fit_lsq,                     &
         D_xc_stiffness_matrix,A_xc_stiffness_matrix,D_xc_stiffness_matrix_lsq,     &
         A_xc_stiffness_matrix_lsq,anisotropy,do_dm_stiffness,simid,J0_matrix,      &
         DM0_mat,DM0_mat_lsq,prn_J0_matrix,Axc_fit_alloy,Dxc_fit_alloy,Tc_alloy,    &
         J0_matrix_alloy,do_ralloy)
      !
      !.. Implicit declarations
      implicit none

      !.. Input variables
      integer, intent(in) :: NA !< Number of atoms in one cell
      integer, intent(in) :: natom !< Number of atoms
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: eta_min !< Minimum  convergence parameters for the stiffness
      integer, intent(in) :: eta_max !< Number of convergence parameters for the stiffness
      integer, intent(in) :: do_anisotropy !< Read anisotropy data (1/0)

      real(dblprec), intent(in) :: A_xc     !< Exchange stiffness constant from rational fit J/m
      real(dblprec), intent(in) :: M_sat    !< Saturation magnetization muB/m^3
      real(dblprec), intent(in) :: Dxc_fit  !< Spin wave stiffnes from rational fit
      real(dblprec), intent(in) :: cell_vol !< Unit cell volume m^3
      real(dblprec), intent(in) :: A_xc_lsq !< Exchange stiffness constant from LSQ fit
      real(dblprec), intent(in) :: D_err_fit   !< Error form spin wave stiffnesrational fit
      real(dblprec), intent(in) :: Dxc_fit_lsq !< Spin wave stiffness LSQ fit
      real(dblprec), dimension(NA,NA), intent(in) :: J0_matrix !< Matrix being used to calculate the Tc-MFA
      real(dblprec), dimension(3,3), intent(in) :: DM0_mat !< DMI spiralization matrix [meVA]
      real(dblprec), dimension(3,3), intent(in) :: DM0_mat_lsq !< DMI spiralization matrix (LSQ fit)  [meVA]
      real(dblprec), dimension(3,3), intent(in) :: D_xc_stiffness_matrix !< Spin wave stiffness tensor rational fit
      real(dblprec), dimension(3,3), intent(in) :: A_xc_stiffness_matrix !< Exchange stiffness constant rational fit
      real(dblprec), dimension(3,3), intent(in) :: D_xc_stiffness_matrix_lsq !< Spin wave stiffness tensor LSQ fit
      real(dblprec), dimension(3,3), intent(in) :: A_xc_stiffness_matrix_lsq !< Exchange stiffness constant LSQ fit
      real(dblprec), dimension(NA,6,Nchmax), intent(in) :: anisotropy !< Input data for anisotropies

      character(len=1), intent(in) :: do_dm_stiffness !< Calculate the DMI spiralization
      character(len=1), intent(in) :: prn_J0_matrix !< Print the full site dependent J0 matrix
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec),dimension(natom),intent(in) :: Axc_fit_alloy !< Exchange stiffness alloy
      real(dblprec),dimension(natom),intent(in) :: Dxc_fit_alloy !< Spin wave stiffness alloy
      real(dblprec),dimension(natom),intent(in) :: Tc_alloy  !< Tc-MFA alloy
      real(dblprec),dimension(na,na,natom),intent(in) :: J0_matrix_alloy !< Exchange matrix alloy
      integer,intent(in) :: do_ralloy

      !.. Local variables
      character(len=30) :: filn

      !.. Executable statements

      ! Open outputfile
      write (filn,'(''micro_ASD.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      ! Printing statements for different micromagnetic variables
      ! .. Basic parameters
      call print_micro_misc(NA,Nchmax,eta_min,eta_max,do_anisotropy,A_xc,M_sat,     &
         cell_vol,anisotropy)
      ! .. Exchange stiffness
      call print_stiffness(A_xc,Dxc_fit,D_err_fit,D_xc_stiffness_matrix,            &
         A_xc_stiffness_matrix,A_xc_lsq,Dxc_fit_lsq,D_xc_stiffness_matrix_lsq,      &
         A_xc_stiffness_matrix_lsq)

      ! .. DMI spiralization
      if (do_dm_stiffness.eq.'Y') then
         call print_dmi_stiffness(M_sat,DM0_mat,DM0_mat_lsq,D_xc_stiffness_matrix,  &
            D_xc_stiffness_matrix_lsq)
      endif

      ! .. J0 matrix printing
      if (prn_J0_matrix.eq.'Y') then
         call print_J0(NA, J0_matrix)
      endif

      call print_J0_vector(NA, J0_matrix)

      if(do_ralloy==1) then
        call print_random_stiffness(natom,na,Axc_fit_alloy,Dxc_fit_alloy,Tc_alloy,  &
         J0_matrix_alloy)
      endif

      close(ofileno)

   end subroutine prn_micro_wrapper

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing the exchange stiffness parameters
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine print_stiffness(A_xc,Dxc_fit,D_err_fit,D_xc_stiffness_matrix,         &
         A_xc_stiffness_matrix,A_xc_lsq,Dxc_fit_lsq,D_xc_stiffness_matrix_lsq,      &
         A_xc_stiffness_matrix_lsq)

      implicit none

      ! .. Input variables
      real(dblprec), intent(in) :: A_xc     !< Exchange stiffness constant from rational fit J/m
      real(dblprec), intent(in) :: Dxc_fit  !< Spin wave stiffnes from rational fit
      real(dblprec), intent(in) :: A_xc_lsq !< Exchange stiffness constant from LSQ fit
      real(dblprec), intent(in) :: D_err_fit !< Error form spin wave stiffnes rational fit
      real(dblprec), intent(in) :: Dxc_fit_lsq !< Spin wave stiffness LSQ fit
      real(dblprec), dimension(3,3), intent(in) :: D_xc_stiffness_matrix !< Spin wave stiffness tensor rational fit
      real(dblprec), dimension(3,3), intent(in) :: A_xc_stiffness_matrix !< Exchange stiffness constant rational fit
      real(dblprec), dimension(3,3), intent(in) :: D_xc_stiffness_matrix_lsq !< Spin wave stiffness tensor LSQ fit
      real(dblprec), dimension(3,3), intent(in) :: A_xc_stiffness_matrix_lsq !< Exchange stiffness constant LSQ fit

      ! Printing of exchange stiffness information obtained from the Jij's
      write (ofileno,'(a)') " "
      write (ofileno,'(2x,a,f10.3,2x,a)') 'Exchange stiffness at T=0 K:',A_xc*1e12,'pJ/m'
      write (ofileno,'(a)') "**************** EXCHANGE STIFFNESS MATRIX [pJ/m] **************"
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') A_xc_stiffness_matrix(1,1)*1e12,A_xc_stiffness_matrix(1,2)*1e12,A_xc_stiffness_matrix(1,3)*1e12
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') A_xc_stiffness_matrix(2,1)*1e12,A_xc_stiffness_matrix(2,2)*1e12,A_xc_stiffness_matrix(2,3)*1e12
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') A_xc_stiffness_matrix(3,1)*1e12,A_xc_stiffness_matrix(3,2)*1e12,A_xc_stiffness_matrix(3,3)*1e12
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)')  " "
      write (ofileno,'(2x,a,f10.3,a,f10.3,2x,a)') 'Spin wave stiffness from Jijs:',Dxc_fit,'±',D_err_fit,'meVÅ^2'
      write (ofileno,'(a)') "**************** SPIN WAVE STIFFNESS MATRIX [meVA^2] ***********"
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') D_xc_stiffness_matrix(1,1),D_xc_stiffness_matrix(1,2),D_xc_stiffness_matrix(1,3)
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') D_xc_stiffness_matrix(2,1),D_xc_stiffness_matrix(2,2),D_xc_stiffness_matrix(2,3)
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') D_xc_stiffness_matrix(3,1),D_xc_stiffness_matrix(3,2),D_xc_stiffness_matrix(3,3)
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)')  " "
      write (ofileno,'(2x,a,f10.3,2x,a)') 'Exchange stiffness at T=0 K (LSQ):',A_xc_lsq*1e12,'pJ/m'
      write (ofileno,'(a)') "**************** EXCHANGE STIFFNESS MATRIX [pJ/m] **************"
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') A_xc_stiffness_matrix_lsq(1,1)*1e12,A_xc_stiffness_matrix_lsq(1,2)*1e12,A_xc_stiffness_matrix_lsq(1,3)*1e12
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') A_xc_stiffness_matrix_lsq(2,1)*1e12,A_xc_stiffness_matrix_lsq(2,2)*1e12,A_xc_stiffness_matrix_lsq(2,3)*1e12
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') A_xc_stiffness_matrix_lsq(3,1)*1e12,A_xc_stiffness_matrix_lsq(3,2)*1e12,A_xc_stiffness_matrix_lsq(3,3)*1e12
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)')  " "
      write (ofileno,'(2x,a,f10.3,2x,a)') 'Spin wave stiffness from Jijs (LSQ):',Dxc_fit_lsq,'meVÅ^2'
      write (ofileno,'(a)') "**************** SPIN WAVE STIFFNESS MATRIX [meVA^2] ***********"
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') D_xc_stiffness_matrix_lsq(1,1),D_xc_stiffness_matrix_lsq(1,2),D_xc_stiffness_matrix_lsq(1,3)
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') D_xc_stiffness_matrix_lsq(2,1),D_xc_stiffness_matrix_lsq(2,2),D_xc_stiffness_matrix_lsq(2,3)
      write (ofileno,'(2x, f14.6,2x,f14.6,2x,f14.6)') D_xc_stiffness_matrix_lsq(3,1),D_xc_stiffness_matrix_lsq(3,2),D_xc_stiffness_matrix_lsq(3,3)
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)')  " "

   end subroutine print_stiffness

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing the exchange DMI stiffness/spiralization
   !> @author
   !> Jonathan Chico
   !> @todo check the implementation of the wavelength making use of the tensorial forms
   !---------------------------------------------------------------------------------
   subroutine print_dmi_stiffness(M_sat,DM0_mat,DM0_mat_lsq,D_xc_stiffness_matrix,  &
         D_xc_stiffness_matrix_lsq)

      use Constants

      implicit none

      ! .. Input variables
      real(dblprec), intent(in) :: M_sat !< Saturation magnetization muB/m^3
      real(dblprec), dimension(3,3), intent(in) :: DM0_mat !< DMI spiralization matrix [meVA]
      real(dblprec), dimension(3,3), intent(in) :: DM0_mat_lsq !< DMI spiralization matrix (LSQ fit) [meVA]
      real(dblprec), dimension(3,3), intent(in) :: D_xc_stiffness_matrix !< Spin wave stiffness tensor rational fit
      real(dblprec), dimension(3,3), intent(in) :: D_xc_stiffness_matrix_lsq !< Spin wave stiffness tensor LSQ fit

      ! .. Local variables
      real(dblprec), dimension(3,3) :: spiral
      real(dblprec), dimension(3,3) :: spiral_lsq

      integer :: ii, jj

      ! Calculating the spiral wavelength from the spin wave stiffness and the
      ! DMI spiralization
      do ii=1,3
         do jj=1,3
            spiral(ii,jj)=4*pi*D_xc_stiffness_matrix(ii,jj)/(DM0_mat(ii,jj)+1.0d-10)
            spiral_lsq(ii,jj)=4*pi*D_xc_stiffness_matrix_lsq(ii,jj)/(DM0_mat_lsq(ii,jj)+1.0d-10)
         enddo
      enddo

      write (ofileno,'(a)') " "
      write (ofileno,'(2x,a)') "DM matrix (q,n)   D_{ij} x r_{ij} "
      write (ofileno,'(a)') "**************** DMI SPIRALIZATION MATRIX [meVA] ***************"
      write (ofileno,'(1x,a)') "            r_x         r_y         r_z"
      write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_x ',DM0_mat(:,1)
      write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_y ',DM0_mat(:,2)
      write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_z ',DM0_mat(:,3)
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)') " "
      write (ofileno,'(a)') "******* DMI SPIRALIZATION MATRIX (LSQ fit) [meVA] **************"
      write (ofileno,'(1x,a)') "            r_x         r_y         r_z"
      write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_x ',DM0_mat_lsq(:,1)
      write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_y ',DM0_mat_lsq(:,2)
      write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_z ',DM0_mat_lsq(:,3)
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)') " "
      write (ofileno,'(a)') "**************** DMI SPIRALIZATION MATRIX [mJ/m^2] *************"
      write (ofileno,'(1x,a)') "            r_x         r_y         r_z"
      write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_x ',DM0_mat(:,1)*1e-10*M_sat/(Joule_ev*4.0_dblprec)
      write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_y ',DM0_mat(:,2)*1e-10*M_sat/(Joule_ev*4.0_dblprec)
      write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_z ',DM0_mat(:,3)*1e-10*M_sat/(Joule_ev*4.0_dblprec)
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)') " "
      !write (ofileno,'(a)') "**************** DMI SPIRALIZATION MATRIX (LSQ FIT) [mJ/m^2] *************"
      !write (ofileno,'(1x,a)') "            r_x         r_y         r_z"
      !write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_x ',DM0_mat_lsq(:,1)*1e-10*M_sat/(Joule_ev*4.0_dblprec)
      !write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_y ',DM0_mat_lsq(:,2)*1e-10*M_sat/(Joule_ev*4.0_dblprec)
      !write (ofileno,'(2x,a,f14.6,2x,f14.6,2x,f14.6)') ' D_z ',DM0_mat_lsq(:,3)*1e-10*M_sat/(Joule_ev*4.0_dblprec)
      !write (ofileno,'(a)') "****************************************************************"
      !write (ofileno,'(a)')  " "
      write (ofileno,'(a)') "**************** SPIRAL WAVELENGTH MATRIX (STIFF/DMI) [A] ******"
      write (ofileno,'(2x,G14.6,2x,G14.6,2x,G14.6)') spiral(1,1),spiral(1,2),spiral(1,3)
      write (ofileno,'(2x,G14.6,2x,G14.6,2x,G14.6)') spiral(2,1),spiral(2,2),spiral(2,3)
      write (ofileno,'(2x,G14.6,2x,G14.6,2x,G14.6)') spiral(3,1),spiral(3,2),spiral(3,3)
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)')  " "
      write (ofileno,'(a)') "*************** SPIRAL WAVELENGTH MATRIX LSQ (STIFF/DMI) [A] ***"
      write (ofileno,'(2x,G14.6,2x,G14.6,2x,G14.6)') spiral_lsq(1,1),spiral_lsq(1,2),spiral_lsq(1,3)
      write (ofileno,'(2x,G14.6,2x,G14.6,2x,G14.6)') spiral_lsq(2,1),spiral_lsq(2,2),spiral_lsq(2,3)
      write (ofileno,'(2x,G14.6,2x,G14.6,2x,G14.6)') spiral_lsq(3,1),spiral_lsq(3,2),spiral_lsq(3,3)
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)')  " "

   end subroutine print_dmi_stiffness

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing micromagnetic information for analysis
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine print_micro_misc(NA,Nchmax,eta_min,eta_max,do_anisotropy,A_xc,M_sat,  &
      cell_vol,anisotropy)

      use Constants

      implicit none

      ! .. Input variables
      integer, intent(in) :: NA      !< Number of atoms in one cell
      integer, intent(in) :: Nchmax  !< Max number of chemical components on each site in cell
      integer, intent(in) :: eta_min !< Minimum  convergence parameters for the stiffness
      integer, intent(in) :: eta_max !< Number of convergence parameters for the stiffness
      integer, intent(in) :: do_anisotropy !< Read anisotropy data (1/0)

      real(dblprec), intent(in) :: A_xc     !< Exchange stiffness constant from rational fit J/m
      real(dblprec), intent(in) :: M_sat    !< Saturation magnetization muB/m^3
      real(dblprec), intent(in) :: cell_vol !< Unit cell volume m^3
      real(dblprec), dimension(NA,6,Nchmax), intent(in) :: anisotropy !< Input data for anisotropies

      ! .. Local variables
      real(dblprec) :: ani_den  !< Anisotropy density MJ/m^3
      real(dblprec) :: DW_width !< Domain wall width m

      ! Calculate micromagnetic parameters dependent on the anisotropy
      if (do_anisotropy.eq.1) then
         ani_den=sum(anisotropy(:,1,1))*mry/cell_vol
         DW_width=sqrt(abs(A_xc)/abs(ani_den))
      else
         ani_den=0.0_dblprec
         DW_width=0.0_dblprec
      endif

      write (ofileno,'(a)') " "
      write (ofileno,'(a)') "**************** MICROMAGNETIC INFORMATION *********************"
      write (ofileno,'(2x,a,G10.3,2x,a)') 'Unit cell volume         :',cell_vol,'m^3'
      write (ofileno,'(2x,a,f10.3,2x,a)') 'Saturation nagnetization :',M_sat*mub/1e6,'MA/m'
      ! If anisotropies are present several parameters can be obtained
      if (do_anisotropy.eq.1) then
         write (ofileno,'(2x,a,f10.3,2x,a)') 'Anisotropy density       :',ani_den/1e6,'MJ/m^3'
         write (ofileno,'(2x,a,f10.3,2x,a)') 'Anisotropy per atom      :',sum(anisotropy(:,1,1)*ry_ev)/NA,'meV'
         write (ofileno,'(2x,a,f10.3,2x,a)') 'Domain wall width        :',DW_width*1e9,'nm'
      endif
      write (ofileno,'(2x,a,i6)') 'ETA MIN :',eta_min
      write (ofileno,'(2x,a,i6)') 'ETA MAX :',eta_max
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)')  " "

   end subroutine print_micro_misc

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing the J0 vector
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine print_J0_vector(NA, J0_matrix)

      use Constants

      implicit none

      ! .. Input variables
      integer, intent(in) :: NA  !< Number of atoms in one cell

      real(dblprec), dimension(NA,NA), intent(in) :: J0_matrix !< Matrix being used to calculate the Tc-MFA

      ! .. Local variables
      integer :: ii,jj

      real(dblprec) :: temp_J0

      write (ofileno,'(a)')  " "
      write (ofileno,'(a)') "**************** J0 VECTOR [meV] *******************************"
      do ii=1, NA
         temp_J0=0.0_dblprec
         do jj=1, NA
            temp_J0=temp_J0+J0_matrix(ii,jj)
         enddo
         write (ofileno,'(2x,i6,2x,f10.3)') ii, temp_J0*ry_ev
      enddo
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)')  " "

   end subroutine print_J0_vector

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing the J0 matrix
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine print_J0(NA, J0_matrix)

      use Constants

      implicit none

      ! .. Input variables
      integer, intent(in) :: NA !< Number of atoms in one cell

      real(dblprec), dimension(NA,NA), intent(in) :: J0_matrix !< Matrix being used to calculate the Tc-MFA

      ! .. Local variables
      integer :: ii,jj

      write (ofileno,'(a)')  " "
      write (ofileno,'(a)') "**************** J0 MATRIX [meV] *******************************"
      do ii=1, NA
         do jj=1, NA
            write (ofileno,'(2x,i6,2x,i6,2x,f10.3)') ii, jj,J0_matrix(ii,jj)*ry_ev
         enddo
      enddo
      write (ofileno,'(a)') "****************************************************************"
      write (ofileno,'(a)')  " "

   end subroutine print_J0

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing site resolved information in random alloys
   !
   !> @author
   !> Lars Bergqvist
   !---------------------------------------------------------------------------------
   subroutine print_random_stiffness(natom,NA,Axc_fit_alloy,Dxc_fit_alloy,Tc_alloy, &
      J0_matrix_alloy)

      use Constants
      implicit none

      ! .. Input variables
      integer, intent(in) :: natom  !< Number of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell

      real(dblprec), dimension(natom), intent(in) :: Axc_fit_alloy !< Exchange stiffness matrix
      real(dblprec), dimension(natom), intent(in) :: Dxc_fit_alloy !< Spin wave  stiffness matrix
      real(dblprec), dimension(natom), intent(in) :: Tc_alloy !< Tc-MFA matrix
      real(dblprec), dimension(NA,NA,natom), intent(in) :: J0_matrix_alloy !< Exchange matrix

      ! .. Local variables
      integer :: ii,jj
      real(dblprec) :: temp_J0

      write (ofileno,'(a)')  " "
      write (ofileno,'(a)') "*** Supercell averaged properties *******************************"
      write (ofileno,'(a)') "*****************************************************************"
      write (ofileno,'(a,f10.4,a)') "Exchange stiffness:   ",sum(Axc_fit_alloy)/natom,"  pJ/m"
      write (ofileno,'(a,f10.4,a)') "Spin wave  stiffness: ",sum(Dxc_fit_alloy)/natom,"  meVÅ^2"
      write (ofileno,'(a,f10.2,a)') "Tc-MFA:               ",sum(Tc_alloy)/natom,"  K"
      write (ofileno,'(a)') "J0 VECTOR [meV]:"
      do ii=1, NA
         temp_J0=0.0_dblprec
         do jj=1, NA
            temp_J0=temp_J0+sum(J0_matrix_alloy(ii,jj,:))/natom
         enddo
         write (ofileno,'(2x,i6,2x,f10.3)') ii, temp_J0*ry_ev
      enddo

      write (ofileno,'(a)')  " "
      write (ofileno,'(a)') "*** Site resolved properties *************************************"
      write (ofileno,'(a)') " atom       A_xc[pJ/m] D_xc[meVÅ^2] Tc-MFA[K]                "
      do ii=1,natom
          write (ofileno,'(2x,i6,2x,f10.4,2x,f10.4,2x,f10.2)') ii,Axc_fit_alloy(ii),Dxc_fit_alloy(ii),Tc_alloy(ii)
      enddo
   end subroutine print_random_stiffness

end module prn_micromagnetic
