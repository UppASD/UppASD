!---------------------------------------------------------------------------
!> @brief
!> Augments module InputData with variables and initialization routines for lattice dynamics
!
!> @author
!> Johan Hellsvik
!---------------------------------------------------------------------------
module LatticeInputData

   use Parameters
   use Profiling

   implicit none

   character(len=1) :: do_ld                            !< Do lattice dynamics ('Y'/'N')
   integer :: initlatt                                  !< Mode of initialization of ionic displacments (1/2/4)
   integer :: lattroteul                                !< Global rotation of ionic displacements
   real(dblprec) , dimension(3) :: lattrotang           !< Euler angles for global rotation of ionic displacements
   character(len=1) :: do_n3                            !< Newton's third law correction of force constant coefficient elements ('Y'/'N')
   character(len=35) :: lattrestartfile                 !< Name of lattice restart file
  !!!SLDTODO Right now only a single iplattdamp value is read. C.f. initial phase for spins
  !!!real(dblprec) :: iplattdamp                          !< Initial phase dissipative damping constant for ionic motion
  real(dblprec), dimension(:), allocatable :: iplattdamp  !< Initial phase dissipative damping constant for ionic motion
   real(dblprec) :: lattdamp                            !< Dissipative damping constant for ionic motion
   character(len=1) :: do_velrsc                        !< Use the Bussi canonical velocity rescaling thermostat ('Y'/'N')
   integer :: velrsc_step                               !< Number of steps between the application of the Bussi thermostat
   real(dblprec) :: velrsc_taut                         !< Relaxation time of the Bussi canonical velocity rescaling thermostat
  character :: do_set_avrgp0                           !< Set average linear momentum to zero
  character :: do_set_avrgu0                           !< Set average displacement to zero
  integer :: imp_max_count                             !< Maximum number of fix-point iterations
  real(dblprec) :: imp_epsilon                         !< Fix-point iteration criteria

   character(len=35) :: phonfile                        !< File name for file with ionic masses and initial u and v
   real(dblprec), &
      dimension(:,:), allocatable :: mion_inp         !< Ionic mass
   real(dblprec), &
      dimension(:,:,:), allocatable :: uvec_inp       !< Initial ionic displacements
   real(dblprec), &
      dimension(:,:,:), allocatable :: vvec_inp       !< Initial ionic velocities


   !LL data
   integer :: do_ll                                     !< Add harmonic lattice forces (LL) term to Hamiltonian (0/1)
   character(len=35) :: llfile                          !< File name for LL data
   real(dblprec), &
      dimension(:,:,:,:), allocatable :: ll_redcoord    !< Neighbour vectors for LL
   real(dblprec), &
      dimension(:,:,:,:,:), allocatable :: ll_inptens !< Coupling tensor for LL
   integer, dimension(:), allocatable :: ll_nn          !< No. shells of neighbours for LL
   integer :: nn_ll_tot                                 !< Calculated number of neighbours with PD interactions
   integer :: max_no_llshells                           !< Actual maximum number of shells for LL interactions


   !LL phonopy data
   integer :: do_ll_phonopy                             !< Read harmonic lattice forces (LL) on phonopy format
   character(len=35) :: ll_phonopyfile                  !< File name for LL-phonopy force constant data
   character(len=35) :: ll_phonopycoordfile             !< File name for LL-phonopy coordinates
   integer :: Natom_phonopy                             !< Number of atoms in phonopy force contant data
   integer :: i0phonopy                                 !< Index of the central atom i
   real(dblprec) :: radius_phonopy                      !< Interaction radius for couplings
   real(dblprec) :: scalefac_phonopy                    !< Scale factor for interactions radius
   integer, dimension(:), allocatable :: atomindex_phonopy  !< Index list of atoms
   real(dblprec), &
      dimension(:,:), allocatable :: ll_coord_phonopy !< Coordinates for 
   real(dblprec), &
      dimension(:,:,:,:), allocatable :: ll_inptens_phonopy !< Coupling tensor for LL phonopy force constants


   !LLL data
   integer :: do_lll                                     !< Add anharmonic lattice forces (LLL) term to Hamiltonian (0/1)
   character(len=35) :: lllfile                          !< File name for LLL data
   real(dblprec), &
      dimension(:,:,:,:), allocatable :: lll_redcoord    !< Neighbour vectors for LLL
   real(dblprec), &
      dimension(:,:,:,:,:), allocatable :: lll_inptens !< Coupling tensor for LLL
   integer, dimension(:), allocatable :: lll_nn          !< No. shells of neighbours for LLL
   integer :: nn_lll_tot                                 !< Calculated number of neighbours with LLL interactions
   integer :: max_no_lllshells                           !< Actual maximum number of shells for LLL interactions


   !LLLL data
   integer :: do_llll                                     !< Add anharmonic lattice forces (LLLL) term to Hamiltonian (0/1)
   character(len=35) :: llllfile                          !< File name for LLLL data
   real(dblprec), &
      dimension(:,:,:,:), allocatable :: llll_redcoord    !< Neighbour vectors for LLLL
   real(dblprec), &
      dimension(:,:,:,:,:), allocatable :: llll_inptens !< Coupling tensor for LLLL
   integer, dimension(:), allocatable :: llll_nn          !< No. shells of neighbours for LLLL
   integer :: nn_llll_tot                                 !< Calculated number of neighbours with LLLL interactions
   integer :: max_no_llllshells                           !< Actual maximum number of shells for LLLL interactions


   !ML data
   integer :: do_ml                                     !< Add spin-lattice coupling (ML) term to Hamiltonian (0/1)
   character(len=35) :: mlfile                          !< File name for ML data
   real(dblprec), &
      dimension(:,:,:,:), allocatable :: ml_redcoord    !< Neighbour vectors for ML
   real(dblprec), &
      dimension(:,:,:,:,:), allocatable :: ml_inptens !< Coupling tensor for ML
   integer, dimension(:), allocatable :: ml_nn          !< No. shells of neighbours for ML
   integer :: nn_ml_tot                                 !< Calculated number of neighbours with ML interactions
   integer :: max_no_mlshells                           !< Actual maximum number of shells for ML interactions


   !MML data
   integer :: do_mml                                     !< Add spin-lattice coupling (MML) term to Hamiltonian (0/1)
   logical:: mml_diag                                    !< Use diagonal part of MML tensor only (T/F)
   logical :: mml_ene_opt                                !< Optimize calculation of GS energy (T/F) 
   character(len=35) :: mmlfile                          !< File name for MML data
   real(dblprec), &
      dimension(:,:,:,:), allocatable :: mml_redcoord  !< Neighbour vectors for MML
   real(dblprec), &
      dimension(:,:,:,:,:), allocatable :: mml_inptens !< Coupling tensor for MML
   integer, dimension(:), allocatable :: mml_nn          !< No. shells of neighbours for MML
   integer :: nn_mml_tot                                 !< Calculated number of neighbours with MML interactions
   integer :: max_no_mmlshells                           !< Actual maximum number of shells for MML interactions
   integer, dimension(:,:), allocatable :: mml_invsym    !< Inversion symmetry of the coupling (1/-1)
   real(dblprec) :: mml_scale                            !< Manual scaling of mml couplings

   !MMLL data
   integer :: do_mmll                                     !< Add spin-lattice coupling (MMLL) term to Hamiltonian (0/1)
   character(len=35) :: mmllfile                          !< File name for MMLL data
   real(dblprec), &
      dimension(:,:,:,:), allocatable :: mmll_redcoord    !< Neighbour vectors for MMLL
   real(dblprec), &
      dimension(:,:,:,:,:), allocatable :: mmll_inptens !< Coupling tensor for MMLL
   integer, dimension(:), allocatable :: mmll_nn          !< No. shells of neighbours for MMLL
   integer :: nn_mmll_tot                                 !< Calculated number of neighbours with MMLL interactions
   integer :: max_no_mmllshells                           !< Actual maximum number of shells for MMLL interactions


   ! Adiabatic phonon Spectra calculation flags
   character(LEN=1) :: do_phonspec              !< Calculate phonon spectra (N/Y)
   character(LEN=1) :: do_phondos               !< Calculate phonon density of states (N/Y/F)
   character(len=35) :: phondosfile             !< Phonon DOS file
   real(dblprec)    :: phondos_sigma            !< Frequency broadening of phonon spectra DOS (in meV)
   integer          :: phondos_freq             !< Number of frequencies of phonon spectra DOS



contains


   !> Sets default values to input variables
   subroutine set_lattinput_defaults()
      !
      implicit none

      do_ld = 'N'
      phonfile = 'phonfile'
      initlatt = 4
      do_n3 ='N'
      lattrestartfile = 'lattrestart.dat'
      lattdamp = 0_dblprec
      if(allocated(iplattdamp)) iplattdamp = 0_dblprec
      do_velrsc = 'N'
      velrsc_step = 100
      velrsc_taut = 0.5_dblprec
    do_set_avrgp0 = 'N';
    do_set_avrgu0 = 'N';
    imp_epsilon = 1e-10;
    imp_max_count = 20;

      !LL data
      llfile = 'llfile'
      do_ll  = 0

      !LL phonopy data
      do_ll_phonopy = 0
      ll_phonopyfile = 'll_phonopyfile'
      ll_phonopycoordfile = 'll_phonopycoordfile'
      Natom_phonopy = 0
      i0phonopy = 0
      radius_phonopy = 1.0_dblprec
      scalefac_phonopy = 1.0_dblprec

      !LLL data
      lllfile = 'lllfile'
      do_lll  = 0

      !LLLL data
      llllfile = 'llllfile'
      do_llll  = 0

      !ML data
      mlfile = 'mlfile'
      do_ml  = 0

      !MML data
      mmlfile = 'mmlfile'
      do_mml  = 0
      mml_diag= .false.
      mml_ene_opt= .false.
      mml_scale=1.0_dblprec

      !MMLL data
      mmllfile = 'mmllfile'
      do_mmll  = 0

      do_phonspec = 'N'
      do_phondos = 'N'

   end subroutine set_lattinput_defaults

   !> Allocate arrays for input data for lattice Hamiltonian
   subroutine allocate_latthamiltonianinput(no_shells, flag) !NA, limit_no_shells, Nchmax, flag)
      use Parameters
      use Profiling
      use InputData, only : NT
      !use LatticeInputData
      implicit none

      integer, intent(in),optional :: no_shells !< Parameter limiting number of exchange coupling shells
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then

         allocate(ll_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(ll_nn))*kind(ll_nn),'ll_nn','allocate_latthamiltonianinput')
         ll_nn=0

         allocate(lll_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(lll_nn))*kind(lll_nn),'lll_nn','allocate_latthamiltonianinput')
         lll_nn=0

         allocate(llll_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(llll_nn))*kind(llll_nn),'llll_nn','allocate_latthamiltonianinput')
         llll_nn=0

         allocate(ml_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(ml_nn))*kind(ml_nn),'ml_nn','allocate_latthamiltonianinput')
         ml_nn=0

         allocate(mml_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(mml_nn))*kind(mml_nn),'mml_nn','allocate_latthamiltonianinput')
         mml_nn=0

         allocate(mmll_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(mmll_nn))*kind(mmll_nn),'mmll_nn','allocate_latthamiltonianinput')
         mmll_nn=0

      else

         i_all=-product(shape(ll_nn))*kind(ll_nn)
         deallocate(ll_nn,stat=i_stat)
         call memocc(i_stat,i_all,'ll_nn','allocate_latthamiltonianinput')
         if (allocated(ll_redcoord)) then
            i_all=-product(shape(ll_redcoord))*kind(ll_redcoord)
            deallocate(ll_redcoord,stat=i_stat)
            call memocc(i_stat,i_all,'ll_redcoord','allocate_latthamiltonianinput')
         endif
         if (allocated(ll_inptens)) then
            i_all=-product(shape(ll_inptens))*kind(ll_inptens)
            deallocate(ll_inptens,stat=i_stat)
            call memocc(i_stat,i_all,'ll_inptens','allocate_latthamiltonianinput')
         endif

         i_all=-product(shape(lll_nn))*kind(lll_nn)
         deallocate(lll_nn,stat=i_stat)
         call memocc(i_stat,i_all,'lll_nn','allocate_latthamiltonianinput')
         if (allocated(lll_redcoord)) then
            i_all=-product(shape(lll_redcoord))*kind(lll_redcoord)
            deallocate(lll_redcoord,stat=i_stat)
            call memocc(i_stat,i_all,'lll_redcoord','allocate_latthamiltonianinput')
         endif
         if (allocated(lll_inptens)) then
            i_all=-product(shape(lll_inptens))*kind(lll_inptens)
            deallocate(lll_inptens,stat=i_stat)
            call memocc(i_stat,i_all,'lll_inptens','allocate_latthamiltonianinput')
         endif

         i_all=-product(shape(llll_nn))*kind(llll_nn)
         deallocate(llll_nn,stat=i_stat)
         call memocc(i_stat,i_all,'llll_nn','allocate_latthamiltonianinput')
         if (allocated(llll_redcoord)) then
            i_all=-product(shape(llll_redcoord))*kind(llll_redcoord)
            deallocate(llll_redcoord,stat=i_stat)
            call memocc(i_stat,i_all,'llll_redcoord','allocate_latthamiltonianinput')
         endif
         if (allocated(llll_inptens)) then
            i_all=-product(shape(llll_inptens))*kind(llll_inptens)
            deallocate(llll_inptens,stat=i_stat)
            call memocc(i_stat,i_all,'llll_inptens','allocate_latthamiltonianinput')
         endif

         i_all=-product(shape(ml_nn))*kind(ml_nn)
         deallocate(ml_nn,stat=i_stat)
         call memocc(i_stat,i_all,'ml_nn','allocate_latthamiltonianinput')
         if (allocated(ml_redcoord)) then
            i_all=-product(shape(ml_redcoord))*kind(ml_redcoord)
            deallocate(ml_redcoord,stat=i_stat)
            call memocc(i_stat,i_all,'ml_redcoord','allocate_latthamiltonianinput')
         endif
         if (allocated(ml_inptens)) then
            i_all=-product(shape(ml_inptens))*kind(ml_inptens)
            deallocate(ml_inptens,stat=i_stat)
            call memocc(i_stat,i_all,'ml_inptens','allocate_latthamiltonianinput')
         endif

         i_all=-product(shape(mml_nn))*kind(mml_nn)
         deallocate(mml_nn,stat=i_stat)
         call memocc(i_stat,i_all,'mml_nn','allocate_latthamiltonianinput')
         if (allocated(mml_redcoord)) then
            i_all=-product(shape(mml_redcoord))*kind(mml_redcoord)
            deallocate(mml_redcoord,stat=i_stat)
            call memocc(i_stat,i_all,'mml_redcoord','allocate_latthamiltonianinput')
         endif
         if (allocated(mml_inptens)) then
            i_all=-product(shape(mml_inptens))*kind(mml_inptens)
            deallocate(mml_inptens,stat=i_stat)
            call memocc(i_stat,i_all,'mml_inptens','allocate_latthamiltonianinput')
         endif

         i_all=-product(shape(mmll_nn))*kind(mmll_nn)
         deallocate(mmll_nn,stat=i_stat)
         call memocc(i_stat,i_all,'mmll_nn','allocate_latthamiltonianinput')
         if (allocated(mmll_redcoord)) then
            i_all=-product(shape(mmll_redcoord))*kind(mmll_redcoord)
            deallocate(mmll_redcoord,stat=i_stat)
            call memocc(i_stat,i_all,'mmll_redcoord','allocate_latthamiltonianinput')
         endif
         if (allocated(mmll_inptens)) then
            i_all=-product(shape(mmll_inptens))*kind(mmll_inptens)
            deallocate(mmll_inptens,stat=i_stat)
            call memocc(i_stat,i_all,'mmll_inptens','allocate_latthamiltonianinput')
         endif

      end if

   end subroutine allocate_latthamiltonianinput



end module LatticeInputData
