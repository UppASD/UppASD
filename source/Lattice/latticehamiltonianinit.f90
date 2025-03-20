!> Routines for mounting the lattice Hamiltonian
!Augments module HamiltonianInit
module LatticeHamiltonianInit
   use Parameters
   use Profiling

!!! TEMPORARY
   use InputData, only : do_hoc_debug
   use LatticeInputData, only : i0phonopy, radius_phonopy, Natom_phonopy, mml_scale, mml_diag
!!! TEMPORARY

   implicit none

   integer, dimension(:,:,:,:), allocatable :: nme !< Neighbour map, elements
   integer, dimension(:,:,:), allocatable :: nm !< Neighbour map
   integer, dimension(:,:), allocatable :: nmdim !< Dimension of neighbour map

   private
   public :: setup_latticehamiltonian
   public :: setup_neighbour_latticehamiltonian, deallocate_nme



contains


   !>  Main routine for mounting the lattice Hamiltonian
   subroutine setup_latticehamiltonian(simid, Natom, &
      NT, NA, N1, N2, N3, &
      C1, C2, C3, BC1, BC2, BC3, &
      atype, anumb, coord, alat, Bas, do_prnstruct, do_sortcoup, sym, do_n3, &
      do_ll, ll_nn, nn_ll_tot, max_no_llshells, ll_listsize, ll_list, ll_tens, ll_redcoord, ll_inptens, &
      do_lll, lll_nn, nn_lll_tot, max_no_lllshells, lll_listsize, lll_list, lll_tens, lll_redcoord, lll_inptens, &
      do_llll, llll_nn, nn_llll_tot, max_no_llllshells, llll_listsize, llll_list, llll_tens, llll_redcoord, llll_inptens, &
      do_ml, ml_nn, nn_ml_tot, max_no_mlshells, ml_listsize, ml_list, ml_tens, ml_redcoord, ml_inptens, &
      do_mml, mml_nn, nn_mml_tot, max_no_mmlshells, mml_listsize, mml_list, mml_tens, mml_redcoord, mml_inptens, mml_invsym, &
      do_mmll, mmll_nn, nn_mmll_tot, max_no_mmllshells, mmll_listsize, mmll_list, mmll_tens, mmll_redcoord, mmll_inptens, &
      do_ralloy, acellnumb, acellnumbrev, achtype, nchmax, &
      natom_full, ammom_inp, atype_ch, asite_ch, achem_ch, &
      lm_listsize, lm_list, lm_tens, &
      lmm_listsize, lmm_list, lmm_tens, &
      llmm_listsize, llmm_list, llmm_tens, &
      do_ll_phonopy, Natom_phonopy, atomindex_phonopy, ll_inptens_phonopy )

      use NeighbourMap, only : setup_nm, setup_nm_nelem
      use LatticeInputData, only : allocate_latthamiltonianinput
      use LatticeHamiltonianData, only : allocate_llhamiltoniandata, allocate_lllhamiltoniandata, &
         allocate_llllhamiltoniandata, allocate_mlhamiltoniandata, allocate_mmlhamiltoniandata, allocate_mmllhamiltoniandata, &
         mml_tens_diag, lmm_tens_diag
      use LatticePrintHamiltonian
      Use Constants

      !.. Implicit declarations
      implicit none

      character(len=8) :: simid                                          !< Name of simulation
      integer, intent(inout) :: Natom                                    !< Number of atoms in system
      integer, intent(in) :: NT                                          !< Number of types of atoms
      integer, intent(in) :: NA                                          !< Number of atoms in one cell
      integer, intent(in) :: N1                                          !< Number of cell repetitions in x direction
      integer, intent(in) :: N2                                          !< Number of cell repetitions in y direction
      integer, intent(in) :: N3                                          !< Number of cell repetitions in z direction

      real(dblprec), dimension(3), intent(in) :: C1                      !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2                      !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3                      !< Third lattice vector
      character(len=1), intent(in) :: BC1                                !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2                                !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3                                !< Boundary conditions in z-direction

      integer, dimension(Natom), intent(inout) :: atype                  !< Type of atom
      integer, dimension(Natom), intent(inout) :: anumb                  !< Atom number in cell
      real(dblprec), dimension(:,:), allocatable, intent(inout) :: coord !< Coordinates of atoms
      real(dblprec) :: alat                                              !< Lattice parameter
      real(dblprec), dimension(3,NA), intent(inout) :: Bas               !< Coordinates for basis atoms
      integer, intent(in) :: do_prnstruct                                !< Print Hamiltonian information (0/1)
      character, intent(in) :: do_sortcoup                               !< Sort the entries of ncoup arrays (Y/N)
      integer, intent(in) :: sym                                         !< Symmetry of system (0-3)
      character(len=1) :: do_n3           !< Newton´s third law correction of force constant coefficient elements ('Y'/'N')

      integer, intent(in) :: do_ll                                       !< Add harmonic (LL) term to the lattice Hamiltonian (0/1)
      integer, dimension(:), intent(in) :: ll_nn                         !< Number of neighbour shells for LL
      integer, intent(inout) :: nn_ll_tot                                !< Calculated number of neighbours with LL interactions
      integer, intent(in) :: max_no_llshells                             !< Limit of number of shells for LL interactions
      ! SLKTODO Ensure proper use of attribute allocatable for these variables
      integer, dimension(:), allocatable, intent(inout) :: ll_listsize          !< Size of neighbour list for LL
      integer, dimension(:,:,:), allocatable, intent(inout) :: ll_list            !< List of neighbours for LL
      real(dblprec), dimension(:,:,:,:), allocatable, intent(inout) :: ll_tens    !< LL coupling tensor
      real(dblprec), dimension(:,:,:,:), intent(in) :: ll_redcoord         !< Coordinates for LL exchange couplings
      real(dblprec), dimension(9,NT,max_no_llshells,NT,NT), intent(in) :: ll_inptens      !< LL couplings
      !real(dblprec), dimension(:,:,:,:,:), intent(in) :: ll_inptens      !< LL couplings

      integer, intent(in) :: do_lll                                      !< Add anharmonic (LLL) term to the lattice Hamiltonian (0/1)
      integer, dimension(:), intent(in) :: lll_nn                        !< Number of neighbour shells for LLL
      integer, intent(inout) :: nn_lll_tot                               !< Calculated number of neighbours with LLL interactions
      integer, intent(in) :: max_no_lllshells                            !< Limit of number of shells for LLL interactions
      integer, dimension(:), allocatable, intent(inout) :: lll_listsize          !< Size of neighbour list for LLL
      integer, dimension(:,:,:), allocatable, intent(inout) :: lll_list            !< List of neighbours for LLL
      real(dblprec), dimension(:,:,:,:,:), allocatable, intent(inout) :: lll_tens    !< LLL coupling tensor
      real(dblprec), dimension(:,:,:,:), intent(in) :: lll_redcoord      !< Coordinates for LLL exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: lll_inptens     !< LLL couplings

      integer, intent(in) :: do_llll                                      !< Add anharmonic (LLLL) term to the lattice Hamiltonian (0/1)
      integer, dimension(:), intent(in) :: llll_nn                        !< Number of neighbour shells for LLLL
      integer, intent(inout) :: nn_llll_tot                               !< Calculated number of neighbours with LLLL interactions
      integer, intent(in) :: max_no_llllshells                            !< Limit of number of shells for LLLL interactions
      integer, dimension(:), allocatable, intent(inout) :: llll_listsize          !< Size of neighbour list for LLLL
      integer, dimension(:,:,:), allocatable, intent(inout) :: llll_list            !< List of neighbours for LLLL
      real(dblprec), dimension(:,:,:,:,:,:), allocatable, intent(inout) :: llll_tens    !< LLLL coupling tensor
      real(dblprec), dimension(:,:,:,:), intent(in) :: llll_redcoord      !< Coordinates for LLLL exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: llll_inptens     !< LLLL couplings

      integer, intent(in) :: do_ml                                      !< Add anharmonic (ML) term to the lattice Hamiltonian (0/1)
      integer, dimension(:), intent(in) :: ml_nn                        !< Number of neighbour shells for ML
      integer, intent(inout) :: nn_ml_tot                               !< Calculated number of neighbours with ML interactions
      integer, intent(in) :: max_no_mlshells                            !< Limit of number of shells for ML interactions
      integer, dimension(:), allocatable, intent(inout) :: ml_listsize          !< Size of neighbour list for ML
      integer, dimension(:,:,:), allocatable, intent(inout) :: ml_list            !< List of neighbours for ML
      real(dblprec), dimension(:,:,:,:), allocatable, intent(inout) :: ml_tens    !< ML coupling tensor
      real(dblprec), dimension(:,:,:,:), intent(in) :: ml_redcoord      !< Coordinates for ML exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: ml_inptens     !< ML couplings

      integer, intent(in) :: do_mml                                      !< Add anharmonic (MML) term to the lattice Hamiltonian (0/1)
      integer, dimension(:), intent(in) :: mml_nn                        !< Number of neighbour shells for MML
      integer, intent(inout) :: nn_mml_tot                               !< Calculated number of neighbours with MML interactions
      integer, intent(in) :: max_no_mmlshells                            !< Limit of number of shells for MML interactions
      integer, dimension(:), allocatable, intent(inout) :: mml_listsize          !< Size of neighbour list for MML
      integer, dimension(:,:,:), allocatable, intent(inout) :: mml_list            !< List of neighbours for MML
      real(dblprec), dimension(:,:,:,:,:), allocatable, intent(inout) :: mml_tens    !< MML coupling tensor
      real(dblprec), dimension(:,:,:,:), intent(in) :: mml_redcoord      !< Coordinates for MML exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: mml_inptens     !< MML couplings
      integer, dimension(:,:) :: mml_invsym                              !< Inversion symmetry of the coupling (1/-1)

      integer, intent(in) :: do_mmll                                      !< Add anharmonic (MMLL) term to the lattice Hamiltonian (0/1)
      integer, dimension(:), intent(in) :: mmll_nn                        !< Number of neighbour shells for MMLL
      integer, intent(inout) :: nn_mmll_tot                               !< Calculated number of neighbours with MMLL interactions
      integer, intent(in) :: max_no_mmllshells                            !< Limit of number of shells for MMLL interactions
      integer, dimension(:), allocatable, intent(inout) :: mmll_listsize          !< Size of neighbour list for MMLL
      integer, dimension(:,:,:), allocatable, intent(inout) :: mmll_list            !< List of neighbours for MMLL
      real(dblprec), dimension(:,:,:,:,:,:), allocatable, intent(inout) :: mmll_tens    !< MMLL coupling tensor
      real(dblprec), dimension(:,:,:,:), intent(in) :: mmll_redcoord      !< Coordinates for MMLL exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: mmll_inptens     !< MMLL couplings

      integer, intent(in) :: do_ralloy                                   !< Random alloy simulation (0/1)
      integer, dimension(:), intent(inout) :: acellnumb                  !< List for translating atom no. in full cell to actual cell
      integer, dimension(:), intent(inout) :: acellnumbrev               !< List for translating atom no. in actual cell to full cell
      integer, dimension(:), intent(inout) :: achtype                    !< Chemical type of atoms (full list)
      integer, intent(in) :: nchmax                                      !< Max number of chemical components on each site in cell

      integer, intent(in) :: Natom_full                                  !< Number of atoms for full system (=Natom if not dilute)
      real(dblprec), dimension(:,:,:), intent(inout) :: ammom_inp        !< Magnetic moment directions from input (for alloys)
      integer, dimension(:), intent(inout) :: atype_ch                   !< Actual type of atom for dilute system
      integer, dimension(:), intent(inout) :: asite_ch                   !< Actual site of atom for dilute system
      integer, dimension(:), intent(inout) :: achem_ch                   !< Chemical type of atoms (reduced list)

      integer, dimension(:), allocatable, intent(inout) :: lm_listsize          !< Size of neighbour list for LM
      integer, dimension(:,:,:), allocatable, intent(inout) :: lm_list            !< List of neighbours for LM
      real(dblprec), dimension(:,:,:,:), allocatable, intent(inout) :: lm_tens    !< LM coupling tensor

      integer, dimension(:), allocatable, intent(inout) :: lmm_listsize          !< Size of neighbour list for LMM
      integer, dimension(:,:,:), allocatable, intent(inout) :: lmm_list            !< List of neighbours for LMM
      real(dblprec), dimension(:,:,:,:,:), allocatable, intent(inout) :: lmm_tens    !< LMM coupling tensor

      integer, dimension(:), allocatable, intent(inout) :: llmm_listsize          !< Size of neighbour list for LLMM
      integer, dimension(:,:,:), allocatable, intent(inout) :: llmm_list            !< List of neighbours for LLMM
      real(dblprec), dimension(:,:,:,:,:,:), allocatable, intent(inout) :: llmm_tens    !< LLMM coupling tensor

      integer, intent(in) :: do_ll_phonopy                             !< Read harmonic lattice forces (LL) on phonopy format
      integer, intent(in) :: Natom_phonopy                             !< Number of atoms in phonopy force contant data
      integer, dimension(:), allocatable, intent(in) :: atomindex_phonopy  !< Index list of atoms
      real(dblprec), &
         dimension(:,:,:,:), allocatable, intent(in) :: ll_inptens_phonopy !< Coupling tensor for LL phonopy force constants

      ! Local variables
      real(dblprec), dimension(9,NT,max_no_llshells,NT,NT,48) :: ll_symtens
      real(dblprec), dimension(27,NT,max_no_lllshells,NT,NT,48) :: lll_symtens
      real(dblprec), dimension(81,NT,max_no_llllshells,NT,NT,48) :: llll_symtens
      real(dblprec), dimension(9,NT,max_no_mlshells,NT,NT,48) :: ml_symtens
      real(dblprec), dimension(27,NT,max_no_mmlshells,NT,NT,48) :: mml_symtens
      real(dblprec), dimension(81,NT,max_no_mmllshells,NT,NT,48) :: mmll_symtens

      integer, dimension(48,max_no_llshells,na) :: ll_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      integer, dimension(48,max_no_lllshells,na) :: lll_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      integer, dimension(48,max_no_llllshells,na) :: llll_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      integer, dimension(48,max_no_mlshells,na) :: ml_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      integer, dimension(48,max_no_mmlshells,na) :: mml_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      integer, dimension(48,max_no_mmllshells,na) :: mmll_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      !integer, dimension(48,max_no_shells,na) :: ll_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      !integer, dimension(max_no_equiv,max_no_shells,na) :: ll_symind  !< Indices for elements of the symmetry degenerate coupling tensor

      integer :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      real(dblprec) :: fc_unitconv_phonopy
      integer :: i, j, ia, ja

      integer :: k
      real(dblprec), dimension(3,3,3) :: tempmmlcoup

      if(do_ll==1) then

         ! Allocate and mount LL Hamiltonian
         ! Makes use of point symmetry operations for the coupling tensor
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for LL interaction'
         call setup_nm_nelem(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            nn_ll_tot, max_no_llshells, max_no_equiv, sym, &
            ll_nn, ll_redcoord, nme, nmdim, 1, &
            do_ralloy, Natom_full, acellnumb, atype_ch, &
            !Nchmax, 9, .true., 2, 1, 1, ll_inptens, ll_symtens, ll_symind)
            Nchmax, 9, .false., 2, 0, 0, ll_inptens, ll_symtens, ll_symind)
            !Nchmax, hdim, do_tens_sym_in, couptensrank, invsym, timesym, couptens, fullcouptens, nm_cell_symind)
            !Nchmax, 9, .true., 2, 1, 1, ll_inptens, ll_symtens)
            !Nchmax, hdim, do_tens_sym, couptensrank, invsym, timesym, couptens, fullcouptens)
         write(*,'(a)') ' done'

         
         if(do_ll_phonopy==1) nn_ll_tot = Natom_phonopy
         !if(do_ll_phonopy==1) nn_ll_tot = Natom

         call allocate_llhamiltoniandata(Natom, nn_ll_tot, 1)

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount LL Hamiltonian'
         call setup_neighbour_latticehamiltonian(Natom, NT, NA, anumb, atype, &
            nn_ll_tot, max_no_equiv, max_no_llshells, 1, &
            9, 1, ll_listsize, ll_nn, ll_list, ll_tens, nme, nmdim, ll_symtens, ll_symind, do_sortcoup, do_n3, &
            !9, 1, ll_listsize, ll_nn, ll_list, ll_tens, nme, nmdim, ll_symtens, do_sortcoup, do_n3, &
            do_ralloy, Natom_full, Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, &
            0, 0)
         write(*,'(a)') ' done'

         ! Print LL interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print LL coordination shells"
            call prn_ll_shells(simid, Natom, NT, NA, N1, N2, N3, atype, &
               max_no_llshells, max_no_equiv, ll_redcoord, ll_nn, nme, nmdim, &
               do_ralloy, Natom_full, Nchmax, acellnumb, acellnumbrev, achtype)
            write(*,'(2x,a)',advance='no') "Print LL interactions"
            call prn_llcoup(Natom, nn_ll_tot, ll_listsize, ll_list, ll_tens, simid)
            write(*,'(a)') ' done'
         end if

         ! Deallocate the large neighbour map.
         call deallocate_nme()

      end if


      ! This overrides the LL couplings that are read and setup by do_ll! 
      if(do_ll_phonopy==1 .and. i0phonopy==0) then

         ! Factor for unit conversion from [ll_inptens_phonopy] = eV / angstrom**2
         ! to [ll_tens] = mRyd / angstrom**2
         fc_unitconv_phonopy = 1000 / Ry_ev
         write(*,*) 'fc_unitconv_phonopy ', fc_unitconv_phonopy
         !fc_unitconv_phonopy = 1000 * mry / (amu * ry_ev)
         !write(*,*) 'fc_unitconv_phonopy ', fc_unitconv_phonopy
         !fc_unitconv_phonopy = ev / amu
         !write(*,*) 'fc_unitconv_phonopy ', fc_unitconv_phonopy

         ! Sorts the atoms according to atomindex_phonopy
         ! and mounts the couplings
         ! SLDTODO Remove the atomindex_phonopy construct altogether
         do i=1,Natom
            ia = atomindex_phonopy(i)
            ll_listsize(ia) = Natom
            do j=1,Natom
               ja = atomindex_phonopy(j)
               if(do_hoc_debug==1) then
                  write(*,*) 'i, j, ia, ja ', i, j, ia, ja
               end if
               ll_list(1,j,i) = j 
               ll_tens(1:3,1:3,j,i) = ll_inptens_phonopy(1:3,1:3,ja,ia) * fc_unitconv_phonopy
            end do
         end do

         ! Print LL phonopy interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print LL phonopy interactions"
            call prn_llphonopycoup(Natom, nn_ll_tot, ll_listsize, ll_list, ll_tens, simid)
            write(*,'(a)') ' done'
         end if

      end if


      if(do_lll==1) then

         ! Allocate and mount LLL Hamiltonian
         ! Makes use of point symmetry operations for the coupling tensor
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for LLL interaction'
         call setup_nm_nelem(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            nn_lll_tot, max_no_lllshells, max_no_equiv, sym, lll_nn, lll_redcoord, nme, nmdim, 2, &
            do_ralloy, Natom_full, acellnumb, atype_ch, &
            Nchmax, 27, .true., 3, 1, 1, lll_inptens, lll_symtens, lll_symind)
            !Nchmax, 27, .true., 3, -1, 1, lll_inptens, lll_symtens)
            !Nchmax, hdim, do_tens_sym, couptensrank, invsym, timesym, couptens, fullcouptens)
         write(*,'(a)') ' done'

         write(*,*) 'nn_lll_tot ', nn_lll_tot
         call allocate_lllhamiltoniandata(Natom, nn_lll_tot, 1)

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount LLL Hamiltonian'
         call setup_neighbour_latticehamiltonian(Natom, NT, NA, anumb, atype, &
            nn_lll_tot, max_no_equiv, max_no_lllshells, 2, &
            27, 1, lll_listsize, lll_nn, lll_list, lll_tens, nme, nmdim, lll_symtens, lll_symind, do_sortcoup, 'N', &
            do_ralloy, Natom_full, Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, &
            0, 0)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nme()

         ! Print LLL interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print LLL interactions"
            call prn_lllcoup(Natom, nn_lll_tot, lll_listsize, lll_list, lll_tens, simid)
            write(*,'(a)') ' done'
         end if

      end if


      if(do_llll==1) then

         ! Allocate and mount LLLL Hamiltonian
         ! Makes use of point symmetry operations for the coupling tensor
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for LLLL interaction'
         call setup_nm_nelem(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            nn_llll_tot, max_no_llllshells, max_no_equiv, sym, llll_nn, llll_redcoord, nme, nmdim, 3, &
            do_ralloy, Natom_full, acellnumb, atype_ch, &
            Nchmax, 81, .true., 4, 1, 1, llll_inptens, llll_symtens, llll_symind)
            !Nchmax, hdim, do_tens_sym, couptensrank, invsym, timesym, couptens, fullcouptens)
         write(*,'(a)') ' done'

         write(*,*) 'nn_llll_tot ', nn_llll_tot
         call allocate_llllhamiltoniandata(Natom, nn_llll_tot, 1)

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount LLLL Hamiltonian'
         call setup_neighbour_latticehamiltonian(Natom, NT, NA, anumb, atype, &
            nn_llll_tot, max_no_equiv, max_no_llllshells, 3, &
            81, 1, llll_listsize, llll_nn, llll_list, llll_tens, nme, nmdim, llll_symtens, llll_symind, do_sortcoup, 'N', &
            do_ralloy, Natom_full, Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, &
            0, 0)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nme()

         ! Print LLLL interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print LLLL interactions"
            call prn_llllcoup(Natom, nn_llll_tot, llll_listsize, llll_list, llll_tens, simid)
            write(*,'(a)') ' done'
         end if

      end if


      if(do_ml==1) then

         ! Allocate and mount ML Hamiltonian
         ! Makes use of point symmetry operations for the coupling tensor
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for ML interaction'
         call setup_nm_nelem(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            nn_ml_tot, max_no_mlshells, max_no_equiv, sym, ml_nn, ml_redcoord, nme, nmdim, 1, &
            do_ralloy, Natom_full, acellnumb, atype_ch, &
            Nchmax, 9, .true., 2, 1, 1, ml_inptens, ml_symtens, ml_symind)
            !Nchmax, 9, .true., 2, -1, -1, ml_inptens, ml_symtens)
            !Nchmax, hdim, do_tens_sym, couptensrank, invsym, timesym, couptens, fullcouptens)
         write(*,'(a)') ' done'

         write(*,*) 'nn_ml_tot ', nn_ml_tot
         call allocate_mlhamiltoniandata(Natom, nn_ml_tot, 1)

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount ML Hamiltonian'
         call setup_neighbour_latticehamiltonian(Natom, NT, NA, anumb, atype, &
            nn_ml_tot, max_no_equiv, max_no_mlshells, 1, &
            9, 1, ml_listsize, ml_nn, ml_list, ml_tens, nme, nmdim, ml_symtens, ml_symind, do_sortcoup, 'N', &
            do_ralloy, Natom_full, Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, &
            1, 0)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nme()

         if(do_hoc_debug==1) then
            write(*,*) 'ml_tens'
            write(*,*) ml_tens
         end if

         call setup_lm_arrays(Natom, max_no_mlshells, nn_ml_tot, ml_listsize, ml_list, ml_tens, lm_listsize, lm_list, lm_tens)

         ! Print ML interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            !write(*,'(2x,a)',advance='no') "Print ML coordination shells"
            !call prn_ml_shells(simid, Natom, NT, NA, N1, N2, N3, atype, &
            !     max_no_llshells, max_no_equiv, ll_redcoord, ll_nn, nme, nmdim, &
            !     do_ralloy, Natom_full, Nchmax, acellnumb, acellnumbrev, achtype)
            write(*,'(2x,a)',advance='no') "Print ML interactions"
            call prn_mlcoup(Natom, nn_ml_tot, ml_listsize, ml_list, ml_tens, simid)
            write(*,'(a)') ' done'
         end if

         ! Print LM interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print LM interactions"
            call prn_lmcoup(Natom, nn_ml_tot, lm_listsize, lm_list, lm_tens, simid)
            write(*,'(a)') ' done'
         end if

      end if


      if(do_mml>=1) then

         ! Allocate and mount MML Hamiltonian
         ! Makes use of point symmetry operations for the coupling tensor
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for MML interaction'
         call setup_nm_nelem(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            nn_mml_tot, max_no_mmlshells, max_no_equiv, sym, mml_nn, mml_redcoord, nme, nmdim, 2, &
            do_ralloy, Natom_full, acellnumb, atype_ch, &
            Nchmax, 27, .true., 3, 1, 1, mml_inptens, mml_symtens, mml_symind)
            !Nchmax, 27, .true., 3, -1, 1, mml_inptens, mml_symtens)
            !Nchmax, hdim, do_tens_sym, couptensrank, invsym, timesym, couptens, fullcouptens)
         write(*,'(a)') ' done'

         write(*,*) 'nn_mml_tot ', nn_mml_tot
         call allocate_mmlhamiltoniandata(Natom, nn_mml_tot, 1)

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount MML Hamiltonian'
         call setup_neighbour_latticehamiltonian(Natom, NT, NA, anumb, atype, &
            nn_mml_tot, max_no_equiv, max_no_mmlshells, 2, &
            27, 1, mml_listsize, mml_nn, mml_list, mml_tens, nme, nmdim, mml_symtens, mml_symind, do_sortcoup, 'N', &
            do_ralloy, Natom_full, Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, &
            1, 1)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nme()

         !if(do_hoc_debug==1) then
         !   write(*,*) 'mml_tens'
         !   write(*,*) mml_tens
         !end if

         !      mml_tens=nint(1.0e3*mml_tens)/1.0e3
         ! Calculates the correct force constant coefficient element 
         ! \Psi_{ii} = -\sum_{j\neq i} Psi_{ij}
         ! as required by Newton´s third law.
         ! SLDTODO Applicable for higher order couplings?
         if(do_n3 == 'Z') then
            !write(*,*) 'Start do_N3'
            do i=1,Natom
               tempmmlcoup = 0.0_dblprec
               !write(*,*) '-----------------'
               do j=1,mml_listsize(i)
                  write(*,*) 'i ', i, ' j ', j
                  if( mml_list(2,j,i) == i ) then
                     !if( mml_list(1,j,i) == i ) then
                     k = j
                     !write(*,*) 'hit ', i, mml_list(1,j,i), j
                  else
                     tempmmlcoup(1:3,1:3,1:3) = tempmmlcoup + mml_tens(1:3,1:3,1:3,j,i)
                     !write(*,*) 'hop ', i, mml_list(1,j,i), j
                  end if
               end do
               !mml_tens(1:3,1:3,1:3,k,i) = -tempmmlcoup(1:3,1:3,1:3)
               mml_tens(1:3,1:3,1:3,k,i) = -0.5*tempmmlcoup(1:3,1:3,1:3)
               !mml_tens(:,k,i) = -0.5*tempmmlcoup
               !mml_tens(:,k,i) = -tempmmlcoup
            end do
         end if

         ! Rescale MML couplings if needed
         if(mml_scale.ne.1.0_dblprec) mml_tens=mml_scale*mml_tens

         ! Setup reduced tensor for diagonal couplings
         if(mml_diag) then
            do i=1,Natom
               do j=1,mml_listsize(i)
                  mml_tens_diag(:,j,i)=mml_tens(1,1,:,j,i) 
               end do
            end do
         end if

         ! Print MML interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print MML interactions"
            call prn_mmlcoup(Natom, nn_mml_tot, mml_listsize, mml_list, mml_tens, simid)
            write(*,'(a)') ' done'
            call chk_mmlcoup(Natom, nn_mml_tot, mml_listsize, mml_list, mml_tens, simid)
         end if

         call setup_lmm_arrays(Natom, max_no_mmlshells, nn_mml_tot, mml_listsize, mml_list, mml_tens, lmm_listsize, lmm_list, lmm_tens)
         !      lmm_tens=nint(1.0e3*lmm_tens)/1.0e3

         ! Print LMM interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print LMM interactions"
            call prn_lmmcoup(Natom, nn_mml_tot, lmm_listsize, lmm_list, lmm_tens, simid)
            write(*,'(a)') ' done'
         end if

         ! Setup reduced tensor for diagonal couplings
         if(mml_diag) then
            do i=1,Natom
               do j=1,mml_listsize(i)
                  lmm_tens_diag(:,j,i)=lmm_tens(:,1,1,j,i) 
               end do
            end do
         end if


      end if


      if(do_mmll==1) then

         ! Allocate and mount MMLL Hamiltonian
         ! Makes use of point symmetry operations for the coupling tensor
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for MMLL interaction'
         call setup_nm_nelem(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            nn_mmll_tot, max_no_mmllshells, max_no_equiv, sym, mmll_nn, mmll_redcoord, nme, nmdim, 3, &
            do_ralloy, Natom_full, acellnumb, atype_ch, &
            Nchmax, 81, .true., 2, 1, 1, mmll_inptens, mmll_symtens, mmll_symind)
         !Nchmax, hdim, do_tens_sym, couptensrank, invsym, timesym, couptens, fullcouptens)
         write(*,'(a)') ' done'

         write(*,*) 'nn_mmll_tot ', nn_mmll_tot
         call allocate_mmllhamiltoniandata(Natom, nn_mmll_tot, 1)

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount MMLL Hamiltonian'
         call setup_neighbour_latticehamiltonian(Natom, NT, NA, anumb, atype, &
            nn_mmll_tot, max_no_equiv, max_no_mmllshells, 3, &
            81, 1, mmll_listsize, mmll_nn, mmll_list, mmll_tens, nme, nmdim, mmll_symtens, mmll_symind, do_sortcoup, 'N', &
            do_ralloy, Natom_full, Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, &
            1, 1)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nme()

         ! Print MMLL interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print MMLL interactions"
            call prn_mmllcoup(Natom, nn_mmll_tot, mmll_listsize, mmll_list, mmll_tens, simid)
            write(*,'(a)') ' done'
         end if

         call setup_llmm_arrays(Natom, max_no_mmllshells, nn_mmll_tot, mmll_listsize, mmll_list, mmll_tens, llmm_listsize, llmm_list, llmm_tens)

         ! Print LLMM interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print LLMM interactions"
            call prn_llmmcoup(Natom, nn_mmll_tot, llmm_listsize, llmm_list, llmm_tens, simid)
            write(*,'(a)') ' done'
         end if

      end if

   end subroutine setup_latticehamiltonian


   !> Mounts terms of the lattice Hamiltonian
   subroutine setup_neighbour_latticehamiltonian(Natom, NT, NA, anumb, atype, &
         max_no_neigh, max_no_equiv, max_no_shells, nelem, &
         hdim, lexp, nlistsize, nn, nlist, ncoup, nm, nmdim, xc, nm_cell_symind, do_sortcoup, do_n3, &
         !hdim, lexp, nlistsize, nn, nlist, ncoup, nm, nmdim, xc, do_sortcoup, do_n3, &
      do_ralloy, Natom_full, Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, &
         lexpi, lexpj)
      !do_ralloy, Natom_full, Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp)
      !
      use Constants
      use Sorting, only : MergeSortIR
      !!! TMP
      use LatticeInputData, only : do_ll_phonopy
      !!! TMP

      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      integer, dimension(Natom), intent(in) :: atype !< Type of atom

      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for the lattice interaction
      integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for the lattice interaction
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for the lattice interaction
      integer, intent(in) :: nelem  !< Number of elements in each coupling (=1 two-site couplings, =2 for three-site couplings etc)

      integer, intent(in) :: hdim  !< Number of elements in Hamiltonian element (scalar or vector)
      integer, intent(in) :: lexpi  !< Power of the spin m_i. Needed for proper rescaling
      integer, intent(in) :: lexpj  !< Power of the spin m_j. Needed for proper rescaling
      integer, intent(in) :: lexp  !< Order of the interaction. Needed for proper rescaling (1=linear,2=quadratic)  <--- presently not used
      integer, dimension(Natom), intent(out):: nlistsize !< Size of neighbour list for the lattice interaction
      integer, dimension(NT), intent(in):: nn !< Number of neighbour shells
      integer, dimension(nelem, max_no_neigh, Natom), intent(out) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(hdim,max_no_neigh,Natom), intent(out) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(Natom, max_no_shells, max_no_equiv, nelem), intent(in) :: nm !< Neighbour map
      integer, dimension(max_no_shells, Natom), intent(in) :: nmdim !< Dimension of neighbour map
      real(dblprec), dimension(hdim, NT, max_no_shells, NT, NT, 48), intent(in) :: xc !< Coupling constants
      !real(dblprec), dimension(hdim, NT, max_no_shells, NT, NT, max_no_equiv), intent(in) :: xc !< Coupling constants
      integer, dimension(48,max_no_shells,na) :: nm_cell_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      !real(dblprec), dimension(hdim, NT, max_no_shells, Nchmax, NT), intent(in) :: xc !< Coupling constants
      character :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
      character(len=1) :: do_n3           !< Newton´s third law correction of force constant coefficient elements ('Y'/'N')

      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)
      real(dblprec), dimension(NA, Nchmax), intent(in) :: ammom_inp !< Magnetic moment directions from input (for alloys)
      !
      integer :: i, j, k, l
      integer :: iatom, jatom, ichtype, jchtype
      real(dblprec) :: fc_unitconv
      real(dblprec), dimension(hdim) :: tempcoup
      integer, dimension(nelem) :: tempn

      integer :: ncount, ncounter
      logical :: exis

      ncoup = 0.0_dblprec

      ! Factor for unit conversion from 
      ! to [ll_tens] = mRyd/angstrom**2
      if(do_ll_phonopy==1 .and. nelem==1 ) then 
         ! In phonopy the harmonic lattice potential is in ev / angstrom**2
         fc_unitconv = 1000.0_dblprec / Ry_ev
      else 
         !The standard harmonic lattice input is in mRyd / angstrom**2
         fc_unitconv = 1.0_dblprec
      end if

      write(*,*) 'fc_unitconv ', fc_unitconv

      iatom=0

      ! Loop over atoms
      do i=1, Natom
         ncount = 1
         ncounter = 1
         ! Shells
         do k=1, Nn(atype(i))
            ! Sites in shell
            if(do_hoc_debug==1) write(*,'(a,i4,a,i4,a,i4,a,i4)') 'i ', i, ' k ', k, ' nmdim ', nmdim(k,i)
            !do j=1, 1  !!! TMP !!!
            do j=1, nmdim(k,i)
               if(do_hoc_debug==1) write(*,'(a,i4,a,i4,a,i4,a,i4,a,i4)') 'i ', i, ' k ', k, ' j ', j, ' nm(i,k,j,1) ', nm(i,k,j,1), ' nm_cell_symind(j,k,atype(i)) ', nm_cell_symind(j,k,atype(i))
               ! Existing coupling
               if (nm(i,k,j,1)>0) then
                  exis=.false.
                  do l=1,ncount-1
                     if(nelem==1) then
                        if(nlist(1,l,i)==nm(i,k,j,1)) exis=.true.
                     elseif(nelem==2) then
                        if(nlist(1,l,i)==nm(i,k,j,1) .and. nlist(2,l,i)==nm(i,k,j,2)) exis=.true.
                     elseif(nelem==3) then
                        if(nlist(1,l,i)==nm(i,k,j,1) .and. nlist(2,l,i)==nm(i,k,j,2) &
                           .and. nlist(3,l,i)==nm(i,k,j,3)) exis=.true.
                     end if
                  end do
                  if(.not.exis) then
                     nlist(1:nelem,ncount,i) = nm(i,k,j,1:nelem)
                     if (do_ralloy==0) then
                        if(do_hoc_debug==1) then
                           !write(*,*) xc(:,atype(i),k,1,atype(nm(i,k,j,1)), j) * fc_unitconv
                           !write(*,'(a,81f10.6)') 'xc_inp ', xc(:,atype(i),k,1,1) * fc_unitconv
                           write(*,'(a,81f10.6)') 'xc_inp ', xc(:,atype(i),k,1,atype(nm(i,k,j,1)), nm_cell_symind(j,k,anumb(i))) * fc_unitconv
                           write(*,*) '  types: ', atype(i) , atype(nm(i,k,j,1)), nm(i,k,j,1), k
                           write(*,'(a,5i4, 9f8.4)') '  XC idx:', atype(i),k,1,atype(nm(i,k,j,1)), nm_cell_symind(j,k,anumb(i)), xc(:,atype(i),k,1,atype(nm(i,k,j,1)), nm_cell_symind(j,k,anumb(i))) 
                           !write(*,'(a,81f10.6)') 'xc_inp ', xc(:,atype(i),k,1,atype(nm(i,k,j,1)), nm_cell_symind(j,k,atype(i))) * fc_unitconv                           
                           !write(*,'(a,81f10.6)') 'xc_inp ', xc(:,atype(i),k,1,atype(nm(i,k,j,1)), j) * fc_unitconv
                           !write(*,'(a,81f10.6)') 'xc_inp ', xc(:,atype(i),k,1,atype(nm(i,k,j,1))) * fc_unitconv
                        end if
                        !ncoup(:,ncount,i) = xc(:,atype(i),k,1,atype(nm(i,k,j,1)), j) * fc_unitconv
                        !ncoup(:,ncount,i) = xc(:,atype(i),k,1,1) * fc_unitconv

                        !ncoup(:,ncount,i) = xc(:,atype(i),k,1,atype(nm(i,k,j,1))) * fc_unitconv

                        ncoup(:,ncount,i) = xc(:,atype(i),k,1,atype(nm(i,k,j,1)), nm_cell_symind(j,k,anumb(i))) * fc_unitconv &
                        !ncoup(:,ncount,i) = xc(:,atype(i),k,1,atype(nm(i,k,j,1)), nm_cell_symind(j,k,atype(i))) * fc_unitconv &                           
                           / ammom_inp(anumb(i),1)**lexpi / ammom_inp(anumb(nm(i,k,j,1)),1)**lexpj

                        !ncoup(:,ncount,i) = xc(:,atype(i),k,1,atype(nm(i,k,j,1)), j) * fc_unitconv &
                        !   / ammom_inp(anumb(i),1)**lexpi / ammom_inp(anumb(nm(i,k,j,1)),1)**lexpj

                        !ncoup(:,ncount,i) = xc(:,atype(i),k,1,atype(nm(i,k,j,1))) * fc_unitconv &
                        !   / ammom_inp(anumb(i),1)**lexpi / ammom_inp(anumb(nm(i,k,j,1)),1)**lexpj

                        ! *** From hamiltonianinit *** 
                        ! coup(:,ncount,i,iconf) = xc(:,atype(i),k,1,1,iconf) * fc2 &
                        !     / ammom_inp(anumb(i),1,iconf)**lexp / ammom_inp(anumb(nm(i,k,j)),1,iconf)**lexp
                        ! *** From hamiltonianinit *** 

                        if(do_hoc_debug==1) then
                           write(*,'(a,i4,a,3i4,a)') 'i ', i, '  nlist(1:nelem,ncount,i)', nlist(1:nelem,ncount,i)!, ' ncoup'
                           write(*,'(a,81f10.6)') 'xc_map ', ncoup(:,ncount,i)
                        end if
                     else
                        !
                        iatom = i
                        jatom = nm(i,k,j,1)
                        ichtype = atype_ch(iatom)
                        jchtype = atype_ch(jatom)
                        ncoup(:,ncount,i) = xc(:,atype_ch(i),k,achem_ch(i),achem_ch(jatom), j)  * fc_unitconv
                        !ncoup(:,ncount,i) = xc(:,atype_ch(i),k,achem_ch(i),achem_ch(jatom))  * fc_unitconv
                     end if
                     ncount = ncount + 1
                  end if
               end if
            end do
         end do
         nlistsize(i) = ncount-1
      end do

      if(do_sortcoup == 'Y') then
         do i=1,Natom
            do j=1,nlistsize(i)
               do k=1,nlistsize(i)-j
                  if ( nlist(1,k,i) .gt. nlist(1,k+1,i) ) then
                     tempn        = nlist(:,k,i)
                     nlist(:,k,i)   = nlist(:,k+1,i)
                     nlist(:,k+1,i) = tempn
                     tempcoup     = ncoup(:,k,i)
                     ncoup(:,k,i)   = ncoup(:,k+1,i)
                     ncoup(:,k+1,i) = tempcoup
                  end if
               end do
            end do
         end do
      end if

      ! Calculates the correct force constant coefficient element 
      ! \Psi_{ii} = -\sum_{j\neq i} Psi_{ij}
      ! as required by Newton´s third law.
      ! SLDTODO Applicable for higher order couplings?
      if(do_n3 == 'Y') then
         !write(*,*) 'Start do_N3'
         do i=1,Natom
            tempcoup = 0.0_dblprec
            !write(*,*) '-----------------'
            do j=1,nlistsize(i)
               if( nlist(1,j,i) == i ) then
                  k = j
                  !write(*,*) 'hit ', i, nlist(1,j,i), j
               else
                  tempcoup = tempcoup + ncoup(:,j,i)
                  !write(*,*) 'hop ', i, nlist(1,j,i), j
               end if
            end do
            ncoup(:,k,i) = -tempcoup
         end do
      end if

   end subroutine setup_neighbour_latticehamiltonian


   !> Deallocate neighbour map array
   subroutine deallocate_nm()

      implicit none

      integer :: i_stat, i_all

      i_all=-product(shape(nm))*kind(nm)
      deallocate(nm,stat=i_stat)
      call memocc(i_stat,i_all,'nm','deallocate_nm')
      i_all=-product(shape(nmdim))*kind(nmdim)
      deallocate(nmdim,stat=i_stat)
      call memocc(i_stat,i_all,'nmdim','deallocate_nm')

   end subroutine deallocate_nm


   !> Deallocate neighbour map array, element
   subroutine deallocate_nme()

      implicit none

      integer :: i_stat, i_all

      i_all=-product(shape(nme))*kind(nme)
      deallocate(nme,stat=i_stat)
      call memocc(i_stat,i_all,'nme','deallocate_nme')
      i_all=-product(shape(nmdim))*kind(nmdim)
      deallocate(nmdim,stat=i_stat)
      call memocc(i_stat,i_all,'nmdim','deallocate_nm')

   end subroutine deallocate_nme


   !> Setup the lm-tensor setup_lm_arrays
   subroutine setup_lm_arrays(Natom, max_no_mlshells, nn_ml_tot, ml_listsize, ml_list, ml_tens, lm_listsize, lm_list, lm_tens)
      !
      use Constants

      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_mlshells                            !< Limit of number of shells for ML interactions
      integer, intent(inout) :: nn_ml_tot                               !< Calculated number of neighbours with ML interactions

      integer, dimension(Natom), intent(in) :: ml_listsize          !< Size of neighbour list for ML
      integer, dimension(1,nn_ml_tot,Natom), intent(in) :: ml_list            !< List of neighbours for ML
      real(dblprec), dimension(3,3,nn_ml_tot,Natom), intent(in) :: ml_tens    !< ML coupling tensor

      integer, dimension(Natom), intent(out) :: lm_listsize          !< Size of neighbour list for LM
      integer, dimension(1,nn_ml_tot,Natom), intent(out) :: lm_list            !< List of neighbours for LM
      real(dblprec), dimension(3,3,nn_ml_tot,Natom), intent(out) :: lm_tens    !< LM coupling tensor

      ! Local scalars
      integer :: i, j, inew
      integer :: ia, ja
      integer, dimension(1) :: tmpneigh

      lm_listsize=0
      lm_list=0
      lm_tens=0_dblprec

      do i=1,Natom
         do j=1,ml_listsize(i)
            inew = ml_list(1,j,i)
            lm_listsize( inew ) = lm_listsize( inew ) + 1
            tmpneigh(1) = i
            !write(23,'(a,3i8,a,4i8)') 'j, i, ml_list(1,j,i) ', j, i, ml_list(1,j,i), '    tmpneigh ', tmpneigh(1)
            !write(23,'(a,4i8,a,4i8)') 'i, j, ml_list(1:2,j,i) ', i, j, ml_list(1:2,j,i), '    tmpneigh ', tmpneigh(1:2)
            !write(23,*) 'i, j, ml_list(1:2,j,i) ', i, j, ml_list(1:2,j,i)
            !write(23,*) 'tmpneigh ', tmpneigh(1:2)
            lm_list(1, lm_listsize(inew), inew) = tmpneigh(1)
            !lm_tens(1:3,1:3,1:3,lm_listsize( inew), inew) = ml_tens(1:3,1:3,1:3,j,i)
            do ja=1,3
               do ia=1,3
                  lm_tens(ja,ia,lm_listsize( inew), inew) = ml_tens(ia,ja,j,i)
                  write(26,*) ja, ia, lm_tens(ja,ia,lm_listsize( inew), inew)
               end do
            end do
            !write(24,'(a,8i8)') 'lm_listsize(inew), inew, lm_list(1:1, lm_listsize(inew), inew) ', &
            !     lm_listsize(inew), inew, lm_list(1:1, lm_listsize(inew), inew) 
            !write(24,*) 'i, j, lm_list(1:2,j,i) ', i, j, lm_list(1:2,j,i)
            !write(23,*) 'lm_list(1:2, lm_listsize(inew), inew) ', inew, lm_listsize(inew)
            !write(23,*) lm_list(1:2, lm_listsize(inew), inew)
         end do
      end do

   end subroutine setup_lm_arrays


   !> Setup the lmm-tensor setup_lmm_arrays
   subroutine setup_lmm_arrays(Natom, max_no_mmlshells, nn_mml_tot, mml_listsize, mml_list, mml_tens, lmm_listsize, lmm_list, lmm_tens)
      !
      use Constants

      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_mmlshells                            !< Limit of number of shells for MML interactions
      integer, intent(inout) :: nn_mml_tot                               !< Calculated number of neighbours with MML interactions

      integer, dimension(Natom), intent(in) :: mml_listsize          !< Size of neighbour list for MML
      integer, dimension(2,nn_mml_tot,Natom), intent(in) :: mml_list            !< List of neighbours for MML
      real(dblprec), dimension(3,3,3,nn_mml_tot,Natom), intent(in) :: mml_tens    !< MML coupling tensor

      integer, dimension(Natom), intent(out) :: lmm_listsize          !< Size of neighbour list for LMM
      integer, dimension(2,nn_mml_tot,Natom), intent(out) :: lmm_list            !< List of neighbours for LMM
      real(dblprec), dimension(3,3,3,nn_mml_tot,Natom), intent(out) :: lmm_tens    !< LMM coupling tensor

      ! Local scalars
      integer :: i, j, inew
      integer :: ia, ja, ka
      integer, dimension(2) :: tmpneigh

      lmm_listsize=0
      lmm_list=0
      lmm_tens=0_dblprec

      do i=1,Natom
         do j=1,mml_listsize(i)
            inew = mml_list(2,j,i)
            lmm_listsize( inew ) = lmm_listsize( inew ) + 1
            tmpneigh(1) = i
            tmpneigh(2) = mml_list(1,j,i)
            !write(23,'(a,4i8,a,4i8)') 'j, i, mml_list(1:2,j,i) ', j, i, mml_list(1:2,j,i), '    tmpneigh ', tmpneigh(1:2)
            !write(23,'(a,4i8,a,4i8)') 'i, j, mml_list(1:2,j,i) ', i, j, mml_list(1:2,j,i), '    tmpneigh ', tmpneigh(1:2)
            !write(23,*) 'i, j, mml_list(1:2,j,i) ', i, j, mml_list(1:2,j,i)
            !write(23,*) 'tmpneigh ', tmpneigh(1:2)
            lmm_list(1:2, lmm_listsize(inew), inew) = tmpneigh(1:2)
            !lmm_tens(1:3,1:3,1:3,lmm_listsize( inew), inew) = mml_tens(1:3,1:3,1:3,j,i)
            do ka=1,3
               do ja=1,3
                  do ia=1,3
                     lmm_tens(ka,ia,ja,lmm_listsize( inew), inew) = mml_tens(ia,ja,ka,j,i)
                     !write(26,*) ka, ia, ja, lmm_tens(ka,ia,ja,lmm_listsize( inew), inew)
                     !write(26,*) ka, ja, ia, lmm_tens(ka,ja,ia,lmm_listsize( inew), inew)
                  end do
               end do
            end do
            !write(24,'(a,8i8)') 'lmm_listsize(inew), inew, lmm_list(1:2, lmm_listsize(inew), inew) ', &
            !     lmm_listsize(inew), inew, lmm_list(1:2, lmm_listsize(inew), inew) 
            !write(24,*) 'i, j, lmm_list(1:2,j,i) ', i, j, lmm_list(1:2,j,i)
            !write(23,*) 'lmm_list(1:2, lmm_listsize(inew), inew) ', inew, lmm_listsize(inew)
            !write(23,*) lmm_list(1:2, lmm_listsize(inew), inew)
         end do
      end do

   end subroutine setup_lmm_arrays


   !> Setup the llmm-tensor setup_llmm_arrays
   subroutine setup_llmm_arrays(Natom, max_no_mmllshells, nn_mmll_tot, mmll_listsize, mmll_list, mmll_tens, llmm_listsize, llmm_list, llmm_tens)
      !
      use Constants

      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_mmllshells                            !< Limit of number of shells for MMLL interactions
      integer, intent(inout) :: nn_mmll_tot                               !< Calculated number of neighbours with MMLL interactions

      integer, dimension(Natom), intent(in) :: mmll_listsize          !< Size of neighbour list for MMLL
      integer, dimension(3,nn_mmll_tot,Natom), intent(in) :: mmll_list            !< List of neighbours for MMLL
      real(dblprec), dimension(3,3,3,3,nn_mmll_tot,Natom), intent(in) :: mmll_tens    !< MMLL coupling tensor

      integer, dimension(Natom), intent(out) :: llmm_listsize          !< Size of neighbour list for LLMM
      integer, dimension(3,nn_mmll_tot,Natom), intent(out) :: llmm_list            !< List of neighbours for LLMM
      real(dblprec), dimension(3,3,3,3,nn_mmll_tot,Natom), intent(out) :: llmm_tens    !< LLMM coupling tensor

      ! Local scalars
      integer :: i, j, inew
      integer :: ia, ja, ka, la
      integer, dimension(3) :: tmpneigh

      llmm_listsize=0
      llmm_list=0
      llmm_tens=0_dblprec

      do i=1,Natom
         do j=1,mmll_listsize(i)
            inew = mmll_list(2,j,i)
            tmpneigh(1) = mmll_list(3,j,i)
            llmm_listsize( inew ) = llmm_listsize( inew ) + 1
            tmpneigh(2) = i
            tmpneigh(3) = mmll_list(1,j,i)
            !write(23,'(a,5i8,a,5i8)') 'j, i, mmll_list(1:3,j,i) ', j, i, mmll_list(1:3,j,i), '    tmpneigh ', tmpneigh(1:3)
            !write(23,'(a,5i8,a,5i8)') 'i, j, mmll_list(1:3,j,i) ', i, j, mmll_list(1:3,j,i), '    tmpneigh ', tmpneigh(1:3)
            !write(23,'(a,8i8)') 'i, j, mmll_list(1:3,j,i) ', i, j, mmll_list(1:3,j,i)
            !write(23,'(a,8i8)') 'tmpneigh ', tmpneigh(1:3)
            llmm_list(1:3, llmm_listsize(inew), inew) = tmpneigh(1:3)
            do la=1,3
               do ka=1,3
                  do ja=1,3
                     do ia=1,3
                        llmm_tens(ka,la,ia,ja,llmm_listsize( inew), inew) = mmll_tens(ia,ja,ka,la,j,i)
                        write(26,*) ka, la, ia, ja, llmm_tens(ka,la,ia,ja,llmm_listsize( inew), inew)
                     end do
                  end do
               end do
            end do
            !write(24,'(a,8i8)') 'llmm_listsize(inew), inew, llmm_list(1:3, llmm_listsize(inew), inew) ', &
            !     llmm_listsize(inew), inew, llmm_list(1:3, llmm_listsize(inew), inew) 
            !write(24,'(a,8i8)') 'i, j, llmm_list(1:3,j,i) ', i, j, llmm_list(1:3,j,i)
            !write(23,*) 'llmm_list(1:2, llmm_listsize(inew), inew) ', inew, llmm_listsize(inew)
            !write(23,*) llmm_list(1:2, llmm_listsize(inew), inew)
         end do
      end do

   end subroutine setup_llmm_arrays


end module LatticeHamiltonianInit
