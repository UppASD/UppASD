!> Data and allocation routines for the lattice Hamiltonian
module LatticeHamiltonianData
   use Parameters
   use Profiling
   !
   implicit none
   !
   !From setup
   integer ::  max_no_llneigh    !< Calculated maximum of neighbours for LL interaction
   integer, dimension(:), allocatable :: ll_listsize      !< Size of neighbour list for LL
   integer, dimension(:,:,:), allocatable :: ll_list        !< List of neighbours for LL 
   real(dblprec), dimension(:,:,:,:), allocatable :: ll_tens    !< Harmonic lattice forces (LL) 

   integer ::  max_no_lllneigh    !< Calculated maximum of neighbours for LLL interaction
   integer, dimension(:), allocatable :: lll_listsize      !< Size of neighbour list for LLL
   integer, dimension(:,:,:), allocatable :: lll_list        !< List of neighbours for LLL 
   real(dblprec), dimension(:,:,:,:,:), allocatable :: lll_tens    !< Anharmonic lattice forces (LLL) 

   integer ::  max_no_llllneigh    !< Calculated maximum of neighbours for LLLL interaction
   integer, dimension(:), allocatable :: llll_listsize      !< Size of neighbour list for LLLL
   integer, dimension(:,:,:), allocatable :: llll_list        !< List of neighbours for LLLL 
   real(dblprec), dimension(:,:,:,:,:,:), allocatable :: llll_tens    !< Quartic lattice forces (LLLL) 

   integer ::  max_no_mlneigh    !< Calculated maximum of neighbours for ML interaction
   integer, dimension(:), allocatable :: ml_listsize      !< Size of neighbour list for ML
   integer, dimension(:,:,:), allocatable :: ml_list        !< List of neighbours for ML 
   real(dblprec), dimension(:,:,:,:), allocatable :: ml_tens    !< spin-lattice coupling (ML)
   integer, dimension(:), allocatable :: lm_listsize      !< Size of neighbour list for LM
   integer, dimension(:,:,:), allocatable :: lm_list        !< List of neighbours for LM 
   real(dblprec), dimension(:,:,:,:), allocatable :: lm_tens    !< spin-lattice coupling (LM)

   integer ::  max_no_mmlneigh    !< Calculated maximum of neighbours for MML interaction
   integer, dimension(:), allocatable :: mml_listsize      !< Size of neighbour list for MML
   integer, dimension(:,:,:), allocatable :: mml_list        !< List of neighbours for MML 
   real(dblprec), dimension(:,:,:,:,:), allocatable :: mml_tens    !< spin-lattice coupling (MML)
   real(dblprec), dimension(:,:,:), allocatable :: mml_tens_diag    !< Diagonal (xc-striction) spin-lattice coupling (MML)
   integer, dimension(:), allocatable :: lmm_listsize      !< Size of neighbour list for LMM
   integer, dimension(:,:,:), allocatable :: lmm_list        !< List of neighbours for LMM 
   real(dblprec), dimension(:,:,:,:,:), allocatable :: lmm_tens    !< spin-lattice coupling (LMM)
   real(dblprec), dimension(:,:,:), allocatable :: lmm_tens_diag    !< spin-lattice coupling (LMM)

   integer ::  max_no_mmllneigh    !< Calculated maximum of neighbours for MMLL interaction
   integer, dimension(:), allocatable :: mmll_listsize      !< Size of neighbour list for MMLL
   integer, dimension(:,:,:), allocatable :: mmll_list        !< List of neighbours for MMLL 
   real(dblprec), dimension(:,:,:,:,:,:), allocatable :: mmll_tens   !< spin-lattice coupling (MMLL)   
   integer, dimension(:), allocatable :: llmm_listsize      !< Size of neighbour list for LLMM
   integer, dimension(:,:,:), allocatable :: llmm_list        !< List of neighbours for LLMM 
   real(dblprec), dimension(:,:,:,:,:,:), allocatable :: llmm_tens    !< spin-lattice coupling (LLMM)

   real(dblprec) :: mm_energy0_const  !< Ground-state energy
   logical       :: mm_energy0_calc = .false. !< Is GS energy calculated?

   public 


contains 


   !> Allocate arrays for LL Hamiltonian
   subroutine allocate_llhamiltoniandata(Natom,nn_ll_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nn_ll_tot !< Calculated number of neighbours with LL interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ll_listsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ll_listsize))*kind(ll_listsize),'ll_listsize','allocate_llhamiltoniandata')
         allocate(ll_list(1,nn_ll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ll_list))*kind(ll_list),'ll_list','allocate_llhamiltoniandata')
         allocate(ll_tens(3,3,nn_ll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ll_tens))*kind(ll_tens),'ll_tens','allocate_llhamiltoniandata')
      else
         i_all=-product(shape(ll_listsize))*kind(ll_listsize)
         deallocate(ll_listsize,stat=i_stat)
         call memocc(i_stat,i_all,'ll_listsize','allocate_llhamiltoniandata')
         i_all=-product(shape(ll_list))*kind(ll_list)
         deallocate(ll_list,stat=i_stat)
         call memocc(i_stat,i_all,'ll_list','allocate_llhamiltoniandata')
         i_all=-product(shape(ll_tens))*kind(ll_tens)
         deallocate(ll_tens,stat=i_stat)
         call memocc(i_stat,i_all,'ll_tens','allocate_llhamiltoniandata')
      end if

   end subroutine allocate_llhamiltoniandata


   !> Allocate arrays for LLL Hamiltonian
   subroutine allocate_lllhamiltoniandata(Natom,nn_lll_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nn_lll_tot !< Calculated number of neighbours with LLL interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(lll_listsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(lll_listsize))*kind(lll_listsize),'lll_listsize','allocate_lllhamiltoniandata')
         allocate(lll_list(2,nn_lll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(lll_list))*kind(lll_list),'lll_list','allocate_lllhamiltoniandata')
         allocate(lll_tens(3,3,3,nn_lll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(lll_tens))*kind(lll_tens),'lll_tens','allocate_lllhamiltoniandata')
      else
         i_all=-product(shape(lll_listsize))*kind(lll_listsize)
         deallocate(lll_listsize,stat=i_stat)
         call memocc(i_stat,i_all,'lll_listsize','allocate_lllhamiltoniandata')
         i_all=-product(shape(lll_list))*kind(lll_list)
         deallocate(lll_list,stat=i_stat)
         call memocc(i_stat,i_all,'lll_list','allocate_lllhamiltoniandata')
         i_all=-product(shape(lll_tens))*kind(lll_tens)
         deallocate(lll_tens,stat=i_stat)
         call memocc(i_stat,i_all,'lll_tens','allocate_lllhamiltoniandata')
      end if

   end subroutine allocate_lllhamiltoniandata


   !> Allocate arrays for LLLL Hamiltonian
   subroutine allocate_llllhamiltoniandata(Natom,nn_llll_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nn_llll_tot !< Calculated number of neighbours with LLLL interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(llll_listsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(llll_listsize))*kind(llll_listsize),'llll_listsize','allocate_llllhamiltoniandata')
         allocate(llll_list(3,nn_llll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(llll_list))*kind(llll_list),'llll_list','allocate_llllhamiltoniandata')
         allocate(llll_tens(3,3,3,3,nn_llll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(llll_tens))*kind(llll_tens),'llll_tens','allocate_llllhamiltoniandata')
      else
         i_all=-product(shape(llll_listsize))*kind(llll_listsize)
         deallocate(llll_listsize,stat=i_stat)
         call memocc(i_stat,i_all,'llll_listsize','allocate_llllhamiltoniandata')
         i_all=-product(shape(llll_list))*kind(llll_list)
         deallocate(llll_list,stat=i_stat)
         call memocc(i_stat,i_all,'llll_list','allocate_llllhamiltoniandata')
         i_all=-product(shape(llll_tens))*kind(llll_tens)
         deallocate(llll_tens,stat=i_stat)
         call memocc(i_stat,i_all,'llll_tens','allocate_llllhamiltoniandata')
      end if

   end subroutine allocate_llllhamiltoniandata


   !> Allocate arrays for ML Hamiltonian
   subroutine allocate_mlhamiltoniandata(Natom,nn_ml_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nn_ml_tot !< Calculated number of neighbours with ML interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ml_listsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ml_listsize))*kind(ml_listsize),'ml_listsize','allocate_mlhamiltoniandata')
         allocate(ml_list(1,nn_ml_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ml_list))*kind(ml_list),'ml_list','allocate_mlhamiltoniandata')
         allocate(ml_tens(3,3,nn_ml_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ml_tens))*kind(ml_tens),'ml_tens','allocate_mlhamiltoniandata')
         allocate(lm_listsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(lm_listsize))*kind(lm_listsize),'lm_listsize','allocate_mlhamiltoniandata')
         allocate(lm_list(1,nn_ml_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(lm_list))*kind(lm_list),'lm_list','allocate_mlhamiltoniandata')
         allocate(lm_tens(3,3,nn_ml_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(lm_tens))*kind(lm_tens),'lm_tens','allocate_mlhamiltoniandata')
      else
         i_all=-product(shape(ml_listsize))*kind(ml_listsize)
         deallocate(ml_listsize,stat=i_stat)
         call memocc(i_stat,i_all,'ml_listsize','allocate_mlhamiltoniandata')
         i_all=-product(shape(ml_list))*kind(ml_list)
         deallocate(ml_list,stat=i_stat)
         call memocc(i_stat,i_all,'ml_list','allocate_mlhamiltoniandata')
         i_all=-product(shape(ml_tens))*kind(ml_tens)
         deallocate(ml_tens,stat=i_stat)
         call memocc(i_stat,i_all,'ml_tens','allocate_mlhamiltoniandata')
         i_all=-product(shape(lm_listsize))*kind(lm_listsize)
         deallocate(lm_listsize,stat=i_stat)
         call memocc(i_stat,i_all,'lm_listsize','allocate_mlhamiltoniandata')
         i_all=-product(shape(lm_list))*kind(lm_list)
         deallocate(lm_list,stat=i_stat)
         call memocc(i_stat,i_all,'lm_list','allocate_mlhamiltoniandata')
         i_all=-product(shape(lm_tens))*kind(lm_tens)
         deallocate(lm_tens,stat=i_stat)
         call memocc(i_stat,i_all,'lm_tens','allocate_mlhamiltoniandata')
      end if

   end subroutine allocate_mlhamiltoniandata


   !> Allocate arrays for MML Hamiltonian
   subroutine allocate_mmlhamiltoniandata(Natom,nn_mml_tot,flag)
      use LatticeInputData, only : mml_diag
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nn_mml_tot !< Calculated number of neighbours with MML interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(mml_listsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(mml_listsize))*kind(mml_listsize),'mml_listsize','allocate_mmlhamiltoniandata')
         allocate(mml_list(2,nn_mml_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(mml_list))*kind(mml_list),'mml_list','allocate_mmlhamiltoniandata')
         allocate(mml_tens(3,3,3,nn_mml_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(mml_tens))*kind(mml_tens),'mml_tens','allocate_mmlhamiltoniandata')
         allocate(lmm_listsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(lmm_listsize))*kind(lmm_listsize),'lmm_listsize','allocate_mmlhamiltoniandata')
         allocate(lmm_list(2,nn_mml_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(lmm_list))*kind(lmm_list),'lmm_list','allocate_mmlhamiltoniandata')
         allocate(lmm_tens(3,3,3,nn_mml_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(lmm_tens))*kind(lmm_tens),'lmm_tens','allocate_mmlhamiltoniandata')
         if(mml_diag) then
            allocate(mml_tens_diag(3,nn_mml_tot,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(mml_tens_diag))*kind(mml_tens_diag),'mml_tens_diag','allocate_mmlhamiltoniandata')
            allocate(lmm_tens_diag(3,nn_mml_tot,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(lmm_tens_diag))*kind(lmm_tens_diag),'lmm_tens_diag','allocate_mmlhamiltoniandata')
         end if

      else
         i_all=-product(shape(mml_listsize))*kind(mml_listsize)
         deallocate(mml_listsize,stat=i_stat)
         call memocc(i_stat,i_all,'mml_listsize','allocate_mmlhamiltoniandata')
         i_all=-product(shape(mml_list))*kind(mml_list)
         deallocate(mml_list,stat=i_stat)
         call memocc(i_stat,i_all,'mml_list','allocate_mmlhamiltoniandata')
         i_all=-product(shape(mml_tens))*kind(mml_tens)
         deallocate(mml_tens,stat=i_stat)
         call memocc(i_stat,i_all,'mml_tens','allocate_mmlhamiltoniandata')
         i_all=-product(shape(lmm_listsize))*kind(lmm_listsize)
         deallocate(lmm_listsize,stat=i_stat)
         call memocc(i_stat,i_all,'lmm_listsize','allocate_mmlhamiltoniandata')
         i_all=-product(shape(lmm_list))*kind(lmm_list)
         deallocate(lmm_list,stat=i_stat)
         call memocc(i_stat,i_all,'lmm_list','allocate_mmlhamiltoniandata')
         i_all=-product(shape(lmm_tens))*kind(lmm_tens)
         deallocate(lmm_tens,stat=i_stat)
         call memocc(i_stat,i_all,'lmm_tens','allocate_mmlhamiltoniandata')
         if(mml_diag) then
            i_all=-product(shape(mml_tens_diag))*kind(mml_tens_diag)
            deallocate(mml_tens_diag,stat=i_stat)
            call memocc(i_stat,i_all,'mml_tens_diag','allocate_mmlhamiltoniandata')
            i_all=-product(shape(lmm_tens_diag))*kind(lmm_tens_diag)
            deallocate(lmm_tens_diag,stat=i_stat)
            call memocc(i_stat,i_all,'lmm_tens_diag','allocate_mmlhamiltoniandata')
         end if
      end if

   end subroutine allocate_mmlhamiltoniandata


   !> Allocate arrays for MMLL Hamiltonian
   subroutine allocate_mmllhamiltoniandata(Natom,nn_mmll_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nn_mmll_tot !< Calculated number of neighbours with MMLL interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(mmll_listsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(mmll_listsize))*kind(mmll_listsize),'mmll_listsize','allocate_mmllhamiltoniandata')
         allocate(mmll_list(3,nn_mmll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(mmll_list))*kind(mmll_list),'mmll_list','allocate_mmllhamiltoniandata')
         allocate(mmll_tens(3,3,3,3,nn_mmll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(mmll_tens))*kind(mmll_tens),'mmll_tens','allocate_mmllhamiltoniandata')
         allocate(llmm_listsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(llmm_listsize))*kind(llmm_listsize),'llmm_listsize','allocate_mmllhamiltoniandata')
         allocate(llmm_list(3,nn_mmll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(llmm_list))*kind(llmm_list),'llmm_list','allocate_mmllhamiltoniandata')
         allocate(llmm_tens(3,3,3,3,nn_mmll_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(llmm_tens))*kind(llmm_tens),'llmm_tens','allocate_mmllhamiltoniandata')
      else
         i_all=-product(shape(mmll_listsize))*kind(mmll_listsize)
         deallocate(mmll_listsize,stat=i_stat)
         call memocc(i_stat,i_all,'mmll_listsize','allocate_mmllhamiltoniandata')
         i_all=-product(shape(mmll_list))*kind(mmll_list)
         deallocate(mmll_list,stat=i_stat)
         call memocc(i_stat,i_all,'mmll_list','allocate_mmllhamiltoniandata')
         i_all=-product(shape(mmll_tens))*kind(mmll_tens)
         deallocate(mmll_tens,stat=i_stat)
         call memocc(i_stat,i_all,'mmll_tens','allocate_mmllhamiltoniandata')
         i_all=-product(shape(llmm_listsize))*kind(llmm_listsize)
         deallocate(llmm_listsize,stat=i_stat)
         call memocc(i_stat,i_all,'llmm_listsize','allocate_mmllhamiltoniandata')
         i_all=-product(shape(llmm_list))*kind(llmm_list)
         deallocate(llmm_list,stat=i_stat)
         call memocc(i_stat,i_all,'llmm_list','allocate_mmllhamiltoniandata')
         i_all=-product(shape(llmm_tens))*kind(llmm_tens)
         deallocate(llmm_tens,stat=i_stat)
         call memocc(i_stat,i_all,'llmm_tens','allocate_mmllhamiltoniandata')
      end if

   end subroutine allocate_mmllhamiltoniandata


end module LatticeHamiltonianData

