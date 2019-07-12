!> Wolff cluster Monte Carlo
!> @author
!> Jonathan Chico
!> @copyright
!> GNU Public License
module MC_WOLFF

   use Parameters
   use Constants
   use Profiling

   implicit none

   integer, dimension(:), allocatable :: wolff_list_neigh_size
   integer, dimension(:,:), allocatable :: wolff_list_neigh

   public :: wolff_list_neigh_size, wolff_list_neigh

contains

   ! This subroutine creates the neighbour list for the wolff algorithm in which the nearest neighbours are needed
   subroutine create_wolff_neigh_list(Natom, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, &
         max_no_neigh, nlistsize, nlist, coord,wolff_list_neigh_size,wolff_list_neigh)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: N1    !< Number of cell repetitions in x direction
      integer, intent(in) :: N2    !< Number of cell repetitions in y direction
      integer, intent(in) :: N3    !< Number of cell repetitions in z direction
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh, Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms

      integer, dimension(:), allocatable :: wolff_list_neigh_size
      integer, dimension(:,:), allocatable :: wolff_list_neigh

      real(dblprec), parameter :: rcutoff=1.05_dblprec

      integer :: iatom, jneigh, jatom, ix, iy, iz, icount
      integer :: i_stat
      integer :: signx, signy, signz
      integer :: perx, pery, perz
      real(dblprec) :: jnorm, jpnorm, jtnorm
      real(dblprec), dimension(3) :: icoord, jcoord, jpcoord, jtcoord
      real(dblprec), dimension(3,-1:1,-1:1,-1:1) :: minmat_coord
      real(dblprec), dimension(-1:1,-1:1,-1:1) :: minmat_norm
      integer, dimension(-1:1,-1:1,-1:1) :: minmat_atom

      ! Allocate storage arrays
      allocate(wolff_list_neigh_size(Natom),stat=i_stat)
      call memocc(i_stat,product(shape(wolff_list_neigh_size))*kind(wolff_list_neigh_size),'wolff_neigh_list_size','create_wolff_neigh_list')
      allocate(wolff_list_neigh(26,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(wolff_list_neigh))*kind(wolff_list_neigh),'wolff_list_neigh','create_wolff_neigh_list')

      ! Prepare for periodic boundaries
      if(BC1=='P') then
         perx=1
      else
         perx=0
      end if
      if(BC2=='P') then
         pery=1
      else
         pery=0
      end if
      if(BC3=='P') then
         perz=1
      else
         perz=0
      end if

      ! Loop over atoms
      do iatom=1, Natom
         icoord=coord(1:3,iatom)

         ! Initialize minima matrix with large entries
         minmat_norm=rcutoff
         minmat_coord=0.0_dblprec
         minmat_atom=0

         ! Loop over neighbours
         do jneigh=1, nlistsize(iatom)
            jatom=nlist(jneigh,iatom)
            jcoord=coord(1:3,jatom)-icoord
            jnorm=sqrt(jcoord(1)**2+jcoord(2)**2+jcoord(3)**2)

            ! Fix for periodicity if present
            jpcoord=jcoord
            jpnorm=jnorm
            do ix=-perx,perx
               do iy=-pery,pery
                  do iz=-perz,perz
                     jtcoord(1)=jcoord(1)+ix*N1*C1(1)+iy*N2*C2(1)+iz*N3*C3(1)
                     jtcoord(2)=jcoord(2)+ix*N1*C1(2)+iy*N2*C2(2)+iz*N3*C3(2)
                     jtcoord(3)=jcoord(3)+ix*N1*C1(3)+iy*N2*C2(3)+iz*N3*C3(3)
                     jtnorm=sqrt(jtcoord(1)**2+jtcoord(2)**2+jtcoord(3)**2)
                     if(jtnorm<jpnorm) then
                        jpcoord=jtcoord
                        jpnorm=jtnorm
                     end if
                  end do
               end do
            end do
            jcoord=jpcoord
            jnorm=jpnorm

            ! Loop over all quadrants and directions to find closest neighbours
            do ix=-1,1
               signx=ix !(-1)**ix
               do iy=-1,1
                  signy=iy !(-1)**iy
                  do iz=-1,1
                     signz=iz !(-1)**iz
                     if(jnorm<minmat_norm(ix,iy,iz)) then
                        if((signx>0.and.jcoord(1)>0.0_dblprec.or.signx==0.and.abs(jcoord(1))<dbl_tolerance &
                           .or.signx<0.and.jcoord(1)<0.0_dblprec) &
                           .and.(signy>0.and.jcoord(2)>0.0_dblprec.or.signy==0 &
                           .and.abs(jcoord(2))<dbl_tolerance.or.signy<0.and.jcoord(2)<0.0_dblprec) &
                           .and.(signz>0.and.jcoord(3)>0.0_dblprec.or.signz==0 &
                           .and.abs(jcoord(3))<dbl_tolerance.or.signz<0.and.jcoord(3)<0.0_dblprec) &
                           ) then
                           minmat_coord(1:3,ix,iy,iz)=jcoord
                           minmat_norm(ix,iy,iz)=jnorm
                           minmat_atom(ix,iy,iz)=jatom
                        end if
                     end if
                  end do
               end do
            end do
         end do

         ! Make neighbour list
         icount=0
         do ix=-1,1
            do iy=-1,1
               do iz=-1,1
                  if(minmat_atom(ix,iy,iz)>0) then
                     icount=icount+1
                     wolff_list_neigh(icount,iatom)=minmat_atom(ix,iy,iz)
                  end if
               end do
            end do
         end do
         wolff_list_neigh_size(iatom)=icount
      end do

   end subroutine create_wolff_neigh_list

   ! This createds the cluster for the Wolff algorithm
   subroutine create_wolff_cluster(k_ensem,Natom,nHam,Mensemble,emom,nlistsize,nlist,&
         Temperature,conf_num,ncoup,max_no_neigh,visited_atoms)

      use RandomNumbers
      implicit none

      integer, intent(in) :: k_ensem
      integer, intent(in) :: Natom
      integer, intent(in) :: nHam
      integer, intent(in) :: Mensemble
      integer, intent(in) :: max_no_neigh

      integer, dimension(Natom), intent(inout) :: visited_atoms
      integer, dimension(:), allocatable :: wolffsnake
      integer, dimension(:), allocatable :: border_atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom

      integer, dimension(nHam), intent(in) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist
      integer, intent(in) :: conf_num !> number of configurations for LSF
      real(dblprec), intent(in) :: Temperature

      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup

      real(dblprec) :: dflip(1)
      integer :: border_index, i_border, snake_head, snake_tail
      integer :: atom_wolff_len, selected_atom,curr_len
      logical :: loop_closed, failure, total_failure,new_border

      total_failure=.false.
      loop_closed=.false.
      failure=.false.

      allocate(wolffsnake(Natom+1))
      allocate(border_atoms(Natom))
      visited_atoms=0

      ! Select the seed atom for the algorithm
      call rng_uniform(dflip,1)
      selected_atom=int(dflip(1)*Natom)+1

      ! Mark the seed atom as visited, later on this array will be used to know which spins will be flipped
      visited_atoms(selected_atom)=1
      ! Increase the size of the number of accepted atoms, as the seed has been accepted
      atom_wolff_len=atom_wolff_len+1
      ! Store the atom in the wolff snake to be able to afterward chek it's neighbours
      wolffsnake(atom_wolff_len)=selected_atom
      ! Initiliaze the border atom, setting up the border index to one and the seed atom as the first border atom
      border_index=1
      border_atoms(1)=selected_atom

      ! Check the neighbours of the seed and see how it goes
      call neighbour_check(k_ensem,Natom,nHam,Mensemble,selected_atom,emom,nlist,nlistsize,max_no_neigh,&
         border_atoms,atom_wolff_len,wolffsnake,border_index,curr_len,new_border,conf_num,ncoup,&
         Temperature,visited_atoms)

      ! Start the snake head at 1
      snake_head=1
      ! start the tail at the number of border atoms accepted from the first atom
      snake_tail=border_index

      ! If no atoms are found acceptable the loop just wont start
      if (.not.new_border) total_failure=.true.

      ! Do the loop intil the cluster is created
      do while (.not. loop_closed .and. .not. total_failure)

         curr_len=0
         ! Do a loop over the border atoms of the previous step
         do i_border=snake_head,snake_tail

            ! Temporarly store the index of the selected atom
            selected_atom=border_atoms(i_border)

            ! Check the neighbours of the recently selected atom and see if one needs to add more
            call neighbour_check(k_ensem,Natom,nHam,Mensemble,selected_atom,emom,nlist,nlistsize,max_no_neigh,&
               border_atoms,atom_wolff_len,wolffsnake,border_index,curr_len,new_border,conf_num,ncoup,&
               Temperature,visited_atoms)

         enddo

         ! If there is no new border atoms added the loop stops
         if (.not.new_border.and.curr_len.eq.0) then
            loop_closed=.true.
            !       write(*,*) 'sum',sum(visited_atoms)
         endif

         ! The new start of the cluster check index
         snake_head=snake_tail
         ! The new end of the cluster check index
         snake_tail=snake_head+curr_len

      enddo

   end subroutine create_wolff_cluster

   ! Check the neighbours to see if they can be added into the cluster
   subroutine neighbour_check(k_ensem,Natom,nHam,Mensemble,selected_atom,emom,nlist,nlistsize,max_no_neigh,&
         border_atoms,atom_wolff_len,wolffsnake,border_index,curr_len,new_border,conf_num,ncoup,&
         Temperature,visited_atoms)

      use RandomNumbers

      implicit none

      integer, intent(in) :: k_ensem
      integer, intent(in) :: Natom
      integer, intent(in) :: nHam
      integer, intent(in) :: Mensemble
      integer, intent(in) :: selected_atom
      integer, intent(in) :: max_no_neigh
      integer, intent(inout) :: atom_wolff_len
      integer, intent(inout) :: border_index
      integer, dimension(Natom), intent(inout) :: visited_atoms
      integer, dimension(Natom), intent(inout) :: border_atoms
      integer, dimension(Natom), intent(inout) :: wolffsnake
      real(dblprec), intent(in) :: Temperature
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
      integer, dimension(nHam), intent(in) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist

      integer, intent(inout) :: curr_len
      logical, intent(out) :: new_border

      integer :: neigh_loop
      real(dblprec) :: flipprob(1)

      new_border=.false.

      ! Do a loop over the neighbours of the selected atom
      do neigh_loop=1,nlistsize(selected_atom)
         ! Check only if it does not belong to the freaking cluster already
         if ( visited_atoms(nlist(neigh_loop,selected_atom)).ne.1) then
            ! Do a check for ISING MODEL first if the spins are parallel
            if (abs((emom(1,selected_atom,k_ensem)*emom(1,nlist(neigh_loop,selected_atom),k_ensem)+&
               emom(2,selected_atom,k_ensem)*emom(2,nlist(neigh_loop,selected_atom),k_ensem)+&
               emom(3,selected_atom,k_ensem)*emom(3,nlist(neigh_loop,selected_atom),k_ensem) )- 1.0_dblprec)<dbl_tolerance) then

               ! Random number to check the probability
               call rng_uniform(flipprob,1)
               ! See if the link is established
               if ( prob_add(Natom,nHam,max_no_neigh,selected_atom,neigh_loop,Temperature,conf_num,ncoup) ) then
                  ! Adding the neighbour to the list of visited atoms
                  visited_atoms(nlist(neigh_loop,selected_atom))=1
                  ! Increasing the size of the list of atoms
                  atom_wolff_len=atom_wolff_len+1
                  ! Storing the atom in the snake of this way one can see the neighbours
                  wolffsnake(atom_wolff_len)=nlist(neigh_loop,selected_atom)
                  ! increase the border index by one
                  border_index=border_index+1
                  ! add the currently selected atom to the border
                  border_atoms(border_index)=nlist(neigh_loop,selected_atom)
                  ! Increase the current length of the added border
                  curr_len=curr_len+1
                  ! As a new border atom that has been accepted has been found the new border flag is true
                  new_border=.true.
               endif

            endif

         endif

      enddo


   end subroutine neighbour_check

   ! This function decides if a link between two atoms is added with a certain probability
   logical function prob_add(Natom,nHam,max_no_neigh,i_atom,j_atom,Temperature,conf_num,ncoup)

      use RandomNumbers
      use Constants
      use HamiltonianData, only : ham

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: nHam
      integer, intent(in) :: max_no_neigh
      integer, intent(in) :: i_atom
      integer, intent(in) :: j_atom
      real(dblprec), intent(in) :: Temperature
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup
      real(dblprec) :: beta
      real(dblprec) :: flipprob(1), prob

      call rng_uniform(flipprob,1)

      ! Right now only nearest neighbour exchange is included
      beta=1.0_dblprec/k_bolt/Temperature
      prob=1.0_dblprec-exp(-2.0_dblprec*ncoup(j_atom,ham%aham(i_atom),1)*beta*mub)

      if (flipprob(1) < prob ) then
         prob_add=.true.
      else
         prob_add=.false.
      endif

   end function prob_add

end module MC_WOLFF
