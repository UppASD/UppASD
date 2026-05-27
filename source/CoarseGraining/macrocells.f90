!-------------------------------------------------------------------------------
! MODULE: MACROCELLS
!
! DESCRIPTION:
!> @brief Used for the subdivision of the lattice into cubic macrocells to simplify the
!> calculation of the dipolar interaction.
!
!> @author
!> Jonathan Chico
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module macrocells

   use Parameters
   use Profiling

   implicit none

   integer :: Num_dip   !< Number of dipolar cells where the macrocell dipole is written
   integer :: Num_macro !< Number of macrocells in the system
   integer :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
   character(len=1) :: do_macro_cells !< Flag to whether perform macrocell decomposition
   character(len=1) :: prn_dip_subset !< Flag to print the macro dipole-dipole field over a subset of cells
   character(len=20) :: dip_file      !< File containing the indexes of the cells where the macro-dipole field will be written
   integer, dimension(:), allocatable :: cell_index      !< Macrocell index for each atom
   integer, dimension(:), allocatable :: dipole_subset   !< List of the cells where the macrocell dipolar interaction is printed
   integer, dimension(:), allocatable :: macro_alistsize !< Number of atoms per macrocell
   integer, dimension(:,:), allocatable :: macro_atom_alist !< List containing the information of which atoms are in a given macrocell
   integer, dimension(:), allocatable :: atom_in_macro_pos !< Position of an atom within its macrocell
   integer :: max_macro_halo_size !< Maximum number of unique halo atoms for any macrocell
   integer :: max_macro_cell_neigh !< Maximum number of neighbouring macrocells for any macrocell
   integer :: max_macro_atom_local_neigh !< Maximum number of valid local neighbours for any atom in any macrocell
   integer, dimension(:), allocatable :: macro_halo_nlistsize !< Number of unique halo atoms per macrocell
   integer, dimension(:,:), allocatable :: macro_halo_to_global !< Local-halo to global atom map (local_halo, macrocell)
   integer, dimension(:,:), allocatable :: macro_atom_to_global !< Local-atom to global atom map (local_atom, macrocell)
   integer, dimension(:,:), allocatable :: macro_atom_local_nlistsize !< Number of neighbours for local atoms (local_atom, macrocell)
   integer, dimension(:,:,:), allocatable :: macro_atom_local_nlist !< Atom neighbour list in macrocell-local halo indices (max_macro_atom_local_neigh, local_atom, macrocell)
   integer, dimension(:), allocatable :: macro_cell_nlistsize !< Number of neighbouring macrocells per macrocell
   integer, dimension(:,:), allocatable :: macro_cell_nlist !< Macrocell neighbour map (local_cell_neigh, macrocell)

   real(dblprec), dimension(:,:), allocatable :: mmom_macro !< Magnitude of the macrocell magnetic moments
   real(dblprec), dimension(:,:), allocatable :: max_coord_macro !< Maximum value of the coordinates per cell
   real(dblprec), dimension(:,:), allocatable :: min_coord_macro !< Minimum value of the coordinates per cell
   real(dblprec), dimension(:,:), allocatable :: mid_coord_macro !< Midpoint of each cell
   real(dblprec), dimension(:,:,:), allocatable :: emom_macro  !< Unit vector of the macrocell magnetic moment
   real(dblprec), dimension(:,:,:), allocatable :: emomM_macro !< The full vector of the macrocell magnetic moment
   
   integer :: block_size_x = 1
   integer :: block_size_y = 1
   integer :: block_size_z = 1

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: print_macrocell_neighbour_map
   !
   !> @brief Print macrocell membership map in row format:
   !> i_cell  n_atoms  atom1 atom2 ...
   !> Here atom* corresponds to atom indices belonging to i_cell.
   !----------------------------------------------------------------------------
   subroutine print_macrocell_neighbour_map(Num_macro, max_num_atom_macro_cell, &
         macro_alistsize, macro_atom_alist, simid)

      implicit none

      integer, intent(in) :: Num_macro
      integer, intent(in) :: max_num_atom_macro_cell
      integer, dimension(Num_macro), intent(in) :: macro_alistsize
      integer, dimension(max_num_atom_macro_cell,Num_macro), intent(in) :: macro_atom_alist
      character(len=8), intent(in) :: simid

      integer :: ofileno, i_cell, j, n_neigh, i_err
      character(len=64) :: out_file

      out_file = 'macro_cell_alist_mapping.' // trim(simid) // '.out'
      open(newunit=ofileno, file=out_file, status='replace', action='write', iostat=i_err)
      if (i_err /= 0) then
         write(*,'(2x,a,1x,a,1x,i0)') 'Warning: could not open macro cell alist mapping file:', trim(out_file), i_err
         return
      end if
      write(ofileno,'(a)') '# i_cell n_atoms atom1 atom2 ...'

      do i_cell = 1, Num_macro
         n_neigh = min(macro_alistsize(i_cell), max_num_atom_macro_cell)
         write(ofileno,'(i8,1x,i8, 10x)', advance='no') i_cell, n_neigh
         do j = 1, n_neigh
            write(ofileno,'(1x,i8)', advance='no') macro_atom_alist(j, i_cell)
         end do
         write(ofileno,*)
      end do

      close(ofileno)
      write(*,'(2x,a,1x,a)') 'Wrote macro cell alist mapping:', trim(out_file)

   end subroutine print_macrocell_neighbour_map

   subroutine print_macro_halo_mapping(Num_macro, max_macro_halo_size, macro_halo_nlistsize, macro_halo_to_global, simid)
      implicit none
      integer, intent(in) :: Num_macro
      integer, intent(in) :: max_macro_halo_size
      integer, dimension(Num_macro), intent(in) :: macro_halo_nlistsize
      integer, dimension(max_macro_halo_size,Num_macro), intent(in) :: macro_halo_to_global
      character(len=8), intent(in) :: simid
      integer :: ofileno, i_cell, j, n_halo, i_err
      character(len=64) :: out_file

      out_file = 'macro_halo_mapping.' // trim(simid) // '.out'
      open(newunit=ofileno, file=out_file, status='replace', action='write', iostat=i_err)
      if (i_err /= 0) then
         write(*,'(2x,a,1x,a,1x,i0)') 'Warning: could not open macro halo mapping file:', trim(out_file), i_err
         return
      end if
      write(ofileno,'(a)') '# i_cell n_halo_atoms atom1 atom2 ...'
      do i_cell = 1, Num_macro
         n_halo = min(macro_halo_nlistsize(i_cell), max_macro_halo_size)
         write(ofileno,'(i8,1x,i8, 10x)', advance='no') i_cell, n_halo
         do j = 1, n_halo
            write(ofileno,'(1x,i8)', advance='no') macro_halo_to_global(j,i_cell)
         end do
         write(ofileno,*)
      end do
      close(ofileno)
      write(*,'(2x,a,1x,a)') 'Wrote macro halo mapping:', trim(out_file)
   end subroutine print_macro_halo_mapping

   subroutine print_macro_macro_mapping(Num_macro, max_macro_cell_neigh, macro_cell_nlistsize, macro_cell_nlist, simid)
      implicit none
      integer, intent(in) :: Num_macro
      integer, intent(in) :: max_macro_cell_neigh
      integer, dimension(Num_macro), intent(in) :: macro_cell_nlistsize
      integer, dimension(max_macro_cell_neigh,Num_macro), intent(in) :: macro_cell_nlist
      character(len=8), intent(in) :: simid
      integer :: ofileno, i_cell, j, n_cell, i_err
      character(len=64) :: out_file

      out_file = 'macro_macro_mapping.' // trim(simid) // '.out'
      open(newunit=ofileno, file=out_file, status='replace', action='write', iostat=i_err)
      if (i_err /= 0) then
         write(*,'(2x,a,1x,a,1x,i0)') 'Warning: could not open macro macro mapping file:', trim(out_file), i_err
         return
      end if
      write(ofileno,'(a)') '#######################################################'
      write(ofileno,'(a,1x,i8)') '# Number of macrocells:', Num_macro
      write(ofileno,'(a,1x,i8)') '# Maximum num of macro-neighbours:', max_macro_cell_neigh
      write(ofileno,'(a)') '#######################################################'
      write(ofileno,'(a)') '#  icell jcell'
      do i_cell = 1, Num_macro
         n_cell = min(macro_cell_nlistsize(i_cell), max_macro_cell_neigh)
         do j = 1, n_cell
            write(ofileno,'(2i8)') i_cell, macro_cell_nlist(j,i_cell)
         end do
      end do
      close(ofileno)
      write(*,'(2x,a,1x,a)') 'Wrote macro macro mapping:', trim(out_file)
   end subroutine print_macro_macro_mapping

   !----------------------------------------------------------------------------
   ! SUBROUTINE: init_macrocell
   !
   !> @brief Initialization routine for the default values for the macrocells.
   !
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine init_macrocell()

      implicit none

      do_macro_cells='N'
      prn_dip_subset='N'
      dip_file='dip_file.dat'
      Num_dip=0
      Num_macro=0
      max_num_atom_macro_cell=0
      max_macro_halo_size=0
      max_macro_cell_neigh=0
      max_macro_atom_local_neigh=0
      block_size_x=1
      block_size_y=1
      block_size_z=1

   end subroutine init_macrocell

   !----------------------------------------------------------------------------
   ! SUBROUTINE: create_macrocell
   !
   !> @brief Routine for the creation of cubic macro cells.
   !> @details It makes use of the geoblocking algorithm implemented by Anders Bergman.
   !> It also creates a series of helper arrays such as to identify which atom belongs
   !> to which macrocell, and which macro cell contains which atoms.
   !> It is intended to be used for the calculation of the macrocell dipolar interaction.
   !
   !> @author Jonathan Chico and Anders Bergman
   !> @todo generalize the routine to handle any kind of cell subdivision. Maybe
   !> some type of Voronoi kind of construction.
   !----------------------------------------------------------------------------
   subroutine create_macrocell(NA,N1,N2,N3,Natom,Mensemble,block_size,coord,&
         Num_macro,max_num_atom_macro_cell,cell_index,macro_alistsize,macro_atom_alist,simid)

      implicit none

      integer, intent(in) :: NA !< Number of atoms in one cell
      integer, intent(in) :: N1 !< Number of cell repetitions in x direction
      integer, intent(in) :: N2 !< Number of cell repetitions in y direction
      integer, intent(in) :: N3 !< Number of cell repetitions in z direction
      integer, intent(in) :: Natom        !< Number of atoms in system`
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      ! .. In/out arguments
      integer, intent(inout) :: Num_macro  !< Number of macrocells in the system
      integer, intent(inout) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(:), allocatable, intent(inout) :: cell_index !< Macrocell index for each atom
      integer, dimension(:), allocatable, intent(inout) :: macro_alistsize  !< Number of atoms per macrocell
      integer, dimension(:,:), allocatable, intent(inout) :: macro_atom_alist !< List containing the information of which atoms are in a given macrocell

      ! .. Local variables
      integer :: dim
      integer :: II1,II2,II3,I0,I1,I2,I3
      integer :: i_stat, i_all
      integer :: ii, kk
      integer :: bnx, bny, bnz, nbins, ibin
      integer :: nblkx, nblky, nblkz
      integer, dimension(:), allocatable :: bin_count, bin_to_macro, atom_bin
      real(dblprec) :: xmin, xmax, ymin, ymax, zmin, zmax
      real(dblprec) :: dx, dy, dz
      logical :: use_regular_partition

      character(len=23) :: output_file
      ! Creation of a file to writeout the geometry of the macro cell
      output_file = 'macro_info.'//simid//'.out'

      ii=0
      kk=0
      dim=0
      ! Calculate the dimensionality of the repetition of the unit cell
      if (N1>1) then
         dim=dim+1
      endif
      if (N2>1) then
         dim=dim+1
      endif
      if (N3>1) then
         dim=dim+1
      endif

      use_regular_partition=(N1*N2*N3*NA==Natom .and. N1>1 .and. N2>1 .and. N3>1)
      if (use_regular_partition) then
         ! Existing fast path for regular supercells.
         nblkx=(N1+block_size_x-1)/block_size_x
         nblky=(N2+block_size_y-1)/block_size_y
         nblkz=(N3+block_size_z-1)/block_size_z
         Num_macro=nblkx*nblky*nblkz
         max_num_atom_macro_cell=NA*(block_size_x*block_size_y*block_size_z)
         call allocate_macrocell(1,Natom,Mensemble)

         print *,'Number of macrocells:',Num_macro
         do II3=0, N3-1, block_size_z
            do II2=0, N2-1, block_size_y
               do II1=0, N1-1, block_size_x
                  kk=kk+1 ! Cell counter
                  do I3=II3,min(II3+block_size_z-1,N3-1)
                     do I2=II2,min(II2+block_size_y-1,N2-1)
                        do I1=II1,min(II1+block_size_x-1,N1-1)
                           do I0=1, NA
                              ii=ii+1 ! Atom counter
                              macro_alistsize(kk)=macro_alistsize(kk)+1
                              cell_index(ii)=kk
                              macro_atom_alist(macro_alistsize(kk),kk)=ii
                              atom_in_macro_pos(ii)=macro_alistsize(kk)
                              max_coord_macro(1,kk)=max(coord(1,ii),max_coord_macro(1,kk))
                              max_coord_macro(2,kk)=max(coord(2,ii),max_coord_macro(2,kk))
                              max_coord_macro(3,kk)=max(coord(3,ii),max_coord_macro(3,kk))
                              min_coord_macro(1,kk)=min(coord(1,ii),min_coord_macro(1,kk))
                              min_coord_macro(2,kk)=min(coord(2,ii),min_coord_macro(2,kk))
                              min_coord_macro(3,kk)=min(coord(3,ii),min_coord_macro(3,kk))
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         if (kk/=Num_macro) then
            write(*,'(2x,a,2(1x,i0))') 'Warning: expected/actual Num_macro mismatch:', Num_macro, kk
            Num_macro=kk
         end if
      else
         ! General path: coordinate-based spatial binning, independent of N1/N2/N3 ordering.
         bnx=max(1,block_size_x)
         bny=max(1,block_size_y)
         bnz=max(1,block_size_z)
         if (bnx==1 .and. bny==1 .and. bnz==1) then
            bnx=max(1,block_size)
            bny=max(1,block_size)
            bnz=max(1,block_size)
         end if
         nbins=bnx*bny*bnz

         allocate(bin_count(nbins),bin_to_macro(nbins),atom_bin(Natom),stat=i_stat)
         bin_count=0
         bin_to_macro=0
         atom_bin=0

         xmin=minval(coord(1,:)); xmax=maxval(coord(1,:))
         ymin=minval(coord(2,:)); ymax=maxval(coord(2,:))
         zmin=minval(coord(3,:)); zmax=maxval(coord(3,:))
         dx=max((xmax-xmin)/real(bnx,dblprec), dbl_tolerance)
         dy=max((ymax-ymin)/real(bny,dblprec), dbl_tolerance)
         dz=max((zmax-zmin)/real(bnz,dblprec), dbl_tolerance)

         do ii=1,Natom
            I1=min(max(1,int((coord(1,ii)-xmin)/dx)+1),bnx)
            I2=min(max(1,int((coord(2,ii)-ymin)/dy)+1),bny)
            I3=min(max(1,int((coord(3,ii)-zmin)/dz)+1),bnz)
            ibin=I1+(I2-1)*bnx+(I3-1)*bnx*bny
            atom_bin(ii)=ibin
            bin_count(ibin)=bin_count(ibin)+1
         end do

         Num_macro=0
         max_num_atom_macro_cell=0
         do ibin=1,nbins
            if (bin_count(ibin)>0) then
               Num_macro=Num_macro+1
               bin_to_macro(ibin)=Num_macro
               max_num_atom_macro_cell=max(max_num_atom_macro_cell,bin_count(ibin))
            end if
         end do
         call allocate_macrocell(1,Natom,Mensemble)
         print *,'Number of macrocells:',Num_macro

         do ii=1,Natom
            kk=bin_to_macro(atom_bin(ii))
            macro_alistsize(kk)=macro_alistsize(kk)+1
            cell_index(ii)=kk
            macro_atom_alist(macro_alistsize(kk),kk)=ii
            atom_in_macro_pos(ii)=macro_alistsize(kk)
            max_coord_macro(1,kk)=max(coord(1,ii),max_coord_macro(1,kk))
            max_coord_macro(2,kk)=max(coord(2,ii),max_coord_macro(2,kk))
            max_coord_macro(3,kk)=max(coord(3,ii),max_coord_macro(3,kk))
            min_coord_macro(1,kk)=min(coord(1,ii),min_coord_macro(1,kk))
            min_coord_macro(2,kk)=min(coord(2,ii),min_coord_macro(2,kk))
            min_coord_macro(3,kk)=min(coord(3,ii),min_coord_macro(3,kk))
         end do

         deallocate(bin_count,bin_to_macro,atom_bin,stat=i_stat)
      end if

      ! Print the midpoint of the macro cells
      open(ofileno,file=output_file)
      do kk=1, Num_macro
         mid_coord_macro(1,kk)=(max_coord_macro(1,kk)+min_coord_macro(1,kk))*0.5_dblprec
         mid_coord_macro(2,kk)=(max_coord_macro(2,kk)+min_coord_macro(2,kk))*0.5_dblprec
         mid_coord_macro(3,kk)=(max_coord_macro(3,kk)+min_coord_macro(3,kk))*0.5_dblprec
         write(ofileno,'(i6,3f16.8)') kk,mid_coord_macro(1,kk),mid_coord_macro(2,kk),mid_coord_macro(3,kk)
      enddo
      close(ofileno)

      ! Deallocate helper arrays of the macrocell coordinates
      i_all=-product(shape(min_coord_macro))*kind(min_coord_macro)
      deallocate(min_coord_macro,stat=i_stat)
      call memocc(i_stat,i_all,'min_coord_macro','create_macrocell')

      i_all=-product(shape(max_coord_macro))*kind(max_coord_macro)
      deallocate(max_coord_macro,stat=i_stat)
      call memocc(i_stat,i_all,'max_coord_macro','create_macrocell')

      i_all=-product(shape(mid_coord_macro))*kind(mid_coord_macro)
      deallocate(mid_coord_macro,stat=i_stat)
      call memocc(i_stat,i_all,'mid_coord_macro','create_macrocell')

      call print_macrocell_neighbour_map(Num_macro, max_num_atom_macro_cell,      &
            macro_alistsize, macro_atom_alist, simid)

   end subroutine create_macrocell

   subroutine build_macro_halo_maps(Natom,mnn,nlist,nlistsize,simid)
      implicit none
      integer, intent(in) :: Natom
      integer, intent(in) :: mnn
      integer, dimension(mnn,Natom), intent(in) :: nlist
      integer, dimension(Natom), intent(in) :: nlistsize
      character(len=8), intent(in) :: simid
      integer :: i, j, a, ii, jj, kk, i_stat, i_all, valid_count
      integer, allocatable :: mark(:), cell_mark(:), cell_touched(:)

      if (.not.allocated(macro_alistsize) .or. .not.allocated(macro_atom_alist)) return

      if (allocated(macro_halo_nlistsize)) then
         i_all=-product(shape(macro_halo_nlistsize))*kind(macro_halo_nlistsize)
         deallocate(macro_halo_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'macro_halo_nlistsize','build_macro_halo_maps')
      end if
      if (allocated(macro_halo_to_global)) then
         i_all=-product(shape(macro_halo_to_global))*kind(macro_halo_to_global)
         deallocate(macro_halo_to_global,stat=i_stat)
         call memocc(i_stat,i_all,'macro_halo_to_global','build_macro_halo_maps')
      end if
      if (allocated(macro_atom_local_nlistsize)) then
         i_all=-product(shape(macro_atom_local_nlistsize))*kind(macro_atom_local_nlistsize)
         deallocate(macro_atom_local_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'macro_atom_local_nlistsize','build_macro_halo_maps')
      end if
      if (allocated(macro_atom_local_nlist)) then
         i_all=-product(shape(macro_atom_local_nlist))*kind(macro_atom_local_nlist)
         deallocate(macro_atom_local_nlist,stat=i_stat)
         call memocc(i_stat,i_all,'macro_atom_local_nlist','build_macro_halo_maps')
      end if
      if (allocated(macro_cell_nlistsize)) then
         i_all=-product(shape(macro_cell_nlistsize))*kind(macro_cell_nlistsize)
         deallocate(macro_cell_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'macro_cell_nlistsize','build_macro_halo_maps')
      end if
      if (allocated(macro_cell_nlist)) then
         i_all=-product(shape(macro_cell_nlist))*kind(macro_cell_nlist)
         deallocate(macro_cell_nlist,stat=i_stat)
         call memocc(i_stat,i_all,'macro_cell_nlist','build_macro_halo_maps')
      end if

      allocate(mark(Natom), cell_mark(Num_macro), cell_touched(Num_macro), stat=i_stat)
      mark=0
      cell_mark=0
      max_macro_halo_size=0
      max_macro_cell_neigh=0

      allocate(macro_halo_nlistsize(Num_macro),stat=i_stat)
      call memocc(i_stat,product(shape(macro_halo_nlistsize))*kind(macro_halo_nlistsize),'macro_halo_nlistsize','build_macro_halo_maps')
      macro_halo_nlistsize=0

      do ii=1,Num_macro
         mark=0
         kk=0
         do a=1,macro_alistsize(ii)
            i=macro_atom_alist(a,ii)
            do j=1,nlistsize(i)
               jj=nlist(j,i)
               if (jj<1 .or. jj>Natom) cycle
               if (mark(jj)==0) then
                  mark(jj)=1
                  kk=kk+1
               end if
            end do
         end do
         macro_halo_nlistsize(ii)=kk
         max_macro_halo_size=max(max_macro_halo_size,kk)
      end do

      allocate(macro_halo_to_global(max(1,max_macro_halo_size),Num_macro),stat=i_stat)
      call memocc(i_stat,product(shape(macro_halo_to_global))*kind(macro_halo_to_global),'macro_halo_to_global','build_macro_halo_maps')
      macro_halo_to_global=0
      max_macro_atom_local_neigh=0
      do ii=1,Num_macro
         do a=1,macro_alistsize(ii)
            i=macro_atom_alist(a,ii)
            valid_count=0
            do j=1,nlistsize(i)
               jj=nlist(j,i)
               if (jj<1 .or. jj>Natom) cycle
               valid_count=valid_count+1
            end do
            max_macro_atom_local_neigh=max(max_macro_atom_local_neigh,valid_count)
         end do
      end do
      allocate(macro_atom_to_global(max_num_atom_macro_cell,Num_macro),stat=i_stat)
      call memocc(i_stat,product(shape(macro_atom_to_global))*kind(macro_atom_to_global),'macro_atom_to_global','build_macro_halo_maps')
      macro_atom_to_global=0
      allocate(macro_atom_local_nlistsize(max_num_atom_macro_cell,Num_macro),stat=i_stat)
      call memocc(i_stat,product(shape(macro_atom_local_nlistsize))*kind(macro_atom_local_nlistsize),'macro_atom_local_nlistsize','build_macro_halo_maps')
      macro_atom_local_nlistsize=0
      allocate(macro_atom_local_nlist(max(1,max_macro_atom_local_neigh),max_num_atom_macro_cell,Num_macro),stat=i_stat)
      call memocc(i_stat,product(shape(macro_atom_local_nlist))*kind(macro_atom_local_nlist),'macro_atom_local_nlist','build_macro_halo_maps')
      macro_atom_local_nlist=0

      do ii=1,Num_macro
         mark=0
         kk=0
         do a=1,macro_alistsize(ii)
            i=macro_atom_alist(a,ii)
            macro_atom_to_global(a,ii)=i
            do j=1,nlistsize(i)
               jj=nlist(j,i)
               if (jj<1 .or. jj>Natom) cycle
               if (mark(jj)==0) then
                  kk=kk+1
                  macro_halo_to_global(kk,ii)=jj
                  mark(jj)=kk
               end if
            end do
         end do
         do a=1,macro_alistsize(ii)
            i=macro_atom_alist(a,ii)
            kk=0
            do j=1,nlistsize(i)
               jj=nlist(j,i)
               if (jj<1 .or. jj>Natom) cycle
               kk=kk+1
               macro_atom_local_nlist(kk,a,ii)=mark(jj)
            end do
            macro_atom_local_nlistsize(a,ii)=kk
         end do

         cell_touched=0
         kk=0
         do j=1,macro_halo_nlistsize(ii)
            jj=macro_halo_to_global(j,ii)
            if (jj<1 .or. jj>Natom) cycle
            a=cell_index(jj)
            if (a<1 .or. a>Num_macro) cycle
            if (a/=ii .and. cell_mark(a)/=ii) then
               cell_mark(a)=ii
               kk=kk+1
               cell_touched(kk)=a
            end if
         end do
         max_macro_cell_neigh=max(max_macro_cell_neigh,kk)
      end do

      allocate(macro_cell_nlistsize(Num_macro),stat=i_stat)
      call memocc(i_stat,product(shape(macro_cell_nlistsize))*kind(macro_cell_nlistsize),'macro_cell_nlistsize','build_macro_halo_maps')
      macro_cell_nlistsize=0
      allocate(macro_cell_nlist(max(1,max_macro_cell_neigh),Num_macro),stat=i_stat)
      call memocc(i_stat,product(shape(macro_cell_nlist))*kind(macro_cell_nlist),'macro_cell_nlist','build_macro_halo_maps')
      macro_cell_nlist=0

      cell_mark=0
      do ii=1,Num_macro
         kk=0
         do j=1,macro_halo_nlistsize(ii)
            jj=macro_halo_to_global(j,ii)
            if (jj<1 .or. jj>Natom) cycle
            a=cell_index(jj)
            if (a<1 .or. a>Num_macro) cycle
            if (a/=ii .and. cell_mark(a)/=ii) then
               cell_mark(a)=ii
               kk=kk+1
               macro_cell_nlist(kk,ii)=a
            end if
         end do
         macro_cell_nlistsize(ii)=kk
      end do

      deallocate(mark,cell_mark,cell_touched,stat=i_stat)

      call print_macro_halo_mapping(Num_macro, max(1,max_macro_halo_size), macro_halo_nlistsize, macro_halo_to_global, simid)
      call print_macro_macro_mapping(Num_macro, max(1,max_macro_cell_neigh), macro_cell_nlistsize, macro_cell_nlist, simid)
   end subroutine build_macro_halo_maps

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_macro_mom
   !> @brief
   !> Calculation of the macro spin moment inside each macrocell.
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine calc_macro_mom(Natom,Num_macro,Mensemble,max_num_atom_macro_cell,&
         macro_alistsize,macro_atom_alist,mmom,emom,emomM,mmom_macro,emom_macro,emomM_macro)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system`
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_alistsize !< Number of atoms per macrocell
      integer, dimension(max_num_atom_macro_cell,Num_macro), intent(in) :: macro_atom_alist !< List containing the information of which atoms are in a given macrocell
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector

      real(dblprec), dimension(Num_macro,Mensemble), intent(inout) :: mmom_macro !< Magnitude of the macrocell magnetic moments
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro !< Unit vector of the macrocell magnetic moment
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro !< The full vector of the macrocell magnetic moment

      ! .. Local variables
      integer :: ii, jj, kk, iatom

      mmom_macro=0.0_dblprec
      emom_macro=0.0_dblprec
      emomM_macro=0.0_dblprec

      !$omp parallel do private(kk,ii,jj,iatom)
      do kk=1,Mensemble
         do ii=1, Num_macro
            do jj=1, macro_alistsize(ii)
               iatom=macro_atom_alist(jj,ii)
               mmom_macro(ii,kk)=mmom_macro(ii,kk)+mmom(iatom,kk)
               emom_macro(:,ii,kk)=emom_macro(:,ii,kk)+emom(:,iatom,kk)
               emomM_macro(:,ii,kk)=emomM_macro(:,ii,kk)+emomM(:,iatom,kk)
            enddo
            emom_macro(:,ii,kk)=emom_macro(:,ii,kk)/norm2(emom_macro(:,ii,kk))
         enddo
      enddo
      !$omp end parallel do

   end subroutine calc_macro_mom

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_trial_macro_mom
   !> @brief Calculation of the trial macrocell moment after change of direction
   !> of a magnetic moment in the Monte Carlo algorithm.
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine calc_trial_macro_mom(kk,iatom,Natom,Mensemble,Num_macro,max_num_atom_macro_cell,&
         macro_alistsize,macro_atom_alist,trialmom,emomM,emomM_macro,macro_mag_trial,macro_trial)

      implicit none

      integer, intent(in) :: kk ! Current ensemble
      integer, intent(in) :: iatom ! Current atom
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_alistsize !< Number of atoms per macrocell
      integer, dimension(max_num_atom_macro_cell,Num_macro), intent(in) :: macro_atom_alist !< List containing the information of which atoms are in a given macrocell
      real(dblprec), dimension(3), intent(in) :: trialmom   !< Trial moment from the montecarlo update
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM !< Current magnetic moment vector
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), intent(out) :: macro_mag_trial   !< Magnitude of the trial macrocell moment
      real(dblprec), dimension(3), intent(out) :: macro_trial  !< Full trial macrocell moment

      ! .. Local variables
      integer :: ii, icell,jatom

      macro_trial=0.0_dblprec

      icell=cell_index(iatom)
      do ii=1,macro_alistsize(icell)
         jatom=macro_atom_alist(ii,icell)
         if (jatom.eq.iatom) then
            macro_trial(:)=macro_trial(:)+trialmom(:)
         else
            macro_trial(:)=macro_trial(:)+emomM(:,jatom,kk)
         endif
      enddo

      macro_mag_trial=norm2(macro_trial)

   end subroutine calc_trial_macro_mom

   !----------------------------------------------------------------------------
   ! SUBROUTINE: read_dipole_subset
   !> @brief subroutine to read the subset of cells over which one will selectively
   !> print the stray fields from all the other cells.
   !> @author Jonatan Chico
   !----------------------------------------------------------------------------
   subroutine read_dipole_subset()

      implicit none

      ! .. Local variables
      integer :: flines,i_stat,iline,itemp,isite

      flines=0
      ! Open file and read the number of lines
      open(ifileno, file=trim(dip_file))
      do
         read(ifileno,*,end=200)  isite
         flines=flines+1
      end do
      200 continue
      ! Set the number of sites where the stray fields are going to be written
      Num_dip=flines

      allocate(dipole_subset(Num_dip),stat=i_stat)
      call memocc(i_stat,product(shape(dipole_subset))*kind(dipole_subset),'dipole_subset','read_dipole_subset')

      rewind(ifileno)
      ! Read the actual site where the dipolar field will be written
      do iline=1,Num_dip
         read(ifileno,*) itemp
         dipole_subset(iline)=itemp
      enddo

   end subroutine read_dipole_subset

   !----------------------------------------------------------------------------
   ! SUBROUTINE: calculate_dipole_subset
   !> @brief routine for writign the dipole-dipole field acting over a subset of
   !> cells in the macrocell dipole approach.
   !> @details The idea is to be able to "crop" a part of a bigger sample and print
   !> the dipolar field acting over that piece due to all the oter macro moments
   !> acting over those cells. This can then be used as an input for an external
   !> magnetic field.
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine calculate_dipole_subset(Natom,Mensemble,Num_macro,max_num_atom_macro_cell,&
      macro_alistsize,macro_atom_alist,emomM_macro,Qdip_macro,simid)

      use Constants

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_alistsize !< Number of atoms per macrocell
      integer, dimension(max_num_atom_macro_cell,Num_macro), intent(in) :: macro_atom_alist !< List containing the information of which atoms are in a given macrocell
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,3,Num_macro,Num_macro), intent(in) :: Qdip_macro !< Matrix for macro spin dipole-dipole interaction
      character(len=8), intent(in) :: simid !< Name of simulation

      ! .. Local variables
      integer :: ii,jj,icell,jcell
      real(dblprec) :: off_energy,fcinv
      character(len=23) :: crop_dip_out
      real(dblprec), dimension(3) :: field

      off_energy=0.0_dblprec

      ! Need to calculate the energy offset produced by the region of interest plus the outer-outer region
      do icell=1,Num_macro
         if (.not. ANY(dipole_subset(:)==icell)) then
            do jj=1,Num_dip
               jcell=dipole_subset(jj)
               off_energy = off_energy - Qdip_macro(1,1,jcell,icell)*emomM_macro(1,icell,1)*emomM_macro(1,jcell,1)
               off_energy = off_energy - Qdip_macro(2,1,jcell,icell)*emomM_macro(1,icell,1)*emomM_macro(2,jcell,1)
               off_energy = off_energy - Qdip_macro(3,1,jcell,icell)*emomM_macro(1,icell,1)*emomM_macro(3,jcell,1)
               off_energy = off_energy - Qdip_macro(1,2,jcell,icell)*emomM_macro(2,icell,1)*emomM_macro(1,jcell,1)
               off_energy = off_energy - Qdip_macro(2,2,jcell,icell)*emomM_macro(2,icell,1)*emomM_macro(2,jcell,1)
               off_energy = off_energy - Qdip_macro(3,2,jcell,icell)*emomM_macro(2,icell,1)*emomM_macro(3,jcell,1)
               off_energy = off_energy - Qdip_macro(1,3,jcell,icell)*emomM_macro(3,icell,1)*emomM_macro(1,jcell,1)
               off_energy = off_energy - Qdip_macro(2,3,jcell,icell)*emomM_macro(3,icell,1)*emomM_macro(2,jcell,1)
               off_energy = off_energy - Qdip_macro(3,3,jcell,icell)*emomM_macro(3,icell,1)*emomM_macro(3,jcell,1)
            enddo
            do jcell=1,Num_macro
               if (.not. ANY(dipole_subset(:)==jcell)) then
                  off_energy = off_energy - Qdip_macro(1,1,jcell,icell)*emomM_macro(1,icell,1)*emomM_macro(1,jcell,1)
                  off_energy = off_energy - Qdip_macro(2,1,jcell,icell)*emomM_macro(1,icell,1)*emomM_macro(2,jcell,1)
                  off_energy = off_energy - Qdip_macro(3,1,jcell,icell)*emomM_macro(1,icell,1)*emomM_macro(3,jcell,1)
                  off_energy = off_energy - Qdip_macro(1,2,jcell,icell)*emomM_macro(2,icell,1)*emomM_macro(1,jcell,1)
                  off_energy = off_energy - Qdip_macro(2,2,jcell,icell)*emomM_macro(2,icell,1)*emomM_macro(2,jcell,1)
                  off_energy = off_energy - Qdip_macro(3,2,jcell,icell)*emomM_macro(2,icell,1)*emomM_macro(3,jcell,1)
                  off_energy = off_energy - Qdip_macro(1,3,jcell,icell)*emomM_macro(3,icell,1)*emomM_macro(1,jcell,1)
                  off_energy = off_energy - Qdip_macro(2,3,jcell,icell)*emomM_macro(3,icell,1)*emomM_macro(2,jcell,1)
                  off_energy = off_energy - Qdip_macro(3,3,jcell,icell)*emomM_macro(3,icell,1)*emomM_macro(3,jcell,1)
               endif
            enddo
         endif
      enddo

      fcinv = mub*0.50_dblprec/mry

      crop_dip_out='crop_field.'//simid//'.out'
      open(ofileno,file=crop_dip_out)
      write(ofileno,'(a,2x,es16.8)') 'Ref. energy', off_energy*fcinv/Natom

      ! Cell of interest
      !$omp parallel do default(shared) schedule(static) private(jcell,icell,jj) reduction(+:field)
      do ii=1,Num_dip
         field=0.0_dblprec
         icell=dipole_subset(ii)
         do jcell=1, Num_macro
            ! If the current cell is not in the considered cells calculate the field
            if (.not. ANY(dipole_subset(:)==jcell)) then
               ! The field generated by ALL the other macrocells acting over the selected cell
               field(:)=field(:)+Qdip_macro(1,:,jcell,icell)*emomM_macro(1,jcell,1)&
               +Qdip_macro(2,:,jcell,icell)*emomM_macro(2,jcell,1)&
               +Qdip_macro(3,:,jcell,icell)*emomM_macro(3,jcell,1)
            endif
         enddo
         ! Print the field acting over each atom in the considered cell
         do jj=1,macro_alistsize(icell)
            write(ofileno,'(i5,x,es16.8,x,es16.8,x,es16.8)') macro_atom_alist(jj,icell), field(1),field(2),field(3)
         enddo
      enddo
      !$omp end parallel do
      close(ofileno)
   end subroutine calculate_dipole_subset

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_macrocell
   !> @brief
   !> Helper routine to allocate/deallocate arrays needed for the macrocell approach.
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_macrocell(flag,Natom,Mensemble)

      implicit none

      integer, intent(in) :: flag
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles

      ! .. Local variables
      integer :: i_stat, i_all

      if (flag>0) then
         allocate(macro_alistsize(Num_macro), stat=i_stat)
         call memocc(i_stat,product(shape(macro_alistsize))*kind(macro_alistsize),'macro_alistsize','allocate_macrocell')
         macro_alistsize=0
         allocate(cell_index(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(cell_index))*kind(cell_index),'cell_index','allocate_macrocell')
         cell_index=0
         allocate(macro_atom_alist(max_num_atom_macro_cell,Num_macro),stat=i_stat)
         call memocc(i_stat,product(shape(macro_atom_alist))*kind(macro_atom_alist),'macro_atom_alist','allocate_macrocell')
         macro_atom_alist=0
         allocate(max_coord_macro(3,Num_macro),stat=i_stat)
         call memocc(i_stat,product(shape(max_coord_macro))*kind(max_coord_macro),'max_coord_macro','allocate_macrocell')
         max_coord_macro=0.0_dblprec
         allocate(min_coord_macro(3,Num_macro),stat=i_stat)
         call memocc(i_stat,product(shape(min_coord_macro))*kind(min_coord_macro),'min_coord_macro','allocate_macrocell')
         min_coord_macro=1.0d9
         allocate(mid_coord_macro(3,Num_macro),stat=i_stat)
         call memocc(i_stat,product(shape(mid_coord_macro))*kind(mid_coord_macro),'mid_coord_macro','allocate_macrocell')
         mid_coord_macro=0.0_dblprec
         allocate(mmom_macro(Num_macro,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mmom_macro))*kind(mmom_macro),'mmom_macro','allocate_macrocell')
         mmom_macro=0.0_dblprec
         allocate(emom_macro(3,Num_macro,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(emom_macro))*kind(emom_macro),'emom_macro','allocate_macrocell')
         emom_macro=0.0_dblprec
         allocate(emomM_macro(3,Num_macro,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(emomM_macro))*kind(emomM_macro),'emomM_macro','allocate_macrocell')
         emomM_macro=0.0_dblprec
         allocate(atom_in_macro_pos(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(atom_in_macro_pos))*kind(atom_in_macro_pos),'atom_in_macro_pos','allocate_macrocell')
         atom_in_macro_pos=0
      else

         i_all=-product(shape(cell_index))*kind(cell_index)
         deallocate(cell_index,stat=i_stat)
         call memocc(i_stat,i_all,'cell_index','allocate_macrocell')

         i_all=-product(shape(macro_atom_alist))*kind(macro_atom_alist)
         deallocate(macro_atom_alist,stat=i_stat)
         call memocc(i_stat,i_all,'macro_atom_alist','allocate_macrocell')

         i_all=-product(shape(macro_alistsize))*kind(macro_alistsize)
         deallocate(macro_alistsize,stat=i_stat)
         call memocc(i_stat,i_all,'macro_alistsize','allocate_macrocell')
         if (allocated(atom_in_macro_pos)) then
            i_all=-product(shape(atom_in_macro_pos))*kind(atom_in_macro_pos)
            deallocate(atom_in_macro_pos,stat=i_stat)
            call memocc(i_stat,i_all,'atom_in_macro_pos','allocate_macrocell')
         end if

         i_all=-product(shape(mmom_macro))*kind(mmom_macro)
         deallocate(mmom_macro,stat=i_stat)
         call memocc(i_stat,i_all,'mmom_macro','allocate_macrocell')

         i_all=-product(shape(emom_macro))*kind(emom_macro)
         deallocate(emom_macro,stat=i_stat)
         call memocc(i_stat,i_all,'emom_macro','allocate_macrocell')

         i_all=-product(shape(emomM_macro))*kind(emomM_macro)
         deallocate(emomM_macro,stat=i_stat)
         call memocc(i_stat,i_all,'emomM_macro','allocate_macrocell')

         if (allocated(dipole_subset)) then
            i_all=-product(shape(dipole_subset))*kind(dipole_subset)
            deallocate(dipole_subset,stat=i_stat)
            call memocc(i_stat,i_all,'dipole_subset','allocate_macrocell')
         endif
         if (allocated(macro_halo_nlistsize)) then
            i_all=-product(shape(macro_halo_nlistsize))*kind(macro_halo_nlistsize)
            deallocate(macro_halo_nlistsize,stat=i_stat)
            call memocc(i_stat,i_all,'macro_halo_nlistsize','allocate_macrocell')
         endif
         if (allocated(macro_halo_to_global)) then
            i_all=-product(shape(macro_halo_to_global))*kind(macro_halo_to_global)
            deallocate(macro_halo_to_global,stat=i_stat)
            call memocc(i_stat,i_all,'macro_halo_to_global','allocate_macrocell')
         endif
         if (allocated(macro_atom_local_nlistsize)) then
            i_all=-product(shape(macro_atom_local_nlistsize))*kind(macro_atom_local_nlistsize)
            deallocate(macro_atom_local_nlistsize,stat=i_stat)
            call memocc(i_stat,i_all,'macro_atom_local_nlistsize','allocate_macrocell')
         endif
         if (allocated(macro_atom_to_global)) then
            i_all=-product(shape(macro_atom_to_global))*kind(macro_atom_to_global)
            deallocate(macro_atom_to_global,stat=i_stat)
            call memocc(i_stat,i_all,'macro_atom_to_global','allocate_macrocell')
         endif
         if (allocated(macro_atom_local_nlist)) then
            i_all=-product(shape(macro_atom_local_nlist))*kind(macro_atom_local_nlist)
            deallocate(macro_atom_local_nlist,stat=i_stat)
            call memocc(i_stat,i_all,'macro_atom_local_nlist','allocate_macrocell')
         endif
         if (allocated(macro_cell_nlistsize)) then
            i_all=-product(shape(macro_cell_nlistsize))*kind(macro_cell_nlistsize)
            deallocate(macro_cell_nlistsize,stat=i_stat)
            call memocc(i_stat,i_all,'macro_cell_nlistsize','allocate_macrocell')
         endif
         if (allocated(macro_cell_nlist)) then
            i_all=-product(shape(macro_cell_nlist))*kind(macro_cell_nlist)
            deallocate(macro_cell_nlist,stat=i_stat)
            call memocc(i_stat,i_all,'macro_cell_nlist','allocate_macrocell')
         endif
      endif

   end subroutine allocate_macrocell

end module macrocells
