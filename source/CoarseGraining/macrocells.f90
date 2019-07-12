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
   integer, dimension(:), allocatable :: macro_nlistsize !< Number of atoms per macrocell
   integer, dimension(:,:), allocatable :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell

   real(dblprec), dimension(:,:), allocatable :: mmom_macro !< Magnitude of the macrocell magnetic moments
   real(dblprec), dimension(:,:), allocatable :: max_coord_macro !< Maximum value of the coordinates per cell
   real(dblprec), dimension(:,:), allocatable :: min_coord_macro !< Minimum value of the coordinates per cell
   real(dblprec), dimension(:,:), allocatable :: mid_coord_macro !< Midpoint of each cell
   real(dblprec), dimension(:,:,:), allocatable :: emom_macro  !< Unit vector of the macrocell magnetic moment
   real(dblprec), dimension(:,:,:), allocatable :: emomM_macro !< The full vector of the macrocell magnetic moment

contains

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
      Num_macro=0
      max_num_atom_macro_cell=0

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
         Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,macro_atom_nlist,simid)

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
      integer, dimension(:), allocatable, intent(inout) :: macro_nlistsize  !< Number of atoms per macrocell
      integer, dimension(:,:), allocatable, intent(inout) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell

      ! .. Local variables
      integer :: dim
      integer :: ii, kk
      integer :: II1,II2,II3,I0,I1,I2,I3
      integer :: i_stat, i_all

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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Notice that the determination of the following parameters assumes a cubic
      ! macrocell, this should be modified for a general shape macrocell
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the number of macro cells
      Num_macro=int(N3*N2*N1/block_size**dim)
      ! Calculate the Maximum number of atoms per macro cell
      max_num_atom_macro_cell=NA*block_size**dim
      call allocate_macrocell(1,Natom,Mensemble)

      ! Create the macrocells lists needed for the macrocell approximation
      do II3=0, N3-1, block_size
         do II2=0, N2-1, block_size
            do II1=0, N1-1, block_size
               kk=kk+1 ! Cell counter
               do I3=II3,min(II3+block_size-1,N3-1)
                  do I2=II2,min(II2+block_size-1,N2-1)
                     do I1=II1,min(II1+block_size-1,N1-1)
                        do I0=1, NA
                           ii=ii+1 ! Atom counter
                           macro_nlistsize(kk)=macro_nlistsize(kk)+1
                           cell_index(ii)=kk
                           macro_atom_nlist(kk,macro_nlistsize(kk))=ii
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

   end subroutine create_macrocell

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_macro_mom
   !> @brief
   !> Calculation of the macro spin moment inside each macrocell.
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine calc_macro_mom(Natom,Num_macro,Mensemble,max_num_atom_macro_cell,&
         macro_nlistsize,macro_atom_nlist,mmom,emom,emomM,mmom_macro,emom_macro,emomM_macro)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system`
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
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
            do jj=1, macro_nlistsize(ii)
               iatom=macro_atom_nlist(ii,jj)
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
         macro_nlistsize,macro_atom_nlist,trialmom,emomM,emomM_macro,macro_mag_trial,macro_trial)

      implicit none

      integer, intent(in) :: kk ! Current ensemble
      integer, intent(in) :: iatom ! Current atom
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
      real(dblprec), dimension(3), intent(in) :: trialmom   !< Trial moment from the montecarlo update
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM !< Current magnetic moment vector
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), intent(out) :: macro_mag_trial   !< Magnitude of the trial macrocell moment
      real(dblprec), dimension(3), intent(out) :: macro_trial  !< Full trial macrocell moment

      ! .. Local variables
      integer :: ii, icell,jatom

      macro_trial=0.0_dblprec

      icell=cell_index(iatom)
      do ii=1,macro_nlistsize(icell)
         jatom=macro_atom_nlist(icell,ii)
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
      macro_nlistsize,macro_atom_nlist,emomM_macro,Qdip_macro,simid)

      use Constants

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
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
         do jj=1,macro_nlistsize(icell)
            write(ofileno,'(i5,x,es16.8,x,es16.8,x,es16.8)') macro_atom_nlist(icell,jj), field(1),field(2),field(3)
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
         allocate(macro_nlistsize(Num_macro), stat=i_stat)
         call memocc(i_stat,product(shape(macro_nlistsize))*kind(macro_nlistsize),'macro_nlistsize','allocate_macrocell')
         macro_nlistsize=0
         allocate(cell_index(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(cell_index))*kind(cell_index),'cell_index','allocate_macrocell')
         cell_index=0
         allocate(macro_atom_nlist(Num_macro,max_num_atom_macro_cell),stat=i_stat)
         call memocc(i_stat,product(shape(macro_atom_nlist))*kind(macro_atom_nlist),'macro_atom_nlist','allocate_macrocell')
         macro_atom_nlist=0
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
      else

         i_all=-product(shape(cell_index))*kind(cell_index)
         deallocate(cell_index,stat=i_stat)
         call memocc(i_stat,i_all,'cell_index','allocate_macrocell')

         i_all=-product(shape(macro_atom_nlist))*kind(macro_atom_nlist)
         deallocate(macro_atom_nlist,stat=i_stat)
         call memocc(i_stat,i_all,'macro_atom_nlist','allocate_macrocell')

         i_all=-product(shape(macro_nlistsize))*kind(macro_nlistsize)
         deallocate(macro_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'macro_nlistsize','allocate_macrocell')

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
      endif

   end subroutine allocate_macrocell

end module macrocells
