!-------------------------------------------------------------------------------
!  MODULE: FixedMom
!> @brief
!> Intended to describe the selective update of magnetic moments.
!> @details This routine contains the necessary data structures such that only
!> specific magnetic moments are updated, while the rest is kept fixed
!> Currently only working with Depondt solver SDEalgh=5
!> @author
!> Jonathan Chico
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module FixedMom
   use Parameters
   use Profiling

   implicit none

   integer :: Nred
   integer, dimension(:), allocatable :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
   integer, dimension(:), allocatable :: fixed_mom_flag !< Site dependent fixed moment flag
   integer, dimension(:,:,:), allocatable :: inp_fixed_mom_flag !< Fixed moments flag defined for atoms in the unit cell
   character(len=1) :: do_fixed_mom !< Do Fixed moment calculation (Y/N)

   public

contains

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: init_fixed_mom
   !> @brief
   !> Set the default values of the fixed moment approach
   !
   !> @author
   !> Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine init_fixed_mom()

      implicit none

      do_fixed_mom='N'
      Nred=0

   end subroutine init_fixed_mom

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: fixed_loadrestart
   !> @brief
   !> Read magnetic moments from file for the case of fixed moments
   !
   !> @author
   !> Jonathan Chico, based on the previously existent loadrestart routine
   !-----------------------------------------------------------------------------
   subroutine fixed_loadrestart(Natom,Mensemble,restartfile,rstep,mmom,emom,emomM,&
         fixed_mom_flag,Nred,ind_mom_flag,ind_list_full)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=1), intent(in) :: ind_mom_flag !< Flag to indicate that there are induced moments being considered
      integer, intent(out) :: rstep !< Starting simulation step
      integer, intent(inout) :: Nred !< Number of moments that can be updated
      integer, dimension(Natom), intent(inout) :: ind_list_full !< Indication of whether a given moment is induced/fixed 1/0
      integer, dimension(Natom), intent(inout) :: fixed_mom_flag !< Site dependent fixed moment flag
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(out) :: mmom !< Magnitude of magnetic moments
      character(len=35), intent(inout) :: restartfile !< File containing restart information

      integer :: i, j, k, l, ios,ired
      logical :: exists

      ired=0
      !.. Executable statements
      inquire(file=restartfile,exist=exists)
      if(exists) then
         open(ifileno,iostat=ios, file=restartfile, status="old")
         read (ifileno,*) rstep
         if (ind_mom_flag=='Y') then
            do i=1,Mensemble
               ired=0
               do j=1, Natom
                  read (ifileno,*) k, l, mmom(j,i), emom(1,j,i), emom(2,j,i), emom(3,j,i),fixed_mom_flag(j),ind_list_full(j)
                  if (fixed_mom_flag(j).eq.1) then
                     ired=ired+1
                  endif
                  emomM(:,j,i)=emom(:,j,i)*mmom(j,i)
               end do
            end do
            close(ifileno)
         else
            do i=1,Mensemble
               ired=0
               do j=1, Natom
                  read (ifileno,*) k, l, mmom(j,i), emom(1,j,i), emom(2,j,i), emom(3,j,i),fixed_mom_flag(j)
                  if (fixed_mom_flag(j).eq.1) then
                     ired=ired+1
                  endif
                  emomM(:,j,i)=emom(:,j,i)*mmom(j,i)
               end do
            end do
            close(ifileno)
         endif
      else
         write(*,*) 'ERROR: Restartfile ',trim(adjustl(restartfile)), ' does not exist.'
         stop
      end if

      Nred=ired
   end subroutine fixed_loadrestart

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: create_aux_fixed_arrays
   !> @brief
   !> Create the necessary auxiliary arrays for the fixed moments calculations
   !
   !> @author
   !> Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine create_aux_fixed_arrays(Natom,Mensemble,initmag,Natom_full,do_ralloy, &
      anumb,achtype,rstep,red_atom_list,Nred,mmom,emom,emomM,restartfile,           &
      ind_mom_flag,ind_list_full,do_mom_legacy)

      use Restart, only: read_mag_conf

      implicit none

      ! .. Input variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: initmag !< Mode of initialization of magnetic moments (1-4)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      integer, dimension(:), allocatable, intent(in) :: achtype !< Chemical type of atoms (full list)
      character(len=1), intent(in) :: ind_mom_flag !< Flag to indicate that there are induced moments being considered
      character(len=1), intent(in) :: do_mom_legacy
      ! .. Output variables
      integer, intent(out) :: rstep !< Starting simulation step
      integer, dimension(:), allocatable, intent(out) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      ! .. In/out variables
      integer, intent(out) :: Nred !< Number of moments that can be updated
      integer, dimension(Natom), intent(inout) :: ind_list_full !< Indication of whether a given moment is induced/fixed 1/0
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom    !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM !< Current magnetic moment vector
      character(len=35), intent(inout) :: restartfile !< File containing restart information

      ! .. Local variables
      integer :: i_stat, i_all
      integer :: ired, iatom

      allocate(fixed_mom_flag(Natom),stat=i_stat)
      call memocc(i_stat,product(shape(fixed_mom_flag))*kind(fixed_mom_flag),'fixed_mom_flag','create_aux_fixed_arrays')
      fixed_mom_flag=0

      ired=0

      ! Fill the flag for fixed magnetic moments for each atom
      ! If the initial magnetizaion is given by the restartfile
      if (initmag.eq.4) then
         if (ind_mom_flag=='Y') then
            call read_mag_conf(Natom,Mensemble,do_mom_legacy,rstep,restartfile,     &
               mmom,emom,emomM,fixed_mom_flag,ind_list_full,Nred)
         else
            call read_mag_conf(Natom,Mensemble,do_mom_legacy,rstep,restartfile, &
            mmom,emom,emomM,fixed_mom_flag,Nred)
         endif
         ! If there is some other restartfile, use the information of the momfile
      else
         ! Do a loop over the atoms in the system
         if (do_ralloy.eq.0) then
            do iatom=1, Natom
               if (inp_fixed_mom_flag(anumb(iatom),1,1).eq.1) then
                  ired=ired+1
                  fixed_mom_flag(iatom)=1
               else
                  fixed_mom_flag(iatom)=0
               endif
            enddo
    else if ((do_ralloy==1.or.do_ralloy==2)) then
            do iatom=1, Natom
               if (inp_fixed_mom_flag(anumb(iatom),achtype(iatom),1).eq.1) then
                  ired=ired+1
                  fixed_mom_flag(iatom)=1
               else
                  fixed_mom_flag(iatom)=0
               endif
            enddo
         endif
         Nred=ired
      endif

      ! allocate the reduced atom list
      allocate(red_atom_list(Nred),stat=i_stat)
      call memocc(i_stat,product(shape(red_atom_list))*kind(red_atom_list),'red_atom_list','create_aux_fixed_arrays')

      ired=0
      do iatom=1,Natom
         if(fixed_mom_flag(iatom).eq.1) then
            ired=ired+1
            red_atom_list(ired)=iatom
         endif
      enddo

      ! Deallocate the helper array containing wheter each site gets updated
      i_all=-product(shape(fixed_mom_flag))*kind(fixed_mom_flag)
      deallocate(fixed_mom_flag,stat=i_stat)
      call memocc(i_stat,i_all,'fixed_mom_flag','create_aux_fixed_arrays')
      !
      if (allocated(inp_fixed_mom_flag)) then
         i_all=-product(shape(inp_fixed_mom_flag))*kind(inp_fixed_mom_flag)
         deallocate(inp_fixed_mom_flag,stat=i_stat)
         call memocc(i_stat,i_all,'inp_fixed_mom_flag','deallocate_rest')
      endif

   end subroutine create_aux_fixed_arrays

end module FixedMom
