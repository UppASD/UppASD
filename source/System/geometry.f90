!------------------------------------------------------------------------------------
! MODULE: Geometry
!> @brief Routines for setting up structural and chemical information about the system
!> @author
!> Anders Bergman, Johan Hellsvik
!> @copyright
!> GNU Public License.
!------------------------------------------------------------------------------------
module Geometry
   use Parameters
   use Profiling
   use table_tet_mesh, only : tet_mesh_wrapper

   implicit none

   private
   public :: setup_geometry,setup_type_and_numb,setup_globcoord,setup_chemicaldata,rescale_lattvec

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: setup_geometry
   !> @brief Wrapper routine for setting up structural and chemical information about the system
   !---------------------------------------------------------------------------------
   subroutine setup_geometry(Natom,NA,N1,N2,N3,Bas,C1,C2,C3,coord,atype,anumb,      &
      atype_inp,anumb_inp,do_prnstruct,do_prn_poscar,tseed,simid,do_ralloy,         &
      Natom_full,Nchmax,Nch,acellnumb,acellnumbrev,achtype,chconc,atype_ch,asite_ch,&
      achem_ch,block_size)

      implicit none

      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: tseed        !< Temperature seed
      integer, intent(in) :: Nchmax       !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: do_prnstruct !< Print Hamiltonian information (0/1)
      integer, intent(in) :: do_prn_poscar    !< Print geometry on POSCAR format (0/1)
      character(len=8), intent(in) :: simid     !< Name of simulation
      integer, dimension(:), intent(in) :: Nch  !< Number of chemical components on each site in cell
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(:,:), intent(in) :: chconc !< Chemical concentration on sites
      ! .. In/Out variables
      integer, intent(inout) :: Natom !< Number of atoms in system
      integer, dimension(Natom_full), intent(inout)   :: atype       !< Type of atom
      integer, dimension(Natom_full), intent(inout)   :: anumb       !< Atom number in cell
      integer, dimension(NA), intent(inout)           :: atype_inp   !< Type of atom from input
      integer, dimension(NA), intent(inout)           :: anumb_inp   !< Atom number in cell from input
      real(dblprec), dimension(3,NA) , intent(inout) :: Bas !< Coordinates for basis atoms
      real(dblprec), dimension(:,:), allocatable, intent(inout) :: coord !< Coordinates of atoms
      ! .. Output variables
      integer, dimension(:), intent(out) :: achtype      !< Chemical type of atoms (full list)
      integer, dimension(:), intent(out) :: atype_ch     !< Actual type of atom for dilute system
      integer, dimension(:), intent(out) :: asite_ch     !< Actual site of atom for dilute system
      integer, dimension(:), intent(out) :: achem_ch     !< Chemical type of atoms (reduced list)
      integer, dimension(:), intent(out) :: acellnumb    !< List for translating atom no. in full cell to actual cell
      integer, dimension(:), intent(out) :: acellnumbrev !< List for translating atom no. in actual cell to full cell

      write (*,'(2x,a)',advance='no') 'Set up types of atoms'
      call setup_type_and_numb(Natom_full,NA,N1,N2,N3,atype,anumb,atype_inp,        &
         anumb_inp,block_size)
      write (*,'(a)') ' done'

      write (*,'(2x,a)',advance='no') 'Set up chemical information of alloy'
      if(do_ralloy==1.or.do_ralloy==2) then
         call setup_chemicaldata(Natom,NA,N1,N2,N3,atype,tseed,do_ralloy,Natom_full,&
            Nchmax,Nch,achtype,acellnumb,acellnumbrev,chconc,atype_ch,asite_ch,     &
            achem_ch)
      else
         Natom=Natom_full
      end if
      write (*,'(a)') ' done'

      write (*,'(2x,a)',advance='no') 'Set up global coordinates'
      call setup_globcoord(Natom,NA,Bas,C1,C2,C3,coord,N1,N2,N3,atype,anumb,        &
         do_prnstruct,do_prn_poscar,simid,do_ralloy,Natom_full,acellnumb,achtype,   &
         atype_ch,asite_ch,achem_ch,block_size)
      write (*,'(a)') ' done'

      ! If the system is 3D calculate the Delaunay tesellation
      if (calc_dim(coord,Natom)=='3D') then
         write(*,'(2x,a)',advance='no') 'Perform Delaunay tesellation'
         call tet_mesh_wrapper(Natom,coord,simid)
         write(*,'(a)') 'done'
      endif

   end subroutine setup_geometry

   !---------------------------------------------------------------------------------
   ! FUNCTION: calc_dim
   !> @brief Calculate the dimensionality of the sample
   !> @details This is a function that determines the dimensionality of the current
   !> sample, it then specifies whether the dimensionality is 1D,2D or 3D
   !---------------------------------------------------------------------------------
   function calc_dim(coord,Natom) result(dimensionality)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atomsd
      ! .. Local variables
      real(dblprec) :: tol
      real(dblprec) :: min_x,min_y,min_z
      real(dblprec) :: max_x,max_y,max_z
      real(dblprec) :: dist_x,dist_y,dist_z
      ! Result
      character(len=3) :: dimensionality 
      
      tol=0.00050_dblprec**2
      ! Find the extremun of the coordinates
      max_x=maxval(coord(1,:))
      min_x=maxval(coord(1,:))
      max_y=maxval(coord(2,:))
      min_y=maxval(coord(2,:))
      max_z=maxval(coord(3,:))
      min_z=maxval(coord(3,:))
      ! Calculate the distance 
      dist_x=abs(max_x-min_x)
      dist_y=abs(max_y-min_y)
      dist_z=abs(max_z-min_z)

      if ((dist_x.le.tol).or.(dist_y.le.tol).or.(dist_z.le.tol)) then
         dimensionality='2D'
         if (((dist_x.le.tol).and.(dist_y.le.tol)).or.                              &
             ((dist_x.le.tol).and.(dist_z.le.tol)).or.                              &
             ((dist_y.le.tol).and.(dist_z.le.tol))) then
            dimensionality='1D'
         endif
      else
         dimensionality='3D'
      endif

   end function calc_dim
   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_type_and_numb
   !> @brief Sets up the type of atoms in the system
   !----------------------------------------------------------------------------
   subroutine setup_type_and_numb(Natom_full,NA,N1,N2,N3,atype,anumb,atype_inp,     &
      anumb_inp,block_size)
      !
      implicit none

      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      integer, intent(in) :: Natom_full   !< Number of atoms in system
      ! .. In/Out variables
      integer, dimension(Natom_full), intent(inout)   :: atype !< Type of atom
      integer, dimension(Natom_full), intent(inout)   :: anumb !< Atom number in cell
      integer, dimension(NA), intent(inout)           :: atype_inp  !< Type of atom from input
      integer, dimension(NA), intent(inout)           :: anumb_inp !< Atom number in cell from input
      !
      integer :: A1
      integer :: i0, i1, i2, i3, ii1, ii2, ii3
      !
      A1=0
      do iI3=0, N3-1, block_size
         do iI2=0, N2-1, block_size
            do iI1=0, N1-1, block_size
               do I3=ii3,min(ii3+block_size-1,N3-1)
                  do I2=ii2,min(ii2+block_size-1,N2-1)
                     do I1=ii1,min(ii1+block_size-1,N1-1)
                        do I0=1, NA
                           A1=A1+1
                           atype(A1)=atype_inp(I0)
                           anumb(A1)=anumb_inp(I0)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine setup_type_and_numb

   !----------------------------------------------------------------------------
   !> @brief Sets up chemical information of the system in case of random alloy
   !> @note IMPORTANT: In case of dilute system, Natom variable is reduced.
   !----------------------------------------------------------------------------
   subroutine setup_chemicaldata(Natom,NA,N1,N2,N3,atype,tseed,do_ralloy,Natom_full,&
      Nchmax,Nch,achtype,acellnumb,acellnumbrev,chconc,atype_ch,asite_ch,achem_ch)
      !
      use RandomNumbers, only : rng_uniform
      !
      implicit none

      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: tseed        !< Temperature seed
      integer, intent(in) :: Nchmax       !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, dimension(NA), intent(in) :: Nch !< Number of chemical components on each site in cell
      real(dblprec), dimension(NA,Nchmax), intent(in) :: chconc !< Chemical concentration on sites
      ! .. Output variables
      integer, intent(out) :: Natom !< Number of atoms in system
      integer, dimension(Natom_full), intent(out) :: achtype      !< Chemical type of atoms (full list)
      integer, dimension(Natom_full), intent(out) :: atype_ch     !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(out) :: asite_ch     !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(out) :: achem_ch     !< Chemical type of atoms (reduced list)
      integer, dimension(Natom_full), intent(out) :: acellnumb    !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full), intent(out) :: acellnumbrev !< List for translating atom no. in actual cell to full cell

      ! .. In/Out variables
      integer, dimension(Natom_full), intent(inout) :: atype !< Type of atom

      integer :: A1
      integer :: i0, i1, i2, i3
      integer :: Ncell, iatom
      integer :: ia, ich, i
      integer :: i_stat, i_all
      integer, dimension(1) :: irn
      integer, dimension(:,:), allocatable :: atoms
      integer, dimension(:,:), allocatable :: qch
      integer, dimension(:), allocatable :: ns, ne
      !integer, dimension(:), allocatable :: atoms2, atoms2T
      real(dblprec), dimension(:), allocatable :: rn

      if (do_ralloy==1.or.do_ralloy==2) then

         ! Same seed as for Monte-Carlo
         Ncell = N1*N2*N3

         ! For rnd planes version
         if  (do_ralloy==2) then
                 Ncell = N3
         endif

         allocate(atoms(Ncell,Na),stat=i_stat)
         call memocc(i_stat,product(shape(atoms))*kind(atoms),'atoms','setup_chemicaldata')
         ! Not used
         !allocate(atoms2(Ncell),stat=i_stat)
         !call memocc(i_stat,product(shape(atoms2))*kind(atoms2),'atoms2','setup_chemicaldata')
         !allocate(atoms2T(Ncell),stat=i_stat)
         !call memocc(i_stat,product(shape(atoms2T))*kind(atoms2T),'atoms2T','setup_chemicaldata')
         allocate(qch(Nchmax,Na),stat=i_stat)
         call memocc(i_stat,product(shape(qch))*kind(qch),'qch','setup_chemicaldata')
         allocate(rn(Ncell),stat=i_stat)
         call memocc(i_stat,product(shape(rn))*kind(rn),'rn','setup_chemicaldata')
         allocate(ns(Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ns))*kind(ns),'ns','setup_chemicaldata')
         allocate(ne(Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ne))*kind(ne),'ne','setup_chemicaldata')
         qch=0
         atoms=0
         achtype=0

         do ia=1,Na
            do ich=1,Nch(ia)
               qch(ich,ia) = nint(chconc(ia,ich)*Ncell)
            end do
         end do

         do ia=1,Na
            call rng_uniform(rn,Ncell)
            !Partitioning of the Ncell random nr:s
            !to sets for each chemical type ich
            ns(1) = 1
            ne(1) = qch(1,ia)
            do ich=2,Nch(ia)
               ns(ich) = ns(ich-1)+qch(ich-1,ia)
               ne(ich) = ne(ich-1)+qch(ich,ia)
            end do
            do i=1,Ncell
               irn = maxloc(rn)
               atoms(i,ia)=irn(1)
               rn(irn)=0.0_dblprec
            end do
            do ich=1,Nch(ia)
               do i=ns(ich),ne(ich)
                  A1=(atoms(i,ia)-1)*Na+ia
                  achtype(A1)=ich
               end do
            end do
         end do
         acellnumb=0
         acellnumbrev=0
         iatom=0
         ! asite_ch and achem_ch contains information of species on each site
         do I3=0, N3-1
            do I2=0, N2-1
               do I1=0, N1-1
                  do I0=1, NA
                     i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                     if (achtype(i) /= 0) then
                        iatom = iatom + 1
                        acellnumb(i) = iatom
                        acellnumbrev(iatom) = i
                        asite_ch(iatom)=I0
                        atype_ch(iatom)=atype(i)
                        if  (do_ralloy==2) then
                                achem_ch(iatom)=achtype(I3)
                        else
                                achem_ch(iatom)=achtype(i)
                        end if
                     end if
                  end do
               end do
            end do
         end do

         ! Natom is reduced in case of dilute system
         Natom = iatom
         ! Deallocate temporary arrays
         i_all=-product(shape(atoms))*kind(atoms)
         deallocate(atoms,stat=i_stat)
         call memocc(i_stat,i_all,'atoms','setup_chemicaldata')
         !i_all=-product(shape(atoms2))*kind(atoms2)
         !deallocate(atoms2,stat=i_stat)
         !call memocc(i_stat,i_all,'atoms2','setup_chemicaldata')
         !i_all=-product(shape(atoms2T))*kind(atoms2T)
         !deallocate(atoms2T,stat=i_stat)
         !call memocc(i_stat,i_all,'atoms2T','setup_chemicaldata')
         i_all=-product(shape(qch))*kind(qch)
         deallocate(qch,stat=i_stat)
         call memocc(i_stat,i_all,'qch','setup_chemicaldata')
         i_all=-product(shape(rn))*kind(rn)
         deallocate(rn,stat=i_stat)
         call memocc(i_stat,i_all,'rn','setup_chemicaldata')
         i_all=-product(shape(ns))*kind(ns)
         deallocate(ns,stat=i_stat)
         call memocc(i_stat,i_all,'ns','setup_chemicaldata')
         i_all=-product(shape(ne))*kind(ne)
         deallocate(ne,stat=i_stat)
         call memocc(i_stat,i_all,'ne','setup_chemicaldata')
      end if

   end subroutine setup_chemicaldata

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_globcoord
   !> @brief Sets up the global coordinates in array coord
   !> @details Translates the basis atoms so that they are always inside the
   !> first supercell.
   !----------------------------------------------------------------------------
   subroutine setup_globcoord(Natom,NA,Bas,C1,C2,C3,coord,N1,N2,N3,atype,anumb,     &
      do_prnstruct,do_prn_poscar,simid,do_ralloy,Natom_full,acellnumb,achtype,      &
      atype_ch,asite_ch,achem_ch,block_size)
      !
      implicit none
      !
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      integer, intent(in) :: do_prnstruct    !< Print Hamiltonian information (0/1)
      integer, intent(in) :: do_prn_poscar   !< Print geometry on POSCAR format (0/1)
      character(len=8), intent(in) :: simid  !< Name of simulation
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      integer, dimension(Natom_full), intent(in) :: acellnumb  !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full), intent(in) :: achtype    !< Chemical type of atoms (full list)
      integer, dimension(Natom_full), intent(in) :: atype_ch   !< Actual type of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: asite_ch   !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: achem_ch   !< Chemical type of atoms (reduced list)
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      ! .. In/Out variables
      integer, intent(inout) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3,NA) , intent(inout) :: Bas !< Coordinates for basis atoms
      real(dblprec), dimension(:,:), allocatable, intent(inout) :: coord !< Coordinates of atoms
      ! .. Local variables
      integer :: i,i_stat
      integer :: i0, i1, i2, i3, ii1, ii2, ii3
      integer :: iatom
      character(len=20) :: filn, filn2
      real(dblprec) :: detmatrix
      real(dblprec), dimension(3) :: icvec,bsf
      real(dblprec), dimension(3,3) :: invmatrix
      !
      
      ! Fold all atoms into first unit cell:
      ! Calculate inverse of basis matrix
      detmatrix=C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+&
         C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+&
         C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1)
      invmatrix=0.0_dblprec
      if (abs(detmatrix)>dbl_tolerance) then
         invmatrix(1,1)=(C2(2)*C3(3)-C3(2)*C2(3))/detmatrix
         invmatrix(1,2)=(C1(3)*C3(2)-C3(3)*C1(2))/detmatrix
         invmatrix(1,3)=(C1(2)*C2(3)-C2(2)*C1(3))/detmatrix
         invmatrix(2,1)=(C2(3)*C3(1)-C3(3)*C2(1))/detmatrix
         invmatrix(2,2)=(C1(1)*C3(3)-C3(1)*C1(3))/detmatrix
         invmatrix(2,3)=(C1(3)*C2(1)-C2(3)*C1(1))/detmatrix
         invmatrix(3,1)=(C2(1)*C3(2)-C3(1)*C2(2))/detmatrix
         invmatrix(3,2)=(C1(2)*C3(1)-C3(2)*C1(1))/detmatrix
         invmatrix(3,3)=(C1(1)*C2(2)-C2(1)*C1(2))/detmatrix
      end if

      do I0=1,NA
         !find coordinate vector in basis coordinates
         icvec(1)=Bas(1,I0)*invmatrix(1,1)+Bas(2,I0)*invmatrix(2,1)+Bas(3,I0)*invmatrix(3,1)
         icvec(2)=Bas(1,I0)*invmatrix(1,2)+Bas(2,I0)*invmatrix(2,2)+Bas(3,I0)*invmatrix(3,2)
         icvec(3)=Bas(1,I0)*invmatrix(1,3)+Bas(2,I0)*invmatrix(2,3)+Bas(3,I0)*invmatrix(3,3)
         ! fold back to original cell
         bsf(1)=floor(icvec(1)+1d-7)
         bsf(2)=floor(icvec(2)+1d-7)
         bsf(3)=floor(icvec(3)+1d-7)
         !
         Bas(1,I0)=Bas(1,I0)-bsf(1)*C1(1)-bsf(2)*C2(1)-bsf(3)*C3(1)
         Bas(2,I0)=Bas(2,I0)-bsf(1)*C1(2)-bsf(2)*C2(2)-bsf(3)*C3(2)
         Bas(3,I0)=Bas(3,I0)-bsf(1)*C1(3)-bsf(2)*C2(3)-bsf(3)*C3(3)
      end do

      ! Open file for writing the coordinates
      if(do_prnstruct==1.or.do_prnstruct==2.or.do_prnstruct==4) then
         write (filn,'(''coord.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn)
      end if
      if(do_prn_poscar==1) then
         filn2='POSCAR'
         open(ofileno2, file=filn2)
         write(ofileno2,'(a)') 'Exported from UppASD'
         write(ofileno2,'(a)') '1.0'
         write(ofileno2,'(3f12.6)') N1*C1(1), N1*C1(2), N1*C1(3)
         write(ofileno2,'(3f12.6)') N2*C2(1), N2*C2(2), N2*C2(3)
         write(ofileno2,'(3f12.6)') N3*C3(1), N3*C3(2), N3*C3(3)
         !write(ofileno2,'(3f12.6)') C1(1), C1(2), C1(3)
         !write(ofileno2,'(3f12.6)') C2(1), C2(2), C2(3)
         !write(ofileno2,'(3f12.6)') C3(1), C3(2), C3(3)
         ! 'Mn' is a dummy label, not used
         write(ofileno2,'(a)') 'Mn'
         write(ofileno2,'(i8)') Natom
         !write(ofileno2,'(a)') 'Direct'
         write(ofileno2,'(a)') 'Cartesian'
      end if

      iatom=0
      ! Allocate coordinate array
      allocate(coord(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(coord))*kind(coord),'coord','setup_globcoord')
      i=0
      do II3=0, N3-1,block_size
         do II2=0, N2-1,block_size
            do II1=0, N1-1,block_size
               do I3=II3, min(II3+block_size-1,N3-1)
                  do I2=II2,min(II2+block_size-1,N2-1)
                     do I1=II1, min(II1+block_size-1,N1-1)
                        do I0=1, NA
                           i=i+1
                           if (do_ralloy==0) then
                              coord(1,i)=I1*C1(1)+I2*C2(1)+I3*C3(1)+Bas(1,I0)
                              coord(2,i)=I1*C1(2)+I2*C2(2)+I3*C3(2)+Bas(2,I0)
                              coord(3,i)=I1*C1(3)+I2*C2(3)+I3*C3(3)+Bas(3,I0)
                              if(do_prnstruct==1.or.do_prnstruct==2.or.do_prnstruct==4) then
                                 write(ofileno,'(i7,3f12.6,2i6)') i,coord(1:3,i),atype(i),anumb(i)
                              end if
                              if(do_prn_poscar==1) then
                                 write(ofileno2,'(3f12.6)') coord(1:3,i)
                              end if
                           else
                              if (achtype(i) /= 0) then
                                 iatom=acellnumb(i)
                                 coord(1,iatom)=I1*C1(1)+I2*C2(1)+I3*C3(1)+Bas(1,I0)
                                 coord(2,iatom)=I1*C1(2)+I2*C2(2)+I3*C3(2)+Bas(2,I0)
                                 coord(3,iatom)=I1*C1(3)+I2*C2(3)+I3*C3(3)+Bas(3,I0)
                                 if(do_prnstruct==1.or.do_prnstruct==2) then
                                    write(ofileno,'(i7,3f12.6,4i7)') iatom,coord(1:3,iatom),&
                                       atype_ch(iatom),asite_ch(iatom),achem_ch(iatom)
                                 end if
                                 if(do_prn_poscar==1) then
                                    write(ofileno2,'(a)') 'Writing to POSCAR for random alloy &
                                         currently not supported.'
                                 end if
                              end if
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      if(do_prnstruct==1.or.do_prnstruct==4) then
         close(ofileno)
      end if
      if(do_prn_poscar==1) close(ofileno2)
    end subroutine setup_globcoord


    subroutine rescale_lattvec(scalefac,C1,C2,C3)
      !
      implicit none
      !
      real(dblprec), intent(in) :: scalefac !< Lattice vector rescaling factor
      real(dblprec), dimension(3) :: C1 !< First lattice vector
      real(dblprec), dimension(3) :: C2 !< Second lattice vector
      real(dblprec), dimension(3) :: C3 !< Third lattice vector

      ! Rescaling of lattice vectors. Intended mainly for debugging
      C1(1:3)=scalefac*C1(1:3)
      C2(1:3)=scalefac*C2(1:3)
      C3(1:3)=scalefac*C3(1:3)
      
    end subroutine rescale_lattvec


   
end module Geometry
