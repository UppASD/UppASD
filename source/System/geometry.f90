!> Routines for setting up structural and chemical information about the system
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module Geometry
  use Parameters
  use Profiling

  implicit none


  private
  public :: setup_geometry

contains


  !> Wrapper routine for setting up structural and chemical information about the system
  subroutine setup_geometry(Natom, NA, N1, N2, N3, Bas, C1, C2, C3, coord, atype, anumb, &
       atype_inp, anumb_inp, do_prnstruct, tseed, simid, do_ralloy, Natom_full, Nchmax, &
       Nch, acellnumb, acellnumbrev, achtype, chconc, atype_ch, asite_ch, achem_ch)

    implicit none

    integer, intent(inout) :: Natom !< Number of atoms in system
    integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    real(dblprec), dimension(3,NA) , intent(inout) :: Bas !< Coordinates for basis atoms
    real(dblprec), dimension(3) , intent(in) :: C1 !< First lattice vector
    real(dblprec), dimension(3) , intent(in) :: C2 !< Second lattice vector
    real(dblprec), dimension(3) , intent(in) :: C3 !< Third lattice vector
    real(dblprec), dimension(:,:),allocatable, intent(inout) :: coord !< Coordinates of atoms
    integer, dimension(Natom_full), intent(inout) :: atype !< Type of atom
    integer, dimension(Natom_full), intent(inout) :: anumb !< Atom number in cell
    integer, dimension(NA), intent(inout) :: atype_inp  !< Type of atom from input
    integer, dimension(NA), intent(inout) :: anumb_inp !< Atom number in cell from input
    integer, intent(in) :: do_prnstruct !< Print Hamiltonian information (0/1)
    integer, intent(in) :: tseed                                !< Temperature seed
    character(len=8) :: simid                                   !< Name of simulation
    integer, intent(in) :: do_ralloy                            !< Random alloy simulation (0/1)
    integer, intent(in) :: Nchmax                               !< Max number of chemical components on each site in cell
    integer, dimension(:), intent(in) :: Nch                    !< Number of chemical components on each site in cell
    integer, dimension(:), intent(out) :: acellnumb             !< List for translating atom no. in full cell to actual cell
    integer, dimension(:), intent(out) :: acellnumbrev          !< List for translating atom no. in actual cell to full cell
    integer, dimension(:), intent(out) :: achtype               !< Chemical type of atoms (full list)
    real(dblprec), dimension(:,:), intent(in) :: chconc         !< Chemical concentration on sites
    integer, dimension(:), intent(out) :: atype_ch              !< Actual type of atom for dilute system
    integer, dimension(:), intent(out) :: asite_ch              !< Actual site of atom for dilute system
    integer, dimension(:), intent(out) :: achem_ch              !< Chemical type of atoms (reduced list)

    write (*,'(2x,a)',advance='no') 'Set up types of atoms'
    call setup_type_and_numb(Natom_full, NA, N1, N2, N3, atype, anumb, atype_inp, anumb_inp)
    write (*,'(a)') ' done'

    write (*,'(2x,a)',advance='no') 'Set up chemical information of alloy'
    if(do_ralloy==1) then
       call setup_chemicaldata(Natom, NA, N1, N2, N3, atype, tseed, &
            do_ralloy, Natom_full, Nchmax, Nch, achtype, acellnumb, acellnumbrev, &
            chconc, atype_ch, asite_ch, achem_ch)
    else
       Natom=Natom_full
    end if
    write (*,'(a)') ' done'

    write (*,'(2x,a)',advance='no') 'Set up global coordinates'
    call setup_globcoord(Natom, NA, Bas, C1, C2, C3, coord, N1, N2, N3, atype, anumb, do_prnstruct, simid, &
         do_ralloy, Natom_full, acellnumb, achtype, atype_ch, asite_ch, achem_ch)
    write (*,'(a)') ' done'


  end subroutine setup_geometry


  !> Sets up the type of atoms in the system
  subroutine setup_type_and_numb(Natom_full, NA, N1, N2, N3, atype, anumb, atype_inp, anumb_inp)
    !
    implicit none

    integer, intent(in) :: Natom_full !< Number of atoms in system
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    integer, dimension(Natom_full), intent(inout) :: atype !< Type of atom
    integer, dimension(Natom_full), intent(inout) :: anumb !< Atom number in cell
    integer, dimension(NA), intent(inout) :: atype_inp  !< Type of atom from input
    integer, dimension(NA), intent(inout) :: anumb_inp !< Atom number in cell from input
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


  !> Sets up chemical information of the system in case of random alloy
  !! IMPORTANT: In case of dilute system, Natom variable is reduced.
  subroutine setup_chemicaldata(Natom, NA, N1, N2, N3, atype, tseed, &
       do_ralloy, Natom_full, Nchmax, Nch, achtype, acellnumb, acellnumbrev, &
       chconc, atype_ch, asite_ch, achem_ch)
    !
    use RandomNumbers, only : rng_uniform
    use Sorting, only : mergesort
    !
    implicit none

    integer, intent(out) :: Natom !< Number of atoms in system
    integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    integer, dimension(Natom_full), intent(inout) :: atype !< Type of atom
    integer, intent(in) :: tseed  !< Temperature seed
    integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
    integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
    integer, dimension(NA), intent(in) :: Nch !< Number of chemical components on each site in cell
    integer, dimension(Natom_full), intent(out) :: achtype !< Chemical type of atoms (full list)
    integer, dimension(Natom_full), intent(out) :: acellnumb !< List for translating atom no. in full cell to actual cell
    integer, dimension(Natom_full), intent(out) :: acellnumbrev !< List for translating atom no. in actual cell to full cell
    real(dblprec), dimension(NA,Nchmax), intent(in) :: chconc !< Chemical concentration on sites
    integer, dimension(Natom_full), intent(out) :: atype_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), intent(out) :: asite_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), intent(out) :: achem_ch !< Chemical type of atoms (reduced list)

    integer :: A1
    integer :: i0, i1, i2, i3
    integer :: Ncell, iatom
    integer :: ia, ich, i
    integer :: i_stat, i_all
    integer, dimension(1) :: irn
    integer, dimension(:,:), allocatable :: atoms
    integer, dimension(:,:), allocatable :: qch
    integer, dimension(:), allocatable :: ns, ne
    integer, dimension(:), allocatable :: atoms2, atoms2T
    real(dblprec), dimension(:), allocatable :: rn

    if (do_ralloy==1) then

       ! Same seed as for Monte-Carlo
       Ncell = N1*N2*N3
       allocate(atoms(Ncell,Na),stat=i_stat)
       call memocc(i_stat,product(shape(atoms))*kind(atoms),'atoms','setup_type_and_numb')
       allocate(atoms2(Ncell),stat=i_stat)
       call memocc(i_stat,product(shape(atoms2))*kind(atoms2),'atoms2','setup_type_and_numb')
       allocate(atoms2T(Ncell),stat=i_stat)
       call memocc(i_stat,product(shape(atoms2T))*kind(atoms2T),'atoms2T','setup_type_and_numb')
       allocate(qch(Nchmax,Na),stat=i_stat)
       call memocc(i_stat,product(shape(qch))*kind(qch),'qch','setup_type_and_numb')
       allocate(rn(Ncell),stat=i_stat)
       call memocc(i_stat,product(shape(rn))*kind(rn),'rn','setup_type_and_numb')
       allocate(ns(Nchmax),stat=i_stat)
       call memocc(i_stat,product(shape(ns))*kind(ns),'ns','setup_type_and_numb')
       allocate(ne(Nchmax),stat=i_stat)
       call memocc(i_stat,product(shape(ne))*kind(ne),'ne','setup_type_and_numb')
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
             rn(irn)=0.0d0
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
                      achem_ch(iatom)=achtype(i)
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
       call memocc(i_stat,i_all,'atoms','setup_type_and_numb')
       i_all=-product(shape(atoms2))*kind(atoms2)
       deallocate(atoms2,stat=i_stat)
       call memocc(i_stat,i_all,'atoms2','setup_type_and_numb')
       i_all=-product(shape(atoms2T))*kind(atoms2T)
       deallocate(atoms2T,stat=i_stat)
       call memocc(i_stat,i_all,'atoms2T','setup_type_and_numb')
       i_all=-product(shape(qch))*kind(qch)
       deallocate(qch,stat=i_stat)
       call memocc(i_stat,i_all,'qch','setup_type_and_numb')
       i_all=-product(shape(rn))*kind(rn)
       deallocate(rn,stat=i_stat)
       call memocc(i_stat,i_all,'rn','setup_type_and_numb')
       i_all=-product(shape(ns))*kind(ns)
       deallocate(ns,stat=i_stat)
       call memocc(i_stat,i_all,'ns','setup_type_and_numb')
       i_all=-product(shape(ne))*kind(ne)
       deallocate(ne,stat=i_stat)
       call memocc(i_stat,i_all,'ne','setup_type_and_numb')
    end if

  end subroutine setup_chemicaldata


  !> Sets up the global coordinates in array coord
  !> Translates the basis atoms so that they are always inside the
  !! first supercell.
  subroutine setup_globcoord(Natom, NA, Bas, C1, C2, C3, coord, N1, N2, N3, atype, anumb, do_prnstruct, simid, &
             do_ralloy, Natom_full, acellnumb, achtype, atype_ch, asite_ch, achem_ch)
    !
    implicit none
    !
    integer, intent(inout) :: Natom !< Number of atoms in system
    integer, intent(in) :: NA  !< Number of atoms in one cell
    real(dblprec), dimension(3,NA) , intent(inout) :: Bas !< Coordinates for basis atoms
    real(dblprec), dimension(3) , intent(in) :: C1 !< First lattice vector
    real(dblprec), dimension(3) , intent(in) :: C2 !< Second lattice vector
    real(dblprec), dimension(3) , intent(in) :: C3 !< Third lattice vector
    real(dblprec), dimension(:,:),allocatable, intent(inout) :: coord !< Coordinates of atoms
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
    integer, intent(in) :: do_prnstruct !< Print Hamiltonian information (0/1)
    character(len=8) :: simid !< Name of simulation
    integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
    integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
    integer, dimension(Natom_full), intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
    integer, dimension(Natom_full), intent(in) :: achtype !< Chemical type of atoms (full list)
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual type of atom for dilute system
    integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)

    integer :: i,i_stat
    integer :: i0, i1, i2, i3, ii1, ii2, ii3
    integer :: iatom
    character(len=20) :: filn
    real(dblprec) :: detmatrix
    real(dblprec), dimension(3) :: icvec,bsf
    real(dblprec), dimension(3,3) :: invmatrix
    !
    ! Fold all atoms into first unit cell:
    ! Calculate inverse of basis matrix
    detmatrix=C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+&
         C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+&
         C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1)
    invmatrix=0.0d0
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
       write (filn,'(''coord.'',a8,''.out'')') simid
       open(ofileno, file=filn)
    end if

    iatom=0
    ! Allocate coordinate array
    allocate(coord(3,Natom),stat=i_stat)
    call memocc(i_stat,product(shape(coord))*kind(coord),'coord','setup_globcoord')
    i=0
    print *,'starting geo',i
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
                               write(ofileno,'(i7,3f12.6,2i6)') i,coord(1:3,i),atype(i), anumb(i)
                            end if
                         else
                            if (achtype(i) /= 0) then
                               iatom=acellnumb(i)
                               coord(1,iatom)=I1*C1(1)+I2*C2(1)+I3*C3(1)+Bas(1,I0)
                               coord(2,iatom)=I1*C1(2)+I2*C2(2)+I3*C3(2)+Bas(2,I0)
                               coord(3,iatom)=I1*C1(3)+I2*C2(3)+I3*C3(3)+Bas(3,I0)
                               if(do_prnstruct==1.or.do_prnstruct==2) then
                                  write(ofileno,'(i7,3f12.6,4i7)') iatom,coord(1:3,iatom),&
                                     atype_ch(iatom), asite_ch(iatom), achem_ch(iatom)
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
    print *,'ending geoblocking',i

    if(do_prnstruct==1.or.do_prnstruct==4) then
       close(ofileno)
    end if
 end subroutine setup_globcoord


end module Geometry
