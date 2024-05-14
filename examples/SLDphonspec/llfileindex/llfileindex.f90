!--------------------------------------------------!---------------------------
!    Program for "llfile"  which is lattice for dynamics from phonopy-VASP
!    By Banasree Sadhukhan, Date : 11/03/2022
!    Revised by Johan Hellsvik, Jan 2023
!--------------------------------------------------!---------------------------

program main
  implicit none
  ! Set array size large enough - replace eventaually with prereading of sizes and dynamic allocation
  integer, parameter :: dim=200000 ! Allocate dimension
  integer, parameter :: dim2=2000 ! Allocate dimension
  real(8) :: phi(9,dim2,dim2)
  integer :: natom_uc, natom_sc! Atoms in unitcell and supercell
  integer :: nx, ny, nz, na
  integer :: ia, iatom
  integer :: jx, jy, jz, ja, jatom
  integer :: x0, y0, z0
  integer :: lx, ly, lz

  ! Read input file specifying the size of the supercell
  ! and the chosen central i site.
  open(unit=1,file="llfileindex.dat")
  read(1,*) nx, ny, nz, na
  read(1,*) x0, y0, z0
  close(1)

  ! Read the force constants
  open(unit=1,file="FORCE_CONSTANTS")
  read(1,*)
  natom_uc = na
  natom_sc = nx*ny*nz*na
  do iatom=1,natom_sc
     do jatom=1,natom_sc
        !write(*,*) iatom, jatom
        read(1,*)
        read(1,*) phi(1:3,iatom,jatom)
        read(1,*) phi(4:6,iatom,jatom)
        read(1,*) phi(7:9,iatom,jatom)
     end do
  end do
  close(1)
  !----------------------------------------------------------------------------------- 
  
  open(unit=1,file="llfile")
  open(unit=2,file="checkllfile")
  ! Loop over the site indices in the Phonopy supercell
  do ia=1,na
     do jz=1,nz
        do jy=1,ny
           do jx=1,nx
              do ja=1,na
                 jatom = (ja-1)*nx*ny*nz + (jz-1)*nx*ny + (jy-1)*nx + jx
                 iatom = (ia-1)*nx*ny*nz + (z0-1)*nx*ny + (y0-1)*nx + x0
                 !write(*,*) iatom, jatom
                 lx = jx - x0
                 ly = jy - y0
                 lz = jz - z0
                 write(1,"(2I4,3I4,9F15.7)") ia, ja, lx, ly, lz, phi(1:9,iatom,jatom)
                 write(2,"(2I4,5I4,9F15.7)") ia, ja, lx, ly, lz, iatom, jatom, phi(1:9,iatom,jatom)
              end do
           end do
        end do
     end do
  end do
  close(1)
  close(2)
  
end program main


