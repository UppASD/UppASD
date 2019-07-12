!-------------------------------------------------------------------------------
! MODULE: ElkGeometry
!> @brief
!> Routines for reading and saving geometry and moment configuration on Elk format
!> @author
!> J. Hellsvik
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module ElkGeometry
   use Parameters
   use Profiling
   use InputData
   use AMS
   use Constants
   use math_functions
   
   implicit none
   public


contains


   !> Prints geometry and spin configuration on ELK format
   !> Uses Cartesian coordinates for the internal coordinates (ELK 'molecule .true.')
   !> @todo Support chemical species. For now 'Mn' is used globally
   subroutine prn_elk_geometry(Natom, Mensemble, simid, mstep, emom, mmom, &
        NA, N1, N2, N3, C1, C2, C3, atype_inp, coord, jfile, maptype, posfiletype, Bas)
     
      ! Why the spurious circular dependencies for the use statements?
      !use InputHandler, only : read_exchange_getMaxNoShells!, read_exchange_getNeighVec
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: mstep !< Current simulation step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector      
      integer, dimension(NA), intent(inout)           :: atype_inp   !< Type of atom from input
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      integer :: maptype !< Format for input data (1=direct format,2=bgfm style)
      character :: posfiletype !< posfile type (D)irect or (C)arteisian coordinates in posfile
      !real(dblprec), dimension(:,:) :: Bas        !< Coordinates for basis atoms
      real(dblprec), dimension(:,:), allocatable :: Bas        !< Coordinates for basis atoms
   
      integer :: i, j, j0
      integer :: flines, no_shells
      integer :: isite, jsite
      character(len=30) :: filn, filn2, jfile
      real(dblprec) :: jcoup

      real(dblprec), dimension(3) :: r_tmp, r_red 
      real(dblprec), dimension(3,Natom) :: mom_rounded
      real(dblprec), dimension(3) :: r_tmp2

      integer :: ichem, iconf
      
      !.. Executable statements

      i=1
      j0=1
      ! Before export, round the Cartesian components of the magnetic moments
      ! to three significant digits
      do j=1, Natom
         mom_rounded(1:3,j) = emom(1:3,j,i)
         mom_rounded(1,j) = nint(1000 * mom_rounded(1,j)) * 0.001_dblprec
         mom_rounded(2,j) = nint(1000 * mom_rounded(2,j)) * 0.001_dblprec
         mom_rounded(3,j) = nint(1000 * mom_rounded(3,j)) * 0.001_dblprec
         mom_rounded(1:3,j) = mom_rounded(1:3,j) / norm2 ( mom_rounded(1:3,j) )
         mom_rounded(1:3,j) = mom_rounded(1:3,j) * mmom(j,i)
      end do

      filn='GEOMETRY.OUT.uppasd'
      open(ofileno, file=filn)

      !write(ofileno,'(a)') 'Exported from UppASD'
      write(ofileno,'(a)') '' 
      write(ofileno,'(a)') 'scale'
      write(ofileno,'(a)') '    1.000'
      write(ofileno,'(a)') '' 
      write(ofileno,'(a)') 'scale1'
      write(ofileno,'(a)') '    1.000'
      write(ofileno,'(a)') '' 
      write(ofileno,'(a)') 'scale2'
      write(ofileno,'(a)') '    1.000'
      write(ofileno,'(a)') '' 
      write(ofileno,'(a)') 'scale3'
      write(ofileno,'(a)') '    1.000'
      write(ofileno,'(a)') '' 
      write(ofileno,'(a)') 'avec'
      write(ofileno,'(4x,3f12.6)') N1*C1(1), N1*C1(2), N1*C1(3)
      write(ofileno,'(4x,3f12.6)') N2*C2(1), N2*C2(2), N2*C2(3)
      write(ofileno,'(4x,3f12.6)') N3*C3(1), N3*C3(2), N3*C3(3)
      write(ofileno,'(a)') ''
      write(ofileno,'(a)') 'molecule'
      write(ofileno,'(a)') '    T'
      write(ofileno,'(a)') ''
      write(ofileno,'(a)') 'atoms'
      write(ofileno,'(4x,i8)') 1
      write(ofileno,'(4x,a)') "'Mn.in'"
      write(ofileno,'(4x,i8)') Natom

      i=1
      do j=1, Natom
         write (ofileno,10004) coord(1:3,j), mom_rounded(1:3,j)
      end do
      
      close(ofileno)

      iconf = 1
      ichem = 1

      ! Write Cartesian positions to posfile_cart
      open(ofileno, file='posfile_cart',status='replace')
      do isite=1,NA
         write(ofileno,'(2i8,5f12.6)') isite, atype_inp(isite), bas(1:3,isite)
      end do
      close(ofileno)
      
      ! Write magnetic moments to momfile_cart
      open(ofileno, file='momfile_cart',status='replace')
      do isite=1,NA
         write(ofileno,'(2i8,5f12.6)') isite, ichem, ammom_inp(isite,ichem,iconf), aemom_inp(1:3,isite,ichem,iconf)
      end do
      close(ofileno)

      ! Get number of exchange interactions
      flines=0
      ! Hard-wire exchange interaction file name for now
      filn='JFILE'
      open(ifileno, file=filn)
      do
         read(ifileno,*,end=200) isite, jsite, r_tmp(1:3), jcoup
         flines=flines+1
      end do
200   continue
      rewind(ifileno)
      write(*,*) 'flines', flines
      
      ! Read the exchange interactions and write them to jfile_cart using
      ! maptype 1, posfiletype C form
      filn2='jfile_cart'
      open(ofileno, file=filn2)

      do j=1,flines
         read(ifileno,*) isite, jsite, r_tmp(1:3), jcoup

         ! Got spurious circular dependencies, copied
         ! lines from read_exchange_getNeighVec as a quick fix
         !call read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

         if(maptype==3) then
            ! Calculate proper neighbour vector (from "RSPt")
            r_tmp2=r_tmp
            r_red(1)=Bas0(1,jsite)-Bas0(1,isite)+C1(1)*r_tmp2(1)+C2(1)*r_tmp2(2)+C3(1)*r_tmp2(3)
            r_red(2)=Bas0(2,jsite)-Bas0(2,isite)+C1(2)*r_tmp2(1)+C2(2)*r_tmp2(2)+C3(2)*r_tmp2(3)
            r_red(3)=Bas0(3,jsite)-Bas0(3,isite)+C1(3)*r_tmp2(1)+C2(3)*r_tmp2(2)+C3(3)*r_tmp2(3)
         elseif(maptype==2) then
            ! Calculate proper neighbour vector (from "bgfm")
            r_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(1)+C2(1)*r_tmp(2)+C3(1)*r_tmp(3)
            r_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(1)+C2(2)*r_tmp(2)+C3(2)*r_tmp(3)
            r_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(1)+C2(3)*r_tmp(2)+C3(3)*r_tmp(3)
         else
            ! Calculates neighbour vectors from direct coordinates or Cartesian
            ! coordinates, corresponding to how the atomic positions are entered
            if (posfiletype=='C') then
               r_red=r_tmp
            elseif (posfiletype=='D') then
               r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
               r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
               r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
            else
               stop 'Only posfiletype = C or D is currently supported'
            endif
         end if

         write(ofileno,'(2i8,5f12.6)') isite, jsite, r_red(1:3), jcoup, norm2(r_red(1:3))
      end do

      close(ifileno)      
      close(ofileno)


10002 format (i8)
10003 format (i8,i8,2x,es16.8,es16.8,es16.8,es16.8)
10004 format (4x,6f12.6)

   end subroutine prn_elk_geometry


   !---------------------------------------------------------------------------
   !> @brief
   !> Read geometry and spin configuration on ELK format
   !> Uses Cartesian coordinates for the internal coordinates (ELK 'molecule .true.')
   !> @todo Support chemical species. For now 'Mn' is used globally
   !
   !> @author
   !> Johan Hellsvik
   !---------------------------------------------------------------------------
   subroutine read_elk_geometry()
      !
      !
      implicit none
      !
      real(dblprec), dimension(:), allocatable :: atype_inp2

      integer :: i, j
      integer :: i_stat, i_err
      integer :: isite, msite, ichem, iconf

      real(dblprec) :: scale, scale1, scale2, scale3
      real(dblprec),dimension(3) :: inpC1, inpC2, inpC3
      real(dblprec), dimension(3,3) :: A, invA, B, T
      
      character(len=30) :: filn

      real(dblprec), dimension(Natom) :: mom_col

      integer :: itype, mtype, iat, msite2, i_all, ifileno2
      real(dblprec), dimension(3) :: tmp, emomref
      
      real(dblprec):: collinear_tolerance, emomcomp, emomcomp2, emomcomp3
      logical :: collinear, ferromagnetic

      !.. Executable statements      
      write(*,'(2x,a)') 'Geometry and initial spin configuration read on ELK format from GEOMETRY.OUT'
      
      iconf = 1
      ichem = 1
      ifileno2 = 99
      collinear_tolerance = 1.0d-6

      filn='GEOMETRY.OUT'
      open(ifileno, file=filn)

      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*) scale
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*) scale1
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*) scale2
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*) scale3
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*) inpC1(1), inpC1(2), inpC1(3)
      read(ifileno,*) inpC2(1), inpC2(2), inpC2(3)
      read(ifileno,*) inpC3(1), inpC3(2), inpC3(3)
      inpC1(1:3) = scale*scale1*inpC1(1:3)
      inpC2(1:3) = scale*scale2*inpC2(1:3)
      inpC3(1:3) = scale*scale3*inpC3(1:3)
      !C1(1:3) = scale*scale1*inpC1(1:3)
      !C2(1:3) = scale*scale2*inpC2(1:3)
      !C3(1:3) = scale*scale3*inpC3(1:3)
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*)
      read(ifileno,*) msite

      ! Set NA to number of atoms in the magnetic cell
      na=msite
      !nt=na
      
      ! Allocate position input arrays
      allocate(bas(3,na),stat=i_stat)
      call memocc(i_stat,product(shape(bas))*kind(bas),'bas','read_elk_geometry')
      allocate(atype_inp(na),stat=i_stat)
      call memocc(i_stat,product(shape(atype_inp))*kind(atype_inp),'atype_inp','read_elk_geometry')
      allocate(anumb_inp(na),stat=i_stat)
      call memocc(i_stat,product(shape(anumb_inp))*kind(anumb_inp),'anumb_inp','read_elk_geometry')

      ! Allocate moment input arrays
      allocate(ammom_inp(na,1,1),stat=i_stat)
      call memocc(i_stat,product(shape(ammom_inp))*kind(ammom_inp),'ammom_inp','read_elk_geometry')
      ammom_inp=0.0_dblprec
      allocate(aemom_inp(3,na,1,1),stat=i_stat)
      call memocc(i_stat,product(shape(aemom_inp))*kind(aemom_inp),'aemom_inp','read_elk_geometry')
      aemom_inp=0.0_dblprec
      allocate(Landeg_ch(na,1,1),stat=i_stat)
      call memocc(i_stat,product(shape(Landeg_ch))*kind(Landeg_ch),'Landeg_ch','read_elk_geometry')
      Landeg_ch=2.0_dblprec
      !Landeg_ch=Landeg_global

      ! Open position file
      open(ifileno2, file=posfile)
      mtype=0
      msite2=0
      ! Pre-read file to get max no. sites, types and chemical species
      do
         read(ifileno2,*,end=200) isite, itype
         msite2=max(msite2,isite)
         mtype=max(mtype,itype)
      end do
      200 continue
      rewind(ifileno2)
      !na=msite

      ! Set NT to number of atoms in the chemical cell
      nt=mtype

      allocate(atype_inp2(na),stat=i_stat)
      call memocc(i_stat,product(shape(atype_inp2))*kind(atype_inp2),'atype_inp2','read_elk_geometry')

      ! Read basis atoms and setup type array
      ! Site, Type, Rx, Ry, Rz
      do iat=1, msite2
         read (ifileno2,*) isite, itype,  tmp
         atype_inp2(isite)=itype
      enddo
      
      close (ifileno2)

      do isite=1,msite
         read (ifileno,*) bas(1:3,isite), aemom_inp(1:3,isite,ichem,iconf)
         ammom_inp(isite,ichem,iconf)=sqrt(aemom_inp(1,isite,ichem,iconf)**2+aemom_inp(2,isite,ichem,iconf)**2+aemom_inp(3,isite,ichem,iconf)**2)
         aemom_inp(1:3,isite,ichem,iconf)=aemom_inp(1:3,isite,ichem,iconf)/ammom_inp(isite,ichem,iconf)
         ! Use the same types as for the atoms in the chemical cell
         atype_inp(isite)=atype_inp2(mod(isite-1,msite2)+1)
         !atype_inp(isite)=isite
         ! Redundant but kept for the time beeing
         anumb_inp(isite)=isite
      end do

      close (ifileno)

      ! Calculate the linear transformation (T) that takes the original lattice vectors (A)
      ! to the new lattice vectors (B)

      ! The original lattice vectors
      A(1,1:3) = C1(1:3)
      A(2,1:3) = C2(1:3)
      A(3,1:3) = C3(1:3)
      ! The new lattice vectors
      B(1,1:3) = inpC1(1:3)
      B(2,1:3) = inpC2(1:3)
      B(3,1:3) = inpC3(1:3)
      ! From B = T*A the linear transformation T is obtained as T=B*invA
      call matinvcramer3(A, invA)
      call matmatmul(3, B, invA, T)
      write(*,'(2x,a)') 'Calculates the linear transformation (T) that takes the original'
      write(*,'(2x,a)') 'lattice vectors (A) to the new lattice vectors (B).'
      write(*,'(2x,a)') 'The elements of the transformation matrix T'
      write(*,'(2x,3f12.6)') T(1,1:3)
      write(*,'(2x,3f12.6)') T(2,1:3)
      write(*,'(2x,3f12.6)') T(3,1:3)
      write(*,'(2x,a)') 'The original lattice vectors'
      write(*,'(2x,3f12.6)') A(1,1:3)
      write(*,'(2x,3f12.6)') A(2,1:3)
      write(*,'(2x,3f12.6)') A(3,1:3)
      write(*,'(2x,a)') 'The new lattice vectors'
      write(*,'(2x,3f12.6)') B(1,1:3)
      write(*,'(2x,3f12.6)') B(2,1:3)
      write(*,'(2x,3f12.6)') B(3,1:3)
      
      ! Replace the original lattice vectors with the new lattice vectors
      C1(1:3)=inpC1(1:3)
      C2(1:3)=inpC2(1:3)
      C3(1:3)=inpC3(1:3)
      
      ! Write Cartesian positions to posfile_cart
      open(ofileno, file='posfile_cart',status='replace')
      do isite=1,msite
         write(ofileno,'(2i8,5f12.6)') isite, atype_inp(isite), bas(1:3,isite)
      end do
      close(ofileno)
      
      ! Write magnetic moments to momfile_cart
      open(ofileno, file='momfile_cart',status='replace')
      do isite=1,msite
         write(ofileno,'(2i8,5f12.6)') isite, ichem, ammom_inp(isite,ichem,iconf), aemom_inp(1:3,isite,ichem,iconf)
      end do
      close(ofileno)
      
      ! Check if the initial spin configuration is collinear
      ! and if it is ferromagnetic
      collinear = .true.
      ferromagnetic = .true.
      emomref(1:3) = aemom_inp(1:3,1,ichem,iconf)
      do isite=1,msite
         emomcomp = emomref(1) * aemom_inp(1,isite,ichem,iconf) &
              + emomref(2) * aemom_inp(2,isite,ichem,iconf) &
              + emomref(3) * aemom_inp(3,isite,ichem,iconf)
         emomcomp2 = 1.0_dblprec - abs(emomcomp)
         emomcomp3 = 1.0_dblprec - emomcomp
         !write(*,'(a,f12.6,a,f12.6,a,f12.6)') 'emomcomp', emomcomp, ' emomcomp2', emomcomp2, ' emomcomp3', emomcomp3
         if( emomcomp2 .gt. collinear_tolerance ) collinear = .false.
         if( emomcomp3 .gt. collinear_tolerance ) ferromagnetic = .false.
      end do
      if( collinear ) then
         if (ferromagnetic ) then
            write(*,'(2x,a)') 'The initial spin configuration is collinear and ferromagnetic.'
         else
            write(*,'(2x,a)') 'The initial spin configuration is collinear antiferromagnetic'
            write(*,'(2x,a)') 'or collinear ferrimagnetic.'
         end if
      else
         write(*,*) 'The initial spin configuration is noncollinear.'
         write(*,*) 'Adiabatic magnon spectra calculation will not be performed.'
         do_ams = 'N'
         do_magdos = 'N'
      end if
      
      i_all=-product(shape(atype_inp2))*kind(atype_inp2)
      deallocate(atype_inp2,stat=i_stat)
      call memocc(i_stat,i_all,'atype_inp2','allocate_initmag')

   end subroutine read_elk_geometry


   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_elkgeometry(ifile)

     use FileParser

     implicit none

     ! ... Formal Arguments ...
     integer, intent(in) :: ifile   !< File to read from
     !
     ! ... Local Variables ...
     character(len=50) :: keyword, cache
     integer :: rd_len, i_err, i_errb
     logical :: comment

     do
10      continue
        ! Read file character for character until first whitespace
        keyword=""
        call bytereader(keyword,rd_len,ifile,i_errb)

        ! converting Capital letters
        call caps2small(keyword)

        ! check for comment markers (currently % and #)
        comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
             (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
             (scan(trim(keyword),'!')==1))

        if (comment) then
           read(ifile,*)
        else
           ! Parse keyword
           keyword=trim(keyword)
           select case(keyword)
           ! This is the flags for the ElkGeometry module

           case('do_prn_elk')
              read(ifile,*,iostat=i_err) do_prn_elk
              if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

           case('do_read_elk')
              read(ifile,*,iostat=i_err) do_read_elk
              if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

           end select
        end if

        ! End of file
        if (i_errb==20) goto 20
        ! End of row
        if (i_errb==10) goto 10
     end do

20   continue

     rewind(ifile)
     return
   end subroutine read_parameters_elkgeometry

   
end module ElkGeometry
