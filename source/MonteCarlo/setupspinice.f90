!> Spin-ice simulation routines
module SetupSpinIce

   use Parameters
   use Profiling
   use Constants
   use SpinIce
   use SpinIceData
   use prn_SpinIce

   implicit none

   integer, dimension(:,:), allocatable :: nm_ver_ver
   integer, dimension(:,:), allocatable :: nm_ver_atom
   integer, dimension(:), allocatable :: nmdim_ver_ver ! Dimension of neighbour map
   integer, dimension(:), allocatable :: nmdim_ver_atom

   private

   public :: setup_ice_neighbours,setup_vertex_geometry
   public :: read_parameters_spinice, spin_ice_init, read_vertices

contains

   !> setup system
   subroutine setup_ice_neighbours(Natom, Mensemble,NT,NA, N1, N2, N3, BC1, BC2, BC3, atype, Bas, sym, simid,coord)

      !.. Implicit declarations
      implicit none

      integer, intent(inout) :: Natom ! Number of atoms in system
      integer, intent(in) :: NT ! Number of types of atoms
      integer, intent(in) :: NA  ! Number of atoms in one cell
      integer, intent(in) :: N1  ! Number of cell repetitions in x direction
      integer, intent(in) :: N2  ! Number of cell repetitions in y direction
      integer, intent(in) :: N3  ! Number of cell repetitions in z direction
      integer, intent(in) :: Mensemble
      character(len=1), intent(in) :: BC1 ! Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 ! Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 ! Boundary conditions in z-direction
      integer, dimension(Natom), intent(inout) :: atype
      real(dblprec), dimension(3,NA), intent(inout) :: Bas !< Coordinates for basis atoms
      integer, intent(in) :: sym ! Symmetry of system (0-3)
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      character(len=8) :: simid ! Name of simulation


      character(len=1) :: neigh_type

      ! This initializes the neighbour lists needed for the Spin Ice system, vertex vertex and certex atom
      call setup_ver_nm(Nvertex, NT_ver, NA_ver, N1, N2, N3, C1_ver, C2_ver, C3_ver, BC1, BC2, BC3, atype_ver, Bas_ver, &
         max_no_equiv_ver, ver_sym, ver_dist_coord, nm_ver_ver, nmdim_ver_ver)
      call setup_ver_atom_nm(Nvertex, NT_ver, NA_ver, NA, N1, N2, N3, C1_ver, C2_ver, C3_ver, BC1, BC2, BC3, atype_ver, Bas_ver, Bas, &
         max_no_equiv_ver_atom, ver_sym, ver_atom_dist_coord, nm_ver_atom, nmdim_ver_atom)

      ! This allocates the arrays that are needed for the neighbour lists
      call allocate_spinice_data(Natom,Mensemble,1)

      ! This need to be called a number of times to make sure that both ver_ver and ver_atom lists are created
      call setup_neighbour_list(Nvertex, max_no_equiv_ver, nlistsize_ver, nlist_ver, nm_ver_ver, nmdim_ver_ver)
      call setup_neighbour_list(Nvertex, max_no_equiv_ver_atom, nlistsize_ver_atom, nlist_ver_atom, nm_ver_atom, nmdim_ver_atom)

      ! This needs to be called two times to make sure that ver-ver and ver-atom are printed in different files
      neigh_type='V'
      call prn_ver_ver(simid, Natom, Nvertex, NT, NT_ver, NA_ver, N1, N2, N3, atype, atype_ver, max_no_equiv_ver, ver_dist_coord,&
         nm_ver_ver, nmdim_ver_ver, neigh_type)
      call prn_ver_neigh(Nvertex, max_no_equiv_ver, nlistsize_ver, nlist_ver, simid, neigh_type)

      neigh_type='A'
      call prn_ver_ver(simid, Natom, Nvertex, NT, NT_ver, NA_ver, N1, N2, N3, atype, atype_ver,max_no_equiv_ver_atom, ver_atom_dist_coord,&
         nm_ver_atom, nmdim_ver_atom,neigh_type)
      call prn_ver_neigh(Nvertex, max_no_equiv_ver_atom, nlistsize_ver_atom, nlist_ver_atom, simid, neigh_type)

      call GEOMETRICAL_ICE(Natom,Nvertex,nlistsize_ver_atom,nlist_ver_atom,max_no_equiv_ver_atom,coord,coord_vertices,vert_ice_coord,simid, &
         BC1,BC2,BC3,ver_atom_dist_coord,NT_ver,atype_ver)

      call deallocate_nm_spin_ice()

   end subroutine setup_ice_neighbours

   !> deallocation of arrays
   subroutine deallocate_nm_spin_ice()

      implicit none

      integer :: i_stat, i_all

      i_all=-product(shape(nm_ver_ver))*kind(nm_ver_ver)
      deallocate(nm_ver_ver,stat=i_stat)
      call memocc(i_stat,i_all,'nm_ver_ver','deallocate_nm_spin_ice')

      i_all=-product(shape(nm_ver_atom))*kind(nm_ver_atom)
      deallocate(nm_ver_atom,stat=i_stat)
      call memocc(i_stat,i_all,'nm_ver_atom','deallocate_nm_spin_ice')

      i_all=-product(shape(nmdim_ver_ver))*kind(nmdim_ver_ver)
      deallocate(nmdim_ver_ver,stat=i_stat)
      call memocc(i_stat,i_all,'nmdim_ver_ver','deallocate_nm_spin_ice')

      i_all=-product(shape(nmdim_ver_atom))*kind(nmdim_ver_atom)
      deallocate(nmdim_ver_atom,stat=i_stat)
      call memocc(i_stat,i_all,'nmdim_ver_atom','deallocate_nm_spin_ice')

   end subroutine deallocate_nm_spin_ice

   !> Initialization of the variables needed for the Spin-Ice measurements
   subroutine spin_ice_init()

      implicit none

      real(dblprec) :: one=1.0_dblprec
      real(dblprec) :: zero=0.0_dblprec

      !Vertices
      NA_ver  = 0
      NT_ver  = 0
      ver_sym = 0
      C1_ver  = (/one,zero,zero/)
      C2_ver  = (/zero,one,zero/)
      C3_ver  = (/zero,zero,one/)
      vertex  = 'vertexfile'

      ver_no          = 0
      mchits          = 0
      vertex_step     = 1000
      vertex_buff     = 10
      prn_vertices    = 'N'
      loop_ave_len    = 0.0_dblprec
      mchits_spin_ice = 0

   end subroutine spin_ice_init

   !> Setup vertices to first supercell.
   subroutine setup_globvertices(Nvertex, NA_ver, Bas_ver, C1_ver, C2_ver, C3_ver, N1, N2, N3, atype_ver, anumb_ver, do_prnstruct, simid,coord_vertices)
      !
      !
      implicit none
      !
      integer, intent(inout) :: Nvertex !< Number of atoms in system
      integer, intent(in) :: NA_ver  !< Number of atoms in one cell
      real(dblprec), dimension(3,NA_ver) , intent(inout) :: Bas_ver !< Coordinates for basis atoms
      real(dblprec), dimension(3) , intent(in) :: C1_ver !< First lattice vector
      real(dblprec), dimension(3) , intent(in) :: C2_ver !< Second lattice vector
      real(dblprec), dimension(3) , intent(in) :: C3_ver !< Third lattice vector
      real(dblprec), dimension(:,:),allocatable,intent(inout) :: coord_vertices !< Coordinates of atoms
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, dimension(Nvertex), intent(in) :: atype_ver !< Type of atom
      integer, dimension(Nvertex), intent(in) :: anumb_ver !< Atom number in cell
      integer, intent(in) :: do_prnstruct !< Print Hamiltonian information (0/1)
      character(len=8) :: simid !< Name of simulation

      integer :: i,i_stat
      integer :: I0, i1, i2, i3
      integer :: iatom
      character(len=20) :: filn
      real(dblprec) :: detmatrix
      real(dblprec), dimension(3) :: icvec, bsf
      real(dblprec), dimension(3,3) :: invmatrix
      !
      ! Fold all atoms into first unit cell:
      ! Calculate inverse of basis matrix

      detmatrix=C1_ver(1)*C2_ver(2)*C3_ver(3)-C1_ver(1)*C2_ver(3)*C3_ver(2)+&
         C1_ver(2)*C2_ver(3)*C3_ver(1)-C1_ver(2)*C2_ver(1)*C3_ver(3)+&
         C1_ver(3)*C2_ver(1)*C3_ver(2)-C1_ver(3)*C2_ver(2)*C3_ver(1)
      invmatrix=0.0_dblprec
      if(abs(detmatrix)>dbl_tolerance) then
         invmatrix(1,1)=(C2_ver(2)*C3_ver(3)-C3_ver(2)*C2_ver(3))/detmatrix
         invmatrix(1,2)=(C1_ver(3)*C3_ver(2)-C3_ver(3)*C1_ver(2))/detmatrix
         invmatrix(1,3)=(C1_ver(2)*C2_ver(3)-C2_ver(2)*C1_ver(3))/detmatrix
         invmatrix(2,1)=(C2_ver(3)*C3_ver(1)-C3_ver(3)*C2_ver(1))/detmatrix
         invmatrix(2,2)=(C1_ver(1)*C3_ver(3)-C3_ver(1)*C1_ver(3))/detmatrix
         invmatrix(2,3)=(C1_ver(3)*C2_ver(1)-C2_ver(3)*C1_ver(1))/detmatrix
         invmatrix(3,1)=(C2_ver(1)*C3_ver(2)-C3_ver(1)*C2_ver(2))/detmatrix
         invmatrix(3,2)=(C1_ver(2)*C3_ver(1)-C3_ver(2)*C1_ver(1))/detmatrix
         invmatrix(3,3)=(C1_ver(1)*C2_ver(2)-C2_ver(1)*C1_ver(2))/detmatrix
      end if

      do I0=1,NA_ver
         !find coordinate vector in basis coordinates
         icvec(1)=Bas_ver(1,I0)*invmatrix(1,1)+Bas_ver(2,I0)*invmatrix(2,1)+Bas_ver(3,I0)*invmatrix(3,1)
         icvec(2)=Bas_ver(1,I0)*invmatrix(1,2)+Bas_ver(2,I0)*invmatrix(2,2)+Bas_ver(3,I0)*invmatrix(3,2)
         icvec(3)=Bas_ver(1,I0)*invmatrix(1,3)+Bas_ver(2,I0)*invmatrix(2,3)+Bas_ver(3,I0)*invmatrix(3,3)
         ! fold back to original cell
         bsf(1)=floor(icvec(1)+1d-7)
         bsf(2)=floor(icvec(2)+1d-7)
         bsf(3)=floor(icvec(3)+1d-7)
         !
         Bas_ver(1,I0)=Bas_ver(1,I0)-bsf(1)*C1_ver(1)-bsf(2)*C2_ver(1)-bsf(3)*C3_ver(1)
         Bas_ver(2,I0)=Bas_ver(2,I0)-bsf(1)*C1_ver(2)-bsf(2)*C2_ver(2)-bsf(3)*C3_ver(2)
         Bas_ver(3,I0)=Bas_ver(3,I0)-bsf(1)*C1_ver(3)-bsf(2)*C2_ver(3)-bsf(3)*C3_ver(3)
      end do
      ! Open file for writing the coordinates
      if(do_prnstruct==4) then
         write (filn,'(''c_ver.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn)
      end if

      iatom=0
      ! Allocate coordinate array
      allocate(coord_vertices(3,Nvertex),stat=i_stat)
      call memocc(i_stat,product(shape(coord_vertices))*kind(coord_vertices),'coord_vertices','setup_vertexcoord')
      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA_ver
                  i=I0+I1*NA_ver+I2*N1*NA_ver+I3*N2*N1*NA_ver
                  coord_vertices(1,i)=I1*C1_ver(1)+I2*C2_ver(1)+I3*C3_ver(1)+Bas_ver(1,I0)
                  coord_vertices(2,i)=I1*C1_ver(2)+I2*C2_ver(2)+I3*C3_ver(2)+Bas_ver(2,I0)
                  coord_vertices(3,i)=I1*C1_ver(3)+I2*C2_ver(3)+I3*C3_ver(3)+Bas_ver(3,I0)
                  if(do_prnstruct==4) then
                     write(ofileno,'(i7,3f12.6,2i6)') i,coord_vertices(1:3,i),atype_ver(i), anumb_ver(i)
                  end if
               end do
            end do
         end do
      end do

      if(do_prnstruct==4) then
         close(ofileno)
      end if
   end subroutine setup_globvertices

   !> Sets up the type of vertices in the system
   subroutine setup_vertex_type_and_numb(Nvertex, NA_ver, N1, N2, N3, atype_ver, anumb_ver, atype_inp_ver, anumb_inp_ver)
      !
      !
      implicit none

      integer, intent(in) :: Nvertex ! Number of atoms in system
      integer, intent(in) :: NA_ver  ! Number of atoms in one cell
      integer, intent(in) :: N1  ! Number of cell repetitions in x direction
      integer, intent(in) :: N2  ! Number of cell repetitions in y direction
      integer, intent(in) :: N3  ! Number of cell repetitions in z direction
      integer, dimension(Nvertex), intent(inout) :: atype_ver ! Type of atom
      integer, dimension(Nvertex), intent(inout) :: anumb_ver ! Atom number in cell
      integer, dimension(NA_ver), intent(inout) :: atype_inp_ver  ! Type of atom from input
      integer, dimension(NA_ver), intent(inout) :: anumb_inp_ver ! Atom number in cell from input
      !
      integer :: A1
      integer :: i0, i1, i2, i3
      !
      !
      A1=0
      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA_ver
                  A1=I0+I1*NA_ver+I2*N1*NA_ver+I3*N2*N1*NA_ver
                  atype_ver(A1)=atype_inp_ver(I0)
                  anumb_ver(A1)=anumb_inp_ver(I0)
               end do
            end do
         end do
      end do
   end subroutine setup_vertex_type_and_numb

   !> Wrapper routine for setting up structural and chemical information about the system
   subroutine setup_vertex_geometry( N1, N2, N3, do_prnstruct, simid)


      !
      implicit none

      integer, intent(in) :: N1  ! Number of cell repetitions in x direction
      integer, intent(in) :: N2  ! Number of cell repetitions in y direction
      integer, intent(in) :: N3  ! Number of cell repetitions in z direction
      integer, intent(in) :: do_prnstruct !< Print Hamiltonian information (0/1)
      character(len=8), intent(in) :: simid !< Name of simulation
      !
      integer :: i_stat
      !
      Nvertex=NA_ver*N1*N2*N3

      allocate(atype_ver(Nvertex),stat=i_stat)
      call memocc(i_stat,product(shape(atype_ver))*kind(atype_ver),'atype_ver','read_initphase')
      allocate(anumb_ver(Nvertex),stat=i_stat)
      call memocc(i_stat,product(shape(anumb_ver))*kind(anumb_ver),'anumb_ver','read_initphase')

      write (*,'(2x,a)',advance='no') 'Set up types of vertices'
      call setup_vertex_type_and_numb(Nvertex, NA_ver, N1, N2, N3, atype_ver, anumb_ver, atype_inp_ver, anumb_inp_ver)
      write (*,'(a)') ' done'

      write (*,'(2x,a)',advance='no') 'Set up global vertex coordinates'
      call setup_globvertices(Nvertex, NA_ver, Bas_ver, C1_ver, C2_ver, C3_ver, N1, N2, N3, atype_ver, anumb_ver, do_prnstruct, simid,coord_vertices)
      write (*,'(a)') ' done'

   end subroutine setup_vertex_geometry


   subroutine GEOMETRICAL_ICE(Natom,Nvertex,nlistsize_ver_atom,nlist_ver_atom,max_no_equiv_ver_atom,coord,coord_vertices,vert_ice_coord,simid, &
         BC1,BC2,BC3,ver_atom_dist_coord,NT_ver,atype_ver)

      implicit none

      integer, intent(in) :: NT_ver
      integer, intent(in) :: Natom
      integer, intent(in) :: Nvertex
      integer, intent(in) :: max_no_equiv_ver_atom
      integer, dimension(Nvertex), intent(in) :: atype_ver ! Type of vertex
      integer, dimension(Nvertex), intent(in) :: nlistsize_ver_atom
      integer, dimension(max_no_equiv_ver_atom,Nvertex) , intent(in) :: nlist_ver_atom

      real(dblprec), dimension(NT_ver,3), intent(in) :: ver_atom_dist_coord
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(3,Nvertex), intent(in) :: coord_vertices

      character(len=8), intent(in) :: simid
      character(len=1), intent(in) :: BC1 ! Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 ! Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 ! Boundary conditions in z-direction

      real(dblprec), dimension(3,Nvertex,max_no_equiv_ver_atom), intent(out) :: vert_ice_coord
      real(dblprec) :: norm,fac, new_norm
      real(dblprec), dimension(3) :: temp

      integer :: vertex, atom

      character(len=30) :: filn

      write(filn,'(''ICE_RULE.'',a,''.out'')') trim(simid)
      open (ofileno,file=filn,position='append')

      do vertex=1, Nvertex

         do atom=1, nlistsize_ver_atom(vertex)

            temp(1) = coord_vertices(1,vertex)-coord(1,nlist_ver_atom(atom,vertex))
            if (BC1=='P') then
               if((abs(temp(1)).gt.abs(ver_atom_dist_coord(atype_ver(vertex),1))).and.temp(1).gt.0.0_dblprec) then
                  temp(1) = -1.0_dblprec*ver_atom_dist_coord(atype_ver(vertex),1)
               else if ((abs(temp(1)).gt.abs(ver_atom_dist_coord(atype_ver(vertex),1))).and.temp(1).lt.0.0_dblprec) then
                  temp(1) = ver_atom_dist_coord(atype_ver(vertex),1)
               endif
            endif
            temp(2) = coord_vertices(2,vertex)-coord(2,nlist_ver_atom(atom,vertex))
            if (BC2=='P') then
               if((abs(temp(2)).gt.abs(ver_atom_dist_coord(atype_ver(vertex),2))).and.temp(2).gt.0.0_dblprec) then
                  temp(2) = -1.0_dblprec*ver_atom_dist_coord(atype_ver(vertex),2)
               else if ((abs(temp(2)).gt.abs(ver_atom_dist_coord(atype_ver(vertex),2))).and.temp(2).lt.0.0_dblprec) then
                  temp(2) = ver_atom_dist_coord(atype_ver(vertex),2)
               endif
            endif

            temp(3) = coord_vertices(3,vertex)-coord(3,nlist_ver_atom(atom,vertex))
            if (BC3=='P') then
               if((abs(temp(3)).gt.abs(ver_atom_dist_coord(atype_ver(vertex),3))).and.temp(3).gt.0.0_dblprec) then
                  temp(3) = -1.0_dblprec*ver_atom_dist_coord(atype_ver(vertex),3)
               else if ((abs(temp(3)).gt.abs(ver_atom_dist_coord(atype_ver(vertex),3))).and.temp(3).lt.0.0_dblprec) then
                  temp(3) = ver_atom_dist_coord(atype_ver(vertex),3)
               endif
            endif

            norm=temp(1)**2+temp(2)**2+temp(3)**2

            !Attention, the sign below has changed since the original implementation.
            fac=1.0_dblprec/sqrt(norm+1d-15)
            vert_ice_coord(1,vertex,atom)=temp(1)*fac
            vert_ice_coord(2,vertex,atom)=temp(2)*fac
            vert_ice_coord(3,vertex,atom)=temp(3)*fac

            new_norm=vert_ice_coord(1,vertex,atom)**2+vert_ice_coord(2,vertex,atom)**2+vert_ice_coord(3,vertex,atom)**2

            write(ofileno,'(i8,i8,3es16.8,es16.9)') vertex, nlist_ver_atom(atom,vertex), vert_ice_coord(1:3,vertex,atom),new_norm

         enddo
      enddo

      close(ofileno)

   end subroutine GEOMETRICAL_ICE


   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_spinice(ifile)
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
         10     continue
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
            !> - simid
         case('vertexfile')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            vertex=adjustl(trim(cache))

         case('prn_vertices') ! Flag for printing the vertices for loop algorithm
            read(ifile,*,iostat=i_err) prn_vertices
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('vertex_step') ! Simulation interval for printing the vertices
            read(ifile,*,iostat=i_err) vertex_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('vertex_buff')  ! Buffer size for the vertex printing
            read(ifile,*,iostat=i_err) vertex_buff
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err


         end select
      end if

      ! End of file
      if (i_errb==20) goto 20
      ! End of row
      if (i_errb==10) goto 10
   end do

   20  continue

   rewind(ifile)
   return
end subroutine read_parameters_spinice

subroutine read_vertices()
   use FileParser

   implicit none

   integer :: i_err,rd_len,i_errb, iat,i_stat
   character(len=50) :: keyword
   logical :: comment

   open(ifileno,file=trim(vertex))
   do
      10     continue
      ! Read file character for character until first whitespace
      keyword=""
      call bytereader(keyword,rd_len,ifileno,i_errb)
      ! converting Capital letters
      call caps2small(keyword)
      ! check for comment markers (currently % and #)
      comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
         (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1)
      ! Parse keyword
      keyword=trim(keyword)
      select case(keyword)

      case('na_ver')

         read(ifileno,*,iostat=i_err) NA_ver
         if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

      case('bas_ver')

         if (NA_ver==0) then
            write(*,*) 'ERROR: ','NA_ver not set before reading ',keyword
         else
            allocate(Bas_ver(3,NA_ver),stat=i_stat)
            call memocc(i_stat,product(shape(Bas_ver))*kind(Bas_ver),'Bas_ver','read_vertex_positions')
            allocate(atype_inp_ver(NA_ver),stat=i_stat)
            call memocc(i_stat,product(shape(atype_inp_ver))*kind(atype_inp_ver),'atype_inp_ver','read_vertex_positions')
            allocate(anumb_inp_ver(NA_ver),stat=i_stat)
            call memocc(i_stat,product(shape(anumb_inp_ver))*kind(anumb_inp_ver),'anumb_inp_ver','read_vertex_positions')

            do iat=1, NA_ver
               read(ifileno,*,iostat=i_err) anumb_inp_ver(iat), atype_inp_ver(iat), Bas_ver(1,iat), Bas_ver(2,iat), Bas_ver(3,iat)
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            end do
            NT_ver=maxval(atype_inp_ver)
         endif

      case('vertex_cell')
         read(ifileno,*,iostat=i_err) C1_ver(1), C1_ver(2), C1_ver(3)
         read(ifileno,*,iostat=i_err) C2_ver(1), C2_ver(2), C2_ver(3)
         read(ifileno,*,iostat=i_err) C3_ver(1), C3_ver(2), C3_ver(3)
         if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

      case('ver_atom')
         allocate(ver_atom_dist_coord(NT_ver,3),stat=i_stat)
         call memocc(i_stat,product(shape(ver_atom_dist_coord))*kind(ver_atom_dist_coord),'ver_atom_dist_coord','read_vertex_positions')

         do iat=1, NT_ver
            read(ifileno,*,iostat=i_err) ver_atom_dist_coord(iat,1), ver_atom_dist_coord(iat,2), ver_atom_dist_coord(iat,3)
         enddo

      case('ver_ver')
         allocate(ver_dist_coord(NT_ver,3),stat=i_stat)
         call memocc(i_stat,product(shape(ver_dist_coord))*kind(ver_dist_coord),'ver_dist_coord','read_vertex_positions')

         do iat=1, NT_ver
            read(ifileno,*,iostat=i_err) ver_dist_coord(iat,1), ver_dist_coord(iat,2), ver_dist_coord(iat,3)
         enddo

      case('ver_sym')
         read(ifileno,*,iostat=i_err) ver_sym
         if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err


      case default
         if(comment.or.len(trim(keyword))==0) then
         else
            print *,"Keyword '",trim(keyword),"' is not recognized"
         end if

      end select
      ! End of file
      if (i_errb==20) goto 20
      ! End of row
      if (i_errb==10) goto 10
   end do
   goto 30
   20  continue
   30  continue
   return

   close(ifileno)

end subroutine read_vertices


end module SetupSpinIce
