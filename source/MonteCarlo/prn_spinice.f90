module prn_spinice
   use SpinIceData


   private

   public :: print_vertices, flush_vertices, buffer_vertex, prn_vertex, prn_mchits, prn_ver_ver, prn_ver_neigh

contains
   !> Wrapper routine for printing the vertex information
   subroutine print_vertices(Natom,sstep,mstep,Mensemble,emom,simid)

      use SpinIce, only : ice_rule

      implicit none

      integer, intent(in) :: Natom         !< Number of atoms in the system
      integer, intent(in) :: sstep         ! Simulation step in logarithmic scale
      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Mensemble     !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom !< Current unit magnetic moment
      character(len=8), intent(in) :: simid             !< Simulation ID

      !  Total trajectory
      if (prn_vertices=='Y') then

         if (mod(sstep-1,vertex_step)==0) then

            ! write step to buffer
            call buffer_vertex(mstep-1,Natom,Mensemble,emom)

            if (bcount_vertex==vertex_buff) then

               ! Write buffer to file
               call prn_vertex(Natom, Mensemble, simid)
               bcount_vertex=1
            endif
         else
            bcount_vertex=bcount_vertex+1
         endif

      endif

   end subroutine print_vertices

   !> Flush vertices measurement, i.e. print to a file if it is the last simulation step
   subroutine flush_vertices(Natom,Mensemble,mcnstep,simid)

      implicit none

      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: Mensemble       !< Number of ensembles
      integer, intent(in) :: mcnstep         !< Total number of simulation steps
      character(len=8), intent(in) :: simid  !< Simulation name ID

      ! All vertices are printed to file in the last iteration
      if (prn_vertices=='Y') then
         ! Write buffer to file
         bcount_vertex=bcount_vertex-1
         call prn_vertex(Natom,Mensemble,simid)
      end if

   end subroutine flush_vertices

   !> Buffer for the vertex measurementm
   subroutine buffer_vertex(mstep,Natom,Mensemble, emom)

      use SpinIce, only : ice_rule

      implicit none

      integer, intent(in) :: mstep         !< Simulation step (used for writing restart vertex)
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom !< Current unit magnetic moment

      !.. Local variables
      integer :: k_ensem, trial_vert
      real(dblprec) :: total_trial,ice_count,ice_count_up,ice_count_down

      do k_ensem=1, Mensemble

         do trial_vert=1, Nvertex
            ice_rule_buffer(trial_vert,k_ensem,bcount_vertex)= &
               ice_rule(Natom,Mensemble,k_ensem,emom,trial_vert,total_trial,ice_count,ice_count_up,ice_count_down)
            total_trial_buffer(trial_vert,k_ensem,bcount_vertex)=total_trial
            ice_count_buffer(trial_vert,k_ensem,bcount_vertex)=ice_count
         end do

      end do

      indxb_vertex(bcount_vertex)=mstep

   end subroutine buffer_vertex

   !> Print the vertex information for Spin Ice systems when using the loop algorithm
   subroutine prn_vertex(Natom, Mensemble, simid)
      !

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation

      !.. Local variables
      character(len=30) :: filn

      integer :: trial_vert,k_ensem,j_buffer

      write (filn,'(''vertex.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")

      do j_buffer=1, bcount_vertex
         do k_ensem=1, Mensemble
            do trial_vert=1, Nvertex
               if (ice_rule_buffer(trial_vert,k_ensem,j_buffer)) then
                  write (ofileno,10002) int(indxb_vertex(j_buffer)), trial_vert, ice_count_buffer(trial_vert,k_ensem,j_buffer), total_trial_buffer(trial_vert,k_ensem,j_buffer),1.0_dblprec
               else
                  write (ofileno,10002) int(indxb_vertex(j_buffer)), trial_vert, ice_count_buffer(trial_vert,k_ensem,j_buffer), total_trial_buffer(trial_vert,k_ensem,j_buffer),-1.0_dblprec
               end if
            end do
         end do
      enddo
      close(ofileno)

      return

      write(*,*) 'Error writing the vertex file'
      10002 format (i8,i8,f14.6,2es16.8)

   end subroutine prn_vertex

   !> Print accepted movements for the Monte Carlo routines, currently functioning only with thr Loop algorithm
   subroutine prn_mchits(Natom,Mensemble,Nvertex,mcnstep,simid)

      implicit none

      integer, intent(in) :: Natom     !< Number of Atoms
      integer, intent(in) :: mcnstep   !< Number of total MC steps
      integer, intent(in) :: Nvertex   !< Number of vertices
      integer, intent(in) :: Mensemble !< Number ,curr_bufferof ensembles
      character(len=8), intent(in) :: simid  !< S,curr_bufferimulation name
      !.. Local Scalar variables
      integer :: i, k, temp_ver
      real(dblprec) :: accept_rate, accept_rate_ice

      character(len=30) :: filn

      temp_ver=0

      write(filn,'(''mchits.'',a,''.out'')') trim(simid)

      open(ofileno,file=filn,position="append")
      do i=1, bcount_vertex
         do k=1, Mensemble
            accept_rate=mchits*1.0_dblprec/(mcnstep*Natom)
            accept_rate_ice=mchits_spin_ice*1.0_dblprec/(mcnstep*Natom)
            write (ofileno,'(i8,i8,es16.8,es16.8,i8,es16.8)') int(indxb_vertex(i)), k, accept_rate, accept_rate_ice, ver_no+temp_ver, loop_ave_len

         enddo
      enddo

      close(ofileno)

      return

      write(*,*) 'Error writing the mchits file'

   end subroutine prn_mchits

   subroutine prn_ver_ver(simid, Natom, Nvertex, NT, NT_ver, NA_ver, N1, N2, N3, atype, atype_ver, max_no_equiv, nncoord, nm, nmdim,neigh_type)

      implicit none

      integer, intent(in) :: Natom   ! Number of atoms in the system
      integer, intent(in) :: Nvertex ! Number of vertex in the system
      integer, intent(in) :: NA_ver  ! Number of vertex in one cell
      integer, intent(in) :: N1  ! Number of cell repetitions in x direction
      integer, intent(in) :: N2  ! Number of cell repetitions in y direction
      integer, intent(in) :: N3  ! Number of cell repetitions in z direction
      integer, dimension(Natom), intent(in) :: atype  ! Type of atoms
      integer, dimension(Nvertex), intent(in) :: atype_ver ! Type of vertex
      integer, intent(in) :: NT     ! Number of types of atoms
      integer, intent(in) :: NT_ver ! Number of types of vertices
      integer, intent(in) :: max_no_equiv ! Calculated maximum of neighbours in one shell for eighter vertex-vertex or vertex-atom
      real(dblprec), dimension(NT_ver,3), intent(in) :: nncoord ! Coordinates of neighbours for vertex-vertex (vertex-atom)
      integer, dimension(Nvertex,max_no_equiv), intent(in) :: nm ! Neighbour map for vertex-vertex (vertex-atom)
      integer, dimension(Nvertex), intent(in) :: nmdim ! Dimension of neighbour map vertex-vertex (vertex-atom)
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: neigh_type ! Type of neighbour map to be done vertex-vertex (vertex-atom)

      integer :: i, j, l, count, i0, i1, i2, i3
      character(len=30) :: filn
      i=0
      I0=0
      !.. Executable statements
      ! print neighbor map
      if ( neigh_type.eq.'V') then
         write (filn,'(''sverver.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn)
         write (ofileno,*) "Data from SpinIce"
      else if (neigh_type.eq.'A') then
         write (filn,'(''sveram.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn)
         write (ofileno,*) "Data from SpinIce"
      end if

      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA_ver
                  i=I0+I1*NA_ver+I2*N1*NA_ver+I3*N2*N1*NA_ver
                  if (i==0) cycle
                  if(neigh_type.eq.'V') then
                     write (ofileno,*) "------------------------------------------------------"
                  else if (neigh_type.eq.'A') then
                     write (ofileno,*) "------------------------------------------------------"
                  endif
                  if(neigh_type.eq.'V') then
                     write (ofileno,40001) i
                  else if(neigh_type.eq.'A') then
                     write (ofileno,40001) i
                  endif
                  if (neigh_type.eq.'V') then
                     write (ofileno,40002) 1,nmdim(i), nncoord(atype_ver(i0),1), nncoord(atype_ver(i0),2), nncoord(atype_ver(i0),3),&
                        sqrt(sum(nncoord(atype_ver(i0),:)**2))
                     write (ofileno,40003)   nm(i,1:nmdim(i))
                     ! new check for self-mapping
                     do l=1,nmdim(i)
                        if (nm(i,l)==i) then
                           write(*,'(1x,a,i6)') 'WARNING: ver_ver entry in neighbour map for vertex',i
                           write (ofileno,'(1x,a)') 'WARNING: ver_ver entry in neighbour map'
                        endif
                     end do
                  else if (neigh_type.eq.'A') then
                     write (ofileno,40428) 1,nmdim(i), nncoord(atype_ver(i0),1),nncoord(atype_ver(i0),2), nncoord(atype_ver(i0),3),&
                        sqrt(sum(nncoord(atype_ver(i0),:)**2))
                     write (ofileno,40003)   nm(i,1:nmdim(i))
                  endif
                  if (neigh_type.eq.'V') then
                     do j=1,NT_ver
                        count=0
                        do l=1,nmdim(i)
                           if (nm(i,l)/=0) then
                              if (atype_ver(nm(i,l))==j) then
                                 count=count+1
                              endif
                           end if
                        end do
                        write (ofileno,40004) j, count
                     end do
                  else if (neigh_type.eq.'A') then
                     do j=1,NT
                        count=0
                        do l=1,nmdim(i)
                           if (nm(i,l)/=0) then
                              if (atype(nm(i,l))==j) then
                                 count=count+1
                              endif
                           end if
                        end do
                        write (ofileno,40429) j, count
                     end do
                  end if
               end do
            end do
         end do
      end do
      if(neigh_type.eq.'V')then
         close(ofileno)
      else if(neigh_type.eq.'A') then
         close(ofileno)
      endif

      40001 format ("Vertex=",i8)
      40002 format ("Shell=",i4,2x,"Number of vertices=",i4,2x,"Shell coordinates:", 4f8.4)
      40003 format ("            ",1X,5I6)
      40004 format ("            Type=",i4,2x,"Number of vertices=",i4)
      40428 format ("Shell=",i4,2x,"Number of atoms=",i4,2x,"Shell coordinates:", 4f8.4)
      40429 format ("            Type=",i4,2x,"Number of atoms=",i4)

   end subroutine prn_ver_ver

   !> Print strength of exchange couplings
   subroutine prn_ver_neigh(Nvertex, max_no_neigh, nlistsize, nlist, simid, neigh_type)

      implicit none

      integer, intent(in) :: Nvertex ! Number of vertices in the system
      integer, intent(in) :: max_no_neigh ! Calculated maximum of neighbours for vertex-vertex (vertex-atom)
      integer, dimension(Nvertex), intent(in) :: nlistsize ! Size of neighbour list of vertex-vertex (vertex-atom)
      integer, dimension(max_no_neigh, Nvertex), intent(in) :: nlist ! Neighbour list for vertex vertex (vertex-atom)
      character(len=8),intent(in) :: simid !< Name of simulation
      character(len=1),intent(in) :: neigh_type ! Type of neighbouring to be done vertex-vertex or vertex atom

      !.. Local variables
      integer :: i
      character(len=20) :: filn

      !.. Executable statements

      if (neigh_type.eq.'V') then
         write (filn,'(''s1vv.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn)

         ! print neighbor list - after sort
         write (ofileno,*) "Sorted data from heisge0"
         do i=1,Nvertex
            write (ofileno,*) "----------------------------------"
            write (ofileno,40001) i,nlistsize(i)
            write (ofileno,40002) nlist(1:nlistsize(i),i)
         end do
         close(ofileno)

         40001  format ("Vertex=",i8,4x,"No vertex neigh=",i7)
         40002  format ("            ",1X,5I6)

      else if (neigh_type.eq.'A') then
         write (filn,'(''s1va.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn)

         ! print neighbor list - after sort
         write (ofileno,*) "Sorted data from heisge0"
         do i=1,Nvertex
            write (ofileno,*) "----------------------------------"
            write (ofileno,40005) i,nlistsize(i)
            write (ofileno,40006) nlist(1:nlistsize(i),i)
         end do
         close(ofileno)

         40005  format ("Vertex=",i8,4x,"No atoms neigh=",i7)
         40006  format ("            ",1X,5I6)

      end if
   end subroutine prn_ver_neigh

end module prn_spinice
