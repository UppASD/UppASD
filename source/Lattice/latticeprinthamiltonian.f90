!> Routines for printing information about the lattice Hamiltonian
!> and the mixed spin-lattice Hamiltonians
module LatticePrintHamiltonian
   use Parameters
   use Profiling
   use Constants


   implicit none
   public


contains


   !> Print structural information about the existing ll-coupling shells
   subroutine prn_ll_shells(simid, Natom, NT, NA, N1, N2, N3,  atype, &
         max_no_shells, max_no_equiv, redcoord, nn, nm, nmdim, &
         do_ralloy, Natom_full, Nchmax, acellnumb, acellnumbrev, achtype)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, intent(in) :: NT !< Number of types of atoms
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      real(dblprec), dimension(NT,max_no_shells,3), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings
      integer, dimension(NT), intent(in) :: nn !< Number of neighbour shells
      integer, dimension(Natom,max_no_shells,max_no_equiv), intent(in) :: nm !< Neighbour map
      integer, dimension(max_no_shells,Natom), intent(in) :: nmdim !< Dimension of neighbour map
      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, dimension(Natom_full), intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full), intent(in) :: acellnumbrev !< List for translating atom no. in actual cell to full cell
      integer, dimension(Natom_full), intent(in) :: achtype !< Chemical type of atoms (full list)

      integer :: i, j, k, l, count, i0, i1, i2, i3
      integer :: iatom, jatom
      character(len=30) :: filn

      !.. Executable statements

      ! print neighbor map
      write (filn,'(''ll_shells.'',a,''.out'')') trim(simid)
      open(18, file=filn)

      write (18,*) "Data from latticehamiltonianinit"
      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA
                  i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                  if (do_ralloy==1) i = acellnumb(i)
                  if (i==0) cycle
                  write (18,*) "------------------------------------------------------"
                  if (do_ralloy==0) then
                     write (18,10001) i
                  else
                     iatom = acellnumbrev(i)
                     write (18,10011) i, atype(i), achtype(i)
                  end if
                  do k=1,nn(atype(i))
                     write(18,10002) k, nmdim(k,i), redcoord(atype(i0),k,1),&
                        redcoord(atype(i0),k,2), redcoord(atype(i0),k,3),&
                        sqrt(sum(redcoord(atype(i0),k,:)**2))
                     write(18,10003)   nm(i,k,1:1) !!! TMP !!!
                     !write(18,10003)   nm(i,k,1:nmdim(k,i))

                     !!!  For lattice dynamics self-mapping is used !!!
                     ! new check for self-mapping
                     !do l=1,nmdim(k,i)
                     !   !print *,'prnB',i,k,l!,shape(nm)
                     !   if (nm(i,k,l)==i) then
                     !      write(*,'(1x,a,i6)') 'WARNING: Jii entry in neighbour map for atom',i
                     !      write(18,'(1x,a)') 'WARNING: Jii entry in neighbour map'
                     !   endif
                     !end do

                     if (do_ralloy==0) then
                        do j=1,NT
                           count=0
                           do l=1,1
                              !do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 if (atype(nm(i,k,l))==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write(18,10004) j, count
                        end do
                     else
                        do j=1,Nchmax
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 jatom = acellnumbrev(nm(i,k,l))
                                 if (achtype(jatom)==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write(18,10014) j, count
                        end do
                     end if
                  end do
               end do
            end do
         end do
      end do
      close(18)

      10001 format ("Atom=",i8)
      10002 format ("Shell=",i4,2x,"Number of atoms=",i4,2x,"Shell coordinates:", 4f8.4)
      10003 format ("            ",1X,5I6)
      10004 format ("            Type=",i4,2x,"Number of atoms=",i4)

      10011 format ("Atom=",i8,2x,"Type=",i4,2x,"Chtype=",i4)
      10014 format ("            Chtype=",i4,2x,"Number of atoms=",i4)

   end subroutine prn_ll_shells


   !> Print directions and strengths of LL couplings
   subroutine prn_llcoup(Natom, max_no_llneigh, ll_listsize, ll_list, ll_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_llneigh !< Calculated number of neighbours with ll interactions
      integer,dimension(Natom), intent(in) :: ll_listsize !< Size of neighbour list for LL
      integer,dimension(max_no_llneigh,Natom), intent(in) :: ll_list   !< List of neighbours for LL
      real(dblprec),dimension(9,max_no_llneigh,Natom), intent(in) :: ll_tens !< LL force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''lldata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted ll couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,ll_listsize(i)
         write (21,10002) ll_list(1:ll_listsize(i),i)
         write (21,10003) ll_tens(1:9,1:ll_listsize(i),i)
      end do
      close(21)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (9es14.6)

   end subroutine prn_llcoup


   !> Print directions and strengths of LL couplings
   subroutine prn_llphonopycoup(Natom, max_no_llneigh, ll_listsize, ll_list, ll_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_llneigh !< Calculated number of neighbours with ll interactions
      integer,dimension(Natom), intent(in) :: ll_listsize !< Size of neighbour list for LL
      integer,dimension(max_no_llneigh,Natom), intent(in) :: ll_list   !< List of neighbours for LL
      real(dblprec),dimension(9,max_no_llneigh,Natom), intent(in) :: ll_tens !< LL force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn, filn2

      !.. Executable statements
      write (filn,'(''llphonopydata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted ll couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,ll_listsize(i)
         write (21,10002) ll_list(1:ll_listsize(i),i)
         write (21,10003) ll_tens(1:9,1:ll_listsize(i),i)
      end do
      close(21)

      !.. Executable statements
      write (filn2,'(''llphonopymat.'',a,''.out'')') trim(simid)
      open(22, file=filn2)

      ! print harmonic interatomic force constants on matrix form
      do i=1,Natom
         do j=1,Natom-1
            write (22,10004,advance='no') ll_tens(1:3,j,i)
         end do
         write (22,10004) ll_tens(1:3,Natom,i)
         do j=1,Natom-1
            write (22,10004,advance='no') ll_tens(4:6,j,i)
         end do
         write (22,10004) ll_tens(4:6,Natom,i)
         do j=1,Natom-1
            write (22,10004,advance='no') ll_tens(7:9,j,i)
         end do
         write (22,10004) ll_tens(7:9,Natom,i)
      end do
      close(22)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (9es14.6)
      10004 format (1000es14.6)

   end subroutine prn_llphonopycoup


   !> Print directions and strengths of LLL couplings
   subroutine prn_lllcoup(Natom, max_no_lllneigh, lll_listsize, lll_list, lll_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_lllneigh !< Calculated number of neighbours with LLL interactions
      integer,dimension(Natom), intent(in) :: lll_listsize !< Size of neighbour list for LLL
      integer,dimension(2,max_no_lllneigh,Natom), intent(in) :: lll_list   !< List of neighbours for LLL
      real(dblprec),dimension(27,max_no_lllneigh,Natom), intent(in) :: lll_tens !< LLL force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''llldata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted lll couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,lll_listsize(i)
         write (21,10002) lll_list(1:2,1:lll_listsize(i),i)
         write (21,10003) lll_tens(1:27,1:lll_listsize(i),i)
      end do
      close(21)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (9es14.6)

   end subroutine prn_lllcoup


   !> Print directions and strengths of LLLL couplings
   subroutine prn_llllcoup(Natom, max_no_llllneigh, llll_listsize, llll_list, llll_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_llllneigh !< Calculated number of neighbours with LLLL interactions
      integer,dimension(Natom), intent(in) :: llll_listsize !< Size of neighbour list for LLLL
      integer,dimension(3,max_no_llllneigh,Natom), intent(in) :: llll_list   !< List of neighbours for LLLL
      real(dblprec),dimension(81,max_no_llllneigh,Natom), intent(in) :: llll_tens !< LLLL force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''lllldata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted llll couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,llll_listsize(i)
         write (21,10002) llll_list(1:3,1:llll_listsize(i),i)
         write (21,10003) llll_tens(1:81,1:llll_listsize(i),i)
      end do
      close(21)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (9es14.6)

   end subroutine prn_llllcoup


   !> Print structural information about the existing ml-coupling shells
   subroutine prn_ml_shells(simid, Natom, NT, NA, N1, N2, N3,  atype, &
         max_no_shells, max_no_equiv, redcoord, nn, nm, nmdim, &
         do_ralloy, Natom_full, Nchmax, acellnumb, acellnumbrev, achtype)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, intent(in) :: NT !< Number of types of atoms
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      real(dblprec), dimension(NT,max_no_shells,3), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings
      integer, dimension(NT), intent(in) :: nn !< Number of neighbour shells
      integer, dimension(Natom,max_no_shells,max_no_equiv), intent(in) :: nm !< Neighbour map
      integer, dimension(max_no_shells,Natom), intent(in) :: nmdim !< Dimension of neighbour map
      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, dimension(Natom_full), intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full), intent(in) :: acellnumbrev !< List for translating atom no. in actual cell to full cell
      integer, dimension(Natom_full), intent(in) :: achtype !< Chemical type of atoms (full list)

      integer :: i, j, k, l, count, i0, i1, i2, i3
      integer :: iatom, jatom
      character(len=30) :: filn

      !.. Executable statements

      ! print neighbor map
      write (filn,'(''ml_shells.'',a,''.out'')') trim(simid)
      open(18, file=filn)

      write (18,*) "Data from latticehamiltonianinit"
      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA
                  i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                  if (do_ralloy==1) i = acellnumb(i)
                  if (i==0) cycle
                  write (18,*) "------------------------------------------------------"
                  if (do_ralloy==0) then
                     write (18,10001) i
                  else
                     iatom = acellnumbrev(i)
                     write (18,10011) i, atype(i), achtype(i)
                  end if
                  do k=1,nn(atype(i))
                     write(18,10002) k, nmdim(k,i), redcoord(atype(i0),k,1),&
                        redcoord(atype(i0),k,2), redcoord(atype(i0),k,3),&
                        sqrt(sum(redcoord(atype(i0),k,:)**2))
                     write(18,10003)   nm(i,k,1:nmdim(k,i))

                     !!!  For lattice dynamics self-mapping is used !!!
                     ! new check for self-mapping
                     !do l=1,nmdim(k,i)
                     !   !print *,'prnB',i,k,l!,shape(nm)
                     !   if (nm(i,k,l)==i) then
                     !      write(*,'(1x,a,i6)') 'WARNING: Jii entry in neighbour map for atom',i
                     !      write(18,'(1x,a)') 'WARNING: Jii entry in neighbour map'
                     !   endif
                     !end do

                     if (do_ralloy==0) then
                        do j=1,NT
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 if (atype(nm(i,k,l))==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write(18,10004) j, count
                        end do
                     else
                        do j=1,Nchmax
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 jatom = acellnumbrev(nm(i,k,l))
                                 if (achtype(jatom)==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write(18,10014) j, count
                        end do
                     end if
                  end do
               end do
            end do
         end do
      end do
      close(18)

      10001 format ("Atom=",i8)
      10002 format ("Shell=",i4,2x,"Number of atoms=",i4,2x,"Shell coordinates:", 4f8.4)
      10003 format ("            ",1X,5I6)
      10004 format ("            Type=",i4,2x,"Number of atoms=",i4)

      10011 format ("Atom=",i8,2x,"Type=",i4,2x,"Chtype=",i4)
      10014 format ("            Chtype=",i4,2x,"Number of atoms=",i4)

   end subroutine prn_ml_shells


   !> Print directions and strengths of ML couplings
   subroutine prn_mlcoup(Natom, max_no_mlneigh, ml_listsize, ml_list, ml_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_mlneigh !< Calculated number of neighbours with ml interactions
      integer,dimension(Natom), intent(in) :: ml_listsize !< Size of neighbour list for ML
      integer,dimension(max_no_mlneigh,Natom), intent(in) :: ml_list   !< List of neighbours for ML
      real(dblprec),dimension(9,max_no_mlneigh,Natom), intent(in) :: ml_tens !< ML force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''mldata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted ml couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,ml_listsize(i)
         write (21,10002) ml_list(1:ml_listsize(i),i)
         write (21,10003) ml_tens(1:9,1:ml_listsize(i),i)
      end do
      close(21)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (9es14.6)

   end subroutine prn_mlcoup


   !> Print directions and strengths of LM couplings
   subroutine prn_lmcoup(Natom, max_no_lmneigh, lm_listsize, lm_list, lm_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_lmneigh !< Calculated number of neighbours with lm interactions
      integer,dimension(Natom), intent(in) :: lm_listsize !< Size of neighbour list for LM
      integer,dimension(max_no_lmneigh,Natom), intent(in) :: lm_list   !< List of neighbours for LM
      real(dblprec),dimension(9,max_no_lmneigh,Natom), intent(in) :: lm_tens !< LM force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''lmdata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted lm couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,lm_listsize(i)
         write (21,10002) lm_list(1:lm_listsize(i),i)
         write (21,10003) lm_tens(1:9,1:lm_listsize(i),i)
      end do
      close(21)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (9es14.6)

   end subroutine prn_lmcoup

  integer function which_neighbour(i,r_vec)
     !
     use SystemData, only: coord
     use InputData, only : Natom
     !.. Implicit declarations
     implicit none
     integer, intent(in) :: i !< Atom to find neighbour for
     real(dblprec),dimension(3), intent(in) :: r_vec !< The neighbour vector
     !
     integer :: j
     real(dblprec),dimension(3) :: r_i, r_j, r_trial
     !
     r_i=coord(:,i)
     !
     which_neighbour=0
     do j=1,Natom
        r_j=coord(:,j)
        call wrap_neighbour_vector(r_i,r_j,r_trial)
        if((r_trial(1)-r_vec(1))**2+(r_trial(2)-r_vec(2))**2+(r_trial(3)-r_vec(3))**2<1.0e-8_dblprec) then
           which_neighbour=j
           exit
        end if
     end do
     return
     !
  end function which_neighbour

  subroutine wrap_neighbour_vector(r_i,r_j,r_nn)
     !
     use InputData, only : BC1,BC2,BC3,N1,N2,N3,C1,C2,C3
     !
     implicit none
     !
     real(dblprec), dimension(3), intent(in) :: r_i   
     real(dblprec), dimension(3), intent(in) :: r_j
     real(dblprec), dimension(3), intent(out) :: r_nn
     !
     integer :: perx, pery, perz, ix, iy, iz
     real(dblprec) :: r_norm, rt_norm
     real(dblprec), dimension(3) :: r_trial
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
     ! Fix for periodicity if present
     r_nn=r_j-r_i
     r_norm=sum(r_nn*r_nn)
     do ix=-perx,perx
        do iy=-pery,pery
           do iz=-perz,perz
              r_trial(1)=r_nn(1)+ix*N1*C1(1)+iy*N2*C2(1)+iz*N3*C3(1)
              r_trial(2)=r_nn(2)+ix*N1*C1(2)+iy*N2*C2(2)+iz*N3*C3(2)
              r_trial(3)=r_nn(3)+ix*N1*C1(3)+iy*N2*C2(3)+iz*N3*C3(3)
              rt_norm=sum(r_trial*r_trial)
              if(rt_norm<r_norm) then
                 r_nn=r_trial
                 r_norm=rt_norm
              end if
           end do
        end do
     end do
  end subroutine wrap_neighbour_vector

  !> Print directions and strengths of MML couplings
  subroutine chk_mmlcoup(Natom, max_no_mmlneigh, mml_listsize, mml_list, mml_tens, simid)
    !
    use SystemData, only : coord
    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: max_no_mmlneigh !< Calculated number of neighbours with MML interactions
    integer,dimension(Natom), intent(in) :: mml_listsize !< Size of neighbour list for MML
    integer,dimension(2,max_no_mmlneigh,Natom), intent(in) :: mml_list   !< List of neighbours for MML
    real(dblprec),dimension(3,3,3,max_no_mmlneigh,Natom), intent(in) :: mml_tens !< MML force constant matrix
    character(len=8) :: simid !< Name of simulation

    !.. Local variables
    integer :: i, j, k, j_idx, al, be, ga, ip, jid, ip_idx
    real(dblprec) :: enorm, tnorm
    real(dblprec), dimension(3) :: r_i

    !.. Executable statements

    ! print neighbor list - after sort
    open(56,file='mmlcheck.out')
    write (*,'(a)',advance='no') "Checking mml couplings.   Error: "
    write(56,'(a)')  " Checking for consistency in mml-mapping according to A_ijk^abg == A_jik^bag  "
    enorm=0
    do i=1,Natom
       r_i=coord(:,i)
       write(56,'(a)') '---------------------------------------------------------------------------'
       write (56,'(a,3i4)') ' Atom: ',i
       do j_idx=1,mml_listsize(i)
          j=mml_list(1,j_idx,i)
          k=mml_list(2,j_idx,i)
          write (56,'(a,3i4,a,3i4,5x,2f10.5)') ' ijk:', i,j,k
          ip_idx=0
          do jid=1,mml_listsize(j)
             if(mml_list(1,jid,j)==i .and. mml_list(2,jid,j)==k) then
                ip_idx=jid
                write (56,'(a,3i4,10x,"|",3i4)')  '   jik found',j,mml_list(1:2,jid,j), j_idx, jid
                exit
             end if
          end do
          if(ip_idx.ne.0) then
             ip=mml_list(1,ip_idx,j)
             tnorm=0.0_dblprec
             do al=1,3
                do be=1,3
                   do ga=1,3
                      if((mml_tens(al,be,ga,j_idx,i)-mml_tens(be,al,ga,ip_idx,j))**2>0) then
                         write(56, '(a,5x,3i2,5x,2f12.6)')  '      warning!  alpha,beta,gamma:', al,be,ga, mml_tens(al,be,ga,j,i),mml_tens(be,al,ga,ip_idx,j)
                      else
                      !  write(56, '(a,2x,2i7,5x,3i2,5x,2f12.6)')  '      ok: ', i,j, al,be,ga, mml_tens(al,be,ga,j,i),mml_tens(be,al,ga,ip_idx,j)
                      end if
                      write(55,'(5i5,2f12.6)') ip_idx,ip, al, be, ga, (mml_tens(al,be,ga,j_idx,i)-mml_tens(be,al,ga,ip_idx,j))**2
                      tnorm=tnorm+(mml_tens(al,be,ga,j_idx,i)-mml_tens(be,al,ga,ip_idx,j))**2
                   end do
                end do
             end do
             if(tnorm==0.0_dblprec) then
                       write(56, '(a,2x,2i7,5x,3i2,5x,2f12.6)')  '      A_ijk^abg == A_jik^bag ok.'
             end if
             enorm=enorm+tnorm
          else
             write (56,'(a,3i4,10x,"|",3i4)')  '   jik missing', i,j,k, j_idx
             enorm=enorm+1.0_dblprec
          end if
       end do
    end do
    print *, enorm
    close(56)


  end subroutine chk_mmlcoup



   !> Print directions and strengths of MML couplings
   subroutine prn_mmlcoup(Natom, max_no_mmlneigh, mml_listsize, mml_list, mml_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_mmlneigh !< Calculated number of neighbours with MML interactions
      integer,dimension(Natom), intent(in) :: mml_listsize !< Size of neighbour list for MML
      integer,dimension(2,max_no_mmlneigh,Natom), intent(in) :: mml_list   !< List of neighbours for MML
      real(dblprec),dimension(27,max_no_mmlneigh,Natom), intent(in) :: mml_tens !< MML force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
    integer :: i,j
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''mmldata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted mml couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,mml_listsize(i)
!      write (21,10002) mml_list(1:2,1:mml_listsize(i),i)
!      write (21,10003) mml_tens(1:27,1:mml_listsize(i),i)
       do j=1,mml_listsize(i)
          write (21,10004) i,mml_list(1:2,j,i),mml_tens(1:27,j,i)
       end do
      end do
      close(21)

10001 format ("Atom=",i8,4x,"No neigh=",i7)
10002 format ("            ",1X,2I6)
!10002 format ("            ",1X,5I6)
10003 format (9es14.6)
10004 format (3x,3i6,27es14.6)

   end subroutine prn_mmlcoup


   !> Print directions and strengths of LMM couplings
   subroutine prn_lmmcoup(Natom, max_no_lmmneigh, lmm_listsize, lmm_list, lmm_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_lmmneigh !< Calculated number of neighbours with LMM interactions
      integer,dimension(Natom), intent(in) :: lmm_listsize !< Size of neighbour list for LMM
      integer,dimension(2,max_no_lmmneigh,Natom), intent(in) :: lmm_list   !< List of neighbours for LMM
      real(dblprec),dimension(27,max_no_lmmneigh,Natom), intent(in) :: lmm_tens !< LMM force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
    integer :: i,j
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''lmmdata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted lmm couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,lmm_listsize(i)
!      write (21,10002) lmm_list(1:2,1:lmm_listsize(i),i)
!      write (21,10003) lmm_tens(1:27,1:lmm_listsize(i),i)
       do j=1,lmm_listsize(i)
          write (21,10004) i,lmm_list(1:2,j,i), lmm_tens(1:27,j,i)
       end do
      end do
      close(21)

10001 format ("Atom=",i8,4x,"No neigh=",i7)
10002 format ("            ",1X,2I6)
!10002 format ("            ",1X,5I6)
10003 format (9es14.6)
10004 format (3x,3i6,27es14.6)

   end subroutine prn_lmmcoup


   !> Print directions and strengths of MMLL couplings
   subroutine prn_mmllcoup(Natom, max_no_mmllneigh, mmll_listsize, mmll_list, mmll_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_mmllneigh !< Calculated number of neighbours with MMLL interactions
      integer,dimension(Natom), intent(in) :: mmll_listsize !< Size of neighbour list for MMLL
      integer,dimension(3,max_no_mmllneigh,Natom), intent(in) :: mmll_list   !< List of neighbours for MMLL
      real(dblprec),dimension(81,max_no_mmllneigh,Natom), intent(in) :: mmll_tens !< MMLL force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''mmlldata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted mmll couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,mmll_listsize(i)
         write (21,10002) mmll_list(1:3,1:mmll_listsize(i),i)
         write (21,10003) mmll_tens(1:81,1:mmll_listsize(i),i)
      end do
      close(21)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,3I6)
      !10002 format ("            ",1X,5I6)
      10003 format (9es14.6)

   end subroutine prn_mmllcoup


   !> Print directions and strengths of LLMM couplings
   subroutine prn_llmmcoup(Natom, max_no_llmmneigh, llmm_listsize, llmm_list, llmm_tens, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_llmmneigh !< Calculated number of neighbours with LLMM interactions
      integer,dimension(Natom), intent(in) :: llmm_listsize !< Size of neighbour list for LLMM
      integer,dimension(3,max_no_llmmneigh,Natom), intent(in) :: llmm_list   !< List of neighbours for LLMM
      real(dblprec),dimension(81,max_no_llmmneigh,Natom), intent(in) :: llmm_tens !< LLMM force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''llmmdata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted llmm couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,llmm_listsize(i)
         write (21,10002) llmm_list(1:3,1:llmm_listsize(i),i)
         write (21,10003) llmm_tens(1:81,1:llmm_listsize(i),i)
      end do
      close(21)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,3I6)
      !10002 format ("            ",1X,5I6)
      10003 format (9es14.6)

   end subroutine prn_llmmcoup


end module LatticePrintHamiltonian
