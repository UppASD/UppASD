!> Calculation of the dynamic matrix and the phonon dispersions
module LatticeSpectra
   !
   ! The phonon dispersions are written to phonspec.SIMID.out, and the
   ! energies are sorted by nsize on each row.
   !
   ! The module contains the public subroutine calculate_phondisp, and uses private subroutines
   ! calc_dynmat_elem, eigenvalue_calculation, and printEnergies. 
   !
   use Parameters
   use Constants
   use LatticeHamiltonianData,  only : ll_list, ll_listsize, ll_tens

   implicit none

   private

   public :: calculate_phondisp


contains


   subroutine calculate_phondisp(simid, Natom, Mensemble, NT, NA, N1, N2, N3, &
         C1, C2, C3, BC1, BC2, BC3, &
         atype, anumb, coord, mioninv, Bas, nq, q, &
         do_phondos, phondosfile, phondos_sigma, phondos_freq)

      character(len=8) :: simid                                          !< Name of simulation
      integer, intent(in) :: Natom                                       !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NT                                          !< Number of types of atoms
      integer, intent(in) :: NA                                          !< Number of atoms in one cell
      integer, intent(in) :: N1                                          !< Number of cell repetitions in x direction
      integer, intent(in) :: N2                                          !< Number of cell repetitions in y direction
      integer, intent(in) :: N3                                          !< Number of cell repetitions in z direction

      real(dblprec), dimension(3), intent(in) :: C1                      !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2                      !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3                      !< Third lattice vector
      character(len=1), intent(in) :: BC1                                !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2                                !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3                                !< Boundary conditions in z-direction

      integer, dimension(Natom), intent(in) :: atype                     !< Type of atom
      integer, dimension(Natom), intent(in) :: anumb                     !< Atom number in cell
      real(dblprec), dimension(3,Natom), intent(in) :: coord             !< Coordinates of atoms
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mioninv   !< Inverse ionic mass
      real(dblprec), dimension(3,NA), intent(in) :: Bas                  !< Coordinates for basis atoms
      integer, intent(in) :: nq                                          !< Number of q-points to sample
      real(dblprec), dimension(3,nq)             :: q                    !< q-points    

      character(len=1), intent(in) :: do_phondos
      character(len=30), intent(in) :: phondosfile
      real(dblprec), intent(in) :: phondos_sigma
      integer, intent(in) :: phondos_freq

      ! local arrays
      integer :: mu
      integer :: nu

      integer :: i
      integer :: I1, I2, I3 !< indices to find suitable cell to perform calculations

      integer :: countstart !< index for last atom "before" the cell where calculations are performed
      real(dblprec), dimension(3) :: zeroVector !(0,0,0) vector  to calculate J(0)

      complex(dblprec), allocatable, dimension(:,:,:) :: A !< the dynamical matrix
      complex(dblprec), allocatable, dimension(:,:,:,:,:) :: Atens !< the dynamical matrix on tensor form
      real(dblprec), allocatable, dimension(:,:) :: wres ! matrix eigenvalues
      complex(dblprec), allocatable, dimension(:,:,:) :: eigv !< matrix eigenvectors
      complex(dblprec),dimension(3,3) :: dynmat_elem 
      real(dblprec), allocatable, dimension(:,:) :: phondos

      ! for output printing
      character(LEN = 30) :: output_file
      character(LEN = 30) :: output_file2
      character(LEN = 30) :: output_file3
      character(LEN = 30) :: output_file4
      integer :: output_file_id,output_file_id2,output_file_id3,output_file_id4

      integer :: ia, ib, ja, jb

      ! if no q allocated then no AMS calculations are performed
      !if (.not. allocated(q)) then 
      !   write(*,*) 'No q-file allocated. No phonon spectra calculation is performed'
      !   goto 10 ! jumps to end of subroutine
      !end if

      ! Add calls to memory-prof
      allocate(A(3*NA,3*NA,nq)) 
      allocate(Atens(3,3,NA,NA,nq)) 
      allocate(eigv(3*NA,3*NA,nq)) 
      allocate(wres(3*NA,nq)) 
      A = 0.0_dblprec
      Atens = 0.0_dblprec
      zeroVector = 0.0_dblprec
      dynmat_elem = 0.0_dblprec
      !phondos=0.0_dblprec

      !write(*,*) 'nq ', nq
      !write(*,*) 'q ', q
      !write(*,*) 'shape ', shape(A)

      output_file = 'phondisp.'//simid//'.out'
      output_file_id = 98
      output_file2 = 'phonfcq.'//simid//'.out'
      output_file_id2 = 99
      output_file3 = 'phonev.'//simid//'.out'
      output_file_id3 = 97
      output_file4 = 'phondos.'//simid//'.out'
      output_file_id4 = 96

      ! since \Phi(q) is depending on positions (scalar product with q-vector), it is desired
      ! to perform calculations on the lattice in the middle of the system, to avoid 
      ! problems caused by periodicy. Calculations are done for the lattice starting with 
      ! atom number countstart+1
      I1 = N1/2
      I2 = N2/2
      I3 = N3/2

      countstart = 0 + I1*NA + I2*N1*NA + I3*N2*N1*NA
      ! this routine might be redundant
      ! Use the mod and the remainder to calculate countstart without if statements
      if (anumb(countstart+1)/=1) then ! could be mod(countstart,NA) ==0 aswell
         do i = 1,NA
            if (anumb(countstart+i+1)==1) then
               countstart = countstart+i
               exit
            else if (anumb(countstart-i+1)==1) then
               countstart = countstart-i
               exit
            end if
         end do
      end if
      print *,' In LatticeSpectra, countstart=',countstart

      ! A-matrix calculation, one for each q
      write (*,'(1x,a)',advance="no") "Mounting the dynamical matrix for phonon spectra calculation"
      A=0.0_dblprec
      do i=1,nq
         !print '(1x,a,i4,3f12.6)','iq',i,q(:,i)
         do mu=1,NA
            do nu=1,NA
               call calc_dynmat_elem(Natom, Mensemble, mu, nu, anumb, q(1:3,i), countstart, dynmat_elem, coord, mioninv)
               Atens(1:3,1:3,mu,nu,i) = -dynmat_elem
               write(401,*) 'Atens', mu, nu, i
               write(401,*) Atens(1:3,1:3,mu,nu,i)
               ia = (mu-1)*3+1
               ib = mu*3
               ja = (nu-1)*3+1
               jb = nu*3
               A(ia:ib,ja:jb,i) = Atens(1:3,1:3,mu,nu,i)
            end do
         end do
      end do
      write(*,'(a)') " done."

      do i=1,1
         do nu=1,NA*3
            write(501,'(100f12.6)') (real(A(mu,nu,i)),mu=1,NA*3)
         end do
         do nu=1,NA*3
            write(502,'(100f12.6)') (aimag(A(mu,nu,i)),mu=1,NA*3)
         end do
      end do

      ! eigenvalues of A (one set per NA) are written to wres
      write (*,'(1x,a)',advance='yes') "Diagonalizing A-matrix for lattice spectra calculation"
      ! Why use only real part of A?
      !A=real(A)
      call eigenvalue_calculation_lapack(A, wres, eigv, 3*na, nq)

      ! Sort eigenvalues and eigenvectors according to magnitude of eigenvalues
      do i=1,nq
         ! Replace with calls to MergeSort subroutines
         call sortEigenVVs(wres(1:3*na,i),eigv(1:3*na,1:3*na,i),3*na)
      end do

      ! print energies to output file
      call printEnergies(output_file, output_file_id, wres, nq, na)
      !call printEnergies(output_file, output_file_id, wres, nq, 3*na)
      !call printEnergies(output_file, output_file_id, wres, nq)
      !call printEnergies(output_file,output_file_id,wres,msat,tcmfa,tcrpa,1)

      ! print eigenvectors to output file
      call printEigVects(output_file3, output_file_id3, eigv, wres, q, nq, na)

      ! print phonon density of states to file
      call phondos_calc(output_file4, output_file_id4, wres, nq, NA, phondos, &
         do_phondos, phondosfile, phondos_sigma, phondos_freq)

      deallocate(A)
      deallocate(Atens)
      deallocate(eigv)
      deallocate(wres)
      write(*,'(a)') " done."
      write(*,'(1x,a)') 'Phonon dispersion calculations done.'
      10  continue
      ! stop

   end subroutine calculate_phondisp


   subroutine calc_dynmat_elem(Natom, Mensemble, mu, nu, anumb, q_vect, countstart, dynmat_elem, coord, mioninv)
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer,intent(in) :: mu
      integer,intent(in) :: nu
      integer, dimension(Natom), intent(in) :: anumb
      real(dblprec),dimension(:),intent(in) :: q_vect
      integer,intent(in) :: countstart
      complex(dblprec),dimension(3,3),intent(out) :: dynmat_elem 
      real(dblprec), dimension(3,Natom), intent(in) :: coord             !< Coordinates of atoms
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mioninv   !< Inverse ionic mass

      real(dblprec) :: qr
      real(dblprec),dimension(3) :: q_vect2pi 
      real(dblprec),dimension(3) :: dist
      integer :: j, mutemp, nutemp
      complex(dblprec) :: i

      i = (0.0_dblprec,1.0_dblprec)

      dynmat_elem = 0.0_dblprec
      q_vect2pi = q_vect*2_dblprec*pi
      !write(*,*), 'q_vect2pi ', q_vect2pi

      mutemp = mu + countstart  
      nutemp = nu + countstart 

      !write(*,*) 'mutemp ', mutemp
      !write(*,*) 'nutemp ', nutemp

      do j=1,ll_listsize(mutemp) 
         !write(401,*) 'mutemp ', mutemp, 'nutemp ', nutemp, ' j ', j, 'll_list(1,j,mutemp)', ll_list(1,j,mutemp)
         !write(*,*) 'shape(anumb) ', shape(anumb), ' shape(ll_list) ', shape(ll_list)
         if (anumb(nutemp)==anumb(ll_list(1,j,mutemp))) then
            !dist(:)=coord(1:3,mutemp)-coord(1:3,nlist(j,mutemp))
            !write(401,*) 'll_tens(1:3,1:3,j,mutemp)'
            !write(401,*) ll_tens(1:3,1:3,j,mutemp)
            !write(401,*) '-----dist------'
            !write(401,*) coord( 1:3 , mutemp )
            !write(401,*) coord( 1:3, ll_list(1,j,mutemp))
            !write(401,*) coord( 1:3 , mutemp ) - coord( 1:3, ll_list(1,j,mutemp) )
            dist(1:3) = coord( 1:3 , mutemp ) - coord( 1:3, ll_list(1,j,mutemp) )
            qr = q_vect2pi(1)*dist(1) + q_vect2pi(2)*dist(2) + q_vect2pi(3)*dist(3)
            !write (401,*) 'qr ', qr
            !dynmat_elem(1:3,1:3) = dynmat_elem(1:3,1:3) + ll_tens(1:3,1:3,j,mutemp) * exp( i * qr ) * &
            !sqrt(mioninv(mutemp,1)) * sqrt(mioninv(ll_list(1,j,mutemp),1)) / ev
            dynmat_elem(1:3,1:3) = dynmat_elem(1:3,1:3) + ll_tens(1:3,1:3,j,mutemp) * exp( i * qr ) * &
               sqrt(mioninv(mutemp,1)) * sqrt(mioninv(ll_list(1,j,mutemp),1))
            !dynmat_elem(1:3,1:3) = dynmat_elem(1:3,1:3) + ll_tens(1:3,1:3,j,mutemp) * exp( i * qr )
            !dynmat_elem(1:3,1:3) = dynmat_elem(1:3,1:3) * sqrt(mioninv(mutemp,1)) * sqrt(mioninv(ll_list(1,j,mutemp),1))
            !write(401,*) 'dynmat_elem(1:3,1:3)'
            !write(401,*) dynmat_elem(1:3,1:3)
         end if
         !dynmat_elem(1:3,1:3) = dynmat_elem(1:3,1:3) * sqrt(mioninv(mutemp,1)) * sqrt(mioninv(ll_list(1,j,mutemp),1))
         !dynmat_elem(1:3,1:3) = dynmat_elem(1:3,1:3) * mioninv(mutemp,1) *mioninv(ll_list(1,j,mutemp),1)
      end do

   end subroutine calc_dynmat_elem


   ! calculates eigenvalues for the matrix (for all k) A(:,:,k)
   subroutine eigenvalue_calculation_lapack(A,wres,eigv,na,nq)

      implicit none
      integer, intent(in) :: na !< Number of atoms in cell
      integer, intent(in) :: nq !< Number of q-points
      complex(dblprec), dimension(na,na,nq),intent(in) :: A !< A-matrix whose eigenvalues are the sought energies 
      complex(dblprec), dimension(na,na,nq),intent(inout) :: eigv !< Eigenvectors from J(q)
      real(dblprec),  dimension(na,nq) :: wres

      integer :: iq
      integer :: lwork, info
      complex(dblprec),dimension(:),allocatable :: work
      complex(dblprec),dimension(:,:),allocatable :: ctemp
      real(dblprec),dimension(:),allocatable :: rwork
      complex(dblprec),dimension(:,:),allocatable :: A_inplace
      complex(dblprec), allocatable, dimension(:) :: cwres
      complex(dblprec), allocatable, dimension(:) :: eig_ave

      lwork=2*na
      allocate(work(lwork))
      allocate(ctemp(na,na))
      allocate(A_inplace(na,na))
      allocate(rwork(2*na))
      allocate(cwres(NA))
      allocate(eig_ave(NA))

      ! eigenvalue calculations performed, energy = abs(real_part +i*imaginary part)
      do iq = 1,nq
         !print *,' AMS: q-vector:',iq
         A_inplace=A(:,:,iq)
         call zgeev('V','V',NA, A(1:NA,1:NA,iq), NA, cwres(1:NA), ctemp, na, eigv(1:NA,1:NA,iq), NA, WORK, LWORK, RWORK, INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in zgeev:',info
         end if
         !wres(1:NA,iq)=sqrt(cwres(1:NA)*conjg(cwres(1:NA)))
         !wres(1:NA,iq)=abs(real(cwres(1:NA)))
         wres(1:NA,iq)=-(real(cwres(1:NA)))
      end do

      deallocate(work)
      deallocate(rwork)
      deallocate(cwres)
      deallocate(ctemp)
      deallocate(A_inplace)
      deallocate(eig_ave)
   end subroutine eigenvalue_calculation_lapack


   subroutine printEnergies(filename, file_id, wres, nq, na)
      !subroutine printEnergies(filename, file_id, wres, msat, tcmfa, tcrpa, flag)
      implicit none

      character(LEN=*),intent(in) :: filename  
      integer,intent(in) :: file_id
      real(dblprec), allocatable, dimension(:,:) :: wres
      !real(dblprec) :: msat, tcmfa, tcrpa
      integer :: nq
      integer :: na

      integer :: i, j

      real(dblprec), dimension(3*na,nq) :: wressqrt
      real(dblprec) :: energyfac

      open(file_id,file=filename)

      !if(flag==1) then
      !   tcmfa=0.0_dblprec ; tcrpa=0.0_dblprec
      !endif
      !k=0 ; r=0

      wressqrt = 0.0_dblprec

      !C.f. VASP
      !FACTOR=SQRT(EVTOJ/((1E-10)**2)/AMTOKG)
      !W*1000*PLANK/EVTOJ/2/PI

      !energyfac = sqrt( ev / (angstrom**2 * amu) ) * 1000 * hbar / (ev * 2 * pi)
      energyfac =sqrt( mRy / (angstrom**2 * amu) ) * (1000 / ev ) * hbar ! / ( 2.0_dblprec * pi))
      write(*,*) 'energyfac ', energyfac

      ! before printing the eigenvalues of each matrix is sorted by size
      do i =1,nq
         !!call sortRealArray(wres(:,i),NA)
         do j=1,3*na
            wressqrt(j,i) = energyfac * sqrt(abs(wres(j,i))) * sign(1.0_dblprec,wres(j,i))
            !wressqrt(j,i) = energyfac * sqrt(wres(j,i))
         end do
         write(file_id,1001) i, wressqrt(1:3*na,i)
         !write(file_id,1001) i, wres(:,i)
         !if (flag==1) then
         !   !if (wres(na,i) >=0.001) then
         !   do ia=1,na
         !      if (wres(ia,i) >= 1.0d-5) then
         !         tcrpa = tcrpa + ( 1._dblprec / wres(ia,i) )
         !         r=r+1
         !      endif
         !   end do
         !   tcmfa = tcmfa + sum(wres(:,i))
         !   k = k + na
         !endif
      end do
      !if (flag==1) then
      !   tcmfa = (msat*tcmfa) / (6*k*k_bolt_ev*1000)
      !   write(*,1002) 'Tc-MFA from AMS:' , tcmfa
      !   tcrpa = ((r*msat) / (6*k_bolt_ev*1000))*(1._dblprec/tcrpa)
      !   write(*,1002) 'Tc-RPA from AMS:' , tcrpa
      !endif

      close(file_id)

      ! "safe" if someone calculates the spectra with many atoms in one cell.
      1001 format (2x,i5,2000es16.8)
      !1001 format (2x,i5,2000f18.12)
      1002 format (2x,a,f10.1)

   end subroutine printEnergies


   subroutine printEigVects(filename,file_id,eigv,wres,q,nq,na)
      implicit none

      integer,intent(in) :: nq
      integer,intent(in) :: na
      integer,intent(in) :: file_id
      complex(dblprec),dimension(na,na,nq), intent(inout) :: eigv
      real(dblprec),dimension(na,nq), intent(inout) :: wres
      real(dblprec),dimension(3,nq), intent(in) :: q
      character(LEN=*),intent(in) :: filename

      integer :: i,j,k

      open(file_id,file=filename)

      !eigv=eigv+(1.0_dblprec,0.0_dblprec)*1d-10
      ! before printing the eigenvalues of each matrix is sorted by size
      do i =1,nq
         ! Ensuring that the "reference" direction is positive (to avoid fluctuations of theta along q)
         !eigv(:,:,i)=eigv(:,:,i)*sign(1.0_dblprec,real(eigv(1,1,i)+1.0d-12))
         write(file_id,1001) "# q-point index ",i, " vector: ", q(1:3,i)
         do j=1,na
            ! Compact format : Value only for each complex number
            write(file_id,1003) wres(j,i),(( real(eigv(k,j,i))),k=1,na)
            !write(file_id,1002) ((real(eigv(k,j,i))),k=1,na)
            !!! Preferred format : Value and argument for each complex number
            !!write(file_id,1002) ((abs(eigv(k,j,i)), atan2(aimag(eigv(k,j,i)),real(eigv(k,j,i)))*180.0_dblprec/acos(-1.0_dblprec)),k=1,na)

            !write(file_id,1002) (eigv(k,j,i),(abs(eigv(k,j,i)), atan2(aimag(eigv(k,j,i)),real(eigv(k,j,i)))*180.0_dblprec/acos(-1.0_dblprec)),k=1,na)
            !write(file_id,1002) ((abs(eigv(j,k,i)), atan2(aimag(eigv(j,k,i)),real(eigv(j,k,i)))*180.0_dblprec/acos(-1.0_dblprec)),k=1,na)
            ! Optional format : Absolute value with the sign of the real part
            !write(file_id,1002) sign(abs(eigv(j,:,i)),real(eigv(j,:,i)))
            ! Optional format : The real part..
            !write(file_id,1002) real(eigv(:,j,i))
         end do
      end do

      close(file_id)

      ! "safe" if someone calculates the spectra with many atoms in one cell.
      1001 format (a,2x,i5,a,3f12.6)
      1002 format (1x,2000f12.6)
      1003 format (1x,f12.6,4x,2000f12.6)
   end subroutine printEigVects


   subroutine sortEigenVVs(array,matrix,nsize)
      ! subroutine to sort the output from cg when calculating eigenvalues,
      implicit none
      integer,intent(in) :: nsize
      real(dblprec),dimension(nsize):: array
      complex(dblprec),dimension(nsize,nsize):: matrix
      integer :: min_index,i,k
      real(dblprec) :: min_value,val_inf
      real(dblprec),dimension(:),allocatable :: sorted_array
      complex(dblprec),dimension(:),allocatable :: temp_row
      complex(dblprec),dimension(:,:),allocatable :: sorted_matrix

      allocate(sorted_array(nsize))
      allocate(temp_row(nsize))
      allocate(sorted_matrix(nsize,nsize))
      sorted_matrix=0.0_dblprec

      val_inf =  1.d10
      min_index = 1
      do k = 1,nsize
         min_value = val_inf
         temp_row=0.0_dblprec
         do i = 1,nsize
            if (array(i)<min_value) then
               min_value = array(i)
               temp_row(1:nsize) = matrix(1:nsize,i)
               min_index = i
            end if
         end do
         array(min_index) = val_inf
         sorted_array(k) = min_value
         sorted_matrix(1:nsize,k) = temp_row(1:nsize)
      end do

      array = sorted_array
      matrix = sorted_matrix

      deallocate(sorted_matrix)
      deallocate(sorted_array)
      deallocate(temp_row)

   end subroutine sortEigenVVs


   subroutine sortRealArray(array,nsize)
      ! subroutine to sort the output from cg when calculating eigenvalues,
      implicit none
      real(dblprec),dimension(:):: array
      integer,intent(in) :: nsize
      integer :: max_index,i,k
      real(dblprec) :: max_value,minus_inf
      real(dblprec),dimension(:),allocatable :: sorted_array

      allocate(sorted_array(nsize))
      minus_inf = -1.d10
      max_index = 1
      do k = 1,nsize
         max_value = minus_inf
         do i = 1,nsize
            if (array(i)>max_value) then
               max_value = array(i)
               max_index = i
            end if
         end do
         array(max_index) = minus_inf
         sorted_array(k) = max_value
      end do

      array(:) = sorted_array

   end subroutine sortRealArray


   !  subroutine printEigenvectors()
   !    ! NOT FINISHED! for printing the eigenvectors
   !    implicit none
   !    integer :: j,i,k,file_id
   !    character(LEN=20) :: filename
   !    filename='ams.eigenvectors.out'
   !
   !    file_id = 33
   !
   !    open(file_id,file=filename)
   !    do k=1,nq
   !       do i = 1,NA
   !          do j = 1,NA
   !             write(file_id,*)
   !          end do
   !       end do
   !    end do
   !    close(file_id)
   !
   !  end subroutine printEigenvectors


   !  subroutine printA(A)
   !    ! obsolete, used for troubleshooting
   !    implicit none
   !    integer :: id =33,i,j
   !    character(LEN=11) :: nom = 'Amatrix.out'
   !    complex(dblprec), allocatable, dimension(:,:,:) :: A 
   !    open(id,file=nom)
   !    do i=1,nq
   !       do j=1,NA
   !          write(id,*) A(j,:,i)
   !       end do
   !       write(id,*)
   !    end do
   !    write(*,*) 'A-matrix printed to Amatrix.out'
   !  end subroutine printA


   subroutine phondos_calc(filename, file_id, wres, nq, NA, phondos, &
         do_phondos, phondosfile, phondos_sigma, phondos_freq)

      implicit none


      character(LEN=*),intent(in) :: filename
      integer,intent(in) :: file_id
      real(dblprec), allocatable, dimension(:,:) :: wres
      integer :: nq                                                      !< Number of q-points to sample
      integer, intent(in) :: NA                                          !< Number of atoms in one cell
      real(dblprec), allocatable,dimension(:,:) :: phondos

      character(len=1), intent(in) :: do_phondos
      character(len=30), intent(in) :: phondosfile
      real(dblprec) :: phondos_sigma
      integer :: phondos_freq

      real(dblprec) :: emin, emax, deltae, fac, dummy
      integer :: i, k, ia, flines
      real(dblprec), allocatable,dimension(:,:) :: mtemp

      if(do_phondos=='A') then

         ! Open input file
         open(401, file=phondosfile)
         flines=0
         ! Pre-read file to get number of lines (frequencies)
         do
            read(401,*,end=200)  dummy,dummy,dummy
            flines=flines+1
         end do
         200    continue
         rewind(401)
         phondos_freq=flines
         allocate(phondos(phondos_freq,2))
         do k=1,phondos_freq
            read(401,*) phondos(k,1),phondos(k,2)
         enddo

         close(401)

      elseif(do_phondos=='Y') then

         ! Calculate AMS m-DOS
         open(file_id,file=filename)

         emin=minval(wres)
         emax=maxval(wres)
         deltae=(emax-emin)/(phondos_freq-1)
         allocate(phondos(phondos_freq,2))
         phondos=0.0_dblprec

         !!call sortRealArray(wres(:,i),NA)

         do i=1,phondos_freq
            phondos(i,1)=emin+(i-1)*deltae
         enddo

         do k=1,phondos_freq
            do i =1,nq     
               do ia=1,na
                  fac=exp(-(wres(ia,i)-phondos(k,1))**2/phondos_sigma)
                  phondos(k,2)=phondos(k,2)+fac
               enddo
            enddo
         enddo

         !Normalize DOS 

         phondos(:,2)=phondos(:,2)/(sum(phondos(:,2))*deltae)

         do k=1,phondos_freq
            write(file_id,1001) phondos(k,1),phondos(k,2),sum(phondos(1:k,2))*deltae
         enddo
      elseif(do_phondos=='N') then
      else    !S/B
         ! Open input file
         open(401, file=phondosfile)
         flines=0
         ! Pre-read file to get number of lines (frequencies)
         do
            read(401,*,end=300)  dummy,dummy,dummy,dummy,dummy
            flines=flines+1
         end do
         300    continue
         rewind(401)
         phondos_freq=flines
         allocate(mtemp(phondos_freq,2))
         allocate(phondos(phondos_freq,2))
         do k=1,phondos_freq
            read(401,*) phondos(k,1),mtemp(k,1),mtemp(k,2)
         enddo
         close(401)

         ! Transversal components (x and y assumed)
         do k=1,phondos_freq
            phondos(k,2)=sqrt(mtemp(k,1)**2+mtemp(k,2)**2)
         enddo
         ! Normalize
         deltae=(phondos(phondos_freq,1)-phondos(1,1))/(phondos_freq-1)
         phondos(:,2)=phondos(:,2)/(sum(phondos(:,2))*deltae)
         deallocate(mtemp)

         !   do k=1,phondos_freq
         !        write(*,1001) phondos(k,1),phondos(k,2),sum(phondos(1:k,2))*deltae
         !   enddo

      endif

      if(do_phondos.ne."N") deallocate(phondos)

      close(file_id)

      1001 format (2x,f18.12,2f18.12)
      1002 format (2x,a,f10.1)


   end subroutine phondos_calc


end module LatticeSpectra

