!> Data and routines for calculating the AMS and S(q,w) for a general non-collinear magnet within the Linear Spin Wave Theory framework
!> Authors
!> Νίκος Ντάλλης
!> Manuel Pereiro
!>
module Diamag_full
  use Parameters
  !use Matfunc
  use Constants
  
  use Profiling
  use AMS,                only :do_magdos, magdosfile, magdos_sigma, &
                                 magdos_freq
  use Qvectors,           only : qpoints, nq, q, qnumber
  use Hamiltoniandata,    only : ham
  use Systemdata,         only : anumb,coord, Landeg
  use Momentdata,         only : mmom, emom, emomM
  use InputData,          only : N1,N2,N3,C1,C2,C3,NA,Natom,simid,hfield,ham_inp                                

  
  !
  implicit none
  !
  !type(ham_t) :: ham
  integer :: nq_ext !< extended nq points
  complex(dblprec), dimension(:,:), allocatable :: h
  integer, save :: counter_dcf 
  character(len=1):: do_diamag_full
  integer::energy_step
  character(len=1):: do_correlation,do_norm_corr
  real(dblprec)::sigma_corr,temp_corr
  real(dblprec),dimension(3) :: q_vector, rot_axis

  ! public subroutines
  public :: setup_finite_hamiltonian_full,read_parameters_diamag_full
  public :: do_diamag_full
  !
 contains

  subroutine setup_diamg_full()
      implicit none
      counter_dcf = 0
      do_diamag_full='N'
      !q_vector=0.000_dblprec
      !rot_axis(1)=0.000_dblprec
      !rot_axis(2)=0.000_dblprec
      !rot_axis(3)=1.000_dblprec
    
  end subroutine setup_diamg_full
  !
  ! Set up the Hamiltonian for first magnetic cell
 subroutine setup_finite_hamiltonian_full(N1,N2,N3,NA,Natom, Mensemble, simid, emomM, mmom)
    !
    use Constants
    !
    implicit none
    ! Global variables
    integer, intent(in) :: N1             !< Number of cell repetitions in x direction
    integer, intent(in) :: N2             !< Number of cell repetitions in y direction
    integer, intent(in) :: N3             !< Number of cell repetitions in z direction
    integer, intent(in) :: NA             !< Number of atoms in one cell
    integer, intent(in) :: Natom          !< Number of atoms in system
    integer, intent(in) :: Mensemble      !< Number of ensembles
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom       !< Current magnetic moment magnitude
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
    character(len=8), intent(in) :: simid !< Name of simulation
    ! Local variables
    integer :: I1,I2,I3                   !< indices to find suitable cell to perform calculations
    integer :: countstart                 !< index for last atom "before" the cell where calculations are performed
    integer :: i_stat, ia, ja, j, hdim, la, i, lwork, info
    real(dblprec), dimension(3) :: zeroVector !(0,0,0) vector  to calculate J(0)
    real(dblprec), dimension(NA) :: eig_val, eig_valdum, eig_valdumzero
    real(dblprec), dimension(2*NA) :: eig_val2
    real(dblprec) :: msat
    real(dblprec), dimension(3,3) :: sumjq_zero
    real(dblprec)  :: emax, emin
    real(dblprec), dimension(:,:), allocatable :: sumjq_tot_zero
    real(dblprec), dimension(NA) :: lande  !< Lande gyromagnetic factor in real units
    real(dblprec), dimension(:,:), allocatable :: eig_val_k
    real(dblprec), dimension(:,:), allocatable :: eig_val_2k
    real(dblprec), dimension(:,:), allocatable :: eig_val_2kred
    real(dblprec), dimension(:,:), allocatable :: eig_val_kred
    real(dblprec), dimension(:,:), allocatable :: jq
    real(dblprec), dimension(3*NA-2) :: rwork
    complex(dblprec) :: udot, ucdot, ucdot_zero, ucdot_minus, ucdot_ani_au,ucdot_ani_av, ucdot_ani_b, vdot
    complex(dblprec), dimension(2*NA,2*NA) :: eig_vec, Tinv_mat
    complex(dblprec), dimension(3) :: ui, vi, uj, vj, ul, vl
    complex(dblprec), dimension(NA,NA) :: dot_mat, A_kdum, A_kdumzero
    complex(dblprec), dimension(:,:,:), allocatable :: A_kjq
    complex(dblprec), dimension(:,:,:), allocatable :: eig_vec_k, Tinv_mat_k
    complex(dblprec), dimension(:,:,:), allocatable :: eig_vec_kred
    complex(dblprec), dimension(:,:,:), allocatable :: A_k
    complex(dblprec), dimension(:,:,:), allocatable :: A_kminus
    complex(dblprec), dimension(:,:,:), allocatable :: B_k
    complex(dblprec), dimension(:,:), allocatable :: C_k
    complex(dblprec), dimension(:,:), allocatable :: A_kzero
    complex(dblprec), dimension(:,:,:), allocatable :: h_k
    complex(dblprec), dimension(:,:,:,:), allocatable :: S_mat
    complex(dblprec), dimension(2*NA-1) :: cwork
    real(dblprec), dimension(3) :: r1,r2,r3,b1,b2,b3,qvector
    real(dblprec) :: c1r1,c2r2,c3r3
    ! for output printing
    character(LEN = 24) :: output_file
    character(LEN = 24) :: output_file1
    character(LEN = 21) :: output_file2
    character(LEN = 20) :: output_file3
    character(LEN = 24) :: output_file4
    character(LEN = 23) :: output_file5
    !
    integer :: output_file_id,output_file_id1,output_file_id2,output_file_id3, output_file_id4, output_file_id5
    !Reciprocal lattice vectors only in direct coordinates (qpoints option D)
    if (qpoints == 'D') then
      !Definition of the ordering wave vector in direct coordinate basis (b1,b2,b3)
      !Calculate reciprocal lattice vectors
      ! r1 = C2xC3
      r1(1)=C2(2)*C3(3)-C2(3)*C3(2)
      r1(2)=C2(3)*C3(1)-C2(1)*C3(3)
      r1(3)=C2(1)*C3(2)-C2(2)*C3(1)
      ! r2 = C3xC1
      r2(1)=C3(2)*C1(3)-C3(3)*C1(2)
      r2(2)=C3(3)*C1(1)-C3(1)*C1(3)
      r2(3)=C3(1)*C1(2)-C3(2)*C1(1)
      ! r3 = C1xC2
      r3(1)=C1(2)*C2(3)-C1(3)*C2(2)
      r3(2)=C1(3)*C2(1)-C1(1)*C2(3)
      r3(3)=C1(1)*C2(2)-C1(2)*C2(1)
      ! cell volume C1*(C2xC3)
      c1r1=C1(1)*r1(1)+C1(2)*r1(2)+C1(3)*r1(3)
      c2r2=C2(1)*r2(1)+C2(2)*r2(2)+C2(3)*r2(3)
      c3r3=C3(1)*r3(1)+C3(2)*r3(2)+C3(3)*r3(3)
      ! b1=r1/(C1*r1)
      b1(1)=r1(1)/c1r1
      b1(2)=r1(2)/c1r1
      b1(3)=r1(3)/c1r1
      ! b2=r2/(C1*r1)
      b2(1)=r2(1)/c2r2
      b2(2)=r2(2)/c2r2
      b2(3)=r2(3)/c2r2
      ! b3=r3/(C1*r1)
      b3(1)=r3(1)/c3r3
      b3(2)=r3(2)/c3r3
      b3(3)=r3(3)/c3r3
      !Ordering wave vector
      qvector(1)=q_vector(1)*b1(1)+q_vector(2)*b2(1)+q_vector(3)*b3(1)
      qvector(2)=q_vector(1)*b1(2)+q_vector(2)*b2(2)+q_vector(3)*b3(2)
      qvector(3)=q_vector(1)*b1(3)+q_vector(2)*b2(3)+q_vector(3)*b3(3)
    else
      write (*,'(1x,a)',advance="yes") "Aborting LSWT-AMS calculation ..."
      write (*,'(1x,a)',advance="yes") "Reason: Only valid for flag qpoints D"
      go to 10
    end if
    ! Defining the number of qpoints
    if (do_correlation == 'Y' .or. (qvector(1) .ne. 0.0_dblprec .or. qvector(2) .ne. 0.0_dblprec .or. qvector(3) .ne. 0.0_dblprec)) then
      nq_ext=nq*3
    else
      nq_ext=nq
    end if
     
    ! dimension of hamiltonian matrix
    hdim=2*NA
    ! Redefinition of Lande factor in real units
    lande(1:NA) = 2.0d0*landeg(1:NA)
   
    !allocate arrays
    allocate(A_k(NA,NA,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(A_k))*kind(A_k),'A_k',' setup_finite_hamiltonian_full')
    allocate(A_kjq(NA,NA,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(A_kjq))*kind(A_kjq),'A_kjq',' setup_finite_hamiltonian_full')
    allocate(A_kzero(NA,NA),stat=i_stat)
    call memocc(i_stat,product(shape(A_kzero))*kind(A_kzero),'A_kzero',' setup_finite_hamiltonian_full')
    allocate(A_kminus(NA,NA,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(A_kminus))*kind(A_kminus),'A_kminus',' setup_finite_hamiltonian_full')
    allocate(B_k(NA,NA,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(B_k))*kind(B_k),'B_k',' setup_finite_hamiltonian_full')
    allocate(C_k(NA,NA),stat=i_stat)
    call memocc(i_stat,product(shape(C_k))*kind(C_k),'C_k',' setup_finite_hamiltonian_full')
    allocate(h_k(hdim,hdim,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(h_k))*kind(h_k),'h_k',' setup_finite_hamiltonian_full')
    allocate(h(hdim,hdim),stat=i_stat)
    call memocc(i_stat,product(shape(h))*kind(h),'h',' setup_finite_hamiltonian_full')  
    allocate(eig_val_k(NA,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(eig_val_k))*kind(eig_val_k),'eig_val_k',' setup_finite_hamiltonian_full')
    allocate(eig_val_2k(2*NA,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(eig_val_2k))*kind(eig_val_2k),'eig_val_2k',' setup_finite_hamiltonian_full')
    allocate(eig_val_2kred(2*NA,nq),stat=i_stat)
    call memocc(i_stat,product(shape(eig_val_2kred))*kind(eig_val_2kred),'eig_val_2kred',' setup_finite_hamiltonian_full')
    allocate(eig_val_kred(NA,nq),stat=i_stat)
    call memocc(i_stat,product(shape(eig_val_kred))*kind(eig_val_kred),'eig_val_kred',' setup_finite_hamiltonian_full')
    allocate(jq(na,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(jq))*kind(jq),'jq',' setup_finite_hamiltonian_full')
    allocate(eig_vec_k(2*NA,2*NA,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(eig_vec_k))*kind(eig_vec_k),'eig_vec_k',' setup_finite_hamiltonian_full')
    allocate(Tinv_mat_k(2*NA,2*NA,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(Tinv_mat_k))*kind(Tinv_mat_k),'Tinv_mat_k',' setup_finite_hamiltonian_full')
    allocate(eig_vec_kred(2*NA,2*NA,nq),stat=i_stat)
    call memocc(i_stat,product(shape(eig_vec_kred))*kind(eig_vec_kred),'eig_vec_kred',' setup_finite_hamiltonian_full')
    allocate(sumjq_tot_zero(NA,nq_ext),stat=i_stat)
    call memocc(i_stat,product(shape(sumjq_tot_zero))*kind(sumjq_tot_zero),'sumjq_tot_zero',' setup_finite_hamiltonian_full')
    if (do_correlation=='Y') then
    allocate(S_mat(nq,3,3,energy_step),stat=i_stat)
    call memocc(i_stat,product(shape(S_mat))*kind(S_mat),'S_mat',' setup_finite_hamiltonian_full')   
    end if
    ! Definition of output files
    output_file = 'eigval-lswt.'//simid//'.out'
    output_file_id = 104
    output_file1 = 'eigvec-lswt.'//simid//'.out'
    output_file_id1 = 105
    output_file2 = 'dos-lswt.'//simid//'.out'
    output_file_id2 = 106
    output_file3 = 'jq-lswt.'//simid//'.out'
    output_file_id3 = 107
    output_file4 = 'dcorrf-real.'//simid//'.out'
    output_file_id4 = 108
    output_file5 = 'dcorrf-img.'//simid//'.out'
    output_file_id5 = 109
    !Initialising arrays and variables
    eig_val_k=0.0d0
    eig_val_2k=0.0d0
    eig_vec_k=0.0d0
    eig_val=0.0d0
    eig_val2=0.0d0
    eig_vec=0.0d0
    msat=0.0d0
    ! since j(q) is depending on positions (scalar product with q-vector), it is desired
    ! to perform calculations on the lattice in the middle of the system, to avoid
    ! problems caused by periodicy. Calculations are done for the lattice starting with
    ! atom number countstart+1
    I1 = N1/2
    I2 = N2/2
    I3 = N3/2
    !
    countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
    zeroVector = 0.0_dblprec
    sumjq_zero = 0.0d0
    ! Calculation of the saturation magnetization
    do i=1,na
      msat=msat+mmom(i+countstart,1)/na
    enddo
    ! Write the initialization of the calculation of H(k)
    if (counter_dcf==0) then
    write (*,'(1x,a)',advance="yes") "Mounting H(k) for non-collinear AMS calculation"
    end if
    ! Toth and Lake (J. Phys.: Condens. Matter 27 (2015) 166002)
    ! Calculation of the anisotropy matrix in the rotating frame
    ! Calculation of A and B matrices
    A_k=0.0d0;A_kminus=0.0d0;B_k=0.0d0;C_k=0.0d0;dot_mat=0.0d0
    do i=1, nq
     do ia=1,NA
        call find_uv(ui,vi,emomM(:,ia,1))
         do ja=1,NA
           call find_uv(uj,vj,emomM(:,ja,1))
           if (ham_inp%do_jtensor==0) then
             ! Phason mode (i or k)
               if (i==1.and.ia==1.and.ja==1.and.counter_dcf==0) then
                 write (*,'(1x,a)',advance="yes") " Calculating the AMS in scalar J representation"
               end if
               udot=sJs(ui,jncoup_q(ia,ja,-q(:,i),qvector,countstart),uj)
               ucdot=sJs(ui,jncoup_q(ia,ja,-q(:,i),qvector,countstart),conjg(uj))
               ucdot_minus=sJs(ui,jncoup_q(ia,ja,q(:,i),qvector,countstart),conjg(uj))
            
               A_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot
               A_kminus(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot_minus
               B_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*udot
               !A_kjq(ia,ja,i)=A_k(ia,ja,i)
                !A_kjq(ia,ja,i)=-2.0d0*A_k(ia,ja,i)
              
           else if (ham_inp%do_jtensor==1) then
             ! Phason mode (i or k)
               if (i==1.and.ia==1.and.ja==1.and.counter_dcf==0) then
                 write (*,'(1x,a)',advance="yes") " Calculating the AMS in tensorial J representation"
               end if
               udot=sJs(ui,jtens_q(ia,ja,-q(:,i),qvector,countstart),uj)
               ucdot=sJs(ui,jtens_q(ia,ja,-q(:,i),qvector,countstart),conjg(uj))
               ucdot_minus=sJs(ui,jtens_q(ia,ja,q(:,i),qvector,countstart),conjg(uj))
              
               A_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot
               A_kminus(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot_minus
               B_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*udot
               
           else
               write (*,'(1x,a)',advance="yes") "  Warning: No A and B AMS matrices implemented for other options of do_jtensor flag"
           end if
         end do
     end do
    end do
    ! Definition for the ordering wave-vector
    if (nq_ext .gt. nq) then
      ! Calculating k+Q and k-Q modes
      do i=nq+1, nq_ext
       do ia=1,NA
          call find_uv(ui,vi,emomM(:,ia,1))
           do ja=1,NA
             call find_uv(uj,vj,emomM(:,ja,1))
             if (ham_inp%do_jtensor==0) then
               ! i+q_vector mode (k+Q)
               if (i .gt. nq .and. i .le. 2*nq) then
                 j=i-nq
                 udot=sJs(ui,jncoup_q(ia,ja,-q(:,j)-qvector(:),qvector,countstart),uj)
                 ucdot=sJs(ui,jncoup_q(ia,ja,-q(:,j)-qvector(:),qvector,countstart),conjg(uj))
                 ucdot_minus=sJs(ui,jncoup_q(ia,ja,q(:,j)+qvector(:),qvector,countstart),conjg(uj))
                
                 A_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot
                 A_kminus(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot_minus
                 B_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*udot
                 
               ! i-q_vector mode (k-Q)
               else
                 j=i-(2*nq)
                 udot=sJs(ui,jncoup_q(ia,ja,-q(:,j)+qvector(:),qvector,countstart),uj)
                 ucdot=sJs(ui,jncoup_q(ia,ja,-q(:,j)+qvector(:),qvector,countstart),conjg(uj))
                 ucdot_minus=sJs(ui,jncoup_q(ia,ja,q(:,j)-qvector(:),qvector,countstart),conjg(uj))
                
                 A_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot
                 A_kminus(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot_minus
                 B_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*udot
                
               end if
             else if (ham_inp%do_jtensor==1) then
               ! i+q_vector mode (k+Q)
               if (i .gt. nq .and. i .le. 2*nq) then
                 j=i-nq
                 udot=sJs(ui,jtens_q(ia,ja,-q(:,j)-qvector(:),qvector,countstart),uj)
                 ucdot=sJs(ui,jtens_q(ia,ja,-q(:,j)-qvector(:),qvector,countstart),conjg(uj))
                 ucdot_minus=sJs(ui,jtens_q(ia,ja,q(:,j)+qvector(:),qvector,countstart),conjg(uj))
                 
                 A_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot
                 A_kminus(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot_minus
                 B_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*udot
                 
               ! i-q_vector mode (k-Q)
               else
                 j=i-(2*nq)
                 udot=sJs(ui,jtens_q(ia,ja,-q(:,j)+qvector(:),qvector,countstart),uj)
                 ucdot=sJs(ui,jtens_q(ia,ja,-q(:,j)+qvector(:),qvector,countstart),conjg(uj))
                 ucdot_minus=sJs(ui,jtens_q(ia,ja,q(:,j)-qvector(:),qvector,countstart),conjg(uj))
                
                 A_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot
                 A_kminus(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot_minus
                 B_k(ia,ja,i)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*udot
                 
               end if
             else
                 write (*,'(1x,a)',advance="yes") "  Warning: No A and B AMS matrices implemented for other options of do_jtensor flag"
             end if
           end do
       end do
      end do
    end if
    !Calculation of A_kzero
     A_kzero=0.0d0
    do ia=1,NA
      call find_uv(ui,vi,emomM(:,ia,1))
          do ja=1,NA
             call find_uv(uj,vj,emomM(:,ja,1))
             if (ham_inp%do_jtensor==0) then
               ucdot_zero=sJs(ui,jncoup_q(ia,ja,zeroVector,qvector,countstart),conjg(uj))
               A_kzero(ia,ja)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot_zero
             else if (ham_inp%do_jtensor==1) then
               ucdot_zero=sJs(ui,jtens_q(ia,ja,zeroVector,qvector,countstart),conjg(uj))
               A_kzero(ia,ja)=0.5d0*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot_zero
             else
               write (*,'(1x,a)',advance="yes") "  Warning: No A_kzero matrices implemented for other options"
             end if
          end do
    end do
    ! Calculation of C matrix
       C_k=0.0d0
       do ia=1,NA
          call find_uv(ui,vi,emomM(:,ia,1))
            do la=1,NA
               call find_uv(ul,vl,emomM(:,la,1))
               if (ham_inp%do_jtensor==0) then
                 vdot=sJs(vi,jncoup_q(ia,la,zeroVector,qvector,countstart),vl)
                 C_k(ia,ia)=C_k(ia,ia)+ mmom(la,1)*vdot
               else if (ham_inp%do_jtensor==1) then
                 vdot=sJs(vi,jtens_q(ia,la,zeroVector,qvector,countstart),vl)
                 C_k(ia,ia)=C_k(ia,ia)+ mmom(la,1)*vdot
               else
                 write (*,'(1x,a)',advance="yes") "  Warning: No C AMS matrices implemented for other options"
               end if
            end do
      end do

     
   ! Setting up the Hamiltonian H(k)
   do i=1,nq_ext
    do ia=1,NA
       do ja=1,NA
          ! Diagonal as-is
          h_k(ia,ja,i)=A_k(ia,ja,i)-C_k(ia,ja)
          ! Off-diagonal as-is
          h_k(ia,ja+NA,i)=B_k(ia,ja,i)
          ! Off-diagonal Transpose and conjugate
          h_k(ia+NA,ja,i)=conjg(B_k(ja,ia,i))
          ! Diagonal Conjugate
          h_k(ia+NA,ja+NA,i)=(A_kminus(ja,ia,i)-C_k(ia,ja)) !conjg(A_kminus(ia,ja,i))-C_k(ia,ja)
       end do
    end do
   end do
   ! Hamiltonian in meV units and the negative sign is for adopting the current
   ! convenction in UppASD for the Heisenberg Hamiltonian defined with a negative sign
   do i=1,nq_ext
     do ia=1,NA
       do ja=1,NA
          ! Diagonal as-is
          h_k(ia,ja,i) = -h_k(ia,ja,i)*mub/mry*ry_ev*lande(ia)
          ! Off-diagonal upper right
          h_k(ia,ja+NA,i) = -h_k(ia,ja+NA,i)*mub/mry*ry_ev*lande(ia)
          ! Off-diagonal lower left
          h_k(ia+NA,ja,i) = -h_k(ia+NA,ja,i)*mub/mry*ry_ev*lande(ia)
          ! Diagonal
          h_k(ia+NA,ja+NA,i) = -h_k(ia+NA,ja+NA,i)*mub/mry*ry_ev*lande(ia)
       enddo
     enddo
   enddo
   ! Add delta to diagonal terms to avoid error in Cholesky decomposition close to zero energy
    do ia=1,2*NA
     h_k(ia,ia,:)= h_k(ia,ia,:)+1.0d-9
     enddo
   ! Diagonalize hamiltonian
   do i=1,nq_ext
     h(:,:) = h_k(:,:,i)
     if (i .eq. nq_ext) then
       call diagonalize_quad_hamiltonian_full(NA,eig_val,eig_val2,eig_vec,Tinv_mat,1,i)
        
       else
       call diagonalize_quad_hamiltonian_full(NA,eig_val,eig_val2,eig_vec,Tinv_mat,0,i)
        
     end if
     eig_val_k(:,i) = eig_val(:)
     eig_val_2k(:,i) = eig_val2(:)
     eig_vec_k(:,:,i) = eig_vec(:,:)
     Tinv_mat_k(:,:,i) = Tinv_mat(:,:)
   end do
   if (counter_dcf .eq. 1) then
    return
   end if
   ! Taken eigenvalues and eigenvectors up to nq with no q-vector
   if (qvector(1) == 0.0_dblprec .and. qvector(2) == 0.0_dblprec .and. qvector(3) == 0.0_dblprec) then
     do i=1,nq
       eig_val_kred(:,i) = eig_val_k(:,i)
       eig_vec_kred(:,:,i) = eig_vec_k(:,:,i)
       eig_val_2kred(:,i) = eig_val_2k(:,i)
     end do
     ! Printing eigenvectors, eigenvalues and DOS
     write(*,'(a)') "  Calculation of ordering temperature:"
     ! Print eigenvalues and eigenvector to output file
     call printeigenvalues(output_file,output_file_id,eig_val_kred,msat,1,0) 
     call printeigenvectors(output_file1,output_file_id1,eig_vec_kred,eig_val_2kred,q,nq,na,0)
     
     write(*,'(a)') "  done."
     ! Print Magnon density of states
     call lswtdos_calc(output_file2,output_file_id2,eig_val_k,0)
     ! Calculate j(q)
     ! See, for example, Eq.(7) in PRB 64, 174402
     lwork=2*NA-1
     ! J(q)
     do i=1,nq
       do ia=1,Na
         do ja=1,Na
           A_kdum(ia,ja)= A_kjq(ia,ja,i)*mmom(ia,1)
         end do
       end do
     call zheev('N','U',NA, A_kdum, NA, eig_valdum, cwork, lwork, rwork, info)
     jq(:,i) = eig_valdum(:)*(mub/mry)*(ry_ev/lande(:))
     end do
     ! J(0)
     do ia=1,Na
       do ja=1,Na
         A_kdumzero(ia,ja)= A_kzero(ia,ja)*mmom(ia,1)
       end do
     end do
     call zheev('N','U',NA, A_kdumzero, NA, eig_valdumzero, cwork, lwork, rwork, info)
     do i= 1,nq
       sumjq_tot_zero(:,i) = jq(:,i)-eig_valdumzero(:)*(mub/mry)*(ry_ev/lande(:))
     end do
     ! Print J(q) and J(q)-J(0)
     call printjq(output_file3,output_file_id3,sumjq_tot_zero,jq,0)
   else
     ! Printing eigenvectors, eigenvalues and DOS
     write(*,'(a)') "  Calculation of ordering temperature:"
     ! Print eigenvalues and eigenvector to output file
     call printeigenvalues(output_file,output_file_id,eig_val_k,msat,1,1)
     call printeigenvectors(output_file1,output_file_id1,eig_vec_k,eig_val_2k,q,nq,na,1)
     !
     write(*,'(a)') "  done."
     ! Print Magnon density of states
     call lswtdos_calc(output_file2,output_file_id2,eig_val_k,1)
     ! Calculate j(q)
     ! See, for example, Eq.(7) in PRB 64, 174402
     lwork=2*NA-1
     ! J(q)
     do i=1,nq_ext
       do ia=1,Na
         do ja=1,Na
           A_kdum(ia,ja)= A_kjq(ia,ja,i)*mmom(ia,1)
         end do
       end do
     call zheev('N','U',NA, A_kdum, NA, eig_valdum, cwork, lwork, rwork, info)
     jq(:,i) = eig_valdum(:)*(mub/mry)*(ry_ev/lande(:))
     end do
     ! J(0)
     do ia=1,Na
       do ja=1,Na
         A_kdumzero(ia,ja)= A_kzero(ia,ja)*mmom(ia,1)
       end do
     end do
     call zheev('N','U',NA, A_kdumzero, NA, eig_valdumzero, cwork, lwork, rwork, info)
     do i= 1,nq_ext
       sumjq_tot_zero(:,i) = jq(:,i)-eig_valdumzero(:)*(mub/mry)*(ry_ev/lande(:))
     end do
     ! Print J(q) and J(q)-J(0)
     call printjq(output_file3,output_file_id3,sumjq_tot_zero,jq,1)
   end if
   ! Calculation of the Dynamical correlation function
   if (do_correlation == 'Y') then
   ! 20 % bigger than maxval(eig_val_k)
   emax=aint(1.20d0*maxval(eig_val_k))
   ! To avoid division by zero in the Bose factor for the Goldstone magnon
   emin=1.0d-12
   call calc_dcf(emax,emin,energy_step,eig_val_2k,Tinv_mat_k,S_mat,qvector,Mensemble)
   
   ! Print dynamical correlation function to file
   call print_dcf(output_file4,output_file_id4,output_file5,output_file_id5,emax,emin,S_mat)
!     ! Calculation of the Powder spin wave spectra
!     if (do_powder_spectra == 'Y') then
!       if (do_correlation == 'N') then
!         write(*,'(a)') "  do_correlation flag is deactivated. Aborting Powder spin wave spectra calculation."
!         go to 10
!       else
!         ! Calculation of the powder spectra
!         call calc_pow_spectra()
!         ! Print Powder spectra to a file
!         call print_pow_spectra()
!       end if
!     end if
   end if
   ! Calculation done
   write(*,'(1x,a)') 'Non-Collinear Linear Spin-Wave Theory AMS Calculation done'
   !
10 continue
   !Deallocate arrays
   deallocate(A_k,stat=i_stat)
   call memocc(i_stat,-product(shape(A_k))*kind(A_k),'A_k',' setup_finite_hamiltonian_full')
   deallocate(A_kjq,stat=i_stat)
   call memocc(i_stat,-product(shape(A_kjq))*kind(A_kjq),'A_kjq',' setup_finite_hamiltonian_full')
   deallocate(A_kzero,stat=i_stat)
   call memocc(i_stat,product(shape(A_kzero))*kind(A_kzero),'A_kzero',' setup_finite_hamiltonian_full')
   deallocate(A_kminus,stat=i_stat)
   call memocc(i_stat,-product(shape(A_kminus))*kind(A_kminus),'A_kminus',' setup_finite_hamiltonian_full')
   deallocate(B_k,stat=i_stat)
   call memocc(i_stat,-product(shape(B_k))*kind(B_k),'B_k',' setup_finite_hamiltonian_full')
   deallocate(C_k,stat=i_stat)
   call memocc(i_stat,-product(shape(C_k))*kind(C_k),'C_k',' setup_finite_hamiltonian_full')
   deallocate(h_k,stat=i_stat)
   call memocc(i_stat,-product(shape(h_k))*kind(h_k),'h_k',' setup_finite_hamiltonian_full')
   deallocate(eig_val_k,stat=i_stat)
   call memocc(i_stat,-product(shape(eig_val_k))*kind(eig_val_k),'eig_val_k',' setup_finite_hamiltonian_full')
   deallocate(eig_val_2k,stat=i_stat)
   call memocc(i_stat,-product(shape(eig_val_2k))*kind(eig_val_2k),'eig_val_2k',' setup_finite_hamiltonian_full')
   deallocate(eig_val_2kred,stat=i_stat)
   call memocc(i_stat,-product(shape(eig_val_2kred))*kind(eig_val_2kred),'eig_val_2kred',' setup_finite_hamiltonian_full')
   deallocate(eig_val_kred,stat=i_stat)
   call memocc(i_stat,-product(shape(eig_val_kred))*kind(eig_val_kred),'eig_val_kred',' setup_finite_hamiltonian_full')
   deallocate(jq,stat=i_stat)
   call memocc(i_stat,-product(shape(jq))*kind(jq),'jq',' setup_finite_hamiltonian_full')
   deallocate(eig_vec_k,stat=i_stat)
   call memocc(i_stat,-product(shape(eig_vec_k))*kind(eig_vec_k),'eig_vec_k',' setup_finite_hamiltonian_full')
   deallocate(Tinv_mat_k,stat=i_stat)
   call memocc(i_stat,-product(shape(Tinv_mat_k))*kind(Tinv_mat_k),'Tinv_mat_k',' setup_finite_hamiltonian_full')
   deallocate(eig_vec_kred,stat=i_stat)
   call memocc(i_stat,-product(shape(eig_vec_kred))*kind(eig_vec_kred),'eig_vec_kred',' setup_finite_hamiltonian_full')
   deallocate(sumjq_tot_zero,stat=i_stat)
   call memocc(i_stat,product(shape(sumjq_tot_zero))*kind(sumjq_tot_zero),'sumjq_tot_zero',' setup_finite_hamiltonian_full')
   if (do_correlation=='Y') then
    deallocate(S_mat,stat=i_stat)
    call memocc(i_stat,-product(shape(S_mat))*kind(S_mat),'S_mat',' setup_finite_hamiltonian_full')  
   end if
  
   !
 end subroutine  setup_finite_hamiltonian_full

  subroutine diagonalize_quad_hamiltonian_full(NA,eig_val,eig_val2,eig_vec,Tinv_mat,flag,flag1)
   !
   use Constants
   !
   implicit none
   ! Global variables
   integer, intent(in) :: NA
   integer, intent(in) :: flag,flag1
   real(dblprec), dimension(NA), intent(out) :: eig_val
   real(dblprec), dimension(2*NA), intent(out) :: eig_val2
   complex(dblprec), dimension(2*NA,2*NA), intent(out) :: eig_vec
   complex(dblprec), dimension(2*NA,2*NA), intent(out) :: Tinv_mat
   ! Local variables
   integer :: i, j, info,lwork,rlwork, hdim, i_stat, ia,step,iia,jja
   real(dblprec), dimension(:), allocatable :: eig_valkgk, eig_valgl
   real(dblprec), dimension(:), allocatable :: rwork
   complex(dblprec) :: cone, czero, im
   complex(dblprec), dimension(:,:), allocatable :: g_mat
   complex(dblprec), dimension(:,:), allocatable :: U_mat
   complex(dblprec), dimension(:,:), allocatable :: Mdum_mat
   complex(dblprec), dimension(:,:), allocatable :: E_mat
   complex(dblprec), dimension(:,:), allocatable :: sqE_mat
   complex(dblprec), dimension(:,:), allocatable :: L_mat
   complex(dblprec), dimension(:,:), allocatable :: K_mat, Kinv_mat
   complex(dblprec), dimension(:,:), allocatable :: W_mat
   complex(dblprec), dimension(:,:), allocatable :: dum_mat, dum_mat1, dum_mat2
   complex(dblprec), dimension(:,:), allocatable :: eig_vecgl, eig_veckgk
   complex(dblprec), dimension(:), allocatable :: cwork
   character(LEN=1) :: unit_trimat
   ! Definition of constants
   hdim=2*NA
   lwork=2*hdim-1
   rlwork=3*hdim-2
   czero=(0.0d0,0.0d0)
   cone=(1.0d0,0.0d0)
   im=(0.0d0,1.0d0)
   ! Allocating arrays
   allocate(g_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(g_mat))*kind(g_mat),'g_mat','diagonalize_quad_hamiltonian_full')
   allocate(U_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(U_mat))*kind(U_mat),'U_mat','diagonalize_quad_hamiltonian_full')
   allocate(K_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(K_mat))*kind(K_mat),'K_mat','diagonalize_quad_hamiltonian_full')
   allocate(Kinv_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(Kinv_mat))*kind(Kinv_mat),'Kinv_mat','diagonalize_quad_hamiltonian_full')
   allocate(Mdum_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(Mdum_mat))*kind(Mdum_mat),'Mdum_mat','diagonalize_quad_hamiltonian_full')
   allocate(E_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(E_mat))*kind(E_mat),'E_mat','diagonalize_quad_hamiltonian_full')
   allocate(sqE_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(sqE_mat))*kind(sqE_mat),'sqE_mat','diagonalize_quad_hamiltonian_full')
   allocate(L_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(L_mat))*kind(L_mat),'L_mat','diagonalize_quad_hamiltonian_full')
   allocate(W_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(W_mat))*kind(W_mat),'W_mat','diagonalize_quad_hamiltonian_full')
   allocate(eig_vecgl(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(eig_vecgl))*kind(eig_vecgl),'eig_vecgl','diagonalize_quad_hamiltonian_full')
   allocate(eig_veckgk(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(eig_veckgk))*kind(eig_veckgk),'eig_veckgk','diagonalize_quad_hamiltonian_full')
   allocate(dum_mat(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(dum_mat))*kind(dum_mat),'dum_mat','diagonalize_quad_hamiltonian_full')
   allocate(dum_mat1(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(dum_mat1))*kind(dum_mat1),'dum_mat1','diagonalize_quad_hamiltonian_full')
   allocate(dum_mat2(hdim,hdim),stat=i_stat)
   call memocc(i_stat,product(shape(dum_mat2))*kind(dum_mat2),'dum_mat2','diagonalize_quad_hamiltonian_full')
   allocate(eig_valkgk(hdim),stat=i_stat)
   call memocc(i_stat,product(shape(eig_valkgk))*kind(eig_valkgk),'eig_valkgk','diagonalize_quad_hamiltonian_full')
   allocate(eig_valgl(hdim),stat=i_stat)
   call memocc(i_stat,product(shape(eig_valgl))*kind(eig_valgl),'eig_valgl','diagonalize_quad_hamiltonian_full')
   allocate(cwork(lwork),stat=i_stat)
   call memocc(i_stat,product(shape(cwork))*kind(cwork),'cwork','diagonalize_quad_hamiltonian_full')
   allocate(rwork(rlwork),stat=i_stat)
   call memocc(i_stat,product(shape(rwork))*kind(rwork),'rwork','diagonalize_quad_hamiltonian_full')
   !Deallocate arrays
   
     ! Definition of g-matrix
     g_mat=(0.0d0,0.0d0)
     do ia=1,NA
        g_mat(ia,ia)=(1.0d0,0.0d0)
        g_mat(ia+NA,ia+NA)=(-1.0d0,0.0d0)
     end do
     ! Initializing arrays
     dum_mat=0.0d0
     dum_mat1=0.0d0
     dum_mat2=0.0d0
     eig_vec=0.0d0
     U_mat=0.0d0
     K_mat=0.0d0
     Kinv_mat=0.0d0
     Tinv_mat=0.0d0
     E_mat=0.0d0
     sqE_mat=0.0d0
     L_mat=0.0d0
     W_mat=0.0d0
     eig_vecgl=0.0d0
     eig_veckgk=0.0d0
     eig_valkgk=0.0d0
     eig_valgl=0.0d0
     cwork=0.0d0
     rwork=0.0d0
     !J. H. P. Colpa, Physica A 93, 327 (1978) and C. Tsallis, Journal of Mathematical Physics 19, 277 (1978)
     ! Redefinition of the hamiltonian
    
     K_mat=h
     ! Cholesky decomposition of h to K'*K
     !call cholesky_dec('U',hdim,K_mat,hdim,info)   !here is routine from MAtfuncIf ypou want this include the Matfunct module
     call zpotrf('U',hdim,K_mat,hdim,info)         ! here is the routine from lapack directly, or MKL
     ! Setting to zero the lower triangular part of K_mat
     do i=1,hdim
       do j=1,hdim
        if ( i .gt. j) then
           K_mat(i,j) = czero
        end if
       enddo
     enddo
     ! Test for non-positive-definiteness of K_mat
     if (info .lt. 0) then
     write (*,"(a,i2,a)") 'The ',info,'-th argument had an illegal value'
     else if (info .gt. 0) then
     write (*,"(a,i1,a)") 'The leading minor in H of order ',info,' is not positive definite,'
     write (*,"(a)") 'and the factorization could not be completed.'
     write (*,"(a)") 'Either the magnetic configuration is wrong'
     write (*,"(a)") 'or there is a goldstone mode at the gamma point (0.00 0.00 0.00).'
     write (*,"(a,i4,a)") 'Problem detected in the',flag1,'-th q point'
     else if (info .eq. 0 .and. flag .eq. 1) then
        write (*,"(a)") '  Cholesky decomposition H(k)=K´K done'
     end if
     !
     ! Compute K*g
     call zgemm('N','N',hdim,hdim,hdim,cone,K_mat,hdim,g_mat,hdim,czero,dum_mat,hdim)
     ! Compute (K*g)*K'
     call zgemm('N','C',hdim,hdim,hdim,cone,dum_mat,hdim,K_mat,hdim,czero,eig_veckgk,hdim)
    

     dum_mat2=eig_veckgk
     ! Eigenvalue solver for KgK'
     call zheev('V','U',hdim,eig_veckgk,hdim,eig_valkgk, cwork, lwork, rwork, info)
     ! Test for diagonalization of  KgK'
     if (info .lt. 0) then
     write (*,"(a,i2,a)") 'The ',info,'-th argument had an illegal value'
     else if (info .gt. 0) then
     write (*,"(a,i1,a)") 'The algorithm failed to converge; ',info,' off-diagonal elements'
     write (*,"(a)") 'of an intermediate tridiagonal form did not converge to zero'
     write (*,"(a)") 'when diagonalizing KgK.'
     else if (info .eq. 0 .and. flag .eq. 1) then
     write (*,"(a)") '  Matrix diagonalization KgK´ done'
     end if
     !
     ! Calculate L
     !
     ! U'*(K*g*K')
     call zgemm('C','N',hdim,hdim,hdim,cone,eig_veckgk,hdim,dum_mat2,hdim,czero,dum_mat1,hdim)
     ! (U'*K*g*K')*U
     call zgemm('N','N',hdim,hdim,hdim,cone,dum_mat1,hdim,eig_veckgk,hdim,czero,L_mat,hdim)
     ! Calculate E=g*L
     call zgemm('N','N',hdim,hdim,hdim,cone,g_mat,hdim,L_mat,hdim,czero,E_mat,hdim)
     ! Calculate K^-1
     Kinv_mat = K_mat
     ! Checking is K_mat is a unit triangular matrix
     do i=1,hdim
      if ( K_mat(i,i)==cone ) then
      !unit triangular matrix
        unit_trimat = 'U'
      else
      !Non-unit triangular matrix
        unit_trimat = 'N'
      end if
     enddo
     call ztrtri('U',unit_trimat,hdim,Kinv_mat,hdim,info)
     ! Compute E^(1/2)
     sqE_mat=sqrt(abs(E_mat))
     ! Calculate K^-1 * U
     call zgemm('N','N',hdim,hdim,hdim,cone,Kinv_mat,hdim,eig_veckgk,hdim,czero,dum_mat,hdim)
     ! Calculate T^-1=(K^-1 * U) * E^1/2
     call zgemm('N','N',hdim,hdim,hdim,cone,dum_mat,hdim,sqE_mat,hdim,czero,Tinv_mat,hdim)
     ! Compute T'^-1 * H
     call zgemm('C','N',hdim,hdim,hdim,cone,Tinv_mat,hdim,h,hdim,czero,dum_mat,hdim)
     ! Calculate (T'^-1 * H) * T^-1
     call zgemm('N','N',hdim,hdim,hdim,cone,dum_mat,hdim,Tinv_mat,hdim,czero,W_mat,hdim)
     ! Eigensolve the problem
     eig_vecgl=W_mat
     call zheev('V','U',hdim,eig_vecgl,hdim,eig_valgl, cwork, lwork, rwork, info)
        ! Test for diagonalization of (T'^-1 * H) * T^-1
     if (info .lt. 0) then
     write (*,"(a,i2,a)") 'The ',info,'-th argument had an illegal value'
     else if (info .gt. 0) then
     write (*,"(a,i1,a)") 'The algorithm failed to converge; ',info,' off-diagonal elements'
     write (*,"(a)") 'of an intermediate tridiagonal form did not converge to zero'
     write (*,"(a)") 'when diagonalizing (T**H)^-1 * H * T^-1.'
     else if (info .eq. 0 .and. flag .eq. 1) then
     write (*,"(a)") '  Diagonalization of the eigenvalue problem (T**H)^-1*H(k)*T^-1 done'
     end if
     
     eig_val2=eig_valgl
    
     ! Sorting out the eigenvalues in increasing order
     call dlasrt( 'I', hdim, eig_valgl, info )
     step=0
     do i = 1,NA
        eig_val(i)=eig_valgl(i+step)
        step=step+1
     enddo
     ! Assign the eigenvectors
     eig_vec=eig_vecgl
   


     deallocate(g_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(g_mat))*kind(g_mat),'g_mat','diagonalize_quad_hamiltonian_full')
     deallocate(U_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(U_mat))*kind(U_mat),'U_mat','diagonalize_quad_hamiltonian_full')
     deallocate(K_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(K_mat))*kind(K_mat),'K_mat','diagonalize_quad_hamiltonian_full')
     deallocate(Kinv_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(Kinv_mat))*kind(Kinv_mat),'Kinv_mat','diagonalize_quad_hamiltonian_full')
     deallocate(Mdum_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(Mdum_mat))*kind(Mdum_mat),'Mdum_mat','diagonalize_quad_hamiltonian_full')
     deallocate(E_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(E_mat))*kind(E_mat),'E_mat','diagonalize_quad_hamiltonian_full')
     deallocate(sqE_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(sqE_mat))*kind(sqE_mat),'sqE_mat','diagonalize_quad_hamiltonian_full')
     deallocate(L_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(L_mat))*kind(L_mat),'L_mat','diagonalize_quad_hamiltonian_full')
     deallocate(W_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(W_mat))*kind(W_mat),'W_mat','diagonalize_quad_hamiltonian_full')
     deallocate(eig_vecgl,stat=i_stat)
     call memocc(i_stat,-product(shape(eig_vecgl))*kind(eig_vecgl),'eig_vecgl','diagonalize_quad_hamiltonian_full')
     deallocate(eig_veckgk,stat=i_stat)
     call memocc(i_stat,-product(shape(eig_veckgk))*kind(eig_veckgk),'eig_veckgk','diagonalize_quad_hamiltonian_full')
     deallocate(dum_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(dum_mat))*kind(dum_mat),'dum_mat','diagonalize_quad_hamiltonian_full')
     deallocate(dum_mat1,stat=i_stat)
     call memocc(i_stat,-product(shape(dum_mat1))*kind(dum_mat1),'dum_mat1','diagonalize_quad_hamiltonian_full')
     deallocate(dum_mat2,stat=i_stat)
     call memocc(i_stat,-product(shape(dum_mat2))*kind(dum_mat2),'dum_mat2','diagonalize_quad_hamiltonian_full')
     deallocate(eig_valkgk,stat=i_stat)
     call memocc(i_stat,-product(shape(eig_valkgk))*kind(eig_valkgk),'eig_valkgk','diagonalize_quad_hamiltonian_full')
     deallocate(eig_valgl,stat=i_stat)
     call memocc(i_stat,-product(shape(eig_valgl))*kind(eig_valgl),'eig_valgl','diagonalize_quad_hamiltonian_full')
     deallocate(cwork,stat=i_stat)
     call memocc(i_stat,-product(shape(cwork))*kind(cwork),'cwork','diagonalize_quad_hamiltonian_full')
     deallocate(rwork,stat=i_stat)
     call memocc(i_stat,-product(shape(rwork))*kind(rwork),'rwork','diagonalize_quad_hamiltonian_full')
     
   !
 end subroutine diagonalize_quad_hamiltonian_full


 subroutine calc_dcf(emax,emin,energy_step,eig_val_2k,Tinv_mat_k,S_mat,qvector,Mensemble)
   !
   use Constants
   !
   implicit none
   ! Global variables
   integer, intent(in) :: energy_step, Mensemble
   real(dblprec), intent(in) :: emax, emin
   real(dblprec), dimension(2*NA,nq_ext), intent(in) :: eig_val_2k
   real(dblprec), dimension(3), intent(in) :: qvector
   complex(dblprec), dimension(2*NA,2*NA,nq_ext), intent(in) :: Tinv_mat_k
   complex(dblprec), allocatable,dimension(:,:,:,:) ::  S_mat
   ! Local variables
   integer :: i, j, ia, ja, alpha, beta, hdim, l, i_stat, step,ii,jj,kk
   real(dblprec) :: a, b
   real(dblprec), dimension(3) :: dist
   real(dblprec), allocatable, dimension(:,:) :: eig_val_2k_redplus, eig_val_2k_redminus
   complex(dblprec) :: i_complex, cone, czero
   complex(dblprec), dimension(3) :: ui, vi, uj, vj
   complex(dblprec), dimension(3,3) :: R1_mat, R1conj_mat
   complex(dblprec), dimension(3,3) :: I_mat, R2_mat, K_mat
   complex(dblprec), allocatable,dimension(:,:,:,:,:) :: Y_mat, Z_mat, V_mat, O_mat
   complex(dblprec), allocatable,dimension(:,:,:,:,:) :: M_mat, Nab_mat
   complex(dblprec), allocatable,dimension(:,:,:,:,:) :: Nab_mat_red
   complex(dblprec), allocatable,dimension(:,:,:,:) :: S_ab, S_ab_plus, S_ab_minus
   complex(dblprec), allocatable,dimension(:,:,:)  :: T_mat
   complex(dblprec), allocatable,dimension(:,:) :: N_mat, Mdum1_mat, Tdum_mat, dum1_mat
   ! Allocating arrays
   
   !allocate(S_mat(nq,3,3,energy_step),stat=i_stat)
   !call memocc(i_stat,product(shape(S_mat))*kind(S_mat),'S_mat','calc_dcf')
   allocate(S_ab(nq,3,3,energy_step),stat=i_stat)
   call memocc(i_stat,product(shape(S_ab))*kind(S_ab),'S_ab','calc_dcf')
   allocate(S_ab_plus(nq,3,3,energy_step),stat=i_stat)
   call memocc(i_stat,product(shape(S_ab_plus))*kind(S_ab_plus),'S_ab_plus','calc_dcf')
   allocate(S_ab_minus(nq,3,3,energy_step),stat=i_stat)
   call memocc(i_stat,product(shape(S_ab_minus))*kind(S_ab_minus),'S_ab_minus','calc_dcf')
   allocate(Y_mat(NA,NA,3,3,nq_ext),stat=i_stat)
   call memocc(i_stat,product(shape(Y_mat))*kind(Y_mat),'Y_mat','calc_dcf')
   allocate(Z_mat(NA,NA,3,3,nq_ext),stat=i_stat)
   call memocc(i_stat,product(shape(Z_mat))*kind(Z_mat),'Z_mat','calc_dcf')
   allocate(V_mat(NA,NA,3,3,nq_ext),stat=i_stat)
   call memocc(i_stat,product(shape(V_mat))*kind(V_mat),'V_mat','calc_dcf')
   allocate(O_mat(NA,NA,3,3,nq_ext),stat=i_stat)
   call memocc(i_stat,product(shape(O_mat))*kind(O_mat),'O_mat','calc_dcf')
   allocate(M_mat(NA*2,NA*2,3,3,nq_ext),stat=i_stat)
   call memocc(i_stat,product(shape(M_mat))*kind(M_mat),'M_mat','calc_dcf')
   allocate(Nab_mat(NA*2,NA*2,3,3,nq_ext),stat=i_stat)
   call memocc(i_stat,product(shape(Nab_mat))*kind(Nab_mat),'Nab_mat','calc_dcf')
   allocate(Nab_mat_red(NA,NA,3,3,nq_ext),stat=i_stat)
   call memocc(i_stat,product(shape(Nab_mat_red))*kind(Nab_mat_red),'Nab_mat_red','calc_dcf')
   allocate(T_mat(NA*2,NA*2,nq_ext),stat=i_stat)
   call memocc(i_stat,product(shape(T_mat))*kind(T_mat),'T_mat','calc_dcf')
   allocate(N_mat(NA*2,NA*2),stat=i_stat)
   call memocc(i_stat,product(shape(N_mat))*kind(N_mat),'N_mat','calc_dcf')
   allocate(Mdum1_mat(NA*2,NA*2),stat=i_stat)
   call memocc(i_stat,product(shape(Mdum1_mat))*kind(Mdum1_mat),'Mdum1_mat','calc_dcf')
   allocate(Tdum_mat(NA*2,NA*2),stat=i_stat)
   call memocc(i_stat,product(shape(Tdum_mat))*kind(Tdum_mat),'Tdum_mat','calc_dcf')
   allocate(dum1_mat(NA*2,NA*2),stat=i_stat)
   call memocc(i_stat,product(shape(dum1_mat))*kind(dum1_mat),'dum1_mat','calc_dcf')
   allocate(eig_val_2k_redplus(2*NA,nq_ext),stat=i_stat)  !!!!!
   call memocc(i_stat,product(shape(eig_val_2k_redplus))*kind(eig_val_2k_redplus),'eig_val_2k_redplus','calc_dcf')
   allocate(eig_val_2k_redminus(2*NA,nq_ext),stat=i_stat)  !!!!!
   call memocc(i_stat,product(shape(eig_val_2k_redminus))*kind(eig_val_2k_redminus),'eig_val_2k_redminus','calc_dcf')
   
   
     ! dimension of hamiltonian matrix
     hdim=2*NA
     ! Initializing arrays
     Y_mat=0.0d0
     Z_mat=0.0d0
     V_mat=0.0d0
     O_mat=0.0d0
     N_mat=0.0d0
     dist=0.0d0
     S_ab=0.0d0
     I_mat=0.0d0
     i_complex = (0.0d0,1.0d0)
     cone=(1.0d0,0.0d0)
     czero=(0.0d0,0.0d0)
     S_ab=0.0d0
     S_ab_plus=0.0d0
     S_ab_minus=0.0d0
     S_mat=0.0d0
     dum1_mat=0.0d0
     T_mat=0.0d0
     Mdum1_mat=0.0d0
     Tdum_mat=0.0d0
     ! Calculation of the Dynamical Correlation Function
     !
     ! Step 1: Calculation of sub-matrices Y, Z, V and O
     do i=1,nq_ext
       do ia=1,NA
          call find_uv(ui,vi,emomM(:,ia,1))
          do ja=1,NA
          call find_uv(uj,vj,emomM(:,ja,1))
            do alpha=1,3
              do beta=1,3
                dist(:)=coord(1:3,ia)-coord(1:3,ja)
                ! Phason mode (i or k)
                if (i .le. nq) then
                  Y_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*(q(1,i)*2_dblprec*pi*dist(1)+q(2,i)*2_dblprec*pi*dist(2)+&
                                              q(3,i)*2_dblprec*pi*dist(3)))*ui(alpha)*conjg(uj(beta))
                  Z_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*(q(1,i)*2_dblprec*pi*dist(1)+q(2,i)*2_dblprec*pi*dist(2)+&
                                              q(3,i)*2_dblprec*pi*dist(3)))*ui(alpha)*uj(beta)
                  V_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*(q(1,i)*2_dblprec*pi*dist(1)+q(2,i)*2_dblprec*pi*dist(2)+&
                                              q(3,i)*2_dblprec*pi*dist(3)))*conjg(ui(alpha))*conjg(uj(beta))
                  O_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*(q(1,i)*2_dblprec*pi*dist(1)+q(2,i)*2_dblprec*pi*dist(2)+&
                                              q(3,i)*2_dblprec*pi*dist(3)))*conjg(ui(alpha))*uj(beta)
                ! i+q_vector mode (k+Q)
                else if (i .gt. nq .and. i .le. 2*nq) then
                  j=i-nq
                  Y_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*((q(1,j)+q_vector(1))*2_dblprec*pi*dist(1)+(q(2,j)+q_vector(2))*2_dblprec*pi*dist(2)+&
                                              (q(3,j)+q_vector(3))*2_dblprec*pi*dist(3)))*ui(alpha)*conjg(uj(beta))
                  Z_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*((q(1,j)+q_vector(1))*2_dblprec*pi*dist(1)+(q(2,j)+q_vector(2))*2_dblprec*pi*dist(2)+&
                                              (q(3,j)+q_vector(3))*2_dblprec*pi*dist(3)))*ui(alpha)*uj(beta)
                  V_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*((q(1,j)+q_vector(1))*2_dblprec*pi*dist(1)+(q(2,j)+q_vector(2))*2_dblprec*pi*dist(2)+&
                                              (q(3,j)+q_vector(3))*2_dblprec*pi*dist(3)))*conjg(ui(alpha))*conjg(uj(beta))
                  O_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*((q(1,j)+q_vector(1))*2_dblprec*pi*dist(1)+(q(2,j)+q_vector(2))*2_dblprec*pi*dist(2)+&
                                              (q(3,j)+q_vector(3))*2_dblprec*pi*dist(3)))*conjg(ui(alpha))*uj(beta)
                ! i-q_vector mode (k-Q)
                else
                  j=i-(2*nq)
                  Y_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*((q(1,j)-q_vector(1))*2_dblprec*pi*dist(1)+(q(2,j)-q_vector(2))*2_dblprec*pi*dist(2)+&
                                              (q(3,j)-q_vector(3))*2_dblprec*pi*dist(3)))*ui(alpha)*conjg(uj(beta))
                  Z_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*((q(1,j)-q_vector(1))*2_dblprec*pi*dist(1)+(q(2,j)-q_vector(2))*2_dblprec*pi*dist(2)+&
                                              (q(3,j)-q_vector(3))*2_dblprec*pi*dist(3)))*ui(alpha)*uj(beta)
                  V_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*((q(1,j)-q_vector(1))*2_dblprec*pi*dist(1)+(q(2,j)-q_vector(2))*2_dblprec*pi*dist(2)+&
                                              (q(3,j)-q_vector(3))*2_dblprec*pi*dist(3)))*conjg(ui(alpha))*conjg(uj(beta))
                  O_mat(ia,ja,alpha,beta,i) = sqrt(mmom(ia,1)*mmom(ja,1))*exp(i_complex*((q(1,j)-q_vector(1))*2_dblprec*pi*dist(1)+(q(2,j)-q_vector(2))*2_dblprec*pi*dist(2)+&
                                              (q(3,j)-q_vector(3))*2_dblprec*pi*dist(3)))*conjg(ui(alpha))*uj(beta)
                end if
              end do
            end do
          end do
       end do
     end do
     ! Step 2: Building the sub-matrices in one single matrix
     do i=1,nq_ext
      do ia=1,NA
        do ja=1,NA
          do alpha=1,3
            do beta=1,3
              ! First block (top-left)
              M_mat(ia,ja,alpha,beta,i)=Y_mat(ia,ja,alpha,beta,i)
              ! Second block (top-right)
              M_mat(ia,ja+NA,alpha,beta,i)=Z_mat(ia,ja,alpha,beta,i)
              ! Third block (bottom-left)
              M_mat(ia+NA,ja,alpha,beta,i)=V_mat(ia,ja,alpha,beta,i)
              ! Fourth block (bottom-right)
              M_mat(ia+NA,ja+NA,alpha,beta,i)=O_mat(ia,ja,alpha,beta,i)
            end do
          end do
        end do
      end do
     end do
     ! Step 3: Calculation of T_mat
     do i=1,nq_ext
       do ia=1,NA
         do ja=1,NA
           ! Transpose is done by changing index
           ! First block (top-left)
           T_mat(ia,ja,i)=conjg(Tinv_mat_k(ja,ia,i))
           ! Second block (top-right)
           T_mat(ia,ja+NA,i)=-conjg(Tinv_mat_k(ja,ia+NA,i))
           ! Third block (bottom-left)
           T_mat(ia+NA,ja,i)=-conjg(Tinv_mat_k(ja+NA,ia,i))
           ! Fourth block (bottom-right)
           T_mat(ia+NA,ja+NA,i)=conjg(Tinv_mat_k(ja+NA,ia+NA,i))
         end do
       end do
     end do
     ! Step 4: Calculation of T**H*M*T for every alpha and beta
     do i=1,nq_ext
       do alpha=1,3
         do beta=1,3
           Mdum1_mat(:,:) = M_mat(:,:,alpha,beta,i)
           Tdum_mat(:,:)=T_mat(:,:,i)
           ! Step 4: Calculation of T**H*M
           call zgemm('C','N',hdim,hdim,hdim,cone,Tdum_mat,hdim,Mdum1_mat,hdim,czero,dum1_mat,hdim)
           ! Setp 5: Calculation of T**H*M*T
           call zgemm('N','N',hdim,hdim,hdim,cone,dum1_mat,hdim,Tdum_mat,hdim,czero,N_mat,hdim)
           Nab_mat(:,:,alpha,beta,i)=N_mat(:,:)
         end do
       end do
     end do
     ! Sorting modes in decreasing order
     step=0
     do ia = 1,NA
       ! for positive q vectors
       Nab_mat_red(ia,ia,:,:,:)=Nab_mat(ia+step,ia+step,:,:,:)
       step=step+1
     enddo
     ! Step 5: Calculation of the dynamical correlation function in the rotating frame
     ! Definition of the energy coefficients
     a=emin-((emax-emin)/(real(energy_step)-1))
     b=(emax-emin)/(real(energy_step)-1)
     ! eigenvalues for positive q vectors
     eig_val_2k_redplus(:,:)=eig_val_2k(:,:)
     ! Calculating the eigenvalues for negative q vectors
     counter_dcf=counter_dcf+1
     do i=1, nq
      q(:,i)=-q(:,i)
     enddo
        ! call  setup_finite_hamiltonian_full(N1,N2,N3,NA,Natom, Mensemble, simid, emomM, mmom)
     ! -q vectors
     eig_val_2k_redminus(:,:)=eig_val_2k(:,:)
     ! Turning q-vectors back to original values
     do i=1, nq
       q(:,i)=-q(:,i)
     enddo
     ! Definition of the dynamical correlation function in the rotating frame S_ab
     do i=1,nq_ext
       do alpha=1,3
         do beta=1,3
           do j=1,energy_step
             ! Phason mode (i or k)
             if (i .le. nq) then
               do ia=1,NA
                  S_ab(i,alpha,beta,j)=S_ab(i,alpha,beta,j)+(1.0d0/(real(NA)))*Nab_mat_red(ia,ia,alpha,beta,i)*&
                  exp(-((a+real(j)*b)-eig_val_2k_redplus(ia,i))**2.0d0/(2.0d0*sigma_corr**2.0d0))*&
                  (1.0d0/(exp(eig_val_2k_redplus(ia,i)/(k_bolt_ev*1000.0d0*temp_corr))-1))+&
                  (1.0d0/(real(NA)))*Nab_mat_red(ia,ia,alpha,beta,i)*&
                  exp(-((a+real(j)*b)+eig_val_2k_redminus(ia,i))**2.0d0/(2.0d0*sigma_corr**2.0d0))*&
                  (1.0d0/(exp(eig_val_2k_redminus(ia,i)/(k_bolt_ev*1000.0d0*temp_corr))-1)+1)
               end do
             ! i+q_vector mode (k+Q)
             else if (i .gt. nq .and. i .le. 2*nq) then
               l=i-nq
               do ia=1,NA
                  S_ab_plus(l,alpha,beta,j)=S_ab_plus(l,alpha,beta,j)+(1.0d0/(real(NA)))*Nab_mat_red(ia,ia,alpha,beta,i)*&
                  exp(-((a+real(j)*b)-eig_val_2k_redplus(ia,i))**2.0d0/(2.0d0*sigma_corr**2.0d0))*&
                  (1.0d0/(exp(eig_val_2k_redplus(ia,i)/(k_bolt_ev*1000.0d0*temp_corr))-1))+&
                  (1.0d0/(real(NA)))*Nab_mat_red(ia,ia,alpha,beta,i)*&
                  exp(-((a+real(j)*b)+eig_val_2k_redminus(ia,i))**2.0d0/(2.0d0*sigma_corr**2.0d0))*&
                  (1.0d0/(exp(eig_val_2k_redminus(ia,i)/(k_bolt_ev*1000.0d0*temp_corr))-1)+1)
               end do
             ! i-q_vector mode (k-Q)
             else
               l=i-(2*nq)
               do ia=1,NA
                  S_ab_minus(l,alpha,beta,j)=S_ab_minus(l,alpha,beta,j)+(1.0d0/(real(NA)))*Nab_mat_red(ia,ia,alpha,beta,i)*&
                  exp(-((a+real(j)*b)-eig_val_2k_redplus(ia,i))**2.0d0/(2.0d0*sigma_corr**2.0d0))*&
                  (1.0d0/(exp(eig_val_2k_redplus(ia,i)/(k_bolt_ev*1000.0d0*temp_corr))-1))+&
                  (1.0d0/(real(NA)))*Nab_mat_red(ia,ia,alpha,beta,i)*&
                  exp(-((a+real(j)*b)+eig_val_2k_redminus(ia,i))**2.0d0/(2.0d0*sigma_corr**2.0d0))*&
                  (1.0d0/(exp(eig_val_2k_redminus(ia,i)/(k_bolt_ev*1000.0d0*temp_corr))-1)+1)
               end do
             end if
           end do
         end do
       end do
     end do
     ! Identity matrix
     I_mat=(0.0d0,0.0d0)
     I_mat(1,1)=(1.0d0,0.0d0)
     I_mat(2,2)=(1.0d0,0.0d0)
     I_mat(3,3)=(1.0d0,0.0d0)
     ! Definition of matrix R2
     R2_mat=cmplx(R_mat(1,qvector,2),0.0d0)
      !R2_mat=(R_mat(1,qvector,2),0.0d0)
     ! Definition of K matrix
     K_mat=cmplx(0.0d0,R_mat(1,qvector,1))
     !K_mat=(0.0d0,R_mat(1,qvector,1))
     ! Definition of matrix R1
     R1_mat=(0.5d0,0.5d0)*(I_mat-i_complex*K_mat-R2_mat)
     ! Defining the complex conjugate of matrix R1
     R1conj_mat=(0.5d0,0.5d0)*(I_mat+i_complex*K_mat-R2_mat)
     ! Definition of the dynamical correlation function in the laboratory frame
    
    
     do i=1,nq
       do j=1,energy_step
         !write(*,*) j
         S_mat(i,:,:,j)=abs(matmul(S_ab(i,:,:,j),R2_mat(:,:))+matmul(S_ab_plus(i,:,:,j),R1_mat(:,:))+matmul(S_ab_minus(i,:,:,j),R1conj_mat(:,:)))   
       end do
     end do
 
     ! Normalising the intensity to 1
     if (do_norm_corr == 'Y') then
       do alpha=1,3
         do beta=1,3
           !Real part of S_mat
           if (maxval(abs(real(S_mat(:,alpha,beta,:)))) /= 0.0d0) then
             S_mat(:,alpha,beta,:)=cmplx(abs(real(S_mat(:,alpha,beta,:)))/maxval(abs(real(S_mat(:,alpha,beta,:)))),aimag(S_mat(:,alpha,beta,:)))
           end if
           !Imaginary part of S_mat
           if (maxval(abs(aimag(S_mat(:,alpha,beta,:)))) /= 0.0d0) then
             S_mat(:,alpha,beta,:)=cmplx(real(S_mat(:,alpha,beta,:)),abs(aimag(S_mat(:,alpha,beta,:)))/maxval(abs(aimag(S_mat(:,alpha,beta,:)))))
           end if
         end do
      end do
     end if
   
     
     deallocate(S_ab,stat=i_stat)
     call memocc(i_stat,-product(shape(S_ab))*kind(S_ab),'S_ab','calc_dcf')
     deallocate(S_ab_plus,stat=i_stat)
     call memocc(i_stat,-product(shape(S_ab_plus))*kind(S_ab_plus),'S_ab_plus','calc_dcf')
     deallocate(S_ab_minus,stat=i_stat)
     call memocc(i_stat,-product(shape(S_ab_minus))*kind(S_ab_minus),'S_ab_minus','calc_dcf')
     deallocate(Y_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(Y_mat))*kind(Y_mat),'Y_mat','calc_dcf')
     deallocate(Z_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(Z_mat))*kind(Z_mat),'Z_mat','calc_dcf')
     deallocate(V_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(V_mat))*kind(V_mat),'V_mat','calc_dcf')
     deallocate(O_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(O_mat))*kind(O_mat),'O_mat','calc_dcf')
     deallocate(M_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(M_mat))*kind(M_mat),'M_mat','calc_dcf')
     deallocate(Nab_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(Nab_mat))*kind(Nab_mat),'Nab_mat','calc_dcf')
     deallocate(Nab_mat_red,stat=i_stat)
     call memocc(i_stat,-product(shape(Nab_mat_red))*kind(Nab_mat_red),'Nab_mat_red','calc_dcf')
     !deallocate(T_mat,stat=i_stat)
     !call memocc(i_stat,-product(shape(T_mat))*kind(T_mat),'T_mat','calc_dcf')
     deallocate(N_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(N_mat))*kind(N_mat),'N_mat','calc_dcf')
     deallocate(Mdum1_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(Mdum1_mat))*kind(Mdum1_mat),'Mdum1_mat','calc_dcf')
     deallocate(Tdum_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(Tdum_mat))*kind(Tdum_mat),'Tdum_mat','calc_dcf')
     deallocate(dum1_mat,stat=i_stat)
     call memocc(i_stat,-product(shape(dum1_mat))*kind(dum1_mat),'dum1_mat','calc_dcf')
     deallocate(eig_val_2k_redplus,stat=i_stat)
     call memocc(i_stat,-product(shape(eig_val_2k_redplus))*kind(eig_val_2k_redplus),'eig_val_2k_redplus','calc_dcf')
     deallocate(eig_val_2k_redminus,stat=i_stat)
     call memocc(i_stat,-product(shape(eig_val_2k_redminus))*kind(eig_val_2k_redminus),'eig_val_2k_redminus','calc_dcf')
    ! deallocate(S_mat,stat=i_stat)
    ! call memocc(i_stat,-product(shape(S_mat))*kind(S_mat),'S_mat',' setup_finite_hamiltonian_full')
  

 end subroutine calc_dcf

 subroutine find_uv(u,v,mom)
   !
   implicit none
   !Global variables
   real(dblprec), dimension(3), intent(in) :: mom
   complex(dblprec), dimension(3), intent(out) :: u
   complex(dblprec), dimension(3), intent(out) :: v
   !Local variables
   real(dblprec) :: mnorm, theta
   real(dblprec), dimension(3,3) :: R_prime, Ku_mat,I_mat
   real(dblprec), dimension(3) :: k_vec, mom_hat
   complex(dblprec) :: im = (0.0d0,1.0d0)
   !
   !
   mnorm=sqrt(sum(mom*mom)+1.0d-15)
   ! Find the rotation matrix R_prime that transforms mom to z-axis.
   ! Done by finding the vector orthogonal to mom and z and then perform a
   ! Rodrigues rotation with the corresponding angle
   !
   ! Unit vector parallel to mom
   mom_hat=mom/mnorm
   ! Angle between mom and z (rotated anticlockwise)
   theta=-acos(mom_hat(3))
   ! k_vec= cross(mom,z) : unit vector defining the rotation axis
   mom_hat=mom_hat+1.0d-15
   k_vec(1)=mom_hat(2)
   k_vec(2)=-mom_hat(1)
   k_vec(3)=0.0d0
   k_vec=k_vec/sqrt(sum(k_vec*k_vec))
   ! cross product matrix (k_vec x mom) also denoted as Ku matrix
   Ku_mat=0.0d0
   Ku_mat(1,2)= -k_vec(3)
   Ku_mat(1,3)=  k_vec(2)
   Ku_mat(2,1)=  k_vec(3)
   Ku_mat(2,3)= -k_vec(1)
   Ku_mat(3,1)= -k_vec(2)
   Ku_mat(3,2)=  k_vec(1)
   ! Identity matrix
   I_mat=0.0d0
   I_mat(1,1)=1.0d0
   I_mat(2,2)=1.0d0
   I_mat(3,3)=1.0d0
   ! Rodrigues rotation matrix
   R_prime=I_mat+sin(theta)*Ku_mat+(1.0d0-cos(theta))*matmul(Ku_mat,Ku_mat)
   ! Definition of u and v vectors
   u=R_prime(:,1)+im*R_prime(:,2)
   v=R_prime(:,3)
   !
 end subroutine find_uv

 function jncoup_q(mu,nu,q_vect,qvector,countstart)
   !
   implicit none
   !Global variables
   integer, intent(in) :: mu
   integer, intent(in) :: nu
   integer, intent(in) :: countstart
   real(dblprec), dimension(:), intent(in) :: q_vect
   real(dblprec), dimension(3), intent(in) :: qvector
   !Local variables
   integer :: j,k, mutemp, nutemp,ii,jj,jat,kk
   real(dblprec), dimension(3) :: q_vect2pi,q_ani
   real(dblprec), dimension(3) :: dist
   real(dblprec), dimension(3,3) :: J_mat, K_mat,Jleftprod_mat, Jrightprod_mat,Kleftprod_mat, Krightprod_mat
   complex(dblprec) :: i
   complex(dblprec), dimension(3,3)  :: jncoup_q,test
   !Initialising variables
   jncoup_q=0.0_dblprec
   q_vect2pi=q_vect*2_dblprec*pi
   mutemp = mu+countstart
   nutemp = nu+countstart
   i = (0.0d0,1.0d0)
   J_mat=0.0d0
   Jleftprod_mat=0.0d0
   Jrightprod_mat=0.0d0
   Kleftprod_mat=0.0d0
   Krightprod_mat=0.0d0
   q_ani(1)=0.0d0
   q_ani(2)=0.0d0
   q_ani(3)=0.0d0
   !Computing Fourier transform of isotropic J
   
   
   kk=0
   do j=1,ham%nlistsize(mutemp)
            
           if (anumb(nutemp)==anumb(ham%nlist(j,mutemp))) then
                dist(:)=coord(1:3,mutemp)-coord(1:3,ham%nlist(j,mutemp))
                do k=1,3
                  J_mat(k,k) = ham%ncoup(j,mutemp,1) 
                end do
          
                 
               if(ham_inp%do_dm==1) then
                   J_mat(1,2) =  -ham%dm_vect(3,j,mutemp)
                   J_mat(1,3) =   ham%dm_vect(2,j,mutemp)
                   J_mat(2,1) =   ham%dm_vect(3,j,mutemp)
                   J_mat(2,3) =  -ham%dm_vect(1,j,mutemp)
                   J_mat(3,1) =  -ham%dm_vect(2,j,mutemp)
                   J_mat(3,2) =   ham%dm_vect(1,j,mutemp)
                endif

                if (ham_inp%do_pd==1) then
                         do ii=1,3
                                  do jj=1,3
                                   if(ii==jj) then
                                      J_mat(ii,jj)=J_mat(ii,jj)+ham%pd_vect(ii,j,mutemp)
                                    else
                                     J_mat(ii,jj)=J_mat(ii,jj)+ham%pd_vect( (3*(ii-1)+jj) ,j,mutemp)
                                    end if
                                  end do
                        end do
                      
                end if
            
                   Jleftprod_mat=matmul(transpose(R_mat(mutemp,qvector,0)),J_mat)
                   Jrightprod_mat=matmul(Jleftprod_mat,R_mat(ham%nlist(j,mutemp),qvector,0))
                   jncoup_q = jncoup_q +1.0d0*Jrightprod_mat*exp(-i*&
                   (q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3)))
           end if
   end do
   
                if ((nutemp==mutemp)  ) then  
                    K_mat=0.0d0
                     if((ham_inp%do_anisotropy==1)  ) then
                               do ii=1,3
                                  do jj=1,3
                                   if(ii==jj) then
                                     K_mat(ii,jj)  = K_mat(ii,jj)  -2.0d0*ham%kaniso(1,mutemp)*(ham%eaniso(ii,mutemp) *ham%eaniso(ii,mutemp) )
                                   else
                                     K_mat(ii,jj)  = K_mat(ii,jj) -2.0d0*ham%kaniso(1,mutemp)*(ham%eaniso(ii,mutemp) *ham%eaniso(jj,mutemp)) 
                                     
                                   endif
                                 enddo
                            enddo
                        if(ham_inp%mult_axis=='Y') then
                              
                                 do ii=1,3
                                       do jj=1,3
                                          if(ii==jj) then
                                          K_mat(ii,jj)  = K_mat(ii,jj)-2.0d0*ham%kaniso_diff(1,mutemp)*(ham%eaniso_diff(ii,mutemp)*ham%eaniso_diff(ii,mutemp))
                                          else
                                         K_mat(ii,jj)  = K_mat(ii,jj) -2.0d0* ham%kaniso_diff(1,mutemp)*(ham%eaniso_diff(ii,mutemp)*ham%eaniso_diff(jj,mutemp))
                                        
                                       endif
                                     enddo
                                 enddo
                           endif
                      ! Here I add the anisotrpy. Actually we do not need the operatios, but I left them for consistency.                    
                        Kleftprod_mat=matmul(transpose(R_mat(mutemp,q_vector,0)),K_mat)  
                        Krightprod_mat=matmul(Kleftprod_mat,R_mat(mutemp,q_vector,0))
                        jncoup_q =  jncoup_q+ Krightprod_mat
                     endif
                endif   
              


   ! End of Fourier transform of J function
 end function jncoup_q

 function jtens_q(mu,nu,q_vect,qvector,countstart)
   !
   implicit none
   !Global variables
   integer, intent(in) :: mu
   integer, intent(in) :: nu
   integer, intent(in) :: countstart
   real(dblprec), dimension(:),intent(in) :: q_vect
   real(dblprec), dimension(3), intent(in) :: qvector
   !Local variables
   integer ::  mutemp, nutemp, j,ii,jj
   real(dblprec), dimension(3) :: q_vect2pi
   real(dblprec), dimension(3) :: dist
   real(dblprec), dimension(3,3) :: Jtensleft_mat, Jtensright_mat
   real(dblprec), dimension(3,3) ::  K_mat,Kleftprod_mat, Krightprod_mat
   complex(dblprec) :: i
   complex(dblprec), dimension(3,3) :: jtens_q
   !Initialising variables
   jtens_q=0.0_dblprec
   q_vect2pi=q_vect*2_dblprec*pi
   mutemp = mu+countstart
   nutemp = nu+countstart
   i = (0.0d0,1.0d0)
   Jtensleft_mat=0.0d0
   Jtensright_mat=0.0d0
    K_mat=0.0d0
    Kleftprod_mat=0.0d0
     Krightprod_mat=0.0d0
   ! Computing Fourier transform of J tensor
   do j=1,ham%nlistsize(mutemp)
           if (anumb(nutemp)==anumb(ham%nlist(j,mutemp))) then
                dist(:)=coord(1:3,mutemp)-coord(1:3,ham%nlist(j,mutemp))
              
                
                Jtensleft_mat=matmul(transpose(R_mat(mutemp,qvector,0)),ham%j_tens(:,:,j,mutemp))
                Jtensright_mat=matmul(Jtensleft_mat,R_mat(ham%nlist(j,mutemp),qvector,0))
                jtens_q = jtens_q + Jtensright_mat*exp(-i* &
                (q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3)))
           end if
   end do
         
    if ((nutemp==mutemp) ) then  
                   
                     if((ham_inp%do_anisotropy==1)  ) then
                               do ii=1,3
                                  do jj=1,3
                                   if(ii==jj) then
                                     K_mat(ii,jj)  = K_mat(ii,jj)  -2.0d0*ham%kaniso(1,mutemp)*(ham%eaniso(ii,mutemp) *ham%eaniso(ii,mutemp) )
                                   else
                                    K_mat(ii,jj)  =K_mat(ii,jj) -2.0d0*ham%kaniso(1,mutemp)*(ham%eaniso(ii,mutemp) *ham%eaniso(jj,mutemp)) 
                                     
                                   endif
                                 enddo
                            enddo
                        if(ham_inp%mult_axis=='Y') then
                              
                                 do ii=1,3
                                       do jj=1,3
                                          if(ii==jj) then
                                          K_mat(ii,jj)  = K_mat(ii,jj)-2.0d0*ham%kaniso_diff(1,mutemp)*(ham%eaniso_diff(ii,mutemp)*ham%eaniso_diff(ii,mutemp))
                                          else
                                         K_mat(ii,jj)  = K_mat(ii,jj) -2.0d0* ham%kaniso_diff(1,mutemp)*(ham%eaniso_diff(ii,mutemp)*ham%eaniso_diff(jj,mutemp))
                                        
                                       endif
                                     enddo
                                 enddo
                           endif
                           
                             
                        Kleftprod_mat=matmul(transpose(R_mat(mutemp,q_vector,0)),K_mat)
                        Krightprod_mat=matmul(Kleftprod_mat,R_mat(mutemp,q_vector,0))
                        jtens_q=  jtens_q+ Krightprod_mat
                     endif
                endif   
       
   ! End of Fourier transform of tensor J function
 end function jtens_q

 function R_mat(atom_num,qvector,flag)
   !
   implicit none
   !Global variables
   integer, intent(in) :: atom_num
   integer, intent(in) :: flag
   real(dblprec), dimension(3), intent(in) :: qvector
   real(dblprec), dimension(3,3) :: R_mat
   !Local variables
   real(dblprec) :: theta
   real(dblprec), dimension(3,3) :: Ku_mat,I_mat
   ! Definition of the rotation axis
   if (q_vector(1)==0.0d0 .and. q_vector(2)==0.0d0 .and. q_vector(3)==0.0d0) then
      rot_axis(:) = emom(:,atom_num,1)
   end if
   ! Find the rotation matrix R performing a
   ! Rodrigues rotation with the corresponding angle (q_vector*R)
   ! Rotation angle
   theta=(qvector(1)*coord(1,atom_num)+qvector(2)*coord(2,atom_num)+qvector(3)*coord(3,atom_num))*pi
   ! Matrix defined by the rotation axis
   Ku_mat=0.0d0
   Ku_mat(1,2)=-rot_axis(3)
   Ku_mat(1,3)= rot_axis(2)
   Ku_mat(2,1)= rot_axis(3)
   Ku_mat(2,3)=-rot_axis(1)
   Ku_mat(3,1)=-rot_axis(2)
   Ku_mat(3,2)= rot_axis(1)
   ! Identity matrix
   I_mat=0.0d0
   I_mat(1,1)=1.0d0
   I_mat(2,2)=1.0d0
   I_mat(3,3)=1.0d0
   if (flag == 0) then
   ! Rodrigues rotation matrix
   R_mat=I_mat+sin(theta)*Ku_mat+(1.0d0-cos(theta))*matmul(Ku_mat,Ku_mat)
   else if (flag == 1) then
   ! Axis of rotation for rotation matrix R1
   R_mat=Ku_mat
   else if (flag == 2) then
   ! Rotation matrix R2
   R_mat(1,1)=rot_axis(1)**2
   R_mat(1,2)=rot_axis(1)*rot_axis(2)
   R_mat(1,3)=rot_axis(1)*rot_axis(3)
   R_mat(2,1)=rot_axis(2)*rot_axis(1)
   R_mat(2,2)=rot_axis(2)**2
   R_mat(2,3)=rot_axis(2)*rot_axis(3)
   R_mat(3,1)=rot_axis(3)*rot_axis(1)
   R_mat(3,2)=rot_axis(3)*rot_axis(2)
   R_mat(3,3)=rot_axis(3)**2
   end if
   !
 end function R_mat

 complex(dblprec) function sJs(u,J,v)
   !
   implicit none
   !
   complex(dblprec), dimension(3), intent(in) :: u
   complex(dblprec), dimension(3,3), intent(in) :: J
   complex(dblprec), dimension(3), intent(in) :: v
   !
   complex(dblprec), dimension(3) :: dt
   !
   dt(1)=J(1,1)*v(1)+J(1,2)*v(2)+J(1,3)*v(3)
   dt(2)=J(2,1)*v(1)+J(2,2)*v(2)+J(2,3)*v(3)
   dt(3)=J(3,1)*v(1)+J(3,2)*v(2)+J(3,3)*v(3)
   !
   sJs=u(1)*dt(1)+u(2)*dt(2)+u(3)*dt(3)
   !
 end function sJs



 subroutine printeigenvalues(filename,file_id,wres,msat,flag,flag1)
  !
  implicit none
  !Global variables
  integer, intent(in) :: file_id, flag, flag1
  character(LEN=*),intent(in) :: filename
  real(dblprec), intent(in) :: msat
  real(dblprec), dimension(Na,nq_ext), intent(in) :: wres
  !Local variables
  integer :: i, k, r, ia, eigv, counter
  real(dblprec) :: tcmfa,tcrpa
  !Open file
  open(file_id,file=filename)
  ! Definition of the maximun value of wres and truncated to the integer
  eigv= int(aint(maxval(abs(real(wres)))))
  if (eigv == 0) then
     eigv = 1
  endif
  ! number of digits
  counter = 0
  do while (eigv .ne. 0)
  eigv=eigv/10
  counter=counter+1
  end do
  ! Redefinition of the counters (keep in mind the number of decimal points in the format 12 decimal points + 1 point=13)
  counter = counter + 13
  !
  if(flag==1) then
    tcmfa=0.d0 ; tcrpa=0.d0
  endif
  k=0 ; r=0
  ! Only printed half of the eigenvalues since that the other half are negative or repeated.
  if(flag1==0) then
    write(file_id,1003) 'q-point','q-point list', 'mode'
    do i =1,nq
         write(file_id,1001) qnumber(i), i, wres(:,i)
         if (flag==1) then
            do ia=1,na
              if (wres(ia,i)>=1.0d-5) then
                  tcrpa=tcrpa+(1.d0/wres(ia,i))
                  r=r+1
              endif
            end do
            tcmfa=tcmfa+sum(wres(:,i))
            k=k+na
         endif
    end do
  else    
      write(file_id,1004) 'q-point','q-point list','Phason mode', 'K+Q mode', 'K-Q mode'    
    do i =1,nq
         write(file_id,1005)  qnumber(i), i, wres(:,i), wres(:,i+nq), wres(:,i+2*nq)
         if (flag==1) then
            do ia=1,na
              if (wres(ia,i)>=1.0d-5) then
                  tcrpa=tcrpa+(1.d0/wres(ia,i))
                  r=r+1
              endif
            end do
            tcmfa=tcmfa+sum(wres(:,i))
            k=k+na
         endif
    end do
  end if
  if (flag==1) then
     tcmfa=(msat*tcmfa)/(6*k*k_bolt_ev*1000)
     write(*,1002) 'Tc-MFA from non-collinear LSWT:' , tcmfa
     tcrpa=((r*msat)/(6*k_bolt_ev*1000))*(1.d0/tcrpa)
     write(*,1002) 'Tc-RPA from non-collinear LSWT:' , tcrpa
  endif
  !

  1001 format (2x,f11.6,2x,i5,2x,9es16.8)
  1002 format (2x,a,f10.1)
  1003 format (2x,a,2x,a,2x,a)
  1004 format (2x,a,2x,a,2x,a,2x,a,2x,a)
  1005 format (2x,f11.6,2x,i5,2x,9es16.8,9es16.8,9es16.8)
  close(file_id)
  !
  end subroutine printeigenvalues

  subroutine printeigenvectors(filename,file_id,eigv,wres,q,nq,na,flag)
   !
   implicit none
   ! Global variables
   integer,intent(in) :: nq
   integer,intent(in) :: na
   integer,intent(in) :: file_id,flag
   real(dblprec),dimension(2*na,nq_ext), intent(in) :: wres
   real(dblprec),dimension(3,nq_ext), intent(in) :: q
   complex(dblprec),dimension(2*na,2*na,nq_ext), intent(in) :: eigv
   character(LEN=*),intent(in) :: filename
   ! Local variables
   integer :: i,j,k,real_eigv,img_eigv, counter_real, counter_img
   character(LEN=2), dimension(2*na,2*na,nq_ext) :: imag_unit
   ! Open file
   open(file_id,file=filename)
   !Definition of variables
   imag_unit = '+i'
   !where(aimag(eigv)<0.0d0) imag_unit = '-i'
   ! Definition of the maximun real and imaginary values of eigv and truncated to the integer
  ! real_eigv= int(aint(maxval(real(eigv))))
  ! img_eigv= int(aint(maxval(aimag(eigv))))
  ! if (real_eigv == 0) then
  !    real_eigv = 1
  ! endif
  ! if (img_eigv == 0) then
  !    img_eigv = 1
  ! endif
   ! number of digits real part
  ! do while (real_eigv .ne. 0)
  ! counter_real = 0
  ! real_eigv=real_eigv/10
  ! counter_real=counter_real+1
  ! end do
   ! number of digits imaginary part
  ! counter_img = 0
  ! do while ( img_eigv .ne. 0)
  ! img_eigv=img_eigv/10
  ! counter_img=counter_img+1
 !  end do
   ! Redefinition of the counters (keep in mind the number of decimal points in the format 3 decimal points + 1 point+ sign=5)
   counter_real =  6
   counter_img = 6
   ! Write file
   if (flag == 0) then
     do i=1,nq
        write(file_id,1001) "# q-point index ",i, " vector: ", q(1:3,i)
        do j=1,2*na
             do k=1,(2*NA)
           ! Compact format : Value only for each complex number
              write(file_id,1002,advance='NO') wres(j,i), real(eigv(k,j,i)),imag_unit(k,j,i),aimag(eigv(k,j,i)) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             end do
        end do
         write(file_id,*)
     end do
   else
     do i=1,nq_ext
        if (i .le. nq) then
        write(file_id,1001) "(Phason mode) q-point index ",i, " vector: ", q(1:3,i)
        else if (i .gt. nq .and. i .le. nq*2 ) then
        write(file_id,1001) "(K+Q mode) q-point index ",i-nq, " vector: ", q(1:3,i)
        else
        write(file_id,1001) "(K-Q mode) q-point index ",i-2*nq, " vector: ", q(1:3,i)
        end if
        do j=1,2*na
           ! Compact format : Value only for each complex number
             do k=1,(2*NA)
              write(file_id,1002,advance='NO') wres(j,i), real(eigv(k,j,i)),imag_unit(k,j,i),aimag(eigv(k,j,i))   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             end do
        end do
         write(file_id,*) 
     end do
   end if
   ! Close file
   close(file_id)
   ! Format file
   1001 format (a,2x,i5,a,3f12.6)
   1002 format (1x,f9.3,2x,f9.3,2x,a2,f9.3,2x)   
   ! end subroutine
 end subroutine printeigenvectors

 subroutine lswtdos_calc(filename,file_id,wres,flag)
  !
  implicit none
  !Global variables
  integer,intent(in) :: file_id,flag
  real(dblprec), dimension(na,nq_ext), intent(in) :: wres
  character(LEN=*),intent(in) :: filename
  !Local variables
  integer :: i, k, ia,i_stat
  real(dblprec), allocatable,dimension(:,:) :: lswtdos, lswtdos1, lswtdos2, lswtdos3
  real(dblprec) :: emin, emin1, emin2, emin3, emax, emax1, emax2, emax3, deltae, deltae1, deltae2, deltae3, fac, fac1, fac2, fac3
  !
    open(file_id,file=filename)
    !
    if (flag==0) then
      emin=minval(wres)
      emax=maxval(wres)
      deltae=(emax-emin)/(magdos_freq-1)
      allocate(lswtdos(magdos_freq,2),stat=i_stat)
      call memocc(i_stat,product(shape(lswtdos))*kind(lswtdos),'lswtdos','lswtdos_calc')
      lswtdos=0.0d0
      ! call sortRealArray(wres(:,i),NA)
      do i=1,magdos_freq
          lswtdos(i,1)=emin+(i-1)*deltae
      enddo
      !
      do k=1,magdos_freq
         do i =1,nq
             do ia=1,na
                fac=exp(-(wres(ia,i)-lswtdos(k,1))**2/magdos_sigma)
                lswtdos(k,2)=lswtdos(k,2)+fac
             enddo
         enddo
      enddo
      !Normalize DOS
      lswtdos(:,2)=lswtdos(:,2)/(sum(lswtdos(:,2))*deltae)
      !
      write(file_id,*) 'Phason mode'
      do k=1,magdos_freq
        write(file_id,1001) lswtdos(k,1),lswtdos(k,2),sum(lswtdos(1:k,2))*deltae
      enddo
    else
      emin1=minval(wres(:,1:nq))
      emax1=maxval(wres(:,1:nq))
      deltae1=(emax1-emin1)/(magdos_freq-1)
      emin2=minval(wres(:,nq+1:nq*2))
      emax2=maxval(wres(:,nq+1:nq*2))
      deltae2=(emax2-emin2)/(magdos_freq-1)
      emin3=minval(wres(:,nq*2+1:nq*3))
      emax3=maxval(wres(:,nq*2+1:nq*3))
      deltae3=(emax3-emin3)/(magdos_freq-1)
      allocate(lswtdos1(magdos_freq,2),stat=i_stat)
      call memocc(i_stat,product(shape(lswtdos1))*kind(lswtdos1),'lswtdos1','lswtdos_calc')
      allocate(lswtdos2(magdos_freq,2),stat=i_stat)
      call memocc(i_stat,product(shape(lswtdos2))*kind(lswtdos2),'lswtdos2','lswtdos_calc')
      allocate(lswtdos3(magdos_freq,2),stat=i_stat)
      call memocc(i_stat,product(shape(lswtdos3))*kind(lswtdos3),'lswtdos3','lswtdos_calc')
      lswtdos=0.0d0
      ! call sortRealArray(wres(:,i),NA)
      do i=1,magdos_freq
          lswtdos1(i,1)=emin1+(i-1)*deltae1
          lswtdos2(i,1)=emin2+(i-1)*deltae2
          lswtdos3(i,1)=emin3+(i-1)*deltae3
      enddo
      !
      do k=1,magdos_freq
         do i=1,nq_ext
             do ia=1,na
               if (i .le. nq) then
                 fac1=exp(-(wres(ia,i)-lswtdos1(k,1))**2/magdos_sigma)
                 lswtdos1(k,2)=lswtdos1(k,2)+fac1
               else if (i .gt. nq .and. i .le. nq*2) then
                 fac2=exp(-(wres(ia,i)-lswtdos2(k,1))**2/magdos_sigma)
                 lswtdos2(k,2)=lswtdos2(k,2)+fac2
               else
                 fac3=exp(-(wres(ia,i)-lswtdos3(k,1))**2/magdos_sigma)
                 lswtdos3(k,2)=lswtdos3(k,2)+fac3
               end if
             enddo
         enddo
      enddo
      !Normalize DOS
      lswtdos1(:,2)=lswtdos1(:,2)/(sum(lswtdos1(:,2))*deltae1)
      lswtdos2(:,2)=lswtdos2(:,2)/(sum(lswtdos2(:,2))*deltae2)
      lswtdos3(:,2)=lswtdos3(:,2)/(sum(lswtdos3(:,2))*deltae3)
      !
      write(file_id,*) 'Phason mode'
      do k=1,magdos_freq
        write(file_id,1001) lswtdos1(k,1),lswtdos1(k,2),sum(lswtdos1(1:k,2))*deltae1
      enddo
      write(file_id,*) 'K+Q mode'
      do k=1,magdos_freq
        write(file_id,1001) lswtdos2(k,1),lswtdos2(k,2),sum(lswtdos2(1:k,2))*deltae2
      enddo
      write(file_id,*) 'K-Q mode'
      do k=1,magdos_freq
        write(file_id,1001) lswtdos3(k,1),lswtdos3(k,2),sum(lswtdos3(1:k,2))*deltae3
      enddo
    end if
    ! deallocate arrays
    if (flag == 0) then
      deallocate(lswtdos)
    else
      deallocate(lswtdos1)
      deallocate(lswtdos2)
      deallocate(lswtdos3)
    end if
  ! Format of variables
  1001 format (2x,f18.12,2f18.12)
  1002 format (2x,a,f10.1)
  ! Close file
  close(file_id)
  !
 end subroutine lswtdos_calc

  subroutine printjq(filename,file_id,jq_zero,jq,flag)
  !
  implicit none
  ! Global variables
  integer,intent(in) :: file_id,flag
  real(dblprec), dimension(na,nq_ext), intent(in) :: jq_zero
  real(dblprec), dimension(na,nq_ext), intent(in) :: jq
  character(LEN=*),intent(in) :: filename
  ! Local variables
  integer :: i, eigv, eigv1, eigv2, eigv3, counter, counter1, counter2, counter3
  ! Open file
  open(file_id,file=filename)
  if (flag==0) then
    ! Definition of the maximun value of jq and truncated to the integer
    eigv= int(aint(maxval(abs(real(jq)))))
    if (eigv == 0) then
       eigv = 1
    endif
    ! number of digits
    counter = 0
    do while (eigv .ne. 0)
    eigv=eigv/10
    counter=counter+1
    end do
    ! Redefinition of the counters (keep in mind the number of decimal points in the format 12 decimal points + 1 point + sign=14)
    counter = counter + 14
    ! Print j(q)
    write(file_id,1002) 'q-point', 'J(q)', 'J(0)'
    do i =1,nq
         write(file_id,1004) i, jq(:,i),jq_zero(:,i)
    end do
  else
    ! Definition of the maximun value of jq and truncated to the integer
    eigv1= int(aint(maxval(abs(real(jq(:,1:nq))))))
    eigv2= int(aint(maxval(abs(real(jq(:,nq+1:2*nq))))))
    eigv3= int(aint(maxval(abs(real(jq(:,2*nq+1:nq*3))))))
    if (eigv1 == 0) then
       eigv1 = 1
    endif
    if (eigv2 == 0) then
       eigv2 = 1
    endif
    if (eigv3 == 0) then
       eigv3 = 1
    endif
    ! number of digits
    counter1 = 0
    do while (eigv1 .ne. 0)
    eigv1=eigv1/10
    counter1=counter1+1
    end do
    counter2 = 0
    do while (eigv2 .ne. 0)
    eigv2=eigv2/10
    counter2=counter2+1
    end do
    counter3 = 0
    do while (eigv3 .ne. 0)
    eigv3=eigv3/10
    counter3=counter3+1
    end do
    ! Redefinition of the counters (keep in mind the number of decimal points in the format 12 decimal points + 1 point + sign=14)
    counter1 = counter1 + 14
    counter2 = counter2 + 14
    counter3 = counter3 + 14
    ! Print j(q)
    write(file_id,1003) 'q-point', 'J(K)', 'J(0)', 'J(K+Q)', 'J(K-Q)'
    do i =1,nq
         write(file_id,1001) i, jq(:,i),jq_zero(:,i), jq(:,i+nq),jq(:,i+2*nq)
    end do
  end if
  !
  1001 format (2x,i5,2x,9es16.8,2x,9es16.8,2x,9es16.8,2x,9es16.8,2x)
  1002 format (2x,a,2x,a,2x,a,2x)
  1003 format (2x,a,5x,a,12x,a,12x,a,12x,a,12x)
  1004 format (2x,i5,2x,9es16.8,2x,9es16.8,2x)
  close(file_id)
  !
  end subroutine printjq

 subroutine print_dcf(filename1,file_id1,filename2,file_id2,emax,emin,S_mat)
   !
   implicit none
   ! Global variables
   integer, intent(in) :: file_id1, file_id2
   real(dblprec), intent(in) :: emax, emin
   complex(dblprec), dimension(nq,3,3,energy_step), intent(in) ::  S_mat
   
   character(LEN=*), intent(in) :: filename1, filename2
   ! Local variables
   integer :: i, j, ii,jj,real_S, img_S, counter_real, counter1_real, counter_img, counter1_img,energy, &
              counter_energy, qpoint, counter_qpoint, alpha, beta
   real :: a,b
   real(dblprec)::ssmr, ssmi,maxr,maxi

   ! Open file
   open(file_id1,file=filename1)
   open(file_id2,file=filename2)
   ! Coeficients of the energy
   a=emin-((emax-emin)/(real(energy_step)-1))
   b=(emax-emin)/(real(energy_step)-1)
   ! Definition of the maximun real and imaginary values of S and truncated to the integer
   real_S = maxval(abs(exponent(real(S_mat))))
   img_S = maxval(abs(exponent(aimag(S_mat))))
   if (real_S == 0) then
      real_S = 1
   endif
   if (img_S == 0) then
      img_S = 1
   endif
   ! number of digits real part
   counter_real = 0
   do while (real_S .ne. 0)
   real_S=real_S/10
   counter_real=counter_real+1
   end do
   ! number of digits imaginary part
   counter_img = 0
   do while ( img_S .ne. 0)
   img_S=img_S/10
   counter_img=counter_img+1
   end do
   ! Counter for the number of positions to be used
   counter1_real = counter_real + 8
   counter1_img = counter_img + 8
   ! number of digits energy_step
   energy=int(emax)
   counter_energy = 0
   do while (energy .ne. 0)
   energy=energy/10
   counter_energy=counter_energy+1
   end do
   ! Redefinition of the counter_energy (sign + dot + 3=5)
   counter_energy=counter_energy + 5
   ! number of digits of qpoints
   qpoint=nq
   counter_qpoint = 0
   do while (qpoint .ne. 0)
   qpoint=qpoint/10
   counter_qpoint=counter_qpoint+1
   end do
   ! Write file
   do alpha=1,3
     do beta=1,3
        if (alpha==1 .and. beta==1) then
          write(file_id1,1001) "S_xx_component"
          write(file_id2,1001) "S_xx_component"
          do j=1,energy_step
            do i=1,nq
                 write(file_id1,1002)  qnumber(i), i, a+real(j)*b, real(S_mat(i,alpha,beta,j))
                 write(file_id2,1003)  qnumber(i), i, a+real(j)*b, aimag(S_mat(i,alpha,beta,j))
                 if (i==nq) then
                   write(file_id1,1001)
                   write(file_id2,1001)
                 end if
            end do
          end do
        else if (alpha==2 .and. beta==2) then
          write(file_id1,1001) "S_yy_component"
          write(file_id2,1001) "S_yy_component"
          do j=1,energy_step
            do i=1,nq
                 write(file_id1,1002)  qnumber(i), i, a+real(j)*b, real(S_mat(i,alpha,beta,j))
                 write(file_id2,1003)  qnumber(i), i, a+real(j)*b, aimag(S_mat(i,alpha,beta,j))
                 if (i==nq) then
                   write(file_id1,1001)
                   write(file_id2,1001)
                 end if
            end do
          end do
        else if (alpha==3 .and. beta==3) then
          write(file_id1,1001) "S_zz_component"
          write(file_id2,1001) "S_zz_component"
          do j=1,energy_step
            do i=1,nq
                 write(file_id1,1002)   qnumber(i), i, a+real(j)*b, real(S_mat(i,alpha,beta,j))
                 write(file_id2,1003)   qnumber(i), i, a+real(j)*b, aimag(S_mat(i,alpha,beta,j))
                 if (i==nq) then
                   write(file_id1,1001)
                   write(file_id2,1001)
                 end if
            end do
          end do
        else if (alpha==1 .and. beta==2) then
          write(file_id1,1001) "S_xy_component"
          write(file_id2,1001) "S_xy_component"
          do j=1,energy_step
            do i=1,nq
                 write(file_id1,1002)  qnumber(i), i, a+real(j)*b, real(S_mat(i,alpha,beta,j))
                 write(file_id2,1003)  qnumber(i), i, a+real(j)*b, aimag(S_mat(i,alpha,beta,j))
                 if (i==nq) then
                   write(file_id1,1001)
                   write(file_id2,1001)
                 end if
            end do
          end do
        else if (alpha==1 .and. beta==3) then
          write(file_id1,1001) "S_xz_component"
          write(file_id2,1001) "S_xz_component"
          do j=1,energy_step
            do i=1,nq
                 write(file_id1,1002)  qnumber(i), i, a+real(j)*b, real(S_mat(i,alpha,beta,j))
                 write(file_id2,1003)  qnumber(i), i, a+real(j)*b, aimag(S_mat(i,alpha,beta,j))
                 if (i==nq) then
                   write(file_id1,1001)
                   write(file_id2,1001)
                 end if
            end do
          end do
        else if (alpha==2 .and. beta==1) then
          write(file_id1,1001) "S_yx_component"
          write(file_id2,1001) "S_yx_component"
          do j=1,energy_step
            do i=1,nq
                 write(file_id1,1002)  qnumber(i), i, a+real(j)*b, real(S_mat(i,alpha,beta,j))
                 write(file_id2,1003)  qnumber(i), i, a+real(j)*b, aimag(S_mat(i,alpha,beta,j))
                 if (i==nq) then
                   write(file_id1,1001)
                   write(file_id2,1001)
                 end if
            end do
          end do
        else if (alpha==2 .and. beta==3) then
          write(file_id1,1001) "S_yz_component"
          write(file_id2,1001) "S_yz_component"
          do j=1,energy_step
            do i=1,nq
                 write(file_id1,1002)  qnumber(i), i, a+real(j)*b, real(S_mat(i,alpha,beta,j))
                 write(file_id2,1003)  qnumber(i), i, a+real(j)*b, aimag(S_mat(i,alpha,beta,j))
                 if (i==nq) then
                   write(file_id1,1001)
                   write(file_id2,1001)
                 end if
            end do
          end do
        else if (alpha==3 .and. beta==1) then
          write(file_id1,1001) "S_zx_component"
          write(file_id2,1001) "S_zx_component"
          do j=1,energy_step
            do i=1,nq
                 write(file_id1,1002)  qnumber(i), i, a+real(j)*b, real(S_mat(i,alpha,beta,j))
                 write(file_id2,1003)  qnumber(i), i, a+real(j)*b, aimag(S_mat(i,alpha,beta,j))
                 if (i==nq) then
                   write(file_id1,1001)
                   write(file_id2,1001)
                 end if
            end do
          end do
        else if (alpha==3 .and. beta==2) then
          write(file_id1,1001) "S_zy_component"
          write(file_id2,1001) "S_zy_component"
            do j=1,energy_step
              do i=1,nq
                 write(file_id1,1002)  qnumber(i), i, a+real(j)*b, real(S_mat(i,alpha,beta,j))
                 write(file_id2,1003)  qnumber(i), i, a+real(j)*b, aimag(S_mat(i,alpha,beta,j))
                 if (i==nq) then
                   write(file_id1,1001)
                   write(file_id2,1001)
                 end if
              end do
           end do
        end if
     end do
   end do


   ! Close files
   close(file_id1)
   close(file_id2)
   ! Format file
   1001 format (a)
   1002 format (1x,f11.6,2x,i8,2x,9es16.8,2x,9es16.8,2x,9es16.8,2x)
   1003 format (1x,f11.6,2x,i8,2x,9es16.8,2x,9es16.8,2x,9es16.8,2x)
   ! end subroutine
 end subroutine print_dcf
 !

  
! subroutine calc_pow_spectra(S_mat)
!   implicit none
!   ! Global variables
!   complex(dblprec), dimension(nq,3,3,energy_step), intent(in) ::  S_mat
!
!   stepsize=(qvecmax-qvecmin)/qvecgrid
!
!   do qvecmod=qvecmin,qvecmax,stepsize
!     call sphere_fibonacci_grid_points (qvecmod,numsphere,qvec_sphere)
!     call calc_dcf(emax,emin,energy_step,eig_val_2k,Tinv_mat_k,S_mat,qvector,0)
!   end do
!     ! Average S_mat
!     if (do_powder_spectra == 'Y') then
!           do alpha=1,3
!            do beta=1,3
!              S_mat_avrg(:,:)=S_mat_avrg(:,:)+S_mat(:,alpha,beta,:)
!            end do
!           end do
!           S_mat_avrg(:,:)=S_mat_avrg(:,:)/9
!     end if
! end subroutine calc_pow_spectra
! !
! subroutine sphere_fibonacci_grid_points (qvecmod,numsphere,xyz)
!   !
!   ! SPHERE_FIBONACCI_GRID_POINTS computes sphere points on a Fibonacci spiral.
!   !
!   !  Parameters:
!   !
!   !    Input, integer  NG, the number of points.
!   !
!   !    Output, real  XYZ(3,NG), the Fibonacci spiral points.
!   !
!   use Constants
!   !
!   implicit none
!   ! Global variables
!   integer, intent(in) :: numsphere !< Number of points over the sphere
!   real(dblprec), intent(in) :: qvecmod !< Module of the qvector (radius of the sphere)
!   real(dblprec), dimension(3,numsphere), intent(out) :: qvec_sphere !< Fibonacci spiral points over the surface of the sphere
!   ! Local variable
!   real (dblprec) :: cphi
!   integer :: i
!   real(dblprec) :: i_r8
!   integer :: j
!   real(dblprec) :: numsphere_r8
!   real(dblprec) :: r8_phi
!   real (dblprec) sphi
!   real (dblprec) theta
!
!   r8_phi = ( 1.0D+00 + sqrt ( 5.0D+00 ) ) / 2.0D+00
!   numsphere_r8 = real ( numsphere, kind = 8 )
!
!   do j = 1, numsphere
!     i_r8 = real ( - numsphere - 1 + 2 * j, kind = 8 )
!     theta = 2.0D+00 * pi * i_r8 / r8_phi
!     sphi = i_r8 / numsphere_r8
!     cphi = sqrt ( ( numsphere_r8 + i_r8 ) * ( numsphere_r8 - i_r8 ) ) / numsphere_r8
!     qvec_sphere(1,j) = qvecmod * cphi * sin ( theta )
!     qvec_sphere(2,j) = qvecmod * cphi * cos ( theta )
!     qvec_sphere(3,j) = qvecmod * sphi
!   end do
!
!   return
! end
!


 subroutine read_parameters_diamag_full(ifile)
      use FileParser
      
      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword,cache, string
      integer :: rd_len,i_err,i,i_stat,i_errb,ii, i_all
      logical :: comment
      real(dblprec) :: tmp
 
      q_vector=0.0_dblprec
      rot_axis(1:2)=0.0_dblprec
      rot_axis(3)=1.0_dblprec
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

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! START OF VARIABLES 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            ! Flags for the Non-Collinear Linear Spin Wave Theory subroutine (diamag_full)
       case('do_diamag_full')
          read(ifile,*,iostat=i_err) do_diamag_full
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
       case('energy_step')
          read(ifile,*,iostat=i_err) energy_step
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
       case('sigma_corr')
          read(ifile,*,iostat=i_err) sigma_corr
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err       
       case('do_correlation')
          read(ifile,*,iostat=i_err) do_correlation
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
       case('do_norm_corr')
          read(ifile,*,iostat=i_err) do_norm_corr
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
       case('q_vector')
          read(ifile,*,iostat=i_err) q_vector
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
       case('rot_axis')
          read(ifile,*,iostat=i_err) rot_axis
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
       case('temp_corr')
          read(ifile,*,iostat=i_err) temp_corr
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

           

            end select
         end if

         ! End of file
         if (i_errb==20) goto 20
         ! End of row
         if (i_errb==10) goto 10
      end do

      20  continue

      return
   end subroutine read_parameters_diamag_full 

end module diamag_full
