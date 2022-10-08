!> Data and routines for calculate frequency based spin correlation function S(w)
!> @authors
!> Anders Bergman, Manuel Pereiro
!> @copyright
!> GNU Public License
!!
module diamag
   use Parameters
   use Profiling
   use Hamiltoniandata,    only : ham
   !
   implicit none
   !
   !
   real(dblprec), dimension(:,:), allocatable :: magham !< M(w)
   complex(dblprec), dimension(:,:,:,:), allocatable :: magmom_qw !< M(w) in reciprocal space
   complex(dblprec), dimension(:,:), allocatable :: A_k
   complex(dblprec), dimension(:,:), allocatable :: A_km
   complex(dblprec), dimension(:,:), allocatable :: B_k
   complex(dblprec), dimension(:,:), allocatable :: C_k
   complex(dblprec), dimension(:,:), allocatable :: h_k

   !
   character(len=1) :: do_diamag    !< Perform frequency based spin-correlation sampling (Y/N/C)
   real(dblprec)    :: diamag_mix   !< Separation between sampling steps
   real(dblprec)    :: diamag_thresh   !< Separation between sampling steps
   integer          :: diamag_niter !< Number of steps to sample

   private
   ! public subroutines
   public :: do_diamag, read_parameters_diamag
   public ::  setup_diamag, setup_finite_hamiltonian, setup_infinite_hamiltonian
   public :: setup_tensor_hamiltonian, setup_altern_hamiltonian

contains

   subroutine setup_diamag()

      implicit none

      do_diamag='Y'
      diamag_niter=1000
      diamag_mix=0.030_dblprec
      diamag_thresh=1.0d-8

   end subroutine setup_diamag

   !> Set up the Hamiltonian for first cell
   subroutine setup_finite_hamiltonian(N1,N2,N3,NT,NA,Natom, nHam, Mensemble, conf_num, simid, emomM, &
         mmom, max_no_neigh, nlistsize, nlist, ncoup, kaniso )

      use Constants
      !
      implicit none
      !
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: NT                                          !< Number of types of atoms
      integer, intent(in) :: NA                                          !< Number of atoms in one cell
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: nHam      !< Number of atoms in Hamiltonian
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector

      !
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(nHam), intent(in) :: nlistsize     !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist       !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup  !< Heisenberg exchange couplings
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso  !< Anisotropy constant
      !
      integer :: i_stat, ia, ja, j, hdim, l, la
      complex(dblprec), dimension(3) :: ui, vi, uj, vj, ul, vl
      complex(dblprec) :: udot, ucdot, vdot
      complex(dblprec), dimension(NA,NA) :: dot_mat
      !
      complex(dblprec), dimension(:,:), allocatable :: eig_vec
      real(dblprec), dimension(:), allocatable :: eig_val
      !
      hdim=2*NA
      allocate(magham(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(magham))*kind(magham),'magham','setup_finite_hamiltonian')
      magham=0.0_dblprec
      allocate(A_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_k))*kind(A_k),'A_k','setup_finite_hamiltonian')
      allocate(B_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(B_k))*kind(B_k),'B_k','setup_finite_hamiltonian')
      allocate(C_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(C_k))*kind(C_k),'C_k','setup_finite_hamiltonian')
      allocate(eig_vec(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_vec))*kind(eig_vec),'eig_vec','setup_finite_hamiltonian')
      allocate(eig_val(hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_val))*kind(eig_val),'eig_val','setup_finite_hamiltonian')

      call setup_diamag()

      !Toth and Lake
      A_k=0.0_dblprec;B_k=0.0_dblprec;C_k=0.0_dblprec;dot_mat=0.0_dblprec
      do ia=1,NA
         call find_uv(ui,vi,emomM(:,ia,1))
         !  print *,'Atom',ia
         do j=1,nlistsize(ia)
            ja=nlist(j,ia)
            call find_uv(uj,vj,emomM(:,ja,1))
            ucdot=ui(1)*conjg(uj(1))+ui(2)*conjg(uj(2))+ui(3)*conjg(uj(3))
            udot=ui(1)*uj(1)+ui(2)*uj(2)+ui(3)*uj(3)
            A_k(ia,ja)=0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ncoup(j,ia,1)*ucdot
            B_k(ia,ja)=0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ncoup(j,ia,1)*udot
         end do
         C_k(ia,ia)=0.0_dblprec
         do l=1,nlistsize(ia)
            la=nlist(l,ia)
            call find_uv(ul,vl,emomM(:,la,1))
            vdot=vi(1)*vl(1)+vi(2)*vl(2)+vi(3)*vl(3)
            C_k(ia,ia)=C_k(ia,ia)+1.0_dblprec*mmom(la,1)*ncoup(l,ia,1)*vdot
         end do
      end do

      allocate(h_k(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(h_k))*kind(h_k),'h_k','setup_finite_hamiltonian')
      h_k=0.0_dblprec
      do ia=1,NA
         do ja=1,NA
            ! Diagonal as-is
            h_k(ia,ja)=A_k(ia,ja)-C_k(ia,ja)
            ! Off-diagonal Transpose and conjugate
            h_k(ia+NA,ja)=conjg(B_k(ja,ia))
            !!!h_k(ia+NA,ja)=B_k(ia,ja)
            ! Off-diagonal as-is
            h_k(ia,ja+NA)=B_k(ia,ja)
            !!!h_k(ia,ja+NA)=conjg(B_k(ja,ia))
            ! Diagonal Conjugate
            h_k(ia+NA,ja+NA)=conjg(A_k(ia,ja))-C_k(ia,ja)
         end do
         ! Add diagonal epsilon to ensure positive definiteness
         !h_k(ia,ia)=h_k(ia,ia)+1.0d-6
      end do

      call diagonalize_quad_hamiltonian(NA,h_k,eig_val,eig_vec)

      return
      !
   end subroutine setup_finite_hamiltonian

   !> Set up the Hamiltonian for first cell
   subroutine setup_infinite_hamiltonian(N1,N2,N3,NT,NA,Natom, nHam, Mensemble, conf_num, simid, emomM, &
         mmom,do_dm)

      use Constants
      use SystemData, only : coord, atype
      use AMS, only : wrap_coord_diff
      use Correlation,        only : q,nq
      !
      implicit none
      !
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: NT                                          !< Number of types of atoms
      integer, intent(in) :: NA                                          !< Number of atoms in one cell
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: nHam      !< Number of atoms in Hamiltonian
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)

      !
      real(dblprec), dimension(3) :: q_vec, q_plus, q_minus, q_0
      real(dblprec), dimension(3) :: dist, z
      real(dblprec), dimension(3) :: dmv
      real(dblprec), dimension(3,3) :: J_n, D_n, R_n, R_m
      complex(dblprec), dimension(3,3) :: J_prime
      !
      integer :: i_stat, ia, ja, j, hdim, l, la, iq, jat, lat
      complex(dblprec), dimension(3) :: ui, vi, uj, vj, ul, vl
      complex(dblprec) :: udot, ucdot, vdot, im
      complex(dblprec), dimension(NA,NA) :: dot_mat
      !
      complex(dblprec), dimension(:,:), allocatable :: eig_vec
      real(dblprec), dimension(:), allocatable :: eig_val
      real(dblprec), dimension(:,:), allocatable :: q_ext
      !
      hdim=2*NA
      allocate(magham(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(magham))*kind(magham),'magham','setup_infinite_hamiltonian')
      allocate(q_ext(3,0:2*nq),stat=i_stat)
      call memocc(i_stat,product(shape(q_ext))*kind(q_ext),'q_ext','setup_infinite_hamiltonian')
      magham=0.0_dblprec
      allocate(A_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_k))*kind(A_k),'A_k','setup_infinite_hamiltonian')
      allocate(A_km(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_km))*kind(A_km),'A_km','setup_infinite_hamiltonian')
      allocate(B_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(B_k))*kind(B_k),'B_k','setup_infinite_hamiltonian')
      allocate(C_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(C_k))*kind(C_k),'C_k','setup_infinite_hamiltonian')
      allocate(eig_vec(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_vec))*kind(eig_vec),'eig_vec','setup_infinite_hamiltonian')
      allocate(eig_val(hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_val))*kind(eig_val),'eig_val','setup_infinite_hamiltonian')
      allocate(h_k(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(h_k))*kind(h_k),'h_k','setup_infinite_hamiltonian')

      call setup_diamag()

      call clone_q(nq,q,q_ext)


      z(1)=0.0_dblprec;z(2)=0.0_dblprec;z(3)=1.0_dblprec


      im=(0.0_dblprec,1.0_dblprec)


      !Toth and Lake
      do iq=1,nq
         q_vec=q(:,iq)*2.0_dblprec*pi
         q_plus=q(:,iq)*2.0_dblprec*pi
         q_minus=-q(:,iq)*2.0_dblprec*pi
         q_0=0.0_dblprec
         !print '(2x,a,3f12.6)' , 'q=',q_vec

         A_k=0.0_dblprec;A_km=0.0_dblprec;B_k=0.0_dblprec;C_k=0.0_dblprec;dot_mat=0.0_dblprec
         do ia=1,NA
            call find_uv(ui,vi,emomM(:,ia,1))
            
            ! Jij exchange
            do j=1,ham%nlistsize(ia)
               ja=ham%nlist(j,ia)
               jat=atype(ja)
               call find_uv(uj,vj,emomM(:,ja,1))
               !print *,'---u_i---',iq
               !print '(6f12.6)',ui
               !print *,'---u_j---',iq
               !print '(6f12.6)',uj
               !print *,'-C(u_j)--',iq
               !print '(6f12.6)',conjg(uj)
               call wrap_coord_diff(Natom,coord,ia,ja,dist)
               !print '(3x,7f12.6)', dist, im,exp(-im*(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))


               !!! ucdot=ui(1)*conjg(uj(1))+ui(2)*conjg(uj(2))+ui(3)*conjg(uj(3))
               !!! udot=ui(1)*uj(1)+ui(2)*uj(2)+ui(3)*uj(3)
               !!! A_k(ia,jat)=A_k(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ham%ncoup(j,ia,1)*ucdot &
               !!!    *exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))
               !!! A_km(ia,jat)=A_km(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ham%ncoup(j,ia,1)*ucdot &
               !!!    *exp(-im *(q_plus(1)*dist(1)+q_plus(2)*dist(2)+q_plus(3)*dist(3)))
               !!! B_k(ia,jat)=B_k(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ham%ncoup(j,ia,1)*udot &
               !!!    *exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))
               ! print '(a,8f12.6)' ,'Scalar: ', ham%ncoup(j,ia,1)*ucdot,ham%ncoup(j,ia,1)*udot, exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))

               !!! print *,'----m_i---'
               !!! print '(6f12.6)', emomM(:,ia,1)
               !!! print *,'----m_ja---'
               !!! print '(6f12.6)', emomM(:,ja,1)
               !!! print *,'---m_jat---'
               !!! print '(6f12.6)', emomM(:,jat,1)
               !!! print *,'----u_i---'
               !!! print '(6f12.6)', ui
               !!! print *,'----u_j---'
               !!! print '(6f12.6)', uj
               !!! print *,'--c(u_j)--'
               !!! print '(6f12.6)', conjg(uj)

                J_n=0.0_dblprec; J_n(1,1)=ham%ncoup(j,ia,1); J_n(2,2)=ham%ncoup(j,ia,1); J_n(3,3)=ham%ncoup(j,ia,1)
                call find_R(R_m,z,emomM(:,ia,1))
                !print *,'----R_m---'
                !print '(3f12.6)', R_m

                call find_R(R_n,emomM(:,jat,1),emomM(:,ja,1))
                !print *,'----R_n---'
                !print '(3f12.6)', R_n
                !call find_R(R_n,emomM(:,ia,1),emomM(:,ja/na+1,1))

                !J_prime=matmul(R_n,J_n)
                !J_prime=matmul(matmul(J_n,R_n),R_m)
                J_prime=J_n
                !print *,'--J_prime--'
                !print '(6f12.6)', J_prime
                !call find_R(R_n,emomM(:,ia,1),emomM(:,ja/na+1,1))

                A_k(ia,jat)=A_k(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,J_prime,conjg(uj)) &
                   *exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))
                A_km(ia,jat)=A_km(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,J_prime,conjg(uj)) &
                   *exp(-im *(q_plus(1)*dist(1)+q_plus(2)*dist(2)+q_plus(3)*dist(3)))
                B_k(ia,jat)=B_k(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,J_prime,uj) &
                   *exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))

                !!! print *,'----A,B----'
                !!! print *,sJs(ui,J_prime,conjg(uj)) 
                !!! print *,sJs(ui,J_prime,uj) 
               !!! ! !print '(a,8f12.6)' ,'Tensor: ', sJs(ui,J_prime,conjg(uj)) , sJs(ui,J_prime,uj), exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))

            end do
            !stop
            ! Jij diagonal
            do l=1,ham%nlistsize(ia)
               la=ham%nlist(l,ia)
               lat=atype(la)
               call find_uv(ul,vl,emomM(:,la,1))

               !vdot=vi(1)*vl(1)+vi(2)*vl(2)+vi(3)*vl(3)
               !C_k(ia,ia)=C_k(ia,ia)+1.0_dblprec*mmom(la,1)*ham%ncoup(l,ia,1)* vdot

               J_n=0.0_dblprec; J_n(1,1)=ham%ncoup(l,ia,1); J_n(2,2)=ham%ncoup(l,ia,1); J_n(3,3)=ham%ncoup(l,ia,1)
               call find_R(R_n,emomM(:,lat,1),emomM(:,la,1))
               !call find_R(R_n,emomM(:,ia,1),emomM(:,la/na+1,1))
               !J_prime=matmul(J_n,R_n)
               !J_prime=matmul(matmul(J_n,R_n),R_m)
               J_prime=J_n

               C_k(ia,ia)=C_k(ia,ia)+1.0_dblprec*mmom(la,1)*sJs(vi,J_prime,vl) 
            end do
            !!! DMI
            ! if(do_dm==1) then
            !    do j=1,ham%dmlistsize(ia)
            !       ja=ham%dmlist(j,ia)
            !       jat=atype(ja)
            !       dmv=ham%dm_vect(:,j,ia)
            !       call find_uv(uj,vj,emomM(:,ja,1))
            !       call wrap_coord_diff(Natom,coord,ia,ja,dist)

            !       !A_k(ia,jat)=A_k(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))* &
            !       !   ( dmv(1)*ui(3)*conjg(uj(2))-dmv(1)*ui(2)*conjg(uj(3)) + &
            !       !   dmv(2)*ui(1)*conjg(uj(3))-dmv(2)*ui(3)*conjg(uj(1)) + &
            !       !   dmv(3)*ui(2)*conjg(uj(1))-dmv(3)*ui(1)*conjg(uj(2)) ) &
            !       !   *exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))
            !       !A_km(ia,jat)=A_km(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))* &
            !       !   ( dmv(1)*ui(3)*conjg(uj(2))-dmv(1)*ui(2)*conjg(uj(3)) + &
            !       !   dmv(2)*ui(1)*conjg(uj(3))-dmv(2)*ui(3)*conjg(uj(1)) + &
            !       !   dmv(3)*ui(2)*conjg(uj(1))-dmv(3)*ui(1)*conjg(uj(2)) ) &
            !       !   *exp(-im *(q_plus(1)*dist(1)+q_plus(2)*dist(2)+q_plus(3)*dist(3)))
            !       !B_k(ia,jat)=B_k(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))* &
            !       !   ( dmv(1)*ui(3)*(uj(2))-dmv(1)*ui(2)*(uj(3)) + &
            !       !   dmv(2)*ui(1)*(uj(3))-dmv(2)*ui(3)*(uj(1)) + &
            !       !   dmv(3)*ui(2)*(uj(1))-dmv(3)*ui(1)*(uj(2)) ) &
            !       !   *exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))

            !       call find_R(R_n,emomM(:,ia,1),emomM(:,ja,1))
            !       D_n=dm2tens(dmv)*0.1
            !       J_prime=matmul(D_n,R_n)
            !       J_prime=D_n


            !    A_k(ia,jat)=A_k(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,J_prime,conjg(uj)) &
            !       *exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))
            !    A_km(ia,jat)=A_km(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,J_prime,conjg(uj)) &
            !       *exp(-im *(q_plus(1)*dist(1)+q_plus(2)*dist(2)+q_plus(3)*dist(3)))
            !    B_k(ia,jat)=B_k(ia,jat)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,J_prime,uj) &
            !       *exp(-im *(q_minus(1)*dist(1)+q_minus(2)*dist(2)+q_minus(3)*dist(3)))

            !       !print '(3x,3i5,2f12.6)', ia,j, ja, 0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ncoup(j,ia,1)*ucdot
            !       !print '(2f12.6)', exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))
            !       !print '(2f12.6)', A_k(ia,jat)
            !    end do
            !    ! DMI diagonal
            !    do l=1,ham%nlistsize(ia)
            !       la=ham%nlist(l,ia)
            !       lat=atype(la)
            !       dmv=ham%dm_vect(:,l,ia)
            !       call find_uv(ul,vl,emomM(:,la,1))
            !       !C_k(ia,ia)=C_k(ia,ia)+1.0_dblprec*mmom(la,1)* &
            !       !   ( dmv(1)*vi(3)*(vl(2))-dmv(1)*vi(2)*(vl(3)) + &
            !       !   dmv(2)*vi(1)*(vl(3))-dmv(2)*vi(3)*(vl(1)) + &
            !       !   dmv(3)*vi(2)*(vl(1))-dmv(3)*vi(1)*(vl(2)) ) 
            !       !if(lat==ia) C_k(ia,ia)=C_k(ia,ia)+1.0_dblprec*mmom(la,1)*ncoup(l,ia,1)*vdot

            !       call find_R(R_n,emomM(:,ia,1),emomM(:,la,1))
            !       D_n=dm2tens(dmv)*0.1
            !       J_prime=matmul(D_n,R_n)
            !       J_prime=D_n

            !       C_k(ia,ia)=C_k(ia,ia)+1.0_dblprec*mmom(la,1)*sJs(vi,J_prime,vl)  

            !    end do
            ! end if
         end do
         !print *,'----Ak---',q_vec(3),sum(ham%ncoup(:,1,1))
         !print '(4f12.6)',A_k
         !print *,'----Bk---',real(A_k(1,1))
         !print '(4f12.6)',B_k
         !print *,'----Ck---', aimag(A_k(1,1))
         !print '(4f12.6)',C_k
         !print *,'---------'

         !stop

         h_k=0.0_dblprec
         do ia=1,NA
            do ja=1,NA
               ! Diagonal as-is
               h_k(ia,ja)=A_k(ia,ja)-C_k(ia,ja)
               ! Off-diagonal Transpose and conjugate
               h_k(ia+NA,ja)=conjg(B_k(ja,ia))
               !!!h_k(ia+NA,ja)=B_k(ia,ja)
               ! Off-diagonal as-is
               h_k(ia,ja+NA)=B_k(ia,ja)
               !!!h_k(ia,ja+NA)=conjg(B_k(ja,ia))
               ! Diagonal Conjugate
               h_k(ia+NA,ja+NA)=conjg(A_km(ia,ja))-C_k(ia,ja)
            end do
            ! Add diagonal epsilon to ensure positive definiteness
            !h_k(ia,ia)=h_k(ia,ia)+1.0d-6
         end do
         !print *,'----h_k--'
         !print '(4f12.6)',h_k

         call diagonalize_quad_hamiltonian(NA,h_k,eig_val,eig_vec)

         write(1000,'(i6,12f12.6)') iq,real(eig_val) *13.605698066_dblprec*4.0_dblprec*2.0_dblprec
         !print '(a,i6,12f12.6)','ev: ',iq,real(eig_val)*13.605698066_dblprec*4.0_dblprec

      end do

      return
      !
   end subroutine setup_infinite_hamiltonian

   !> Set up the Hamiltonian for first cell
   subroutine setup_tensor_hamiltonian(N1,N2,N3,NT,NA,Natom, Mensemble, simid, emomM, mmom,do_dm)

      use Constants
      use SystemData, only : coord, atype
      use AMS, only : wrap_coord_diff
      use Correlation,        only : q,nq
      !
      implicit none
      !
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)

      !
      real(dblprec), dimension(3) :: q_vec
      real(dblprec), dimension(3) :: dist
      real(dblprec), dimension(3) :: dmv
      real(dblprec), dimension(3) :: z
      !
      integer :: i_stat, ia, ja, j, hdim, l, la, iq, jat, lat, nq_ext
      complex(dblprec), dimension(3) :: ui, vi, uj, vj, ul, vl
      complex(dblprec) :: udot, ucdot, vdot, im
      complex(dblprec), dimension(NA,NA) :: dot_mat
      !
      complex(dblprec), dimension(:,:), allocatable :: eig_vec
      real(dblprec), dimension(:), allocatable :: eig_val
      complex(dblprec), dimension(:,:,:,:,:), allocatable ::   jtens_q  !< FT of Exchange tensor
      real(dblprec), dimension(:,:), allocatable :: q_ext
      !
      hdim=2*NA
      nq_ext=2*nq
      allocate(jtens_q(3,3,NA,NA,0:nq_ext),stat=i_stat)
      call memocc(i_stat,product(shape(jtens_q))*kind(jtens_q),'jtens_q','setup_tensor_hamiltonian')
      !
      allocate(magham(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(magham))*kind(magham),'magham','setup_tensor_hamiltonian')
      allocate(q_ext(3,0:nq_ext),stat=i_stat)
      call memocc(i_stat,product(shape(q_ext))*kind(q_ext),'q_ext','setup_tensor_hamiltonian')
      magham=0.0_dblprec
      allocate(A_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_k))*kind(A_k),'A_k','setup_tensor_hamiltonian')
      allocate(A_km(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_km))*kind(A_km),'A_km','setup_tensor_hamiltonian')
      allocate(B_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(B_k))*kind(B_k),'B_k','setup_tensor_hamiltonian')
      allocate(C_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(C_k))*kind(C_k),'C_k','setup_tensor_hamiltonian')
      allocate(eig_vec(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_vec))*kind(eig_vec),'eig_vec','setup_tensor_hamiltonian')
      allocate(eig_val(hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_val))*kind(eig_val),'eig_val','setup_tensor_hamiltonian')
      allocate(h_k(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(h_k))*kind(h_k),'h_k','setup_tensor_hamiltonian')

      call setup_diamag()

      im=(0.0_dblprec,1.0_dblprec)

      call clone_q(nq,q,q_ext)
      !print *,' test q: '
      !print '(3f10.3)',q_ext(:,0:min(3,nq))

      !print *, 'setup_Jtens_q ', nq_ext,nq
      !print '(3f12.5)',q_ext

      z(1)=0.0_dblprec;z(2)=0.0_dblprec;z(3)=1.0_dblprec

      call setup_Jtens_q(Natom,Mensemble,NA,emomM,do_dm,q_ext,nq_ext,Jtens_q)

      !print *,' some q: '
      !print '(3f10.3)',q_ext(:,0:min(3,nq))
      !Toth and Lake
      do iq=1,nq
         !q_vec=q(:,iq)*2.0_dblprec*pi
         !print '(2x,a,3f12.6)' , 'q=',q_ext(:,iq)*2.0_dblprec*pi
         !print '(2x,a,3f12.6)' , 'q_min=',q_ext(:,iq+nq)*2.0_dblprec*pi

         A_k=0.0_dblprec;A_km=0.0_dblprec;B_k=0.0_dblprec;C_k=0.0_dblprec;dot_mat=0.0_dblprec
         do ia=1,NA
            !call find_uv(ui,vi,z)
            call find_uv(ui,vi,emomM(:,ia,1))
            do ja=1,NA
               !call find_uv(uj,vj,z)
               call find_uv(uj,vj,emomM(:,ja,1))
         !   
               !!!   print *,'---u_i---',ia
               !!!   print '(6f12.6)',ui
               !!!   print *,'---v_i---',ia
               !!!   print '(6f12.6)',vi
               !!!   print *,'---u_j---',ja
               !!!   print '(6f12.6)',uj
               !!!   print *,'-conj(u_j)--',ja
               !!!   print '(6f12.6)',conjg(uj)
               !!!   print *,'---v_j---',ja
               !!!   print '(6f12.6)',vj
               !!!   print *,'----Jq---',ia,ja,iq
               !!!   print '(6f12.6)',Jtens_q(:,:,ia,ja,iq)
               !!!   print *,'---J-q---',ia,ja,iq+nq
               !!!   print '(6f12.6)',Jtens_q(:,:,ia,ja,iq+nq)
               !!!   print *,'---A_ij--',ia,ja
               !!!   print '(6f12.6)',sJs(ui,Jtens_q(:,:,ia,ja,iq+nq),conjg(uj))
               !!!   print *,'---Ak_ij--',ia,ja
               !!!   print '(6f12.6)',sJs(ui,Jtens_q(:,:,ia,ja,iq),conjg(uj))
               !!!   print *,'---B_ij--',ja,ia
               !!!   print '(6f12.6)',sJs(ui,Jtens_q(:,:,ia,ja,iq),uj)
               !!!   print *,'---C_ij--',ja,ia
               !!!   print '(6f12.6)',sJs(vi,Jtens_q(:,:,ia,ja,0),vj)
               !!!   print *,'------------------------------------------------'
               A_k(ia,ja) =A_k(ia,ja) +0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,Jtens_q(:,:,ia,ja,iq+nq),conjg(uj))
               A_km(ia,ja)=A_km(ia,ja)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,Jtens_q(:,:,ia,ja,iq),conjg(uj))
               B_k(ia,ja) =B_k(ia,ja) +0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,Jtens_q(:,:,ia,ja,iq+nq),uj)
               C_k(ia,ia) =C_k(ia,ia) +1.0_dblprec*mmom(ja,1)*sJs(vi,Jtens_q(:,:,ia,ja,0),vj)
            end do
         end do

         !!! print *,'---Ak---',q_vec(3)
         !!! print '(4f12.6)',A_k
         !!! print *,'---Akm--',sum(ham%ncoup(:,1,1))
         !!! print '(4f12.6)',A_km
         !!! print *,'----Bk---',real(A_k(1,1))
         !!! print '(4f12.6)',B_k
         !!! print *,'----Ck---', aimag(A_k(1,1))
         !!! print '(4f12.6)',C_k
         !!! print *,'---------'

         h_k=0.0_dblprec
         do ia=1,NA
            do ja=1,NA
               ! Diagonal as-is
               h_k(ia,ja)=A_k(ia,ja)-C_k(ia,ja)
               ! Off-diagonal Transpose and conjugate
               h_k(ia+NA,ja)=conjg(B_k(ja,ia))
               !!!h_k(ia+NA,ja)=B_k(ia,ja)
               ! Off-diagonal as-is
               h_k(ia,ja+NA)=B_k(ia,ja)
               !!!h_k(ia,ja+NA)=conjg(B_k(ja,ia))
               ! Diagonal transpose
               h_k(ia+NA,ja+NA)=(A_km(ja,ia))-C_k(ia,ja)
               ! Diagonal Conjugate
               !h_k(ia+NA,ja+NA)=conjg(A_km(ia,ja))-C_k(ia,ja)
            end do
            !write(201,'(10f18.12)') real(A_k(ia,:))
            !write(202,'(10f18.12)') real(B_k(ia,:))
            !write(203,'(10f18.12)') real(C_k(ia,:))
            !write(301,'(10f18.12)') aimag(A_k(ia,:))
            !write(302,'(10f18.12)') aimag(B_k(ia,:))
            !write(303,'(10f18.12)') aimag(C_k(ia,:))
            ! Add diagonal epsilon to ensure positive definiteness
            !h_k(ia,ia)=h_k(ia,ia)+1.0d-5
         end do
         write(4000,'(100f18.10)') h_k

         !print *,'---h_k---'
         !print '(4f12.6)',h_k

         call diagonalize_quad_hamiltonian(NA,h_k,eig_val,eig_vec)

         write(1000,'(i6,12f12.6)') iq,real(eig_val)*13.605698066_dblprec*4.0_dblprec*2.0_dblprec
         !print '(a,i6,12f12.6)','ev: ',iq,real(eig_val)*13.605698066_dblprec*4.0_dblprec

      end do

      return
      !
   end subroutine setup_tensor_hamiltonian

   subroutine diagonalize_quad_hamiltonian(NA,h_in,eig_val,eig_vec)
      !
      use Constants
      !
      implicit none
      !
      !integer, intent(in) :: Natom     !< Number of atoms in system
      !integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA
      complex(dblprec), dimension(2*NA,2*NA) :: h_in
      complex(dblprec), dimension(2*NA,2*NA) :: eig_vec
      real(dblprec), dimension(2*NA) :: eig_val
      !real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
      !integer, intent(in) :: Mensemble !< Number of ensembles
      !real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      !
      !
      complex(dblprec), dimension(:,:), allocatable :: g_mat
      complex(dblprec), dimension(:,:), allocatable :: U_mat
      complex(dblprec), dimension(:,:), allocatable :: T_mat
      complex(dblprec), dimension(:,:), allocatable :: E_mat
      complex(dblprec), dimension(:,:), allocatable :: sqE_mat
      complex(dblprec), dimension(:,:), allocatable :: L_mat
      complex(dblprec), dimension(:,:), allocatable :: K_mat
      complex(dblprec), dimension(:,:), allocatable :: iK_mat
      complex(dblprec), dimension(:,:), allocatable :: dum_mat, dum_mat2
      complex(dblprec), dimension(:,:), allocatable :: x_mat
      complex(dblprec), dimension(:,:,:), allocatable :: bigS
      complex(dblprec), dimension(:), allocatable :: cwork
      real(dblprec), dimension(:), allocatable :: rwork
      !
      integer :: info, lwork, hdim, i_stat, ia, ja

      complex(dblprec) :: cone, czero, fcinv, dia_eps, im

      complex(dblprec), dimension(3) :: ul, vl
      !
      !
      hdim=2*NA
      czero=(0.0_dblprec,0.0_dblprec)
      cone=(1.0_dblprec,0.0_dblprec)
      im=(0.0_dblprec,1.0_dblprec)
      dia_eps=1.0d-12
      !
      fcinv = 0.5_dblprec*mub/mry
      !
      allocate(g_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(g_mat))*kind(g_mat),'g_mat','diagonalize_quad_hamiltonian')
      allocate(U_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(U_mat))*kind(U_mat),'U_mat','diagonalize_quad_hamiltonian')
      allocate(K_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(K_mat))*kind(K_mat),'K_mat','diagonalize_quad_hamiltonian')
      allocate(iK_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(iK_mat))*kind(iK_mat),'iK_mat','diagonalize_quad_hamiltonian')
      allocate(T_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(T_mat))*kind(T_mat),'T_mat','diagonalize_quad_hamiltonian')
      allocate(E_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(E_mat))*kind(E_mat),'E_mat','diagonalize_quad_hamiltonian')
      allocate(sqE_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(sqE_mat))*kind(sqE_mat),'sqE_mat','diagonalize_quad_hamiltonian')
      allocate(L_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(L_mat))*kind(L_mat),'L_mat','diagonalize_quad_hamiltonian')
      allocate(dum_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(dum_mat))*kind(dum_mat),'dum_mat','diagonalize_quad_hamiltonian')
      allocate(dum_mat2(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(dum_mat2))*kind(dum_mat2),'dum_mat2','diagonalize_quad_hamiltonian')
      allocate(x_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(x_mat))*kind(x_mat),'x_mat','diagonalize_quad_hamiltonian')
      allocate(bigS(3,NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(bigS))*kind(bigS),'bigS','diagonalize_quad_hamiltonian')
      !
      g_mat=(0.0_dblprec,0.0_dblprec)
      do ia=1,NA
         g_mat(ia,ia)=(1.0_dblprec,0.0_dblprec)
         g_mat(ia+NA,ia+NA)=(-1.0_dblprec,0.0_dblprec)
      end do
      !
      ! Cholesky decomposition of h to K'*K
      K_mat=-h_in*fcinv
      !print *,'-----real-----'
      !print '(6f20.10)',-real(h_in)
      !print *,'-----imag-----'
      !print '(6f20.10)',-aimag(h_in)
      !K_mat=-h_k*fcinv
      lwork=4*NA-1
      allocate(cwork(lwork))
      allocate(rwork(6*NA-2))
      eig_vec=K_mat
      call zheev('V','U',hdim,eig_vec,hdim,eig_val, cwork, lwork, rwork, info)
      write(3000,'(100g20.6)') eig_val
      !print *,'eig_val int',maxval(eig_val),minval(eig_val)
      ! Add offset to ensure positive definiteness, if needed.
      !dia_eps=0.0_dblprec
      !if(minval(eig_val)<0.0_dblprec) dia_eps=dia_eps-minval(eig_val)
      deallocate(cwork)
      deallocate(rwork)
      ! Add eps to diagonal to ensure positive definiteness
      do ia=1,hdim
         K_mat(ia,ia)=K_mat(ia,ia)+2.0_dblprec*dia_eps
      end do
      ! Cholesky
      call zpotrf('U',hdim,K_mat,hdim,info)
      if(info==0) then  ! Positive-definit matrix, Colpa diagonalization ok
         do ia=1,hdim
            do ja=ia+1,hdim
               K_mat(ja,ia)=0.0_dblprec
            end do
         end do
         call zgemm('N','C',hdim,hdim,hdim,cone,g_mat,hdim,K_mat,hdim,czero,dum_mat,hdim)
         call zgemm('N','N',hdim,hdim,hdim,cone,K_mat,hdim,dum_mat,hdim,czero,eig_vec,hdim)
         allocate(cwork(lwork))
         allocate(rwork(6*NA-2))
         ! Eigenvaluesolver for HgH'
         call zheev('V','U', hdim, eig_vec, hdim, eig_val, cwork, lwork, rwork, info)
         deallocate(cwork)
         deallocate(rwork)
      else
         print *,' Warning in diamag: non-positive definite matrix in zpotrf', info
      end if
      !print *,'zheeev',info
      call shuffle_eig(eig_val,eig_vec,hdim)
      !do ia=1,hdim
      !   print '(100f10.4)',real(eig_val(ia))
      !   !print '(100f10.4)',real(eig_vec(ia,:))
      !end do
      !ABs
      x_mat=eig_vec
      !do ia=1,NA
      !   do ja=1,NA
      !      call find_uv(ul,vl,emomM(:,ja,1))
      !      bigS(:,ia,ja)=sqrt(mmom(ia,1))*sqrt(0.5_dblprec)*(conjg(ul)*x_mat(ia,ja)+ul*x_mat(ia+NA,ja))+vl*(mmom(ia,1)-x_mat(ia+NA,ja)*x_mat(ia,ja))
      !   end do
      !end do

      ! Calculate L
      call zgemm('C','N',hdim,hdim,hdim,cone,K_mat,hdim,eig_vec,hdim,czero,dum_mat,hdim)
      call zgemm('N','N',hdim,hdim,hdim,cone,g_mat,hdim,dum_mat,hdim,czero,dum_mat2,hdim)
      call zgemm('N','N',hdim,hdim,hdim,cone,K_mat,hdim,dum_mat2,hdim,czero,dum_mat,hdim)
      call zgemm('C','N',hdim,hdim,hdim,cone,eig_vec,hdim,dum_mat,hdim,czero,L_mat,hdim)
      !
      ! Eigensolve E=g*L
      !print *,'zgemm'
      call zgemm('N','N',hdim,hdim,hdim,cone,g_mat,hdim,L_mat,hdim,czero,E_mat,hdim)
      ! Calculate K^-1
      iK_mat=K_mat
      call ztrtri('U','N',hdim,iK_mat,hdim,info)
      !print *,'ztrtri',info
      sqE_mat=sqrt(abs(E_mat))
      !print *,'zgemm'
      call zgemm('N','N',hdim,hdim,hdim,cone,eig_vec,hdim,sqE_mat,hdim,czero,dum_mat,hdim)
      !print *,'zgemm'
      call zgemm('N','N',hdim,hdim,hdim,cone,iK_mat,hdim,dum_mat,hdim,czero,T_mat,hdim)

      x_mat=0.0_dblprec
      do ia=1,hdim
         !call zgemv('N',hdim,hdim,cone,T_mat,hdim,eig_vec(ia,1:hdim),1,cone,x_mat(ia,1:hdim),1)
         !call zgemv('N',hdim,hdim,cone,T_mat,hdim,eig_vec(ia,1:hdim),1,cone,x_mat(1:hdim,ia),1)
         !call zgemv('N',hdim,hdim,cone,T_mat,hdim,eig_vec(1:hdim,ia),1,cone,x_mat(1:hdim,ia),1)
      end do
      x_mat=T_mat
      !print *,'Re x_mat'
      !open(ofileno,file='diaval.simid.out')
      do ia=1,hdim
         write(2000   , '(100f10.4)')eig_val(ia) *13.605698066_dblprec*4.0_dblprec*2.0_dblprec,real(x_mat(:,ia))
      end do
      !close(ofileno)
      !do ia=1,NA
      !   do ja=1,NA
      !      call find_uv(ul,vl,emomM(:,ja,1))
      !      bigS(:,ia,ja)=sqrt(mmom(ia,1))*sqrt(0.5_dblprec)*(conjg(ul)*x_mat(ia,ja)+ul*conjg(x_mat(ia,ja)))+vl*(mmom(ia,1)-conjg(x_mat(ia,ja))*x_mat(ia,ja))
      !   end do
      !end do
      return
      !
   end subroutine diagonalize_quad_hamiltonian


   !> Set up the Hamiltonian for first cell
   subroutine diagonalize_finite_hamiltonian(NA)

      use Constants
      !
      implicit none
      !
      integer, intent(in) :: NA                                          !< Number of atoms in one cell
      !
      integer :: i_stat
      integer :: lwork,info, hdim


      real(dblprec), dimension(:,:), allocatable :: Hprime, eig_vec
      real(dblprec), dimension(:), allocatable :: eig_val, work
      !

      allocate(Hprime(3*NA,3*NA),stat=i_stat)
      call memocc(i_stat,product(shape(Hprime))*kind(Hprime),'Hprime','diagonalize_finite_hamiltonian')
      !
      allocate(eig_vec(3*NA,3*NA),stat=i_stat)
      call memocc(i_stat,product(shape(eig_vec))*kind(eig_vec),'eig_vec','diagonalize_finite_hamiltonian')
      !!!
      allocate(eig_val(3*NA),stat=i_stat)
      call memocc(i_stat,product(shape(eig_val))*kind(eig_val),'eig_val','diagonalize_finite_hamiltonian')
      !!!
      lwork=NA*NA


      hdim=2*NA
      eig_vec= magham
      lwork=4*3*NA
      allocate(work(lwork),stat=i_stat)
      call memocc(i_stat,product(shape(work))*kind(work),'work','diagonalize_finite_hamiltonian')
      !
      call dsyev('V','L',hdim,eig_vec,hdim,eig_val,work,lwork,info)
      if(info.ne.0) then
         print '(2x,a,i4)', 'Problem in dsyev (diagonalize_finite_hamiltonian):',info
      end if
      !
      magham=magham
      eig_val=eig_val
      print *,'Eigenmoments 1'
      print '(1x,36f12.6)',eig_vec(:,1)
      !
      call find_eigenvector(hdim,eig_vec,eig_val)
      return
      !
   end subroutine diagonalize_finite_hamiltonian

   subroutine find_eigenvector(hdim,eig_vec,eig_val)

      use Constants
      !
      implicit none
      !
      integer, intent(in) :: hdim                   !< Number of atoms in one cell
      real(dblprec), dimension(hdim,hdim), intent(in) :: eig_vec
      real(dblprec), dimension(hdim), intent(in) :: eig_val
      !
      integer :: i_stat, ia, ivec, iter, i

      real(dblprec) :: enorm
      real(dblprec) :: alpha
      real(dblprec), dimension(:,:), allocatable :: Hprime, moments, err_norm
      real(dblprec), dimension(:), allocatable :: mom_vec, error
      !

      allocate(Hprime(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(Hprime))*kind(Hprime),'Hprime','find_eigenvector')
      !
      allocate(moments(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(moments))*kind(moments),'moments','find_eigenvector')
      !
      allocate(err_norm(hdim,diamag_niter),stat=i_stat)
      call memocc(i_stat,product(shape(err_norm))*kind(err_norm),'err_norm','find_eigenvector')
      !
      allocate(mom_vec(hdim),stat=i_stat)
      call memocc(i_stat,product(shape(mom_vec))*kind(mom_vec),'mom_vec','find_eigenvector')
      !
      allocate(error(hdim),stat=i_stat)
      call memocc(i_stat,product(shape(error))*kind(error),'error','find_eigenvector')

      do ivec=1,hdim
         Hprime=magham/maxval(abs(eig_val))
         alpha=1.0_dblprec
         do ia=1,hdim
            Hprime(ia,ia)=Hprime(ia,ia)-eig_val(ivec)/maxval(abs(eig_val))
         end do
         mom_vec=eig_vec(:,ivec);
         print *,'Intensities:'
         print '(1x,36f8.4)',(sqrt(sum(mom_vec((i-1)*3+1:i*3)**2)),i=1,hdim/3)
         call normalize_moments(mom_vec,hdim,1)
         print *,'IVEC: ',ivec
         print *,'MoMenTs:',eig_val(ivec)
         print '(1x,36f12.6)',(mom_vec)
         call dgemv('N',hdim,hdim,alpha,magham,hdim,eig_vec(:,ivec),1,0.0_dblprec,error,1)
         print *,'magham * vec :'
         print '(1x,36g12.4)',error
         call dgemv('N',hdim,hdim,alpha,Hprime,hdim,eig_vec(:,ivec),1,0.0_dblprec,error,1)
         print *,'Hprime * vec :'
         print '(1x,36g12.4)',error
         iter=1;
         enorm=1.0;
         do while (enorm>diamag_thresh.and.iter<=diamag_niter)
            call dgemv('N',hdim,hdim,alpha,Hprime,hdim,mom_vec,1,0.0_dblprec,error,1)
            enorm=sqrt(sum(error*error))
            print *,'Hprime * vec :',enorm
            print '(1x,36g12.4)',error
            if(sum(error*error)>1.0_dblprec*hdim/3.0_dblprec) call normalize_moments(error,hdim,1)
            mom_vec=mom_vec-1.00*error
            call normalize_moments(mom_vec,hdim,1)
            print *,'New trial vec:'
            print '(1x,36g12.4)',mom_vec
            err_norm(ivec,iter)=enorm;
            iter=iter+1;
         end do
         call dgemv('N',hdim,hdim,1.0_dblprec,Hprime,hdim,mom_vec,1,0.0_dblprec,error,1)
         print *,'Moments:',eig_val(ivec)
         print '(1x,36f12.6)',(mom_vec)
         print *,'Magham * vec :'
         print '(1x,36g12.4)',error
         call dgemv('N',hdim,hdim,alpha,Hprime,hdim,mom_vec,1,0.0_dblprec,error,1)
         print *,'Hprime * vec :'
         print '(1x,36g12.4)',error
         err_norm(ivec,diamag_niter)=enorm;
         moments(:,ivec)=mom_vec
      end do

      print *,'errors'
      print '(1x,9g12.4)',err_norm(:,diamag_niter)
      !
      call moments_to_angles(moments,hdim,eig_val)
      !
      return
      !
   end subroutine find_eigenvector

   subroutine moments_to_angles(moments,ndim,eig_val)
      !
      !
      implicit none
      !
      !
      integer, intent(in) :: ndim
      real(dblprec), dimension(ndim,ndim), intent(in) :: moments
      real(dblprec), dimension(ndim), intent(in) :: eig_val
      !
      real(dblprec) :: pi, mi_norm, mj_norm, mdot
      real(dblprec), dimension(ndim/3,ndim) :: angles
      real(dblprec), dimension(3) :: mom_i, mom_j

      integer :: i,j, idx

      angles=0.0_dblprec
      print *,'Angles:'
      pi=acos(-1.0_dblprec)
      do j=1,ndim
         idx=0
         mom_i(1:3)=moments(1:3,j)
         mi_norm=sqrt(sum(mom_i*mom_i))+1.0d-15
         do i=1,ndim,3
            mom_j(1:3)=moments(i:i+2,j)
            mj_norm=sqrt(sum(mom_j*mom_j))+1.0d-15
            mdot=sum(mom_i*mom_j)/(mi_norm*mj_norm)
            idx=idx+1
            angles(idx,j)=acos(mdot)*180.0_dblprec/pi
         end do
      end do
      !do j=1,ndim,3
      !   print '(1x,20f10.4)', eig_val(j), angles(:,j)
      !end do
      !
      return
      !
   end subroutine moments_to_angles

   subroutine normalize_moments(moments,nrow,ncol)
      !
      implicit none
      !
      !
      integer, intent(in) :: nrow
      integer, intent(in) :: ncol
      real(dblprec), dimension(ncol,nrow), intent(inout) :: moments
      !
      real(dblprec) :: mnorm
      integer :: i,j
      !
      do i=1,nrow,3
         do j=1,ncol
            mnorm=sum(moments(j,i:i+2)**2)**0.5_dblprec+1.0d-15
            moments(j,i+0)=moments(j,i+0)/mnorm
            moments(j,i+1)=moments(j,i+1)/mnorm
            moments(j,i+2)=moments(j,i+2)/mnorm
         end do
      end do
      !
      return
      !
   end subroutine normalize_moments

   subroutine find_uv(u,v,mom)
      !
      implicit none
      !
      complex(dblprec), dimension(3), intent(out) :: u
      complex(dblprec), dimension(3), intent(out) :: v
      real(dblprec), dimension(3), intent(in) :: mom
      !

      real(dblprec) :: mnorm
      complex(dblprec) :: im = (0.0_dblprec,1.0_dblprec)
      real(dblprec), dimension(3,3) :: R_prime
      real(dblprec), dimension(3) :: mom_hat
      !
      !
      mnorm=sqrt(sum(mom*mom)+1.0d-15)
      ! Find the rotation matrix R_prime that transforms mom to z-axis.
      ! Done by finding the vector orthogonal to mom and z and then perform a
      ! Rodrigues rotation with the corresponding angle
      !
      ! Unit vector parallel to mom
      mom_hat=mom/mnorm
      ! Follow Toth-Lake recipe for R'
      R_prime(3,:)=mom_hat
      if(abs(mom_hat(3))==1) then  !if m=00+-1 then R_2=0+-10
         R_prime(2,1)=0.0_dblprec;R_prime(2,2)=mom_hat(3);R_prime(2,3)=0.0_dblprec
      else  !else R_2=001 x m
         !R_prime(2,1)=-mom_hat(2)
         !R_prime(2,2)=mom_hat(3)
         !!R_prime(2,1)=mom_hat(2)
         !!R_prime(2,2)=-mom_hat(3)
         !R_prime(2,3)=0.0_dblprec
         R_prime(2,1)=mom_hat(2)*1.0_dblprec-0.0_dblprec*mom_hat(3)
         R_prime(2,2)=mom_hat(3)*0.0_dblprec-1.0_dblprec*mom_hat(1)
         R_prime(2,3)=mom_hat(1)*0.0_dblprec-0.0_dblprec*mom_hat(2)
         mnorm=sqrt(sum(R_prime(2,:)*R_prime(2,:)))
         R_prime(2,:)=R_prime(2,:)/mnorm
      end if  !R_1 = R_2 x R_3
      R_prime(1,1)=R_prime(2,2)*R_prime(3,3)-R_prime(2,3)*R_prime(3,2)
      R_prime(1,2)=R_prime(2,3)*R_prime(3,1)-R_prime(2,1)*R_prime(3,3)
      R_prime(1,3)=R_prime(2,1)*R_prime(3,2)-R_prime(2,2)*R_prime(3,1)
      mnorm=sqrt(sum(R_prime(1,:)*R_prime(1,:)))
      R_prime(1,:)=R_prime(1,:)/mnorm

      u=R_prime(1,:)+im*R_prime(2,:)
      v=R_prime(3,:)
      !
      return
      !
   end subroutine find_uv

   subroutine find_R(R,mom_0,mom_n)
      !
      implicit none
      !
      real(dblprec), dimension(3,3), intent(out) :: R
      real(dblprec), dimension(3), intent(in) :: mom_0
      real(dblprec), dimension(3), intent(in) :: mom_n
      !

      real(dblprec) :: mnorm
      real(dblprec) :: angle
      real(dblprec) :: mdot
      complex(dblprec) :: im = (0.0_dblprec,1.0_dblprec)
      real(dblprec), dimension(3,3) :: K, K2, one
      real(dblprec), dimension(3) :: mom_hat_0, mom_hat_n
      real(dblprec), dimension(3) :: k_hat
      real(dblprec), dimension(3) :: x,y
      !
      !
      ! Find the rotation matrix R that transforms mom_n to mom_0
      ! Done by finding the vector orthogonal to mom_n and mom_0 and then perform a
      ! Rodrigues rotation with the corresponding angle
      !
      x(1)=1.0_dblprec;x(2)=0.0_dblprec;x(3)=0.0_dblprec
      y(1)=0.0_dblprec;y(2)=1.0_dblprec;y(3)=0.0_dblprec
      K=0.0_dblprec
      one=0.0_dblprec
      one(1,1)=1.0_dblprec;one(2,2)=1.0_dblprec;one(3,3)=1.0_dblprec
      !
      ! Unit vector parallel to mom
      mnorm=sqrt(sum(mom_0*mom_0)+1.0d-15)
      mom_hat_0=mom_0/mnorm
      mnorm=sqrt(sum(mom_n*mom_n)+1.0d-15)
      mom_hat_n=mom_n/mnorm
      mdot=sum(mom_hat_0*mom_hat_n)
      angle=acos(mdot)

      ! Find perpendicular vector from cross-product
      !print '(a,3f12.6)', 'k_hat',k_hat
      !print '(a,f12.6)', 'angle: ',angle*180/3.1415
      !
      ! If m_0 and m_n are parallel, then use other approach
      if(1.0_dblprec-abs(mdot)<1.0e-9_dblprec) then
         ! Two additional checks: first if m_0 // x
         if(1.0_dblprec-abs(sum(mom_hat_0*x))<1.0e-9_dblprec) then
            k_hat(1)=mom_hat_0(2)*y(3)-mom_hat_0(3)*y(2)
            k_hat(2)=mom_hat_0(3)*y(1)-mom_hat_0(1)*y(3)
            k_hat(3)=mom_hat_0(1)*y(2)-mom_hat_0(2)*y(1)
         else
            k_hat(1)=mom_hat_0(2)*x(3)-mom_hat_0(3)*x(2)
            k_hat(2)=mom_hat_0(3)*x(1)-mom_hat_0(1)*x(3)
            k_hat(3)=mom_hat_0(1)*x(2)-mom_hat_0(2)*x(1)
         end if
         !!! ! Parallel mom_0, mom_n
         !!! if(abs(angle)<1.0_dblprec) then
         !!!    R=one
         !!! else
         !!!    R=-one
         !!! end if
      else
         ! Not parallel
         k_hat(1)=mom_hat_0(2)*mom_hat_n(3)-mom_hat_0(3)*mom_hat_n(2)
         k_hat(2)=mom_hat_0(3)*mom_hat_n(1)-mom_hat_0(1)*mom_hat_n(3)
         k_hat(3)=mom_hat_0(1)*mom_hat_n(2)-mom_hat_0(2)*mom_hat_n(1)
      end if
      k_hat=k_hat/sqrt(sum(k_hat*k_hat)+1.0e-12)
      K(2,1)= k_hat(3);K(1,2)=-k_hat(3)
      K(3,1)=-k_hat(2);K(1,3)= k_hat(2)
      K(2,3)=-k_hat(1);K(3,2)= k_hat(1)
      K2=matmul(K,K)
      R=one+K*sin(angle)+(1.0_dblprec-cos(angle))*K2
      !print *,'---I---'
      !print '(3f12.6)',one
      !print *,'---K---'
      !print '(3f12.6)',K
      !print *,'---K2--'
      !print '(3f12.6)',K2
      !print *,'---R---'
      !print '(3f12.6)',R
      !
      return
      !
   end subroutine find_R


   subroutine find_O(O,mom_0)
      !
      implicit none
      !
      complex(dblprec), dimension(3,3), intent(out) :: O
      real(dblprec), dimension(3), intent(in) :: mom_0
      !
      real(dblprec) :: theta, phi
      real(dblprec), dimension(3,3) :: Oz,Oy
      !
      !
      theta=acos(mom_0(1))
      phi=acos(mom_0(3))
      Oz=0.0_dblprec
      Oy=0.0_dblprec
      Oz(1,1)=cos(phi);Oz(2,2)=cos(phi);Oz(3,3)=1.0_dblprec
      Oz(2,1)=-sin(phi);Oz(1,2)=sin(phi)
      Oy(1,1)=cos(theta);Oy(2,2)=1.0_dblprec;Oy(3,3)=cos(theta)
      Oy(3,1)=-sin(theta);Oy(1,3)=sin(theta)
      O=matmul(Oz,Oy)
      !
      !
      return
      !
   end subroutine find_O

   subroutine get_M(M,mom_0)
      !
      implicit none
      !
      complex(dblprec), dimension(3,3), intent(out) :: M
      real(dblprec), intent(in) :: mom_0
      !
      complex(dblprec) :: im
      !
      im=(0.0_dblprec,1.0_dblprec)
      !
      M=0.0_dblprec
      M(1,1)=1.0_dblprec;M(1,2)=1.0_dblprec
      M(2,1)=-im ;M(1,2)=im
      M(3,3)=sqrt(2.0_dblprec/mom_0)
      M=sqrt(mom_0/2.0_dblprec)
      !
      !
      return
      !
   end subroutine get_M

  subroutine get_Hij(Hij,Jij,emom_i,emom_j,mmom_i,mmom_j,i,j)
      !
      implicit none
      !
      complex(dblprec), dimension(2,2), intent(out) :: Hij
      real(dblprec), dimension(3,3), intent(in) :: Jij
      real(dblprec), dimension(3), intent(in) :: emom_i
      real(dblprec), dimension(3), intent(in) :: emom_j
      real(dblprec), intent(in) :: mmom_i
      real(dblprec), intent(in) :: mmom_j
      integer, intent(in) :: i
      integer, intent(in) :: j
      !
      complex(dblprec), dimension(3,3) :: Jbar, M_i, M_j
      complex(dblprec), dimension(3,3) :: O_i, O_j
      complex(dblprec), dimension(3,3) :: Mprod
      !
      call find_O(O_i,emom_i)
      call find_O(O_j,emom_j)
      call get_M(M_i,mmom_i)
      call get_M(M_j,mmom_i)
      !
      Mprod=M_j
      Mprod=matmul(O_j,Mprod)
      Mprod=matmul(Jij,Mprod)
      Mprod=matmul(transpose(O_i),Mprod)
      Mprod=matmul(transpose(conjg(M_i)),Mprod)
      Hij=Mprod(1:2,1:2)
      if(i==j) then
         Hij(1,1)=Hij(1,1)-Mprod(3,3)
         Hij(2,2)=Hij(2,2)-Mprod(3,3)
      end if
      !
      !
      return
      !
     end subroutine get_Hij
   !!! !> Set up the Hamiltonian for first cell
   !!! subroutine setup_tensor_hamiltonian(N1,N2,N3,NT,NA,Natom, Mensemble, conf_num, simid, emomM, &
   !!!       mmom, max_no_neigh, nlistsize, nlist, jtens)

   !!!    use Constants
   !!!    !
   !!!    implicit none
   !!!    !
   !!!    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
   !!!    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
   !!!    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
   !!!    integer, intent(in) :: NT                                          !< Number of types of atoms
   !!!    integer, intent(in) :: NA                                          !< Number of atoms in one cell
   !!!    integer, intent(in) :: Natom     !< Number of atoms in system
   !!!    integer, intent(in) :: Mensemble !< Number of ensembles
   !!!    integer, intent(in) :: conf_num  !< number of configurations for LSF
   !!!    character(len=8), intent(in) :: simid !< Name of simulation
   !!!    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
   !!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector

   !!!    !
   !!!    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
   !!!    integer, dimension(Natom), intent(in) :: nlistsize     !< Size of neighbour list for Heisenberg exchange couplings
   !!!    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist       !< Neighbour list for Heisenberg exchange couplings
   !!!    real(dblprec), dimension(3,3,max_no_neigh,Natom), intent(in) :: jtens  !< Heisenberg exchange couplings in tensor form
   !!!    !
   !!!    integer :: i_stat, ia, ja, j, hdim, l, la
   !!!    complex(dblprec), dimension(3) :: ui, vi, uj, vj, ul, vl
   !!!    complex(dblprec) :: udot, ucdot, vdot
   !!!    complex(dblprec), dimension(NA,NA) :: dot_mat
   !!!    !
   !!!    !
   !!!    hdim=2*NA
   !!!    allocate(magham(hdim,hdim),stat=i_stat)
   !!!    call memocc(i_stat,product(shape(magham))*kind(magham),'magham','setup_finite_hamiltonian')
   !!!    magham=0.0_dblprec
   !!!    allocate(A_k(NA,NA),stat=i_stat)
   !!!    call memocc(i_stat,product(shape(A_k))*kind(A_k),'A_k','setup_finite_hamiltonian')
   !!!    allocate(B_k(NA,NA),stat=i_stat)
   !!!    call memocc(i_stat,product(shape(B_k))*kind(B_k),'B_k','setup_finite_hamiltonian')
   !!!    allocate(C_k(NA,NA),stat=i_stat)
   !!!    call memocc(i_stat,product(shape(C_k))*kind(C_k),'C_k','setup_finite_hamiltonian')

   !!!    call setup_diamag()

   !!!    !Toth and Lake
   !!!    A_k=0.0_dblprec;B_k=0.0_dblprec;C_k=0.0_dblprec;dot_mat=0.0_dblprec
   !!!    do ia=1,NA
   !!!       call find_uv(ui,vi,emomM(:,ia,1))
   !!!       do j=1,nlistsize(ia)
   !!!          ja=nlist(j,ia)
   !!!          call find_uv(uj,vj,emomM(:,ja,1))
   !!!          udot=sJs(ui,jtens(:,:,j,ia),uj)
   !!!          ucdot=sJs(ui,jtens(:,:,j,ia),conjg(uj))
   !!!          A_k(ia,ja)=0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ucdot
   !!!          B_k(ia,ja)=0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*udot
   !!!          dot_mat(ia,ja)=udot
   !!!       end do
   !!!       !
   !!!       C_k(ia,ia)=0.0_dblprec
   !!!       do l=1,nlistsize(ia)
   !!!          la=nlist(l,ia)
   !!!          call find_uv(ul,vl,emomM(:,la,1))
   !!!          vdot=sJs(vi,jtens(:,:,l,ia),vl)
   !!!          C_k(ia,ia)=C_k(ia,ia)+1.0_dblprec*mmom(la,1)*vdot
   !!!       end do
   !!!    end do

   !!!    allocate(h_k(hdim,hdim),stat=i_stat)
   !!!    call memocc(i_stat,product(shape(h_k))*kind(h_k),'h_k','setup_finite_hamiltonian')

   !!!    do ia=1,NA
   !!!       do ja=1,NA
   !!!          ! Diagonal as-is
   !!!          h_k(ia,ja)=A_k(ia,ja)-C_k(ia,ja)
   !!!          ! Off-diagonal Transpose and conjugate
   !!!          h_k(ia+NA,ja)=conjg(B_k(ja,ia))
   !!!          ! Off-diagonal as-is
   !!!          h_k(ia,ja+NA)=B_k(ia,ja)
   !!!          ! Diagonal Conjugate
   !!!          h_k(ia+NA,ja+NA)=conjg(A_k(ia,ja))-C_k(ia,ja)
   !!!       end do
   !!!       ! Add diagonal epsilon to ensure positive definiteness
   !!!       !h_k(ia,ia)=h_k(ia,ia)+1.0d-6
   !!!    end do

   !!!    !call diagonalize_quad_hamiltonian(Natom,Mensemble,NA,mmom,emomM)

   !!!    return
   !!!    !
   !!! end subroutine setup_tensor_hamiltonian

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
      return
   end function sJs

   subroutine shuffle_eig(eval,evec,ndim)
      !
      implicit none
      !
      integer, intent(in) :: ndim
      real(dblprec), dimension(ndim), intent(inout) :: eval
      complex(dblprec), dimension(ndim,ndim), intent(inout) :: evec
      !
      integer :: i, nhit, thit
      complex(dblprec), dimension(:,:), allocatable :: mtemp
      complex(dblprec), dimension(:), allocatable :: dmtemp
      real(dblprec), dimension(:), allocatable :: vtemp
      real(dblprec) :: dtemp
      !
      allocate(mtemp(ndim,ndim))
      mtemp=(0.0_dblprec,0.0_dblprec)
      allocate(dmtemp(ndim))
      allocate(vtemp(ndim))
      vtemp=0.0_dblprec
      !
      vtemp=eval
      mtemp=evec
      nhit=1
      i=1
      thit=0
      ! Stupid simple sort but small data set so ok.
      do while(nhit==1)
         nhit=0
         do i=1,ndim-1
            if(real(vtemp(i+1))>real(vtemp(i))) then
               nhit=1
               thit=thit+1
               dtemp=vtemp(i+1)
               vtemp(i+1)=vtemp(i)
               vtemp(i)=dtemp
               dmtemp=mtemp(:,i+1)
               mtemp(:,i+1)=mtemp(:,i)
               mtemp(:,i)=dmtemp
            end if
         end do
      end do
      eval=vtemp
      evec=mtemp
      !
      deallocate(mtemp)
      deallocate(dmtemp)
      deallocate(vtemp)
      !
   end subroutine shuffle_eig

   !> Set up the Hamiltonian for first cell
   subroutine setup_Jtens_q(Natom,Mensemble,NA,emomM,do_dm,q,nq,Jtens_q)

      use Constants
      use SystemData, only : coord, atype
      use AMS, only : wrap_coord_diff
      !
      implicit none
      !
      integer, intent(in) :: Natom   !< Number of atoms
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA      !< Number of types
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: nq      !< Number of q-points
      real(dblprec), dimension(3,0:nq), intent(in) :: q !< Qpoints 
      complex(dblprec), dimension(3,3,na,na,0:nq), intent(out) :: Jtens_q

      !
      real(dblprec), dimension(3) :: q_vec
      real(dblprec), dimension(3) :: dist
      real(dblprec), dimension(3) :: dmv
      real(dblprec), dimension(3) :: z
      real(dblprec), dimension(3,3) :: R_n, J_n, D_n, R_m
      complex(dblprec)  :: FTfac, im
      integer :: iq, ia, ja, j, jat, ih
      !
      !real(dblprec), dimension(3,3,na,na,0:nq) :: Jtens_q
      !
      im=(0.0_dblprec,1.0_dblprec)
      Jtens_q=0.0_dblprec
      z(1)=0.0_dblprec;z(2)=0.0_dblprec;z(3)=1.0_dblprec
      !
      !!!$omp parallel do default(shared) private(iq,q_vec,ia,ja,j,jat,FTfac,dmv)
      do iq=0,nq
         ! Ensure that iq=0 corresponds to gamma
         if(iq>0) then
            q_vec=q(:,iq)*2.0_dblprec*pi
         else
            q_vec=0.0_dblprec
         end if
         !print '(2x,a,i5,3f12.6)' , 'q=',iq, q_vec

         do ia=1,NA
            ih=ham%aHam(ia)
            ! Jij exchange
            do j=1,ham%nlistsize(ih)
               ja=ham%nlist(j,ia)
               jat=atype(ja)
               call wrap_coord_diff(Natom,coord,ia,ja,dist)
               !print '(a,3f12.8)',' |DIST--> ',dist

               !call find_R(R_n,z,emomM(:,ja,1))
               !call find_R(R_m,emomM(:,jat,1),z)
               call find_R(R_n,emomM(:,jat,1),emomM(:,ja,1))
               !call find_R(R_n,emomM(:,ia,1),emomM(:,ja/natom+1,1))
               !print *,'  -----R_n----- ',ia,jat
               !print '(3f10.4)', R_n
               !print *,' --->',ia,j,jat
               J_n=0.0_dblprec
               J_n(1,1)=ham%ncoup(j,ih,1)
               J_n(2,2)=ham%ncoup(j,ih,1)
               J_n(3,3)=ham%ncoup(j,ih,1)
               !print *,'  -----J_n----- ',ia,jat
               !print '(3f10.4)', matmul(J_n,R_n)

              ! print '(a,3f10.4)', 'q_vec:',q_vec(1:3)
              ! print '(a,3f10.4)', 'dist:',dist(1),dist(2),dist(3)
              ! print '(a,3f10.4)', 'q dot d',q_vec(1)*dist(1),q_vec(2)*dist(2),q_vec(3)*dist(3)
              ! print '(a,3f10.4)', 'q dot d',q_vec(3),dist(3),q_vec(3)*dist(3)
               FTfac=exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))
            !print *,' ..Jtens_q... ', FTfac

               Jtens_q(:,:,ia,jat,iq)=Jtens_q(:,:,ia,jat,iq)+matmul(J_n,R_n)*FTfac
               !Jtens_q(:,:,ia,jat,iq)=Jtens_q(:,:,ia,jat,iq)+J_n*FTfac
            !print '(6f10.4)', Jtens_q(:,:,ia,jat,iq)
            end do
            ! DMI
            if(do_dm==1) then
               do j=1,ham%dmlistsize(ih)
                  ja=ham%dmlist(j,ia)
                  jat=atype(ja)
                  dmv=ham%dm_vect(:,j,ih)

                  call find_R(R_n,emomM(:,jat,1),emomM(:,ja,1))
                  D_n=dm2tens(dmv)

                  FTfac=exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))
                  Jtens_q(:,:,ia,jat,iq)=Jtens_q(:,:,ia,jat,iq)+matmul(D_n,R_n)*FTfac
               end do
            end if
            !print *,'--Jtens-final--',iq
            !print '(9f10.4)', Jtens_q(:,:,ia,1,iq)
            !print *,'-            -'
            !print *,'-            -'
            !print *,'-            -'
            !print *,'-            -'
         end do


         !!! print *,'--Jtens-final-real-',iq
         !!! do ia=1,NA
         !!!    do ja=1,NA
         !!!       print *,ia,ja
         !!!       print '(3f10.4)', real(Jtens_q(:,:,ia,ja,iq))
         !!!    end do
         !!! end do

         !!! print *,'--Jtens-final-complex-',iq
         !!! do ia=1,NA
         !!!    do ja=1,NA
         !!!       print *,ia,ja
         !!!       print '(3f10.4)', aimag(Jtens_q(:,:,ia,ja,iq))
         !!!    end do
         !!! end do

         do ia=1,NA
            do ja=1,NA
               Jtens_q(:,:,ia,ja,iq)=0.5_dblprec*(Jtens_q(:,:,ia,ja,iq)+transpose(conjg(Jtens_q(:,:,ja,ia,iq))))
               Jtens_q(:,:,ja,ia,iq)=transpose(conjg(Jtens_q(:,:,ia,ja,iq)))
            end do
         end do


      end do
      !!!$omp end parallel do

      return
      !
   end subroutine setup_Jtens_q

   function dm2tens(dmvec) result(dtens)
      !
      implicit none
      !
      real(dblprec), dimension(3), intent(in) :: dmvec !< a DM vector
      !
      real(dblprec), dimension(3,3) :: dtens
      !

      dtens(1,1)= 0.0_dblprec
      dtens(2,1)= dmvec(3)
      dtens(3,1)=-dmvec(2)
      dtens(1,2)=-dmvec(3)
      dtens(2,2)= 0.0_dblprec
      dtens(3,2)= dmvec(1)
      dtens(1,3)=-dmvec(2)
      dtens(2,3)= dmvec(1)
      dtens(3,3)= 0.0_dblprec

      return
      !
   end function dm2tens

   subroutine clone_q(nq,q,q_ext)
      !
      implicit none
      !
      integer, intent(in) :: nq !< No. qpoints from input
      real(dblprec), dimension(3,nq), intent(in) :: q !< qpoints from input
      real(dblprec), dimension(3,0:2*nq), intent(out) :: q_ext !< Extended set of q-points
      !
      integer :: iq
      !

      q_ext(:,0)=0.0_dblprec
      do iq=1,nq
         q_ext(:,iq)=q(:,iq)
         !print '(a,i4,3f10.6)', 'clone_q:',iq,q(:,iq)
         q_ext(:,iq+nq)=-q(:,iq)
      end do

      return
   end subroutine clone_q

   !---------------------------------------------------------------------------
   subroutine read_parameters_diamag(ifile)
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
            ! START OF VARIABLES FOR FREQ SPIN CORRELATION
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            case('do_diamag') ! Perform frequency based spin-correlation sampling
               read(ifile,*,iostat=i_err) do_diamag
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! END OF LEGACY VARIABLES
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            case default
               if(len(trim(keyword))>0) then
                  read(ifile,*)
               end if

            end select
         end if

         ! End of file
         if (i_errb==20) goto 20
         ! End of row
         if (i_errb==10) goto 10
      end do

      20  continue

      return
   end subroutine read_parameters_diamag 

   !> Set up the Hamiltonian for first cell
   subroutine setup_altern_hamiltonian(N1,N2,N3,NT,NA,Natom, nHam, Mensemble, conf_num, simid, emomM, &
         mmom,do_dm)

      use Constants
      use SystemData, only : coord, atype
      use AMS, only : wrap_coord_diff
      use Correlation,        only : q,nq
      !
      implicit none
      !
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: NT                                          !< Number of types of atoms
      integer, intent(in) :: NA                                          !< Number of atoms in one cell
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: nHam      !< Number of atoms in Hamiltonian
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)

      !
      real(dblprec), dimension(3) :: q_vec, q_plus, q_minus, q_0
      real(dblprec), dimension(3) :: dist, z
      real(dblprec), dimension(3) :: dmv
      real(dblprec), dimension(3,3) :: J_n, D_n, R_n, R_m
      real(dblprec), dimension(3,3) :: O_i, O_j
      complex(dblprec), dimension(3,3) :: J_prime
      complex(dblprec), dimension(3,3) :: M_i, M_j
      complex(dblprec), dimension(2,2) :: H_ij
      !
      integer :: i_stat, ia, ja, j, hdim, l, la, iq, jat, lat
      complex(dblprec), dimension(3) :: ui, vi, uj, vj, ul, vl
      complex(dblprec) :: udot, ucdot, vdot, im
      complex(dblprec), dimension(NA,NA) :: dot_mat
      !
      complex(dblprec), dimension(:,:), allocatable :: eig_vec
      real(dblprec), dimension(:), allocatable :: eig_val
      real(dblprec), dimension(:,:), allocatable :: q_ext
      !
      hdim=2*NA
      allocate(eig_vec(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_vec))*kind(eig_vec),'eig_vec','setup_infinite_hamiltonian')
      allocate(eig_val(hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_val))*kind(eig_val),'eig_val','setup_infinite_hamiltonian')
      allocate(h_k(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(h_k))*kind(h_k),'h_k','setup_infinite_hamiltonian')

      call setup_diamag()

      !call clone_q(nq,q,q_ext)


      im=(0.0_dblprec,1.0_dblprec)


      !Santos, Santos, and Santos
      do iq=1,nq
         q_vec=q(:,iq)*2.0_dblprec*pi
         !q_plus=q(:,iq)*2.0_dblprec*pi
         !q_minus=-q(:,iq)*2.0_dblprec*pi
         !q_0=0.0_dblprec
         !print '(2x,a,3f12.6)' , 'q=',q_vec

         h_k=0.0_dblprec
         do ia=1,NA
            
            ! Jij exchange
            do j=1,ham%nlistsize(ia)
               ja=ham%nlist(j,ia)
               jat=atype(ja)
               
               call wrap_coord_diff(Natom,coord,ia,ja,dist)


               J_n=0.0_dblprec; J_n(1,1)=ham%ncoup(j,ia,1); J_n(2,2)=ham%ncoup(j,ia,1); J_n(3,3)=ham%ncoup(j,ia,1)

               call get_Hij(H_ij,J_n,emomM(:,ia,1),emomM(:,ja,1),mmom(ia,1),mmom(ja,1),ia,jat)

               J_prime=J_n

               h_k(ia,jat)      =h_k(ia,jat)      +H_ij(1,1)*exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))
               h_k(ia+NA,jat)   =h_k(ia+NA,jat)   +H_ij(2,1)*exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))
               h_k(ia,jat+NA)   =h_k(ia,jat+NA)   +H_ij(1,2)*exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))
               h_k(ia+NA,jat+NA)=h_k(ia+NA,jat+NA)+H_ij(2,2)*exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))

            end do
         end do

         call diagonalize_altern_hamiltonian(NA,h_k,eig_val,eig_vec)

         write(1000,'(i6,12f12.6)') iq,real(eig_val)*13.605698066_dblprec*4.0_dblprec
         !print '(a,i6,12f12.6)','ev: ',iq,real(eig_val)*13.605698066_dblprec*4.0_dblprec

      end do

      return
      !
   end subroutine setup_altern_hamiltonian

   subroutine diagonalize_altern_hamiltonian(NA,h_in,eig_val,eig_vec)
      !
      use Constants
      !
      implicit none
      !
      integer, intent(in) :: NA
      complex(dblprec), dimension(2*NA,2*NA) :: h_in
      complex(dblprec), dimension(2*NA,2*NA) :: eig_vec
      real(dblprec), dimension(2*NA) :: eig_val
      !
      !
      complex(dblprec), dimension(:,:), allocatable :: g_mat
      complex(dblprec), dimension(:,:), allocatable :: U_mat
      complex(dblprec), dimension(:,:), allocatable :: D_mat
      complex(dblprec), dimension(:,:), allocatable :: H_mat
      complex(dblprec), dimension(:,:), allocatable :: L_eig
      complex(dblprec), dimension(:,:), allocatable :: R_eig
      complex(dblprec), dimension(:), allocatable :: cwork
      real(dblprec), dimension(:), allocatable :: rwork
      !
      integer :: info, lwork, hdim, i_stat, ia, ja

      complex(dblprec) :: cone, czero, fcinv, dia_eps, im

      complex(dblprec), dimension(3) :: ul, vl
      !
      !
      hdim=2*NA
      czero=(0.0_dblprec,0.0_dblprec)
      cone=(1.0_dblprec,0.0_dblprec)
      im=(0.0_dblprec,1.0_dblprec)
      dia_eps=1.0d-12
      !
      fcinv = 0.5_dblprec*mub/mry
      !
      allocate(g_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(g_mat))*kind(g_mat),'g_mat','diagonalize_altern_hamiltonian')
      allocate(R_eig(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(R_eig))*kind(R_eig),'R_eig','diagonalize_altern_hamiltonian')
      allocate(L_eig(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(L_eig))*kind(L_eig),'K_mat','diagonalize_altern_hamiltonian')
      allocate(D_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(D_mat))*kind(D_mat),'D_mat','diagonalize_altern_hamiltonian')
      allocate(H_mat(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(H_mat))*kind(H_mat),'H_mat','diagonalize_altern_hamiltonian')
      !
      g_mat=(0.0_dblprec,0.0_dblprec)
      do ia=1,NA
         g_mat(ia,ia)=(-1.0_dblprec,0.0_dblprec)
         g_mat(ia+NA,ia+NA)=(1.0_dblprec,0.0_dblprec)
      end do
      !
      H_mat=matmul(g_mat,h_in)
      D_mat=H_mat
      print *,'---real h_in---------'
      print '(6f22.16)',real(h_in)
      print *,'---imag h_in---------'
      print '(6f22.16)',aimag(h_in)
      !
      ! General eigenvalue problem

      lwork=6*hdim-1
      allocate(cwork(lwork))
      allocate(rwork(2*hdim))
      call zgeev('V','V',hdim,D_mat,hdim,eig_val,L_eig,hdim,R_eig,hdim,cwork,lwork,rwork,info)
      deallocate(cwork)
      deallocate(rwork)



      do ia=1,hdim
         write(2000   , '(100f10.4)')  eig_val(ia) !*13.605698066_dblprec*4.0_dblprec
      end do
      !close(ofileno)
      !do ia=1,NA
      !   do ja=1,NA
      !      call find_uv(ul,vl,emomM(:,ja,1))
      !      bigS(:,ia,ja)=sqrt(mmom(ia,1))*sqrt(0.5_dblprec)*(conjg(ul)*x_mat(ia,ja)+ul*conjg(x_mat(ia,ja)))+vl*(mmom(ia,1)-conjg(x_mat(ia,ja))*x_mat(ia,ja))
      !   end do
      !end do
      return
      !
   end subroutine diagonalize_altern_hamiltonian



end module diamag
