!> Data and routines for calculate frequency based spin correlation function S(w)
!> @authors
!> Anders Bergman, Manuel Pereiro
!> @copyright
!> GNU Public License
!!
module diamag_variants
   use Parameters
   use Profiling
   use Hamiltoniandata,    only : ham
   use diamag
   !
   implicit none
   !
   !
   private
   ! public subroutines
   public ::  setup_finite_hamiltonian, setup_infinite_hamiltonian
   public ::  setup_altern_hamiltonian

contains

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
      integer :: i_stat, ia, ja, j, hdim, l, la, iq, jat, lat, iv
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


end module diamag_variants
