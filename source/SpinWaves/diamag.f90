!SpinWaves/diamag.f90!> Data and routines for calculate frequency based spin correlation function S(w)
!> @authors
!> Anders Bergman, Manuel Pereiro
!> @copyright
!> GNU Public License
!!
module diamag
   use Parameters
   use Profiling
   use Hamiltoniandata,    only : ham
   use InputData,   only : ham_inp
   !
   implicit none
   !
   !
   real(dblprec), dimension(:,:), allocatable :: magham !< M(w)
   complex(dblprec), dimension(:,:), allocatable :: A_k
   complex(dblprec), dimension(:,:), allocatable :: A_km
   complex(dblprec), dimension(:,:), allocatable :: B_k
   complex(dblprec), dimension(:,:), allocatable :: C_k
   complex(dblprec), dimension(:,:), allocatable :: h_k
   complex(dblprec), dimension(:,:,:,:,:), allocatable :: S_prime
   complex(dblprec), dimension(:,:,:), allocatable :: ektij
   !
   real(dblprec), dimension(:,:), allocatable :: nc_eval_q   !< Eigenvalues from NC-AMS
   real(dblprec), dimension(:,:), allocatable :: nc_eval_qchern   !< Eigenvalues from NC-AMS
   complex(dblprec),  dimension(:,:,:), allocatable :: nc_evec_qchern   !< Eigenvalues from NC-AMS in Chern
   real(dblprec), dimension(:,:,:), allocatable :: nc_evec_q   !< Eigenvalues from NC-AMS
   !
   character(len=1) :: do_diamag       !< Perform frequency based spin-correlation sampling (Y/N/C)
   real(dblprec)    :: diamag_mix      !< Separation between sampling steps
   real(dblprec)    :: diamag_thresh   !< Separation between sampling steps
   real(dblprec)    :: diamag_eps      !< Diagonal offset for positive definitive diagonalization
   integer          :: diamag_niter    !< Number of steps to sample
   integer          :: diamag_nfreq    !< Number of frequencies w in S(k,w)
   real(dblprec), dimension(3)    :: diamag_qvect=0.0_dblprec !< Spin spiral ordering vector
   real(dblprec), dimension(3)    :: diamag_nvect !< Spin spiral pitch vector

   !


   private
   ! public subroutines
   public :: do_diamag, read_parameters_diamag,clone_q,diagonalize_quad_hamiltonian,&
             find_uv,setup_ektij,setup_jtens2_q,setup_jtens_q,sJs
   public :: setup_tensor_hamiltonian, nc_evec_qchern, nc_eval_qchern
   public :: diamag_qvect, nc_eval_q, nc_evec_q

contains

   subroutine setup_diamag()

      implicit none

      do_diamag='N'
      diamag_niter=1000
      diamag_mix=0.030_dblprec
      diamag_thresh=1.0d-8
      diamag_nfreq=200
      !diamag_eps=-1.0_dblprec
      !diamag_qvect=0.0_dblprec
      !diamag_nvect(1:2)=0.0_dblprec;diamag_nvect(3)=1.0_dblprec

   end subroutine setup_diamag

   !> Set up the Hamiltonian for first cell
   subroutine setup_tensor_hamiltonian(NA,Natom, Mensemble, simid, emomM, mmom, q_vect, nq,flag)

      use Constants
      use AMS, only : magdos_calc, printEnergies
      !
      implicit none
      !
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: flag !< Type of calculation (0=NC-AMS, 1=Chern number)
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      real(dblprec), dimension(3,nq), intent(in) :: q_vect !< Qpoints
      integer, intent(in) :: nq    !< Number of qpoints
      !
      real(dblprec), dimension(3) :: z
      !
      integer :: i_stat, ia, ja, hdim, iq, nq_ext
      integer :: alfa, beta, i_all
      complex(dblprec), dimension(3) :: ui, vi, uj, vj
      complex(dblprec) :: im
      complex(dblprec), dimension(NA,NA) :: dot_mat
      !
      complex(dblprec), dimension(:,:), allocatable :: eig_vec
      real(dblprec), dimension(:), allocatable :: eig_val
      complex(dblprec), dimension(:,:,:,:,:), allocatable ::   jtens_q  !< FT of Exchange tensor
      real(dblprec), dimension(:,:), allocatable :: q_ext
      !
      real(dblprec) :: msat,tcmfa,tcrpa
      !
      character(LEN = 25) :: ncams_file
      !
      if (flag == 0) then
        print '(1x,a)', 'Calculating LSWT magnon dispersions'
      end if

      ! Hamiltonian dimension = 4x number of atoms
      hdim=2*NA
      ! Extended q-mesh to capture -q and phasons
      nq_ext=nq*6
      allocate(magham(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(magham))*kind(magham),'magham','setup_tensor_hamiltonian')
      magham=0.0_dblprec
      allocate(q_ext(3,0:nq_ext),stat=i_stat)
      call memocc(i_stat,product(shape(q_ext))*kind(q_ext),'q_ext','setup_tensor_hamiltonian')
      allocate(jtens_q(3,3,NA,NA,0:nq_ext),stat=i_stat)
      call memocc(i_stat,product(shape(jtens_q))*kind(jtens_q),'jtens_q','setup_tensor_hamiltonian')
      !
      !
      allocate(ektij(NA,NA,0:nq_ext),stat=i_stat)
      call memocc(i_stat,product(shape(ektij))*kind(ektij),'ektij','setup_tensor_hamiltonian')
      allocate(S_prime(hdim,hdim,3,3,nq_ext),stat=i_stat)
      call memocc(i_stat,product(shape(S_prime))*kind(S_prime),'S_prime','setup_tensor_hamiltonian')
      S_prime=0.0_dblprec
      !
      allocate(A_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_k))*kind(A_k),'A_k','setup_tensor_hamiltonian')
      allocate(A_km(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_km))*kind(A_km),'A_km','setup_tensor_hamiltonian')
      allocate(B_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(B_k))*kind(B_k),'B_k','setup_tensor_hamiltonian')
      allocate(C_k(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(C_k))*kind(C_k),'C_k','setup_tensor_hamiltonian')
      !
      allocate(h_k(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(h_k))*kind(h_k),'h_k','setup_tensor_hamiltonian')
      if (flag == 0) then
        allocate(nc_eval_q(hdim,nq_ext),stat=i_stat)
        call memocc(i_stat,product(shape(nc_eval_q))*kind(nc_eval_q),'nc_eval_q','setup_tensor_hamiltonian')
        allocate(nc_evec_q(hdim,hdim,nq_ext),stat=i_stat)
        call memocc(i_stat,product(shape(nc_evec_q))*kind(nc_evec_q),'nc_evec_q','setup_tensor_hamiltonian')
      else
        allocate(nc_eval_qchern(hdim,nq_ext),stat=i_stat)
        call memocc(i_stat,product(shape(nc_eval_qchern))*kind(nc_eval_qchern),'nc_eval_qchern','setup_tensor_hamiltonian')
        allocate(nc_evec_qchern(hdim,hdim,nq_ext),stat=i_stat)
        call memocc(i_stat,product(shape(nc_evec_qchern))*kind(nc_evec_qchern),'nc_evec_qchern','setup_tensor_hamiltonian')
      end if
         !
      allocate(eig_vec(hdim,hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_vec))*kind(eig_vec),'eig_vec','setup_tensor_hamiltonian')
      allocate(eig_val(hdim),stat=i_stat)
      call memocc(i_stat,product(shape(eig_val))*kind(eig_val),'eig_val','setup_tensor_hamiltonian')

      !call setup_diamag()

      im=(0.0_dblprec,1.0_dblprec)

      call clone_q(nq,q_vect,nq_ext,q_ext,diamag_qvect)

      z(1)=0.0_dblprec;z(2)=0.0_dblprec;z(3)=1.0_dblprec

      ! Create J(q) tensort
      if (ham_inp%do_jtensor==1) then
         call setup_Jtens2_q(Natom,Mensemble,NA,emomM,q_ext,nq_ext,Jtens_q)
      else
         call setup_Jtens_q(Natom,Mensemble,NA,emomM,q_ext,nq_ext,Jtens_q)
      end if

      ! Create array of exp(i*k*(ri-rj))
      call setup_ektij(Natom,NA,q_ext,nq_ext,ektij)

      !Toth and Lake loop over k-vectors 
      !print *,nq,nq_ext
      do iq=1,3*nq
         !do iq=1,nq_ext

         ! Sub-matrices A(k),B(k),C
         A_k=0.0_dblprec;A_km=0.0_dblprec;B_k=0.0_dblprec;C_k=0.0_dblprec;dot_mat=0.0_dblprec
         do ia=1,NA
            !call find_uv(ui,vi,z)
            call find_uv(ui,vi,emomM(:,ia,1))
            do ja=1,NA
               !call find_uv(uj,vj,z)
               call find_uv(uj,vj,emomM(:,ja,1))

               ! iq+3*nq -> negative q
               A_k(ia,ja) =A_k(ia,ja) +0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,Jtens_q(:,:,ia,ja,iq+3*nq),conjg(uj))
               A_km(ia,ja)=A_km(ia,ja)+0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,Jtens_q(:,:,ia,ja,iq),conjg(uj))
               B_k(ia,ja) =B_k(ia,ja) +0.5_dblprec*sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*sJs(ui,Jtens_q(:,:,ia,ja,iq+3*nq),uj)
               C_k(ia,ia) =C_k(ia,ia) +1.0_dblprec*mmom(ja,1)*sJs(vi,Jtens_q(:,:,ia,ja,0),vj)

               ! Also mount S´ for later correlation function use 
               do alfa=1,3
                  do beta=1,3
                     ! Y(alfa,beta,k)
                     S_prime(ia,ja,alfa,beta,iq)= &
                        sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ui(alfa)*conjg(uj(beta))*ektij(ia,ja,iq)
                     ! Z(alfa,beta,k)
                     S_prime(ia,ja+NA,alfa,beta,iq)= &
                        sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*ui(alfa)*uj(beta)*ektij(ia,ja,iq)
                     ! V(alfa,beta,k)
                     S_prime(ia+NA,ja,alfa,beta,iq)= &
                        sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*conjg(ui(alfa))*conjg(uj(beta))*ektij(ia,ja,iq)
                     ! W(alfa,beta,k)
                     S_prime(ia+NA,ja+NA,alfa,beta,iq)= &
                        sqrt(mmom(ia,1))*sqrt(mmom(ja,1))*conjg(ui(alfa))*uj(beta)*ektij(ia,ja,iq)
                  end do
               end do

            end do
         end do

         ! Hamiltonian h = [A(k)-C B(k) ; B´(k) A(-k)+C ]
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
            ! Add diagonal epsilon to ensure positive definiteness
            !h_k(ia,ia)=h_k(ia,ia)+1.0d-5
         end do

         !!! write(2001,*) iq
         !!! write(2001,*) 'A'
         !!! write(2001,'(12f12.6)')A_k
         !!! write(2001,*) 'B'
         !!! write(2001,'(12f12.6)')B_k
         !!! write(2001,*)'C'
         !!! write(2001,'(12f12.6)')C_k
         !!! write(2001,*)'h'
         !!! write(2001,'(24f12.6)')h_k

         !if (iq<=nq) then 
         !   eig_fcinv = 0.5_dblprec*mub/mry
         !   eig_A=A_k-C_k
         !   eig_lwork=NA*NA !hdim*hdim
         !   call zheev('V','U',NA,eig_A,NA,eig_W, eig_work, eig_lwork, eig_rwork, eig_info)
         !   !call zheev('V','U',hdim,eig_A,hdim,eig_W, eig_work, eig_lwork, eig_rwork, eig_info)

         !   !e_rot=f_cross_product(emomM(:,ia,1),q(:,iq))
         !   !e_norm = norm2(e_rot)
         !   !if(e_norm>0.0_dblprec) e_rot=e_rot/e_norm

         !   write(4200, '(i5,3x,3f16.8,4x,12f20.12)'),iq,q(:,iq),minval(-eig_W)*eig_fcinv*ry_ev*4.0_dblprec*2.0_dblprec
         !   !write(4200, '(i5,3x,3f16.8,4x,3f16.8,12f20.12)'),iq,q(:,iq),e_rot,minval(-eig_W)*eig_fcinv*ry_ev*4.0_dblprec*2.0_dblprec
         !end if

         ! Diagonalize Hamiltonian 
         call diagonalize_quad_hamiltonian(NA,h_k,eig_val,eig_vec,iq,nq_ext,S_prime)

         ! Store eigenvalues and vectors (eigenvalues in meV)
         if (flag==0) then
           nc_eval_q(:,iq)=real(eig_val)*ry_ev*4.0_dblprec
           nc_evec_q(:,:,iq)=abs(eig_vec)
         else
           nc_eval_qchern(:,iq)=real(eig_val)*ry_ev*4.0_dblprec
           nc_evec_qchern(:,:,iq)=eig_vec
         end if
      end do

      if (flag==0) then
        ncams_file = 'ncams.'//trim(simid)//'.out'
        tcmfa=0.0_dblprec
        tcrpa=0.0_dblprec
        msat=sum(mmom(:,1),1)/natom
        call printEnergies(ncams_file,nc_eval_q(1:NA,1:nq),msat,tcmfa,tcrpa,NA,2)
        !if(norm2(diamag_qvect)>1.0e-12_dblprec) then
        ncams_file = 'ncams+q.'//trim(simid)//'.out'
        call printEnergies(ncams_file,nc_eval_q(1:NA,nq+1:2*nq),msat,tcmfa,tcrpa,NA,2)
        ncams_file = 'ncams-q.'//trim(simid)//'.out'
        call printEnergies(ncams_file,nc_eval_q(1:NA,2*nq+1:3*nq),msat,tcmfa,tcrpa,NA,2)
        !end if

        !!! write (filn,'(''ncams_evec.'',a,''.out'')') trim(simid)
        !!! open(ofileno,file=filn, position='append')
        !!! do iq=1,nq
        !!!    do iv=1,hdim
        !!!       write(ofileno,'(2i6,100f18.8)') iq,iv,nc_evec_q(:,iv,iq)
        !!!    end do
        !!! end do
        !!! close(ofileno)
        !   subroutine magdos_calc(filename,wres,magdos,iat)
        !if (do_magdos=='Y') then
        !   magdos_file = 'ncmagdos.'//simid//'.out'
        !   call magdos_calc(magdos_file,nc_eval_q(1:na,1:nq),na)
        !   !call magdos_calc(magdos_file,nc_eval_q(1:na,:),nc_magdos,na,nq)
        !end if


        call diamag_sqw(NA,nq,nq_ext,diamag_nfreq,q_vect,nc_eval_q,S_prime,simid,diamag_nvect)

      end if
      
      i_all=-product(shape(magham))*kind(magham)
      deallocate(magham,stat=i_stat)
      call memocc(i_stat,i_all,'magham','setup_tensor_hamiltonian')
      i_all=-product(shape(q_ext))*kind(q_ext)
      deallocate(q_ext,stat=i_stat)
      call memocc(i_stat,i_all,'q_ext','setup_tensor_hamiltonian')

      i_all=-product(shape(jtens_q))*kind(jtens_q)
      deallocate(jtens_q,stat=i_stat)
      call memocc(i_stat,i_all,'jtens_q','setup_tensor_hamiltonian')

      i_all=-product(shape(ektij))*kind(ektij)
      deallocate(ektij,stat=i_stat)
      call memocc(i_stat,i_all,'ektij','setup_tensor_hamiltonian')
      i_all=-product(shape(S_prime))*kind(S_prime)
      deallocate(S_prime,stat=i_stat)
      call memocc(i_stat,i_all,'S_prime','setup_tensor_hamiltonian')

      i_all=-product(shape(A_k))*kind(A_k)
      deallocate(A_k,stat=i_stat)
      call memocc(i_stat,i_all,'A_k','setup_tensor_hamiltonian')
      i_all=-product(shape(A_km))*kind(A_km)
      deallocate(A_km,stat=i_stat)
      call memocc(i_stat,i_all,'A_km','setup_tensor_hamiltonian')
      i_all=-product(shape(B_k))*kind(B_k)
      deallocate(B_k,stat=i_stat)
      call memocc(i_stat,i_all,'B_k','setup_tensor_hamiltonian')
      i_all=-product(shape(C_k))*kind(C_k)
      deallocate(C_k,stat=i_stat)
      call memocc(i_stat,i_all,'C_k','setup_tensor_hamiltonian')
      i_all=-product(shape(h_k))*kind(h_k)
      deallocate(h_k,stat=i_stat)
      call memocc(i_stat,i_all,'h_k','setup_tensor_hamiltonian')
      !
      if (flag == 0) then
      i_all=-product(shape(nc_eval_q))*kind(nc_eval_q)
      deallocate(nc_eval_q,stat=i_stat)
      call memocc(i_stat,i_all,'nc_eval_q','setup_tensor_hamiltonian')
      i_all=-product(shape(nc_evec_q))*kind(nc_evec_q)
      deallocate(nc_evec_q,stat=i_stat)
      call memocc(i_stat,i_all,'nc_evec_q','setup_tensor_hamiltonian')
      end if
      !
      i_all=-product(shape(eig_vec))*kind(eig_vec)
      deallocate(eig_vec,stat=i_stat)
      call memocc(i_stat,i_all,'eig_vec','setup_tensor_hamiltonian')
      i_all=-product(shape(eig_val))*kind(eig_val)
      deallocate(eig_val,stat=i_stat)
      call memocc(i_stat,i_all,'eig_val','setup_tensor_hamiltonian')

      if (flag == 0) then
      print '(1x,a)', 'LSWT done.'
      end if
      

      return
      !
   end subroutine setup_tensor_hamiltonian


   subroutine diagonalize_quad_hamiltonian(NA,h_in,eig_val,eig_vec,iq,nq_ext,S_prime)
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
      integer, intent(in)  :: iq    !< Current q-point index
      integer, intent(in)  :: nq_ext    !< Number of qpoints
      complex(dblprec), dimension(2*NA,2*NA,3,3,nq_ext), intent(inout) :: S_prime
      !
      !real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
      !integer, intent(in) :: Mensemble !< Number of ensembles
      !real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      !
      !
      complex(dblprec), dimension(2*NA,2*NA) :: g_mat
      complex(dblprec), dimension(2*NA,2*NA) :: E_mat
      complex(dblprec), dimension(2*NA,2*NA) :: sqE_mat
      complex(dblprec), dimension(2*NA,2*NA) :: L_mat
      complex(dblprec), dimension(2*NA,2*NA) :: T_mat
      complex(dblprec), dimension(2*NA,2*NA) :: K_mat
      complex(dblprec), dimension(2*NA,2*NA) :: KgK_mat
      complex(dblprec), dimension(2*NA,2*NA) :: iK_mat
      complex(dblprec), dimension(2*NA,2*NA) :: dum_mat
      complex(dblprec), dimension(2*NA,2*NA) :: x_mat
      complex(dblprec), dimension(:), allocatable :: cwork
      real(dblprec), dimension(:), allocatable :: rwork
      !
      integer :: info, lwork, hdim, ia, ja
      integer :: alfa, beta

      complex(dblprec) :: cone, czero, fcinv, im, dia_eps

      !  print *,'In diagonalize quad Hamiltonian'
      !  print '(12f10.6)',real(h_in)
      !
      !
      hdim=2*NA
      czero=(0.0_dblprec,0.0_dblprec)
      cone=(1.0_dblprec,0.0_dblprec)
      im=(0.0_dblprec,1.0_dblprec)
      ! Add offset to ensure positive definiteness, if needed.
      ! Can be controlled by input parameter `nc_eps`
      if (diamag_eps>-1.0_dblprec) then
         dia_eps=diamag_eps
      else
         dia_eps=1.0d-6
      end if
      !
      fcinv = 0.5_dblprec*mub/mry
      !
      T_mat=(0.0_dblprec,0.0_dblprec)
      !
      g_mat=(0.0_dblprec,0.0_dblprec)
      do ia=1,NA
         g_mat(ia,ia)=(1.0_dblprec,0.0_dblprec)
         g_mat(ia+NA,ia+NA)=(-1.0_dblprec,0.0_dblprec)
      end do
      !
      ! Cholesky decomposition of h to K´*K
      !!! print *,'h'
      !!! print '(8f12.6)',h_in
      !!! print *,'Re h'
      !!! print '(4f12.6)',real(h_in)
      !!! print *,'Im h'
      !!! print '(4f12.6)',aimag(h_in)
      K_mat=-h_in*fcinv
      lwork=4*NA-1
      allocate(cwork(lwork))
      allocate(rwork(6*NA-2))
      eig_vec=K_mat
      ! Pre-diagonalization to estimate eigenvalues (not used for final eigenvalues)
      call zheev('V','U',hdim,eig_vec,hdim,eig_val, cwork, lwork, rwork, info)
      !dia_eps=0.0_dblprec
      !dia_eps=1.0e-6_dblprec
      !if(minval(eig_val)<0.0_dblprec) print *,'zheev',info,dia_eps,minval(eig_val)
      !if(minval(eig_val)<0.0_dblprec) print *,'zheev',real(K_mat(1,1)-K_mat(1,2))
      !if(minval(eig_val)<0.0_dblprec) print '(2f10.6)',real(K_mat)
      !if(minval(eig_val)<0.0_dblprec) print *,'-----------------'
      !if(minval(eig_val)<0.0_dblprec) print '(2f10.6)',aimag(K_mat)
      !if(minval(eig_val)<0.0_dblprec) print *,'-----------------'
      if(minval(eig_val)<0.0_dblprec) dia_eps=dia_eps-minval(eig_val)*1.0_dblprec
      deallocate(cwork)
      deallocate(rwork)

      ! Add eps to diagonal to ensure positive definiteness
      do ia=1,hdim
         !K_mat(ia,ia)=K_mat(ia,ia)+dia_eps
         K_mat(ia,ia)=K_mat(ia,ia)+2.0_dblprec*dia_eps
      end do

      ! Cholesky
      call zpotrf('U',hdim,K_mat,hdim,info)
      ! Setting to zero the lower triangular part of K_mat
      do ia=1,hdim
         do ja=1,hdim
            if ( ia .gt. ja) then
               K_mat(ia,ja) = czero
            end if
         enddo
      enddo

      if(info==0) then  ! Positive-definit matrix, Colpa diagonalization ok
         do ia=1,hdim
            do ja=ia+1,hdim
               K_mat(ja,ia)=0.0_dblprec
            end do
         end do
         call zgemm('N','C',hdim,hdim,hdim,cone,g_mat,hdim,K_mat,hdim,czero,dum_mat,hdim)
         call zgemm('N','N',hdim,hdim,hdim,cone,K_mat,hdim,dum_mat,hdim,czero,eig_vec,hdim)
      else
         print *,' Warning in diamag: non-positive definite matrix in zpotrf', iq, info
         ! print '(12f10.6)',real(h_in)
         call fallback_bosonic_diag(NA, h_in, eig_val, eig_vec, iq)
         ! print *,'eig_val',eig_val
      end if

      allocate(cwork(lwork))
      allocate(rwork(6*NA-2))
      ! Symmetrization of K (from SpinW code, not TothLake)
      ! do ia=1,hdim
      !    do ja=1,hdim
      !       eig_vec(ia,ja)=0.5_dblprec*(eig_vec(ia,ja)+conjg(eig_vec(ja,ia)))
      !       eig_vec(ja,ia)=eig_vec(ia,ja)
      !    end do
      ! end do
      ! Eigenvaluesolver for HgH´
      KgK_mat=eig_vec
      call zheev('V','U', hdim, eig_vec, hdim, eig_val, cwork, lwork, rwork, info)
      deallocate(cwork)
      deallocate(rwork)
      !!!      else
      !!!         print *,' Warning in diamag: non-positive definite matrix in zpotrf', iq, info
      !!!      end if
      !print *,'zheeev',info
      call shuffle_eig(eig_val,eig_vec,hdim)
      !eig_val=eig_val/(ry_ev*4.0_dblprec)
      !call dlasrt( 'I', hdim, eig_val, info )
      !
      !ABs
      x_mat=eig_vec
      !do ia=1,NA
      !   do ja=1,NA
      !      call find_uv(ul,vl,emomM(:,ja,1))
      !      bigS(:,ia,ja)=sqrt(mmom(ia,1))*sqrt(0.5_dblprec)*(conjg(ul)*x_mat(ia,ja)+ul*x_mat(ia+NA,ja))+vl*(mmom(ia,1)-x_mat(ia+NA,ja)*x_mat(ia,ja))
      !   end do
      !end do

      ! Calculate L
      call zgemm('C','N',hdim,hdim,hdim,cone,eig_vec,hdim,KgK_mat,hdim,czero,dum_mat,hdim)
      call zgemm('N','N',hdim,hdim,hdim,cone,dum_mat,hdim,eig_vec,hdim,czero,L_mat,hdim)
      
      ! Eigensolve E=g*L
      call zgemm('N','N',hdim,hdim,hdim,cone,g_mat,hdim,L_mat,hdim,czero,E_mat,hdim)
      ! Calculate K^-1
      iK_mat=K_mat
      call ztrtri('U','N',hdim,iK_mat,hdim,info)
      ! E^1/2
      sqE_mat=sqrt(abs(E_mat))
      ! U*E^1/2
      call zgemm('N','N',hdim,hdim,hdim,cone,eig_vec,hdim,sqE_mat,hdim,czero,dum_mat,hdim)
      ! T=K^-1*U*E^1/2
      call zgemm('N','N',hdim,hdim,hdim,cone,iK_mat,hdim,dum_mat,hdim,czero,T_mat,hdim)

      x_mat=0.0_dblprec
      do ia=1,hdim
         call zgemv('N',hdim,hdim,cone,T_mat,hdim,eig_vec(ia,1:hdim),1,cone,x_mat(ia,1:hdim),1)
         call zgemv('N',hdim,hdim,cone,T_mat,hdim,eig_vec(ia,1:hdim),1,cone,x_mat(1:hdim,ia),1)
         call zgemv('N',hdim,hdim,cone,T_mat,hdim,eig_vec(1:hdim,ia),1,cone,x_mat(1:hdim,ia),1)
      end do
      !x_mat=T_mat

      do alfa=1,3
         do beta=1,3
            ! Calculate T´ * [Y Z ; V W ] * T
            call zgemm('N','N',hdim,hdim,hdim,cone,S_prime(:,:,alfa,beta,iq),hdim,T_mat,hdim,czero,dum_mat,hdim)
            call zgemm('C','N',hdim,hdim,hdim,cone,T_mat,hdim,dum_mat,hdim,czero,S_prime(:,:,alfa,beta,iq),hdim)
         end do
      end do
      !
      !print *,'Re x_mat'
      !open(ofileno,file='diaval.simid.out')
      !do ia=1,hdim
      !   write(2000   , '(100f10.4)') eig_val(ia) *13.605698066_dblprec*4.0_dblprec*2.0_dblprec,real(x_mat(:,ia))
      !end do
      !close(ofileno)
      !do ia=1,NA
      !   do ja=1,NA
      !      call find_uv(ul,vl,emomM(:,ja,1))
      !      bigS(:,ia,ja)=sqrt(mmom(ia,1))*sqrt(0.5_dblprec)*(conjg(ul)*x_mat(ia,ja)+ul*conjg(x_mat(ia,ja)))+vl*(mmom(ia,1)-conjg(x_mat(ia,ja))*x_mat(ia,ja))
      !   end do
      !end do
      !
      return
      !
   end subroutine diagonalize_quad_hamiltonian

   subroutine fallback_bosonic_diag(NA, h_in, eig_val, eig_vec, iq)
      use Constants
      implicit none

      integer, intent(in) :: NA, iq
      integer :: hdim, ia, info, lwork
      complex(dblprec), intent(in)  :: h_in(2*NA, 2*NA)
      real(dblprec),    intent(out) :: eig_val(2*NA)
      complex(dblprec), intent(out) :: eig_vec(2*NA, 2*NA)

      complex(dblprec), allocatable :: alpha_mat(:), beta_mat(:), work(:)
      complex(dblprec), allocatable :: eta(:,:), A_copy(:,:)
      real(dblprec),    allocatable :: rwork(:)
      complex(dblprec) :: czero
      complex(dblprec) :: dummy_vl(1,1)

      czero = (0.0_dblprec, 0.0_dblprec)
      hdim = 2 * NA

      !— allocate everything —
      allocate(alpha_mat(hdim), beta_mat(hdim))
      allocate(eta(hdim, hdim), A_copy(hdim, hdim))
      allocate(work(8 * hdim))
      allocate(rwork(hdim))

      eig_vec = czero
      eig_val = 0.0_dblprec

      !— build η = diag(1, …, 1, –1, …, –1) —
      eta = czero
      do ia = 1, NA
         eta(ia,   ia)       = (1.0_dblprec, 0.0_dblprec)
         eta(NA+ia, NA+ia)   = (-1.0_dblprec, 0.0_dblprec)
      end do

      !— copy H into a working array —
      A_copy = h_in

      !- Pretty-print A_copy
      ! print *, 'A_copy'
      ! do ia = 1, hdim
      !    print '(20f12.6)', h_in(ia, 1:hdim)
      ! end do

      !- Add small diagonal offset to ensure positive definiteness
      !— (can be controlled by input parameter `nc_eps`) —
      if (diamag_eps > -1.0_dblprec) then
         do ia = 1, hdim
            A_copy(ia, ia) = A_copy(ia, ia) + diamag_eps
         end do
      else
         do ia = 1, hdim
            A_copy(ia, ia) = A_copy(ia, ia) + 1.0d-6
         end do
      end if

      !— solve A_copy * v = ω * η * v  via ZGGEV —
      call zggev( 'N', 'V', hdim,                       &
                  A_copy, hdim,                         &
                  eta,   hdim,                         &
                  alpha_mat, beta_mat,                 &
                  dummy_vl, 1,                         & ! VL unused
                  eig_vec, hdim,                       & ! VR returned
                  work, 8*hdim,                        & ! complex workspace
                  rwork,                                & ! <<<<< required real workspace
                  info )

      if (info /= 0) then
         print *, 'ERROR: ZGGEV failed at iq=', iq, ' info=', info
         eig_val = 0.0_dblprec
         eig_vec = czero
      else
         do ia = 1, hdim
            if (abs(beta_mat(ia)) > 1d-12) then
               eig_val(ia) = real(alpha_mat(ia) / beta_mat(ia))
            else
               eig_val(ia) = 0.0_dblprec
            end if
         end do
      end if

      !— clean up —
      deallocate(alpha_mat, beta_mat, eta, A_copy, work, rwork)

   end subroutine fallback_bosonic_diag


   subroutine diamag_sqw(NA,nq,nq_ext,nw,q,nc_eval_q,S_prime,simid,n_vec)
      !
      use Constants 
      !
      implicit none
      !
      integer, intent(in) :: NA      !< Number of atoms in cell
      integer, intent(in) :: nq      !< Number of k-vectors
      integer, intent(in) :: nq_ext  !< Extended number of k-vectors
      integer, intent(in) :: nw      !< Number of freqiencies
      real(dblprec), dimension(3,0:nq), intent(in) :: q !< Qpoints 
      real(dblprec), dimension(2*NA,nq_ext), intent(in)  :: nc_eval_q !< Eigenvalues
      complex(dblprec), dimension(2*NA,2*NA,3,3,nq_ext), intent(in)  :: S_prime !< Dynamical matrix
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3), intent(in) :: n_vec  !< Spin spiral pitch vector
      !
      complex(dblprec), dimension(:,:,:,:), allocatable :: sqw_prime !< Dynamical correlation function (primed)
      complex(dblprec), dimension(:,:,:,:), allocatable :: sqw_diamag !< Dynamical correlation function
      complex(dblprec), dimension(:,:), allocatable :: sqw_ortho    !< Orthogonal projection of DCF
      real(dblprec), dimension(:), allocatable :: hbarw !< Energies for corresponding w
      real(dblprec), dimension(:), allocatable :: bdist !< BE distributio density
      !
      real(dblprec) :: qnorm, qnorm2, polfac
      real(dblprec), dimension(3,3) ::  unit3 
      complex(dblprec), dimension(3,3) ::  R1, R2, eye, skew
      complex(dblprec), dimension(:,:), allocatable :: sqwint_diamag
      !
      integer :: iq, ia, iw, hdim
      integer :: i_all, i_stat
      integer :: alfa, beta
      real(dblprec) :: eig_max, delta_e
      real(dblprec) :: dtemp !< Temperature in correlation function
      character(len=30) :: filn
      complex(dblprec) :: cone = (0.0_dblprec,1.0_dblprec)

      hdim=2*NA
      !dtemp=1.0e-3_dblprec
      dtemp=5.0_dblprec

      ! Allocate dynamical correlation function matrix
      allocate(sqw_prime(3,3,0:nw,nq_ext),stat=i_stat)
      call memocc(i_stat,product(shape(sqw_prime))*kind(sqw_prime),'sqw_prime','diamag_sqw')
      sqw_prime=0.0_dblprec
      allocate(sqw_diamag(3,3,0:nw,nq),stat=i_stat)
      call memocc(i_stat,product(shape(sqw_diamag))*kind(sqw_diamag),'sqw_diamag','diamag_sqw')
      sqw_diamag=0.0_dblprec
      allocate(sqw_ortho(0:nw,nq),stat=i_stat)
      call memocc(i_stat,product(shape(sqw_ortho))*kind(sqw_ortho),'sqw_ortho','diamag_sqw')
      sqw_diamag=0.0_dblprec

      ! Allocate energy array
      allocate(hbarw(0:nw),stat=i_stat)
      call memocc(i_stat,product(shape(hbarw))*kind(hbarw),'hbarw','diamag_sqw')
      allocate(bdist(0:nw),stat=i_stat)
      call memocc(i_stat,product(shape(bdist))*kind(bdist),'bdist','diamag_sqw')
      allocate(sqwint_diamag(0:nw,nq),stat=i_stat)
      call memocc(i_stat,product(shape(sqwint_diamag))*kind(sqwint_diamag),'sqwint_diamag','diamag_sqw')

      ! Fill energy array
      eig_max=maxval(nc_eval_q)
      delta_e=1.1_dblprec*eig_max/(1.0_dblprec*nw)
      do iw=0,nw
         hbarw(iw)=iw*delta_e
         bdist(iw)=1.0_dblprec/(exp(1.0e-3_dblprec*hbarw(iw)/(k_bolt_ev*dtemp))-1.0_dblprec)
         !write (999,'(i5,2g20.10)') iw, hbarw(iw), bdist(iw)
      end do

      ! Construct rotation matrices R1, R2
      R2(1,1)=n_vec(1)*n_vec(1); R2(1,2)=n_vec(1)*n_vec(2); R2(1,3)=n_vec(1)*n_vec(3)
      R2(2,1)=n_vec(2)*n_vec(1); R2(2,2)=n_vec(2)*n_vec(2); R2(2,3)=n_vec(2)*n_vec(3)
      R2(3,1)=n_vec(3)*n_vec(1); R2(3,2)=n_vec(3)*n_vec(2); R2(3,3)=n_vec(3)*n_vec(3)
      eye=0.0_dblprec;eye(1,1)=1.0_dblprec;eye(2,2)=1.0_dblprec;eye(3,3)=1.0_dblprec
      skew=0.0_dblprec;
      skew(1,2)=-n_vec(3);skew(2,1)= n_vec(3)
      skew(1,3)= n_vec(2);skew(3,1)=-n_vec(2)
      skew(2,3)=-n_vec(1);skew(3,2)= n_vec(1)
      skew=-skew
      R1=eye-cone*skew-R2
      R1=0.5_dblprec*R1

      sqw_prime=0.0_dblprec

      ! Loop over q-vectors
      do iq=1,nq_ext/2
         !write(999,'(9f18.8)') abs(S_prime(1,1,:,:,iq))
         do alfa=1,3
            do beta=1,3
               do iw=0,nw
                  do ia=1,na
                     sqw_prime(alfa,beta,iw,iq) = sqw_prime(alfa,beta,iw,iq) &
                        + S_prime(ia,ia,alfa,beta,iq) &    ! S_prime
                        * real(exp(-0.1_dblprec*(hbarw(iw)-(nc_eval_q(ia,iq)))**2)) &   ! delta function (check broadening, now 10)
                        * (0*bdist(iw)+0.5_dblprec*(1.0_dblprec+1.0_dblprec))  ! occupation (Check gii)
                  end do
                  !!Negative eigenvalues 
                  do ia=na+1,2*na
                     sqw_prime(alfa,beta,iw,iq)=sqw_prime(alfa,beta,iw,iq) &
                        + S_prime(ia,ia,alfa,beta,iq) &    ! S_prime
                        * real(exp(-0.1_dblprec*(hbarw(iw)+(nc_eval_q(ia,iq)))**2)) &   ! delta function (check broadening, now 10)
                        * (0*bdist(iw)+0.5_dblprec*(1.0_dblprec-1.0_dblprec))  ! occupation (Check gii)
                  end do
               end do
            end do
         end do
      end do
      !
      do iq=1,nq
         do iw=0,nw
            sqw_diamag(:,:,iw,iq)=matmul(sqw_prime(:,:,iw,iq),R2) + &
               matmul(sqw_prime(:,:,iw,iq+nq),(R1)) + &
               matmul(sqw_prime(:,:,iw,iq+2*nq),conjg(R1))
            !            sqw_diamag(:,:,iw,iq)=sqw_prime(:,:,iw,iq)
         end do
      end do
      !
      write (filn,'(''ncsqw.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position='append')
      do iq=1,nq
         do iw=1,nw
            !write(ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, abs(sqw_diamag(:,:,iw,iq))
            write(ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, real(sqw_diamag(:,:,iw,iq))
         end do
      end do
      close(ofileno)
      !
      !!! write (filn,'(''ncsqw_ortho.'',a,''.out'')') trim(simid)
      !!! open(ofileno,file=filn, position='append')
      !!! do iq=1,nq
      !!!    do iw=1,nw
      !!!       sqw_ortho(iw,iq)=0.0_dblprec
      !!!       do alfa=1,3
      !!!          do beta=1,3
      !!!             sqw_ortho(iw,iq) = sqw_ortho(iw,iq) + (1.0_dblprec-diamag_qvect(alfa)*diamag_qvect(beta)*sqw_diamag(alfa,beta,iw,iq)
      !!!          end do
      !!!       end do
      !!!       !write(ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, abs(sqw_diamag(:,:,iw,iq))
      !!!       write(ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, real(sqw_ortho(iw,iq)) , aimag(sqw_ortho(iw,iq))
      !!!    end do
      !!! end do
      !!! close(ofileno)
      !      Calculate the scattering intensity from the tensorial S(q,w)
      ! Write its real, imaginary, and absolute value to file.
      write (filn,'(''ncsqw_intensity.'',a,''.out'')') trim(simid)
      unit3 = 0.0_dblprec
      unit3(1,1) = 1.0_dblprec;  unit3(2,2) = 1.0_dblprec; unit3(3,3) = 1.0_dblprec
      open(ofileno, file=filn)
      do iq=1,Nq
         qnorm = norm2(q(1:3,iq))
         qnorm2 = qnorm**2+1.0e-12_dblprec
         do iw=0,nw
            do alfa=1,3
               do beta=1,3
                  polfac = unit3(alfa,beta) - q(alfa,iq) * q(beta,iq) / qnorm2
                  sqwint_diamag(iw,iq) = sqwint_diamag(iw,iq) + polfac * sqw_diamag(alfa,beta,iw,iq)
               end do
            end do
            write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, real(sqwint_diamag(iw,iq)), &
               aimag(sqwint_diamag(iw,iq)), abs(sqwint_diamag(iw,iq))
         end do
      end do
      close(ofileno)


      i_all=-product(shape(bdist))*kind(bdist)
      deallocate(bdist,stat=i_stat)
      call memocc(i_stat,i_all,'bdist','diamag_sqw')
      i_all=-product(shape(hbarw))*kind(hbarw)
      deallocate(hbarw,stat=i_stat)
      call memocc(i_stat,i_all,'hbarw','diamag_sqw')

      i_all=-product(shape(sqwint_diamag))*kind(sqwint_diamag)
      deallocate(sqwint_diamag,stat=i_stat)
      call memocc(i_stat,i_all,'sqwint_diamag','diamag_sqw')
      i_all=-product(shape(sqw_ortho))*kind(sqw_ortho)
      deallocate(sqw_ortho,stat=i_stat)
      call memocc(i_stat,i_all,'sqw_ortho','diamag_sqw')
      i_all=-product(shape(sqw_diamag))*kind(sqw_diamag)
      deallocate(sqw_diamag,stat=i_stat)
      call memocc(i_stat,i_all,'sqw_diamag','diamag_sqw')
      i_all=-product(shape(sqw_prime))*kind(sqw_prime)
      deallocate(sqw_prime,stat=i_stat)
      call memocc(i_stat,i_all,'sqw_prime','diamag_sqw')

      return
      !
      10005 format (i7,3f10.6,2x,i7,2x,9es16.8)
      !
   end subroutine diamag_sqw


   !!!    subroutine moments_to_angles(moments,ndim,eig_val)
   !!!       !
   !!!       !
   !!!       implicit none
   !!!       !
   !!!       !
   !!!       integer, intent(in) :: ndim
   !!!       real(dblprec), dimension(ndim,ndim), intent(in) :: moments
   !!!       real(dblprec), dimension(ndim), intent(in) :: eig_val
   !!!       !
   !!!       real(dblprec) :: pi, mi_norm, mj_norm, mdot
   !!!       real(dblprec), dimension(ndim/3,ndim) :: angles
   !!!       real(dblprec), dimension(3) :: mom_i, mom_j
   !!! 
   !!!       integer :: i,j, idx
   !!! 
   !!!       angles=0.0_dblprec
   !!!       print *,'Angles:'
   !!!       pi=acos(-1.0_dblprec)
   !!!       do j=1,ndim
   !!!          idx=0
   !!!          mom_i(1:3)=moments(1:3,j)
   !!!          mi_norm=sqrt(sum(mom_i*mom_i))+1.0d-15
   !!!          do i=1,ndim,3
   !!!             mom_j(1:3)=moments(i:i+2,j)
   !!!             mj_norm=sqrt(sum(mom_j*mom_j))+1.0d-15
   !!!             mdot=sum(mom_i*mom_j)/(mi_norm*mj_norm)
   !!!             idx=idx+1
   !!!             angles(idx,j)=acos(mdot)*180.0_dblprec/pi
   !!!          end do
   !!!       end do
   !!!       !do j=1,ndim,3
   !!!       !   print '(1x,20f10.4)', eig_val(j), angles(:,j)
   !!!       !end do
   !!!       !
   !!!       return
   !!!       !
   !!!    end subroutine moments_to_angles

   !!!    subroutine normalize_moments(moments,nrow,ncol)
   !!!       !
   !!!       implicit none
   !!!       !
   !!!       !
   !!!       integer, intent(in) :: nrow
   !!!       integer, intent(in) :: ncol
   !!!       real(dblprec), dimension(ncol,nrow), intent(inout) :: moments
   !!!       !
   !!!       real(dblprec) :: mnorm
   !!!       integer :: i,j
   !!!       !
   !!!       do i=1,nrow,3
   !!!          do j=1,ncol
   !!!             mnorm=sum(moments(j,i:i+2)**2)**0.5_dblprec+1.0d-15
   !!!             moments(j,i+0)=moments(j,i+0)/mnorm
   !!!             moments(j,i+1)=moments(j,i+1)/mnorm
   !!!             moments(j,i+2)=moments(j,i+2)/mnorm
   !!!          end do
   !!!       end do
   !!!       !
   !!!       return
   !!!       !
   !!!    end subroutine normalize_moments

   subroutine find_uv(u,v,mom)
      !
      use math_functions, only : f_cross_product
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
      real(dblprec), dimension(3), parameter :: x_hat = (/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
      real(dblprec), dimension(3), parameter :: y_hat = (/0.0_dblprec,1.0_dblprec,0.0_dblprec/)
      real(dblprec), dimension(3), parameter :: z_hat = (/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
      !
      mnorm=sqrt(sum(mom*mom)+1.0d-15)
      ! Find the rotation matrix R_prime that transforms mom to z-axis.
      ! Done by finding the vector orthogonal to mom and z and then perform a
      ! Rodrigues rotation with the corresponding angle
      !
      ! Unit vector parallel to mom
      mom_hat=mom/mnorm
      !mnorm=1.0_dblprec
      !mom_hat=z_hat
      ! Follow Toth-Lake recipe for R´
      R_prime(3,:)=mom_hat
      !if(abs(mom_hat(3))==1.0_dblprec) then  !if m=00+-1 then R_2=0+-10
      if(abs(mom_hat(3)-1.0_dblprec)<0.9999_dblprec) then  !if m=00+-1 then R_2 = x^  cross m (not R_2=0+-10)

         R_prime(2,:)=f_cross_product(x_hat,mom_hat)

         !!! !R_prime(3,1)=0.0_dblprec;R_prime(3,2)=0.0_dblprec;R_prime(3,3)=1.0_dblprec
         !!! !R_prime(2,1)=0.0_dblprec;R_prime(2,2)=1.0_dblprec;R_prime(2,3)=0.0_dblprec
      else  !else R_2=001 x m
         R_prime(2,:)=f_cross_product(z_hat,mom_hat)
         !!! !R_prime(2,1)=-mom_hat(2)
         !!! !R_prime(2,2)=mom_hat(3)
         !!! !!R_prime(2,1)=mom_hat(2)
         !!! !!R_prime(2,2)=-mom_hat(3)
         !!! !R_prime(2,3)=0.0_dblprec
         !!! R_prime(2,1)=mom_hat(2)*1.0_dblprec-0.0_dblprec*mom_hat(3)
         !!! R_prime(2,2)=mom_hat(3)*0.0_dblprec-1.0_dblprec*mom_hat(1)
         !!! R_prime(2,3)=mom_hat(1)*0.0_dblprec-0.0_dblprec*mom_hat(2)
         mnorm=sqrt(sum(R_prime(2,:)*R_prime(2,:)))
         R_prime(2,:)=R_prime(2,:)/mnorm
      end if  !R_1 = R_2 x R_3
      R_prime(1,:)=f_cross_product(R_prime(2,:),R_prime(3,:))
      !!! !R_prime(1,1)=R_prime(2,2)*R_prime(3,3)-R_prime(2,3)*R_prime(3,2)
      !!! !R_prime(1,2)=R_prime(2,3)*R_prime(3,1)-R_prime(2,1)*R_prime(3,3)
      !!! !R_prime(1,3)=R_prime(2,1)*R_prime(3,2)-R_prime(2,2)*R_prime(3,1)
      mnorm=sqrt(sum(R_prime(1,:)*R_prime(1,:)))
      R_prime(1,:)=R_prime(1,:)/mnorm

      u=R_prime(1,:)+im*R_prime(2,:)
      v=R_prime(3,:)
      !
      return
      !
   end subroutine find_uv

   subroutine find_R(R,mom_0,mom_n,e_rot)
      !
      implicit none
      !
      real(dblprec), dimension(3,3), intent(out) :: R
      real(dblprec), dimension(3), intent(in) :: mom_0
      real(dblprec), dimension(3), intent(in) :: mom_n
      real(dblprec), dimension(3), intent(in), optional :: e_rot
      !

      real(dblprec) :: mnorm
      real(dblprec) :: e_norm
      real(dblprec) :: angle
      real(dblprec) :: mdot
      real(dblprec), dimension(3,3) :: K, K2, one
      real(dblprec), dimension(3) :: mom_hat_0, mom_hat_n
      real(dblprec), dimension(3) :: k_hat
      real(dblprec), dimension(3) :: x,y,z
      !
      !
      ! Find the rotation matrix R that transforms mom_n to mom_0
      ! Done by finding the vector orthogonal to mom_n and mom_0 and then perform a
      ! Rodrigues rotation with the corresponding angle
      !
      x(1)=1.0_dblprec;x(2)=0.0_dblprec;x(3)=0.0_dblprec
      y(1)=0.0_dblprec;y(2)=1.0_dblprec;y(3)=0.0_dblprec
      z(1)=0.0_dblprec;z(2)=0.0_dblprec;z(3)=1.0_dblprec
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
      if (mdot>1.0_dblprec) then
         mdot=1.0_dblprec
      else if (mdot<-1.0_dblprec) then
         mdot=-1.0_dblprec
      end if

      angle=acos(mdot)

      ! Find perpendicular vector from cross-product
      ! print '(a,3f18.8)', 'k_hat',k_hat
      ! print '(a,f18.8)', 'angle: ',angle*180/3.1415
      !
      if (present(e_rot)) then
         e_norm=norm2(e_rot)
         if (e_norm>0.0_dblprec) then
            k_hat=e_rot/e_norm
         else
            k_hat=z
         end if

      else


         ! If m_0 and m_n are parallel, then use other approach
         if(1.0_dblprec-abs(mdot)<1.0e-9_dblprec) then
            ! Two additional checks: first if m_0 // x
            if(1.0_dblprec-abs(sum(mom_hat_0*x))<1.0e-5_dblprec) then   ! x
               k_hat(1)=mom_hat_0(2)*y(3)-mom_hat_0(3)*y(2)
               k_hat(2)=mom_hat_0(3)*y(1)-mom_hat_0(1)*y(3)
               k_hat(3)=mom_hat_0(1)*y(2)-mom_hat_0(2)*y(1)
            else if(1.0_dblprec-abs(sum(mom_hat_0*y))<1.0e-5_dblprec) then  ! y
               k_hat(1)=mom_hat_0(2)*z(3)-mom_hat_0(3)*z(2)
               k_hat(2)=mom_hat_0(3)*z(1)-mom_hat_0(1)*z(3)
               k_hat(3)=mom_hat_0(1)*z(2)-mom_hat_0(2)*z(1)
            else  ! z
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
      end if
      k_hat=k_hat/sqrt(sum(k_hat*k_hat)+1.0e-12)
      K(2,1)= k_hat(3);K(1,2)=-k_hat(3)
      K(3,1)=-k_hat(2);K(1,3)= k_hat(2)
      K(2,3)=-k_hat(1);K(3,2)= k_hat(1)
      K2=matmul(K,K)
      R=one+K*sin(angle)+(1.0_dblprec-cos(angle))*K2
      !print *,'---I---'
      !print '(3f18.8)',one
      !print *,'---K---'
      !print '(3f18.8)',K
      !print *,'---K2--'
      !print '(3f18.8)',K2
      !print *,'---R---'
      !print '(3f18.8)',R
      !
      !print '(3f18.8)',k_hat
      return
      !
   end subroutine find_R


   !!!    subroutine find_O(O,mom_0)
   !!!       !
   !!!       implicit none
   !!!       !
   !!!       complex(dblprec), dimension(3,3), intent(out) :: O
   !!!       real(dblprec), dimension(3), intent(in) :: mom_0
   !!!       !
   !!!       real(dblprec) :: theta, phi
   !!!       real(dblprec), dimension(3,3) :: Oz,Oy
   !!!       !
   !!!       !
   !!!       theta=acos(mom_0(1))
   !!!       phi=acos(mom_0(3))
   !!!       Oz=0.0_dblprec
   !!!       Oy=0.0_dblprec
   !!!       Oz(1,1)=cos(phi);Oz(2,2)=cos(phi);Oz(3,3)=1.0_dblprec
   !!!       Oz(2,1)=-sin(phi);Oz(1,2)=sin(phi)
   !!!       Oy(1,1)=cos(theta);Oy(2,2)=1.0_dblprec;Oy(3,3)=cos(theta)
   !!!       Oy(3,1)=-sin(theta);Oy(1,3)=sin(theta)
   !!!       O=matmul(Oz,Oy)
   !!!       !
   !!!       !
   !!!       return
   !!!       !
   !!!    end subroutine find_O

   !!!    subroutine get_M(M,mom_0)
   !!!       !
   !!!       implicit none
   !!!       !
   !!!       complex(dblprec), dimension(3,3), intent(out) :: M
   !!!       real(dblprec), intent(in) :: mom_0
   !!!       !
   !!!       complex(dblprec) :: im
   !!!       !
   !!!       im=(0.0_dblprec,1.0_dblprec)
   !!!       !
   !!!       M=0.0_dblprec
   !!!       M(1,1)=1.0_dblprec;M(1,2)=1.0_dblprec
   !!!       M(2,1)=-im ;M(1,2)=im
   !!!       M(3,3)=sqrt(2.0_dblprec/mom_0)
   !!!       M=sqrt(mom_0/2.0_dblprec)
   !!!       !
   !!!       !
   !!!       return
   !!!       !
   !!!    end subroutine get_M

   !!!   subroutine get_Hij(Hij,Jij,emom_i,emom_j,mmom_i,mmom_j,i,j)
   !!!       !
   !!!       implicit none
   !!!       !
   !!!       complex(dblprec), dimension(2,2), intent(out) :: Hij
   !!!       real(dblprec), dimension(3,3), intent(in) :: Jij
   !!!       real(dblprec), dimension(3), intent(in) :: emom_i
   !!!       real(dblprec), dimension(3), intent(in) :: emom_j
   !!!       real(dblprec), intent(in) :: mmom_i
   !!!       real(dblprec), intent(in) :: mmom_j
   !!!       integer, intent(in) :: i
   !!!       integer, intent(in) :: j
   !!!       !
   !!!       complex(dblprec), dimension(3,3) :: Jbar, M_i, M_j
   !!!       complex(dblprec), dimension(3,3) :: O_i, O_j
   !!!       complex(dblprec), dimension(3,3) :: Mprod
   !!!       !
   !!!       call find_O(O_i,emom_i)
   !!!       call find_O(O_j,emom_j)
   !!!       call get_M(M_i,mmom_i)
   !!!       call get_M(M_j,mmom_i)
   !!!       !
   !!!       Mprod=M_j
   !!!       Mprod=matmul(O_j,Mprod)
   !!!       Mprod=matmul(Jij,Mprod)
   !!!       Mprod=matmul(transpose(O_i),Mprod)
   !!!       Mprod=matmul(transpose(conjg(M_i)),Mprod)
   !!!       Hij=Mprod(1:2,1:2)
   !!!       if(i==j) then
   !!!          Hij(1,1)=Hij(1,1)-Mprod(3,3)
   !!!          Hij(2,2)=Hij(2,2)-Mprod(3,3)
   !!!       end if
   !!!       !
   !!!       !
   !!!       return
   !!!       !
   !!!      end subroutine get_Hij
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
   subroutine setup_ektij(Natom,NA,q,nq,ektij)

      use Constants
      use SystemData, only : coord
      use Math_functions, only : f_wrap_coord_diff
      !
      implicit none
      !
      integer, intent(in) :: Natom   !< Number of atoms
      integer, intent(in) :: NA      !< Number of types
      integer, intent(in) :: nq      !< Number of q-points
      real(dblprec), dimension(3,0:nq), intent(in) :: q !< Qpoints 
      complex(dblprec), dimension(na,na,0:nq), intent(out) :: ektij

      !
      real(dblprec), dimension(3) :: q_vec
      real(dblprec), dimension(3) :: dist
      complex(dblprec)  :: FTfac, im
      integer :: iq, ia, ja
      !
      !real(dblprec), dimension(3,3,na,na,0:nq) :: Jtens_q
      !
      im=(0.0_dblprec,1.0_dblprec)
      !
      do iq=0,nq
         ! Ensure that iq=0 corresponds to gamma
         if(iq>0) then
            q_vec=q(:,iq)*2.0_dblprec*pi
         else
            q_vec=0.0_dblprec
         end if

         do ia=1,NA
            ! Jij exchange
            do ja=1,NA
               call f_wrap_coord_diff(Natom,coord,ia,ja,dist)

               FTfac=exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))
               ektij(ia,ja,iq)=FTfac

            end do
         end do

      end do

      return
      !
   end subroutine setup_ektij

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
      dtens(1,3)= dmvec(2)
      dtens(2,3)=-dmvec(1)
      dtens(3,3)= 0.0_dblprec

      return
      !
   end function dm2tens

   function sa2tens(savec) result(satens)
      !
      implicit none
      !
      real(dblprec), dimension(3), intent(in) :: savec !< an SA vector
      !
      real(dblprec), dimension(3,3) :: satens
      !

      satens(1,1)= 0.0_dblprec
      satens(2,1)= savec(3)
      satens(3,1)= savec(2)
      satens(1,2)= savec(3)
      satens(2,2)= 0.0_dblprec
      satens(3,2)= savec(1)
      satens(1,3)= savec(2)
      satens(2,3)= savec(1)
      satens(3,3)= 0.0_dblprec

      return
      !
   end function sa2tens

   function pd2tens(pdvec) result(pdtens)
      !
      implicit none
      !
      real(dblprec), dimension(9), intent(in) :: pdvec !< an SA vector
      !
      real(dblprec), dimension(3,3) :: pdtens
      !

      pdtens(1,1)= pdvec(1)
      pdtens(2,1)= pdvec(2)
      pdtens(3,1)= pdvec(3)
      pdtens(1,2)= pdvec(4)
      pdtens(2,2)= pdvec(5)
      pdtens(3,2)= pdvec(6)
      pdtens(1,3)= pdvec(7)
      pdtens(2,3)= pdvec(8)
      pdtens(3,3)= pdvec(9)

      return
      !
   end function pd2tens

   function ani2tens(anivec) result(anitens)
      !
      implicit none
      !
      real(dblprec), dimension(3), intent(in) :: anivec !< an SA vector
      !
      real(dblprec), dimension(3,3) :: anitens
      !

      anitens(1,:)= anivec(1)*anivec(:)
      anitens(2,:)= anivec(2)*anivec(:)
      anitens(3,:)= anivec(3)*anivec(:)

      return
      !
   end function ani2tens

   subroutine clone_q(nq,q,nq_ext,q_ext,diamag_qvect)
      !
      implicit none
      !
      integer, intent(in) :: nq !< No. qpoints from input
      real(dblprec), dimension(3,nq), intent(in) :: q !< qpoints from input
      integer, intent(in) :: nq_ext !< Extended no. qpoints 
      real(dblprec), dimension(3,0:nq_ext), intent(out) :: q_ext !< Extended set of q-points
      real(dblprec), dimension(3), intent(in) :: diamag_qvect !< Spin spiral q-vector
      !
      integer :: iq
      !

      q_ext(:,0)=0.0_dblprec
      if (norm2(diamag_qvect)>0.0_dblprec) print '(2x,a,3f12.6)','Ordering vector: ',diamag_qvect
      do iq=1,nq
         ! Q = q, q+q_0, q-q0
         q_ext(:,iq+0*nq)=q(:,iq)
         q_ext(:,iq+1*nq)=q(:,iq)+diamag_qvect
         q_ext(:,iq+2*nq)=q(:,iq)-diamag_qvect
         ! Positive and negative q
         q_ext(:,iq+3*nq)=-q_ext(:,iq+0*nq)
         q_ext(:,iq+4*nq)=-q_ext(:,iq+1*nq)
         q_ext(:,iq+5*nq)=-q_ext(:,iq+2*nq)
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
      character(len=50) :: keyword
      integer :: rd_len,i_err,i_errb
      logical :: comment

      call setup_diamag()
      diamag_nvect(1:2)=0.0_dblprec;diamag_nvect(3)=1.0_dblprec

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

         case('nc_qvect') ! Perform frequency based spin-correlation sampling
            read(ifile,*,iostat=i_err) diamag_qvect
            if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

         case('nc_nvect') ! Perform frequency based spin-correlation sampling
            read(ifile,*,iostat=i_err) diamag_nvect
            if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

         case('nc_eps') ! Perform frequency based spin-correlation sampling
            read(ifile,*,iostat=i_err) diamag_eps
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
subroutine setup_Jtens_q(Natom,Mensemble,NA,emomM,q,nq,Jtens_q)

   use Constants
   use SystemData, only : coord
   use Math_functions, only : f_wrap_coord_diff
   !
   implicit none
   !
   integer, intent(in) :: Natom   !< Number of atoms
   integer, intent(in) :: Mensemble !< Number of ensembles
   integer, intent(in) :: NA      !< Number of types
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
   integer, intent(in) :: nq      !< Number of q-points
   real(dblprec), dimension(3,0:nq), intent(in) :: q !< Qpoints 
   complex(dblprec), dimension(3,3,na,na,0:nq), intent(out) :: Jtens_q

   !
   real(dblprec), dimension(3) :: q_vec
   real(dblprec), dimension(3) :: dist
   real(dblprec), dimension(3) :: dmv, cmv, amv
   real(dblprec), dimension(9) :: pmv
   real(dblprec), dimension(3) :: z
   real(dblprec), dimension(3,3) :: R_n, J_n, D_n, C_n, A_n, P_n
   complex(dblprec)  :: FTfac, im
   integer :: iq, ia, ja, j, jat, ih
   !
   !real(dblprec), dimension(3,3,na,na,0:nq) :: Jtens_q
   !
   im=(0.0_dblprec,1.0_dblprec)
   Jtens_q=0.0_dblprec
   z(1)=0.0_dblprec;z(2)=0.0_dblprec;z(3)=1.0_dblprec
   !
   !!!$omp parallel do default(shared) private(iq,q_vec,ia,ja,j,jat,FTfac,dmv,cmv,pmv,amv, J_n, D_n, C_n, P_n, A_n)
   do iq=0,nq
      ! Ensure that iq=0 corresponds to gamma
      if(iq>0) then
         q_vec=q(:,iq)*2.0_dblprec*pi
      else
         q_vec=0.0_dblprec
      end if
      !!! print '(2x,a,i5,3f18.8)' , 'q=',iq, q_vec

      do ia=1,NA
         ih=ham%aHam(ia)
         ! Jij, Dij, Cij exchange (pair interactions)
         do j=1,ham%nlistsize(ih)
            ja=ham%nlist(j,ia)
            !!! jat=atype(ja)
            jat=mod(ja-1,NA)+1

            call f_wrap_coord_diff(Natom,coord,ia,ja,dist)
            call find_R(R_n,emomM(:,jat,1),emomM(:,ja,1))
            !write(10000,'(a,2i6,3f10.4)')'------',ia,ja,q(:,iq)
            !write(10000,'(3f12.6)')real(R_n)


            J_n=0.0_dblprec
            J_n(1,1)=ham%ncoup(j,ih,1)
            J_n(2,2)=ham%ncoup(j,ih,1)
            J_n(3,3)=ham%ncoup(j,ih,1)
            if (ham_inp%do_dm==1) then
               dmv=ham%dm_vect(:,j,ih)
               D_n=dm2tens(dmv)
               J_n=J_n-D_n
            end if
            if (ham_inp%do_sa==1) then
               cmv=-ham%sa_vect(:,j,ih)
               C_n=sa2tens(cmv)
               J_n=J_n+C_n
            end if
            if (ham_inp%do_pd==1) then
               pmv=-ham%pd_vect(:,j,ih)
               P_n=pd2tens(pmv)
               J_n=J_n+P_n
            end if

            FTfac=exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))

            !print '(a,2i6,3f10.4)','----------------',ia,jat,q(:,iq)
            !print '(3f12.6)',real(R_n)
            !print '(a,2i6,3f10.4)','------'
            !print '(3f12.6)',real(J_n)
            !print '(a,2i6,3f10.4)','------'
            !print '(3f12.6)',real(matmul(J_n,R_n))

            Jtens_q(:,:,ia,jat,iq)=Jtens_q(:,:,ia,jat,iq)+matmul(J_n,R_n)*FTfac
         end do
         ! Anisotropies (on-site interactions)
         ! Anisotropy
         if (ham_inp%do_anisotropy==1) then
            A_n = 0.0_dblprec
            if (ham%taniso(ia)==1) then
               ! Uniaxial anisotropy
               !call uniaxial_anisotropy_field(i, k, beff_s,Natom,Mensemble,ham_inp%mult_axis,emomM)
               amv = ham%eaniso(:,ia)
               A_n = - 2.0_dblprec * ham%kaniso(1,ia)*ani2tens(amv)
               if (ham_inp%mult_axis=='Y') then
                  amv = ham%eaniso_diff(:,ia)
                  A_n = A_n - 2.0_dblprec * ham%kaniso_diff(1,ia) * ani2tens(amv)
               end if
            elseif (ham%taniso(ia)==2) then
               ! Cubic anisotropy
               !call cubic_anisotropy_field(i, k, beff_s,Natom,Mensemble,ham_inp%mult_axis,emomM)
               print *,'cub not supported yet'
            elseif (ham%taniso(ia)==7)then
               ! Uniaxial and cubic anisotropy
               !call uniaxial_anisotropy_field(i, k, beff_s,Natom,Mensemble,ham_inp%mult_axis,emomM)
               !tfield=0.0_dblprec
               !call cubic_anisotropy_field(i, k, tfield,Natom,Mensemble,ham_inp%mult_axis,emomM)
               !beff_q=beff_q+tfield*ham%sb(i)
               print *,'both not supported yet'
            endif
            Jtens_q(:,:,ia,ia,iq)=Jtens_q(:,:,ia,ia,iq)+A_n !*FTfac
         endif
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


!> Set up the Hamiltonian for first cell
subroutine setup_Jtens2_q(Natom,Mensemble,NA,emomM,q,nq,Jtens_q)

   use Constants
   use SystemData, only : coord
   use Math_functions, only : f_wrap_coord_diff
   !
   implicit none
   !
   integer, intent(in) :: Natom   !< Number of atoms
   integer, intent(in) :: Mensemble !< Number of ensembles
   integer, intent(in) :: NA      !< Number of types
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
   integer, intent(in) :: nq      !< Number of q-points
   real(dblprec), dimension(3,0:nq), intent(in) :: q !< Qpoints 
   complex(dblprec), dimension(3,3,na,na,0:nq), intent(out) :: Jtens_q

   !
   real(dblprec), dimension(3) :: q_vec
   real(dblprec), dimension(3) :: dist
   real(dblprec), dimension(3) :: amv
   real(dblprec), dimension(3) :: z
   real(dblprec), dimension(3,3) :: R_n, J_n, A_n
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

      do ia=1,NA
         ih=ham%aHam(ia)
         ! Jij exchange
         do j=1,ham%nlistsize(ih)
            ja=ham%nlist(j,ia)

            jat=mod(ja-1,NA)+1
            call f_wrap_coord_diff(Natom,coord,ia,ja,dist)

            call find_R(R_n,emomM(:,jat,1),emomM(:,ja,1))

            J_n=ham%j_tens(:,:,j,ih)

            FTfac=exp(-im *(q_vec(1)*dist(1)+q_vec(2)*dist(2)+q_vec(3)*dist(3)))

            Jtens_q(:,:,ia,jat,iq)=Jtens_q(:,:,ia,jat,iq)+matmul(J_n,R_n)*FTfac

         end do
         ! Anisotropies (on-site interactions)
         ! Anisotropy
         if (ham_inp%do_anisotropy==1) then
            A_n = 0.0_dblprec
            if (ham%taniso(ia)==1) then
               ! Uniaxial anisotropy
               !call uniaxial_anisotropy_field(i, k, beff_s,Natom,Mensemble,ham_inp%mult_axis,emomM)
               amv = ham%eaniso(:,ia)
               !A_n = + 4.0_dblprec * ham%kaniso(1,ia)*ani2tens(amv)
               A_n = - 2.0_dblprec * ham%kaniso(1,ia)*ani2tens(amv)
               if (ham_inp%mult_axis=='Y') then
                  amv = ham%eaniso_diff(:,ia)
                  A_n = A_n - 2.0_dblprec * ham%kaniso_diff(1,ia) * ani2tens(amv)
               end if
            elseif (ham%taniso(ia)==2) then
               ! Cubic anisotropy
               !call cubic_anisotropy_field(i, k, beff_s,Natom,Mensemble,ham_inp%mult_axis,emomM)
               print *,'cub not supported yet'
            elseif (ham%taniso(ia)==7)then
               ! Uniaxial and cubic anisotropy
               !call uniaxial_anisotropy_field(i, k, beff_s,Natom,Mensemble,ham_inp%mult_axis,emomM)
               !tfield=0.0_dblprec
               !call cubic_anisotropy_field(i, k, tfield,Natom,Mensemble,ham_inp%mult_axis,emomM)
               !beff_q=beff_q+tfield*ham%sb(i)
               print *,'both not supported yet'
            endif
            Jtens_q(:,:,ia,ia,iq)=Jtens_q(:,:,ia,ia,iq)+A_n !*FTfac
         endif
      end do


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
end subroutine setup_Jtens2_q

end module diamag

