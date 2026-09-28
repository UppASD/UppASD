program test_sa_neighbor_lists
   ! Regression for an SA bond whose neighbour list differs from exchange.
   ! The FM moments are collinear and the SA tensor is symmetric and traceless.
   use Parameters, only : dblprec
   use HamiltonianData, only : ham
   use InputData, only : ham_inp, BC1, BC2, BC3, C1, C2, C3
   use SystemData, only : coord
   use diamag, only : setup_Jtens_q
   implicit none

   integer, parameter :: natom=3, na=1, nens=1, nq=1
   real(dblprec) :: emomM(3,natom,nens), q(3,0:nq)
   complex(dblprec) :: Jtens_q(3,3,na,na,0:nq)
   real(dblprec), parameter :: tol=1.0e-10_dblprec

   allocate(coord(3,natom))
   coord(:,1)=(/0.0_dblprec,0.0_dblprec,0.0_dblprec/)
   coord(:,2)=(/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
   coord(:,3)=(/2.0_dblprec,0.0_dblprec,0.0_dblprec/)
   C1=(/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
   C2=(/0.0_dblprec,1.0_dblprec,0.0_dblprec/)
   C3=(/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
   BC1='0'
   BC2='0'
   BC3='0'

   allocate(ham%aHam(na), ham%nlistsize(na), ham%nlist(1,natom), ham%ncoup(1,na,nens))
   allocate(ham%salistsize(na), ham%salist(1,natom), ham%sa_vect(3,1,na))
   ham%aHam=1
   ham%nlistsize=1
   ham%nlist=0
   ham%nlist(1,1)=2
   ham%ncoup=0.0_dblprec
   ham%salistsize=1
   ham%salist=0
   ham%salist(1,1)=3
   ! sa2tens maps this vector to a symmetric, traceless yz tensor.
   ham%sa_vect=0.0_dblprec
   ham%sa_vect(1,1,1)=1.0_dblprec

   ham_inp%do_dm=0
   ham_inp%do_sa=1
   ham_inp%do_pd=0
   ham_inp%do_anisotropy=0
   emomM=0.0_dblprec
   emomM(3,:,:)=1.0_dblprec
   q=0.0_dblprec
   q(1,1)=0.25_dblprec

   call setup_Jtens_q(natom,nens,na,emomM,q,nq,Jtens_q)

   ! The SA bond is the second neighbour, so q=1/4 gives exp(-i*pi)=-1.
   ! With cmv=-sa_vect, the (2,3) tensor element is -1 before the phase.
   if (abs(real(Jtens_q(2,3,1,1,1))-1.0_dblprec)>tol .or. &
       abs(aimag(Jtens_q(2,3,1,1,1)))>tol) then
      error stop 'SA neighbour-list regression failed'
   end if
end program test_sa_neighbor_lists
