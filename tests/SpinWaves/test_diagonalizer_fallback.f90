program test_diagonalizer_fallback
   use Parameters, only : dblprec
   use diamag, only : setup_diamag, diagonalize_quad_hamiltonian
   implicit none

   integer, parameter :: na=1, hdim=2, nq_ext=1
   complex(dblprec) :: h_in(hdim,hdim), eig_vec(hdim,hdim)
   complex(dblprec) :: S_prime(hdim,hdim,3,3,nq_ext)
   real(dblprec) :: eig_val(hdim), para_err

   call setup_diamag()
   h_in=(0.0_dblprec,0.0_dblprec)
   h_in(1,1)=(-1.0_dblprec,0.0_dblprec)
   h_in(2,2)=(-1.0_dblprec,0.0_dblprec)
   S_prime=(0.0_dblprec,0.0_dblprec)

   call diagonalize_quad_hamiltonian(na,h_in,eig_val,eig_vec,1,nq_ext,S_prime, &
      .false.,para_err,.true.)

   if (para_err>1.0e-10_dblprec) error stop 'fallback test: paraunitarity failure'
   if (eig_val(1)<=0.0_dblprec .or. eig_val(2)>=0.0_dblprec) then
      error stop 'fallback test: wrong boson branch ordering'
   end if
end program test_diagonalizer_fallback
