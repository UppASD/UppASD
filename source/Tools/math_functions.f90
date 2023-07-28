!------------------------------------------------------------------------------------
! MODULE: math_functions
!> @brief Collection of mathematical functions and subroutines.
!> @details This contains a series of mathematical functions and subroutines which are
!> used in several places in the code, the idea is for the code to be more homogeneous
!> and to ensure that all the definitions are the same throughout the package.
!> @author Jonathan Chico, Anders Bergman, Lars Bergqvist, Johan Hellsvik
!------------------------------------------------------------------------------------
module math_functions

   use Parameters

   implicit none

   public

contains
  !!!  !---------------------------------------------------------------------------------
  !!!  ! FUNCTION: f_dist
  !!!  !> @brief Calculates the norm of the difference of two vectors
  !!!  !> @author Jonathan Chico
  !!!  !---------------------------------------------------------------------------------
  !!!  function f_dist(A_in,B_in) result(dist)

  !!!     implicit none

  !!!     real(dblprec), dimension(3), intent(in) :: A_in
  !!!     real(dblprec), dimension(3), intent(in) :: B_in
  !!!     real(dblprec) :: dist

  !!!     dist = norm2(A_in-B_in)

  !!!  end function f_dist

   !---------------------------------------------------------------------------------
   ! FUNCTION: f_volume
   !> @brief Calculates the cross product between two 3D vectors
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   function f_cross_product(A_in,B_in) result(AcrossB)

      implicit none

      real(dblprec), dimension(3), intent(in) :: A_in
      real(dblprec), dimension(3), intent(in) :: B_in
      real(dblprec), dimension(3) :: AcrossB

      AcrossB(1)=A_in(2)*B_in(3)-A_in(3)*B_in(2)
      AcrossB(2)=A_in(3)*B_in(1)-A_in(1)*B_in(3)
      AcrossB(3)=A_in(1)*B_in(2)-A_in(2)*B_in(1)

   end function f_cross_product

   !---------------------------------------------------------------------------------
   ! FUNCTION: f_volume
   !> @brief Calculates the volume apnned by 3 vectors
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   function f_volume(A_in,B_in,C_in) result(volume)

      implicit none

      real(dblprec), dimension(3), intent(in) :: A_in
      real(dblprec), dimension(3), intent(in) :: B_in
      real(dblprec), dimension(3), intent(in) :: C_in
      real(dblprec) :: volume

      real(dblprec), dimension(3) :: AcrossB

      AcrossB=f_cross_product(A_in,B_in)

      volume = sum(C_in(:)*AcrossB(:))

   end function f_volume

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: f_reciprocal_lattice
   !> @brief Calculation of the reciprocal lattice vectors
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine f_reciprocal_lattice(C1,C2,C3,B1,B2,B3)

      implicit none

      real(dblprec), dimension(3), intent(in) :: C1
      real(dblprec), dimension(3), intent(in) :: C2
      real(dblprec), dimension(3), intent(in) :: C3
      real(dblprec), dimension(3), intent(out) :: B1
      real(dblprec), dimension(3), intent(out) :: B2
      real(dblprec), dimension(3), intent(out) :: B3

      real(dblprec) :: cell_vol

      cell_vol = f_volume(C1,C2,C3)

      B1= f_cross_product(C2,C3)/cell_vol
      B2= f_cross_product(C3,C1)/cell_vol
      B3= f_cross_product(C1,C2)/cell_vol

   end subroutine f_reciprocal_lattice

   !---------------------------------------------------------------------------------
   ! FUNCTION: f_normalize_vec
   !> @brief Returns a normalized vector
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   function f_normalize_vec(vec,vec_length) result(norm_vec)

      implicit none

      integer, intent(in) :: vec_length
      real(dblprec), dimension(vec_length), intent(in) :: vec

      real(dblprec), dimension(vec_length) :: norm_vec

      norm_vec(:) = vec(:)/norm2(vec)

   end function f_normalize_vec

   !-----------------------------------------------------------------------------
   !  FUNCTION: f_logfun
   !> @brief
   !> Function to write the measurements in logarithmic scale
   !-----------------------------------------------------------------------------
   integer function f_logfun(mstep)
      implicit none
      !
      integer, intent(in) :: mstep !< Current simulation step
      !
      integer :: sstep,pstep,tstep
      !
      if(mstep==0) then
         sstep=0
      else
         tstep=int(log(real(mstep))*10) 
         pstep=int(log(real(mstep-1))*10)
         if(tstep>pstep) then
            sstep=tstep
         else
            sstep=-1
         end if
      end if
      f_logfun=sstep
   end function f_logfun

   !--------------------------------------------!
   ! @date 2014/09/01 - Thomas Nystrand
   ! - Moved to separate function
   !--------------------------------------------!
   function f_logstep(mstep,logsamp) result(sstep)
      character(len=1)   :: logsamp
      integer,intent(in) :: mstep

      integer            :: sstep

      if(logsamp=='Y') then
         sstep=f_logfun(mstep)
         if(sstep<0) return
      else
         sstep=mstep
      end if
   end function f_logstep


   
   !  Subroutines moved from module Neighbourmap
   
   !  Used ???
   !  Multiplication of square matrices
   subroutine matmatmul(N, A, B, C)

      implicit none

      integer, intent(in) :: N !< Dimension of matrices
      real(dblprec), dimension(N,N),intent(in) :: A  !< first matrix
      real(dblprec), dimension(N,N),intent(in) :: B  !< second matrix
      real(dblprec), dimension(N,N),intent(out) :: C  !< the product of the two matrices

      !.. Scalar variables
      integer :: i,j,k

      C = 0.0_dblprec

      do i=1,N
         do j=1,N
            do k=1,N
               C(i,j) = C(i,j) + A(i,k) * B(k,j)
            end do
         end do
      end do

   end subroutine matmatmul

   !  Used ???
   !  Multiplication of a square matrices with a vector
   subroutine matvecmul(N, A, B, C)

      implicit none

      integer, intent(in) :: N !< Dimension of matrices
      real(dblprec), dimension(N,N),intent(in) :: A  !< matrix
      real(dblprec), dimension(N),intent(in) :: B  !< vector
      real(dblprec), dimension(N),intent(out) :: C  !< the product of the matrix and the vector

      !.. Scalar variables
      integer :: i,j

      C = 0.0_dblprec

      do i=1,N
         do j=1,N
            C(i) = C(i) + A(i,j) * B(j)
         end do
      end do

   end subroutine matvecmul

   ! Inversion of 3 x 3 matrix with Cramer's rule
   subroutine matinvcramer3(A, B)

      implicit none

      real(dblprec), dimension(3,3),intent(in) :: A  !< Input matrix
      real(dblprec), dimension(3,3),intent(out) :: B  !< The inverse of the input matrix

      !.. Scalar variables
      real(dblprec), dimension(3) :: a1, a2, a3
      real(dblprec), dimension(3) :: b1, b2, b3
      real(dblprec) :: detA, invdetA

      a1(1:3) = A(1:3,1)
      a2(1:3) = A(1:3,2)
      a3(1:3) = A(1:3,3)

      b1(1:3) = crossproduct(a2,a3)
      b2(1:3) = crossproduct(a3,a1)
      b3(1:3) = crossproduct(a1,a2)

      detA = dot_product(a1,b1)
      invdetA = 1.0_dblprec / detA

      B(1,1:3) = invdetA * b1(1:3)
      B(2,1:3) = invdetA * b2(1:3)
      B(3,1:3) = invdetA * b3(1:3)

   end subroutine matinvcramer3

   ! Cross product of two 3-vectors A and B
   function crossproduct(A, B) result(C)

      implicit none

      real(dblprec), dimension(3),intent(in) :: A  !< Left factor
      real(dblprec), dimension(3),intent(in) :: B  !< Right factor

      real(dblprec), dimension(3) :: C  !< The cross product of A and B

      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)

   end function crossproduct

   ! The determinant of a  3 x 3 matrix A
   function det3(A) result(C)

      real(dblprec), dimension(3,3),intent(in) :: A  !< Input matrix

      real(dblprec) :: C  !< The determinant of A


      !.. Scalar variables
      real(dblprec), dimension(3) :: a1, a2, a3
      real(dblprec), dimension(3) :: b1

      a1(1:3) = A(1:3,1)
      a2(1:3) = A(1:3,2)
      a3(1:3) = A(1:3,3)

      b1(1:3) = crossproduct(a2,a3)

      C = dot_product(a1,b1)

   end function det3


   ! Transformation of a rank-1 tensor (vector). Passive rotation
   function transt1(A, R) result(B)

      implicit none

      real(dblprec), dimension(3),intent(in) :: A  !< input tensor
      real(dblprec), dimension(3,3),intent(in) :: R  !< symmetry element
      real(dblprec), dimension(3) :: B  !< output tensor

      !.. Scalar variables
      integer :: i
      integer :: p

      B = 0.0_dblprec
      do i=1,3
         do p=1,3
            B(i) = B(i) + R(i,p) * A(p)
         end do
      end do

   end function transt1

   ! Transformation of a rank-2 tensor.  Passive rotation
   function transt2(Avec, R) result(Bvec)

      real(dblprec), dimension(9),intent(in) :: Avec  !< input tensor (flattened)
      real(dblprec), dimension(3,3),intent(in) :: R  !< symmetry element
      real(dblprec), dimension(9) :: Bvec  !< output tensor (flattened)

      ! Scalar variables
      integer :: i,j
      integer :: p,q

      ! Array variables
      real(dblprec), dimension(3,3) :: A, B

      A = reshape( Avec, (/3,3/) )
      B = 0.0_dblprec

      do i=1,3
         do j=1,3
            do p=1,3
               do q=1,3
                  B(i,j) = B(i,j) + R(i,p) * R(j,q) * A(p,q)
               end do
            end do
         end do
      end do

      Bvec = reshape( B, (/9/) )

   end function transt2

   ! Transformation of a rank-3 tensor.  Passive rotation
   function transt3(Avec, T) result(Bvec)

      real(dblprec), dimension(27),intent(in) :: Avec  !< input tensor (flattened)
      real(dblprec), dimension(3,3),intent(in) :: T  !< symmetry element
      real(dblprec), dimension(27) :: Bvec  !< output tensor (flattened)

      ! Scalar variables
      integer :: i,j,k
      integer :: p,q,r

      ! Array variables
      real(dblprec), dimension(3,3,3) :: A, B

      real(dblprec), dimension(3,3) :: Tinv

      A = reshape( Avec, (/3,3,3/) )
      B = 0.0_dblprec
      Tinv = -T

      do i=1,3
         do j=1,3
            do k=1,3
               do p=1,3
                  do q=1,3
                     do r=1,3
                        B(i,j,k) = B(i,j,k) + T(i,p) * T(j,q) * T(k,r) * A(p,q,r)
                     end do
                  end do
               end do
            end do
         end do
      end do

      Bvec = reshape( B, (/27/) )

   end function transt3

   ! Transformation of a rank-4 tensor.  Passive rotation
   function transt4(Avec, T) result(Bvec)

      real(dblprec), dimension(81),intent(in) :: Avec  !< input tensor (flattened)
      real(dblprec), dimension(3,3),intent(in) :: T  !< symmetry element
      real(dblprec), dimension(81) :: Bvec  !< output tensor (flattened)

      ! Scalar variables
      integer :: i,j,k,l
      integer :: p,q,r,s

      ! Array variables
      real(dblprec), dimension(3,3,3,3) :: A, B

      A = reshape( Avec, (/3,3,3,3/) )
      B = 0.0_dblprec

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do p=1,3
                     do q=1,3
                        do r=1,3
                           do s=1,3
                              B(i,j,k,l) = B(i,j,k,l) + T(i,p) * T(j,q) * T(k,r) * T(l,s) * A(p,q,r,s)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      Bvec = reshape( B, (/81/) )

   end function transt4


   !------------------------------------------------------------------------------------
   ! SUBROUTINE: gramms
   !> @brief Gramm-Schmidt orthogonalization
   !------------------------------------------------------------------------------------
   subroutine gramms(v_in,m_out,ldim)
      !
      integer, intent(in) :: ldim
      real(dblprec), dimension(3,ldim), intent(in) :: v_in
      real(dblprec), dimension(3,3,ldim), intent(out) :: m_out
      !
      real(dblprec), dimension(3) :: v_trial
      real(dblprec)  :: v_in_norm
      integer :: i
      !
      do i=1,ldim
         v_in_norm=sqrt(v_in(1,i)*v_in(1,i)+v_in(2,i)*v_in(2,i)+v_in(3,i)*v_in(3,i))
         ! Set trial vector to x
         v_trial=(/1.0_dblprec, 0.0_dblprec, 0.0_dblprec/)
         v_trial=v_trial/sqrt(sum(v_trial*v_trial))
         ! If input vector is close to x, set trial vector to y..
         if(abs(v_in(1,i))/v_in_norm>0.99995_dblprec) then
            v_trial=(/0.0_dblprec, 1.0_dblprec, 0.0_dblprec/)
            v_trial=v_trial/sqrt(sum(v_trial*v_trial))
         end if
         ! Keep input vector as third output vector
         m_out(1:3,3,i)=v_in(1:3,i)/v_in_norm
         ! Get first orthogonal vector through cross product
         m_out(1,1,i)=m_out(2,3,i)*v_trial(3)-m_out(3,3,i)*v_trial(2)
         m_out(2,1,i)=m_out(3,3,i)*v_trial(1)-m_out(1,3,i)*v_trial(3)
         m_out(3,1,i)=m_out(1,3,i)*v_trial(2)-m_out(2,3,i)*v_trial(1)
         ! Normalize
         m_out(:,1,i)=m_out(:,1,i) / sqrt(sum(m_out(:,1,i)*m_out(:,1,i)))
         ! Get second orthogonal vector by cross product
         m_out(1,2,i)=m_out(2,3,i)*m_out(3,1,i)-m_out(3,3,i)*m_out(2,1,i)
         m_out(2,2,i)=m_out(3,3,i)*m_out(1,1,i)-m_out(1,3,i)*m_out(3,1,i)
         m_out(3,2,i)=m_out(1,3,i)*m_out(2,1,i)-m_out(2,3,i)*m_out(1,1,i)
         ! Normalize for safety..
         m_out(:,2,i)=m_out(:,2,i) / sqrt(sum(m_out(:,2,i)*m_out(:,2,i)))
         ! Get second orthogonal vector by cross product
         m_out(1,1,i)=m_out(2,2,i)*m_out(3,3,i)-m_out(3,2,i)*m_out(2,3,i)
         m_out(2,1,i)=m_out(3,2,i)*m_out(1,3,i)-m_out(1,2,i)*m_out(3,3,i)
         m_out(3,1,i)=m_out(1,2,i)*m_out(2,3,i)-m_out(2,2,i)*m_out(1,3,i)
         ! Normalize for safety..
         m_out(:,1,i)=m_out(:,1,i) / sqrt(sum(m_out(:,1,i)*m_out(:,1,i)))
      end do
      !
   end subroutine gramms


   !------------------------------------------------------------------------------------
   ! SUBROUTINE: gramms2
   !> @brief Gramm-Schmidt orthogonalization alternative
   !------------------------------------------------------------------------------------
   subroutine gramms2(v_in,v_qt,m_out,ldim)
      !
      integer, intent(in) :: ldim
      real(dblprec), dimension(3,ldim), intent(in) :: v_in
      real(dblprec), dimension(3,ldim), intent(in) :: v_qt
      real(dblprec), dimension(3,ldim), intent(out) :: m_out
      !
      real(dblprec), dimension(3) :: v_perp1, v_perp2
      integer :: i
      !
      do i=1,ldim
         ! Get first orthogonal vector by cross product
         v_perp1(1)=v_qt(2,i)*v_in(3,i)-v_qt(3,i)*v_in(2,i)
         v_perp1(2)=v_qt(3,i)*v_in(1,i)-v_qt(2,i)*v_in(3,i)
         v_perp1(3)=v_qt(1,i)*v_in(2,i)-v_qt(1,i)*v_in(1,i)
         ! Get second orthogonal vector by cross product
         v_perp2(1)=v_perp1(2)*v_in(3,i)-v_perp1(3)*v_in(2,i)
         v_perp2(2)=v_perp1(3)*v_in(1,i)-v_perp1(2)*v_in(3,i)
         v_perp2(3)=v_perp1(1)*v_in(2,i)-v_perp1(1)*v_in(1,i)
         ! Calculate projections
         m_out(1,i)=sum(v_perp1*v_in(:,i))
         m_out(2,i)=sum(v_perp2*v_in(:,i))
         m_out(3,i)=sum(v_qt(:,i)*v_in(:,i))
      end do
      !
   end subroutine gramms2

   function f_interp_search(key,arr,n) result(i)
      !
      implicit none
      !
      real(dblprec), intent(in) :: key
      integer, intent(in) :: n
      real(dblprec), dimension(n) :: arr
      !
      integer :: i
      integer :: idum
      integer :: low, high
      !
      low=1;high=n
      if (key<=arr(low)) then
         i=low
         return
      end if
      if (key>=arr(high)) then
         i=high
         return
      end if

      do while (low<n)
         idum=int(low+(((key-arr(low))*(high-low))/(arr(high)-arr(low))))
         if (key>=arr(idum).and.key<arr(idum+1)) then
             i = idum
             return
          end if

         if (key<arr(idum)) then
            high=idum-1
         else
            low=idum+1
         end if

      end do

      i=idum
      return


   end function f_interp_search


   function f_binary_search(key,arr,n) result(i)
      !
      implicit none
      !
      real(dblprec), intent(in) :: key
      integer, intent(in) :: n
      real(dblprec), dimension(n) :: arr
      !
      integer :: i
      integer :: idum
      integer :: low, high
      !
      low=1;high=n
      if (key<=arr(low)) then
         i=low
         return
      end if
      if (key>=arr(high)) then
         i=high
         return
      end if

      do while (low<n)
         idum=(low+high)/2
         if (key>=arr(idum).and.key<arr(idum+1)) then
             i = idum
             return
          end if

         if (key<arr(idum)) then
            high=idum
         else
            low=idum
         end if

      end do

      i=idum
      return


   end function f_binary_search


   function f_joint_search(key,arr,n) result(i)
      !
      implicit none
      !
      real(dblprec), intent(in) :: key
      integer, intent(in) :: n
      real(dblprec), dimension(n) :: arr
      !
      integer :: i
      integer :: idum
      integer :: low, high
      !
      integer :: safety
      !
      low=1;high=n
      if (key<=arr(low)) then
         i=low
         return
      end if
      if (key>=arr(high)) then
         i=high
         return
      end if

      !First one interpolation step
      idum=int(low+(((key-arr(low))*(high-low))/(arr(high)-arr(low))))

      if (key>=arr(idum).and.key<arr(idum+1)) then
         i = idum
         return
      end if

      if (key<arr(idum)) then
         high=idum-1
      else
         low=idum+1
      end if

      !Then continue with binary search
      safety=0
      do while (low<n.and.safety<2*n)
         safety=safety+1
         idum=(low+high)/2
         if (key>=arr(idum).and.key<arr(idum+1)) then
            i = idum
            return
         end if

         if (key<arr(idum)) then
            high=idum
         else
            low=idum
         end if

      end do

      i=idum
      return


   end function f_joint_search

   function f_interp_1d(x_val,x_arr,y_arr,x_len) result(y_val)
      !
      !
      implicit none
      !
      real(dblprec), intent(in) :: x_val
      integer, intent(in) :: x_len
      real(dblprec), dimension(x_len), intent(in) :: x_arr
      real(dblprec), dimension(x_len), intent(in) :: y_arr
      !
      real(dblprec) :: y_val
      !
      real(dblprec) :: y1,y2,x1,x2,frac
      integer :: xi1, xi2 

      if (x_val <= x_arr(1)) then
         y_val = y_arr(1)
      else if (x_val >= x_arr(x_len)) then
         y_val = y_arr(x_len)
      else

         ! A full search is slower but safer
         xi1 = 1
         do while (x_arr(xi1)<x_val)
            xi1 = xi1 + 1
         end do
         ! The binary search sometimes get stuck
         !xi1 = f_joint_search(x_val,x_arr,x_len)

         xi2 = xi1 + 1

         x1 = x_arr(xi1)
         x2 = x_arr(xi2)

         y1 = y_arr(xi1)
         y2 = y_arr(xi2)

         frac = (x_val - x1) / (x2 - x1)
         y_val = y1 + frac *  (y2-y1)
      end if

      return

   end function f_interp_1d



subroutine f_wrap_coord_diff(Natom,coord,i_atom,j_atom,cdiff)
   !
   use InputData, only : BC1,BC2,BC3,C1,C2,C3,N1,N2,N3
   !
   implicit none
   !
   integer, intent(in) :: Natom
   real(dblprec), dimension(3,Natom), intent(in) :: coord
   integer, intent(in) :: i_atom
   integer, intent(in) :: j_atom
   real(dblprec), dimension(3), intent(out) :: cdiff
   !
   real(dblprec), dimension(3) :: odiff, oshift, mdiff
   integer :: x,y,z
   integer :: xmin,xmax,ymin,ymax,zmin,zmax
   !
   odiff=coord(:,j_atom) - coord(:,i_atom)
   !odiff=coord(:,i_atom) - coord(:,j_atom)
   !onorm=norm2(odiff)
   !
   xmax=0;xmin=0;ymax=0;ymin=0;zmax=0;zmin=0
   if(BC1=='P') then
      xmax=1
      xmin=-1
   end if
   if(BC2=='P') then
      ymax=1
      ymin=-1
   end if
   if(BC3=='P') then
      zmax=1
      zmin=-1
   end if
   !
   mdiff=odiff
   do z=zmin,zmax
      do y=ymin,ymax
         do x=xmin,xmax
            oshift = odiff + x*(N1)*C1 + y*(N2)*C2 + z*(N3)*C3
            !oshift = odiff + x*(N1-1)*C1 + y*(N2-1)*C2 + z*(N3-1)*C3
            if(norm2(oshift)<norm2(mdiff))  mdiff = oshift
         end do
      end do
   end do
   !
   !print '(2i6,6f10.6)',i_atom,j_atom,mdiff!, oshift
   cdiff=mdiff
   return
   !
end subroutine f_wrap_coord_diff

subroutine f_get_periodic_shifts(nshifts,shift_arr)
!subroutine f_wrap_coord_diff(Natom,coord,i_atom,j_atom,cdiff)
   !
   use InputData, only : BC1,BC2,BC3,C1,C2,C3,N1,N2,N3
   !
   implicit none
   !
   integer, intent(out) :: nshifts
   real(dblprec), dimension(3,27), intent(out) :: shift_arr
   !
   integer :: x,y,z
   integer :: xmin,xmax,ymin,ymax,zmin,zmax
   !
   xmax=0;xmin=0;ymax=0;ymin=0;zmax=0;zmin=0
   if(BC1=='P') then
      xmax=1
      xmin=-1
   end if
   if(BC2=='P') then
      ymax=1
      ymin=-1
   end if
   if(BC3=='P') then
      zmax=1
      zmin=-1
   end if
   !
   nshifts=0
   do z=zmin,zmax
      do y=ymin,ymax
         do x=xmin,xmax
            nshifts=nshifts+1
            shift_arr(:,nshifts) =  x*(N1)*C1 + y*(N2)*C2 + z*(N3)*C3
         end do
      end do
   end do
   !
   return
   !
end subroutine f_get_periodic_shifts

!---------------------------------------------------------------------------------
! FUNCTION: f_norms
!> @brief Calculate norms for 2d-array
!> @author Anders Bergman
!---------------------------------------------------------------------------------
function f_norms(dim_1, dim_2, array) result (norms)
      implicit none


      integer, intent(in) :: dim_1        !< Dimension to take the norm over
      integer, intent(in) :: dim_2        !< Dimension of entries 
      real(dblprec), dimension(dim_1, dim_2), intent(in) :: array
      real(dblprec), dimension(dim_2) :: norms

      integer :: i_idx

      !$omp parallel do default(shared) private(i_idx)
      do i_idx = 1, dim_2
            norms(i_idx) = sum(array(:, i_idx) * array(:, i_idx))
            norms(i_idx) = sqrt(norms(i_idx))
      end do
      !$omp end parallel do

end function f_norms

!---------------------------------------------------------------------------------
! FUNCTION: f_max_norm
!> @brief Finds largest norm in array of vectors
!> @author Anders Bergman
!---------------------------------------------------------------------------------
function f_max_norm(dim_1, dim_2, array) result (m_norm)
      implicit none


      integer, intent(in) :: dim_1        !< Dimension to take the norm over
      integer, intent(in) :: dim_2        !< Dimension of entries 
      real(dblprec), dimension(dim_1, dim_2), intent(in) :: array
      real(dblprec) :: m_norm

      integer :: i_idx
      real(dblprec) :: t_norm
      real(dblprec), dimension(dim_2) :: t_array

      do i_idx = 1, dim_2
            t_array(i_idx) = sum(array(:, i_idx) * array(:, i_idx))
      end do

      m_norm = maxval(t_array)
      m_norm = sqrt(m_norm)


end function f_max_norm

end module math_functions
