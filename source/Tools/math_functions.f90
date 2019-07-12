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

   
end module math_functions
