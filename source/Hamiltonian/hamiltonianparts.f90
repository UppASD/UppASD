!> Calculate effective field by applying the derivative of the Hamiltonian
!> \details The effective field, \f$\mathbf{B}_i\f$, on an atom \f$\textit{i}\f$, is calculated from 
!> \f$ \mathbf{B}_i=-\frac{\partial \mathbf{H}}{\partial \mathbf{m}_i},\f$ where primarily the part of 
!> the Hamiltonian, \f$\mathbf{H}\f$, which represents interatomic exchange interactions, 
!> \f$\mathbf{H}_\mathrm{ex}\f$, are considered. For this we use the classical Heisenberg Hamiltonian, 
!> \f$ \mathbf{H}_\mathrm{ex}=-\frac{1}{2}\sum_{i\neq j}J_{ij}\mathbf{m}_i\cdot\mathbf{m}_j,\f$ where 
!> \f$i\f$ and \f$j\f$ are atomic indices and \f$J_{ij}\f$ is the strength of the exchange interaction, 
!> which is calculated from first principles theory.
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module HamiltonianParts
   use Parameters
   use Profiling
   use HamiltonianData
   use MomentData, only : emomM

   implicit none


   private


contains


   subroutine heisenberg_field(Natom, Mensemble, i, k, field, emomM)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      !
      integer :: j

      !Exchange term
      do j=1,nlistsize(i)
         field = field + ncoup(j,i)*emomM(1:3,nlist(j,i),k) 
      end do
   end subroutine heisenberg_field


   subroutine dzyaloshinskii_moriya_field(Natom, Mensemble, i, k, field)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      !
      integer :: j

      ! Dzyaloshinskii_moriya term
      do j=1,dmlistsize(i)
         field(1) = field(1) - dm_vect(3,j,i)*emomM(2,dmlist(j,i),k) +&
            dm_vect(2,j,i)*emomM(3,dmlist(j,i),k)
         field(2) = field(2) - dm_vect(1,j,i)*emomM(3,dmlist(j,i),k) +&
            dm_vect(3,j,i)*emomM(1,dmlist(j,i),k)
         field(3) = field(3) - dm_vect(2,j,i)*emomM(1,dmlist(j,i),k) +&
            dm_vect(1,j,i)*emomM(2,dmlist(j,i),k)
      end do

   end subroutine dzyaloshinskii_moriya_field

   subroutine pseudo_dipolar_field(Natom, Mensemble, i, k, field)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      !
      integer :: j

      ! Pseudo-Dipolar term
      do j=1,pdlistsize(i)
         field(1) = field(1) + pd_vect(1,j,i)*emomM(1,pdlist(j,i),k) +&
            pd_vect(4,j,i)*emomM(2,pdlist(j,i),k) +&
            pd_vect(5,j,i)*emomM(3,pdlist(j,i),k)
         field(2) = field(2) + pd_vect(4,j,i)*emomM(1,pdlist(j,i),k) +&
            pd_vect(2,j,i)*emomM(2,pdlist(j,i),k) +&
            pd_vect(6,j,i)*emomM(3,pdlist(j,i),k)
         field(3) = field(3) + pd_vect(5,j,i)*emomM(1,pdlist(j,i),k) +&
            pd_vect(6,j,i)*emomM(2,pdlist(j,i),k) +&
            pd_vect(3,j,i)*emomM(3,pdlist(j,i),k)
      end do

   end subroutine pseudo_dipolar_field

   subroutine dzyaloshinskii_moriya_bq_field(Natom, Mensemble, i, k, field)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      !
      integer :: j
      real(dblprec), dimension(3) :: dot !< Work array

      ! BIQDM term
      do j=1,biqdmlistsize(i)
         dot(1) = emomM(1,i,k)*emomM(2,biqdmlist(j,i),k)-&
            emomM(2,i,k)*emomM(1,biqdmlist(j,i),k)
         dot(2) = emomM(2,i,k)*emomM(3,biqdmlist(j,i),k)-&
            emomM(3,i,k)*emomM(2,biqdmlist(j,i),k)
         dot(3) = emomM(3,i,k)*emomM(1,biqdmlist(j,i),k)-&
            emomM(1,i,k)*emomM(3,biqdmlist(j,i),k)
         field(1) = field(1) + 2.0d0*biqdm_vect(1,j,i)*(&
            dot(1)*emomM(3,biqdmlist(j,i),k)-&
            dot(2)*emomM(2,biqdmlist(j,i),k))
         field(2) = field(2) + 2.0d0*biqdm_vect(1,j,i)*(&
            dot(2)*emomM(1,biqdmlist(j,i),k)-&
            dot(3)*emomM(3,biqdmlist(j,i),k))
         field(3) = field(3) + 2.0d0*biqdm_vect(1,j,i)*(&
            dot(3)*emomM(2,biqdmlist(j,i),k)-&
            dot(1)*emomM(1,biqdmlist(j,i),k))
      end do
   end subroutine dzyaloshinskii_moriya_bq_field

   subroutine biquadratic_field(Natom, Mensemble, i, k, field)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      !
      integer :: j
      real(dblprec) :: dot

      ! Biquadratic exchange term 
      do j=1,bqlistsize(i)
         dot=emomM(1,bqlist(j,i),k)*emomM(1,i,k)+&
            emomM(2,bqlist(j,i),k)*emomM(2,i,k)+&
            emomM(3,bqlist(j,i),k)*emomM(3,i,k)
         field = field + 2.0d0*j_bq(j,i)*dot*emomM(1:3,bqlist(j,i),k)
      end do
   end subroutine biquadratic_field

   subroutine dipolar_field(Natom, Mensemble, i, k, field)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      !
      integer :: j
      real(dblprec) :: dot

      ! Dipolar term
      do j=1,Natom
         field(1) = field(1) + Qdip(1,1,j,i)*emomM(1,j,k) + Qdip(2,1,j,i)*emomM(2,j,k) + Qdip(3,1,j,i)*emomM(3,j,k)
         field(2) = field(2) + Qdip(1,2,j,i)*emomM(1,j,k) + Qdip(2,2,j,i)*emomM(2,j,k) + Qdip(3,2,j,i)*emomM(3,j,k)
         field(3) = field(3) + Qdip(1,3,j,i)*emomM(1,j,k) + Qdip(2,3,j,i)*emomM(2,j,k) + Qdip(3,3,j,i)*emomM(3,j,k)
      end do
   end subroutine dipolar_field

   subroutine uniaxial_anisotropy_field(Natom, Mensemble, i, k, field,emomM)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      !
      real(dblprec) :: tt1,tt2,tt3

      ! Uniaxial anisotropy
      ! cos(theta)
      tt1=emomM(1,i,k)*eaniso(1,i)+emomM(2,i,k)*eaniso(2,i)+emomM(3,i,k)*eaniso(3,i)
      ! k1 + 2*k2*sin^2(theta) = k1 + 2*k2*(1-cos^2(theta))
      tt2=kaniso(1,i)+2.0d0*kaniso(2,i)*(1.0d0-tt1*tt1)
      ! 2 * cos(theta)* [k1 + 2*k2*sin^2(theta)]
      tt3= 2.0d0*tt1*tt2

      field  = field - tt3*eaniso(1:3,i)

   end subroutine uniaxial_anisotropy_field

   subroutine cubic_anisotropy_field(Natom, Mensemble, i, k, field)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      !

      ! Cubic anisotropy 
      field(1) = field(1)  & 
         + 2.0d0*kaniso(1,i)*emomM(1,i,k)*(emomM(2,i,k)**2+emomM(3,i,k)**2) &
         + 2.0d0*kaniso(2,i)+emomM(1,i,k)*emomM(2,i,k)**2*emomM(3,i,k)**2
      field(2) = field(2)  &
         + 2.0d0*kaniso(1,i)*emomM(2,i,k)*(emomM(3,i,k)**2+emomM(1,i,k)**2) &
         + 2.0d0*kaniso(2,i)+emomM(2,i,k)*emomM(3,i,k)**2*emomM(1,i,k)**2
      field(3) = field(3)  &
         + 2.0d0*kaniso(1,i)*emomM(3,i,k)*(emomM(1,i,k)**2+emomM(2,i,k)**2) &
         + 2.0d0*kaniso(2,i)+emomM(3,i,k)*emomM(1,i,k)**2*emomM(2,i,k)**2

   end subroutine cubic_anisotropy_field

   subroutine tensor_heisenberg_field(Natom, Mensemble, i, k, field)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: i !< Atom to calculate effective field for
      integer, intent(in) :: k !< Current ensemble
      real(dblprec), dimension(3), intent(inout) :: field !< Effective field
      !
      integer :: j

      !Exchange term
      do j=1,nlistsize(i)
         field(1) = field(1) + j_tens(1,1,j,i)*emomM(1,nlist(j,i),k) + j_tens(1,2,j,i)*emomM(2,nlist(j,i),k) & 
            + j_tens(1,3,j,i)*emomM(3,nlist(j,i),k)   
         field(2) = field(2) + j_tens(2,1,j,i)*emomM(1,nlist(j,i),k) + j_tens(2,2,j,i)*emomM(2,nlist(j,i),k) & 
            + j_tens(2,3,j,i)*emomM(3,nlist(j,i),k)   
         field(3) = field(3) + j_tens(3,1,j,i)*emomM(1,nlist(j,i),k) + j_tens(3,2,j,i)*emomM(2,nlist(j,i),k) & 
            + j_tens(3,3,j,i)*emomM(3,nlist(j,i),k)   
      end do
   end subroutine tensor_heisenberg_field


end module HamiltonianParts
