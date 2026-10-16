! SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
! SPDX-FileCopyrightText: 2026 UppASD contributors
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Reference model of the UppASD multiscale interface interpolation.
!! atomInterpolation and the body of multiscaleInterpolateInterfaces are
!! transcribed from thirdparty/UppASD/source/Multiscale/multiscaleinterpolation.f90
!! (authors Edgar Mendez, Nikos Ntallis, Manuel Pereiro).
!! The transcription replaces the MomentData and interfaceInterpolation module
!! globals with dummy arguments.
!! Deviates from UppASD: neighbours are read from copies of emom and emom2 to avoid
!! a race between interpolated atoms.
module multiscale_model
   use iso_c_binding, only: c_int, c_double
   implicit none
   private
   public :: model_multiscale_interpolate_interfaces

   integer, parameter :: dblprec = c_double

   type LocalInterpolationInfo
      integer(c_int) :: nrInterpAtoms
      integer(c_int), dimension(:), pointer :: indices
      integer(c_int), dimension(:), pointer :: firstNeighbour
      real(dblprec), dimension(:), pointer :: weights
      integer(c_int), dimension(:), pointer :: neighbours
   end type LocalInterpolationInfo

contains

   function atomInterpolation(interp,atom,ensemble,emom) result(v)
      implicit none
      type(LocalInterpolationInfo), intent(in) :: interp
      integer, intent(in) :: atom
      integer, intent(in) :: ensemble
      real(dblprec), dimension(:,:,:), intent(in) :: emom
      real(dblprec), dimension(3) :: v

      integer :: i
      integer :: index

      if (associated(interp%indices)) then
         index = interp%indices(atom)
         if (index .ne. 0) then
            v = (/ 0, 0, 0 /)
            do i=interp%firstNeighbour(index), &
                 interp%firstNeighbour(index+1)-1
               v = v + emom(:,interp%neighbours(i),ensemble) * interp%weights(i)
            end do
            v = v / sqrt(sum(v**2))
         else
            v = emom(:,atom,ensemble)
         end if
      end if
   end function atomInterpolation

   !> Applies multiscaleInterpolateInterfaces in place to caller-owned arrays.
   !! All index arrays are 1-based and laid out as in LocalInterpolationInfo.
   !! @param[in] atom_count Number of atoms (Natom)
   !! @param[in] ensemble_count Number of ensembles (Mensemble)
   !! @param[in] row_count Number of interpolated atoms (nrInterpAtoms)
   !! @param[in] weight_count Number of weights and neighbours
   !! @param[in] indices Row of each atom, 0 if not interpolated. Shape (atom_count)
   !! @param[in] first_neighbour Row offsets into weights and neighbours. Shape (row_count+1)
   !! @param[in] neighbours Neighbour atom of each weight. Shape (weight_count)
   !! @param[in] weights Interpolation weights. Shape (weight_count)
   !! @param[in] mmom Moment magnitudes. Shape (atom_count, ensemble_count)
   !! @param[inout] emom Unit moments. Shape (3, atom_count, ensemble_count)
   !! @param[inout] emom2 Unit moments. Shape (3, atom_count, ensemble_count)
   !! @param[inout] emomM Moments. Written only at interpolated atoms. Shape (3, atom_count, ensemble_count)
   subroutine model_multiscale_interpolate_interfaces(atom_count, ensemble_count, &
         row_count, weight_count, indices, first_neighbour, neighbours, weights, &
         mmom, emom, emom2, emomM) bind(C, name="model_multiscale_interpolate_interfaces")
      implicit none
      integer(c_int), value :: atom_count, ensemble_count, row_count, weight_count
      integer(c_int), intent(in), target :: indices(atom_count)
      integer(c_int), intent(in), target :: first_neighbour(row_count+1)
      integer(c_int), intent(in), target :: neighbours(weight_count)
      real(dblprec), intent(in), target :: weights(weight_count)
      real(dblprec), intent(in) :: mmom(atom_count, ensemble_count)
      real(dblprec), intent(inout) :: emom(3, atom_count, ensemble_count)
      real(dblprec), intent(inout) :: emom2(3, atom_count, ensemble_count)
      real(dblprec), intent(inout) :: emomM(3, atom_count, ensemble_count)

      type(LocalInterpolationInfo) :: interfaceInterpolation
      integer :: atom,index,ens
      real(dblprec), allocatable :: emom_in(:,:,:), emom2_in(:,:,:)

      emom_in = emom
      emom2_in = emom2

      interfaceInterpolation%nrInterpAtoms = row_count
      interfaceInterpolation%indices => indices
      interfaceInterpolation%firstNeighbour => first_neighbour
      interfaceInterpolation%weights => weights
      interfaceInterpolation%neighbours => neighbours

      !$omp parallel do private(atom,index,ens)
      do atom=1,ubound(interfaceInterpolation%indices, 1)
         index = interfaceInterpolation%indices(atom)
         if (index .ne. 0) then
            do ens=1,ubound(emom,3)
               emom(:,atom,ens) = &
                    atomInterpolation(interfaceInterpolation,atom,ens,emom_in)
               emom2(:,atom,ens) = &
                    atomInterpolation(interfaceInterpolation,atom,ens,emom2_in)
               emomM(:,atom,ens) = emom(:,atom,ens) * mmom(atom,ens)
            end do
         end if
      end do
      !$omp end parallel do
   end subroutine model_multiscale_interpolate_interfaces

end module multiscale_model
