!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro
 module MultiscaleInterpolation
  use Parameters
  
  ! Here is an example of use of the LocalInterpolationInfo structure (serializes it)  
  ! do i=1,Natoms
  !    if(lii%indices(i) .ne. 0) then
  !       print *,"ATOM", i
  !       do j=lii%firstNeighbourIndex(lii%indices(i)),&
  !            lii%firstNeighbourIndex(lii%indices(i)+1)-1
  !          print *,"NEIGHBOUR", lii%neighbours(j), lii%weights(j)
  !       end do
  !    endif
  ! end do

  !! Represents the parameters required for an interpolation over a set of atoms of the form
  !! normalized(Σ(w_i m_i))  where w_i=weights, and i iterates over a set of neighbours
  type LocalInterpolationInfo
     integer ::  nrInterpAtoms !> Number of atoms affected
     !> index on firstNeighbour corresponding to each atom.
     !! indices(i) contains 0 if the atom is not affected by the interpolation.
     integer, dimension(:), pointer :: indices
     !> index of the first neighbour in weights and neighbours
     integer, dimension(:), pointer :: firstNeighbour
     !> Per-neighbour coefficient
     real(dblprec), dimension(:), pointer :: weights
     !> Atom indices for neighbours participating in the interpolation
     integer, dimension(:), pointer :: neighbours
  end type LocalInterpolationInfo

  type(LocalInterpolationInfo) interfaceInterpolation
contains
  
  !> Initializes a new LocalInterpolationInfo structure
  !! Don´t use uninitialized structures, as non-nullified pointers inside them
  !! will cause crashes when deallocating or running solvers
  subroutine newLocalInterpolationInfo(localInterp)
    implicit none
    type(LocalInterpolationInfo), intent(inout) :: localInterp
    
    localInterp%nrInterpAtoms = 0
    nullify(localInterp%indices)
    nullify(localInterp%firstNeighbour)
    nullify(localInterp%weights)
    nullify(localInterp%neighbours)
  end subroutine newLocalInterpolationInfo

  !! Releases all allocated pointers in a LocalInterpolationInfo structure
  !! and resets it. 
  subroutine deleteLocalInterpolationInfo(localInterp)
    use Profiling
  implicit none
    type(LocalInterpolationInfo), intent(inout) :: localInterp
    integer :: i_stat,var_size
    
    localInterp%nrInterpAtoms = 0
    if(associated(localInterp%indices)) then
       var_size = product(shape(localInterp%indices))*kind(localInterp%indices)
       deallocate(localInterp%indices, stat=i_stat)
       call memocc(i_stat, -var_size, &
            'LocalInterpolationInfo%indices','deleteDampingBand')
       nullify(localInterp%indices)
    end if
    if(associated(localInterp%firstNeighbour)) then
       var_size = product(shape(localInterp%firstNeighbour))*kind(localInterp%firstNeighbour)
       deallocate(localInterp%firstNeighbour, stat=i_stat)
       call memocc(i_stat, -var_size, &
            'LocalInterpolationInfo%firstNeighbour','deleteDampingBand')
       nullify(localInterp%firstNeighbour)
    end if
    if(associated(localInterp%weights)) then
       var_size = product(shape(localInterp%weights))*kind(localInterp%weights)
       deallocate(localInterp%weights, stat=i_stat)
       call memocc(i_stat, -var_size, &
            'LocalInterpolationInfo%weights','deleteDampingBand')
       nullify(localInterp%weights)
    end if
    if(associated(localInterp%neighbours)) then
       var_size = product(shape(localInterp%neighbours))*kind(localInterp%neighbours)
       deallocate(localInterp%neighbours, stat=i_stat)
       call memocc(i_stat, -var_size, &
            'LocalInterpolationInfo%neighbours','deleteDampingBand')
       nullify(localInterp%neighbours)
    end if
  end subroutine deleteLocalInterpolationInfo

    !> Evaluates the interpolation
  !! (Σ(w_i m_i))  where w_i=weights
  !! The norm is also interpolated
  function denormalInterpolation(interp,atom,ensemble,arr) result(v)
    implicit none
    type(LocalInterpolationInfo), intent(in) :: interp
    integer, intent(in) :: atom !< atom index
    integer, intent(in) :: ensemble !< current ensemble
    real(dblprec), dimension(:,:,:), intent(in) :: arr   !< Current unit moment vector
    real(dblprec), dimension(3) :: v,t
    real(dblprec) :: norm, div
    integer :: i
    integer :: index

    if (associated(interp%indices)) then
       index = interp%indices(atom)
       if (index .ne. 0) then
          v = (/ 0, 0, 0 /)
          div = 0
          norm = 0
          do i=interp%firstNeighbour(index), &               
               interp%firstNeighbour(index+1)-1
             t = arr(:,interp%neighbours(i),ensemble) * interp%weights(i)
             v = v + t 
             div = div + interp%weights(i)
             norm = norm + sqrt(sum(t**2))
          end do
          !n0 = sqrt(sum(v**2))
          !if (n0 > 1d-10) &
          !     v = (v / n0)! * (norm/div)
       else
          v = arr(:,atom,ensemble)
       end if
    end if
  end function denormalInterpolation

  
  !> Evaluates the spherical interpolation
  !! (Σ(w_i m_i))  where w_i=weights
  !! The norm is also interpolated
  function denormalSlerp(interp,atom,ensemble,arr) result(v)
    implicit none
    type(LocalInterpolationInfo), intent(in) :: interp
    integer, intent(in) :: atom !< atom index
    integer, intent(in) :: ensemble !< current ensemble
    real(dblprec), dimension(:,:,:), intent(in) :: arr   !< Current unit moment vector
    real(dblprec), dimension(3) :: v,cart,pol,acc
    real(dblprec) :: div
    integer :: i
    integer :: index

    if (associated(interp%indices)) then
       index = interp%indices(atom)
       if (index .ne. 0) then
          acc = (/ 0, 0, 0 /)
          div = 0
          do i=interp%firstNeighbour(index), &               
               interp%firstNeighbour(index+1)-1
             cart = arr(:,interp%neighbours(i),ensemble) 
             pol(1) = sqrt(sum(cart**2))
             pol(2) = atan2(cart(2),cart(1))
             pol(3) = acos(cart(3)/pol(1))             
             acc = acc + pol * interp%weights(i) 
             div = div + interp%weights(i) 
          end do
          acc = acc / div
          v(1) = acc(1) * cos(acc(2)) * sin(acc(3))
          v(2) = acc(1) * sin(acc(2)) * sin(acc(3))
          v(3) = acc(1) * cos(acc(3))
       else
          v = arr(:,atom,ensemble)
       end if
    end if
  end function denormalSlerp
  
  
  !> Evaluates the interpolation
  !! normalized(Σ(w_i m_i))  where w_i=weights
  !! and i iterates over the neighbours of atom defined by the interpolation
  !! if the atom is not interpolated, returns its moment
  function atomInterpolation(interp,atom,ensemble,emom) result(v)
    implicit none
    type(LocalInterpolationInfo), intent(in) :: interp
    integer, intent(in) :: atom !< atom index
    integer, intent(in) :: ensemble !< current ensemble
    real(dblprec), dimension(:,:,:), intent(in) :: emom   !< Current unit moment vector
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

  !> determines wether or not an atom is interpolated
  function isInterpolated(interp,atom) result(v)
    implicit none
    type(LocalInterpolationInfo), intent(in) :: interp
    integer, intent(in) :: atom !< atom index
    logical :: v
    v = .false.
    if (associated(interp%indices)) then
       v = interp%indices(atom) .ne. 0
    endif

  end function isInterpolated


  subroutine multiscaleInterpolateInterfaces()
    use MomentData, only: emom,emom2, emomM, mmom
    implicit none
    integer :: atom,index,ens

    !$omp parallel do private(atom,index,ens)
    do atom=1,ubound(interfaceInterpolation%indices, 1) 
       index = interfaceInterpolation%indices(atom)
       if (index .ne. 0) then
          do ens=1,ubound(emom,3)
             emom(:,atom,ens) = &
                  atomInterpolation(interfaceInterpolation,atom,ens,emom)
             emom2(:,atom,ens) = &
                  atomInterpolation(interfaceInterpolation,atom,ens,emom2)
             emomM(:,atom,ens) = emom(:,atom,ens) * mmom(atom,ens)
          end do
       end if
    end do
    !$omp end parallel do
  end subroutine multiscaleInterpolateInterfaces

  
   
  subroutine multiscaleInterpolateArray(array)
    implicit none
    real(dblprec), dimension(:,:,:), intent(inout) :: array
    integer :: atom,index,ens    
    
    !$omp parallel do private(atom,index,ens)
    do atom=1,ubound(interfaceInterpolation%indices, 1) 
       index = interfaceInterpolation%indices(atom)
       if (index .ne. 0) then
          do ens=1,ubound(array,3)
             array(:,atom,ens) = &
                  !denormalSlerp(interfaceInterpolation,atom,ens,array)
                  denormalInterpolation(interfaceInterpolation,atom,ens,array)
          end do
       end if
    end do
    !$omp end parallel do
  end subroutine multiscaleInterpolateArray
  
  
  subroutine printInterpolationData(interp)
    implicit none
    type(LocalInterpolationInfo), intent(in) :: interp
    integer :: i,j
    do i=1,ubound(interp%indices,1)
       if(interp%indices(i) .ne. 0) then
          print *,"ATOM", i
          do j=interp%firstNeighbour(interp%indices(i)),&
               interp%firstNeighbour(interp%indices(i)+1)-1
             print *,"NEIGHBOUR", interp%neighbours(j), interp%weights(j)
          end do
       endif
    end do
  end subroutine printInterpolationData
 
  
end module MultiscaleInterpolation
