!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

 module MultiscaleGradients
  use MultiscaleInterpolation

  !! The calculation of the gradients is designed to be
  !!  fundamentally equivalent to an interpolation of local moments.
  !! For this reason here we make use of the interpolation functionality thought for
  !!  the atomistic-boundary interface.
  type(LocalInterpolationInfo), target :: gradientsIp

  private
  public initializeMultiscaleGradients, calculateMGradients
contains


  subroutine initializeMultiscaleGradients(gradMatrix, nAtom)
    use SparseMatrix, only : SpMatrix    
    implicit none
    type(SpMatrix), intent(in) :: gradMatrix
    integer, intent(in) :: nAtom

    integer :: atom0, origin, nLinked, nLinks, nextLinkSetIdx, link
    !Alias for sorter name
    type(LocalInterpolationInfo), pointer :: lii
    lii => gradientsIp

    call newLocalInterpolationInfo(gradientsIp)

    nLinked = countDifferentPositive(gradMatrix%row%values(1:gradMatrix%row%length))
    nLinks  = gradMatrix%row%length

    allocate(lii%indices(nAtom))
    allocate(lii%firstNeighbour(nLinked+1))
    allocate(lii%neighbours(nLinks))
    allocate(lii%weights(nLinks))

    lii%indices = 0
    lii%nrInterpAtoms = nLinked
    
    atom0 = -1
    nextLinkSetIdx = 1
    do link=1,nLinks
       origin = gradMatrix%row%values(link)
       if (origin /= atom0) then
          atom0 = origin
          lii%indices(origin) = nextLinkSetIdx
          lii%firstNeighbour(nextLinkSetIdx) = link 
          nextLinkSetIdx = nextLinkSetIdx + 1
       end if       
    end do
    
    lii%firstNeighbour(nLinked+1) = nLinks
    lii%neighbours = gradMatrix%col%values(1:gradMatrix%col%length)
    lii%weights = gradMatrix%entries%values(1:gradMatrix%entries%length)
    
  contains
    integer pure function countDifferentPositive(array)
      implicit none
      integer, dimension(:), intent(in) :: array
      integer :: i,previous,count
      
      count = 0
      previous = -1
      do i = 1,ubound(array,1)
         if (array(i) /= previous) then
            count = count + 1
            previous = array(i)
         end if
      end do
      countDifferentPositive = count
    end function countDifferentPositive
  end subroutine initializeMultiscaleGradients


  !> Evaluates the interpolation
  !! (Σ(w_i m_i))  where w_i=weights
  function gradientInterpolation(interp,atom,ensemble,arr) result(v)
    implicit none
    type(LocalInterpolationInfo), intent(in) :: interp
    integer, intent(in) :: atom !< atom index
    integer, intent(in) :: ensemble !< current ensemble
    real(dblprec), dimension(:,:,:), intent(in) :: arr   !< Current unit moment vector
    real(dblprec), dimension(3) :: v,t
    integer :: i
    integer :: index
    if (associated(interp%indices)) then
       index = interp%indices(atom)
       if (index .ne. 0) then
          v = (/ 0, 0, 0 /)
          do i=interp%firstNeighbour(index), &               
               interp%firstNeighbour(index+1)-1
             t = arr(:,interp%neighbours(i),ensemble) * interp%weights(i)
             v = v + t 
          end do          
       else
          v = arr(:,atom,ensemble)
       end if
    end if
  end function gradientInterpolation

  
  subroutine calculateMGradients(array,output)
    implicit none
    real(dblprec), dimension(:,:,:), intent(in) :: array
    real(dblprec), dimension(:,:,:), intent(inout) :: output
    integer :: atom,index,ens
    !$omp parallel do private(atom,index,ens)
    do atom=1,ubound(interfaceInterpolation%indices, 1) 
       do ens=1,ubound(array,3)
          output(:,atom,ens) = &
               gradientInterpolation(gradientsIp,atom,ens,array)
       end do
    end do
    !$omp end parallel do
  end subroutine calculateMGradients
  

  
end module MultiscaleGradients
  
