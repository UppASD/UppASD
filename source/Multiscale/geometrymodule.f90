!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module GeometryModule
    use Parameters
    use DynamicArray
    use KdTreeModule
implicit none

    type SpaceStruct
        integer :: spatDimension
        real(dblprec), dimension(3) :: universeSize
        logical, dimension(3) :: periodicBoundary
    end type SpaceStruct

    logical*1, dimension(:), allocatable :: kdtree_already_added

public distToClosest, getNeighbours, getDistance, getDistanceSquared, &
       getDirectionalVector, SpaceStruct, finalizeGeometry

private
contains

subroutine finalizeGeometry()
  implicit none
  if(allocated(kdtree_already_added)) then
     deallocate(kdtree_already_added)
  end if
end subroutine

!> Returns the indices of all the points that are within a given radius of a point.
!! @param[in] space
!! @param[in] positions The positions of the points that should be considered
!! @param[in] tree A tree that have the indices of the points in positions that are of interest
!! @param[in] point The point that the other points should be close to
!! @param[in] radius The radius
!! @param[in,out] neighIndex The indices of the points that are close.
subroutine getNeighbours(space, positions, tree, point, radius, neighIndex)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(KdTree), intent(in) :: tree
    real(dblprec), dimension(3), intent(in) :: point
    real(dblprec), intent(in) :: radius
    type(DynArrayInt), intent(inout) :: neighIndex
    
    real(dblprec), dimension(3) :: pointTmp
    integer, dimension(3) :: offset_width
    integer :: i, j, k
    integer, dimension(3), parameter :: factors = (/ 0, 1, -1 /)
    
    if (.not. allocated(kdtree_already_added)) then
        allocate(kdtree_already_added(ubound(positions, 2)))
        kdtree_already_added = .false.
    elseif (ubound(kdtree_already_added, 1) < ubound(positions, 2)) then
        deallocate(kdtree_already_added)
        allocate(kdtree_already_added(ubound(positions, 2)))
        kdtree_already_added = .false.
    end if
    
    offset_width = 1
    do i = 1, space%spatDimension
        if (space%periodicBoundary(i)) offset_width(i) = 3
    enddo
    
    call clearArray(neighIndex)
    
    do i = 1, offset_width(1)
        pointTmp(1) = point(1) + factors(i)*space%universeSize(1)
        do j = 1, offset_width(2)
            pointTmp(2) = point(2) + factors(j)*space%universeSize(2)
            do k = 1, offset_width(3)
                pointTmp(3) = point(3) + factors(k)*space%universeSize(3)
                
                call getNeighboursKdTree(space, positions, tree, tree%root, pointTmp, radius, neighIndex)
            enddo
        enddo 
    enddo
    do i = 1, neighIndex%length
        kdtree_already_added(neighIndex%values(i)) = .false.
    end do
end subroutine getNeighbours

recursive subroutine getNeighboursKdTree(space, positions, tree, node, point, radius, neighIndex)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(KdTree), intent(in) :: tree
    type(KdTreeNode), intent(in) :: node
    real(dblprec), dimension(3), intent(in) :: point
    real(dblprec), intent(in) :: radius
    type(DynArrayInt), intent(inout) :: neighIndex
    
    integer :: i
    real(dblprec) :: distSquared
    real(dblprec) :: radiusSquared
    
    if (isLeaf(node)) then
        radiusSquared = radius**2
        do i = node%startPoint, node%endPoint
            distSquared = getDistanceSquared(space, point, positions(:, tree%indices(i))) 
            if (distSquared <= radiusSquared .and. .not. kdtree_already_added(tree%indices(i))) then
                call addEntry(neighIndex, tree%indices(i))
                kdtree_already_added(tree%indices(i)) = .true.
            endif
        enddo
    else
        if (node%split >= point(node%splitDim) - radius - 1d-6) then
            call getNeighboursKdTree(space, positions, tree, node%leftNode, point, radius, neighIndex)
        end if
        if (node%split <= point(node%splitDim) + radius + 1d-6) then
            call getNeighboursKdTree(space, positions, tree, node%rightNode, point, radius, neighIndex)
        end if
    endif
end subroutine

!> The distance between point and the closest point in
subroutine distToClosest(space, positions, tree, point, distance, closest)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(KdTree), intent(in) :: tree
    real(dblprec), dimension(3), intent(in) :: point
    real(dblprec), intent(out), optional :: distance
    integer, intent(out), optional  :: closest

    real(dblprec) :: minDist
    integer :: closestAtom
    real(dblprec), dimension(3) :: pointTmp
    integer, dimension(3) :: offset_width
    integer :: i, j, k
    integer, dimension(3), parameter :: factors = (/ 0, 1, -1 /)
    
    offset_width = 1
    do i = 1, space%spatDimension
        if (space%periodicBoundary(i)) offset_width(i) = 3
    enddo
    
    minDist = huge(minDist)
    closestAtom = -1
    do i = 1, offset_width(1)
        pointTmp(1) = point(1) + factors(i)*space%universeSize(1)
        do j = 1, offset_width(2)
            pointTmp(2) = point(2) + factors(j)*space%universeSize(2)
            do k = 1, offset_width(3)
                pointTmp(3) = point(3) + factors(k)*space%universeSize(3)
                
                call distSquaredToClosestKdTree(positions, tree, tree%root, pointTmp, minDist, closestAtom)
            enddo
        enddo
    enddo

    if (present(distance)) then
        distance = sqrt(minDist)
    endif
    if (present(closest)) then
        closest = closestAtom
    endif
end subroutine distToClosest

!! This function does not care about periodic boundary conditions
recursive subroutine distSquaredToClosestKdTree(positions, tree, node, point, distance, closest)
implicit none
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(KdTree), intent(in) :: tree
    type(KdTreeNode), intent(in) :: node
    real(dblprec), dimension(3), intent(in) :: point
    real(dblprec), intent(inout) :: distance
    integer, intent(inout) :: closest
    
    integer :: i, j
    real(dblprec) :: tmpDistance

    if (isLeaf(node)) then
        do j = node%startPoint, node%endPoint
            i = tree%indices(j)
            tmpDistance = sum((point - positions(:, i))**2)
            if (tmpDistance < distance) then
                distance = tmpDistance
                closest = i
            endif
        end do
    else
        if (point(node%splitDim) < node%split) then
            call distSquaredToClosestKdTree(positions, tree, node%leftNode, point, distance, closest)
            if (distance > (node%split - point(node%splitDim))**2) then
                call distSquaredToClosestKdTree(positions, tree, node%rightNode, point, distance, closest)
            endif
        else
            call distSquaredToClosestKdTree(positions, tree, node%rightNode, point, distance, closest)
            if (distance > (point(node%splitDim) - node%split)**2) then
                call distSquaredToClosestKdTree(positions, tree, node%leftNode, point, distance, closest)
            endif
        endif
    end if
end subroutine distSquaredToClosestKdTree

!> This function returns the distance between point1 and point2 in the
!! geometry specified by geom. The function must satisfy that the
!! distance is equal to the 2-norm of the directional vector returned by
!! get DirectionalVector
pure real(dblprec) function getDistance(space, point1, point2) result(distance)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(3), intent(in) :: point1
    real(dblprec), dimension(3), intent(in) :: point2

    distance = sqrt(getDistanceSquared(space, point1, point2))
end function getDistance

pure real(dblprec) function getDistanceSquared(space, point1, point2) result(distance)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(3), intent(in) :: point1
    real(dblprec), dimension(3), intent(in) :: point2

    real(dblprec), dimension(3) :: dirVect
    call getDirectionalVector(space, point1, point2, dirVect)
    distance = sum(dirVect(1:space%spatDimension)**2)
end function getDistanceSquared

!> Returns shortest vector from point1 to point2.
!! Considers periodicity 
!! The vector may point in the opposite direction to where the point is.
pure subroutine getDirectionalVector(space, point1, point2, dirVect)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(3), intent(in) :: point1
    real(dblprec), dimension(3), intent(in) :: point2
    real(dblprec), dimension(3), intent(out) :: dirVect

    integer :: i
    real(dblprec) :: d2

    dirVect = 0.0_dblprec
    do i=1,space%spatDimension
       if (space%periodicBoundary(i)) then
          dirVect(i) = (point2(i) - point1(i))
          d2 = (point2(i) - point1(i) + space%universeSize(i))
          if(abs(d2) < abs(dirVect(i))) then
             dirVect(i) = d2
          else
             d2 = (point2(i) - point1(i) - space%universeSize(i))
             if(abs(d2) < abs(dirVect(i))) then
                dirVect(i) = d2
             endif
          endif
       else
          dirVect(i) = (point2(i) - point1(i))
       endif
    enddo
end subroutine getDirectionalVector

end module GeometryModule
