!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module KdTreeModule
    use Parameters
implicit none
    integer, parameter :: MAX_INDICES_IN_LEAF = 10
    
    type KdTreeNode
        real(dblprec) :: split
        type(KdTreeNode), pointer :: leftNode
        type(KdTreeNode), pointer :: rightNode
        integer :: startPoint
        integer :: endPoint
        integer :: splitDim
    end type KdTreeNode
    
    type KdTree
        integer, dimension(:), allocatable :: indices
        type(KdTreeNode) :: root
    end type KdTree

public buildKdTree, deallocTree, isLeaf, KdTree, &
       KdTreeNode
    
private
contains

pure logical function isLeaf(node)
implicit none
    type(KdTreeNode), intent(in) :: node

    isLeaf = node%startPoint > 0
end function

subroutine buildKdTree(tree, points, limitIndices)
implicit none
    type(KdTree), intent(inout) :: tree
    real(dblprec), dimension(:, :), intent(in) :: points
    integer, dimension(:), intent(in), optional :: limitIndices
    
    integer :: i
    
    if (present(limitIndices)) then
        allocate(tree%indices(ubound(limitIndices, 1)))
        tree%indices = limitIndices
    else
        tree%indices = (/ (i, i = 1, ubound(points, 2)) /)
    endif
    
    call newNode(tree%root)
    call buildKdTreeRecursive(tree%root, points, tree%indices, 1, ubound(tree%indices, 1))
end subroutine buildKdTree

pure subroutine newNode(node)
implicit none
    type(KdTreeNode), intent(inout) :: node
    
    nullify(node%leftNode)
    nullify(node%rightNode)
    node%split = 0
    node%startPoint = 0
    node%endPoint = 0
end subroutine newNode

subroutine deallocTree(tree)
implicit none
    type(KdTree), intent(inout) :: tree
    
    deallocate(tree%indices)
    call deallocKdTreeNode(tree%root)
end subroutine

recursive subroutine deallocKdTreeNode(node)
implicit none
    type(KdTreeNode), intent(inout) :: node
    
    if (associated(node%leftNode)) then
        call deallocKdTreeNode(node%leftNode)
        deallocate(node%leftNode)
        nullify(node%leftNode)
    endif
    if (associated(node%rightNode)) then
        call deallocKdTreeNode(node%rightNode)
        deallocate(node%rightNode)
        nullify(node%rightNode)
    endif
    node%startPoint = 0
    node%endPoint = 0
end subroutine deallocKdTreeNode

recursive subroutine buildKdTreeRecursive(node, positions, indices, startIndex, endIndex)
implicit none
    type(KdTreeNode), intent(inout) :: node
    real(dblprec), dimension(:, :), intent(in) :: positions
    integer, dimension(:), intent(inout) :: indices
    integer, intent(in) :: startIndex
    integer, intent(in) :: endIndex
    
    integer :: rightStart
    
    if (startIndex + MAX_INDICES_IN_LEAF > endIndex) then
        node%startPoint = startIndex
        node%endPoint = endIndex
        return
    else
        call computeSplit(node%split, node%splitDim, positions, indices, startIndex, endIndex)
        call partitionPoints(node%splitDim, node%split, positions, indices, startIndex, endIndex, rightStart)
        allocate(node%leftNode)
        allocate(node%rightNode)
        call newNode(node%leftNode)
        call newNode(node%rightNode)
        call buildKdTreeRecursive(node%leftNode, positions, indices, startIndex, rightStart - 1)
        call buildKdTreeRecursive(node%rightNode, positions, indices, rightStart, endIndex)
    end if
end subroutine buildKdTreeRecursive

subroutine computeSplit(split, splitDim, positions, indices, startIndex, endIndex)
implicit none
    real(dblprec), intent(out) :: split
    integer, intent(out) :: splitDim
    real(dblprec), dimension(:, :), intent(in) :: positions
    integer, dimension(:), intent(in) :: indices
    integer, intent(in) :: startIndex
    integer, intent(in) :: endIndex
    
    real(dblprec), dimension(3) :: posMean
    real(dblprec), dimension(3) :: sqrPosSum
    real(dblprec), dimension(3) :: variance
    integer :: i
    
    posMean = 0
    sqrPosSum = 0
    do i = startIndex, endIndex
        posMean = posMean + positions(:, indices(i))
        sqrPosSum = sqrPosSum + positions(:, indices(i))**2
    enddo
    posMean = posMean/(endIndex - startIndex + 1)
    variance = sqrPosSum/(endIndex - startIndex + 1) - posMean**2
    splitDim = maxloc(variance, 1)
    split = posMean(splitDim)
end subroutine computeSplit

subroutine partitionPoints(splitDim, split, positions, indices, startIndex, endIndex, rightStart)
implicit none
    integer, intent(in) :: splitDim
    real(dblprec), intent(in) :: split
    real(dblprec), dimension(:, :), intent(in) :: positions
    integer, dimension(:), intent(inout) :: indices
    integer, intent(in) :: startIndex
    integer, intent(in) :: endIndex
    integer, intent(out) :: rightStart
    
    integer :: i
    integer :: tmpIndex
    
    rightStart = startIndex
    do i = startIndex, endIndex
        if (positions(splitDim, indices(i)) < split) then
            tmpIndex = indices(i)
            indices(i) = indices(rightStart)
            indices(rightStart) = tmpIndex
            rightStart = rightStart + 1
        endif
    enddo

end subroutine partitionPoints

end module KdTreeModule
