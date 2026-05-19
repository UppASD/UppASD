!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module areaCoefficients
  use ShapeModule, only: BoxShape
  use GeometryModule
  use Parameters
  use DynamicArray
  use SortModule
  use KdTreeModule
  use FiniteDifference
  implicit none

  integer :: count_dtc
  integer, parameter :: CACHE_SIZE = 5000
  real(dblprec), dimension(3,CACHE_SIZE) :: closest_cache_point
  integer, dimension(CACHE_SIZE) :: closest_cache_value


  abstract interface
     subroutine SubboxPeriodicHandler(subbox, wrapped)
       import :: BoxShape
       type(BoxShape), intent(in) :: subbox
       ! True if the box is consequence of periodic boundary condition wrapping
       logical,intent(in) :: wrapped 
     end subroutine SubboxPeriodicHandler
  end interface
  
  abstract interface
     subroutine SubboxHandler(subbox, isAtomistic)
       import :: BoxShape
       type(BoxShape), intent(in) :: subbox
       logical,intent(in) :: isAtomistic
     end subroutine SubboxHandler
  end interface

public iterativeAreaCoefficients, boxesAreaCoefficients
  
private
contains

  !> For each atom around the given box, approximates the area that is
  !! closest to that atom than all the others.
  !! @param mesh Finite difference mesh drawing.
  !! @param box Region where the area is calculated, atoms outside may be considered.
  !! @param latSp Lattice spacing.
  !! @param positions Coordinates of atoms.
  !! @param tree  KdTree used to find atoms in positions.
  !! @param[out] atoms indices of the atoms considered for the calculation.
  !! @param[out] areas area affected by each atom listed in atoms (same order)
  subroutine boxesAreaCoefficients(mesh, box, latSp, positions, tree, atoms, areas)
    implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    type(BoxShape), intent(in) :: box
    real(dblprec),  intent(in) :: latSp
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(KdTree),   intent(in)    :: tree
    type(DynArrayInt),  intent(inout) :: atoms
    type(DynArrayReal), intent(inout) :: areas

    type(BoxShape) :: exp_box, atom_bbox
    integer, pointer,dimension(:) :: atomistic,continuous
    integer :: first_atom, first_cont
    first_cont=-1
    
    exp_box%corner = box%corner - mesh%boxSize*5d-1 - 1d-10
    exp_box%sizes  = box%sizes + mesh%boxSize + 2d-10
   
    call clearArray(atoms)
    call clearArray(areas)
    
    call getNeighbours(mesh%space,positions,tree, &
         exp_box%corner + exp_box%sizes/2.0_dblprec, &
         sqrt(real(mesh%space%spatDimension))/2.0_dblprec * &
         maxval(exp_box%sizes(1:mesh%space%spatDimension)), &
         atoms)

 
    call splitIndicesByDomain(mesh,positions,atoms, atomistic,continuous, &
         first_atom, first_cont)
    
    call ensureAllocLength(areas,atoms%length)
    areas%values = 0
    areas%length = atoms%length

    call splitBoxByDomain(mesh,box,boxesByDomain)

  contains

    !> Gives two pointers, one containing the set of atom indices
    !! in the atomistic domain and the other in the continuous.
    !! For that purpose the indices are reordered in atoms,
    !! so that all atoms come after first_atom and all continuous nodes
    !! are after first_cont.
    subroutine splitIndicesByDomain(mesh,positions,atoms,&
         atomistic,continuous, first_atom, first_cont)
      implicit none
      type(FiniteDiffMesh), intent(in) :: mesh
      real(dblprec), dimension(:, :), intent(in) :: positions
      type(DynArrayInt),  intent(inout) :: atoms
      integer, dimension(:), pointer, intent(inout) :: atomistic, continuous
      integer, intent(inout) :: first_atom, first_cont
      integer :: i,j, swap_tmp

      j = atoms%length
      i = 1
      do while(i<j)
         if(isPointInsideAtomisticBox(mesh,positions(:,atoms%values(i)))) then
            swap_tmp = atoms%values(i)
            atoms%values(i) = atoms%values(j)
            atoms%values(j) = swap_tmp
            j = j - 1
         else
            i = i + 1
         end if
      end do
      if(isPointInsideAtomisticBox(mesh,positions(:,atoms%values(i)))) then
         first_atom = i
         first_cont = 1
         continuous => atoms%values(1:i-1)
         atomistic => atoms%values(i:atoms%length)
      else
         first_atom = i+1
         first_cont = 1
         continuous => atoms%values(1:i)
         atomistic => atoms%values(i+1:atoms%length)
      end if      
    end subroutine splitIndicesByDomain

    !! Callback that receives each subbox
    subroutine boxesByDomain(box,isAtom)
      implicit none      
      type(BoxShape), intent(in) :: box
      logical, intent(in) :: isAtom
      
      integer :: i, off
      integer, pointer,dimension(:) :: atom_set

      if(isAtom) then
         atom_bbox%sizes = latSp
         off = first_atom - 1 
         atom_set => atomistic
      else
         atom_bbox%sizes = mesh%boxSize
         off = first_cont - 1
         atom_set => continuous
      end if
      do i=1,ubound(atom_set,1)
         atom_bbox%corner = positions(:,atom_set(i)) - atom_bbox%sizes*5d-1
         areas%values(i+off) = areas%values(i+off) +&
              intersectVolumeBC(atom_bbox,box,mesh%space)
      end do
    end subroutine boxesByDomain
    
  end subroutine boxesAreaCoefficients


  
  !> Given a box, splits it in subboxes such that none of them is affected by
  !! periodic boundary conditions.
  subroutine splitBoxByPbc(space, box, handler)
    implicit none
    type(SpaceStruct), intent(in) :: space
    type(BoxShape), intent(in) :: box
    procedure(SubboxPeriodicHandler) :: handler
    
    type(BoxShape) :: spaceBB, trimmedBox, offset, candidate
    integer :: i,j,k, dim
    integer, dimension(2,3) :: ranks
    logical, dimension(3) :: periodic

    dim = space%spatDimension
    ! Space bounding box
    spaceBB%corner = 0
    spaceBB%sizes = space%universeSize
    ! Trimming the box so it is never larger than the universe prevents self-intersection,
    ! which would cause two calls to the handler on overlapping boxes.
    trimmedBox = box
    trimmedBox%sizes(1:dim) = min(trimmedBox%sizes(1:dim),space%universeSize(1:dim))

    ! We´ll offset the box on each direction with periodicity, one universe size below and above.
    ! Here we find the offsets we´ll apply, in terms of universe sizes, for each direction.
    ! ranks(1,:) contains the lower-bound offsets and ranks(2,:) the upper bounds.
    periodic = space%periodicBoundary .and. (abs(space%universeSize) >1d-10)    
    ranks(1,:) = merge((/-1,-1,-1/),(/0,0,0/), periodic)
    ranks(2,:) = merge((/1,1,1/),(/0,0,0/), periodic)

    do i=ranks(1,1),ranks(2,1)
       do j=ranks(1,2),ranks(2,2)
          do k=ranks(1,3),ranks(2,3)
             offset = trimmedBox
             offset%corner = offset%corner &
                  +((/i,j,k/)*space%universeSize)
             ! Intersect the offset box with spaceBB, if the intersection is nonempty notify.
             candidate = intersectBoxes(offset,spaceBB,dim)
             if (all(abs(candidate%sizes(1:dim)) > 1d-10)) then
                call handler(candidate, any((/i,j,k/)/=0))
             end if
          end do
       end do
    end do
    
  end subroutine splitBoxByPbc
  
  
  !> Given a mesh and a box, splits the box in subboxes in such a way that each 
  !! subbox is either totally inside the atomistic or totally inside the continuum
  !! domain. It also deals with boundary conditions, guaranteeing that each subbox
  !! does not intersect a boundary and is fully inside the universe.
  !! @parameter mesh 
  !! @parameter box Box to split
  !! @parameter handler Function called for every subbox.
  subroutine splitBoxByDomain(mesh, box, handler)
    implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    type(BoxShape), intent(in) :: box
    procedure(SubboxHandler) :: handler

    call splitBoxByPbc(mesh%space, box, nonperiodic)
  contains
    subroutine nonperiodic(box, wrapped)
      implicit none
      type(BoxShape),intent(in) :: box      
      logical, intent(in) :: wrapped

      integer :: i,j,k, dim 
      integer,dimension(3) :: nbox
      type(BoxShape) :: fdbox, candidate
      real(dblprec), dimension(3) :: corner

      if(1==0) print *,wrapped !Suppress waring for unused parameter
      
      dim = mesh%space%spatDimension
      
      ! Span of fdiff boxes of box
      nbox = 0
      nbox(1:dim) = floor(box%sizes(1:dim) / mesh%boxSize(1:dim)) + 1
      ! Top-left corner of the first fdiff box
      corner = 0 
      corner(1:dim) = floor(box%corner(1:dim) / mesh%boxSize(1:dim)) * mesh%boxSize(1:dim)
      
      fdbox%sizes = mesh%boxSize
      do i=0,nbox(1)
         do j=0,nbox(2)
            do k=0,nbox(3)
               fdbox%corner = corner + mesh%boxSize*(/i,j,k/) ! Shift to the current fdiff box
               candidate = intersectBoxes(box,fdbox,dim)

               if (all(abs(candidate%sizes(1:dim)) > 1d-10)) then
                  call handler(candidate, &
                       isPointInsideAtomisticBox(mesh,fdbox%corner + 0.5*fdbox%sizes))
               end if               
            end do
         end do
      end do
      
    end subroutine nonperiodic
  end subroutine splitBoxByDomain

  !> Given two boxes calculates the box result of their intersection.
  !! If the boxes do not intersect, a box of size 0,0,0 and corner 0,0,0 is
  !! returned instead. Spatial dimension is taken in account.
  !! @param box1 
  !! @param box2 
  !! @param spatDimension 
  function intersectBoxes(box1,box2,spatDimension) result(inter)
  implicit none
    type(BoxShape),intent(in) :: box1,box2
    integer, intent(in) :: spatDimension
    type(BoxShape) :: inter
    real(dblprec),dimension(2,3) :: box1c,box2c
    integer :: i

    box1c(1,:) = box1%corner
    box1c(2,:) = box1%corner+box1%sizes
    box2c(1,:) = box2%corner
    box2c(2,:) = box2%corner+box2%sizes

    inter%corner=0
    inter%sizes=0
    
    do i=1,spatDimension
       if ((box1c(1,i) < box2c(2,i)) .and. (box2c(1,i) < box1c(2,i))) then
          inter%corner(i) = max(box1c(1,i),box2c(1,i))
          inter%sizes(i)  = min(box1c(2,i),box2c(2,i)) - inter%corner(i)
       else
          inter%corner(i) = 0
          inter%sizes(i)  = 0
          return 
       end if
    end do    
    
  end function intersectBoxes


  
  
  !> Finds the volume of the intersection between box1 and box2
  !! Takes care of periodic boundary conditions
  !! @param box1
  !! @param box2
  !! @param space  
  function intersectVolumeBC( box1, box2, space) result(volume)
  implicit none
    type(BoxShape),intent(in) :: box1,box2
    type(SpaceStruct), intent(in) :: space
    real(dblprec) :: volume

    type(BoxShape) :: box_test
    integer, dimension(3) :: offsetWidth, offset
    integer :: i,j,k

    volume = 0
    
    ! Wiggle box1 and add the total volume of intersections with box2
    offsetWidth = merge((/1,1,1/),(/0,0,0/), space%periodicBoundary)
    box_test%sizes = box1%sizes
    do i = -1*offsetWidth(1),offsetWidth(1)
       do j = -1*offsetWidth(2),offsetWidth(2)
          do k = -1*offsetWidth(3),offsetWidth(3)
             offset = (/i,j,k/)
             box_test%corner = box1%corner + space%universeSize * offset
             volume = volume + intersectVolume(box_test,box2,space%spatDimension)
          end do
       end do
    end do   
  end function intersectVolumeBC

  !! Calculates the volume of the intersection bewteen box1 and box2
  !! Does NOT take care of periodic bc
  function intersectVolume(box1, box2, spatDimension) result(volume)
  implicit none
    type(BoxShape),intent(in) :: box1,box2
    integer, intent(in) :: spatDimension
    real(dblprec) :: volume
    type(BoxShape) :: intersection

    intersection = intersectBoxes(box1,box2, spatDimension)
    volume = product(intersection%sizes(1:spatDimension))
  end function intersectVolume  
  
  !> For each atom around the given box, approximates the area that is
  !! closest to that atom than all the others.
  !! This method is slower than boxesAreaCoefficients but more accurate for 2D/3D
  !! @param box Region where the area is calculated
  !! @param expansion Atoms that are outside the box but closer than expansion will be taken into account too
  !! @param geometry
  !! @param limit Depth limit in iterative process. Larguer means slower but more accurate.
  !! @param atoms[out] indices of the atoms considered for the calculation.
  !! @param areas[out] area affected by each atom listed in atoms (same order)
  subroutine iterativeAreaCoefficients(box,expansion,positions,space,tree,limit, atoms,areas)
  implicit none     
    type(BoxShape), intent(in) :: box
    real(dblprec), dimension(3), intent(in) :: expansion
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(SpaceStruct), intent(in) :: space
    type(KdTree), intent(in) :: tree
    integer,intent(in) :: limit
    type(DynArrayInt),intent(inout) :: atoms
    type(DynArrayReal),intent(inout) :: areas
    
    type(KdTree) :: atomTree

    type(BoxShape) :: exp_box, sane_box
    exp_box%corner = box%corner - expansion
    exp_box%sizes  = box%sizes + 2*expansion
   
    call clearArray(atoms)

    call getNeighbours(space,positions,tree, &
         exp_box%corner + exp_box%sizes/2.0_dblprec, &
         sqrt(real(space%spatDimension))/2.0_dblprec * maxval(exp_box%sizes(1:space%spatDimension)), &
         atoms)
    call sort(atoms%values(1:atoms%length))

    sane_box = box
    if (space%spatDimension < 3) then
       sane_box%corner(3) = 0
       sane_box%sizes(3) = 0
       if (space%spatDimension < 2) then
          sane_box%corner(2) = 0
          sane_box%sizes(2) = 0
       end if
    end if

    !! Reset the cache
    count_dtc = 0
    closest_cache_point = -1
    
    call ensureAllocLength(areas,atoms%length)
    areas%values = 0
    areas%length = atoms%length
    if(atoms%length .gt. 0) then
       call buildKdTree(atomTree, positions, atoms%values(1:atoms%length))
       call iterativeAreaCoefficientsAux(sane_box,space, positions, &
            atomTree,limit, atoms, areas%values)
       call deallocTree(atomTree)
    end if

    
  contains
    !> Recursively find the area affected by each atom inside the given box.
    !! atoms is expected to contain the indices of all possible atoms in geometry
    !!  that intersect the box.
    !! atom_accum should be 0-initialized at the beginning,
    !!  and will contain the area affected by each atom in atoms when the subroutine completes
    !! limit is the depth  iterative method stop splitting the box and
    !!  approximates the area affected by each atom instead.
    recursive subroutine iterativeAreaCoefficientsAux(box, space, positions, tree, limit, atoms, atom_accum)
      implicit none
      type(BoxShape), intent(in)             :: box
      type(SpaceStruct), intent(in) :: space
      real(dblprec), dimension(:, :), intent(in) :: positions
      type(KdTree), intent(in) :: tree
      integer, intent(in)                  :: limit
      type(DynArrayInt), intent(in) :: atoms
      real(dblprec), dimension(:), intent(inout) :: atom_accum      

      logical :: unique, last_iter
      real(dblprec) :: area
      integer :: atom_idx, i, searchIndex
      integer, dimension(8) :: closest_atoms
      real(dblprec), dimension(3) :: offset, half_sizes
      integer, dimension(3) :: mask
      type(BoxShape) :: subbox
      integer :: count

      count =  closestAtomsToCorners(box, space, positions, tree, closest_atoms)
      unique = all(closest_atoms(2:count) .eq. closest_atoms(1))
      last_iter = unique .or. (limit .eq. 0)

      if(last_iter) then
         area = boxVolume(box,space%spatDimension) / count
         do atom_idx=1,count
            searchIndex = searchSortedArray(closest_atoms(atom_idx), atoms%values(1:atoms%length))
            atom_accum(searchIndex) = atom_accum(searchIndex) + area
         end do
      else
         half_sizes = box%sizes/2
         do i=0,(2**(space%spatDimension)-1)
            mask=(/ibits(i,0,1), ibits(i,1,1), ibits(i,2,1)/)
            offset = half_sizes*mask
            subbox%corner = box%corner + offset
            subbox%sizes = half_sizes
            call iterativeAreaCoefficientsAux(subbox, space, positions, tree, &
                 limit-1, atoms, atom_accum)
         end do
      end if

    end subroutine iterativeAreaCoefficientsAux

  end subroutine  iterativeAreaCoefficients

  !> given a box and a set of atoms in a geometry, find the atoms that are closer to each of the corners.
  !! The indices of the closest atoms are returned in the last parameter,
  !! that must be a preallocated integer array of length 8.
  function closestAtomsToCorners(box,space,positions,tree,closest_atoms) result(count)
    use GeometryModule
    implicit none
    type(BoxShape), intent(in)     :: box
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(KdTree), intent(in) :: tree
    integer, dimension(8), intent(out) :: closest_atoms

    integer                 :: atom, i, j

    integer :: count
    real(dblprec), dimension(3)::corner
    real(dblprec), dimension(8) :: distanceSquared
    real(dblprec) :: tmpDistance

    count = 2**space%spatDimension
    do i=1,count
       corner = boxCorner(box,i)
       call cachedClosest(space, corner, positions, tree, atom)
       closest_atoms(i) = atom
       distanceSquared(i) = getDistanceSquared(space, positions(:, atom), corner)
    end do
    
    corner = boxCorner(box, 1)
    do j = 2, count
        if (closest_atoms(j) == closest_atoms(1)) cycle
        tmpDistance = getDistanceSquared(space, positions(:, closest_atoms(j)), corner)
        if (abs(tmpDistance - distanceSquared(1)) < 1d-5) then
            distanceSquared(1) = tmpDistance
            closest_atoms(1) = closest_atoms(j)
        end if
    end do
    do i = 1, count - 1
        do j = i + 1, count
            if (closest_atoms(j) == closest_atoms(i)) cycle
            corner = boxCorner(box, j)
            tmpDistance = getDistanceSquared(space, positions(:, closest_atoms(i)), corner)
            if (abs(tmpDistance - distanceSquared(j)) < 1d-5) then
                distanceSquared(j) = tmpDistance
                closest_atoms(j) = closest_atoms(i)
            end if
        end do
    end do
  end function  closestAtomsToCorners

  !! Cached version of distToClosest
  !! Depending on the access pattern, effective calls to distToClosest are
  !! close to halved with this cache.
  subroutine cachedClosest(space, point, positions, tree, closest)
    implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(3), intent(in) :: point
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(KdTree), intent(in) :: tree
    integer, intent(out), optional  :: closest

    integer :: hash, index
    integer(kind=4) :: tmp1,tmp2,tmp3 
    
    tmp1 = floor(huge(tmp1) * (point(1) / space%universeSize(1)))
    tmp2 = floor(huge(tmp2) * (point(2) / space%universeSize(2)))
    tmp3 = floor(huge(tmp3) * (point(3) / space%universeSize(3)))
    hash = &
         xor(tmp1 * 1024, xor(tmp2 * 2048, tmp3 * 256)) + &
         xor(tmp1 / 1024, xor(tmp2 / 2048, tmp3 / 256))
    
    index = modulo(hash,CACHE_SIZE) + 1

    if(all(closest_cache_point(1:3,index) .eq. point)) then
       closest = closest_cache_value(index)
    else
       count_dtc = count_dtc + 1
       call distToClosest(space, positions, tree, point, closest=closest)
       closest_cache_value(index) = closest
       closest_cache_point(1:3,index) = point
    end if

  end subroutine cachedClosest
  
  !> Gives the nth corner of a box.
  !! The corners are sorted in such a way that A is the first, B is the last.
  !! In general the corner I has its nth coordinate as B if the bit n-1 of the number I is 1
  !! @param box The box whose corner is needed
  !! @param n   The index of the corner, the expected value ranges from 1 to 2**dims
  function boxCorner(box, n) result (corner)
    type(BoxShape), intent(in) :: box
    Integer, intent(in) :: n
    real(dblprec), dimension(3) :: corner

    Integer :: i
    Integer :: mask(3)
    i = n-1
    mask = (/ibits(i,0,1), ibits(i,1,1), ibits(i,2,1)/)
    corner = box%corner + (box%sizes * mask)
  end function boxCorner

  !> volume of a box, given dimensions
  function boxVolume(rgn,dims) result(area)
    implicit none
    type(BoxShape), intent(in)     :: rgn
    integer, intent(in) :: dims
    real(dblprec) :: area
    area = product(rgn%sizes(1:dims))
  end function boxVolume
  
end module areaCoefficients
