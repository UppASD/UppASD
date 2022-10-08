!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module DynamicMultiscale 
  use parameters
  use ShapeModule, only : BoxShape
  implicit none

  integer, parameter :: DR_CONT     = 0
  integer, parameter :: DR_FOCUS    = 1
  integer, parameter :: DR_ENVELOPE = 2
  
  
  type DynamicState
     integer, dimension(:,:,:), allocatable :: dynamicRegions 
     
     !WIP real(dblprec) :: atomisticMaxVolume  ! Max volume covered by atomistic regions
     real(dblprec) :: spawnGradientThreshold   ! Minimal gradient needed to spawn a new atomistic region
     real(dblprec) :: despawnGradientThreshold ! Minimal gradient needed to keep an atomistic region (must be less [Half?] than spawnGradientThreshold)
     real(dblprec),dimension(3) :: focusSize   ! Size of the box that covers a focused region
     real(dblprec),dimension(3) :: envelopeSize! Size of the region surrounding the focus
     
     real(dblprec), dimension(:), allocatable  :: gradients !< Gradients     
     integer :: lastUpdated ! Timestep the dynamic state was checked for last time     
  end type DynamicState

  
  private
  public DynamicState, &
       initializeDynamicState, updateDynamicState

contains

  subroutine initializeDynamicState( &
       resolution, &
       spawnGradientThreshold, despawnGradientThreshold, &
       focusSize, envelopeSize, &
       state)
    implicit none
    integer, dimension(3) :: resolution
    real(dblprec), intent(in) :: spawnGradientThreshold
    real(dblprec), optional,  intent(in) :: despawnGradientThreshold
    real(dblprec),dimension(3) :: focusSize
    real(dblprec),dimension(3) :: envelopeSize   
    type(DynamicState), intent(inout) :: state


    allocate(state%dynamicRegions(resolution(1),resolution(2),resolution(2)))
    state%dynamicRegions = DR_CONT 
    state%spawnGradientThreshold = spawnGradientThreshold
    state%despawnGradientThreshold = spawnGradientThreshold * 5d-1
    if (present(despawnGradientThreshold)) state%despawnGradientThreshold = despawnGradientThreshold
    state%focusSize = focusSize
    state%envelopeSize = envelopeSize
    state%lastUpdated = -1
  end subroutine initializeDynamicState

  subroutine updateDynamicState(coords, moments, natom, nlist, nlistsize, mensemble, state)
    implicit none
    real(dblprec), dimension(:,:), intent(inout)     :: coords          !< coordinates
    real(dblprec), dimension(:,:,:), intent(inout)   :: moments         !< Normalized moments
    integer, intent(in)                              :: natom,mensemble !< number of atoms and ensembles
    integer, dimension(:,:), allocatable, intent(in) :: nlist           !< Neighbour list 
    integer, dimension(:), allocatable, intent(in)   :: nlistsize       !< Size of neighbour list
    type(DynamicState), intent(inout)                :: state

    integer, parameter :: CONT = 0
    integer, parameter :: NEEDS_ATOM = 1
    integer, parameter :: ACCEPTS_ATOM = 2

    integer, dimension(3) :: res
    
    integer, dimension(&
         ubound(state%dynamicRegions,1),&
         ubound(state%dynamicRegions,2),&
         ubound(state%dynamicRegions,3)) :: currentState

    integer :: atom, i,j,k

    currentState = state%dynamicRegions
    
    
    call calculateGradients(coords, moments, natom, mensemble, &
         nlist, nlistsize, state%gradients)
    do atom = 1, natom       
       !i = floor(coords(1,atom) * ubound(state%dynamicRegions,1) / universeSize(1))
       !j = floor(coords(2,atom) * ubound(state%dynamicRegions,2) / universeSize(2))
       !k = floor(coords(3,atom) * ubound(state%dynamicRegions,3) / universeSize(3))
!
!       if (currentState(i,j,k) == CONT .and. &
!            state%gradients(atom) > state%spawnGradientThreshold) then
!          currentState(i,j,k) = NEEDS_ATOM
!          
!          
!       end if
       
       !if (state%gradients > state%spawnGradientThreshold) then
          
       !end if
    end do
    
  end subroutine updateDynamicState
    
  
  subroutine findInterpolationCoefficients(oldCoordTree, newCoord, &
       lookupRadius, elements)    
    use KdTreeModule
    implicit none
    type(KdTree), intent(in) :: oldCoordTree
    real(dblprec), dimension(:,:), intent(in) :: newCoord
    real(dblprec), intent(in) :: lookupRadius
    real(dblprec), dimension(:), allocatable :: elements

    
    
  end subroutine findInterpolationCoefficients


  
  
  subroutine calculateGradients (&
       coord, emom, natom, mensemble, &
       nlist, nlistsize,&
       gradients)
    implicit none
    real(dblprec), dimension(:,:), intent(in)   :: coord     !< Coordinates
    real(dblprec), dimension(:,:,:), intent(in) :: emom      !< normalized moments
    integer, intent(in)                         :: natom     !< atom count
    integer, intent(in)                         :: mensemble !< atom count
    integer, dimension(:,:), intent(in)         :: nlist     !< Neighbour list 
    integer, dimension(:),   intent(in)         :: nlistsize !< Size of neighbour list
    real(dblprec), dimension(:), allocatable, intent(inout) :: gradients !< Gradients, Reallocated if needed

    real(dblprec),dimension(3) :: dif
    real(dblprec) :: rij2
    integer :: i,j,k

    if (.not. allocated(gradients)) then
       allocate(gradients(natom))
    elseif (ubound(gradients,1) < natom) then
       deallocate(gradients)
       allocate(gradients(natom))
    end if
    
    do i=1, natom
       dif = 0
       do k=1, mensemble
          do j=1,nlistsize(i)
             rij2 = sum( (coord(:,nlist(j,i))-coord(:,i)) ** 2)
             dif = dif + (emom(1:3, nlist(j,i),k) - emom(1:3, i,k))**2/rij2
          end do
       end do
       gradients(i) = norm2(dif/mensemble)
    end do
    
  end subroutine calculateGradients

end module DynamicMultiscale
