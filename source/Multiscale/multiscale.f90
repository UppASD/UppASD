!> Fortran frontend for the tool
!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro
module Multiscale
  use Anisotropy
  use AtomGenerator
  use AtomLinker
  use Configuration
  use DampingBand
  use FiniteDifference
  use GeometryModule
  use InteractionInput
  use KDTreeModule
  use Parameters
  use ShapeModule
  use SparseMatrix
  
  type MultiscaleRegions
     ! Indices of atoms that belongs to the fully coarse-grained region
     integer, dimension(:), allocatable :: coarseGrained
     ! Indices of atoms that belong to the partially coarse-grained region
     integer, dimension(:), allocatable :: partiallyCoarseGrained
     ! Indices of atoms that are unaffected by coarse-graining
     integer, dimension(:), allocatable :: realAtoms
     ! Indices of atoms in the damping band
     integer, dimension(:), allocatable :: dampingAtoms
     ! Indices outside the damping band
     integer, dimension(:), allocatable :: nonDampingAtoms 
  end type MultiscaleRegions


  type AtomAnisotropies
     integer, pointer, dimension(:,:) :: anisotropyTypes ! size (atoms,axes)
     real(dblprec), pointer, dimension(:,:,:) :: anisotropyKs ! size (2,atoms,axes)
     real(dblprec), pointer, dimension(:,:,:) :: anisotropyE  ! size (3,atoms,axes)
     real(dblprec), pointer, dimension(:,:)   :: anisotropyRatios ! size(atoms,axes)
  end type AtomAnisotropies
  
  type MultiscaleSetup
     real(dblprec), allocatable, dimension(:,:) :: moments

     real(dblprec), allocatable, dimension(:,:) :: positions

     type(SpMatrix) :: atomsExchange
     type(SpMatrix) :: atomsDm

     type(SpMatrix) :: interpolationWeights

     type(SpMatrix) :: dampingBandWeights
     real(dblprec), allocatable, dimension(:) :: dampingBandAttenuation
     
     type(AtomAnisotropies) :: anisotropies
     type(SpMatrix)         :: gradientLinks
  end type MultiscaleSetup

  real(dblprec), dimension(:,:,:,:), allocatable :: multiscaleBackbuffer
  integer :: multiscaleBackbufferHead


  private
  public &
       AtomAnisotropies, &
       MultiscaleOptions, MultiscaleSetup, MultiscaleRegions, &
       runMultiscaleSetup, createAtoms, setupLinks, &
       setupInterpolationWeights, &
       setupDampingBand, setupAnisotropies, createRegions, &
       multiscaleBackbuffer, multiscaleBackbufferHead, &
       deallocateRegions, deallocateSetup
  
contains

  
  !> Prepare a multiscale setup.
  subroutine runMultiscaleSetup(options, setup,atomRegions)
    use SortModule, only: sort
    implicit none
    type(MultiscaleOptions), intent(inout) :: options 
    type(MultiscaleSetup), intent(out) :: setup
    type(MultiscaleRegions), intent(out) :: atomRegions

    type(FiniteDiffMesh) :: mesh
    integer, dimension(:, :, :), allocatable :: finiteDiffIndices

    type(KdTree) :: realAtomsIndices
    type(KdTree) :: paddingAtomsIndices
    type(KdTree) :: nonDampingIndices
    type(AtomSetupInfo) :: atomSetup
    integer :: nrOfFiniteDiff


    
    ! Set the sizes outside the dimensionality to 1, to avoid numerical problems with 0's
    options%space%universeSize((options%space%spatDimension+1):3) = 1d0;

    call createFiniteDiffMesh(options%space, options%finiteDiffBoxes, &
         options%atomisticShapes, mesh)

    call createAtoms(options, mesh, setup%positions, atomSetup, setup%moments,&
         realAtomsIndices, paddingAtomsIndices, &
         finiteDiffIndices, nrOfFiniteDiff)


    call setupAnisotropies(options, setup%positions, setup%anisotropies)

    print *, 'Creating the regions...'
    call createRegions(options%space, setup%positions, realAtomsIndices%indices, paddingAtomsIndices, options%coarseGrainedWidth, &
         options%partCoarseGrainedWidth, options%dampingbandWidth, atomRegions)

    call setupLinks(options,mesh,setup%positions,atomSetup,finiteDiffIndices, &
         atomRegions, setup%atomsExchange, setup%atomsDm)
    call deallocateAtomSetupInfo(atomSetup)

    call setupInterpolationWeights(options,mesh, setup%positions, realAtomsIndices, &
         paddingAtomsIndices%indices, finiteDiffIndices, &
         setup%interpolationWeights)   

    call buildKdTree(nonDampingIndices, setup%positions, atomRegions%nonDampingAtoms)
    call setupDampingBand(mesh,options,setup%positions, realAtomsIndices,&
         paddingAtomsIndices,atomRegions%dampingAtoms, &
         nonDampingIndices, finiteDiffIndices, &
         setup%dampingBandWeights, setup%dampingBandAttenuation  )
    call deallocTree(nonDampingIndices)  

    call setupGradients(options, mesh, finiteDiffIndices, &
         realAtomsIndices, paddingAtomsIndices, setup%positions, setup%gradientLinks)

    call sort(atomRegions%coarseGrained)
    call sort(atomRegions%partiallyCoarseGrained)
    call sort(atomRegions%realAtoms)
    call sort(atomRegions%dampingAtoms)
    call sort(atomRegions%nonDampingAtoms)
    
    call sortMatrixByRowAndColumn(setup%interpolationWeights)
    call sortMatrixByRow(setup%atomsExchange)
    call sortMatrixByRowAndColumn(setup%atomsDm)
    call sortMatrixByRow(setup%dampingBandWeights)
    
    deallocate(finiteDiffIndices)

    call deallocTree(realAtomsIndices)
    call deallocTree(paddingAtomsIndices)
    call deallocateFiniteDiffMesh(mesh)
    call finalizeGeometry()

  end subroutine runMultiscaleSetup
  
  !> deallocate an atomRegions structure
  subroutine deallocateRegions(atomRegions)
    implicit none
    type(MultiscaleRegions), intent(inout) :: atomRegions

    deallocate(atomRegions%partiallyCoarseGrained)
    deallocate(atomRegions%coarseGrained)
    deallocate(atomRegions%realAtoms)
    deallocate(atomRegions%dampingAtoms)
    deallocate(atomRegions%nonDampingAtoms)

  end subroutine deallocateRegions

  !> deallocate multiscaleSetup structure
  subroutine deallocateSetup(setup)
    implicit none
    type(MultiscaleSetup), intent(inout) :: setup

    if(allocated(setup%moments)) &
         deallocate(setup%moments)
    if(allocated(setup%positions)) &
         deallocate(setup%positions)
    call deallocSpMatrix(setup%atomsExchange)
    call deallocSpMatrix(setup%atomsDm)
    call deallocSpMatrix(setup%interpolationWeights)
    call deallocSpMatrix(setup%dampingBandWeights)
    if(allocated(setup%dampingBandAttenuation)) &
         deallocate(setup%dampingBandAttenuation)
    call deallocateAnisotropies(setup%anisotropies)
    call deallocSpMatrix(setup%gradientLinks)
  end subroutine deallocateSetup

  !> Deallocate atom anisotropies
  subroutine deallocateAnisotropies(anisotropies)
    implicit none
    type(AtomAnisotropies), intent(inout) :: anisotropies

    if (associated(anisotropies%anisotropyTypes)) &
         deallocate(anisotropies%anisotropyTypes)
    if(associated(anisotropies%anisotropyKs)) & 
         deallocate(anisotropies%anisotropyKs)
    if(associated(anisotropies%anisotropyE)) & 
         deallocate(anisotropies%anisotropyE)
    if(associated(anisotropies%anisotropyRatios)) & 
         deallocate(anisotropies%anisotropyRatios)
        
  end subroutine deallocateAnisotropies
  
  !> Given a finite difference mesh and a MultiscaleOptions,
  !! builds the corresponfing real and padding atoms.
  !! This is done iterating over regions in the mesh that are atomistic and adding
  !! atoms as specified by the unit cell, wherever the holeShapes are not present.
  !! Padding atoms are placed based on their distance to the atomistic region.
  !! @param[in] opts Options, indicates the unit cell, holes and padding width
  !! @param[in] mesh Finite difference mesh
  !! @param[out] atomPositions Atom coordinates, dimensions (3, Natoms)
  !! @param[out] atomSetup Information about the atom types and reference atoms in the UC
  !! @param[out] atomMoments Moments from the input file
  !! @param[out] realAtomIndices k-d tree indexing all atoms without coarse graining.
  !! @param[out] paddingAtomIndices k-d tree indexing all atoms affected by the damping band.
  !! @param[out] finiteDiffIndices Integer indices of finite-difference nodes by their i,j,k indices. A tree is not needed as nodes are placed in order in a grid.
  !! @param[out] nrOfFiniteDiff Number of finite difference nodes.
  subroutine createAtoms(opts, mesh, atomPositions, atomSetup, atomMoments, &
       realAtomIndices, paddingAtomIndices, finiteDiffIndices, nrOfFiniteDiff)
    use MomentRegions
    implicit none

    type(MultiscaleOptions), intent(in) :: opts
    type(FiniteDiffMesh), intent(in) :: mesh
    real(dblprec), dimension(:, :), allocatable, intent(out) :: atomPositions
    type(AtomSetupInfo), intent(out) :: atomSetup
    real(dblprec), dimension(:, :), allocatable, intent(out) :: atomMoments
    type(KdTree), intent(out) :: realAtomIndices
    type(KdTree), intent(out) :: paddingAtomIndices
    integer, dimension(:, :, :), allocatable, intent(out) :: finiteDiffIndices
    integer, intent(out) :: nrOfFiniteDiff

    type(AtomLinkedList) :: generatedAtomList
    real(dblprec), dimension(:, :), allocatable :: realAtomPositions
    integer, dimension(:), allocatable :: tmpPaddingIndices
    integer :: nrOfAtomsGenerated
    integer :: nrOfPaddingAtoms
    integer :: i
    type(MomentList), pointer :: p
    
    call allocateAtomLinkedList(generatedAtomList)

    print *, 'Generating atoms...'
    call generateAtoms(mesh, &
         opts%unitCell, opts%holeShapes, opts%atomZones, &
         generatedAtomList, nrOfAtomsGenerated)

    allocate(realAtomPositions(3, nrOfAtomsGenerated))
    call extractPositions(generatedAtomList, realAtomPositions)  

    call buildKdTree(realAtomIndices, realAtomPositions)


    print *, 'Generating padding atoms...'
    call generatePaddingAtoms(mesh, opts%unitCell, realAtomPositions, &
         realAtomIndices, opts%paddingWidth, generatedAtomList, nrOfPaddingAtoms)
    deallocate(realAtomPositions)


    tmpPaddingIndices = (/ (i, i = nrOfAtomsGenerated + 1, nrOfAtomsGenerated + nrOfPaddingAtoms) /)

    call allocateAtomSetupInfo(atomSetup, nrOfAtomsGenerated + nrOfPaddingAtoms)
    call extractTypes(generatedAtomList, atomSetup%atomTypes)
    call extractFromUnitcellLocation(generatedAtomList, atomSetup%fromUnitcellLocation)
   
    call generateFiniteDifferenceAtoms(mesh, nrOfAtomsGenerated + nrOfPaddingAtoms, &
         opts%continuumMomentMagnitude, finiteDiffIndices, generatedAtomList, nrOfFiniteDiff)
    allocate(atomPositions(3, nrOfAtomsGenerated + nrOfPaddingAtoms + nrOfFiniteDiff))
    allocate(atomMoments(3, nrOfAtomsGenerated + nrOfPaddingAtoms + nrOfFiniteDiff))

    call extractPositions(generatedAtomList, atomPositions)
    call extractMoments(generatedAtomList, atomMoments)
    do i=1,ubound(atomPositions,2)
       p => momentFromPoint(atomPositions(:,i),opts%momentRegions)
       if (associated(p)) then
          atomMoments(:,i) = p%parameters%magnitude * p%parameters%direction
       end if
    end do
   
    call buildKdTree(paddingAtomIndices, atomPositions, tmpPaddingIndices)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(tmpPaddingIndices)

    print *, 'Number of atoms generated: ', nrOfAtomsGenerated
    print *, 'Number of padding atoms generated: ', nrOfPaddingAtoms

    call deallocateAtomList(generatedAtomList)
  end subroutine createAtoms

  !> Calculates interactions between atoms and creates releated configuration files.
  !! @param[in] opts Options, indicates space features and exchange parameters for different regions.
  !! @param[in] mesh Finite differences mesh
  !! @param[in] atomPositions Atom coordinates, dimensions (3, Natoms)
  !! @param[in] atomSetup Information about the atom types and reference atoms in the UC
  !! @param[in] finiteDiffIndices
  !! @param[in] atomRegions Delimiters for different regions inside realAtoms
  !! @param[out] atomsExchange Sparse matrix defining exchanges between pairs of atoms.
  subroutine setupLinks(opts, mesh, atomPositions, atomSetup, finiteDiffIndices, atomRegions, atomsExchange, atomsDm)
    use AtomLinker
    use DM
    implicit none

    type(MultiscaleOptions) , intent(inout)                  :: opts
    type(FiniteDiffMesh)    , intent(in)                     :: mesh
    real(dblprec)           , intent(in), dimension(:, :)    :: atomPositions
    type(AtomSetupInfo)     , intent(in)                     :: atomSetup
    integer                 , intent(in), dimension(:, :, :) :: finiteDiffIndices
    type(MultiscaleRegions) , intent(in)                     :: atomRegions    
    type(SpMatrix)          , intent(out)                    :: atomsExchange
    type(SpMatrix)          , intent(out)                    :: atomsDm

    integer :: i, dim
    real(dblprec) :: realMaxRadius, coarseMaxRadius

    type(KdTree) :: totalTree, realTree, partTree, coarseTree
    integer, dimension(:), allocatable :: atomIndices
    integer, dimension(3) :: maxNrOfUnitcellLocations
 
    dim = opts%space%spatDimension
    maxNrOfUnitcellLocations = 1
    maxNrOfUnitcellLocations(1:dim) = &
         floor(mesh%space%universeSize(1:dim)/opts%unitcell%size(1:dim) + 0.5d0)
    allocate(atomIndices(atomSetup%nrOfAtoms))
    atomIndices = (/ (i, i = 1, atomSetup%nrOfAtoms) /)
    call buildKdTree(totalTree, atomPositions, atomIndices)
    deallocate(atomIndices)

    call buildKdTree(realTree, atomPositions, atomRegions%realAtoms)
    call buildKdTree(partTree, atomPositions, atomRegions%partiallyCoarseGrained)
    call buildKdTree(coarseTree, atomPositions, atomRegions%coarseGrained)

    call allocSpMatrix(atomsExchange)
    call allocSpMatrix(atomsDm)
    print *, 'Linking the finite difference atoms'

    call createFiniteDiffLinks(mesh, &
         opts%continuumExchangeCoef, opts%continuumDm, &
         finiteDiffIndices, atomsExchange, atomsDm)
    
    print *, 'Creating links between the atoms...'
    realMaxRadius = &
         getMaxInteractionRadius(opts%realExchange, opts%unitcell%size) &
         + opts%linkErrorTolerance
    
    coarseMaxRadius = &
         getMaxInteractionRadius(opts%coarseExchange, opts%unitcell%size) &
         + opts%linkErrorTolerance

    call createExchangeMatrix(opts%space, atomPositions, opts%linkErrorTolerance,&
         totalTree, realTree, partTree, coarseTree, &
         realExchangeLaw, realMaxRadius, coarseExchangeLaw, coarseMaxRadius, &
         atomsExchange)

    realMaxRadius = &
         getMaxInteractionRadius(opts%realDm, opts%unitcell%size) &
         + opts%linkErrorTolerance
    
    coarseMaxRadius = &
         getMaxInteractionRadius(opts%coarseDm, opts%unitcell%size) &
         + opts%linkErrorTolerance
    
    call createDmMatrix(opts%space, atomPositions,opts%linkErrorTolerance,&
         totalTree, realTree, partTree, coarseTree, &
         realDmLaw, realMaxRadius, coarseDmLaw, coarseMaxRadius, &
         atomsDm)

    
    call deallocTree(totalTree)
    call deallocTree(realTree)
    call deallocTree(partTree)
    call deallocTree(coarseTree)

  contains

    pure real(dblprec) function realExchangeLaw(atomI, atomJ)
      use InteractionInput
      integer, intent(in) :: atomI, atomJ

      real(dblprec), dimension(1) :: iv

      iv = getInteractionValue(opts%realExchange, opts%space, &
           maxNrOfUnitcellLocations, atomPositions, atomSetup, atomI, atomJ, &
           opts%linkErrorTolerance)
      realExchangeLaw = iv(1)
      
    end function realExchangeLaw
    pure real(dblprec) function coarseExchangeLaw(atomI, atomJ)
      use InteractionInput
      integer, intent(in) :: atomI, atomJ

      real(dblprec), dimension(1) :: iv

      iv = getInteractionValue(opts%coarseExchange, opts%space, &
           maxNrOfUnitcellLocations, atomPositions, atomSetup, atomI, atomJ, &
           opts%linkErrorTolerance)
      coarseExchangeLaw = iv(1)
      
    end function coarseExchangeLaw
    
    pure function realDmLaw(atomI, atomJ) result(v)
      use InteractionInput
      integer, intent(in) :: atomI, atomJ
      real(dblprec), dimension(3) :: v

      v = getInteractionValue(opts%realDm, opts%space, &
           maxNrOfUnitcellLocations, atomPositions, atomSetup, atomI, atomJ, &
           opts%linkErrorTolerance)
            
    end function realDmLaw
    pure function coarseDmLaw(atomI, atomJ) result(v)
      use InteractionInput
      integer, intent(in) :: atomI, atomJ
      real(dblprec), dimension(3) :: v

      v = getInteractionValue(opts%coarseDm, opts%space, &
           maxNrOfUnitcellLocations, atomPositions, atomSetup, atomI, atomJ, &
           opts%linkErrorTolerance)
      
    end function coarseDmLaw
    
    
  end subroutine setupLinks

  !> Calculates weights used to calculate two interpolations in the interface between atomistic and continuum domains. All results are stored in files.
  !! These interpolations are used to compute the magnetic moment of padding atoms and interface nodes.
  !! @param[in]  opts Options, Provides space configuration and lattice spacing.
  !! @param[in]  mesh Finite differences mesh.
  !! @param[in]  atomPositions Atom coordinates, dimensions (3, Natoms)
  !! @param[in]  realAtomIndices k-d tree indexing real atoms.
  !! @param[in]  paddingIndices Integer indices of atoms in the padding zone.
  !! @param[in]  finiteDiffIndices 
  !! @param[out] interpolationWeights calculated interpolation weights.
  subroutine setupInterpolationWeights(opts, mesh, atomPositions, realAtomIndices,&
       paddingIndices, finiteDiffIndices, interpolationWeights)
    use DynamicArray
    use FiniteDifference
    implicit none
    type(MultiscaleOptions), intent(inout) :: opts
    type(FiniteDiffMesh), intent(in) :: mesh
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(KdTree), intent(in) :: realAtomIndices
    integer, dimension(:), intent(in) :: paddingIndices
    integer, dimension(:, :, :), allocatable, intent(in) :: finiteDiffIndices
    type(SpMatrix), intent(out) :: interpolationWeights

    print *, 'Creating interpolation weights for the finite difference nodes...'
    call allocSpMatrix(interpolationWeights)    
    call createFiniteDiffInterpolationLinks(mesh, opts%space, atomPositions,&
         realAtomIndices, mesh%boxSize, opts%atomLatSp, finiteDiffIndices,&
         interpolationWeights)
    print *, 'Creating interpolation weights for the padding atoms...'
    call createPaddingInterpolationWeights(atomPositions, paddingIndices, mesh, finiteDiffIndices,&
         interpolationWeights)
    
    call removeZeros(interpolationWeights)

  end subroutine setupInterpolationWeights
  
  subroutine setupGradients(opts, mesh, finiteDiffIndices, realAtomIndices, paddingIndices, positions, indices)
    use GradientIndices
    implicit none
    type(MultiscaleOptions), intent(inout) :: opts
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, dimension(:, :, :), intent(in) :: finiteDiffIndices
    type(KdTree), intent(in) :: realAtomIndices
    type(KdTree), intent(in) :: paddingIndices
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(SpMatrix), intent(out) :: indices
    
    !! Allocate output
    call allocSpMatrix(indices)

    if (sum(abs(opts%sttVector)) > 0) then
       print *,'Calculating STT coefficient vectors.'
       call generateDirectionalDerivativeLinks(&
            realAtomIndices, paddingIndices, positions, opts%atomLatSp, mesh, &
            finiteDiffIndices, opts%sttWindowSize, opts%sttVector, &
            indices)
    end if
  end subroutine setupGradients
  

  !> Calculates coefficients used in the damping band filter. These parameters are stored as files.
  !! @param[in]  mesh Finite difference mesh
  !! @param[in]  opts Options, used for space features, lattice spacing and damping band width.
  !! @param[in]  atomPositions Atom coordinates, dimensions (3, Natoms)
  !! @param[in]  realAtomIndices k-d tree indexing real atoms.
  !! @param[in]  paddingIndices k-d of atoms in the padding zone.
  !! @param[in]  dampingIndices Integer indices of atoms affected by the damping band.
  !! @param[in]  nonDampingindices Integer indices of atoms not affected by the damping band.
  !! @param[in]  finiteDiffIndices indices of finite-difference nodes.
  !! @param[out] dampingAvCoeffs Averaging weights for the damping band
  !! @param[out] dampingSpCoeffs Position-based damping band attenuation coefficients.
  subroutine setupDampingBand(mesh, opts, atomsPositions, realAtomIndices, &
       paddingIndices, dampingIndices, nonDampingIndices, finiteDiffIndices, &
       dampingAvCoeffs, dampingSpCoeffs)
    use DynamicArray
    implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    type(MultiscaleOptions), intent(inout) :: opts
    real(dblprec), dimension(:, :), intent(in) :: atomsPositions
    type(KdTree), intent(in) :: realAtomIndices
    type(KdTree), intent(in) :: paddingIndices
    integer, dimension(:), intent(in) :: dampingIndices
    type(KdTree), intent(in) :: nonDampingIndices
    integer, dimension(:, :, :), intent(in) :: finiteDiffIndices
    type(SpMatrix), intent(out) :: dampingAvCoeffs
    real(dblprec), dimension(:), allocatable, intent(out) :: dampingSpCoeffs

    type(DynArrayInt) :: finiteIndices
    integer, dimension(:), allocatable :: allIndices
    type(KdTree) :: totalTree

    if(opts%dampingBandStrength < 1d-16) then
       print *, 'Damping band strength is zero; no damping band.'
       allocate(dampingSpCoeffs(0))
       call allocSpMatrix(dampingAvCoeffs)
    else    
       print *, 'Calculating distance coefficients for damping band...'
       call dampingPositionalCoefficients(dampingIndices, nonDampingIndices, &
            opts%space, atomsPositions, paddingIndices, opts%atomLatSp, &
            opts%dampingBandStrength, opts%dampingbandWidth, dampingSpCoeffs)
       
       print *, 'Calculating averaging coefficients for the damping band...'
       
       call newArray(finiteIndices)
       call getNonInterpolationIndices(finiteDiffIndices, finiteIndices)
       allocate(allIndices(finiteIndices%length + ubound(realAtomIndices%indices, 1)))
       allIndices(1:finiteIndices%length) = finiteIndices%values(1:finiteIndices%length)
       allIndices((finiteIndices%length + 1):(ubound(allIndices, 1))) = realAtomIndices%indices
       
       call buildKdTree(totalTree, atomsPositions, allIndices)
       
       call allocSpMatrix(dampingAvCoeffs)
       call dampingAreaCoeff(mesh,opts%space, atomsPositions, totalTree, &
            dampingIndices, opts%atomLatSp, opts%windowSize, dampingAvCoeffs)
       print *, 'Multiscale setup ready.'
       
       deallocate(allIndices)
       call deallocArray(finiteIndices)
       call deallocTree(totalTree)
    end if
  end subroutine setupDampingBand

  !> Calculates per-atom anisotropy parameters given the atom postions and the
  !! anisotropy regions in the options.
  !! @param[in]  opts Options, Provides space configuration and lattice spacing
  !! @param[in]  atomPositions Atom coordinates, dimensions (3, Natoms)
  !! @param[out] anisotropies Per-atom anisotropy parameters.
  subroutine setupAnisotropies(opts, atomPositions, anisotropies)
    use ShapeModule
    use Anisotropy  
    implicit none
    type(MultiscaleOptions), intent(inout) :: opts 
    real(dblprec), dimension(:, :), allocatable,intent(in) :: atomPositions
    type(AtomAnisotropies), intent(out) :: anisotropies
    
    type(AnisotropyList), pointer :: p 
    integer :: i, axes, nAtoms, axis

    nAtoms = ubound(atomPositions,2)
    axes = countAnisotropyAxes(opts%anisotropies)

    nullify(anisotropies%anisotropyTypes)
    nullify(anisotropies%anisotropyKs)
    nullify(anisotropies%anisotropyE)
    nullify(anisotropies%anisotropyRatios)
    
    if (axes > 0 .and. nAtoms > 0) then
       allocate(anisotropies%anisotropyTypes(nAtoms,axes))
       anisotropies%anisotropyTypes = 0
       allocate(anisotropies%anisotropyKs(2,nAtoms,axes))
       anisotropies%anisotropyKs = 0
       allocate(anisotropies%anisotropyE(3,nAtoms,axes))
       !! UppASD divides over the norm of e, cannot be 0
       anisotropies%anisotropyE(1,:,:) = 1 
       anisotropies%anisotropyE(2,:,:) = 0
       anisotropies%anisotropyE(3,:,:) = 0
       allocate(anisotropies%anisotropyRatios(nAtoms,axes))
       anisotropies%anisotropyRatios = 0
       
       do i = 1, nAtoms
          p => anisotropyFromPoint(atomPositions(:,i),opts%anisotropies)
          if(associated(p)) then
             do axis=1, p%axes
                anisotropies%anisotropyTypes(i,axis) = p%parameters(axis)%type
                anisotropies%anisotropyKs(1,i,axis) = p%parameters(axis)%K1
                anisotropies%anisotropyKs(2,i,axis) = p%parameters(axis)%K2
                anisotropies%anisotropyE(:,i,axis) = p%parameters(axis)%e
                anisotropies%anisotropyRatios(i,axis) = p%parameters(axis)%ratio
             end do
          end if
       end do
    end if
  end subroutine setupAnisotropies
     
  !> Split the atomistic zone into regions
  !! @param[in]  space Space structure
  !! @param[in]  atomPositions Atom coordinates, dimensions (3, Natoms)
  !! @param[in]  realAtomIndices integer indices real atoms.
  !! @param[in]  paddingAtoms k-d tree indexing atoms in the padding zone.
  !! @param[in]  coarseGrainedWidth width in spatial units of the fully coarse-grained atoms
  !! @param[in]  partiallyCoarsegrainedwidth width in spatial units of the partially coarse-grained atoms
  !! @param[in]  dampingBandWidth width of the damping band in spatial units
  !! @param[in,out] atomRegions Regions structure that will be populated with the indices of different regions of the setup.
  subroutine createRegions(space, atomPositions, realAtomIndices, paddingAtoms, coarseGrainedWidth, &
       partiallyCoarseGrainedWidth, dampingbandWidth, atomRegions)
    implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    integer, dimension(:), intent(in) :: realAtomIndices
    type(KdTree), intent(in) :: paddingAtoms
    real(dblprec), intent(in) :: coarseGrainedWidth
    real(dblprec), intent(in) :: partiallyCoarseGrainedWidth
    real(dblprec), intent(in) :: dampingbandWidth
    type(MultiscaleRegions), intent(inout) :: atomRegions

    integer, parameter :: COARSE_GRAINED = 0
    integer, parameter :: PART_COARSE_GRAINED = 1
    integer, parameter :: REAL_ATOM = 2
    integer, parameter :: DAMPING_ATOM = 3

    integer, dimension(ubound(realAtomIndices, 1)) :: atomRegionId
    integer :: nrTotalAtoms, nrCoarseGrained, nrPartCG, nrRealAtom, nrDamping
    real(dblprec) :: distance
    integer :: i
    nrTotalAtoms = ubound(realAtomIndices, 1)
    nrCoarseGrained = 0
    nrPartCG = 0
    nrRealAtom = 0
    nrDamping = 0
    atomRegionId = 0
    do i = 1, nrTotalAtoms
       call distToClosest(space, atomPositions, paddingAtoms, atomPositions(:, realAtomIndices(i)), &
            distance = distance)
       if (distance < coarseGrainedWidth) then
          atomRegionId(i) = ibset(atomRegionId(i), COARSE_GRAINED)
          nrCoarseGrained = nrCoarseGrained + 1
       elseif (distance < coarseGrainedWidth + partiallyCoarseGrainedWidth) then
          atomRegionId(i) = ibset(atomRegionId(i), PART_COARSE_GRAINED)
          nrPartCG = nrPartCG + 1
       else
          atomRegionId(i) =  ibset(atomRegionId(i), REAL_ATOM)
          nrRealAtom = nrRealAtom + 1
       endif
       if (distance < dampingbandWidth) then
          atomRegionId(i) = ibset(atomRegionId(i), DAMPING_ATOM)
          nrDamping = nrDamping + 1
       endif
    enddo

    allocate(atomRegions%coarseGrained(nrCoarseGrained))
    allocate(atomRegions%partiallyCoarseGrained(nrPartCG))
    allocate(atomRegions%realAtoms(nrRealAtom))
    allocate(atomRegions%dampingAtoms(nrDamping))
    allocate(atomRegions%nonDampingAtoms(nrTotalAtoms-nrDamping))

    nrCoarseGrained = 1
    nrPartCG = 1
    nrRealAtom = 1
    nrDamping = 1
    do i = 1, nrTotalAtoms
       if (btest(atomRegionId(i), COARSE_GRAINED)) then
          atomRegions%coarseGrained(nrCoarseGrained) = realAtomIndices(i)
          nrCoarseGrained = nrCoarseGrained + 1
       elseif (btest(atomRegionId(i), PART_COARSE_GRAINED)) then
          atomRegions%partiallyCoarseGrained(nrPartCG) = realAtomIndices(i)
          nrPartCG = nrPartCG + 1
       elseif (btest(atomRegionId(i), REAL_ATOM)) then
          atomRegions%realAtoms(nrRealAtom) = realAtomIndices(i)
          nrRealAtom = nrRealAtom + 1
       endif
       if (btest(atomRegionId(i), DAMPING_ATOM)) then
          atomRegions%dampingAtoms(nrDamping) = realAtomIndices(i)
          nrDamping = nrDamping + 1
       else
          atomRegions%nonDampingAtoms(i - nrDamping + 1) = realAtomIndices(i)
       endif
    enddo
  end subroutine createRegions
end module Multiscale
