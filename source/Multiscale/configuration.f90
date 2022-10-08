!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module Configuration
  use Anisotropy
  use AtomGenerator
  use FileInput
  use FiniteDifference
  use GeometryModule
  use InteractionInput
  use MomentRegions
  use Parameters
  use ShapeModule
  use SparseMatrix


  type MultiscaleOptions
     type(SpaceStruct)                     :: space
     real(dblprec)                         :: dampingbandWidth
     real(dblprec)                         :: coarseGrainedWidth
     real(dblprec)                         :: partCoarseGrainedWidth
     real(dblprec)                         :: paddingWidth
     integer              , dimension(3)   :: finiteDiffBoxes
     type(ShapeList)                       :: atomisticShapes
     type(ShapeList)                       :: holeShapes
     type(AtomCell)                        :: unitCell
     type(InteractionInfo)                 :: realExchange
     type(InteractionInfo)                 :: coarseExchange
     type(InteractionInfo)                 :: realDm
     type(InteractionInfo)                 :: coarseDm     
     real(dblprec)        , dimension(3,3) :: continuumDm
     real(dblprec)                         :: atomLatSp
     real(dblprec)        , dimension(3)   :: continuumExchangeCoef
     real(dblprec)        , dimension(3)   :: windowSize
     real(dblprec)                         :: continuumMomentMagnitude
     type(AtomZoneList)   , pointer        :: atomZones
     type(AnisotropyList) , pointer        :: anisotropies
     type(MomentList)     , pointer        :: momentRegions
     real(dblprec)                         :: linkErrorTolerance
     real(dblprec)                         :: dampingBandStrength
     real(dblprec), dimension(3)           :: sttVector
     real(dblprec), dimension(3)           :: sttWindowSize
  end type MultiscaleOptions


  private
  public MultiscaleOptions, readMultiscaleOptions, deallocateOptions

  type(InputFields) :: fields
contains

  subroutine readMultiscaleOptions(configFile, options)
    implicit none
    character(len=*), intent(in) :: configFile
    type(MultiscaleOptions),intent(out) :: options

    integer :: ierr
    character(len=512) :: errMsg

    call allocateOptions(options)
    call defaultOptions(options)

    call allocateInputFields(fields)
    call setup_inputfields(fields,options, configFile,ierr,errMsg)

    call deallocateInputFields(fields)

    if (ierr /= 0) then
       print '(a)', '------------- ERROR ------------'
       print '(a)', trim(adjustl(errMsg))
       stop 2 ! ENOENT
    else
       call validate_inputedInformation(options)
       call print_inputedInformation(options)
    end if

  end subroutine readMultiscaleOptions

  subroutine defaultOptions(opts)
    implicit none
    type(MultiscaleOptions), intent(inout) :: opts

    opts%dampingbandWidth       = 0.0_dblprec;
    opts%coarseGrainedWidth     = 0.0_dblprec;
    opts%partCoarseGrainedWidth = 0.0_dblprec;
    opts%paddingWidth           = 0.0_dblprec;
    opts%unitCell%size          = 0.0_dblprec
    opts%unitcell%nrOfAtoms     = 0
    opts%finiteDiffBoxes        = 0;
    opts%atomLatSp              = 1.0_dblprec
    opts%continuumExchangeCoef  = 0.0_dblprec
    opts%continuumDm            = 0.0_dblprec
    opts%linkErrorTolerance     = 1d-5
    opts%dampingBandStrength    = 0.0_dblprec
    opts%sttVector              = (/0.0_dblprec,0.0_dblprec,0.0_dblprec/)
    opts%sttWindowSize          = (/0.0_dblprec,0.0_dblprec,0.0_dblprec/)
    nullify(opts%atomZones)
    nullify(opts%anisotropies)
    nullify(opts%momentRegions)

    call initialiseInteractionInfo(opts%realExchange)
    call initialiseInteractionInfo(opts%coarseExchange)
    call initialiseInteractionInfo(opts%realDm)
    call initialiseInteractionInfo(opts%coarseDm)

  end subroutine defaultOptions


  !> Allocate the elements of a MultiscaleOptions structure
  subroutine allocateOptions(options)
    implicit none
    type(MultiscaleOptions), intent(inout) :: options
    nullify (options%atomZones)
    nullify (options%anisotropies)
    nullify (options%momentRegions)
    call allocateShapeList(options%atomisticShapes)
    call allocateShapeList(options%holeShapes)
  end subroutine allocateOptions

  !> Deallocate all elements of a MultiscaleOptions structure.
  !! @param[inout] opts Options struct to deallocate
  subroutine deallocateOptions(opts)
    implicit none
    type(MultiscaleOptions), intent(inout) :: opts

    call deallocateInteractionInfo(opts%realExchange)
    call deallocateInteractionInfo(opts%coarseExchange)
    call deallocateInteractionInfo(opts%realDm)
    call deallocateInteractionInfo(opts%coarseDm)
    call deallocateShapeList(opts%atomisticShapes)
    call deallocateShapeList(opts%holeShapes)
    call deallocateAtomCell(opts%unitCell)
    call deallocateZoneList(opts%atomZones)
    call deallocateAnisotropyList(opts%anisotropies)
    call deallocateMomentList(opts%momentRegions)
  end subroutine deallocateOptions



  subroutine print_inputedInformation(opts)
    use ShapeModule
    implicit none
    type(MultiscaleOptions), intent(in) :: opts
    print *, 'dimension: ', opts%space%spatDimension
    print *, 'universe size: ', opts%space%universeSize
    print *, 'periodic boundaries: ', opts%space%periodicBoundary
    print *, 'damping band width: ', opts%dampingbandWidth
    print *, 'coarse-grained width: ', opts%coarseGrainedWidth
    print *, 'part coarse-grained width: ', opts%partCoarseGrainedWidth
    print *, 'padding width: ', opts%paddingWidth
    print *, 'finite-diff boxes: ', opts%finiteDiffBoxes
    print *, 'max interaction dist: ', opts%atomLatSp
    print *, '----------------------'
    
    print *, 'Continuum exchange: ', opts%continuumExchangeCoef
    print *, 'Real exchange: '
    call printInteractionInfo(opts%realExchange)
    print *, '---------------------'
    print *, 'Coarse exchange: '
    call printInteractionInfo(opts%coarseExchange)
    print *, '----------------------'
    print *, ''

    print *, 'Continuum DM: '
    print *, '  x: ', opts%continuumDM(:,1)
    print *, '  y: ', opts%continuumDM(:,2)
    print *, '  z: ', opts%continuumDM(:,3)
    print *, 'Real DM: '
    call printInteractionInfo(opts%realDm)
    print *, '---------------------'
    print *, 'Coarse DM: '
    call printInteractionInfo(opts%coarseDm)
    print *, '----------------------'
    print *, ''

    print *, '---------------------'
    print *, 'Atomistic shapes: '
    call printShapeList(opts%atomisticShapes)
    print *, '---------------------'
    print *, 'Hole shapes:'
    call printShapeList(opts%holeShapes)
    print *, '---------------------'
    print *, 'Unit cell: '
    print *, 'Size: ', opts%unitCell%size
    print *, 'Number of atoms: ', ubound(opts%unitCell%chemicalTypes, 1)
    print *, '---------------------'
  end subroutine print_inputedInformation


  subroutine validate_inputedInformation(options)
    implicit none
    type(MultiscaleOptions), intent(inout) :: options

    logical :: failed

    failed = .false.
    if (options%space%spatDimension < 1 .or. options%space%spatDimension > 3) then
       print '(a,i4,a)', 'The dimension cant be', options%space%spatDimension, ' it must be either 1, 2 or 3.'
       failed = .true.
    end if
    if (any(options%space%universeSize < 0.0)) then
       print '(a)', 'The size of the universe must be non negative.'
       failed = .true.
    end if

    if (any(options%unitCell%size < 0.0)) then
       print '(a)', 'The unitcell size must be non negative.'
       failed = .true.
    end if

    if (any(options%finiteDiffBoxes <= 0)) then
       print '(a)', 'The number of finite difference boxes must be positive.'
       failed = .true.
    end if

    if (options%atomLatSp <= 0) then
       print '(a)', 'The lattice spacing must be positive.'
       failed = .true.
    end if

    if (getMaxInteractionRadius(options%realExchange, options%unitcell%size) >= options%partCoarseGrainedWidth) then
       print '(a)', 'WARNING: The exchange radius is larger than the partially coarse grained width.'
    else if (getMaxInteractionRadius(options%coarseExchange, options%unitcell%size) &
         > getMaxInteractionRadius(options%realExchange, options%unitcell%size)) then
       print '(a)', 'WARNING: The coarse exchange radius is larger than the real exchange radius.'
    end if

    if (getMaxInteractionRadius(options%coarseExchange, options%unitcell%size) > options%paddingWidth) then
       print '(a)', 'WARNING: The padding width might be too small, it is smaller than the exchange radius in the coarse grained.'
    end if

    if (maxval(abs(options%continuumDM)) == 0    &
         .and. allocated(options%coarseDM%interactionValue) &
         .and. size(options%coarseDM%interactionValue) > 0) then

       print '(a)', 'WARNING: There is a coarse-grained DM but no continuum DM. Did you forget to specify one?'
    end if
    if (maxval(abs(options%continuumDM)) /= 0    &
         .and. (.not. allocated(options%coarseDM%interactionValue) &
         .or. size(options%coarseDM%interactionValue) == 0)) then

       print '(a)', 'WARNING: There is a continuum DM but no coarse-grained DM. Did you forget to specify one?'
    end if

    if (maxval(abs(options%sttVector)) > 0.0 .and. &
         maxval(abs(options%sttWindowSize)) < options%atomLatSp) then
       print '(a)', 'WARNING: Setting stt_window_size automatically since the user did not specify one'
       options%sttWindowSize = 4.0 * options%atomLatSp
    end if
    
    if (failed) stop 22
  end subroutine validate_inputedInformation

  !> Initialises the list of fields accepted in the configuration file.
  !! Associates a MultiscaleOptions structure to the fields.
  !! @param[inout] fields Fields structure to be initialised. Must be preallocated.
  !! @param[in] opts Pointer. Options structure whose members will be associated to fields of the configuration. 
  subroutine setup_inputfields(fields, opts, configFile, ierr, errMsg)
    use InteractionInput
    implicit none
    type(InputFields), intent(inout) :: fields
    type(MultiscaleOptions), intent(inout) :: opts

    character(len=*), intent(in) :: configFile
    integer,intent (out) :: ierr
    character(len=*), intent(out) :: errMsg
    
    call addField(fields, 'dimension', opts%space%spatDimension, .true., &
         "Syntax: dimension (1 | 2 | 3)" // NEW_LINE('A') // "Specifies the dimensionality of the simulation.")
    call addField(fields, 'damping_band_width', opts%dampingbandWidth, .true., &
         "Syntax: damping_band_width <nonnegative real>" // NEW_LINE('A') // "Specifies the width of the damping band," // &
         " the width is measured from the padding atoms.")
    call addField(fields, 'coarse_grained_width', opts%coarseGrainedWidth, .true.,&
         "Syntax: coarse_grained_width <nonnegative real>" // NEW_LINE('A') // &
         "Specifies the width of the coarse grained region, the width is measured from the padding atoms.")
    call addField(fields, 'part_coarse_grained_width', opts%partCoarseGrainedWidth, .true., &
         "Syntax: part_coarse_grained_width <nonnegative real>" // NEW_LINE('A') // &
         "Specifies the width of the partially coarse grained region, the width is measured from the padding atoms.")
    call addField(fields, 'padding_width', opts%paddingWidth, .true., &
         "Syntax: padding_width <nonnegative real>" // NEW_LINE('A') // &
         "Specifies the width of the padding region, the width is measured from the closest atom.")
    call addField(fields, 'periodic_boundary', opts%space%periodicBoundary, .true., &
         "Syntax: periodic_boundary (T | F) (T | F) (T | F)" // NEW_LINE('A') // &
         "Specifies if the boundary should be periodic or not, "// &
         "T indicates a periodic boundary, F indicates a vaccuum boundary.")
    call addField(fields, 'universe_size', opts%space%universeSize, .true., &
         "Syntax: universe_size <x> <y> <z>" // NEW_LINE('A') // &
         "Specifies the size of the space that the simulation takes place in.")
    call addField(fields, 'finitediff_boxes', opts%finiteDiffBoxes, .true., &
         "Syntax: finitediff_boxes <x> <y> <z>" // NEW_LINE('A') // &
         "The number of finite difference boxes in the direction of each dimension. "//&
         "Boxes in three dimensions are boxes that have eight finite difference grid points as its corners.")
    call addField(fields, 'atomistic_shape', atomisticShapeListParser, .true., .false., &
         "Syntax: atomistic_shape [file <filename>]" // NEW_LINE('A') // &
         "The shapes that defines the atomistic region.")
    call addField(fields, 'hole_shape', holeShapeListParser, .false., .false., &
         "Syntax: hole_shape [file <filename>]" // NEW_LINE('A') // &
         "Defines holes in the atomistic region.")
    call addField(fields, 'unitcell_atoms', unitcellParser, .true., .false., &
         "Syntax: unitcell_atoms <x> <y> <z>" // NEW_LINE('A') // &
         "The atoms that are placed in the unit cell.")
    call addField(fields, 'exchange_atoms', realExchangeParser, .true., .false., &
         "Syntax: exchange_atoms [unitcell] [file <filename>]" // NEW_LINE('A') // &
         "The exchange between the atoms.")
    call addField(fields, 'exchange_coarse', coarseExchangeParser, .true., .false., &
         "Syntax: exchange_coarse [unitcell] [file <filename>]" // NEW_LINE('A') // &
         "The exchange between the coarse grained atoms.")

    call addField(fields, 'dm_atoms', realDmParser, .false., .false., &
         "Syntax: dm_atoms [unitcell] [file <filename>]" // NEW_LINE('A') // &
         "The Dzyaloshinsky-Moriya interaction vectors between the atoms.")
    call addField(fields, 'dm_coarse', coarseDmParser, .false., .false., &
         "Syntax: dm_coarse [unitcell] [file <filename>]" // NEW_LINE('A') // &
         "The Dzyaloshinsky-Moriya interaction vectors between the coarse grained atoms.")

    call addField(fields, 'atom_lattice_spacing', opts%atomLatSp, .true., &
         "Syntax: atom_lattice_spacing <postive real>" // NEW_LINE('A') // &
         "This parameter will be use to approximate the side of the Wigner-Seitz cell around each atom." // NEW_LINE('A') // &
         "Each atom will be enclosed by a box of side the value of this parameter." // NEW_LINE('A') // &
         "The approximation is used to calculate weights for averaging padding atoms and interface nodes" // NEW_LINE('A') // &
         "and in the coefficients used in the damping band. The optimal value depends on the lattice.")
    call addField(fields, 'continuum_exchange_coef', continuumExchangeParser, .true., .false., &
         "Syntax: continuum_exchange_coef <real> <real> <real>" // NEW_LINE('A') // &
         "Specifies the continuum exchange coefficient. Equal components can be omitted.")
    call addField(fields, 'damping_band_window_size', opts%windowSize, .true., &
         "Syntax: damping_band_window_size <real> <real> <real>" // NEW_LINE('a') // &
         "Specifies the window size for the dampingband.")
    call addField(fields, 'continuum_moment_magnitude', opts%continuumMomentMagnitude, .true., &
         "Syntax: continuum_moment_magnitude <positive real>" // NEW_LINE('a') // &
         "Specifies the magnitude of the moments in the continuum region")

    call addField(fields, 'anisotropy', anisotropyParser, .false., .true., &
         "Syntax: anisotropy (uniaxial|cubic|both <ratio>) K1 K2 ex ey ez [file <filename>]"&
         // NEW_LINE('A') // &
         "Specifies regions that are affected by the given anisotropy.")

    call addField(fields, 'moments', momentRegionsParser, .false., .true., &
         "Syntax: moments Magnitude mx my mz [file <filename>]"&
         // NEW_LINE('A') // &
         "Specifies regions that will have the specified magnetic moment. "&
         // NEW_LINE('A') // &
         "These regions take precedence over the moment specified in the unit cell.")

    call addField(fields, 'link_error_tolerance', opts%linkErrorTolerance, .false., &
         "Syntax: link_error_tolerance <positive real>" &
         // NEW_LINE('A') // &
         "Optional. Specifies the error tolerance when linking two atoms using a directional vector." &
         // NEW_LINE('A') // &
         "Defaults to 1d-5, but should be adjusted if you require other precission.")

    call addField(fields, 'continuum_dm', continuumDmParser, .false., .false., &
         "Syntax: continuum_dm dx1 dx2 dx3 / dy1 dy2 dy3 / dz1 dz2 dz3" &
         // NEW_LINE('A') // &
         "Optional. Specifies the Dzyaloshinsky-Moriya interaction vectors for the continuum." &
         // NEW_LINE('A') // &
         "If a vector is not specified, its value defaults to 0.")
         
    call addField(fields, 'damping_band_strength', opts%dampingBandStrength, .false., &
         "Syntax: damping_band_strength <positive real>" &
         // NEW_LINE('A') // &
         "Optional. Specifies the strength parameter of the damping band." &
         // NEW_LINE('A') // &
         "When 0, no damping band coefficients are calculated. Defaults to 0.")

    call addField(fields, 'define_zone', atomZoneParser, .false., .true., &
         "Syntax: define_zone <zone id> [file <filename>]" // NEW_LINE('A') // &
         "Defines regions that, in the event of being intersected by atomistic regions," &
         // NEW_LINE('A') // &
         "are filled with the atoms labeled by the given zone id." &
         // NEW_LINE('A') // &
         "Special care must be put so that these regions don't contact the continuum," &
         // NEW_LINE('A') // &
         "since atoms in the continuum are assumed to belong to the default zone.")

    call addField(fields, 'stt_vector', opts%sttVector, .false., &
         "Syntax: stt_vector <x> <y> <z>" // NEW_LINE('A') // &
         "Specifies the direction and strength of the STT.")
    call addField(fields, 'stt_window_size', opts%sttWindowSize, .false., &
         "Syntax: stt_window_size <x> <y> <z>" // NEW_LINE('A') // &
         "Specifies the size of the window used to approximate the directional" &
         // NEW_LINE('A') // &
         "gradients for atomistic spin transfer torque. " &
         // NEW_LINE('A') // &
         "When unspecified 4*lattice_spacing is used." &
         // NEW_LINE('A') // &
         "NaN weights could be introduced due to having a too small window.")
    
    call loadInputFile(trim(configFile), fields, ierr, errMsg)

  contains

    subroutine unitcellParser(fData)
      use MultiscaleFileParser
      implicit none
      type(FileData), intent(inout) :: fData

      call parseAtomCell(fData, opts%unitCell)
    end subroutine unitcellParser

    subroutine atomZoneParser(fData)
      use MultiscaleFileParser
      use AtomGenerator, only: parseZoneList
      implicit none
      type(FileData), intent(inout) :: fData

      call parseZoneList(fData, opts%atomZones)
    end subroutine atomZoneParser
    
    subroutine atomisticShapeListParser(fData)
      use MultiscaleFileParser
      implicit none
      type(FileData), intent(inout) :: fData

      call parseShapeList(fData, opts%atomisticShapes)
    end subroutine atomisticShapeListParser

    subroutine holeShapeListParser(fData)
      use MultiscaleFileParser
      implicit none
      type(FileData), intent(inout) :: fData

      call parseShapeList(fData, opts%holeShapes)
    end subroutine holeShapeListParser

    subroutine realExchangeParser(fData)
      use MultiscaleFileParser
      implicit none
      type(FileData), intent(inout) :: fData

      call parseInteractionList(fData, 1, opts%realExchange) 
    end subroutine realExchangeParser

    subroutine coarseExchangeParser(fData)
      use MultiscaleFileParser
      implicit none
      type(FileData), intent(inout) :: fData

      call parseInteractionList(fData, 1, opts%coarseExchange) 
    end subroutine coarseExchangeParser


    subroutine realDmParser(fData)
      use MultiscaleFileParser
      implicit none
      type(FileData), intent(inout) :: fData

      call parseInteractionList(fData, 3, opts%realDm) 
    end subroutine realDmParser

    subroutine coarseDmParser(fData)
      use MultiscaleFileParser
      implicit none
      type(FileData), intent(inout) :: fData

      call parseInteractionList(fData, 3, opts%coarseDm) 
    end subroutine coarseDmParser



    subroutine anisotropyParser(fData)
      use MultiscaleFileParser
      use Anisotropy
      implicit none
      type(FileData), intent(inout) :: fData

      call parseAnisotropyList(fData, opts%anisotropies) 
    end subroutine anisotropyParser

    subroutine momentRegionsParser(fData)
      use MultiscaleFileParser
      use MomentRegions
      implicit none
      type(FileData), intent(inout) :: fData

      call parseMomentList(fData, opts%momentRegions) 
    end subroutine momentRegionsParser


    subroutine continuumExchangeParser(fData)
      use MultiscaleFileParser
      implicit none
      type(FileData), intent(inout) :: fData
      
      integer :: component, i
      
      do component = 1,3
         if (isAtEndOfLine(fData)) then
            if (component .eq. 1) then               
               call createErrorMsg(fData,1,'At least one exchange component must be specified')
               return
            else
               do i = component,3
                  opts%continuumExchangeCoef(i) = &
                       opts%continuumExchangeCoef(component-1)
               end do
               call readNextLine(fData)
               return
            end if
         else
            call parseReal(fData, opts%continuumExchangeCoef(component))
            if (fData%ierr /= 0) return
            call readNextWord(fData)
         end if
      end do      
      call readNextLine(fData)

    end subroutine continuumExchangeParser
    
    subroutine continuumDmParser(fData)
      use MultiscaleFileParser
      implicit none
      type(FileData), intent(inout) :: fData

      integer :: vector, component

      do vector = 1, 3
         do component = 1,3
            if (isAtEndOfLine(fData)) then
               call createErrorMsg(fData,1,'Incomplete DM vector')
               return
            else
               call parseReal(fData, opts%continuumDm(component,vector))
               if (fData%ierr /= 0) return
               call readNextWord(fData)
            end if
         end do
         if (isAtEndOfLine(fData)) then
            call readNextLine(fData)
            return
         end if
         call skipOptionalDelimiter(fData)
      end do
      
      call readNextLine(fData)

    end subroutine continuumDmParser

  end subroutine setup_inputfields
end module Configuration
