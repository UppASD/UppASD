  !> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro
  module MultiscaleDampingBand
  use Parameters
  use Profiling
  use MultiscaleInterpolation

  type DampingBandData
     type(LocalInterpolationInfo)             :: interpolation
     real(dblprec), dimension(:), pointer     :: coefficients
     real(dblprec), dimension(:,:,:), pointer :: preinterpolation
     logical                                  :: enable
  end type DampingBandData

  type(DampingBandData) :: dampingBand
  
contains 


  !! Update copy of emom used for calculating local moment averages
  subroutine update_backbuffer(dband)
    use MomentData, only : emom
    use Multiscale, only : multiscaleBackbuffer, multiscaleBackbufferHead
    implicit none
    type(DampingBandData), intent(in) :: dband
    integer :: j,k
    if(dband%enable) then
       multiscaleBackbufferHead = multiscaleBackbufferHead+1
       if (multiscaleBackbufferHead > ubound(multiscaleBackbuffer,4)) &
            multiscaleBackbufferHead = 1

       !! Parallel data copy is only faster if caches happen to help,
       !!  maybe this should not be parallelised
       !$omp parallel do default(shared) private(k,j) collapse(2)
       do k=1, ubound(emom,3)
          do j=1, ubound(emom,2)
             multiscaleBackbuffer(:,j,k,multiscaleBackbufferHead) =&
                  emom(:,j,k)
          end do
       end do
       !$omp end parallel do
    end if
  end subroutine update_backbuffer

  
  !! Initialize damping band 
  subroutine initDampingBand(dband)      
    implicit none
    type(DampingBandData),intent(inout) :: dband

    dband%enable = .false.
    nullify(dband%coefficients)
    nullify(dband%preinterpolation)
    call newLocalInterpolationInfo(dband%interpolation)

  end subroutine initDampingBand

  !! Deinitialize damping band 
  !! when needed deallocates memory.
  subroutine deallocateDampingBand(dband)
    implicit none
    type(DampingBandData),intent(inout) :: dband
    integer :: var_size, i_stat

    dband%enable = .false.
    if(associated(dband%coefficients)) then
       var_size = product(shape(dband%coefficients))*kind(dband%coefficients)
       deallocate(dband%coefficients, stat=i_stat)
       call memocc(i_stat, -var_size, &
            'DampingBandData%coefficients','deallocateDampingBand')
       nullify(dband%coefficients)
    end if
    if(associated(dband%preinterpolation)) then
       var_size = product(shape(dband%preinterpolation))*kind(dband%preinterpolation)
       deallocate(dband%preinterpolation, stat=i_stat)
       call memocc(i_stat, -var_size, &
            'DampingBandData%preinterpolation','deallocateDampingBand')
       nullify(dband%coefficients)
    end if
    call deleteLocalInterpolationInfo(dband%interpolation)  
  end subroutine deallocateDampingBand

  
  !! Calculate the damping vector used for each atom in the damping band
  !! This is done by applying a local interpolation of the moments of the atoms
  !! in each damping atom.
  !! ToDo: avoid iteration over all atoms
  subroutine dampingBandPreinterpolation(dband,emom)
    implicit none
    type(DampingBandData),intent(inout) :: dband
    real(dblprec),dimension(:,:,:),intent(in) :: emom
    integer :: dampingBandIndex,j,k

    if (dband%enable) then
       do k=1,ubound(dband%preinterpolation,3)
          do j=1,ubound(dband%interpolation%indices,1)
             dampingBandIndex = dband%interpolation%indices(j)
             if (dampingBandIndex .ne. 0) then
                dband%preinterpolation(:,dampingBandIndex,k) =&
                     atomInterpolation(dband%interpolation,j,k,emom)
             end if
          end do
       end do
    end if
  end subroutine dampingBandPreinterpolation

  
  !! Calculate the damping vector used in the corrector step.
  !! The damping vectors are the mean of the previous damping vectors and
  !! the damping vectors calculated in for the predictor step.
  !! ToDo: avoid iteration over all atoms
  subroutine dampingBandCorrPreinterpolation(dband,emom2)
    implicit none
    type(DampingBandData),intent(inout) :: dband
    real(dblprec),dimension(:,:,:),intent(in) :: emom2
    integer :: dampingBandIndex,j,k
    if (dband%enable) then
       do k=1,ubound(dband%preinterpolation,3)
          do j=1,ubound(dband%preinterpolation,2)
             dampingBandIndex = dband%interpolation%indices(j)
             if (dampingBandIndex .ne. 0) then
                dband%preinterpolation(:,dampingBandIndex,k) = &
                     atomInterpolation(dband%interpolation,j,k,emom2)
             end if
          end do
       end do
    end if
  end subroutine dampingBandCorrPreinterpolation


  !! Calculate the damping vector used in the corrector step.
  !! This function averages the predictor interpolation with the current.
  !! It should be used for methods that predict half timestep ahead.
  !! The damping vectors are the mean of the previous damping vectors and
  !! the damping vectors calculated in for the predictor step.
  !! ToDo: avoid iteration over all atoms
  subroutine dampingBandCorrPreinterpolationAvg(dband,emom2)
    implicit none
    type(DampingBandData),intent(inout) :: dband
    real(dblprec),dimension(:,:,:),intent(in) :: emom2
    integer :: dampingBandIndex,j,k
    if (dband%enable) then
       do k=1,ubound(dband%preinterpolation,3)
          do j=1,ubound(dband%preinterpolation,2)
             dampingBandIndex = dband%interpolation%indices(j)
             if (dampingBandIndex .ne. 0) then
                dband%preinterpolation(:,dampingBandIndex,k) = &
                     (dband%preinterpolation(:,dampingBandIndex,k) + &
                     atomInterpolation(dband%interpolation,j,k,emom2))*5d-1
             end if
          end do
       end do
    end if
  end subroutine dampingBandCorrPreinterpolationAvg

                  

  
end module MultiscaleDampingBand
