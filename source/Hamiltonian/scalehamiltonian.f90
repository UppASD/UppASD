module ScaleHamiltonian

    use Parameters
    use Constants

    implicit none

    real(dblprec), dimension(:,:,:), allocatable :: ncoup_orig !< Original coupling constant before scaling
    real(dblprec), dimension(:,:,:,:), allocatable :: j_tens_orig !< Original coupling constant before scaling
    real(dblprec), dimension(:,:,:), allocatable :: dm_vect_orig !< Original DMI vectors before scaling
    real(dblprec), dimension(:), allocatable :: jscale_factor !< Scaling matrix for Hamiltonian
    real(dblprec), dimension(:,:), allocatable :: dmiscale_factor !< Component-wise scaling factors for DMI (3 x N)
    logical :: jscaling_flag  !< Flag to indicate if the scaling matrix is used for J_ij
    logical :: jscaling_dynamic  !< Flag to indicate dynamic scaling
    logical :: jscaling_dmi  !< Flag to indicate if DMI should be scaled as well
    integer :: jscaling_natoms !< Number of atoms in the scaling file
    integer, dimension(:), allocatable :: jscaling_atomlist !< List of atoms in the scaling file
    character(len=30) :: jscaling_file !< Name of the scaling file
    real(dblprec) :: jscaling_freq !< Frequency of dynamic scaling (in THz)
    real(dblprec) :: jscaling_phase !< Phase of dynamic scaling (in radians)
    real(dblprec) :: jscaling_prefactor !< Prefactor for dynamic scaling

    contains

    subroutine read_jscaling_file()
        use InputData, only: Natom
    
        !< Read scaling parameters from the input file
        integer :: i, j, n, m, idum, i_err
        real(dblprec) :: temp

        if (.not. jscaling_dynamic) then
            ! Read scaling file
            open(unit=ifileno2, file=jscaling_file, status='old', action='read')
            jscaling_natoms = 0
            do
                read(ifileno2, *, iostat=i_err) temp
                if (i_err /= 0) exit
                jscaling_natoms = jscaling_natoms + 1 
            end do
            rewind(ifileno2)
            ! Initialize scaling matrix
            allocate(jscale_factor(jscaling_natoms))
            jscale_factor = 1.0_dblprec
            if (jscaling_dmi) then
                allocate(dmiscale_factor(3,jscaling_natoms))
                dmiscale_factor = 1.0_dblprec
            end if

            ! Read scaling file
            open(unit=ifileno2, file=jscaling_file, status='old', action='read')
            if (.not. jscaling_dmi) then
                do i = 1, jscaling_natoms
                    read(ifileno2, *) jscale_factor(i)
                end do
            else
                ! When jscaling_dmi is true we expect each line to contain three DMI component scale factors
                ! Format assumption (static mode): dm_x dm_y dm_z  (one row per atom)
                do i = 1, jscaling_natoms
                    read(ifileno2, *) jscale_factor(i), dmiscale_factor(1,i), dmiscale_factor(2,i), dmiscale_factor(3,i)
                end do
            end if
            close(ifileno2)
        else
            ! Read scaling file
            open(unit=ifileno2, file=jscaling_file, status='old', action='read')
            jscaling_natoms = 0
            do
                read(ifileno2, *, iostat=i_err) idum
                if (i_err /= 0) exit
                jscaling_natoms = jscaling_natoms + 1 
            end do
            rewind(ifileno2)
            print *,'Jscaling file read with ', jscaling_natoms, ' atoms'

            ! Initialize scaling matrix
            allocate(jscaling_atomlist(jscaling_natoms))
            allocate(jscale_factor(jscaling_natoms))
            jscale_factor = 1.0_dblprec
            if (jscaling_dmi) then
                allocate(dmiscale_factor(3,jscaling_natoms))
                dmiscale_factor = 1.0_dblprec
            end if
            ! Read scaling file
            if (.not. jscaling_dmi) then
                do i = 1, jscaling_natoms
                    read(ifileno2, *) jscaling_atomlist(i), jscale_factor(i)
                end do
            else
                ! When jscaling_dmi is true we expect each line to contain: id dm_x dm_y dm_z jscale
                ! Format assumption (dynamic mode): id dm_x dm_y dm_z jscale
                do i = 1, jscaling_natoms
                    read(ifileno2, *) jscaling_atomlist(i), jscale_factor(i), dmiscale_factor(1,i), dmiscale_factor(2,i), dmiscale_factor(3,i)
                end do
            end if
            close(ifileno2)
        end if

    end subroutine read_jscaling_file
    
    subroutine apply_local_jscaling()
        use InputData, only: Natom, do_reduced, ham_inp
        use HamiltonianData, only : ham
    
        implicit none
        integer :: i_atom, j_neigh, j_atom, i_neigh
    
        real(dblprec), dimension(:,:), allocatable :: temp_ncoup

        if (do_reduced == 'Y') then
            write(*,*) 'Applying J scaling to reduced Hamiltonian not possible for `do_reduced` = Y'
            return
        end if

        if (ham_inp%do_jtensor/=1) then
            if (.not. allocated(ncoup_orig)) then
                allocate(ncoup_orig(size(ham%ncoup, 1), size(ham%ncoup, 2), size(ham%ncoup, 3)))
                ncoup_orig = ham%ncoup  ! Store original coupling matrix
            end if
    
            ! Allocate temporary array for scaled coupling matrix
            ! allocate(temp_ncoup(size(ham%ncoup, 1), size(ham%ncoup, 2)))
            ! temp_ncoup = ham%ncoup(:,:, 1)  ! Copy original coupling matrix
            ham%ncoup = ncoup_orig  ! Reset to original coupling matrix
    
            ! Apply scaling to the coupling matrix
            ! First scale all couplings Jij for each given atom
            do i_atom = 1, Natom
                ham%ncoup(:, i_atom, 1) = ham%ncoup(:, i_atom, 1) * jscale_factor(i_atom)
                ! Then loop over the neighbors to apply the scaling
                ! to the J_ji couplings
                do j_neigh=1,ham%nlistsize(i_atom)
                    j_atom = ham%nlist(j_neigh, i_atom)
                    do i_neigh=1,ham%nlistsize(j_atom)
                        if (i_atom == ham%nlist(i_neigh, j_atom)) then
                            ham%ncoup(i_neigh, j_atom, 1) = ham%ncoup(i_neigh, j_atom, 1) * jscale_factor(i_atom)
                        end if
                    end do
                end do
            end do
    
            ! deallocate(temp_ncoup)
        else ! Tensor Jij
            if (.not. allocated(j_tens_orig)) then
                allocate(j_tens_orig(size(ham%j_tens, 1), size(ham%j_tens, 2), size(ham%j_tens, 3), size(ham%j_tens, 4)))
                j_tens_orig = ham%j_tens  ! Store original tensor coupling matrix
            end if

            ham%j_tens = j_tens_orig  ! Reset to original tensor coupling matrix

            ! Apply scaling to the tensor coupling matrix
            ! First scale all couplings Jij for each given atom
            do i_atom = 1, Natom
                ham%j_tens(:, :, :, i_atom) = ham%j_tens(:, :, :, i_atom) * jscale_factor(i_atom)
                ! Then loop over the neighbors to apply the scaling
                ! to the J_ji couplings
                do j_neigh=1,ham%nlistsize(i_atom)
                    j_atom = ham%nlist(j_neigh, i_atom)
                    do i_neigh=1,ham%nlistsize(j_atom)
                        if (i_atom == ham%nlist(i_neigh, j_atom)) then
                            ham%j_tens(:, :, i_neigh, j_atom) = ham%j_tens(:, :, i_neigh, j_atom) * jscale_factor(i_atom)
                        end if
                    end do
                end do
            end do
    
            ! deallocate(temp_ncoup)
        end if

            ! Apply scaling to DMI vectors if requested
            if (jscaling_dmi .and. ham_inp%do_dm==1 .and. allocated(ham%dm_vect)) then
                if (.not. allocated(dm_vect_orig)) then
                    allocate(dm_vect_orig(size(ham%dm_vect,1), size(ham%dm_vect,2), size(ham%dm_vect,3)))
                    dm_vect_orig = ham%dm_vect  ! Store original DMI vectors
                end if

                ham%dm_vect = dm_vect_orig  ! Reset to original DMI vectors

                ! Loop over atoms and their DM neighbours and apply local scaling
                do i_atom = 1, Natom
                    do i_neigh = 1, ham%dmlistsize(i_atom)
                        j_atom = ham%dmlist(i_neigh, i_atom)
                        ! Scale the D_ij vector corresponding to (i_atom -> j_atom)
                        if (allocated(dmiscale_factor)) then
                            ham%dm_vect(1, i_neigh, i_atom) = ham%dm_vect(1, i_neigh, i_atom) * dmiscale_factor(1,i_atom)
                            ham%dm_vect(2, i_neigh, i_atom) = ham%dm_vect(2, i_neigh, i_atom) * dmiscale_factor(2,i_atom)
                            ham%dm_vect(3, i_neigh, i_atom) = ham%dm_vect(3, i_neigh, i_atom) * dmiscale_factor(3,i_atom)
                        else
                            ham%dm_vect(:, i_neigh, i_atom) = ham%dm_vect(:, i_neigh, i_atom) * jscale_factor(i_atom)
                        end if
                        ! Find the corresponding reverse entry D_ji and scale it
                        do j_neigh = 1, ham%dmlistsize(j_atom)
                            if (i_atom == ham%dmlist(j_neigh, j_atom)) then
                                if (allocated(dmiscale_factor)) then
                                    ham%dm_vect(1, j_neigh, j_atom) = ham%dm_vect(1, j_neigh, j_atom) * dmiscale_factor(1,i_atom)
                                    ham%dm_vect(2, j_neigh, j_atom) = ham%dm_vect(2, j_neigh, j_atom) * dmiscale_factor(2,i_atom)
                                    ham%dm_vect(3, j_neigh, j_atom) = ham%dm_vect(3, j_neigh, j_atom) * dmiscale_factor(3,i_atom)
                                else
                                    ham%dm_vect(:, j_neigh, j_atom) = ham%dm_vect(:, j_neigh, j_atom) * jscale_factor(i_atom)
                                end if
                            end if
                        end do
                    end do
                end do
            end if
    
    end subroutine apply_local_jscaling


    subroutine apply_dynamic_jscaling(time)
        use InputData, only: Natom, do_reduced, ham_inp
        use HamiltonianData, only : ham
    
        implicit none
    
        real(dblprec), intent(in) :: time
    
    
        integer :: i_atom, j_neigh, j_atom, i_neigh, idx
        real(dblprec) :: scale_factor

        if (do_reduced == 'Y') then
            write(*,*) 'Applying J scaling to reduced Hamiltonian not possible for `do_reduced` = Y'
            return
        end if

        ! real(dblprec), dimension(:,:), allocatable :: temp_ncoup
    
        if (ham_inp%do_jtensor/=1) then
    
        if (.not. allocated(ncoup_orig)) then
            !print *,'Storing copy of ham%ncoup, shape:', shape(ham%ncoup)
            allocate(ncoup_orig(size(ham%ncoup, 1), size(ham%ncoup, 2), size(ham%ncoup, 3)))
            !print *,'ncoup_orig stored, shape:', shape(ncoup_orig)
            ncoup_orig = ham%ncoup  ! Store original coupling matrix
        end if
    
        ! Allocate temporary array for scaled coupling matrix
        ! allocate(temp_ncoup(size(ham%ncoup, 1), size(ham%ncoup, 2)))
        ! temp_ncoup = ham%ncoup(:,:, 1)  ! Copy original coupling matrix
        scale_factor = jscaling_prefactor * sin(2.0_dblprec*pi*jscaling_freq * time + jscaling_phase)
        write(1000,*) 'Time: ', time, ' Scale factor: ', scale_factor, jscaling_natoms
        ! print *, 'Time: ', time, ' Scale factor: ', scale_factor, jscaling_natoms
        ! Apply scaling to the coupling matrix
        ! First scale all couplings Jij for each given atom
        do idx = 1, jscaling_natoms
            ! print *,idx, 'Applying J scaling to atom:', jscaling_atomlist(idx), scale_factor
            i_atom = jscaling_atomlist(idx)
            ham%ncoup(:, i_atom, 1) = ncoup_orig(:, i_atom, 1) * ( 1.0_dblprec + scale_factor * jscale_factor(idx))
            ! print *, 'Scaled J for atom ', i_atom, ' to: ', ham%ncoup(:, i_atom, 1)
            ! Then loop over the neighbors to apply the scaling
            ! to the J_ji couplings
            do j_neigh=1,ham%nlistsize(i_atom)
                j_atom = ham%nlist(j_neigh, i_atom)
                do i_neigh=1,ham%nlistsize(j_atom)
                    if (i_atom == ham%nlist(i_neigh, j_atom)) then
                                    ham%ncoup(:, i_atom, 1) = ncoup_orig(:, i_atom, 1) * ( 1.0_dblprec + scale_factor * jscale_factor(idx))
                    end if
                end do
            end do
        end do
        ! Dynamic scaling for DMI
        if (jscaling_dmi .and. ham_inp%do_dm==1 .and. allocated(ham%dm_vect)) then
            if (.not. allocated(dm_vect_orig)) then
                allocate(dm_vect_orig(size(ham%dm_vect,1), size(ham%dm_vect,2), size(ham%dm_vect,3)))
                dm_vect_orig = ham%dm_vect  ! Store original DMI vectors
            end if

            do idx = 1, jscaling_natoms
                i_atom = jscaling_atomlist(idx)
                do i_neigh = 1, ham%dmlistsize(i_atom)
                    j_atom = ham%dmlist(i_neigh, i_atom)
                    if (allocated(dmiscale_factor)) then
                        ! print *, 'Scaling DMI for atom ', i_atom, ' with factors: ', dmiscale_factor(:,idx), ' and jscale: ', jscale_factor(idx)
                        ham%dm_vect(1, i_neigh, i_atom) = dm_vect_orig(1, i_neigh, i_atom) * ( 1.0_dblprec + scale_factor * dmiscale_factor(1,idx))
                        ham%dm_vect(2, i_neigh, i_atom) = dm_vect_orig(2, i_neigh, i_atom) * ( 1.0_dblprec + scale_factor * dmiscale_factor(2,idx))
                        ham%dm_vect(3, i_neigh, i_atom) = dm_vect_orig(3, i_neigh, i_atom) * ( 1.0_dblprec + scale_factor * dmiscale_factor(3,idx))
                        ! print *, 'Scaled DMI for atom ', i_atom, ' to: ', ham%dm_vect(:, i_neigh, i_atom)
                    else
                        ham%dm_vect(:, i_neigh, i_atom) = dm_vect_orig(:, i_neigh, i_atom) * ( 1.0_dblprec + scale_factor * jscale_factor(idx))
                    end if
                    ! Find and scale reverse entry
                    do j_neigh = 1, ham%dmlistsize(j_atom)
                        if (i_atom == ham%dmlist(j_neigh, j_atom)) then
                            if (allocated(dmiscale_factor)) then
                                ! ham%dm_vect(1, j_neigh, j_atom) = dm_vect_orig(1, j_neigh, j_atom) * ( 1.0_dblprec + scale_factor * dmiscale_factor(1,i_atom))
                                ! ham%dm_vect(2, j_neigh, j_atom) = dm_vect_orig(2, j_neigh, j_atom) * ( 1.0_dblprec + scale_factor * dmiscale_factor(2,i_atom))
                                ! ham%dm_vect(3, j_neigh, j_atom) = dm_vect_orig(3, j_neigh, j_atom) * ( 1.0_dblprec + scale_factor * dmiscale_factor(3,i_atom))
                                ham%dm_vect(1, j_neigh, j_atom) = dm_vect_orig(1, j_neigh, j_atom) * ( 1.0_dblprec + scale_factor * dmiscale_factor(1,idx))
                                ham%dm_vect(2, j_neigh, j_atom) = dm_vect_orig(2, j_neigh, j_atom) * ( 1.0_dblprec + scale_factor * dmiscale_factor(2,idx))
                                ham%dm_vect(3, j_neigh, j_atom) = dm_vect_orig(3, j_neigh, j_atom) * ( 1.0_dblprec + scale_factor * dmiscale_factor(3,idx))
                            else
                                ham%dm_vect(:, j_neigh, j_atom) = dm_vect_orig(:, j_neigh, j_atom) * ( 1.0_dblprec + scale_factor * jscale_factor(idx))
                            end if
                        end if
                    end do
                end do
            end do
        end if
        else

        if (.not. allocated(j_tens_orig)) then
            allocate(j_tens_orig(size(ham%j_tens, 1), size(ham%j_tens, 2), size(ham%j_tens, 3), size(ham%j_tens, 4)))
            j_tens_orig = ham%j_tens  ! Store original tensor coupling matrix
        end if

        scale_factor = 1.0_dblprec + jscaling_prefactor * sin(2.0_dblprec*pi*jscaling_freq * time + jscaling_phase)
        write(1000,*) 'Time: ', time, ' Scale factor: ', scale_factor, jscaling_natoms
        ! Apply scaling to the tensor coupling matrix
        ! First scale all couplings Jij for each given atom
        do idx = 1, jscaling_natoms
            i_atom = jscaling_atomlist(idx)
            ham%j_tens(:, :, :, i_atom) = j_tens_orig(:, :, :, i_atom) * scale_factor
            ! Then loop over the neighbors to apply the scaling
            ! to the J_ji couplings
            do j_neigh=1,ham%nlistsize(i_atom)
            j_atom = ham%nlist(j_neigh, i_atom)
            do i_neigh=1,ham%nlistsize(j_atom)
                if (i_atom == ham%nlist(i_neigh, j_atom)) then
                ham%j_tens(:, :, i_neigh, j_atom) = j_tens_orig(:, :, i_neigh, j_atom) * scale_factor
                end if
            end do
            end do
        end do
        end if
    
    end subroutine apply_dynamic_jscaling
    
    subroutine apply_jscaling()
        implicit none
        !< Read and apply the J scaling factors from the file
        ! write(*, '(a,a)', advance='no') 'Reading J local scaling factors from file: ', trim(jscaling_file)
        ! call read_jscaling_file()
        ! write(*,'(a)') ' done'
        write(*,'(a)', advance='no') 'Applying local J scaling factors to Hamiltonian...'
        call apply_local_jscaling()
        write(*,'(a)') ' done'
    end subroutine apply_jscaling
    
    subroutine jscaling_init()
       !
       implicit none
       !
        jscaling_flag = .false.
        jscaling_dynamic = .false.
        jscaling_dmi = .false.
        jscaling_file = 'jscaling.dat'
        jscaling_freq = 0.0
        jscaling_phase = 0.0

    end subroutine jscaling_init
   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_jscaling(ifile)

      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword, cache
      integer :: rd_len, i_err, i_errb, i_stat, ii
      logical :: comment

      print *,'Reading j-scaling data '
      do
         10     continue
         ! Read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,ifile,i_errb)

         ! converting Capital letters
         call caps2small(keyword)

         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
            (scan(trim(keyword),'!')==1))

         if (comment) then
            read(ifile,*)
         else
            ! Parse keyword
            keyword=trim(keyword)
            select case(keyword)

        case('jscaling_flag')
           read(ifile,*,iostat=i_err) jscaling_flag
           print *,'jscaling_flag:', jscaling_flag
           if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

        case('jscaling_dynamic')
           read(ifile,*,iostat=i_err) jscaling_dynamic
           print *,'jscaling_dynamic:', jscaling_dynamic
           if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('jscaling_dmi')
              read(ifile,*,iostat=i_err) jscaling_dmi
              print *,'jscaling_dmi:', jscaling_dmi
              if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

        case('jscaling_prefactor')
           read(ifile,*,iostat=i_err) jscaling_prefactor
           print *,'jscaling_prefactor:', jscaling_prefactor
           if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

        case('jscaling_freq')
           read(ifile,*,iostat=i_err) jscaling_freq
           print *,'jscaling_freq:', jscaling_freq
           if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

        case('jscaling_phase')
           read(ifile,*,iostat=i_err) jscaling_phase
           print *,'jscaling_phase:', jscaling_phase
           if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

        case('jscaling_file')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            jscaling_file=trim(adjustl(cache))
            print *,'jscaling_file:', trim(jscaling_file)

            call read_jscaling_file()

         end select
      end if

      ! End of file
      if (i_errb==20) goto 20
      ! End of row
      if (i_errb==10) goto 10
   end do

   20  continue

      print *,'Finished reading j-scaling data '
      print *,'jscaling_file:', trim(jscaling_file)
      print *,'jscaling_prefactor:', jscaling_prefactor
      print *,'jscaling_freq:', jscaling_freq
      print *,'jscaling_phase:', jscaling_phase
   rewind(ifile)
   return
end subroutine read_parameters_jscaling

end module ScaleHamiltonian
