!-------------------------------------------------------------------------------
! MODULE: de_qminimizer
!> @brief
!> Differential Evolution optimization for spin-spiral energy minimization
!> @author
!> Generated for Anders Bergman's qminimizer
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module de_minimizer
   use Parameters
   use Profiling
   use FieldData
   use HamiltonianActions
   use RandomNumbers, only: rng_uniform, rng_gaussianP
   
   implicit none
   
   private
   
   ! DE algorithm parameters
   integer :: de_pop_size = 50          !< Population size (NP), typically 10*n_params
   integer :: de_max_gen = 200          !< Maximum generations
   real(dblprec) :: de_F = 0.8_dblprec  !< Mutation factor (0.5-1.0)
   real(dblprec) :: de_CR = 0.9_dblprec !< Crossover probability (0.8-0.95)
   real(dblprec) :: de_tol = 1.0e-6_dblprec !< Convergence tolerance
   integer :: de_strategy = 1           !< DE strategy: 1=rand/1/bin, 2=best/1/bin
   character(len=1) :: de_constraint_mode = 'P' !< 'P'rojection or 'N'ormalization only
   logical :: de_qxy_plane = .false.    !< Constrain q-vectors to x-y plane (qz = 0)
   logical :: de_skyrmion_mode = .false. !< Skyrmion mode: optimize 1 q, generate 3 by 120° rotation
   
   ! Gradient refinement parameters
   logical :: de_use_local_search = .true. !< Use gradient-based local refinement
   integer :: de_local_steps = 20       !< Maximum local search steps per refinement
   integer :: de_local_freq = 10        !< Apply local search every N generations
   real(dblprec) :: de_local_tol = 1.0e-8_dblprec !< Local search tolerance
   real(dblprec) :: de_local_step = 0.01_dblprec  !< Initial step size for gradient descent
   
   ! Problem parameters
   integer :: de_n_qvec = 1             !< Number of q-vectors to optimize
   integer :: de_n_params               !< Total number of parameters
   logical :: de_verbose = .true.       !< Print progress information
   
   ! Best solution tracking
   real(dblprec), dimension(:), allocatable :: de_best_params
   real(dblprec) :: de_best_energy
   
   public :: de_minimize_qvectors, read_parameters_de, de_qminimizer_init
   public :: de_pop_size, de_max_gen, de_F, de_CR, de_n_qvec, de_skyrmion_mode
   
contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_minimize_qvectors
   !> @brief Main Differential Evolution minimization routine
   !> @author Generated
   !-----------------------------------------------------------------------------
   subroutine de_minimize_qvectors(Natom, Mensemble, NA, coord, emomM, mmom, &
      hfield, OPT_flag, max_no_constellations, maxNoConstl, unitCellType, &
      constlNCoup, constellations, constellationsNeighType, Num_macro, &
      cell_index, emomM_macro, macro_nlistsize, simid, q_result, s_result, &
      n_vec_result, energy_result)
      
      use InputData, only: N1, N2, N3
      use MomentData, only: emom
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: Natom, Mensemble, NA
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom
      real(dblprec), dimension(3), intent(in) :: hfield
      logical, intent(in) :: OPT_flag
      integer, intent(in) :: max_no_constellations
      integer, dimension(Mensemble), intent(in) :: maxNoConstl
      integer, dimension(Natom, Mensemble), intent(in) :: unitCellType
      real(dblprec), dimension(ham%max_no_neigh, max_no_constellations, Mensemble), intent(in) :: constlNCoup
      real(dblprec), dimension(3, max_no_constellations, Mensemble), intent(in) :: constellations
      integer, dimension(ham%max_no_neigh, max_no_constellations, Mensemble), intent(in) :: constellationsNeighType
      integer, intent(in) :: Num_macro
      integer, dimension(Natom), intent(in) :: cell_index
      real(dblprec), dimension(3, Num_macro, Mensemble), intent(in) :: emomM_macro
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize
      character(len=8), intent(in) :: simid
      
      ! Output parameters
      real(dblprec), dimension(3, 3), intent(out) :: q_result
      real(dblprec), dimension(3, 3), intent(out) :: s_result  ! Multiple s-vectors for skyrmion mode
      real(dblprec), dimension(3), intent(out) :: n_vec_result
      real(dblprec), intent(out) :: energy_result
      
      ! Local variables
      real(dblprec), dimension(:,:), allocatable :: population  ! (n_params, pop_size)
      real(dblprec), dimension(:), allocatable :: energies      ! (pop_size)
      real(dblprec), dimension(:), allocatable :: trial_params  ! (n_params)
      real(dblprec), dimension(:), allocatable :: mutant_params ! (n_params)
      real(dblprec), dimension(3, 3) :: q_temp, s_temp  ! Temporary arrays for extracting results
      integer :: gen, i, j, best_idx
      real(dblprec) :: trial_energy, energy_std, energy_range
      integer :: i_stat, i_all
      character(len=30) :: filn
      logical :: converged
      integer :: n_qvec_actual  ! Actual number of q-vectors (3 in skyrmion mode, de_n_qvec otherwise)
      
      ! Set number of parameters
      ! Skyrmion mode: 3 (q-vector) + 2 (s spherical) + 2 (n spherical) = 7
      !                (qz constrained to 0 if de_qxy_plane = T)
      ! Normal mode: 3*n_qvec (q-vectors) + 4 (s and n_vec in spherical coords)
      if (de_skyrmion_mode) then
         de_n_params = 3 + 2 + 2  ! q, s_spherical, n_spherical
         n_qvec_actual = 3        ! Generate 3 q-vectors for energy evaluation
      else
         de_n_params = 3 * de_n_qvec + 4
         n_qvec_actual = de_n_qvec
      end if
      
      if (de_verbose) then
         write(*,'(a)') repeat('=',80)
         write(*,'(a)') ' Differential Evolution Spin-Spiral Minimization'
         write(*,'(a)') repeat('=',80)
         write(*,'(a,i4)')    '  Population size:     ', de_pop_size
         write(*,'(a,i4)')    '  Max generations:     ', de_max_gen
         write(*,'(a,f6.3)')  '  Mutation factor (F): ', de_F
         write(*,'(a,f6.3)')  '  Crossover prob (CR): ', de_CR
         write(*,'(a,i4)')    '  Number of q-vectors: ', n_qvec_actual
         write(*,'(a,i4)')    '  Total parameters:    ', de_n_params
         if (de_skyrmion_mode) then
            write(*,'(a)') '  Mode: Skyrmion (3 q-vectors with 120 deg rotation)'
         end if
         if (de_qxy_plane) then
            write(*,'(a)') '  Q-vector constraint: x-y plane (qz = 0)'
         end if
         write(*,'(a)') repeat('=',80)
      end if
      
      ! Allocate arrays
      allocate(population(de_n_params, de_pop_size), stat=i_stat)
      call memocc(i_stat, product(shape(population))*kind(population), &
                  'population', 'de_minimize_qvectors')
      allocate(energies(de_pop_size), stat=i_stat)
      call memocc(i_stat, product(shape(energies))*kind(energies), &
                  'energies', 'de_minimize_qvectors')
      allocate(trial_params(de_n_params), stat=i_stat)
      call memocc(i_stat, product(shape(trial_params))*kind(trial_params), &
                  'trial_params', 'de_minimize_qvectors')
      allocate(mutant_params(de_n_params), stat=i_stat)
      call memocc(i_stat, product(shape(mutant_params))*kind(mutant_params), &
                  'mutant_params', 'de_minimize_qvectors')
      allocate(de_best_params(de_n_params), stat=i_stat)
      call memocc(i_stat, product(shape(de_best_params))*kind(de_best_params), &
                  'de_best_params', 'de_minimize_qvectors')
      
      ! Initialize population
      call de_initialize_population(population)
      
      ! Evaluate initial population
      do i = 1, de_pop_size
         energies(i) = de_evaluate_energy(population(:,i), Natom, Mensemble, NA, &
            coord, emomM, mmom, hfield, OPT_flag, max_no_constellations, &
            maxNoConstl, unitCellType, constlNCoup, constellations, &
            constellationsNeighType, Num_macro, cell_index, emomM_macro, &
            macro_nlistsize, N1, N2, N3)
      end do
      
      ! Find initial best
      best_idx = minloc(energies, dim=1)
      de_best_energy = energies(best_idx)
      de_best_params = population(:, best_idx)
      
      ! Open output file
      write(filn,'(''de_optimization.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")
      write(ofileno,'(a)') "#  Gen    Best_Energy(mRy)    Avg_Energy(mRy)    Energy_Std    "
      
      if (de_verbose) then
         write(*,'(a,f16.8)') '  Initial best energy: ', de_best_energy
      end if
      
      ! Main DE loop
      converged = .false.
      do gen = 1, de_max_gen
         
         ! Evolve population
         do i = 1, de_pop_size
            
            ! Mutation and crossover
            call de_mutate_and_crossover(population, i, mutant_params, trial_params)
            
            ! Evaluate trial solution
            trial_energy = de_evaluate_energy(trial_params, Natom, Mensemble, NA, &
               coord, emomM, mmom, hfield, OPT_flag, max_no_constellations, &
               maxNoConstl, unitCellType, constlNCoup, constellations, &
               constellationsNeighType, Num_macro, cell_index, emomM_macro, &
               macro_nlistsize, N1, N2, N3)
            
            ! Selection: greedy
            if (trial_energy < energies(i)) then
               population(:, i) = trial_params
               energies(i) = trial_energy
               
               ! Update global best
               if (trial_energy < de_best_energy) then
                  de_best_energy = trial_energy
                  de_best_params = trial_params
               end if
            end if
            
         end do
         
         ! Local refinement of best solution every de_local_freq generations
         if (de_use_local_search .and. mod(gen, de_local_freq) == 0) then
            call de_local_refinement(de_best_params, Natom, Mensemble, NA, coord, emomM, mmom, &
               hfield, OPT_flag, max_no_constellations, maxNoConstl, unitCellType, &
               constlNCoup, constellations, constellationsNeighType, Num_macro, &
               cell_index, emomM_macro, macro_nlistsize, N1, N2, N3, de_best_energy)
            
            ! Update population with refined best
            best_idx = minloc(energies, 1)
            population(:, best_idx) = de_best_params
            energies(best_idx) = de_best_energy
         end if
         
         ! Statistics
         energy_std = sqrt(sum((energies - sum(energies)/de_pop_size)**2) / de_pop_size)
         energy_range = maxval(energies) - minval(energies)
         
         ! Output progress
         if (mod(gen, 10) == 0 .or. gen == 1) then
            write(ofileno, '(i6,4es20.10)') gen, de_best_energy, &
               sum(energies)/de_pop_size, energy_std
            
            if (de_verbose) then
               write(*,'(a,i5,a,f16.8,a,f16.8,a,es12.4)') &
                  '  Gen ', gen, '  Best: ', de_best_energy, &
                  '  Avg: ', sum(energies)/de_pop_size, '  Std: ', energy_std
            end if
         end if
         
         ! Check convergence
         if (energy_std < de_tol .and. gen > 20) then
            converged = .true.
            if (de_verbose) then
               write(*,'(a,i5)') '  Converged at generation: ', gen
            end if
            exit
         end if
         
      end do
      
      close(ofileno)
      
      ! Extract results
      call de_params_to_vectors(de_best_params, q_temp, s_temp, n_vec_result)
      q_result = q_temp
      s_result = s_temp  ! Return all s-vectors (3 for skyrmion, 1 for normal)
      energy_result = de_best_energy
      
      ! IMPORTANT: Set up final spin configuration with best parameters
      ! This ensures emomM and emom are correctly initialized for the rest of the simulation
      energy_result = de_evaluate_energy(de_best_params, Natom, Mensemble, NA, &
         coord, emomM, mmom, hfield, OPT_flag, max_no_constellations, maxNoConstl, &
         unitCellType, constlNCoup, constellations, constellationsNeighType, &
         Num_macro, cell_index, emomM_macro, macro_nlistsize, N1, N2, N3)
      
      ! Print final results
      if (de_verbose) then
         write(*,'(a)') repeat('=',80)
         write(*,'(a)') ' Optimization Results'
         write(*,'(a)') repeat('=',80)
         write(*,'(a,f16.8,a)') '  Final best energy: ', de_best_energy, ' mRy'
         do i = 1, n_qvec_actual
            write(*,'(a,i2,a,3f12.6)') '  q-vector ', i, ': ', q_result(:,i)
         end do
         if (de_skyrmion_mode) then
            ! Print all 3 s-vectors for skyrmion mode
            do i = 1, 3
               write(*,'(a,i2,a,3f12.6)') '  s-vector ', i, ': ', s_result(:,i)
            end do
         else
            write(*,'(a,3f12.6)') '  s-vector:  ', s_result(:,1)
         end if
         write(*,'(a,3f12.6)') '  n-vector:  ', n_vec_result
         write(*,'(a,f12.6)') '  s·n check: ', dot_product(s_result(:,1), n_vec_result)
         write(*,'(a)') repeat('=',80)
      end if
      
      ! Deallocate
      i_all = -product(shape(population))*kind(population)
      deallocate(population, stat=i_stat)
      call memocc(i_stat, i_all, 'population', 'de_minimize_qvectors')
      
      i_all = -product(shape(energies))*kind(energies)
      deallocate(energies, stat=i_stat)
      call memocc(i_stat, i_all, 'energies', 'de_minimize_qvectors')
      
      i_all = -product(shape(trial_params))*kind(trial_params)
      deallocate(trial_params, stat=i_stat)
      call memocc(i_stat, i_all, 'trial_params', 'de_minimize_qvectors')
      
      i_all = -product(shape(mutant_params))*kind(mutant_params)
      deallocate(mutant_params, stat=i_stat)
      call memocc(i_stat, i_all, 'mutant_params', 'de_minimize_qvectors')
      
   end subroutine de_minimize_qvectors
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_initialize_population
   !> @brief Initialize population with random valid parameters
   !-----------------------------------------------------------------------------
   subroutine de_initialize_population(population)
      implicit none
      real(dblprec), dimension(:,:), intent(out) :: population
      integer :: i, j, n_qvec_params, n_qvec_optimize
      real(dblprec), dimension(3) :: rnd
      real(dblprec) :: theta, phi
      
      ! In skyrmion mode, only optimize first q-vector
      if (de_skyrmion_mode) then
         n_qvec_optimize = 1
         n_qvec_params = 3
      else
         n_qvec_optimize = de_n_qvec
         n_qvec_params = 3 * de_n_qvec
      end if
      
      do i = 1, de_pop_size
         ! Initialize q-vectors: uniform in unit ball
         do j = 1, n_qvec_optimize
            call rng_uniform(rnd, 3)
            ! Use rejection sampling for uniform distribution in ball
            do while (sum(rnd**2) > 1.0_dblprec)
               call rng_uniform(rnd, 3)
               rnd = 2.0_dblprec * rnd - 1.0_dblprec
            end do
            population(3*(j-1)+1:3*j, i) = rnd * 0.8_dblprec  ! Start with smaller q
            
            ! If x-y plane constraint is active, set qz = 0
            if (de_qxy_plane) then
               population(3*j, i) = 0.0_dblprec
            end if
         end do
         
         ! Initialize s-vector in spherical coordinates (theta_s, phi_s)
         call rng_uniform(rnd, 2)
         population(n_qvec_params+1, i) = rnd(1) * 3.14159265358979_dblprec  ! theta: [0, pi]
         population(n_qvec_params+2, i) = rnd(2) * 2.0_dblprec * 3.14159265358979_dblprec  ! phi: [0, 2pi]
         
         ! Initialize n_vec in spherical coordinates, ensuring it's not parallel to s
         call rng_uniform(rnd, 2)
         theta = rnd(1) * 3.14159265358979_dblprec
         phi = rnd(2) * 2.0_dblprec * 3.14159265358979_dblprec
         ! Add offset to avoid parallelism
         theta = theta + 0.5_dblprec
         if (theta > 3.14159265358979_dblprec) theta = theta - 3.14159265358979_dblprec
         population(n_qvec_params+3, i) = theta
         population(n_qvec_params+4, i) = phi
      end do
      
   end subroutine de_initialize_population
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_mutate_and_crossover
   !> @brief DE mutation and crossover operations
   !-----------------------------------------------------------------------------
   subroutine de_mutate_and_crossover(population, target_idx, mutant, trial)
      implicit none
      real(dblprec), dimension(:,:), intent(in) :: population
      integer, intent(in) :: target_idx
      real(dblprec), dimension(:), intent(out) :: mutant
      real(dblprec), dimension(:), intent(out) :: trial
      
      integer :: a, b, c, j, r_idx
      real(dblprec), dimension(1) :: rnd
      integer, dimension(3) :: indices
      logical :: valid
      
      ! Select three distinct random individuals (different from target)
      valid = .false.
      do while (.not. valid)
         call de_select_random_indices(de_pop_size, target_idx, indices)
         a = indices(1)
         b = indices(2)
         c = indices(3)
         valid = (a /= b) .and. (b /= c) .and. (a /= c)
      end do
      
      ! Mutation: mutant = a + F * (b - c)
      if (de_strategy == 1) then
         ! DE/rand/1
         mutant = population(:, a) + de_F * (population(:, b) - population(:, c))
      else if (de_strategy == 2) then
         ! DE/best/1
         mutant = de_best_params + de_F * (population(:, b) - population(:, c))
      end if
      
      ! Apply constraints to mutant
      call de_apply_constraints(mutant)
      
      ! Crossover: binomial
      call rng_uniform(rnd, 1)
      r_idx = int(rnd(1) * de_n_params) + 1  ! Random parameter index
      
      do j = 1, de_n_params
         call rng_uniform(rnd, 1)
         if (rnd(1) < de_CR .or. j == r_idx) then
            trial(j) = mutant(j)
         else
            trial(j) = population(j, target_idx)
         end if
      end do
      
      ! Apply constraints to trial
      call de_apply_constraints(trial)
      
   end subroutine de_mutate_and_crossover
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_select_random_indices
   !> @brief Select random distinct indices
   !-----------------------------------------------------------------------------
   subroutine de_select_random_indices(pop_size, exclude, indices)
      implicit none
      integer, intent(in) :: pop_size, exclude
      integer, dimension(3), intent(out) :: indices
      real(dblprec), dimension(1) :: rnd
      integer :: i
      
      do i = 1, 3
         call rng_uniform(rnd, 1)
         indices(i) = int(rnd(1) * pop_size) + 1
         do while (indices(i) == exclude .or. any(indices(1:i-1) == indices(i)))
            call rng_uniform(rnd, 1)
            indices(i) = int(rnd(1) * pop_size) + 1
         end do
      end do
      
   end subroutine de_select_random_indices
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_apply_constraints
   !> @brief Apply constraints to parameter vector
   !-----------------------------------------------------------------------------
   subroutine de_apply_constraints(params)
      implicit none
      real(dblprec), dimension(:), intent(inout) :: params
      integer :: i, n_qvec_params, n_qvec_optimize
      real(dblprec) :: qnorm
      real(dblprec), dimension(3) :: s_vec, n_vec, s_perp
      real(dblprec) :: pi, dot_sn
      
      pi = 4.0_dblprec * atan(1.0_dblprec)
      
      ! In skyrmion mode, only constrain the first q-vector
      if (de_skyrmion_mode) then
         n_qvec_optimize = 1
         n_qvec_params = 3
      else
         n_qvec_optimize = de_n_qvec
         n_qvec_params = 3 * de_n_qvec
      end if
      
      ! Constrain q-vectors to unit ball: |q| <= 1
      do i = 1, n_qvec_optimize
         ! If x-y plane constraint is active, enforce qz = 0
         if (de_qxy_plane) then
            params(3*i) = 0.0_dblprec
         end if
         
         qnorm = sqrt(sum(params(3*(i-1)+1:3*i)**2))
         if (qnorm > 1.0_dblprec) then
            params(3*(i-1)+1:3*i) = params(3*(i-1)+1:3*i) / qnorm
         end if
      end do
      
      ! Wrap spherical coordinates
      ! s-vector angles
      if (params(n_qvec_params+1) < 0.0_dblprec) params(n_qvec_params+1) = 0.0_dblprec
      if (params(n_qvec_params+1) > pi) params(n_qvec_params+1) = pi
      params(n_qvec_params+2) = modulo(params(n_qvec_params+2), 2.0_dblprec*pi)
      
      ! n_vec angles
      if (params(n_qvec_params+3) < 0.0_dblprec) params(n_qvec_params+3) = 0.0_dblprec
      if (params(n_qvec_params+3) > pi) params(n_qvec_params+3) = pi
      params(n_qvec_params+4) = modulo(params(n_qvec_params+4), 2.0_dblprec*pi)
      
      ! Optional: Ensure s and n_vec are not too parallel
      if (de_constraint_mode == 'P') then
         call de_spherical_to_cartesian(params(n_qvec_params+1), params(n_qvec_params+2), s_vec)
         call de_spherical_to_cartesian(params(n_qvec_params+3), params(n_qvec_params+4), n_vec)
         
         dot_sn = abs(dot_product(s_vec, n_vec))
         if (dot_sn > 0.95_dblprec) then
            ! Make n_vec perpendicular to s_vec using Gram-Schmidt
            n_vec = n_vec - dot_product(n_vec, s_vec) * s_vec
            n_vec = n_vec / sqrt(sum(n_vec**2))
            ! Convert back to spherical
            call de_cartesian_to_spherical(n_vec, params(n_qvec_params+3), params(n_qvec_params+4))
         end if
      end if
      
   end subroutine de_apply_constraints
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_spherical_to_cartesian
   !> @brief Convert spherical to Cartesian coordinates
   !-----------------------------------------------------------------------------
   subroutine de_spherical_to_cartesian(theta, phi, vec)
      implicit none
      real(dblprec), intent(in) :: theta, phi
      real(dblprec), dimension(3), intent(out) :: vec
      
      vec(1) = sin(theta) * cos(phi)
      vec(2) = sin(theta) * sin(phi)
      vec(3) = cos(theta)
      
   end subroutine de_spherical_to_cartesian
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_cartesian_to_spherical
   !> @brief Convert Cartesian to spherical coordinates
   !-----------------------------------------------------------------------------
   subroutine de_cartesian_to_spherical(vec, theta, phi)
      implicit none
      real(dblprec), dimension(3), intent(in) :: vec
      real(dblprec), intent(out) :: theta, phi
      real(dblprec) :: r
      
      r = sqrt(sum(vec**2))
      theta = acos(vec(3) / (r + 1.0e-12_dblprec))
      phi = atan2(vec(2), vec(1))
      if (phi < 0.0_dblprec) phi = phi + 2.0_dblprec * 4.0_dblprec * atan(1.0_dblprec)
      
   end subroutine de_cartesian_to_spherical
   
   !-----------------------------------------------------------------------------
   ! FUNCTION: de_evaluate_energy
   !> @brief Evaluate energy for given parameter vector
   !-----------------------------------------------------------------------------
   function de_evaluate_energy(params, Natom, Mensemble, NA, coord, emomM, mmom, &
      hfield, OPT_flag, max_no_constellations, maxNoConstl, unitCellType, &
      constlNCoup, constellations, constellationsNeighType, Num_macro, &
      cell_index, emomM_macro, macro_nlistsize, N1, N2, N3) result(energy)
      
      use Math_functions, only: f_wrap_coord_diff
      use MomentData, only: emom
      
      implicit none
      
      real(dblprec), dimension(:), intent(in) :: params
      integer, intent(in) :: Natom, Mensemble, NA, N1, N2, N3
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom
      real(dblprec), dimension(3), intent(in) :: hfield
      logical, intent(in) :: OPT_flag
      integer, intent(in) :: max_no_constellations
      integer, dimension(Mensemble), intent(in) :: maxNoConstl
      integer, dimension(Natom, Mensemble), intent(in) :: unitCellType
      real(dblprec), dimension(ham%max_no_neigh, max_no_constellations, Mensemble), intent(in) :: constlNCoup
      real(dblprec), dimension(3, max_no_constellations, Mensemble), intent(in) :: constellations
      integer, dimension(ham%max_no_neigh, max_no_constellations, Mensemble), intent(in) :: constellationsNeighType
      integer, intent(in) :: Num_macro
      integer, dimension(Natom), intent(in) :: cell_index
      real(dblprec), dimension(3, Num_macro, Mensemble), intent(in) :: emomM_macro
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize
      real(dblprec) :: energy
      
      ! Local variables
      real(dblprec), dimension(3, 3) :: q_vecs  ! Max 3 q-vectors (for skyrmion mode)
      real(dblprec), dimension(3, 3) :: s_vecs  ! Multiple s-vectors for skyrmion mode
      real(dblprec), dimension(3) :: n_vec, m_j, srvec
      integer :: ia, k, iq, countstart, I1, I2, I3, n_qvec_use
      real(dblprec) :: pi, qr
      
      pi = 4.0_dblprec * atan(1.0_dblprec)
      
      ! Determine actual number of q-vectors to use
      if (de_skyrmion_mode) then
         n_qvec_use = 3
      else
         n_qvec_use = de_n_qvec
      end if
      
      ! Convert parameters to vectors
      call de_params_to_vectors(params, q_vecs, s_vecs, n_vec)
      
      ! Set up external field
      k = 1  ! Use first ensemble
      do ia = 1, Natom
         external_field(1:3,ia,k) = hfield
         beff(1:3,ia,k) = 0.0_dblprec
         beff1(1:3,ia,k) = 0.0_dblprec
         beff2(1:3,ia,k) = 0.0_dblprec
      end do
      
      ! Reference atom (center of system)
      I1 = N1/2
      I2 = N2/2
      I3 = N3/2
      countstart = 0 + I1*NA + I2*N1*NA + I3*N2*N1*NA
      
      ! Set up spin configuration
      do ia = 1, Natom
         call f_wrap_coord_diff(Natom, coord, ia, countstart+1, srvec)
         
         m_j = 0.0_dblprec
         do iq = 1, n_qvec_use
            qr = q_vecs(1,iq)*srvec(1) + q_vecs(2,iq)*srvec(2) + q_vecs(3,iq)*srvec(3)
            qr = 2.0_dblprec * pi * qr
            ! Spin spiral: m = n*cos(q·r) + s*sin(q·r)
            ! Each q-vector uses its own s-vector (important for skyrmions)
            m_j = m_j + n_vec * cos(qr) + s_vecs(:, iq) * sin(qr)
         end do
         
         ! Normalize
         m_j = m_j / sqrt(sum(m_j**2) + 1.0e-12_dblprec)
         emomM(1:3,ia,k) = m_j * mmom(ia,k)
         emom(1:3,ia,k) = m_j
      end do
      
      ! Calculate energy
      energy = 0.0_dblprec
      call effective_field(Natom, Mensemble, 1, Natom, emomM, mmom, &
         external_field, time_external_field, beff, beff1, beff2, &
         OPT_flag, max_no_constellations, maxNoConstl, unitCellType, &
         constlNCoup, constellations, constellationsNeighType, &
         energy, Num_macro, cell_index, emomM_macro, macro_nlistsize, &
         NA, N1, N2, N3)
      
      energy = energy / Natom
      
   end function de_evaluate_energy
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_params_to_vectors
   !> @brief Convert parameter array to q, s, n_vec
   !> In skyrmion mode, generates 3 q-vectors with 120-degree rotation
   !-----------------------------------------------------------------------------
   subroutine de_params_to_vectors(params, q_vecs, s_vecs, n_vec)
      use Depondt, only: rodmat
      implicit none
      real(dblprec), dimension(:), intent(in) :: params
      real(dblprec), dimension(3, 3), intent(out) :: q_vecs  ! Fixed size for max q-vectors
      real(dblprec), dimension(3, 3), intent(out) :: s_vecs  ! Multiple s-vectors for skyrmion mode
      real(dblprec), dimension(3), intent(out) :: n_vec
      integer :: i, n_qvec_params, qq
      real(dblprec) :: pi, theta
      real(dblprec), dimension(3,3) :: R_mat
      
      pi = 4.0_dblprec * atan(1.0_dblprec)
      
      if (de_skyrmion_mode) then
         ! Skyrmion mode: optimize single q, s, n_vec and generate 3q by rotation
         ! Following sweep_q3 approach exactly
         n_qvec_params = 3  ! q-vector params
         
         ! Extract first q-vector
         q_vecs(:, 1) = params(1:3)
         
         ! Extract s-vector from spherical coordinates
         call de_spherical_to_cartesian(params(4), params(5), s_vecs(:, 1))
         
         ! Extract n_vec from spherical coordinates
         call de_spherical_to_cartesian(params(6), params(7), n_vec)
         
         ! Generate q2, q3 and s2, s3 by 120° and 240° rotations around z-axis
         ! This matches sweep_q3 exactly
         do qq = 2, 3
            theta = 2.0_dblprec * pi / 3.0_dblprec * real(qq-1, dblprec)
            call rodmat([0.0_dblprec, 0.0_dblprec, 1.0_dblprec], theta, R_mat)
            q_vecs(:, qq) = matmul(R_mat, q_vecs(:, 1))
            s_vecs(:, qq) = matmul(R_mat, s_vecs(:, 1))
         end do
         
      else
         ! Normal mode: independent q-vectors
         n_qvec_params = 3 * de_n_qvec
         
         ! Extract q-vectors
         do i = 1, de_n_qvec
            q_vecs(:, i) = params(3*(i-1)+1:3*i)
         end do
         
         ! Extract s-vector from spherical coordinates (same for all q in normal mode)
         call de_spherical_to_cartesian(params(n_qvec_params+1), params(n_qvec_params+2), s_vecs(:, 1))
         do i = 2, 3
            s_vecs(:, i) = s_vecs(:, 1)
         end do
         
         ! Extract n_vec from spherical coordinates
         call de_spherical_to_cartesian(params(n_qvec_params+3), params(n_qvec_params+4), n_vec)
      end if
      
   end subroutine de_params_to_vectors
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: read_parameters_de
   !> @brief Read DE parameters from input file
   !-----------------------------------------------------------------------------
   subroutine read_parameters_de(ifile)
      use FileParser
      implicit none
      
      integer, intent(in) :: ifile
      character(len=50) :: keyword
      integer :: rd_len, i_err, i_errb
      logical :: comment
      
      do
         10 continue
         keyword = ""
         call bytereader(keyword, rd_len, ifile, i_errb)
         call caps2small(keyword)
         
         comment = (scan(trim(keyword),'%')==1) .or. (scan(trim(keyword),'#')==1) .or. &
                   (scan(trim(keyword),'*')==1) .or. (scan(trim(keyword),'=')==1) .or. &
                   (scan(trim(keyword),'!')==1)
         
         if (comment) then
            read(ifile,*)
         else
            keyword = trim(keyword)
            select case(keyword)
            case('de_pop_size')
               read(ifile,*,iostat=i_err) de_pop_size
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_max_gen')
               read(ifile,*,iostat=i_err) de_max_gen
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_f')
               read(ifile,*,iostat=i_err) de_F
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_cr')
               read(ifile,*,iostat=i_err) de_CR
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_tol')
               read(ifile,*,iostat=i_err) de_tol
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_strategy')
               read(ifile,*,iostat=i_err) de_strategy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_n_qvec')
               read(ifile,*,iostat=i_err) de_n_qvec
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_constraint_mode')
               read(ifile,*,iostat=i_err) de_constraint_mode
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_qxy_plane')
               read(ifile,*,iostat=i_err) de_qxy_plane
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_skyrmion_mode')
               read(ifile,*,iostat=i_err) de_skyrmion_mode
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_verbose')
               read(ifile,*,iostat=i_err) de_verbose
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_use_local_search')
               read(ifile,*,iostat=i_err) de_use_local_search
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_local_steps')
               read(ifile,*,iostat=i_err) de_local_steps
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_local_freq')
               read(ifile,*,iostat=i_err) de_local_freq
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_local_tol')
               read(ifile,*,iostat=i_err) de_local_tol
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            case('de_local_step')
               read(ifile,*,iostat=i_err) de_local_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            end select
         end if
         
         if (i_errb==20) goto 20
         if (i_errb==10) goto 10
      end do
      
      20 continue
      rewind(ifile)
      
   end subroutine read_parameters_de
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_qminimizer_init
   !> @brief Initialize DE parameters with defaults
   !-----------------------------------------------------------------------------
   subroutine de_qminimizer_init()
      implicit none
      
      de_pop_size = 50
      de_max_gen = 200
      de_F = 0.8_dblprec
      de_CR = 0.9_dblprec
      de_tol = 1.0e-6_dblprec
      de_strategy = 1
      de_n_qvec = 1
      de_constraint_mode = 'P'
      de_qxy_plane = .false.
      de_skyrmion_mode = .false.
      de_verbose = .true.
      de_best_energy = 1.0e10_dblprec
      
   end subroutine de_qminimizer_init
   
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_setup_final_configuration
   !> @brief Set up final spin configuration with optimized parameters
   !-----------------------------------------------------------------------------
   subroutine de_setup_final_configuration(Natom, Mensemble, NA, coord, emomM, &
      mmom, q_vecs, s_vec, n_vec_in, N1, N2, N3)
      
      use Math_functions, only: f_wrap_coord_diff
      use MomentData, only: emom
      
      implicit none
      
      integer, intent(in) :: Natom, Mensemble, NA, N1, N2, N3
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom
      real(dblprec), dimension(3, de_n_qvec), intent(in) :: q_vecs
      real(dblprec), dimension(3), intent(in) :: s_vec, n_vec_in
      
      integer :: ia, k, iq, countstart, I1, I2, I3
      real(dblprec) :: pi, qr
      real(dblprec), dimension(3) :: m_j, srvec
      
      pi = 4.0_dblprec * atan(1.0_dblprec)
      
      ! Reference atom
      I1 = N1/2
      I2 = N2/2
      I3 = N3/2
      countstart = 0 + I1*NA + I2*N1*NA + I3*N2*N1*NA
      
      do k = 1, Mensemble
         do ia = 1, Natom
            call f_wrap_coord_diff(Natom, coord, ia, countstart+1, srvec)
            
            m_j = 0.0_dblprec
            do iq = 1, de_n_qvec
               qr = q_vecs(1,iq)*srvec(1) + q_vecs(2,iq)*srvec(2) + q_vecs(3,iq)*srvec(3)
               qr = 2.0_dblprec * pi * qr
               m_j = m_j + n_vec_in * cos(qr) + s_vec * sin(qr)
            end do
            
            ! Normalize
            m_j = m_j / sqrt(sum(m_j**2) + 1.0e-12_dblprec)
            emomM(1:3,ia,k) = m_j * mmom(ia,k)
            emom(1:3,ia,k) = m_j
         end do
      end do
      
   end subroutine de_setup_final_configuration

      !-----------------------------------------------------------------------------
   ! SUBROUTINE: de_local_refinement  
   !> @brief Gradient-based local refinement using numerical gradients
   !> @details Uses finite differences to compute gradients and performs
   !>          gradient descent with adaptive step size
   !-----------------------------------------------------------------------------
   subroutine de_local_refinement(params, Natom, Mensemble, NA, coord, emomM, mmom, &
      hfield, OPT_flag, max_no_constellations, maxNoConstl, unitCellType, &
      constlNCoup, constellations, constellationsNeighType, Num_macro, &
      cell_index, emomM_macro, macro_nlistsize, N1, N2, N3, energy)
      
      implicit none
      
      real(dblprec), dimension(:), intent(inout) :: params
      integer, intent(in) :: Natom, Mensemble, NA, N1, N2, N3
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom
      real(dblprec), dimension(3), intent(in) :: hfield
      logical, intent(in) :: OPT_flag
      integer, intent(in) :: max_no_constellations
      integer, dimension(Mensemble), intent(in) :: maxNoConstl
      integer, dimension(Natom, Mensemble), intent(in) :: unitCellType
      real(dblprec), dimension(ham%max_no_neigh, max_no_constellations, Mensemble), intent(in) :: constlNCoup
      real(dblprec), dimension(3, max_no_constellations, Mensemble), intent(in) :: constellations
      integer, dimension(ham%max_no_neigh, max_no_constellations, Mensemble), intent(in) :: constellationsNeighType
      integer, intent(in) :: Num_macro
      integer, dimension(Natom), intent(in) :: cell_index
      real(dblprec), dimension(3, Num_macro, Mensemble), intent(in) :: emomM_macro
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize
      real(dblprec), intent(inout) :: energy
      
      ! Local variables
      real(dblprec), dimension(:), allocatable :: gradient, params_new, params_old
      real(dblprec) :: energy_new, energy_old, step_size, grad_norm
      real(dblprec) :: improvement
      integer :: iter, i
      real(dblprec), parameter :: eps = 1.0e-6_dblprec
      
      allocate(gradient(de_n_params))
      allocate(params_new(de_n_params))
      allocate(params_old(de_n_params))
      
      params_old = params
      energy_old = energy
      step_size = de_local_step
      
      if (de_verbose) write(*,'(a)') '  Starting gradient-based local refinement...'
      
      do iter = 1, de_local_steps
         
         ! Compute numerical gradient using finite differences
         do i = 1, de_n_params
            params_new = params
            params_new(i) = params(i) + eps
            call de_apply_constraints(params_new)
            
            energy_new = de_evaluate_energy(params_new, Natom, Mensemble, NA, coord, emomM, mmom, &
               hfield, OPT_flag, max_no_constellations, maxNoConstl, unitCellType, &
               constlNCoup, constellations, constellationsNeighType, Num_macro, &
               cell_index, emomM_macro, macro_nlistsize, N1, N2, N3)
            
            gradient(i) = (energy_new - energy) / eps
         end do
         
         grad_norm = sqrt(sum(gradient**2))
         
         ! Gradient descent step
         params_new = params - step_size * gradient / (grad_norm + 1.0e-12_dblprec)
         call de_apply_constraints(params_new)
         
         energy_new = de_evaluate_energy(params_new, Natom, Mensemble, NA, coord, emomM, mmom, &
            hfield, OPT_flag, max_no_constellations, maxNoConstl, unitCellType, &
            constlNCoup, constellations, constellationsNeighType, Num_macro, &
            cell_index, emomM_macro, macro_nlistsize, N1, N2, N3)
         
         improvement = energy - energy_new
         
         ! Adaptive step size
         if (improvement > 0.0_dblprec) then
            params = params_new
            energy = energy_new
            step_size = min(step_size * 1.2_dblprec, 0.1_dblprec)
            
            if (mod(iter, 10) == 0 .and. de_verbose) then
               write(*,'(a,i4,a,f16.8,a,es12.4)') '    Iter ', iter, &
                  '  Energy: ', energy, '  |grad|: ', grad_norm
            end if
            
            ! Check convergence
            if (grad_norm < de_local_tol) then
               if (de_verbose) write(*,'(a,i4)') '    Converged at iteration: ', iter
               exit
            end if
         else
            step_size = step_size * 0.5_dblprec
            if (step_size < 1.0e-8_dblprec) then
               if (de_verbose) write(*,'(a)') '    Step size too small, stopping.'
               exit
            end if
         end if
         
      end do
      
      if (de_verbose) then
         write(*,'(a,f16.8,a)') '  Local refinement complete. Final energy: ', energy, ' mRy'
         write(*,'(a,f16.8,a)') '  Energy improvement: ', (energy_old - energy)*1000.0_dblprec, ' mRy'
      end if
      
      deallocate(gradient, params_new, params_old)
      
   end subroutine de_local_refinement

end module de_minimizer