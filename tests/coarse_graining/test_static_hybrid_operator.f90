program test_static_hybrid_operator

   use Parameters, only : dblprec
   use BlockTopology
   use stiffness, only : coarse_material_type, COARSE_MUB_SI
   use CoarseTensorOperator
   use SmoothProjectedOperator
   use StaticHybridOperator

   implicit none

   integer :: failures
   real(dblprec), parameter :: pi = acos(-1.0_dblprec)
   real(dblprec), parameter :: mub_si = COARSE_MUB_SI

   failures = 0
   call test_limiting_masks_and_ownership()
   call test_anisotropic_skew_ownership_non_overlap()
   call test_buffer_patch_and_derivatives()
   call test_texture_translation_and_refinement()

   if (failures /= 0) then
      write(*,'(a,i0)') 'static hybrid operator tests failed: ',failures
      stop 1
   end if
   write(*,'(a)') 'static hybrid operator tests passed'

contains

   subroutine test_limiting_masks_and_ownership()
      type(block_topology_type) :: topology
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      type(static_hybrid_operator_type) :: hybrid
      type(static_hybrid_energy_type) :: energy
      type(coarse_energy_terms_type) :: coarse_energy
      integer, parameter :: n = 12
      integer :: atom, bond, status
      character(len=512) :: message
      logical :: mask(n)
      integer :: bonds(2,n)
      real(dblprec) :: displacement(3,n), matrix(3,3,n), coordinate(3,n)
      real(dblprec) :: fine(3,n), coarse(3,1,n), fine_field(3,n)
      real(dblprec) :: coarse_field(3,1,n), baseline_field(3,n)
      real(dblprec) :: direct_field(3,n), external(3,n)
      real(dblprec) :: onsite_energy(n), onsite_field(3,n), baseline_energy

      call make_chain_fixture(n,1,topology,tensor,projection,bonds,displacement, &
         matrix,coordinate,status,message)
      call check(status == STATIC_HYBRID_OK,'limiting-mask fixture builds: '//trim(message))
      do atom = 1,n
         fine(:,atom) = (/cos(0.31_dblprec*atom),sin(0.31_dblprec*atom), &
            0.2_dblprec*cos(0.17_dblprec*atom)/)
         fine(:,atom) = fine(:,atom)/norm3(fine(:,atom))
         coarse(:,1,atom) = fine(:,atom)
         external(:,atom) = (/0.02_dblprec,-0.01_dblprec,0.03_dblprec/) * &
            (1.0_dblprec+0.01_dblprec*atom)
         onsite_energy(atom) = -mub_si*1.7_dblprec* &
            dot_product(external(:,atom),fine(:,atom))
         onsite_field(:,atom) = external(:,atom)
      end do

      mask = .true.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,1,status,message)
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,fine_field, &
         coarse_field,energy,status,message,onsite_energy,onsite_field,external)
      call atomistic_baseline(fine,bonds,matrix,external,baseline_energy,baseline_field)
      call check(status == STATIC_HYBRID_OK,'all-fine dispatch evaluates: '//trim(message))
      call check_close(energy%total_j,baseline_energy,2.0d-14, &
         'all-fine dispatch equals baseline atomistic energy')
      call check(maxval(abs(fine_field-baseline_field)) < 2.0d-13* &
         max(1.0_dblprec,maxval(abs(baseline_field))), &
         'all-fine dispatch equals baseline atomistic field')
      call check(maxval(abs(coarse_field)) == 0.0_dblprec .and. &
         count(hybrid%atomistic_bond_owner) == n, &
         'all-fine dispatch has no coarse owner and owns every atomistic bond once')

      mask = .false.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,1,status,message)
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,fine_field, &
         coarse_field,energy,status,message,onsite_energy,onsite_field,external)
      call evaluate_coarse_tensor_operator(tensor,coarse(:,1,:),direct_field, &
         coarse_energy,status,message,external_field_t=external)
      call check(status == COARSE_TENSOR_OK,'direct all-coarse reference evaluates')
      call check_close(energy%total_j,coarse_energy%total_j,2.0d-14, &
         'all-coarse hybrid dispatch equals accepted tensor energy')
      call check(maxval(abs(coarse_field(:,1,:)-direct_field)) < 2.0d-13* &
         max(1.0_dblprec,maxval(abs(direct_field))) .and. &
         maxval(abs(fine_field)) == 0.0_dblprec, &
         'all-coarse hybrid dispatch equals accepted tensor field')
      call check(count(hybrid%atomistic_bond_owner) == 0, &
         'all-coarse dispatch owns no atomistic bond')

      mask = .false.
      mask(5) = .true.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,0,status,message)
      do bond = 1,n
         atom = bonds(1,bond)
         call check(hybrid%atomistic_bond_owner(bond) .neqv. &
            (hybrid%coarse_block(atom) .and. hybrid%coarse_block(bonds(2,bond))), &
            'each mixed exchange/DMI bond has exactly one energy owner')
      end do
   end subroutine test_limiting_masks_and_ownership

   !> RCG-05F: re-runs test_limiting_masks_and_ownership's own three
   !> ownership-non-overlap checks (all-fine/all-coarse/mixed) on a
   !> genuinely anisotropic-block-width (width=3, so buffer_width_blocks
   !> differs between the chain axis and the other two -- see
   !> statichybridoperator.f90:180-187) AND non-orthogonal (skew) cell
   !> fixture, rather than the width=1 cubic one
   !> test_limiting_masks_and_ownership itself uses. Confirms short-range
   !> (atomistic_bond_owner) and on-site (coarse_block/onsite_owner)
   !> ownership -- coarsetensoroperator.f90:478-479,499-500 -- remain
   !> non-overlapping (never both, never neither) under that geometry, not
   !> only the cubic ones this invariant was presumably last checked
   !> against.
   subroutine test_anisotropic_skew_ownership_non_overlap()
      type(block_topology_type) :: topology
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      type(static_hybrid_operator_type) :: hybrid
      type(static_hybrid_energy_type) :: energy
      type(coarse_energy_terms_type) :: coarse_energy
      ! make_chain_fixture's first argument is an ATOM count (like
      ! production's ncell), not a block count: nblocks = n/width, exactly
      ! like ncell/block_size_x in production. test_limiting_masks_and_
      ! ownership's own n=12,width=1 fixture never exposed this distinction
      ! (atoms==blocks there); width=3 here makes it load-bearing.
      integer, parameter :: n = 30, width = 3, nblocks = n/width
      integer :: atom, block, bond, status
      character(len=512) :: message
      logical :: mask(nblocks)
      integer :: bonds(2,n)
      real(dblprec) :: displacement(3,n), matrix(3,3,n), coordinate(3,n)
      real(dblprec) :: fine(3,n), coarse(3,1,nblocks), fine_field(3,n)
      real(dblprec) :: coarse_field(3,1,nblocks), baseline_field(3,n)
      real(dblprec) :: direct_field(3,nblocks)
      real(dblprec) :: external_atom(3,n), external_block(3,nblocks)
      real(dblprec) :: onsite_energy(n), onsite_field(3,n), baseline_energy
      real(dblprec), parameter :: a = 2.0d-10
      real(dblprec) :: skew(3,3)

      ! Same isolate-one-variable-at-a-time skew as
      ! tests/coarse_graining/e2e/ownership_dipole_skew/README.md and
      ! test_coarse_tensor_operator.f90's test_dipole_unmasked_and_exactly_
      ! once: row 1 stays (a,0,0) so this chain's own bond geometry is
      ! unaffected; rows 2/3 carry the non-orthogonality.
      skew(:,1) = (/a,0.0_dblprec,0.0_dblprec/)
      skew(:,2) = (/0.35_dblprec*a,0.9_dblprec*a,0.0_dblprec/)
      skew(:,3) = (/0.2_dblprec*a,0.1_dblprec*a,1.15_dblprec*a/)

      call make_chain_fixture(n,width,topology,tensor,projection,bonds,displacement, &
         matrix,coordinate,status,message,cell_override=skew)
      call check(status == STATIC_HYBRID_OK, &
         'anisotropic/skew ownership fixture builds: '//trim(message))
      if (status /= STATIC_HYBRID_OK) return
      call check(topology%n_spatial_blocks == nblocks, &
         'anisotropic/skew fixture has the expected atom/block split (n/width)')
      do atom = 1,n
         fine(:,atom) = (/cos(0.23_dblprec*atom),sin(0.23_dblprec*atom), &
            0.2_dblprec*cos(0.13_dblprec*atom)/)
         fine(:,atom) = fine(:,atom)/norm3(fine(:,atom))
         external_atom(:,atom) = (/0.02_dblprec,-0.01_dblprec,0.03_dblprec/) * &
            (1.0_dblprec+0.01_dblprec*atom)
         onsite_energy(atom) = -mub_si*1.7_dblprec* &
            dot_product(external_atom(:,atom),fine(:,atom))
         onsite_field(:,atom) = external_atom(:,atom)
      end do
      do block = 1,nblocks
         coarse(:,1,block) = (/cos(0.19_dblprec*block),sin(0.19_dblprec*block), &
            0.2_dblprec*cos(0.11_dblprec*block)/)
         coarse(:,1,block) = coarse(:,1,block)/norm3(coarse(:,1,block))
         external_block(:,block) = (/0.015_dblprec,-0.008_dblprec,0.02_dblprec/) * &
            (1.0_dblprec+0.01_dblprec*block)
      end do

      mask = .true.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,1,status,message)
      call check(status == STATIC_HYBRID_OK, &
         'anisotropic/skew all-fine setup: '//trim(message))
      call check(hybrid%buffer_width_blocks(1) /= hybrid%buffer_width_blocks(2) .or. &
         hybrid%buffer_width_blocks(1) /= hybrid%buffer_width_blocks(3), &
         'fixture is genuinely anisotropic: buffer_width_blocks is not the same on every axis')
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,fine_field, &
         coarse_field,energy,status,message,onsite_energy,onsite_field,external_block)
      call atomistic_baseline(fine,bonds,matrix,external_atom,baseline_energy,baseline_field)
      call check(status == STATIC_HYBRID_OK, &
         'anisotropic/skew all-fine dispatch evaluates: '//trim(message))
      call check_close(energy%total_j,baseline_energy,2.0d-14, &
         'anisotropic/skew all-fine dispatch equals baseline atomistic energy')
      call check(maxval(abs(fine_field-baseline_field)) < 2.0d-13* &
         max(1.0_dblprec,maxval(abs(baseline_field))), &
         'anisotropic/skew all-fine dispatch equals baseline atomistic field')
      call check(maxval(abs(coarse_field)) == 0.0_dblprec .and. &
         count(hybrid%atomistic_bond_owner) == n, &
         'anisotropic/skew all-fine dispatch has no coarse owner and owns every atomistic bond once')

      mask = .false.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,1,status,message)
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,fine_field, &
         coarse_field,energy,status,message,onsite_energy,onsite_field,external_block)
      call evaluate_coarse_tensor_operator(tensor,coarse(:,1,:),direct_field, &
         coarse_energy,status,message,external_field_t=external_block)
      call check(status == COARSE_TENSOR_OK,'anisotropic/skew direct all-coarse reference evaluates')
      call check_close(energy%total_j,coarse_energy%total_j,2.0d-14, &
         'anisotropic/skew all-coarse hybrid dispatch equals accepted tensor energy')
      call check(maxval(abs(coarse_field(:,1,:)-direct_field)) < 2.0d-13* &
         max(1.0_dblprec,maxval(abs(direct_field))) .and. &
         maxval(abs(fine_field)) == 0.0_dblprec, &
         'anisotropic/skew all-coarse hybrid dispatch equals accepted tensor field')
      call check(count(hybrid%atomistic_bond_owner) == 0, &
         'anisotropic/skew all-coarse dispatch owns no atomistic bond')

      mask = .false.
      mask(3) = .true.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,0,status,message)
      call check(status == STATIC_HYBRID_OK, &
         'anisotropic/skew single-seed fixture builds: '//trim(message))
      do bond = 1,n
         call check(hybrid%atomistic_bond_owner(bond) .neqv. &
            ((.not. hybrid%atomistic_atom(bonds(1,bond))) .and. &
             (.not. hybrid%atomistic_atom(bonds(2,bond)))), &
            'anisotropic/skew: each mixed exchange/DMI bond has exactly one energy owner')
      end do
      do block = 1,nblocks
         call check(hybrid%atomistic_block(block) .neqv. hybrid%coarse_block(block), &
            'anisotropic/skew: every block is atomistic xor coarse, never both or neither')
      end do
      do atom = 1,n
         call check(hybrid%atomistic_atom(atom) .eqv. &
            hybrid%atomistic_block(topology%atom_to_block(atom)), &
            'anisotropic/skew: per-atom ownership matches its owning block exactly')
      end do
   end subroutine test_anisotropic_skew_ownership_non_overlap

   subroutine test_buffer_patch_and_derivatives()
      type(block_topology_type) :: topology
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      type(static_hybrid_operator_type) :: hybrid
      type(static_hybrid_energy_type) :: energy, uniform_energy, opposite_chirality_energy
      integer, parameter :: n = 32
      integer :: atom, block, status
      character(len=512) :: message
      logical :: mask(n)
      integer :: bonds(2,n)
      real(dblprec), parameter :: q = 2.0_dblprec*pi/real(n,dblprec)
      real(dblprec) :: displacement(3,n), matrix(3,3,n), coordinate(3,n)
      real(dblprec) :: fine(3,n), coarse(3,1,n), fine_field(3,n)
      real(dblprec) :: coarse_field(3,1,n), tangent(3)
      real(dblprec) :: baseline_field(3,n), baseline_energy, baseline_uniform
      real(dblprec) :: numerical, analytic

      call make_chain_fixture(n,1,topology,tensor,projection,bonds,displacement, &
         matrix,coordinate,status,message,dmi_bond=0.17d-21)
      mask = .false.
      mask(8) = .true.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,1,status,message)
      call check(status == STATIC_HYBRID_OK,'mixed buffer setup: '//trim(message))
      call check(all(hybrid%buffer_width_blocks == (/2,2,2/)) .and. &
         hybrid%maximum_interaction_radius_m > 0.0_dblprec, &
         'buffer width is interaction-radius ceiling plus safety dilation')
      call check(count(hybrid%block_state == STATIC_HYBRID_FINE) == 1 .and. &
         count(hybrid%block_state == STATIC_HYBRID_BUFFER) == 4, &
         'periodic dilation produces distinct static fine and buffer ownership')
      do atom = 1,n
         if (hybrid%atomistic_bond_owner(atom)) then
            call check(hybrid%atomistic_atom(bonds(1,atom)) .or. &
               hybrid%atomistic_atom(bonds(2,atom)), &
               'every owned cross-interface interaction touches the atomistic buffer')
         end if
      end do

      fine = 0.0_dblprec
      coarse = 0.0_dblprec
      fine(1,:) = 1.0_dblprec
      coarse(1,1,:) = 1.0_dblprec
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,fine_field, &
         coarse_field,uniform_energy,status,message)
      call check(status == STATIC_HYBRID_OK,'uniform mixed patch evaluates')
      call check(max_tangent_atoms(fine,fine_field,hybrid%atomistic_atom) < 2.0d-13 .and. &
         max_tangent_blocks(coarse,coarse_field,hybrid%coarse_block) < 2.0d-13, &
         'uniform interface patch has no spurious torque')

      do atom = 1,n
         fine(:,atom) = (/cos(q*real(atom-1,dblprec)), &
            sin(q*real(atom-1,dblprec)),0.0_dblprec/)
         coarse(:,1,atom) = fine(:,atom)
      end do
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,fine_field, &
         coarse_field,energy,status,message)
      call atomistic_baseline(fine,bonds,matrix,baseline_energy=baseline_energy, &
         baseline_field=baseline_field)
      fine = 0.0_dblprec
      fine(1,:) = 1.0_dblprec
      call atomistic_baseline(fine,bonds,matrix,baseline_energy=baseline_uniform, &
         baseline_field=baseline_field)
      call check_close(energy%total_j-uniform_energy%total_j, &
         baseline_energy-baseline_uniform,3.0d-14, &
         'constant long-wave spiral passes the block-one interface patch test')

      ! DMI-HYBRID-CROSSING operator fixture: the fine/buffer bonds and the
      ! coarse interior use the same directed DMI energy, so +D_zx must prefer
      ! the negative-q member of this otherwise exchange-degenerate pair.
      do atom = 1,n
         fine(:,atom) = (/cos(-q*real(atom-1,dblprec)), &
            sin(-q*real(atom-1,dblprec)),0.0_dblprec/)
         coarse(:,1,atom) = fine(:,atom)
      end do
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,fine_field, &
         coarse_field,opposite_chirality_energy,status,message)
      call check(energy%total_j > opposite_chirality_energy%total_j, &
         'mixed-resolution interface preserves the DMI-preferred negative-q chirality')

      do atom = 1,n
         fine(:,atom) = (/cos(q*real(atom-1,dblprec)+ &
            0.07_dblprec*sin(2.0_dblprec*q*real(atom-1,dblprec))), &
            sin(q*real(atom-1,dblprec)+ &
            0.07_dblprec*sin(2.0_dblprec*q*real(atom-1,dblprec))),0.0_dblprec/)
         coarse(:,1,atom) = fine(:,atom)
      end do
      atom = 8
      tangent = (/-fine(2,atom),fine(1,atom),0.0_dblprec/)
      call hybrid_finite_difference(hybrid,fine,coarse,matrix,.true.,atom,tangent,numerical)
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,fine_field, &
         coarse_field,energy,status,message)
      analytic = -mub_si*1.7_dblprec*dot_product(fine_field(:,atom),tangent)
      call check_close(numerical,analytic,1.0d-6, &
         'fine-side mixed interface field is the atomistic energy derivative')

      block = 11
      call check(hybrid%coarse_block(block),'selected derivative block is coarse-side ghost owner')
      tangent = (/-coarse(2,1,block),coarse(1,1,block),0.0_dblprec/)
      call hybrid_finite_difference(hybrid,fine,coarse,matrix,.false.,block,tangent,numerical)
      analytic = -mub_si*tensor%block_moment_mub* &
         dot_product(coarse_field(:,1,block),tangent)
      call check_close(numerical,analytic,5.0d-7, &
         'coarse ghost reaction uses the normalized prolongation adjoint')
   end subroutine test_buffer_patch_and_derivatives

   subroutine test_texture_translation_and_refinement()
      real(dblprec) :: wall_error(2), skyrmion_error(2), refinement_error(2)

      call translated_wall_error(19.5_dblprec,wall_error(1))
      call translated_wall_error(25.5_dblprec,wall_error(2))
      call check(maxval(wall_error) < 2.0d-13, &
         'domain-wall pair crosses a fixed interface without measurable block-one pinning')

      call translated_skyrmion_error(7.5_dblprec,skyrmion_error(1))
      call translated_skyrmion_error(10.5_dblprec,skyrmion_error(2))
      call check(maxval(skyrmion_error) < 3.0d-13, &
         'translated skyrmion remains energy-stable across the fixed interface')

      call spiral_refinement_error(4,refinement_error(1))
      call spiral_refinement_error(2,refinement_error(2))
      call check(refinement_error(2) < 0.35_dblprec*refinement_error(1), &
         'static interface error decreases under block refinement')
   end subroutine test_texture_translation_and_refinement

   subroutine translated_wall_error(center,error)
      real(dblprec), intent(in) :: center
      real(dblprec), intent(out) :: error
      type(block_topology_type) :: topology
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      type(static_hybrid_operator_type) :: hybrid
      type(static_hybrid_energy_type) :: energy, uniform_energy
      integer, parameter :: n = 64
      integer :: atom, status
      character(len=512) :: message
      logical :: mask(n)
      integer :: bonds(2,n)
      real(dblprec) :: displacement(3,n), matrix(3,3,n), coordinate(3,n)
      real(dblprec) :: fine(3,n), coarse(3,1,n), ff(3,n), cf(3,1,n)
      real(dblprec) :: baseline, baseline0, unused(3,n), phase, profile

      call make_chain_fixture(n,1,topology,tensor,projection,bonds,displacement, &
         matrix,coordinate,status,message)
      mask = .false.
      mask(16:22) = .true.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,1,status,message)
      do atom = 1,n
         phase = 2.0_dblprec*pi*(real(atom-1,dblprec)-center)/real(n,dblprec)
         profile = tanh(sin(phase)/0.22_dblprec)
         fine(:,atom) = (/sqrt(max(0.0_dblprec,1.0_dblprec-profile*profile)), &
            0.0_dblprec,profile/)
         coarse(:,1,atom) = fine(:,atom)
      end do
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,ff,cf,energy,status,message)
      call atomistic_baseline(fine,bonds,matrix,baseline_energy=baseline,baseline_field=unused)
      fine = 0.0_dblprec
      coarse = 0.0_dblprec
      fine(3,:) = 1.0_dblprec
      coarse(3,1,:) = 1.0_dblprec
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,ff,cf,uniform_energy,status,message)
      call atomistic_baseline(fine,bonds,matrix,baseline_energy=baseline0,baseline_field=unused)
      error = relative_error(energy%total_j-uniform_energy%total_j,baseline-baseline0)
   end subroutine translated_wall_error

   subroutine translated_skyrmion_error(center_x,error)
      real(dblprec), intent(in) :: center_x
      real(dblprec), intent(out) :: error
      type(block_topology_type) :: topology
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      type(static_hybrid_operator_type) :: hybrid
      type(static_hybrid_energy_type) :: energy, uniform_energy
      integer, parameter :: nx = 20, ny = 16, natom = nx*ny, nbond = 2*natom
      integer :: atom, x, y, status
      character(len=512) :: message
      logical :: mask(natom)
      integer :: bonds(2,nbond)
      real(dblprec) :: displacement(3,nbond), matrix(3,3,nbond), coordinate(3,natom)
      real(dblprec) :: fine(3,natom), coarse(3,1,natom), ff(3,natom), cf(3,1,natom)
      real(dblprec) :: baseline, baseline0, unused(3,natom), dx, dy, radius, theta, phi

      call make_square_fixture(nx,ny,topology,tensor,projection,bonds,displacement, &
         matrix,coordinate,status,message)
      mask = .false.
      do atom = 1,natom
         x = modulo(atom-1,nx)
         if (x >= 5 .and. x <= 8) mask(atom) = .true.
      end do
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,1,status,message)
      do atom = 1,natom
         x = modulo(atom-1,nx)
         y = (atom-1)/nx
         dx = real(x,dblprec)-center_x
         dy = real(y,dblprec)-7.5_dblprec
         radius = sqrt(dx*dx+dy*dy)
         theta = pi*max(0.0_dblprec,1.0_dblprec-radius/4.5_dblprec)
         phi = atan2(dy,dx)+0.5_dblprec*pi
         fine(:,atom) = (/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)
         coarse(:,1,atom) = fine(:,atom)
      end do
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,ff,cf,energy,status,message)
      call atomistic_baseline(fine,bonds,matrix,baseline_energy=baseline,baseline_field=unused)
      fine = 0.0_dblprec
      coarse = 0.0_dblprec
      fine(3,:) = 1.0_dblprec
      coarse(3,1,:) = 1.0_dblprec
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,ff,cf,uniform_energy,status,message)
      call atomistic_baseline(fine,bonds,matrix,baseline_energy=baseline0,baseline_field=unused)
      error = relative_error(energy%total_j-uniform_energy%total_j,baseline-baseline0)
   end subroutine translated_skyrmion_error

   subroutine spiral_refinement_error(width,error)
      integer, intent(in) :: width
      real(dblprec), intent(out) :: error
      type(block_topology_type) :: topology
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      type(static_hybrid_operator_type) :: hybrid
      type(static_hybrid_energy_type) :: energy, uniform_energy
      integer, parameter :: n = 64
      integer :: atom, block, status
      character(len=512) :: message
      logical, allocatable :: mask(:)
      integer :: bonds(2,n)
      real(dblprec) :: displacement(3,n), matrix(3,3,n), coordinate(3,n)
      real(dblprec) :: fine(3,n), ff(3,n), unused(3,n)
      real(dblprec), allocatable :: coarse(:,:,:), cf(:,:,:)
      real(dblprec) :: baseline, baseline0, q

      call make_chain_fixture(n,width,topology,tensor,projection,bonds,displacement, &
         matrix,coordinate,status,message)
      allocate(mask(topology%n_spatial_blocks),coarse(3,1,topology%n_spatial_blocks), &
         cf(3,1,topology%n_spatial_blocks))
      mask = .false.
      mask(4:min(6,size(mask))) = .true.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,1,status,message)
      q = 2.0_dblprec*pi/real(n,dblprec)
      do atom = 1,n
         fine(:,atom) = (/cos(q*(coordinate(1,atom)+0.5_dblprec)*real(width,dblprec)), &
            sin(q*(coordinate(1,atom)+0.5_dblprec)*real(width,dblprec)),0.0_dblprec/)
      end do
      do block = 1,topology%n_spatial_blocks
         coarse(:,1,block) = (/cos(q*real(width,dblprec)*(real(block,dblprec)-0.5_dblprec)), &
            sin(q*real(width,dblprec)*(real(block,dblprec)-0.5_dblprec)),0.0_dblprec/)
      end do
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,ff,cf,energy,status,message)
      call atomistic_baseline(fine,bonds,matrix,baseline_energy=baseline,baseline_field=unused)
      fine = 0.0_dblprec
      coarse = 0.0_dblprec
      fine(1,:) = 1.0_dblprec
      coarse(1,1,:) = 1.0_dblprec
      call evaluate_static_hybrid_operator(hybrid,fine,coarse,matrix,ff,cf,uniform_energy,status,message)
      call atomistic_baseline(fine,bonds,matrix,baseline_energy=baseline0,baseline_field=unused)
      error = relative_error(energy%total_j-uniform_energy%total_j,baseline-baseline0)
   end subroutine spiral_refinement_error

   subroutine make_chain_fixture(ncell,width,topology,tensor,projection,bonds, &
         displacement,matrix,coordinate,status,message,dmi_bond,cell_override)
      integer, intent(in) :: ncell, width
      type(block_topology_type), intent(out) :: topology
      type(coarse_tensor_operator_type), intent(out) :: tensor
      type(smooth_projected_operator_type), intent(out) :: projection
      integer, intent(out) :: bonds(2,ncell)
      real(dblprec), intent(out) :: displacement(3,ncell), matrix(3,3,ncell)
      real(dblprec), intent(out) :: coordinate(3,ncell)
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      real(dblprec), intent(in), optional :: dmi_bond
      ! RCG-05F: an optional non-cubic (and, if row 1 is left as (a,0,0),
      ! genuinely skew) atom cell override, so callers can exercise
      ! rebuild_static_hybrid_ownership's inverse_block_m1 metric
      ! (statichybridoperator.f90:180-187) on a non-orthogonal cell without
      ! a second fixture builder. Row 1 must stay (a,0,0) for the chain's
      ! own bond `displacement`/`coordinate` construction below to remain
      ! physically correct; only rows 2/3 (never used by this 1D chain's
      ! own geometry) are meaningful to skew.
      real(dblprec), intent(in), optional :: cell_override(3,3)
      type(coarse_material_type) :: material
      type(coarse_operator_options_type) :: options
      integer :: atom
      real(dblprec), parameter :: a = 2.0d-10, pair_j = 3.0d-21
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3), dbond, moment(ncell)

      dbond = 0.0_dblprec
      if (present(dmi_bond)) dbond = dmi_bond
      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      if (present(cell_override)) cell = cell_override
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,1,(/ncell,1,1/), &
         ncell,(/width,1,1/),cell,(/1/),status,message)
      if (status /= BLOCK_TOPOLOGY_OK) then
         status = STATIC_HYBRID_INVALID_TOPOLOGY
         return
      end if
      exchange = 0.0_dblprec
      exchange(1,1) = pair_j/(2.0_dblprec*a)
      dmi = 0.0_dblprec
      dmi(3,1) = dbond/(a*a)
      call make_material(cell,exchange,dmi,material)
      call setup_coarse_tensor_operator(tensor,topology,material,options,status,message)
      if (status /= COARSE_TENSOR_OK) return
      coordinate = 0.0_dblprec
      moment = 1.7_dblprec
      do atom = 1,ncell
         coordinate(1,atom) = (real(atom,dblprec)-0.5_dblprec)/real(width,dblprec)-0.5_dblprec
         bonds(:,atom) = (/atom,modulo(atom,ncell)+1/)
         displacement(:,atom) = (/a,0.0_dblprec,0.0_dblprec/)
         matrix(:,:,atom) = 0.0_dblprec
         matrix(1,1,atom) = pair_j
         matrix(2,2,atom) = pair_j
         matrix(3,3,atom) = pair_j
         matrix(1,2,atom) = -dbond
         matrix(2,1,atom) = dbond
      end do
      call setup_smooth_projected_operator(projection,topology,coordinate,moment,status,message)
      if (status == SMOOTH_PROJECTED_OK) status = STATIC_HYBRID_OK
   end subroutine make_chain_fixture

   subroutine make_square_fixture(nx,ny,topology,tensor,projection,bonds, &
         displacement,matrix,coordinate,status,message)
      integer, intent(in) :: nx, ny
      type(block_topology_type), intent(out) :: topology
      type(coarse_tensor_operator_type), intent(out) :: tensor
      type(smooth_projected_operator_type), intent(out) :: projection
      integer, intent(out) :: bonds(:,:)
      real(dblprec), intent(out) :: displacement(:,:), matrix(:,:,:), coordinate(:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      type(coarse_material_type) :: material
      type(coarse_operator_options_type) :: options
      integer :: atom, x, y, right, up, bond
      real(dblprec), parameter :: a = 2.0d-10, pair_j = 3.0d-21
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3), moment(nx*ny)

      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,1,(/nx,ny,1/), &
         nx*ny,(/1,1,1/),cell,(/1/),status,message)
      exchange = 0.0_dblprec
      exchange(1,1) = pair_j/(2.0_dblprec*a)
      exchange(2,2) = exchange(1,1)
      dmi = 0.0_dblprec
      call make_material(cell,exchange,dmi,material)
      call setup_coarse_tensor_operator(tensor,topology,material,options,status,message)
      coordinate = 0.0_dblprec
      moment = 1.7_dblprec
      bond = 0
      do y = 0,ny-1
         do x = 0,nx-1
            atom = 1+x+nx*y
            coordinate(:,atom) = (/real(x,dblprec),real(y,dblprec),0.0_dblprec/)
            right = 1+modulo(x+1,nx)+nx*y
            up = 1+x+nx*modulo(y+1,ny)
            bond = bond+1
            bonds(:,bond) = (/atom,right/)
            displacement(:,bond) = (/a,0.0_dblprec,0.0_dblprec/)
            matrix(:,:,bond) = 0.0_dblprec
            matrix(1,1,bond) = pair_j
            matrix(2,2,bond) = pair_j
            matrix(3,3,bond) = pair_j
            bond = bond+1
            bonds(:,bond) = (/atom,up/)
            displacement(:,bond) = (/0.0_dblprec,a,0.0_dblprec/)
            matrix(:,:,bond) = matrix(:,:,bond-1)
         end do
      end do
      call setup_smooth_projected_operator(projection,topology,coordinate,moment,status,message)
      if (status == SMOOTH_PROJECTED_OK) status = STATIC_HYBRID_OK
   end subroutine make_square_fixture

   subroutine make_material(cell,exchange,dmi,material)
      real(dblprec), intent(in) :: cell(3,3), exchange(3,3), dmi(3,3)
      type(coarse_material_type), intent(out) :: material

      material%ready = .true.
      material%n_basis = 1
      material%n_channels = 1
      material%cell_volume_m3 = cell(1,1)*cell(2,2)*cell(3,3)
      allocate(material%basis_to_channel(1),material%basis_moment_mub(1), &
         material%channel_moment_mub(1),material%channel_gamma(1), &
         material%channel_damping(1),material%exchange_stiffness(3,3,1,1), &
         material%spiralization(3,3,1,1))
      material%basis_to_channel = 1
      material%basis_moment_mub = 1.7_dblprec
      material%channel_moment_mub = 1.7_dblprec
      material%channel_gamma = 1.76d11
      material%channel_damping = 0.01_dblprec
      material%exchange_stiffness(:,:,1,1) = exchange
      material%spiralization(:,:,1,1) = dmi
      material%diagnostics%channel_dynamics_parameters_validated = .true.
      material%diagnostics%small_q_energy_validated = .true.
   end subroutine make_material

   subroutine atomistic_baseline(direction,bonds,matrix,external,baseline_energy,baseline_field)
      real(dblprec), intent(in) :: direction(:,:), matrix(:,:,:)
      integer, intent(in) :: bonds(:,:)
      real(dblprec), intent(in), optional :: external(:,:)
      real(dblprec), intent(out) :: baseline_energy, baseline_field(:,:)
      integer :: atom, bond, i, j

      baseline_energy = 0.0_dblprec
      baseline_field = 0.0_dblprec
      do bond = 1,size(bonds,2)
         i = bonds(1,bond)
         j = bonds(2,bond)
         baseline_energy = baseline_energy - dot_product(direction(:,i), &
            matmul(matrix(:,:,bond),direction(:,j)))
         baseline_field(:,i) = baseline_field(:,i) + &
            matmul(matrix(:,:,bond),direction(:,j))/(mub_si*1.7_dblprec)
         baseline_field(:,j) = baseline_field(:,j) + &
            matmul(transpose(matrix(:,:,bond)),direction(:,i))/(mub_si*1.7_dblprec)
      end do
      if (present(external)) then
         do atom = 1,size(direction,2)
            baseline_energy = baseline_energy - mub_si*1.7_dblprec* &
               dot_product(external(:,atom),direction(:,atom))
            baseline_field(:,atom) = baseline_field(:,atom)+external(:,atom)
         end do
      end if
   end subroutine atomistic_baseline

   subroutine hybrid_finite_difference(hybrid,fine,coarse,matrix,is_fine,index,tangent,numerical)
      type(static_hybrid_operator_type), intent(in) :: hybrid
      real(dblprec), intent(in) :: fine(:,:), coarse(:,:,:), matrix(:,:,:), tangent(3)
      logical, intent(in) :: is_fine
      integer, intent(in) :: index
      real(dblprec), intent(out) :: numerical
      type(static_hybrid_energy_type) :: plus_energy, minus_energy
      integer :: status
      character(len=512) :: message
      real(dblprec), parameter :: epsilon_fd = 2.0d-7
      real(dblprec) :: plus_fine(3,size(fine,2)), minus_fine(3,size(fine,2))
      real(dblprec) :: plus_coarse(3,1,size(coarse,3)), minus_coarse(3,1,size(coarse,3))
      real(dblprec) :: ff(3,size(fine,2)), cf(3,1,size(coarse,3))

      plus_fine = fine
      minus_fine = fine
      plus_coarse = coarse
      minus_coarse = coarse
      if (is_fine) then
         plus_fine(:,index) = plus_fine(:,index)+epsilon_fd*tangent
         minus_fine(:,index) = minus_fine(:,index)-epsilon_fd*tangent
         plus_fine(:,index) = plus_fine(:,index)/norm3(plus_fine(:,index))
         minus_fine(:,index) = minus_fine(:,index)/norm3(minus_fine(:,index))
      else
         plus_coarse(:,1,index) = plus_coarse(:,1,index)+epsilon_fd*tangent
         minus_coarse(:,1,index) = minus_coarse(:,1,index)-epsilon_fd*tangent
         plus_coarse(:,1,index) = plus_coarse(:,1,index)/norm3(plus_coarse(:,1,index))
         minus_coarse(:,1,index) = minus_coarse(:,1,index)/norm3(minus_coarse(:,1,index))
      end if
      call evaluate_static_hybrid_operator(hybrid,plus_fine,plus_coarse,matrix,ff,cf, &
         plus_energy,status,message)
      call evaluate_static_hybrid_operator(hybrid,minus_fine,minus_coarse,matrix,ff,cf, &
         minus_energy,status,message)
      numerical = (plus_energy%total_j-minus_energy%total_j)/(2.0_dblprec*epsilon_fd)
   end subroutine hybrid_finite_difference

   real(dblprec) function max_tangent_atoms(direction,field,active)
      real(dblprec), intent(in) :: direction(:,:), field(:,:)
      logical, intent(in) :: active(:)
      integer :: atom
      real(dblprec) :: tangent(3)

      max_tangent_atoms = 0.0_dblprec
      do atom = 1,size(active)
         if (.not. active(atom)) cycle
         tangent = field(:,atom)-direction(:,atom)*dot_product(direction(:,atom),field(:,atom))
         max_tangent_atoms = max(max_tangent_atoms,norm3(tangent))
      end do
   end function max_tangent_atoms

   real(dblprec) function max_tangent_blocks(direction,field,active)
      real(dblprec), intent(in) :: direction(:,:,:), field(:,:,:)
      logical, intent(in) :: active(:)
      integer :: block
      real(dblprec) :: tangent(3)

      max_tangent_blocks = 0.0_dblprec
      do block = 1,size(active)
         if (.not. active(block)) cycle
         tangent = field(:,1,block)-direction(:,1,block)* &
            dot_product(direction(:,1,block),field(:,1,block))
         max_tangent_blocks = max(max_tangent_blocks,norm3(tangent))
      end do
   end function max_tangent_blocks

   pure real(dblprec) function relative_error(actual,expected)
      real(dblprec), intent(in) :: actual, expected
      relative_error = abs(actual-expected)/max(abs(actual),abs(expected),tiny(1.0_dblprec))
   end function relative_error

   pure real(dblprec) function norm3(vector)
      real(dblprec), intent(in) :: vector(3)
      norm3 = sqrt(sum(vector*vector))
   end function norm3

   subroutine check(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         failures = failures+1
         write(*,'(a)') 'FAIL: '//trim(message)
      end if
   end subroutine check

   subroutine check_close(actual,expected,relative_tolerance,message)
      real(dblprec), intent(in) :: actual, expected, relative_tolerance
      character(len=*), intent(in) :: message
      real(dblprec) :: scale

      scale = max(abs(expected),abs(actual),tiny(1.0_dblprec))
      call check(abs(actual-expected) <= relative_tolerance*scale,message)
      if (abs(actual-expected) > relative_tolerance*scale) then
         write(*,'(a,2es24.15)') '  actual/expected: ',actual,expected
      end if
   end subroutine check_close

end program test_static_hybrid_operator
