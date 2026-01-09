!-------------------------------------------------------------------------------
! MODULE: SpinTorques
!> @brief
!> Routines for calculating \f$\frac{\partial\mathbf{m}}{\partial\mathbf{r}}\f$.
!> Needed for generalized spin-torque term
!> \f$\left(\mathbf{m} \times \left(\mathbf{m}\times\frac{\partial\mathbf{m}}{\partial \mathbf{r}}\right)\right)\f$
!> @details For calculating the spin transfer torques, we use the standard adiabatic
!> and non-adiabatic terms as introduced in the LLG equation by Zhang & Li. PRL 93, 127204 (2004).
!> The terms are rewritten to suit the LL equations used in UppASD and here we
!> use the same formulas as Schieback et. al, Eur. Phys. J. B 59, 429, (2007)
!> Currently the torques are only calculated as their respective fields
!> (i.e. missing the preceeding \f$\mathbf{m} \times\f$) since that is taken care of in the Depondt solver
!>
!> === CANONICAL CURRENT DENSITY APPROACH ===
!> The module uses a canonical current density representation where all currents are
!> internally stored and used in physical units (A/m²). This provides:
!> - Single source of truth: jdens_cell (uniform) or jdens_site (site-resolved)
!> - Consistent handling across Zhang-Li, Slonczewski, and SHE torques
!> - Clear precedence: direct jdens input > legacy jvec conversion
!> - Transparent unit conversions with detailed logging
!>
!> === USER INPUT OPTIONS ===
!> Users can specify currents in three ways (in order of precedence):
!> 1. Direct physical current density (RECOMMENDED):
!>    jdens  1.0e12  0.0  0.0    ! A/m² - used directly, no conversion
!>
!> 2. Legacy dimensionless current (backward compatible):
!>    jvec  1.0  0.0  0.0         ! code units - converted to A/m² automatically
!>    jsite Y                     ! optional: site-resolved via jvecfile
!>
!> 3. Site-resolved physical current (future extension):
!>    Currently implemented via legacy jvecfile + automatic conversion
!>
!> === SHE SIGMA DETERMINATION ===
!> For Spin Hall Effect torques, the spin accumulation direction (sigma) is determined by:
!> 1. If she_n_vec provided: sigma = n × j (proper physics for interface SHE)
!>    she_n_vec  0.0  0.0  1.0    ! surface normal
!> 2. If she_n_vec not provided: sigma = ĵ (legacy behavior, useful for testing)
!>
!> All current density resolution is handled by resolve_current_density() which
!> prints comprehensive diagnostics about the current source, conversions applied,
!> and final magnitudes used.
!>
!> @author
!> Anders Bergman
!> @notes
!> J. Chico
!> - Added the SHE torque and the calculation of the current density
!> - Added a general model for SOT
!> A. Bergman
!> - Implemented canonical current density approach (2026)
!> - Added resolve_current_density() for single source of truth
!> - Updated all STT mechanisms to use canonical jdens
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module SpinTorques

   use Parameters
   use Profiling

   implicit none

   !Spin-transfer torque inputs
   character(len=1) :: STT             !< Treat spin transfer torque? (Y/N)
   character(len=1) :: jsite           !< Treat site dependent jvec
   character(len=1) :: do_she          !< Treat SHE spin transfer torque
   character(len=1) :: do_sot          !< Treat SHE spin transfer torque
   character(len=1) :: sot_site_pol    !< Treat site dependent jvec
   character(len=35) :: jvecfile       !< File name for the site dependent jvec
   character(len=35) :: sot_site_file  !< File name for the site dependent SOT polarization
   real(dblprec) :: adibeta      !< Adiabacity parameter for STT
   real(dblprec) :: spin_pol     !< Spin polarization
   real(dblprec) :: SHE_angle    !< Spin Hall angle
   real(dblprec), dimension(3) :: she_n_vec !< Surface/interface normal for SHE (for sigma = n x j)
   real(dblprec) :: sot_field    !< SOT field-like torque strength (Tesla)
   real(dblprec) :: thick_ferro  !< Thickness of ferromagnetic layer (t_f/alat)
   real(dblprec) :: sot_damping  !< SOT damping-like torque strength (Tesla)
   real(dblprec) :: jscale       !< Scale factor for local current density
   ! === LEGACY INPUT VARIABLES (kept for backward compatibility) ===
   real(dblprec), dimension(3) :: jvec !< Legacy spin current vector (dimensionless code units) - converted to jdens
   real(dblprec), dimension(3) :: jdens !< Current density input/output (A/m^2) - may be user input or converted from jvec
   logical :: jdens_user_provided !< Flag: was jdens provided by user in input file? (tracked before set_curr_density conversion)
   real(dblprec) :: stt_dens_conv !< Conversion factor: jdens = stt_dens_conv * jvec (A/m^2 per code unit)
   ! === CANONICAL CURRENT DENSITY VARIABLES (single source of truth, used internally) ===
   logical :: have_jdens !< Flag: was jdens provided directly by user (vs converted from jvec)?
   logical :: use_jdens_site !< Flag: using site-resolved current density (jdens_site) vs uniform (jdens_cell)?
   real(dblprec), dimension(3) :: jdens_cell !< Canonical cell-uniform current density (A/m^2) - ALL torques use this
   real(dblprec), dimension(:,:), allocatable :: jdens_site !< Canonical site-resolved current density (3,Natom) (A/m^2) - used if available
   real(dblprec) :: b_rt_fac      !< Prefactor for Slonczewski STT
   real(dblprec) :: b_zhang_li_fac !< Prefactor for Zhang-Li STT (atomistic Schieback formula with alat conversion)
   real(dblprec) :: sot_rt_fac      !< Prefactor for Spin Hall Torque (revamped SHE-torque)
   real(dblprec), dimension(3) :: she_sigma_vec !< Polarization vector for SHE
   real(dblprec), dimension(3) :: sot_pol_vec  !< Polarization vector for SOT
   real(dblprec), dimension(:,:), allocatable :: sitenatomjvec !< Site dependent spin current vector
   real(dblprec), dimension(:,:), allocatable :: sitenatom_stt_pol
   real(dblprec), dimension(:), allocatable :: sitenatom_stt_jcur
   real(dblprec), dimension(:,:), allocatable :: sitenatom_sot_pol

   !Spin-transfer torque data arrays
   real(dblprec), dimension(:), allocatable  :: stt_prefac !< Prefactor for the STT
   real(dblprec), dimension(:,:,:), allocatable :: dmomdr  !< Current magnetic moment vector
   real(dblprec), dimension(:,:,:), allocatable :: btorque !< Spin transfer torque
   real(dblprec), dimension(:,:,:), allocatable :: she_btorque !< SHE spin transfer torque
   real(dblprec), dimension(:,:,:), allocatable :: sot_btorque !< SOT spin transfer torque

   public

contains

   !---------------------------------------------------------------------------
   !> @brief
   !> Wrapper for the actual calculation of the spin torques
   !> @details This routine calculates all enabled spin-transfer torques using
   !> the canonical current density (jdens_cell or jdens_site) as the single
   !> source of truth. The routine:
   !> - For Zhang-Li (stt='A'): extracts direction for gradient, applies magnitude
   !> - For Slonczewski (stt='S'): uses canonical magnitude with polarization from sitenatom_stt_pol
   !> - For SHE (do_she='Y'): uses sigma from she_sigma_vec (computed from n×j or ĵ)
   !> - For SOT (do_sot='Y'): uses explicit sot_field/sot_damping parameters (Tesla)
   !>   NOTE: SOT is phenomenological and independent of current density
   !>
   !> IMPORTANT: This routine assumes resolve_current_density() has already been
   !> called during initialization to populate jdens_cell/jdens_site from user inputs.
   !>
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine calculate_spintorques(Natom, Mensemble,lambda1_array,emom,mmom)
      !
      use Gradients, only : differentiate_moments
      use Constants, only : gama
      !
      implicit none
      !
      !.. Input variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom  !< Current magnetic moment vector$
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmom  !< Current magnetic moment vector$
      !
      ! .. Local variables
      integer :: j, k
      integer :: iatom
      real(dblprec) :: jmag, jdir_mag
      real(dblprec), dimension(3) :: jvec_i, jdir_i
      real(dblprec), parameter :: eps = 1.0e-15_dblprec
      real(dblprec) :: jmag_loc
      real(dblprec), dimension(3) :: jtmp
      real(dblprec), parameter :: eps_j = 1.0e-30_dblprec  ! tighter than eps for robust |j| test
      
      ! Guard: if no STT/SHE/SOT is enabled, exit early
      if (stt /= 'A' .and. stt /= 'S' .and. stt /= 'F' .and. do_she /= 'Y' .and. do_sot /= 'Y') then
         return
      endif
      
      if (stt == 'A') then
         ! Zhang–Li gradient STT (Schieback / Zhang–Li atomistic form)
         !
         ! We feed the gradient kernel with the drift-velocity-like coefficient in lattice units:
         !
         !   u_latt(i) = (b_zhang_li_fac / mmom_ref(i)) * jdens(i)
         !
         ! where b_zhang_li_fac must be defined consistently with your coordinate convention:
         !
         !   If coord is in lattice units (dimensionless) and b_zhang_li_fac was built as
         !     b_zhang_li_fac = (P*g/(2e)) * (V_at/alat)
         !   then this is correct as written (NO further /alat here).
         !
         ! The gradient routine returns dmomdr ≈ (u_latt · ∇_latt) m, with units 1/time.
      
         ! Build sitenatomjvec(:,iatom) = u_latt(:,iatom) (NOT normalized)
         if (use_jdens_site) then
            do iatom = 1, Natom
               jmag_loc = norm2(jdens_site(:,iatom))
               if (jmag_loc > eps_j .and. mmom(iatom,1) > eps) then
                  sitenatomjvec(:,iatom) = b_zhang_li_fac * jdens_site(:,iatom) / mmom(iatom,1)
               else
                  sitenatomjvec(:,iatom) = 0.0_dblprec
               endif
            enddo
         else
            jmag_loc = norm2(jdens_cell)
            if (jmag_loc > eps_j) then
               jtmp = jdens_cell
            else
               jtmp = 0.0_dblprec
            endif
         
            do iatom = 1, Natom
               if (mmom(iatom,1) > eps) then
                  sitenatomjvec(:,iatom) = b_zhang_li_fac * jtmp / mmom(iatom,1)
               else
                  sitenatomjvec(:,iatom) = 0.0_dblprec
               endif
            enddo
         endif
      
         ! Compute directional derivative using u_latt
         call differentiate_moments(Natom, Mensemble, emom, dmomdr, sitenatomjvec)
      
         ! Apply the Zhang–Li algebra. dmomdr already includes the physical amplitude.
         !$omp parallel do default(shared) private(j,k)
         do j = 1, Natom
            stt_prefac(j) = -(1.0_dblprec + adibeta * lambda1_array(j))
            do k = 1, Mensemble
               btorque(1,j,k) = (lambda1_array(j) - adibeta) * dmomdr(1,j,k)
               btorque(2,j,k) = (lambda1_array(j) - adibeta) * dmomdr(2,j,k)
               btorque(3,j,k) = (lambda1_array(j) - adibeta) * dmomdr(3,j,k)
            enddo
         enddo
         !$omp end parallel do
      
         ! Adds the m × dmomdr term multiplied by stt_prefac(j)
         call mom_cross_dmomdr(Natom, Mensemble, emom)
         
         ! Correct for gamma: The Zhang-Li prefactor uses the electron g-factor (g_e_abs)
         ! without the gyromagnetic ratio (gamma). However, depondt.f90 multiplies all
         ! effective fields by gamma during time evolution. To avoid double-counting gamma,
         ! we divide btorque by gamma here. This ensures consistent scaling with the
         ! Hamiltonian field (beff), which is NOT multiplied by gamma in its definition.
         !$omp parallel do default(shared) private(j,k) collapse(2)
         do k = 1, Mensemble
            do j = 1, Natom
               btorque(:,j,k) = btorque(:,j,k) / gama
            enddo
         enddo
         !$omp end parallel do
         ! external_field=external_field+btorque
      else if(stt=='S') then
         ! Slonczewski torque a la Evans (J. Phys.: Condens. Matter 35 025801 (2023))
         call slonczewski_field(Natom, Mensemble,emom)
         !call mom_cross_mfixed(Natom, Mensemble,emom)
         ! external_field=external_field+btorque
      else if(stt=='F') then
         ! Fixed layer (Older Slonczewski implementation)
         call mom_cross_mfixed(Natom, Mensemble,emom)
         ! external_field=external_field+btorque
      end if

      ! Calculates the spin hall effect generated spin transfer torque
      if (do_she=='Y') then
         call SHE_torque(Natom,Mensemble,lambda1_array,emom,mmom)
      endif

      ! Calculates the general field-like and damping-like terms that describe the SOT
      if (do_sot=='Y') then
         call SOT_torque(Natom,Mensemble,lambda1_array,emom)
      endif

   end subroutine calculate_spintorques

   !---------------------------------------------------------------------------
   !> @brief
   !> Calculate
   !> \f$\left(\mathbf{m}\times\left(\mathbf{u}\cdot\frac{\partial}{\partial\mathbf{r}}\right)\mathbf{m}\right)\f$
   !> (which then ends up as one part of the spin transfer torque) for atomic damping dependence
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine mom_cross_dmomdr(Natom, Mensemble,emomM)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector

      integer :: iatom, k

      !$omp parallel do default(shared) private(iatom,k)
      do iatom=1, Natom
         do k=1, Mensemble
            btorque(1,iatom,k)=btorque(1,iatom,k)+stt_prefac(iatom)*(emomM(2,iatom,k)*dmomdr(3,iatom,k) &
               -emomM(3,iatom,k)*dmomdr(2,iatom,k))
            btorque(2,iatom,k)=btorque(2,iatom,k)+stt_prefac(iatom)*(emomM(3,iatom,k)*dmomdr(1,iatom,k) &
               -emomM(1,iatom,k)*dmomdr(3,iatom,k))
            btorque(3,iatom,k)=btorque(3,iatom,k)+stt_prefac(iatom)*(emomM(1,iatom,k)*dmomdr(2,iatom,k) &
               -emomM(2,iatom,k)*dmomdr(1,iatom,k))
         end do
      end do
      !$omp end parallel do

   end subroutine mom_cross_dmomdr

   !---------------------------------------------------------------------------
   !> @brief
   !> Calculates the spin transfer torque for currents passing through a fixed ferromagnetic layer
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine mom_cross_mfixed(Natom, Mensemble,emomM)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector

      integer :: iatom, k

      btorque=0.0_dblprec
      do k=1, Mensemble
         do iatom=1, Natom
            btorque(1,iatom,k)=emomM(2,iatom,k)*sitenatomjvec(3,iatom)-emomM(3,iatom,k)*sitenatomjvec(2,iatom)
            btorque(2,iatom,k)=emomM(3,iatom,k)*sitenatomjvec(1,iatom)-emomM(1,iatom,k)*sitenatomjvec(3,iatom)
            btorque(3,iatom,k)=emomM(1,iatom,k)*sitenatomjvec(2,iatom)-emomM(2,iatom,k)*sitenatomjvec(1,iatom)
         end do
      end do

   end subroutine mom_cross_mfixed
   !---------------------------------------------------------------------------
   !> @brief
   !> Calculates the spin transfer torque for currents passing through a fixed ferromagnetic layer
   !> Updated formalism according to J. Phys.: Condens. Matter 35 (2023) 025801
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine slonczewski_field(Natom, Mensemble,emom)
      use math_functions, only : f_cross_product
      use damping, only : lambda1_array
      use MomentData, only : mmom
      !B = BSTT (p−α m×p)+BSTT (m×p+α p)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom  !< Current magnetic moment vector

      real(dblprec) :: stt_asym, stt_pfac, stt_dot, jmag_i
      integer :: iatom, k
      real(dblprec), parameter :: eps = 1.0e-15_dblprec

      stt_asym = spin_pol**2
      btorque=0.0_dblprec
      do k=1, Mensemble
         do iatom=1, Natom
            ! Get current magnitude from canonical source
            if (use_jdens_site) then
               jmag_i = norm2(jdens_site(:,iatom))
            else
               jmag_i = norm2(jdens_cell)
            endif
            
            stt_dot = dot_product(emom(:,iatom,k), sitenatom_stt_pol(:,iatom))
            ! Prefactor involves density, physical constants, and area/length
            ! Here we need to divide with local moment
            ! Use canonical current magnitude (already in A/m^2)
            stt_pfac = b_rt_fac*jmag_i*spin_pol/(1.0_dblprec+stt_asym*stt_dot)/mmom(iatom,k)
            ! First add precessional contribution (B^S_P * (p - alpha m x p ))
            btorque(:,iatom,k) = btorque(:,iatom,k) + stt_pfac * adibeta * ( &
               sitenatom_stt_pol(:,iatom) &
               - lambda1_array(iatom) * f_cross_product(emom(:,iatom,k),sitenatom_stt_pol(:,iatom)))
            ! Then add damping contribution (B^S_R * (m x p + alpha p))
            btorque(:,iatom,k) = btorque(:,iatom,k) + stt_pfac * ( &
               f_cross_product(emom(:,iatom,k),sitenatom_stt_pol(:,iatom)) + &
               lambda1_array(iatom) * sitenatom_stt_pol(:,iatom) )

         end do
      end do

   end subroutine slonczewski_field

   !-----------------------------------------------------------------------------
   !> @brief
   !> Calculate the SHE torque generated by a current passing through a non magnetic material with SOC
   !> Note: Revamped formalism according to Meo 2023 (J. Phys. Condens. Matter 35 215801)
   !
   !> @author
   !> Jonathan Chico, Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine SHE_torque(Natom,Mensemble,lambda1_array,emom,mmom)
      use math_functions, only : f_cross_product

      implicit none

      ! .. Input variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom  !< Current magnetic moment vector$
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmom  !< Current magnetic moment vector$

      ! ... Local variables
      integer :: iatom, k
      real(dblprec) :: she_fact

      she_btorque=0.0_dblprec
      ! Factor for the strenght of the spin hall torque
      if (thick_ferro == 0.0_dblprec) then
         write(*,*) 'ERROR: Thickness of ferromagnetic layer (thick_ferro) must be non-zero.'
         stop
      endif
      ! Previous factor:
      ! she_fact=she_angle/(spin_pol*thick_ferro)
      ! New factor according to Meo 2023 below (see also set_curr_density)

      she_btorque = 0.0_dblprec
      !$omp parallel do default(shared) private(iatom,k)
      do k=1, Mensemble
         do iatom=1, Natom
            she_fact = sot_rt_fac / mmom(iatom,k)
            ! Precession part (B^S_P * (p - alpha m x p ))
            she_btorque(:,iatom,k) = she_btorque(:,iatom,k) + she_fact * adibeta * ( &
               she_sigma_vec - lambda1_array(iatom) * f_cross_product(emom(:,iatom,k),she_sigma_vec) &
               )
            ! Then add damping contribution (B^S_R * (m x p + alpha p))
            she_btorque(:,iatom,k) = she_btorque(:,iatom,k) + she_fact * ( &
               f_cross_product(emom(:,iatom,k),she_sigma_vec) + &
               lambda1_array(iatom) * she_sigma_vec )
            ! she_btorque(1,iatom,k) =-she_fact*sitenatomjvec(1,iatom)*emom(3,iatom,k)-lambda1_array(iatom)*mmom(iatom,k)*sitenatomjvec(2,iatom)
            ! she_btorque(2,iatom,k) =-she_fact*sitenatomjvec(2,iatom)*emom(3,iatom,k)+lambda1_array(iatom)*mmom(iatom,k)*sitenatomjvec(1,iatom)
            ! she_btorque(3,iatom,k) = she_fact*sitenatomjvec(2,iatom)*emom(2,iatom,k)+she_fact*sitenatomjvec(1,iatom)*emom(1,iatom,k)
            ! print '(g12.4, 3g12.4)', she_fact, she_btorque(:,iatom,k)
            ! print '(g12.4, 3g12.4)', she_fact, she_sigma_vec
         enddo
      enddo
      !$omp end parallel do

   end subroutine SHE_torque

   !----------------------------------------------------------------------------
   ! SUBROUTINE
   !> @brief Calculates the general phenomenological form of SOT for the LLG
   !> @details Computes spin-orbit torque using phenomenological field-like and
   !> damping-like components. This is INDEPENDENT of current density (unlike STT/SHE)
   !> and uses explicit effective field strengths.
   !>
   !> TORQUE FORM:
   !> The SOT can always be written as:
   !> - Field-like term: τ_FL * m × P
   !> - Damping-like term: τ_DL * m × (m × P)
   !>
   !> where:
   !> - P = sitenatom_sot_pol (polarization vector, dimensionless, normalized)
   !> - τ_FL = sot_field (field-like strength, Tesla)
   !> - τ_DL = sot_damping (damping-like strength, Tesla)
   !>
   !> IMPLEMENTATION (including α-dependent cross-terms):
   !> sot_btorque = -(τ_FL - α*τ_DL) * P - (τ_DL + α*τ_FL) * (m × P)
   !>
   !> where α = lambda1_array (Gilbert damping parameter, dimensionless)
   !>
   !> UNITS:
   !> - sot_field: Tesla (effective magnetic field)
   !> - sot_damping: Tesla (effective magnetic field)
   !> - sot_btorque: Tesla (added to effective field in LLG)
   !>
   !> USER INPUT:
   !> sot_field    0.01   # Tesla - field-like SOT strength
   !> sot_damping  0.05   # Tesla - damping-like SOT strength
   !> sot_pol_vec  1.0  0.0  0.0   # Polarization direction (normalized automatically)
   !>
   !> POLARIZATION SOURCE:
   !> - Uniform: sot_pol_vec (module variable)
   !> - Site-resolved: sitenatom_sot_pol (from sot_site_file)
   !>
   !> NOTE: This is a phenomenological model. For microscopic SOT from current density
   !> (e.g., spin Hall effect), use do_she='Y' instead, which derives torques from
   !> canonical jdens and computes proper σ = n × j.
   !>
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine SOT_torque(Natom,Mensemble,lambda1_array,emom)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom  !< Current magnetic moment vector$

      integer :: iatom, ens
      sot_btorque=0.0_dblprec
      !$omp parallel do default(shared) private(iatom,ens)
      do ens=1,Mensemble
         do iatom=1,Natom
            !
            sot_btorque(1,iatom,ens)=-(sot_field-lambda1_array(iatom)*sot_damping)*sitenatom_sot_pol(1,iatom)&
            -(sot_damping+lambda1_array(iatom)*sot_field)*&
            (emom(2,iatom,ens)*sitenatom_sot_pol(3,iatom)-emom(3,iatom,ens)*sitenatom_sot_pol(2,iatom))
            !
            sot_btorque(2,iatom,ens)=-(sot_field-lambda1_array(iatom)*sot_damping)*sitenatom_sot_pol(2,iatom)&
            -(sot_damping+lambda1_array(iatom)*sot_field)*&
            (emom(3,iatom,ens)*sitenatom_sot_pol(1,iatom)-emom(1,iatom,ens)*sitenatom_sot_pol(2,iatom))
            !
            sot_btorque(3,iatom,ens)=-(sot_field-lambda1_array(iatom)*sot_damping)*sitenatom_sot_pol(3,iatom)&
            -(sot_damping+lambda1_array(iatom)*sot_field)*&
            (emom(1,iatom,ens)*sitenatom_sot_pol(2,iatom)-emom(2,iatom,ens)*sitenatom_sot_pol(1,iatom))
         enddo
      enddo
      !$omp end parallel do

   end subroutine SOT_torque

   !----------------------------------------------------------------------------
   !> @brief
   !> Allocation and deallocation of the arrays for the STT calculation
   !
   !> @author
   !> Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_stt_data(Natom,Mensemble,flag)

      implicit none

      integer, intent(in) :: flag
      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble

      integer :: i_stat, i_all

      if (flag==1) then
         ! Allocate the SHE torque field
         if (do_she=='Y') then
            allocate(she_btorque(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(she_btorque))*kind(she_btorque),'she_btorque','allocate_stt_data')
            she_btorque=0.0_dblprec
         endif

         ! Always allocate btorque (used by prn_fields even when STT is disabled)
         allocate(btorque(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(btorque))*kind(btorque),'btorque','allocate_stt_data')
         btorque=0.0_dblprec
         
         if (stt=='A') then
            allocate(stt_prefac(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(stt_prefac))*kind(stt_prefac),'stt_prefac','allocate_stt_sata')
            stt_prefac=0.0_dblprec
            allocate(dmomdr(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(dmomdr))*kind(dmomdr),'dmomdr','allocate_stt_data')
            dmomdr=0.0_dblprec
         endif
         if (do_sot=='Y') then
            allocate(sot_btorque(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(sot_btorque))*kind(sot_btorque),'sot_btorque','allocate_stt_data')
            sot_btorque=0.0_dblprec
         endif

      else
         ! Always deallocate btorque (always allocated in flag==1 block)
         i_all=-product(shape(btorque))*kind(btorque)
         deallocate(btorque,stat=i_stat)
         call memocc(i_stat,i_all,'btorque','allocate_stt_data')
         
         if (do_she=='Y') then
            i_all=-product(shape(she_btorque))*kind(she_btorque)
            deallocate(she_btorque,stat=i_stat)
            call memocc(i_stat,i_all,'she_btorque','allocate_stt_data')
         endif
         if (stt=='A') then
            i_all=-product(shape(stt_prefac))*kind(stt_prefac)
            deallocate(stt_prefac,stat=i_stat)
            call memocc(i_stat,i_all,'stt_prefac','allocate_stt_data')
            i_all=-product(shape(dmomdr))*kind(dmomdr)
            deallocate(dmomdr,stat=i_stat)
            call memocc(i_stat,i_all,'dmomdr','allocate_stt_data')
         endif
         if(do_sot=='Y') then
            i_all=-product(shape(sot_btorque))*kind(sot_btorque)
            deallocate(sot_btorque,stat=i_stat)
            call memocc(i_stat,i_all,'sot_btorque','allocate_stt_data')
         endif
      endif

   end subroutine allocate_stt_data

   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_stt(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword,cache
      integer :: rd_len, i_err, i_errb
      logical :: comment

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

            case('stt')
               read(ifile,*,iostat=i_err) stt
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('jvec')
               read(ifile,*,iostat=i_err) jvec
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('jdens')
               read(ifile,*,iostat=i_err) jdens
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ! Track that jdens was provided by user (before set_curr_density might convert jvec)
               if (norm2(jdens) > 1.0e-15_dblprec) then
                  jdens_user_provided = .true.
               endif

            case('jsite')
               read(ifile,*,iostat=i_err) jsite
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('jscale')
               read(ifile,*,iostat=i_err) jscale
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('jvecfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               jvecfile=adjustl(trim(cache))

            case('adibeta')
               read(ifile,*,iostat=i_err) adibeta
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_she') ! Consider the SHE generated torque
               read(ifile,*,iostat=i_err) do_she
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('spin_pol') ! Spin polarization
               read(ifile,*,iostat=i_err) spin_pol
               if(i_err/=0) write(*,*) 'ERROR; Reading ',trim(keyword),' data',i_err

            case('she_angle') ! Spin Hall Angle
               read(ifile,*,iostat=i_err) she_angle
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('she_n_vec') ! Surface/interface normal for SHE
               read(ifile,*,iostat=i_err) she_n_vec
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('thick_ferro') ! Thickness of the ferromagnetic layer
               read(ifile,*,iostat=i_err) thick_ferro
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('do_sot') ! Consider the general SOT torque
               read(ifile,*,iostat=i_err) do_sot
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sot_pol_vec')
               read(ifile,*,iostat=i_err) sot_pol_vec
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sot_field')
               read(ifile,*,iostat=i_err) sot_field
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sot_damping')
               read(ifile,*,iostat=i_err) sot_damping
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sot_site_pol')
               read(ifile,*,iostat=i_err) sot_site_pol
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sot_site_file')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               sot_site_file=adjustl(trim(cache))
            end select
         end if

         ! End of file
         if (i_errb==20) goto 20
         ! End of row
         if (i_errb==10) goto 10
      end do

      20  continue

      rewind(ifile)
      return
   end subroutine read_parameters_stt

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: init_stt
   !> @brief Set default values for STT calculations
   !-----------------------------------------------------------------------------
   subroutine init_stt()
      !
      implicit none

      !Spin transfer torque
      STT           = "N"
      jsite         = "N"
      do_she        = "N"
      do_sot        = "N"
      adibeta       = 0.0_dblprec
      spin_pol      = 0.0_dblprec
      jvecfile      = 'jvecfile'
      sot_field     = 0.0_dblprec
      she_angle     = 0.0_dblprec
      thick_ferro   = 0.0_dblprec
      sot_damping   = 0.0_dblprec
      sot_site_pol  = "N"
      sot_site_file = 'site_pol'
      jvec          = (/0.0_dblprec,0.0_dblprec,0.0_dblprec/)
      jdens         = (/0.0_dblprec,0.0_dblprec,0.0_dblprec/)
      jscale        = 1.0_dblprec
      jdens_user_provided = .false.
      ! Track whether jdens was provided by user (before set_curr_density converts jvec)
      jdens_user_provided = .false.
      ! Canonical current variables
      have_jdens    = .false.
      use_jdens_site = .false.
      jdens_cell    = (/0.0_dblprec,0.0_dblprec,0.0_dblprec/)
      she_n_vec     = (/0.0_dblprec,0.0_dblprec,0.0_dblprec/)

   end subroutine init_stt

   !---------------------------------------------------------------------------
   !> @brief
   !> Read site-dependent currents from file or set uniform values
   !> @details Populates sitenatomjvec(3,Natom) with dimensionless current vectors.
   !>
   !> BEHAVIOR:
   !> - If jsite='Y': reads site-resolved vectors from jvecfile
   !>   File format: site_index  jx  jy  jz  (one line per site)
   !> - If jsite='N': sets all sites to uniform jvec value
   !>
   !> LEGACY ARRAYS POPULATED (for backward compatibility):
   !> - sitenatomjvec(3,Natom): raw dimensionless current vectors
   !> - sitenatom_stt_pol(3,Natom): normalized polarization directions
   !> - sitenatom_stt_jcur(Natom): current magnitudes in A/m² (using stt_dens_conv)
   !>
   !> IMPORTANT: These legacy arrays are NOT the canonical source of truth.
   !> The canonical current density is established by resolve_current_density()
   !> which converts sitenatomjvec to jdens_site (physical A/m²) and uses that
   !> for all torque calculations. The legacy arrays are kept only for:
   !> - Backward compatibility with fixed-layer mode (stt='F')
   !> - Polarization direction for Slonczewski (sitenatom_stt_pol)
   !>
   !> CALL ORDER:
   !> Must be called after set_curr_density() (to have stt_dens_conv)
   !> and before resolve_current_density() (which converts sitenatomjvec → jdens_site).
   !>
   !> @author
   !> Jonathan Chico, Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_jvecfile(Natom)

      implicit none

      integer, intent(in) :: Natom

      integer :: i,flines, isite, i_stat
      real(dblprec) :: pnorm

      allocate(sitenatomjvec(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(sitenatomjvec))*kind(sitenatomjvec),'sitenatomjvec','read_jvecfile')

      if (jsite=='Y') then
         open(ifileno, file=trim(jvecfile))
         flines=0
         ! Pre-read file to get number of lines
         do
            read(ifileno,*,end=200)  isite
            flines=flines+1
         end do

         200 continue

         rewind(ifileno)

         write(*,'(2x,a)') 'Reading site dependent currents'

         ! If the size of the file is NATOM then there is no problem
         if ( Natom.eq.flines ) then

            do i=1, flines
               read(ifileno,*) isite, sitenatomjvec(1,isite), sitenatomjvec(2,isite), sitenatomjvec(3,isite)
            end do
         else
            write(*,*) 'WARNING: Size of the SITEATOMJVEC is not NATOM'
            do i=1, flines
               read(ifileno,*) isite, sitenatomjvec(1,isite), sitenatomjvec(2,isite), sitenatomjvec(3,isite)
            end do

         end if

         close(ifileno)

         ! If there is no site dependent current set it to be constant
      else
         do i=1, Natom
            sitenatomjvec(1,i)=jvec(1)
            sitenatomjvec(2,i)=jvec(2)
            sitenatomjvec(3,i)=jvec(3)
         enddo
      endif

      ! From jvec, calculate normalized polarization vector and magnitude
      allocate(sitenatom_stt_pol(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(sitenatom_stt_pol))*kind(sitenatom_stt_pol),'sitenatom_stt_pol','read_jvecfile')
      allocate(sitenatom_stt_jcur(Natom),stat=i_stat)
      call memocc(i_stat,product(shape(sitenatom_stt_jcur))*kind(sitenatom_stt_jcur),'sitenatom_stt_jcur','read_jvecfile')

      do i=1, Natom
         pnorm = norm2(sitenatomjvec(:,i))
         sitenatom_stt_jcur(i) = pnorm * stt_dens_conv
         sitenatom_stt_pol(:,i) = sitenatomjvec(:,i) / (pnorm + 1.0e-15_dblprec)
      end do
   end subroutine read_jvecfile

   !---------------------------------------------------------------------------
   !> @brief
   !> Read site dependent currents
   !
   !> @author
   !> Jonathan Chico, Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_sot_pol_site(Natom)

      implicit none

      integer, intent(in) :: Natom

      integer :: i,flines, isite, i_stat

      allocate(sitenatom_sot_pol(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(sitenatom_sot_pol))*kind(sitenatom_sot_pol),'sitenatom_sot_pol','read_sot_pol_site')

      if (sot_site_pol=='Y') then
         open(ifileno, file=trim(sot_site_file))
         flines=0
         ! Pre-read file to get number of lines
         do
            read(ifileno,*,end=200)  isite
            flines=flines+1
         end do

         200 continue

         rewind(ifileno)

         write(*,'(2x,a)') 'Reading site dependent polarizations'

         ! If the size of the file is NATOM then there is no problem
         if ( Natom.eq.flines ) then

            do i=1, flines
               read(ifileno,*) isite, sitenatom_sot_pol(1,isite), sitenatom_sot_pol(2,isite), sitenatom_sot_pol(3,isite)
            end do
         else
            write(*,*) 'WARNING: Size of the sitenatom_sot_pol is not NATOM'
            do i=1, flines
               read(ifileno,*) isite, sitenatom_sot_pol(1,isite), sitenatom_sot_pol(2,isite), sitenatom_sot_pol(3,isite)
            end do

         end if

         close(ifileno)

         ! If there is no site dependent polarization set it to be constant
      else
         do i=1, Natom
            sitenatom_sot_pol(1,i)=sot_pol_vec(1)
            sitenatom_sot_pol(2,i)=sot_pol_vec(2)
            sitenatom_sot_pol(3,i)=sot_pol_vec(3)
         enddo
      endif

   end subroutine read_sot_pol_site
   !---------------------------------------------------------------------------
   !> @brief
   !> Resolve canonical current density from user inputs (single source of truth)
   !> @details This routine establishes the canonical current density representation
   !> used by all subsequent torque calculations. It implements clear precedence rules:
   !>
   !> PRECEDENCE (highest to lowest):
   !> 1. jdens (direct physical input in A/m²):
   !>    - Used directly as jdens_cell
   !>    - have_jdens = .true.
   !>    - If jvec also provided: warning issued, jvec ignored
   !>
   !> 2. jvec or jvecfile (legacy dimensionless input):
   !>    - Converted to jdens via: jdens = stt_dens_conv * jvec
   !>    - If jsite='Y': creates site-resolved jdens_site(3,Natom)
   !>    - If jsite='N': creates uniform jdens_cell(3)
   !>    - have_jdens = .false. (indicates conversion was used)
   !>    - Clear note printed about conversion factor and source
   !>
   !> 3. No current:
   !>    - jdens_cell = 0, all torques disabled
   !>    - Warning printed if STT/SHE are enabled
   !>
   !> OUTPUT:
   !> - jdens_cell(3): canonical uniform current (A/m²) - always set
   !> - jdens_site(3,Natom): canonical site-resolved (A/m²) - allocated if jsite='Y'
   !> - use_jdens_site: flag indicating which to use
   !> - Comprehensive printed summary with magnitudes, conversions, warnings
   !>
   !> REQUIREMENTS:
   !> Must be called after:
   !> - set_curr_density() (to compute stt_dens_conv)
   !> - read_jvecfile() (to populate sitenatomjvec if jsite='Y')
   !> Before:
   !> - Any call to calculate_spintorques()
   !>
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine resolve_current_density(Natom)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system

      integer :: i
      real(dblprec) :: jmag, jmag_min, jmag_max
      real(dblprec), parameter :: eps = 1.0e-15_dblprec

      ! Initialize canonical current state
      jdens_cell = 0.0_dblprec
      use_jdens_site = .false.
      if (allocated(jdens_site)) deallocate(jdens_site)

      write(*,'(1x,a)') '========================================'
      write(*,'(1x,a)') 'CURRENT DENSITY RESOLUTION'
      write(*,'(1x,a)') '========================================'

      ! Priority 1: User provided jdens directly (physical A/m^2)
      if (jdens_user_provided) then
         jdens_cell = jdens
         have_jdens = .true.
         write(*,'(1x,a)') 'Current source: jdens (user-provided physical current density)'
         write(*,'(1x,a,3(ES14.6,1x),a)') '  jdens (A/m^2): ', jdens_cell(1), jdens_cell(2), jdens_cell(3)
         write(*,'(1x,a,ES14.6,a)') '  |jdens|: ', norm2(jdens_cell), ' A/m^2'
         
         ! Check if both jdens and jvec were provided (warn about precedence)
         ! if (norm2(jvec) > eps) then
         !    write(*,'(1x,a)') 'WARNING: Both jdens and legacy jvec provided; using jdens and ignoring jvec'
         ! endif

      ! Priority 2: Legacy jvec or jvecfile (convert to physical jdens)
      else if (norm2(jvec) > eps .or. jsite == 'Y') then
         have_jdens = .false.
         
         ! Site-resolved legacy current
         if (jsite == 'Y') then
            ! sitenatomjvec already populated by read_jvecfile
            ! Convert to canonical site-resolved jdens
            allocate(jdens_site(3,Natom))
            jdens_site = 0.0_dblprec
            jmag_min = huge(1.0_dblprec)
            jmag_max = 0.0_dblprec
            
            do i=1, Natom
               jdens_site(:,i) = sitenatomjvec(:,i) * stt_dens_conv
               jmag = norm2(jdens_site(:,i))
               if (jmag > eps) then
                  jmag_min = min(jmag_min, jmag)
                  jmag_max = max(jmag_max, jmag)
               endif
            end do
            
            ! Set cell average for convenience
            jdens_cell = sum(jdens_site, dim=2) / real(Natom, dblprec)
            use_jdens_site = .true.
            
            write(*,'(1x,a)') 'Current source: legacy jvecfile (converted to physical jdens)'
            write(*,'(1x,a,ES14.6,a)') '  Conversion factor (stt_dens_conv): ', stt_dens_conv, ' A/m^2 per code unit'
            write(*,'(1x,a)') '  Site-resolved current: active'
            write(*,'(1x,a,ES14.6,a)') '  Min |jdens_site|: ', jmag_min, ' A/m^2'
            write(*,'(1x,a,ES14.6,a)') '  Max |jdens_site|: ', jmag_max, ' A/m^2'
            write(*,'(1x,a,3(ES14.6,1x),a)') '  Average jdens: ', jdens_cell(1), jdens_cell(2), jdens_cell(3), ' A/m^2'
            write(*,'(1x,a,ES14.6,a)') '  |<jdens>|: ', norm2(jdens_cell), ' A/m^2'
            
         ! Cell-uniform legacy current
         else
            jdens_cell = jvec * stt_dens_conv
            write(*,'(1x,a)') 'Current source: legacy jvec (converted to physical jdens)'
            write(*,'(1x,a,ES14.6,a)') '  Conversion factor (stt_dens_conv): ', stt_dens_conv, ' A/m^2 based on `alat`.'
            write(*,'(1x,a,3(ES14.6,1x),a)') '  jvec (code units): ', jvec(1), jvec(2), jvec(3)
            write(*,'(1x,a,3(ES14.6,1x),a)') '  jdens (A/m^2): ', jdens_cell(1), jdens_cell(2), jdens_cell(3)
            write(*,'(1x,a,ES14.6,a)') '  |jdens|: ', norm2(jdens_cell), ' A/m^2'
         endif

      ! No current provided
      else
         write(*,'(1x,a)') 'Current source: none'
         write(*,'(1x,a)') '  All current-driven torques will be zero'
         have_jdens = .false.
      endif

      ! Sanity check: warn about suspicious magnitudes
      jmag = norm2(jdens_cell)
      if (jmag > eps) then
         if (jmag < 1.0e6_dblprec) then
            write(*,'(1x,a)') 'WARNING: |jdens| < 10^6 A/m^2 - check units (seems low for STT)'
         else if (jmag > 1.0e14_dblprec) then
            write(*,'(1x,a)') 'WARNING: |jdens| > 10^14 A/m^2 - check units (seems very high)'
         endif
      endif

      write(*,'(1x,a)') '========================================'

   end subroutine resolve_current_density
   !---------------------------------------------------------------------------
   !> @brief
   !> Calculate conversion factors and prefactors for spin-transfer torques
   !> @details This routine computes fundamental quantities needed for STT/SHE:
   !>
   !> CONVERSION FACTOR:
   !> - stt_dens_conv: converts dimensionless jvec to physical jdens (A/m²)
   !>   Formula: stt_dens_conv = e * M_tot * gamma * a / (V_cell * P)
   !>   where M_tot = total magnetic moment per cell
   !>         a = lattice constant
   !>         V_cell = unit cell volume
   !>         P = spin_pol (spin polarization)
   !>
   !> TORQUE PREFACTORS:
   !> - b_rt_fac: Slonczewski STT prefactor (Meo 2023 formalism)
   !>   Includes: ℏ, spin polarization, area, depth, and physical constants
   !>   Applied per site in slonczewski_field() after dividing by local moment
   !>
   !> - sot_rt_fac: SHE torque prefactor
   !>   Includes: ℏ, spin Hall angle, area, depth, current magnitude
   !>   Applied per site in SHE_torque() after dividing by local moment
   !>
   !> SHE SIGMA DETERMINATION:
   !> Computes spin accumulation direction for SHE based on:
   !> - If she_n_vec provided: sigma = (n × j) / |n × j|  [proper interface physics]
   !> - If she_n_vec not provided: sigma = j / |j|        [legacy/testing mode]
   !> Guards against parallel n||j (undefined sigma, torque disabled)
   !>
   !> LEGACY COMPATIBILITY:
   !> Handles jvec <-> jdens conversion with clear precedence:
   !> - If both provided: jdens wins, jvec updated for consistency
   !> - If only jvec: converted to jdens
   !> - If only jdens: jvec updated (for any legacy code references)
   !>
   !> NOTE: This routine does NOT establish the canonical current density.
   !> That is done by resolve_current_density() which must be called after this.
   !>
   !> @author
   !> Jonathan Chico, Anders Bergman
   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine set_curr_density(NA,Natom,Nchmax,conf_num,alat,spin_pol,C1,C2,C3,jvec,ammom_inp)
      use InputData,   only : N3
      use Constants
      use math_functions, only : f_cross_product
      implicit none

      ! .. Input variables
      integer, intent(in) :: NA
      integer, intent(in) :: Natom
      integer, intent(in) :: Nchmax
      integer, intent(in) :: conf_num
      real(dblprec), intent(inout) :: alat
      real(dblprec), intent(inout) :: spin_pol
      real(dblprec), dimension(3), intent(in) :: C1, C2, C3
      real(dblprec), dimension(3), intent(inout) :: jvec
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp

      ! .. Local variables
      real(dblprec) :: cell_vol          ! [m^3]
      real(dblprec) :: V_at              ! [m^3] per atom in cell (magnetic atoms counted by NA)
      real(dblprec) :: total_mom         ! [mu_B] sum of moments in the reference cell
      real(dblprec) :: xy_area           ! [m^2]
      real(dblprec) :: jmag              ! [A/m^2]
      real(dblprec), dimension(3) :: jdir
      real(dblprec), dimension(3) :: n_norm, sigma_raw
      real(dblprec) :: detC              ! dimensionless determinant of C1,C2,C3
      real(dblprec) :: g_factor          ! dimensionless electron g-factor (derived)
      real(dblprec), parameter :: eps = 1.0e-15_dblprec

      !-----------------------------------------------------------------------
      ! 0) Defaults / sanity
      !-----------------------------------------------------------------------
      if (alat <= 0.0_dblprec .or. abs(alat-1.0_dblprec) < eps) then
         alat = 2.856e-10_dblprec
         write(*,'(1x,a,2x,ES14.6,1x,a)') 'No lattice constant given, assuming BCC Fe lattice constant:', alat, 'm'
      endif

      if (spin_pol == 0.0_dblprec) then
         spin_pol = 1.0_dblprec
         write(*,'(1x,a,2x,F10.6)') 'No polarization set, assuming 100%:', spin_pol
      endif

      !-----------------------------------------------------------------------
      ! 1) Cell volume and "atomic volume" per simulated atom (within NA)
      !    Here C1,C2,C3 are assumed dimensionless lattice vectors; alat sets meters.
      !-----------------------------------------------------------------------
      detC = (C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2) + &
              C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3) + &
              C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1))

      cell_vol = detC * alat**3         ! [m^3]
      V_at     = cell_vol / real(NA, dblprec)

      ! Total moment in the reference cell (user indicates ammom_inp(:,1,1) are magnitudes in mu_B)
      total_mom = sum(ammom_inp(:,1,1)) ! [mu_B]

      !-----------------------------------------------------------------------
      ! 2) Legacy jvec -> physical jdens conversion factor (kept as in your design)
      !
      ! NOTE: This conversion is code-convention dependent. Keeping your existing
      ! formula is fine if it reproduces legacy behavior. It is not used for new
      ! Zhang–Li physics; Zhang–Li uses physical jdens directly.
      !-----------------------------------------------------------------------
      stt_dens_conv = ev*total_mom*gama*alat/(cell_vol*spin_pol)

      !-----------------------------------------------------------------------
      ! 3) Interface area (for SHE/SOT-style "slab" models): area spanned by C1,C2
      !-----------------------------------------------------------------------
      xy_area = alat**2 * norm2(f_cross_product(C1,C2))   ! [m^2]

      !-----------------------------------------------------------------------
      ! 4) Slonczewski prefactor (Meo-style): ħ P /(2e) * (A_per_atom) / μB
      !    Downstream you typically multiply by |j| and divide by local mmom.
      !-----------------------------------------------------------------------
      b_rt_fac = hbar * spin_pol / (2.0_dblprec*ev) * xy_area / real(NA*N3, dblprec) / mub

      !-----------------------------------------------------------------------
      ! 5) Zhang–Li (Schieback / Zhang–Li) prefactor for lattice-unit gradients
      !
      ! Schieback: u = (P g μB)/(2 e Ms) * j
      ! Atomistic: Ms ≈ (μ/ V_at), μ = mmom*μB  => μB cancels:
      !     u_phys = (P g /(2e)) * (V_at/mmom) * j
      !
      ! Your coordinates are in lattice units, so the stencil computes ∂/∂(lattice-unit).
      ! Therefore we must use u_latt = u_phys / alat.
      !
      ! We build a factor b_zhang_li_fac such that:
      !     u_latt = (b_zhang_li_fac / mmom) * jdens
      !
      ! Formula (using electron g-factor g_e_abs from constants):
      !   b_zhang_li_fac = (P * g_e_abs / (2e)) * (V_at/alat)
      !
      ! NOTE: This prefactor does NOT include gamma. After torques are computed and stored
      ! in btorque(:,:,:), they will be divided by gamma (see below) to correct for the fact
      ! that depondt.f90 will multiply them by gamma during time evolution. This ensures
      ! consistent scaling across all field components (Hamiltonian beff, STT btorque, etc.).
      !
      ! Electron g-factor (dimensionless):
      !   g_e_abs = 2.00231930436182 (CODATA 2018 value)
      !-----------------------------------------------------------------------
      b_zhang_li_fac = spin_pol * g_e_abs / (2.0_dblprec*ev) * (V_at/alat)


      ! Units: (1/C) * (m^3/m) = m^2/C, so (b_zhang_li_fac * jdens) -> 1/s

      !-----------------------------------------------------------------------
      ! 6) Canonical jdens selection / legacy interoperability
      !-----------------------------------------------------------------------
      if (jdens_user_provided) then
         if (norm2(jvec) > eps) then
            write(*,'(1x,a)') 'NOTE: Both jvec and jdens are set; jdens takes precedence.'
         endif
         ! Keep legacy-consistency for any remaining code using jvec
         jvec = jdens / stt_dens_conv
      else if (norm2(jvec) > eps) then
         jdens = jvec * stt_dens_conv
      else
         jdens = 0.0_dblprec
      endif

      !-----------------------------------------------------------------------
      ! 7) Print summary
      !-----------------------------------------------------------------------
      write(*,'(1x,a,ES14.6,1x,a)') 'Cell volume:', cell_vol, 'm^3'
      write(*,'(1x,a,ES14.6,1x,a)') 'Atomic volume V_at:', V_at, 'm^3'
      write(*,'(1x,a,ES14.6)')      'Sum of moments in ref. cell (total_mom) [mu_B]:', total_mom
      write(*,'(1x,a,ES14.6,1x,a)') 'Legacy conversion (stt_dens_conv):', stt_dens_conv, '(A/m^2) per legacy jvec-unit'
      write(*,'(1x,a,ES14.6,1x,a)') 'Interface area (C1xC2):', xy_area, 'm^2'
      write(*,'(1x,a,ES14.6)')      'Slonczewski prefactor (b_rt_fac):', b_rt_fac
      write(*,'(1x,a,ES14.6,1x,a)') 'Zhang–Li lattice prefactor (b_zhang_li_fac):', b_zhang_li_fac, 'm^2/C'
      write(*,'(1x,a,3(ES14.6,1x),a,ES14.6)') 'jdens [A/m^2]:', jdens(1), jdens(2), jdens(3), ' |j|=', norm2(jdens)

      !-----------------------------------------------------------------------
      ! 8) SHE sigma determination (unchanged logic; uses physical jdens)
      !-----------------------------------------------------------------------
      jmag = norm2(jdens)

      if (jmag <= eps) then
         she_sigma_vec = 0.0_dblprec
         sot_rt_fac    = 0.0_dblprec
         if (do_she == 'Y') write(*,'(1x,a)') 'WARNING: SHE enabled but |jdens| ~ 0; SHE torque will be zero'
      else
         jdir = jdens / jmag

         if (norm2(she_n_vec) > eps) then
            n_norm = she_n_vec / norm2(she_n_vec)
            sigma_raw = f_cross_product(n_norm, jdir)

            if (norm2(sigma_raw) <= eps) then
               she_sigma_vec = 0.0_dblprec
               sot_rt_fac    = 0.0_dblprec
               if (do_she == 'Y') write(*,'(1x,a)') 'WARNING: SHE: she_n_vec parallel to jdens; sigma undefined; torque disabled'
            else
               she_sigma_vec = sigma_raw / norm2(sigma_raw)
               sot_rt_fac = hbar * SHE_angle / (2.0_dblprec*ev) * xy_area / real(NA*N3, dblprec) / mub * jmag
               if (do_she == 'Y') then
                  write(*,'(1x,a)') 'SHE sigma mode: sigma = n × j'
                  write(*,'(1x,a,3(G14.6,1x))') '  n (normalized): ', n_norm
                  write(*,'(1x,a,3(G14.6,1x))') '  j direction:    ', jdir
                  write(*,'(1x,a,3(G14.6,1x))') '  sigma:          ', she_sigma_vec
                  write(*,'(1x,a,ES14.6)') '  SHE strength (sot_rt_fac): ', sot_rt_fac
               endif
            endif
         else
            she_sigma_vec = jdir
            sot_rt_fac = hbar * SHE_angle / (2.0_dblprec*ev) * xy_area / real(NA*N3, dblprec) / mub * jmag
            if (do_she == 'Y') then
               write(*,'(1x,a)') 'SHE sigma mode: sigma = jdir (no she_n_vec provided)'
               write(*,'(1x,a,3(G14.6,1x))') '  sigma (= j direction): ', she_sigma_vec
               write(*,'(1x,a,ES14.6)') '  SHE strength (sot_rt_fac): ', sot_rt_fac
            endif
         endif
      endif

   end subroutine set_curr_density


end module SpinTorques
