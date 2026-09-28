!------------------------------------------------------------------------------------
!> @brief Trajectory orbital angular momentum in the pyswatter convention.
!>
!> The transverse field is `psi = m_x + i*m_y` in one global frame fixed at
!> initialization.  This is the Holstein--Primakoff annihilation convention;
!> positive lambda_L therefore denotes magnon OAM along +z.  It is the complex
!> conjugate of the m_x - i*m_y convention used in parts of the literature.
!------------------------------------------------------------------------------------
module orbital_angular_momentum

   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
   use Parameters
   use Profiling
   use InputData, only : Landeg_glob
   use Mesh2D, only : nsimp, simp, tri_area, site_area, site_wsum, grad_b, grad_c, &
      site_tri_ptr, site_tri_idx

   implicit none

   private

   character(len=1), public :: do_oam = 'N' !< Deprecated trajectory OAM alias
   character(len=1), public :: do_oam_traj = 'N' !< Enable trajectory OAM
   integer, public :: oam_step_traj = 100       !< Sampling interval
   integer, public :: oam_buff_traj = 10        !< Number of buffered rows
   real(dblprec), public :: oam_origin(3) = 0.0_dblprec !< Fixed origin
   logical, public :: oam_origin_set = .false.  !< Whether an origin was supplied
   character(len=4), public :: oam_weight = 'site' !< Site or area weighting
   real(dblprec), public :: oam_sigma_max = 0.6_dblprec !< Centroid spread guard
   real(dblprec), public :: oam_gfactor = 0.0_dblprec !< Optional g factor
   integer, allocatable, public :: oam_sublattice(:) !< Optional unit-cell sublattice list
   logical, public :: oam_sublattice_set = .false. !< Whether a sublattice list was supplied
   real(dblprec), public :: oam_lambda_centroid_sum = 0.0_dblprec
   integer, public :: oam_lambda_centroid_count = 0

   public :: oam_defaults, oam_init, oam_sample, oam_flush

   logical :: oam_active = .false.
   logical :: oam_initialized = .false.
   logical :: oam_header_written = .false.
   logical :: oam_norm_warning = .false.
   integer :: oam_natom = 0
   integer :: oam_mensemble = 0
   integer :: oam_n1 = 0
   integer :: oam_n2 = 0
   integer :: oam_na = 1
   integer :: oam_ncolumns = 9
   integer :: oam_nsubblocks = 0
   integer :: oam_rstep = 0
   integer :: oam_nbuffer = 0
   real(dblprec) :: oam_g = 2.0_dblprec
   real(dblprec) :: oam_frame(3,3) = 0.0_dblprec
   real(dblprec) :: oam_origin_xy(2) = 0.0_dblprec
   real(dblprec) :: oam_cell(2,2) = 0.0_dblprec
   real(dblprec) :: oam_inv_cell(2,2) = 0.0_dblprec
   real(dblprec) :: oam_half_short = 0.0_dblprec
   logical :: oam_periodic(2) = .false.
   character(len=8) :: oam_simid = ''

   real(dblprec), allocatable :: oam_coord(:,:)
   complex(dblprec), allocatable :: oam_psi(:)
   complex(dblprec), allocatable :: oam_fem_dx(:)
   complex(dblprec), allocatable :: oam_fem_dy(:)
   complex(dblprec), allocatable :: oam_site_dx(:)
   complex(dblprec), allocatable :: oam_site_dy(:)
   integer, allocatable :: oam_step_buffer(:)
   real(dblprec), allocatable :: oam_row_buffer(:,:)

contains

   !---------------------------------------------------------------------------------
   !> @brief Set defaults for trajectory OAM input.
   !---------------------------------------------------------------------------------
   subroutine oam_defaults()
      do_oam_traj = 'N'
      do_oam = 'N'
      oam_step_traj = 100
      oam_buff_traj = 10
      oam_origin = 0.0_dblprec
      oam_origin_set = .false.
      oam_weight = 'site'
      oam_sigma_max = 0.6_dblprec
      oam_gfactor = 0.0_dblprec
      oam_sublattice_set = .false.
      oam_lambda_centroid_sum = 0.0_dblprec
      oam_lambda_centroid_count = 0
      if (allocated(oam_sublattice)) deallocate(oam_sublattice)
   end subroutine oam_defaults

   !---------------------------------------------------------------------------------
   !> @brief Initialize the global frame and trajectory OAM work arrays.
   !---------------------------------------------------------------------------------
   subroutine oam_init(Natom,Mensemble,NA,N1,N2,coord,C1,C2,BC1,BC2,emom,simid,rstep)

      integer, intent(in) :: Natom, Mensemble, NA, N1, N2, rstep
      real(dblprec), intent(in) :: coord(3,Natom), C1(3), C2(3)
      character(len=1), intent(in) :: BC1, BC2
      real(dblprec), intent(in) :: emom(3,Natom,Mensemble)
      character(len=8), intent(in) :: simid

      integer :: i, j, k, i_stat
      real(dblprec) :: mean_m(3), mean_norm, alignment, det_cell, seed_dot
      real(dblprec) :: seed(3), ex(3), ey(3), ez(3)

      if (do_oam_traj /= 'Y') return
      if (Natom < 1 .or. Mensemble < 1 .or. NA < 1 .or. N1 < 1 .or. N2 < 1) then
         write(*,'(1x,a)') 'Trajectory OAM disabled: invalid system dimensions.'
         return
      end if
      if (oam_step_traj < 1) oam_step_traj = 1
      if (oam_buff_traj < 1) oam_buff_traj = 1
      if (trim(adjustl(oam_weight)) /= 'site' .and. trim(adjustl(oam_weight)) /= 'area') then
         write(*,'(1x,a,a)') 'Trajectory OAM: unknown oam_weight ',trim(oam_weight)
         oam_weight = 'site'
      end if
      if (oam_sublattice_set) then
         if (.not.allocated(oam_sublattice) .or. size(oam_sublattice)<1) then
            write(*,'(1x,a)') 'Trajectory OAM disabled: oam_sublattice is empty.'
            return
         end if
         do i=1,size(oam_sublattice)
            if (oam_sublattice(i)<1 .or. oam_sublattice(i)>NA) then
               write(*,'(1x,a,i0,a,i0)') 'Trajectory OAM disabled: oam_sublattice entry ', &
                  oam_sublattice(i), ' is outside 1..', NA
               return
            end if
            do j=1,i-1
               if (oam_sublattice(i)==oam_sublattice(j)) then
                  write(*,'(1x,a)') 'Trajectory OAM disabled: oam_sublattice contains duplicates.'
                  return
               end if
            end do
         end do
         oam_nsubblocks = 0
      else
         oam_nsubblocks = 0
         if (NA > 1) oam_nsubblocks = NA
      end if

      mean_m = 0.0_dblprec
      do k=1,Mensemble
         do i=1,Natom
            mean_m = mean_m + emom(:,i,k)
         end do
      end do
      mean_m = mean_m / real(Natom*Mensemble,dblprec)
      mean_norm = sqrt(dot_product(mean_m,mean_m))
      if (mean_norm <= 1.0e-14_dblprec) then
         write(*,'(1x,a)') 'Trajectory OAM disabled: the initial average moment is zero.'
         return
      end if
      ez = mean_m / mean_norm
      alignment = 1.0_dblprec
      do k=1,Mensemble
         do i=1,Natom
            alignment = min(alignment,dot_product(emom(:,i,k),ez))
         end do
      end do
      if (alignment < 0.9_dblprec) then
         write(*,'(1x,a,f10.6,a)') 'Trajectory OAM disabled: non-collinear initial state (minimum alignment ', &
            alignment, ' < 0.9).'
         return
      end if

      seed = (/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
      seed_dot = abs(dot_product(seed,ez))
      if (seed_dot >= 0.9_dblprec) seed = (/0.0_dblprec,1.0_dblprec,0.0_dblprec/)
      ex = seed - dot_product(seed,ez)*ez
      ex = ex / sqrt(dot_product(ex,ex))
      ey = cross_product(ez,ex)
      oam_frame(:,1) = ex
      oam_frame(:,2) = ey
      oam_frame(:,3) = ez

      oam_natom = Natom
      oam_mensemble = Mensemble
      oam_n1 = N1
      oam_n2 = N2
      oam_na = NA
      oam_ncolumns = 9 + 3*oam_nsubblocks
      oam_rstep = rstep
      oam_simid = simid
      oam_periodic = (/BC1=='P',BC2=='P'/)
      oam_origin_xy = 0.0_dblprec
      allocate(oam_coord(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(oam_coord))*kind(oam_coord),'oam_coord','oam_init')
      oam_coord = coord
      if (oam_origin_set) then
         oam_origin_xy = oam_origin(1:2)
      else
         oam_origin_xy(1) = sum(coord(1,:))/real(Natom,dblprec)
         oam_origin_xy(2) = sum(coord(2,:))/real(Natom,dblprec)
      end if

      oam_cell(:,1) = real(N1,dblprec)*C1(1:2)
      oam_cell(:,2) = real(N2,dblprec)*C2(1:2)
      det_cell = oam_cell(1,1)*oam_cell(2,2)-oam_cell(2,1)*oam_cell(1,2)
      if (abs(det_cell) <= 1.0e-14_dblprec) then
         write(*,'(1x,a)') 'Trajectory OAM disabled: the in-plane cell is singular.'
         call oam_release()
         return
      end if
      oam_inv_cell(1,1) = oam_cell(2,2)/det_cell
      oam_inv_cell(1,2) = -oam_cell(1,2)/det_cell
      oam_inv_cell(2,1) = -oam_cell(2,1)/det_cell
      oam_inv_cell(2,2) = oam_cell(1,1)/det_cell
      oam_half_short = 0.5_dblprec*min(sqrt(dot_product(oam_cell(:,1),oam_cell(:,1))), &
         sqrt(dot_product(oam_cell(:,2),oam_cell(:,2))))

      if (oam_gfactor > 0.0_dblprec) then
         oam_g = oam_gfactor
      else if (Landeg_glob > 0.0_dblprec) then
         oam_g = Landeg_glob
      else
         oam_g = 2.0_dblprec
      end if

      allocate(oam_psi(Natom),stat=i_stat)
      call memocc(i_stat,product(shape(oam_psi))*kind(oam_psi),'oam_psi','oam_init')
      allocate(oam_site_dx(Natom),stat=i_stat)
      call memocc(i_stat,product(shape(oam_site_dx))*kind(oam_site_dx),'oam_site_dx','oam_init')
      allocate(oam_site_dy(Natom),stat=i_stat)
      call memocc(i_stat,product(shape(oam_site_dy))*kind(oam_site_dy),'oam_site_dy','oam_init')
      allocate(oam_fem_dx(max(1,nsimp)),stat=i_stat)
      call memocc(i_stat,product(shape(oam_fem_dx))*kind(oam_fem_dx),'oam_fem_dx','oam_init')
      allocate(oam_fem_dy(max(1,nsimp)),stat=i_stat)
      call memocc(i_stat,product(shape(oam_fem_dy))*kind(oam_fem_dy),'oam_fem_dy','oam_init')
      allocate(oam_step_buffer(oam_buff_traj),stat=i_stat)
      call memocc(i_stat,product(shape(oam_step_buffer))*kind(oam_step_buffer),'oam_step_buffer','oam_init')
      allocate(oam_row_buffer(oam_ncolumns,oam_buff_traj),stat=i_stat)
      call memocc(i_stat,product(shape(oam_row_buffer))*kind(oam_row_buffer),'oam_row_buffer','oam_init')

      oam_nbuffer = 0
      oam_header_written = .false.
      oam_norm_warning = .false.
      oam_lambda_centroid_sum = 0.0_dblprec
      oam_lambda_centroid_count = 0
      oam_active = .true.
      oam_initialized = .true.
   end subroutine oam_init

   !---------------------------------------------------------------------------------
   !> @brief Sample the trajectory OAM and buffer one output row when scheduled.
   !---------------------------------------------------------------------------------
   subroutine oam_sample(mstep,emom,mmom)

      integer, intent(in) :: mstep
      real(dblprec), intent(in) :: emom(3,oam_natom,oam_mensemble)
      real(dblprec), intent(in) :: mmom(oam_natom,oam_mensemble)

      integer :: k, j
      integer :: nfinite(oam_ncolumns)
      real(dblprec) :: row(oam_ncolumns), row_sum(oam_ncolumns)

      if (.not.oam_active) return
      if (mod(mstep-oam_rstep-1,oam_step_traj) /= 0) return

      row_sum = 0.0_dblprec
      nfinite = 0
      do k=1,oam_mensemble
         call oam_evaluate_ensemble(k,emom,mmom,row)
         do j=1,9
            if (ieee_is_finite(row(j))) then
               row_sum(j) = row_sum(j) + row(j)
               nfinite(j) = nfinite(j) + 1
            end if
         end do
         do j=10,oam_ncolumns
            if (ieee_is_finite(row(j))) then
               row_sum(j) = row_sum(j) + row(j)
               nfinite(j) = nfinite(j) + 1
            end if
         end do
      end do
      row = ieee_value(0.0_dblprec,ieee_quiet_nan)
      do j=1,oam_ncolumns
         if (nfinite(j) > 0) row(j) = row_sum(j)/real(nfinite(j),dblprec)
      end do
      if (ieee_is_finite(row(2))) then
         oam_lambda_centroid_sum = oam_lambda_centroid_sum + row(2)
         oam_lambda_centroid_count = oam_lambda_centroid_count + 1
      end if
      call oam_buffer_row(mstep,row)
   end subroutine oam_sample

   !---------------------------------------------------------------------------------
   !> @brief Flush buffered trajectory OAM rows and release its work arrays.
   !---------------------------------------------------------------------------------
   subroutine oam_flush()
      if (.not.oam_initialized) return
      call oam_write_buffer()
      call oam_release()
      oam_active = .false.
      oam_initialized = .false.
   end subroutine oam_flush

   subroutine oam_evaluate_ensemble(k,emom,mmom,row)
      integer, intent(in) :: k
      real(dblprec), intent(in) :: emom(3,oam_natom,oam_mensemble)
      real(dblprec), intent(in) :: mmom(oam_natom,oam_mensemble)
      real(dblprec), intent(out) :: row(:)

      integer :: i, ip, itri, isub, offset
      real(dblprec) :: nm, lambda_origin, lambda_centroid, rx, ry, sigma
      logical :: centroid_valid, norm_valid

      do i=1,oam_natom
         oam_psi(i) = cmplx(dot_product(emom(:,i,k),oam_frame(:,1)), &
            dot_product(emom(:,i,k),oam_frame(:,2)),dblprec)
      end do

      do itri=1,nsimp
         oam_fem_dx(itri) = grad_b(1,itri)*oam_psi(simp(1,itri)) + &
            grad_b(2,itri)*oam_psi(simp(2,itri)) + grad_b(3,itri)*oam_psi(simp(3,itri))
         oam_fem_dy(itri) = grad_c(1,itri)*oam_psi(simp(1,itri)) + &
            grad_c(2,itri)*oam_psi(simp(2,itri)) + grad_c(3,itri)*oam_psi(simp(3,itri))
      end do

      oam_site_dx = (0.0_dblprec,0.0_dblprec)
      oam_site_dy = (0.0_dblprec,0.0_dblprec)
      !$omp parallel do default(shared) private(i,ip,itri) schedule(static)
      do i=1,oam_natom
         if (site_wsum(i) > 0.0_dblprec) then
            do ip=site_tri_ptr(i),site_tri_ptr(i+1)-1
               itri = site_tri_idx(ip)
               oam_site_dx(i) = oam_site_dx(i) + tri_area(itri)*oam_fem_dx(itri)
               oam_site_dy(i) = oam_site_dy(i) + tri_area(itri)*oam_fem_dy(itri)
            end do
            oam_site_dx(i) = oam_site_dx(i)/site_wsum(i)
            oam_site_dy(i) = oam_site_dy(i)/site_wsum(i)
         end if
      end do
      !$omp end parallel do

      row = ieee_value(0.0_dblprec,ieee_quiet_nan)
      call oam_group_metrics(k,emom,mmom,0,nm,lambda_origin,lambda_centroid,rx,ry,sigma, &
         centroid_valid,norm_valid)
      row(1) = lambda_origin
      row(2) = lambda_centroid
      row(3) = nm
      row(5) = nm
      row(7) = rx
      row(8) = ry
      row(9) = sigma
      if (norm_valid) then
         if (centroid_valid) then
            row(4) = nm*lambda_centroid
            row(6) = row(5)+row(4)
         end if
      end if

      do isub=1,oam_nsubblocks
         offset = 9 + 3*(isub-1)
         call oam_group_metrics(k,emom,mmom,isub,nm,lambda_origin,lambda_centroid,rx,ry,sigma, &
            centroid_valid,norm_valid)
         row(offset+1) = lambda_origin
         row(offset+2) = lambda_centroid
         row(offset+3) = nm
      end do
   end subroutine oam_evaluate_ensemble

   subroutine oam_group_metrics(k,emom,mmom,group,nm,lambda_origin,lambda_centroid,rx,ry,sigma, &
      centroid_valid,norm_valid)
      integer, intent(in) :: k, group
      real(dblprec), intent(in) :: emom(3,oam_natom,oam_mensemble)
      real(dblprec), intent(in) :: mmom(oam_natom,oam_mensemble)
      real(dblprec), intent(out) :: nm, lambda_origin, lambda_centroid, rx, ry, sigma
      logical, intent(out) :: centroid_valid, norm_valid

      integer :: i
      real(dblprec) :: weight, norm, psi_weight, angle
      real(dblprec) :: reduced(2), sin_sum(2), cos_sum(2), d(2), lever(2), lever_centroid(2)
      real(dblprec) :: ell, spread_limit

      nm = 0.0_dblprec
      norm = 0.0_dblprec
      do i=1,oam_natom
         if (.not.oam_site_selected(i,group)) cycle
         nm = nm + mmom(i,k)/oam_g*(1.0_dblprec-dot_product(emom(:,i,k),oam_frame(:,3)))
         if (site_wsum(i) > 0.0_dblprec) then
            weight = oam_site_weight(i)
            norm = norm + abs(oam_psi(i))**2*weight
         end if
      end do

      lambda_origin = ieee_value(0.0_dblprec,ieee_quiet_nan)
      lambda_centroid = ieee_value(0.0_dblprec,ieee_quiet_nan)
      rx = ieee_value(0.0_dblprec,ieee_quiet_nan)
      ry = ieee_value(0.0_dblprec,ieee_quiet_nan)
      sigma = ieee_value(0.0_dblprec,ieee_quiet_nan)
      centroid_valid = .false.
      norm_valid = norm >= 1.0e-14_dblprec
      if (.not.norm_valid) then
         if (group==0 .and. .not.oam_norm_warning) then
            write(*,'(1x,a)') 'WARNING: trajectory OAM norm is below 1e-14; writing NaN.'
            oam_norm_warning = .true.
         end if
         return
      end if

      sin_sum = 0.0_dblprec
      cos_sum = 0.0_dblprec
      rx = 0.0_dblprec
      ry = 0.0_dblprec
      do i=1,oam_natom
         if (.not.oam_site_selected(i,group)) cycle
         if (site_wsum(i) <= 0.0_dblprec) cycle
         weight = oam_site_weight(i)
         psi_weight = abs(oam_psi(i))**2*weight
         reduced = matmul(oam_inv_cell,oam_coord(1:2,i))
         if (is_periodic(1)) then
            angle = 2.0_dblprec*acos(-1.0_dblprec)*reduced(1)
            sin_sum(1) = sin_sum(1) + psi_weight*sin(angle)
            cos_sum(1) = cos_sum(1) + psi_weight*cos(angle)
         else
            rx = rx + psi_weight*reduced(1)
         end if
         if (is_periodic(2)) then
            angle = 2.0_dblprec*acos(-1.0_dblprec)*reduced(2)
            sin_sum(2) = sin_sum(2) + psi_weight*sin(angle)
            cos_sum(2) = cos_sum(2) + psi_weight*cos(angle)
         else
            ry = ry + psi_weight*reduced(2)
         end if
      end do
      do i=1,2
         if (is_periodic(i)) then
            angle = atan2(sin_sum(i),cos_sum(i))/(2.0_dblprec*acos(-1.0_dblprec))
            reduced(i) = modulo(angle,1.0_dblprec)
         else if (i==1) then
            reduced(i) = rx/norm
         else
            reduced(i) = ry/norm
         end if
      end do
      rx = oam_cell(1,1)*reduced(1)+oam_cell(1,2)*reduced(2)
      ry = oam_cell(2,1)*reduced(1)+oam_cell(2,2)*reduced(2)

      sigma = 0.0_dblprec
      lambda_origin = 0.0_dblprec
      lambda_centroid = 0.0_dblprec
      do i=1,oam_natom
         if (.not.oam_site_selected(i,group)) cycle
         if (site_wsum(i) <= 0.0_dblprec) cycle
         weight = oam_site_weight(i)
         psi_weight = abs(oam_psi(i))**2*weight
         lever = oam_coord(1:2,i)-oam_origin_xy
         ell = aimag(conjg(oam_psi(i))*(lever(1)*oam_site_dy(i)-lever(2)*oam_site_dx(i)))
         lambda_origin = lambda_origin + ell*weight
         d = oam_coord(1:2,i)-(/rx,ry/)
         reduced = matmul(oam_inv_cell,d)
         if (is_periodic(1)) reduced(1)=reduced(1)-real(nint(reduced(1)),dblprec)
         if (is_periodic(2)) reduced(2)=reduced(2)-real(nint(reduced(2)),dblprec)
         lever_centroid = matmul(oam_cell,reduced)
         ell = aimag(conjg(oam_psi(i))*(lever_centroid(1)*oam_site_dy(i)- &
            lever_centroid(2)*oam_site_dx(i)))
         lambda_centroid = lambda_centroid + ell*weight
         sigma = sigma + psi_weight*dot_product(lever_centroid,lever_centroid)
      end do
      lambda_origin = lambda_origin/norm
      lambda_centroid = lambda_centroid/norm
      sigma = sqrt(sigma/norm)
      spread_limit = oam_sigma_max*oam_half_short
      centroid_valid = sigma <= spread_limit
      if (.not.centroid_valid) lambda_centroid = ieee_value(0.0_dblprec,ieee_quiet_nan)
   end subroutine oam_group_metrics

   logical function oam_site_selected(i,group)
      integer, intent(in) :: i, group
      integer :: isub, j

      isub = mod(i-1,oam_na)+1
      if (group > 0) then
         oam_site_selected = isub == group
      else if (.not.oam_sublattice_set) then
         oam_site_selected = .true.
      else
         oam_site_selected = .false.
         do j=1,size(oam_sublattice)
            if (isub == oam_sublattice(j)) then
               oam_site_selected = .true.
               exit
            end if
         end do
      end if
   end function oam_site_selected

   subroutine oam_buffer_row(mstep,row)
      integer, intent(in) :: mstep
      real(dblprec), intent(in) :: row(:)
      oam_nbuffer = oam_nbuffer + 1
      oam_step_buffer(oam_nbuffer) = mstep
      oam_row_buffer(:,oam_nbuffer) = row
      if (.not.oam_header_written) call oam_write_header()
      if (oam_nbuffer >= oam_buff_traj) call oam_write_buffer()
   end subroutine oam_buffer_row

   subroutine oam_write_header()
      character(len=1024) :: filn, columns
      integer :: ios, isub

      write(filn,'("oam_traj.",a,".out")') trim(oam_simid)
      open(unit=ofileno,file=trim(filn),status='replace',action='write',iostat=ios)
      if (ios /= 0) error stop 'Trajectory OAM: unable to open output file'
      write(ofileno,'(a)') '# psi = m_x + i*m_y in the global frame fixed at oam_init; lambda_L > 0 means magnon OAM along +z.'
      write(ofileno,'(a,a)') '# oam_weight = ',trim(oam_weight)
      write(ofileno,'(a,es24.16)') '# g = ',oam_g
      write(ofileno,'(a,3(es24.16,1x))') '# origin = ',oam_origin_xy(1),oam_origin_xy(2),oam_origin(3)
      write(ofileno,'(a)') '# lambda_L_centroid is referenced to the |psi|^2 centroid R; this removes the drift term ' // &
         '(R x P)_z but not the envelope winding l.'
      columns = '# step lambda_L_origin lambda_L_centroid N_m Lz_tot_hbar dSz_hbar balance R_x R_y sigma_psi'
      do isub=1,oam_nsubblocks
         write(columns(len_trim(columns)+1:),'(a,i0,a,i0,a,i0)') ' lambda_L_origin_s',isub, &
            ' lambda_L_centroid_s',isub,' N_m_s',isub
      end do
      write(ofileno,'(a)') trim(columns)
      close(ofileno)
      oam_header_written = .true.
   end subroutine oam_write_header

   subroutine oam_write_buffer()
      character(len=256) :: filn, fmt
      integer :: i, ios

      if (oam_nbuffer == 0 .or. .not.oam_header_written) return
      write(filn,'("oam_traj.",a,".out")') trim(oam_simid)
      open(unit=ofileno,file=trim(filn),position='append',action='write',iostat=ios)
      if (ios /= 0) error stop 'Trajectory OAM: unable to append output file'
      write(fmt,'("(i8,1x,",i0,"(es24.16,1x))")') oam_ncolumns
      do i=1,oam_nbuffer
         write(ofileno,fmt) oam_step_buffer(i),oam_row_buffer(:,i)
      end do
      close(ofileno)
      oam_nbuffer = 0
   end subroutine oam_write_buffer

   subroutine oam_release()
      integer :: i_stat, i_all

      if (allocated(oam_coord)) then
         i_all=-product(shape(oam_coord))*kind(oam_coord)
         deallocate(oam_coord,stat=i_stat)
         call memocc(i_stat,i_all,'oam_coord','oam_release')
      end if
      if (allocated(oam_psi)) then
         i_all=-product(shape(oam_psi))*kind(oam_psi)
         deallocate(oam_psi,stat=i_stat)
         call memocc(i_stat,i_all,'oam_psi','oam_release')
      end if
      if (allocated(oam_fem_dx)) then
         i_all=-product(shape(oam_fem_dx))*kind(oam_fem_dx)
         deallocate(oam_fem_dx,stat=i_stat)
         call memocc(i_stat,i_all,'oam_fem_dx','oam_release')
      end if
      if (allocated(oam_fem_dy)) then
         i_all=-product(shape(oam_fem_dy))*kind(oam_fem_dy)
         deallocate(oam_fem_dy,stat=i_stat)
         call memocc(i_stat,i_all,'oam_fem_dy','oam_release')
      end if
      if (allocated(oam_site_dx)) then
         i_all=-product(shape(oam_site_dx))*kind(oam_site_dx)
         deallocate(oam_site_dx,stat=i_stat)
         call memocc(i_stat,i_all,'oam_site_dx','oam_release')
      end if
      if (allocated(oam_site_dy)) then
         i_all=-product(shape(oam_site_dy))*kind(oam_site_dy)
         deallocate(oam_site_dy,stat=i_stat)
         call memocc(i_stat,i_all,'oam_site_dy','oam_release')
      end if
      if (allocated(oam_step_buffer)) then
         i_all=-product(shape(oam_step_buffer))*kind(oam_step_buffer)
         deallocate(oam_step_buffer,stat=i_stat)
         call memocc(i_stat,i_all,'oam_step_buffer','oam_release')
      end if
      if (allocated(oam_row_buffer)) then
         i_all=-product(shape(oam_row_buffer))*kind(oam_row_buffer)
         deallocate(oam_row_buffer,stat=i_stat)
         call memocc(i_stat,i_all,'oam_row_buffer','oam_release')
      end if
   end subroutine oam_release

   real(dblprec) function oam_site_weight(i) result(weight)
      integer, intent(in) :: i
      if (trim(adjustl(oam_weight)) == 'area') then
         weight = site_area(i)
      else
         weight = 1.0_dblprec
      end if
   end function oam_site_weight

   logical function is_periodic(axis)
      integer, intent(in) :: axis
      if (axis < 1 .or. axis > 2) then
         is_periodic = .false.
      else
         is_periodic = oam_periodic(axis)
      end if
   end function is_periodic

   function cross_product(a,b) result(c)
      real(dblprec), intent(in) :: a(3), b(3)
      real(dblprec) :: c(3)
      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)
   end function cross_product

end module orbital_angular_momentum
