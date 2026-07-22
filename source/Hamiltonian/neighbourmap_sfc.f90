!-------------------------------------------------------------------------------
! MODULE: NeighbourMapSFC
!> @brief General cell-linked neighbour mapping for UppASD.
!>
!> @details
!> This module provides setup_nm_sfc as a reliable optional replacement for
!> setup_nm.  The historical module/routine name is retained for compatibility,
!> but the fragile Morton-index-window heuristic has been replaced by a
!> deterministic spatial cell index in supercell fractional coordinates.
!>
!> Main properties:
!>   * Supports arbitrary (including non-orthogonal) lattice vectors.
!>   * Supports periodic, open, and mixed boundary conditions.
!>   * Supports interactions with periodic images of the same atom.
!>   * Uses symmetry-expanded neighbour vectors exactly as setup_nm does.
!>   * Optionally enforces target atom types through nntype.
!>   * Requires exactly one matching atom for every generated neighbour.
!>   * Fails explicitly on missing or ambiguous mappings.
!>   * Expected O(N + N_interactions) scaling at fixed density.
!>   * OpenMP parallelization over central atoms.
!>
!> Random-alloy/dilute indexing is intentionally rejected in this implementation.
!> The ordinary setup_nm path should be used when do_ralloy == 1.
!-------------------------------------------------------------------------------
module NeighbourMapSFC
   use Parameters
   use Profiling
   !$ use omp_lib

   implicit none
   private

   public :: setup_nm_sfc

   integer, dimension(:), allocatable :: nsym
   real(dblprec), dimension(:,:,:,:), allocatable :: sym_mats
   integer, dimension(:,:), allocatable :: nmdimt

   ! Mapping and symmetry-vector tolerances are intentionally independent.
   ! Both are squared only at the point of comparison.
   real(dblprec), parameter :: map_tol_default = 1.0d-5
   real(dblprec), parameter :: sym_tol_default = 1.0d-8

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_nm_sfc
   !> Construct a full neighbour map using a periodic cell-linked spatial index.
   !----------------------------------------------------------------------------
   subroutine setup_nm_sfc(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, &
      atype, coord, max_no_neigh, max_no_shells, max_no_equiv, sym, nn, redcoord,   &
      nm, nmdim, do_ralloy, Natom_full, acellnumb, atype_ch, nntype)

      implicit none

      integer, intent(in) :: NT, NA, Natom
      integer, intent(in) :: N1, N2, N3
      integer, intent(in) :: Natom_full
      integer, intent(in) :: max_no_shells
      integer, intent(in) :: sym
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: atype
      real(dblprec), dimension(3), intent(in) :: C1, C2, C3
      character(len=1), intent(in) :: BC1, BC2, BC3
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NT,max_no_shells,3), intent(in) :: redcoord

      integer, intent(out) :: max_no_neigh
      integer, intent(out) :: max_no_equiv
      integer, dimension(:,:,:), allocatable, intent(out) :: nm
      integer, dimension(:,:), allocatable, intent(out) :: nmdim

      integer, intent(in), optional :: do_ralloy
      integer, dimension(Natom_full), optional, intent(in) :: acellnumb
      integer, dimension(Natom_full), optional, intent(in) :: atype_ch
      integer, dimension(NT,max_no_shells), optional, intent(in) :: nntype

      integer :: ralloy
      integer :: i_stat, i_all
      integer :: nelem
      integer :: iat, jat, itype, jtype, ishell, inei
      integer :: counter
      integer :: num_threads

      real(dblprec), dimension(3) :: S1, S2, S3
      real(dblprec), dimension(3,3) :: inv_super
      real(dblprec) :: det_super
      logical, dimension(3) :: periodic

      real(dblprec), dimension(:,:), allocatable :: frac
      real(dblprec), dimension(3) :: frac_min, frac_max, frac_range
      real(dblprec), dimension(3) :: atom_frac_min, atom_frac_max
      real(dblprec), dimension(3) :: recip_norm
      integer, dimension(3) :: nbin
      integer :: ncells
      integer, dimension(:), allocatable :: head, next_atom

      real(dblprec), dimension(:,:,:,:), allocatable :: nncoord

      real(dblprec), dimension(3) :: target_vec, target_pos
      real(dblprec), dimension(3) :: target_frac, target_frac_raw
      real(dblprec), dimension(3) :: dcart, dfrac, residual
      real(dblprec) :: res2, best_res2
      integer :: best_jat, nmatch
      integer, dimension(3) :: ibase, ib
      integer :: dx, dy, dz, cid, nbid
      integer, dimension(27) :: bin_ids
      integer :: n_bin_ids, ibin
      logical :: typematch
      logical :: onsite_target, outside_open

      integer :: fail_count
      integer :: fail_iat, fail_shell, fail_inei, fail_nmatch
      real(dblprec) :: fail_best_res
      real(dblprec), dimension(3) :: fail_target
      integer :: ambiguous_count

      real(dblprec), parameter :: map_tol  = map_tol_default
      real(dblprec), parameter :: map_tol2 = map_tol_default*map_tol_default

      write(*,'(2x,a)') 'Setting up cell-linked SFC neighbour mapping...'

      ! Validate the mode before allocating module state.
      ralloy = 0
      if (present(do_ralloy)) ralloy = do_ralloy

      if (ralloy /= 0 .and. ralloy /= 1) then
         error stop 'setup_nm_sfc: do_ralloy must be 0 or 1'
      end if

      if (ralloy == 1) then
         write(*,'(2x,a)') &
            'setup_nm_sfc: random-alloy/dilute mapping is not supported.'
         write(*,'(2x,a)') &
            'Use the ordinary setup_nm routine for do_ralloy = 1.'
         error stop 'setup_nm_sfc: unsupported random-alloy workflow'
      end if

      if (Natom <= 0 .or. NA <= 0 .or. NT <= 0) then
         error stop 'setup_nm_sfc: invalid atom/type dimensions'
      end if
      if (N1 <= 0 .or. N2 <= 0 .or. N3 <= 0) then
         error stop 'setup_nm_sfc: supercell repetitions must be positive'
      end if
      if (maxval(atype) > NT .or. minval(atype) < 1) then
         error stop 'setup_nm_sfc: atype contains an invalid type index'
      end if
      if (maxval(nn) > max_no_shells .or. minval(nn) < 0) then
         error stop 'setup_nm_sfc: nn is inconsistent with max_no_shells'
      end if

      periodic(1) = (BC1 == 'P')
      periodic(2) = (BC2 == 'P')
      periodic(3) = (BC3 == 'P')

      nelem = 1

      !$ num_threads = omp_get_max_threads()
      !$ write(*,'(4x,a,i0,a)') 'OpenMP enabled with ', num_threads, ' threads'

      ! Build the simulation-supercell vectors.  Periodic wrapping must use the
      ! full simulated box, not the primitive-cell vectors.
      S1 = real(N1,dblprec)*C1
      S2 = real(N2,dblprec)*C2
      S3 = real(N3,dblprec)*C3

      call invert_lattice(S1,S2,S3,inv_super,det_super)

      ! Symmetry setup and expansion.
      allocate(nsym(NT),stat=i_stat)
      call memocc(i_stat,product(shape(nsym))*kind(nsym),'nsym','setup_nm_sfc')
      if (i_stat /= 0) error stop 'setup_nm_sfc: allocation failure for nsym'
      nsym = 0

      allocate(nmdimt(maxval(nn),NT),stat=i_stat)
      call memocc(i_stat,product(shape(nmdimt))*kind(nmdimt),'nmdimt','setup_nm_sfc')
      if (i_stat /= 0) error stop 'setup_nm_sfc: allocation failure for nmdimt'
      nmdimt = 0

      call get_symops(sym,NT)

      max_no_equiv = maxval(nsym)
      if (max_no_equiv <= 0) then
         error stop 'setup_nm_sfc: no symmetry operations generated'
      end if

      allocate(nncoord(3,max_no_equiv,max_no_shells,NT),stat=i_stat)
      call memocc(i_stat,product(shape(nncoord))*kind(nncoord), &
                  'nncoord','setup_nm_sfc')
      if (i_stat /= 0) error stop 'setup_nm_sfc: allocation failure for nncoord'
      nncoord = 0.0_dblprec

      call get_fullnnlist(NT,nn,nelem,max_no_shells,redcoord, &
                          max_no_equiv,nncoord)

      allocate(nm(Natom,maxval(nn),max_no_equiv),stat=i_stat)
      call memocc(i_stat,product(shape(nm))*kind(nm),'nm','setup_nm_sfc')
      if (i_stat /= 0) error stop 'setup_nm_sfc: allocation failure for nm'
      nm = 0

      allocate(nmdim(maxval(nn),Natom),stat=i_stat)
      call memocc(i_stat,product(shape(nmdim))*kind(nmdim), &
                  'nmdim','setup_nm_sfc')
      if (i_stat /= 0) error stop 'setup_nm_sfc: allocation failure for nmdim'
      nmdim = 0

      ! Convert all explicit atom coordinates to supercell fractional
      ! coordinates. Periodic components are canonicalized to [0,1).
      allocate(frac(3,Natom),stat=i_stat)
      if (i_stat /= 0) error stop 'setup_nm_sfc: allocation failure for frac'

      do iat = 1, Natom
         call cart_to_frac(coord(:,iat),inv_super,frac(:,iat))
         call canonicalize_fractional(frac(:,iat),periodic)
      end do

      ! Define the fractional domain used by the spatial bins.
      atom_frac_min = minval(frac,dim=2)
      atom_frac_max = maxval(frac,dim=2)
      frac_min = atom_frac_min
      frac_max = atom_frac_max

      do cid = 1, 3
         if (periodic(cid)) then
            frac_min(cid) = 0.0_dblprec
            frac_max(cid) = 1.0_dblprec
         else
            if (frac_max(cid)-frac_min(cid) < map_tol) then
               frac_min(cid) = frac_min(cid)-0.5_dblprec
               frac_max(cid) = frac_max(cid)+0.5_dblprec
            else
               frac_min(cid) = frac_min(cid)-2.0_dblprec*map_tol
               frac_max(cid) = frac_max(cid)+2.0_dblprec*map_tol
            end if
         end if
      end do
      frac_range = frac_max-frac_min

      ! recip_norm(k) bounds the change in fractional component k caused by a
      ! Cartesian displacement of length map_tol.
      recip_norm(1) = sqrt(sum(inv_super(:,1)**2))
      recip_norm(2) = sqrt(sum(inv_super(:,2)**2))
      recip_norm(3) = sqrt(sum(inv_super(:,3)**2))

      call choose_bins(Natom,frac_range,recip_norm,map_tol,nbin)
      ncells = nbin(1)*nbin(2)*nbin(3)

      allocate(head(ncells),stat=i_stat)
      if (i_stat /= 0) error stop 'setup_nm_sfc: allocation failure for head'
      allocate(next_atom(Natom),stat=i_stat)
      if (i_stat /= 0) error stop 'setup_nm_sfc: allocation failure for next_atom'

      head = 0
      next_atom = 0

      ! Build linked lists for all spatial cells.
      do iat = 1, Natom
         call fractional_bin(frac(:,iat),frac_min,frac_range,nbin,periodic,ib)
         cid = linear_bin(ib,nbin)
         next_atom(iat) = head(cid)
         head(cid) = iat
      end do

      write(*,'(4x,a,3(i0,1x),a,i0)') &
         'Spatial bins: ',nbin(1),nbin(2),nbin(3),' total=',ncells
      write(*,'(4x,a,es12.4)') 'Cartesian mapping tolerance: ',map_tol

      fail_count = 0
      fail_iat = 0
      fail_shell = 0
      fail_inei = 0
      fail_nmatch = 0
      fail_best_res = huge(1.0_dblprec)
      fail_target = 0.0_dblprec
      ambiguous_count = 0

      ! Each central atom owns a disjoint nm/nmdim slice, so this loop is safe
      ! to parallelize.  The linked spatial index is read-only.
      !$omp parallel do default(shared) schedule(dynamic,16)                         &
      !$omp private(iat,itype,ishell,inei,counter,target_vec,target_pos,target_frac, &
      !$omp         target_frac_raw,ibase,bin_ids,n_bin_ids,dx,dy,dz,ib,nbid,ibin,  &
      !$omp         jat,jtype,typematch,onsite_target,outside_open,dcart,dfrac,      &
      !$omp         residual,res2,best_res2,                                         &
      !$omp         best_jat,nmatch)
      do iat = 1, Natom
         itype = atype(iat)

         do ishell = 1, nn(itype)
            counter = 0

            do inei = 1, nmdimt(ishell,itype)
               target_vec = nncoord(:,inei,ishell,itype)
               target_pos = coord(:,iat)+target_vec
               call cart_to_frac(target_pos,inv_super,target_frac_raw)
               target_frac = target_frac_raw
               call canonicalize_fractional(target_frac,periodic)

               ! In an open direction, a symmetry-generated neighbour may lie
               ! outside the explicitly simulated box.  This is a legitimate
               ! boundary omission and is represented by nm=0, as in setup_nm.
               outside_open = .false.
               if (.not.periodic(1)) outside_open = outside_open .or. &
                  target_frac_raw(1) < atom_frac_min(1)-map_tol*recip_norm(1) .or. &
                  target_frac_raw(1) > atom_frac_max(1)+map_tol*recip_norm(1)
               if (.not.periodic(2)) outside_open = outside_open .or. &
                  target_frac_raw(2) < atom_frac_min(2)-map_tol*recip_norm(2) .or. &
                  target_frac_raw(2) > atom_frac_max(2)+map_tol*recip_norm(2)
               if (.not.periodic(3)) outside_open = outside_open .or. &
                  target_frac_raw(3) < atom_frac_min(3)-map_tol*recip_norm(3) .or. &
                  target_frac_raw(3) > atom_frac_max(3)+map_tol*recip_norm(3)

               call fractional_bin(target_frac,frac_min,frac_range,nbin, &
                                   periodic,ibase)

               ! Search the target bin and its immediate neighbours. choose_bins
               ! guarantees that one bin width is at least twice the maximum
               ! fractional displacement associated with map_tol.
               n_bin_ids = 0
               do dz = -1, 1
                  do dy = -1, 1
                     do dx = -1, 1
                        ib = ibase + [dx,dy,dz]
                        if (.not. valid_or_wrap_bin(ib,nbin,periodic)) cycle
                        nbid = linear_bin(ib,nbin)
                        call append_unique_bin(nbid,bin_ids,n_bin_ids)
                     end do
                  end do
               end do

               best_res2 = huge(1.0_dblprec)
               best_jat = 0
               nmatch = 0
               onsite_target = (dot_product(target_vec,target_vec) < map_tol2)

               do ibin = 1, n_bin_ids
                  jat = head(bin_ids(ibin))
                  do while (jat /= 0)
                     ! Exclude only the true zero-displacement on-site target.
                     ! A translated periodic image of the same atom is valid.
                     if (onsite_target .and. jat == iat) then
                        jat = next_atom(jat)
                        cycle
                     end if

                     jtype = atype(jat)
                     typematch = .true.
                     if (present(nntype)) then
                        typematch = (jtype == nntype(itype,ishell))
                     end if

                     if (typematch) then
                        dcart = coord(:,jat)-target_pos
                        call cart_to_frac(dcart,inv_super,dfrac)

                        if (periodic(1)) dfrac(1) = dfrac(1)-anint(dfrac(1))
                        if (periodic(2)) dfrac(2) = dfrac(2)-anint(dfrac(2))
                        if (periodic(3)) dfrac(3) = dfrac(3)-anint(dfrac(3))

                        residual = dfrac(1)*S1+dfrac(2)*S2+dfrac(3)*S3
                        res2 = dot_product(residual,residual)

                        if (res2 < map_tol2) then
                           nmatch = nmatch+1
                           if (res2 < best_res2) then
                              best_res2 = res2
                              best_jat = jat
                           end if
                        end if
                     end if

                     jat = next_atom(jat)
                  end do
               end do

               if (nmatch == 1) then
                  counter = counter+1
                  nm(iat,ishell,inei) = best_jat
               else if (nmatch == 0 .and. outside_open) then
                  ! Legitimate missing neighbour beyond an open boundary.
                  nm(iat,ishell,inei) = 0
               else
                  !$omp critical(setup_nm_sfc_failure)
                  fail_count = fail_count+1
                  if (fail_iat == 0) then
                     fail_iat = iat
                     fail_shell = ishell
                     fail_inei = inei
                     fail_nmatch = nmatch
                     fail_best_res = sqrt(best_res2)
                     fail_target = target_pos
                  end if
                  if (nmatch > 1) ambiguous_count = ambiguous_count+1
                  !$omp end critical(setup_nm_sfc_failure)
               end if
            end do

            ! Match setup_nm semantics: nmdim is the full symmetry-expanded
            ! shell dimension. Entries outside open boundaries remain zero.
            nmdim(ishell,iat) = nmdimt(ishell,itype)
         end do
      end do
      !$omp end parallel do

      if (fail_count > 0) then
         write(*,'(2x,a,i0)') 'setup_nm_sfc: failed mappings = ',fail_count
         write(*,'(2x,a,i0)') 'setup_nm_sfc: ambiguous mappings = ',ambiguous_count
         write(*,'(2x,a,3i10)') 'First failure iat, shell, inei = ', &
                                fail_iat,fail_shell,fail_inei
         write(*,'(2x,a,i10)') 'Number of matching atoms = ',fail_nmatch
         write(*,'(2x,a,3es24.14)') 'Target Cartesian position = ',fail_target
         if (fail_best_res >= 0.0_dblprec .and. &
             fail_best_res < huge(1.0_dblprec)/2.0_dblprec) then
            write(*,'(2x,a,es24.14)') 'Best residual = ',fail_best_res
         end if
         error stop 'setup_nm_sfc: incomplete or ambiguous neighbour mapping'
      end if

      ! Calculate the actual maximum total neighbour count per atom.
      max_no_neigh = 1
      do iat = 1, Natom
         counter = 0
         do ishell = 1, nn(atype(iat))
            counter = counter+nmdim(ishell,iat)
         end do
         max_no_neigh = max(max_no_neigh,counter)
      end do

      write(*,'(4x,a,i0)') 'Maximum neighbours per atom: ',max_no_neigh
      write(*,'(4x,a,i0)') 'Maximum equivalent neighbours: ',max_no_equiv
      write(*,'(2x,a)') 'Cell-linked SFC neighbour mapping completed successfully.'

      ! Temporary storage cleanup. Output arrays nm and nmdim remain allocated.
      deallocate(head,stat=i_stat)
      deallocate(next_atom,stat=i_stat)
      deallocate(frac,stat=i_stat)

      i_all = -product(shape(nncoord))*kind(nncoord)
      deallocate(nncoord,stat=i_stat)
      call memocc(i_stat,i_all,'nncoord','setup_nm_sfc')

      i_all = -product(shape(nmdimt))*kind(nmdimt)
      deallocate(nmdimt,stat=i_stat)
      call memocc(i_stat,i_all,'nmdimt','setup_nm_sfc')

      i_all = -product(shape(nsym))*kind(nsym)
      deallocate(nsym,stat=i_stat)
      call memocc(i_stat,i_all,'nsym','setup_nm_sfc')

      i_all = -product(shape(sym_mats))*kind(sym_mats)
      deallocate(sym_mats,stat=i_stat)
      call memocc(i_stat,i_all,'sym_mats','setup_nm_sfc')

   end subroutine setup_nm_sfc


   !----------------------------------------------------------------------------
   ! SUBROUTINE: invert_lattice
   !> Invert a lattice-vector matrix using the convention employed by setup_nm.
   !----------------------------------------------------------------------------
   subroutine invert_lattice(C1,C2,C3,invmatrix,detmatrix)
      implicit none
      real(dblprec), dimension(3), intent(in) :: C1,C2,C3
      real(dblprec), dimension(3,3), intent(out) :: invmatrix
      real(dblprec), intent(out) :: detmatrix

      detmatrix = C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+ &
                  C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+ &
                  C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1)

      if (abs(detmatrix) <= dbl_tolerance) then
         error stop 'NeighbourMapSFC: singular lattice matrix'
      end if

      invmatrix(1,1)=(C2(2)*C3(3)-C3(2)*C2(3))/detmatrix
      invmatrix(1,2)=(C1(3)*C3(2)-C3(3)*C1(2))/detmatrix
      invmatrix(1,3)=(C1(2)*C2(3)-C2(2)*C1(3))/detmatrix
      invmatrix(2,1)=(C2(3)*C3(1)-C3(3)*C2(1))/detmatrix
      invmatrix(2,2)=(C1(1)*C3(3)-C3(1)*C1(3))/detmatrix
      invmatrix(2,3)=(C1(3)*C2(1)-C2(3)*C1(1))/detmatrix
      invmatrix(3,1)=(C2(1)*C3(2)-C3(1)*C2(2))/detmatrix
      invmatrix(3,2)=(C1(2)*C3(1)-C3(2)*C1(1))/detmatrix
      invmatrix(3,3)=(C1(1)*C2(2)-C2(1)*C1(2))/detmatrix
   end subroutine invert_lattice


   !----------------------------------------------------------------------------
   ! SUBROUTINE: cart_to_frac
   !----------------------------------------------------------------------------
   pure subroutine cart_to_frac(cart,invmatrix,frac)
      implicit none
      real(dblprec), dimension(3), intent(in) :: cart
      real(dblprec), dimension(3,3), intent(in) :: invmatrix
      real(dblprec), dimension(3), intent(out) :: frac

      frac(1)=cart(1)*invmatrix(1,1)+cart(2)*invmatrix(2,1)+ &
              cart(3)*invmatrix(3,1)
      frac(2)=cart(1)*invmatrix(1,2)+cart(2)*invmatrix(2,2)+ &
              cart(3)*invmatrix(3,2)
      frac(3)=cart(1)*invmatrix(1,3)+cart(2)*invmatrix(2,3)+ &
              cart(3)*invmatrix(3,3)
   end subroutine cart_to_frac


   !----------------------------------------------------------------------------
   ! SUBROUTINE: canonicalize_fractional
   !----------------------------------------------------------------------------
   pure subroutine canonicalize_fractional(frac,periodic)
      implicit none
      real(dblprec), dimension(3), intent(inout) :: frac
      logical, dimension(3), intent(in) :: periodic
      integer :: k

      do k=1,3
         if (periodic(k)) then
            frac(k)=frac(k)-floor(frac(k))
            if (frac(k) >= 1.0_dblprec-10.0_dblprec*epsilon(1.0_dblprec)) then
               frac(k)=0.0_dblprec
            end if
         end if
      end do
   end subroutine canonicalize_fractional


   !----------------------------------------------------------------------------
   ! SUBROUTINE: choose_bins
   !> Select approximately O(N) bins while ensuring that each fractional bin
   !> width is at least twice the map-tolerance bound in that component.
   !----------------------------------------------------------------------------
   subroutine choose_bins(Natom,frac_range,recip_norm,map_tol,nbin)
      implicit none
      integer, intent(in) :: Natom
      real(dblprec), dimension(3), intent(in) :: frac_range,recip_norm
      real(dblprec), intent(in) :: map_tol
      integer, dimension(3), intent(out) :: nbin

      real(dblprec) :: geom, base
      real(dblprec), dimension(3) :: scaled
      integer, dimension(3) :: max_by_tol
      integer :: k
      integer(kind=8) :: prod_bins, max_total

      geom=(max(frac_range(1),tiny(1.0_dblprec))* &
            max(frac_range(2),tiny(1.0_dblprec))* &
            max(frac_range(3),tiny(1.0_dblprec)))**(1.0_dblprec/3.0_dblprec)

      base=max(1.0_dblprec,real(Natom,dblprec)**(1.0_dblprec/3.0_dblprec))

      do k=1,3
         scaled(k)=base*frac_range(k)/geom
         nbin(k)=max(1,nint(scaled(k)))

         if (recip_norm(k)*map_tol > tiny(1.0_dblprec)) then
            max_by_tol(k)=max(1,int(frac_range(k)/(2.0_dblprec*recip_norm(k)*map_tol)))
            nbin(k)=min(nbin(k),max_by_tol(k))
         end if
      end do

      ! Bound spatial-index memory. Reducing bin counts only increases bin width
      ! and therefore preserves the immediate-neighbour-bin completeness proof.
      max_total=max(64_8,8_8*int(Natom,kind=8))
      prod_bins=int(nbin(1),8)*int(nbin(2),8)*int(nbin(3),8)

      do while (prod_bins > max_total)
         k=maxloc(nbin,dim=1)
         nbin(k)=max(1,(nbin(k)+1)/2)
         prod_bins=int(nbin(1),8)*int(nbin(2),8)*int(nbin(3),8)
      end do
   end subroutine choose_bins


   !----------------------------------------------------------------------------
   ! SUBROUTINE: fractional_bin
   !----------------------------------------------------------------------------
   pure subroutine fractional_bin(frac,frac_min,frac_range,nbin,periodic,ib)
      implicit none
      real(dblprec), dimension(3), intent(in) :: frac,frac_min,frac_range
      integer, dimension(3), intent(in) :: nbin
      logical, dimension(3), intent(in) :: periodic
      integer, dimension(3), intent(out) :: ib

      real(dblprec) :: x
      integer :: k

      do k=1,3
         x=(frac(k)-frac_min(k))/frac_range(k)
         if (periodic(k)) x=x-floor(x)
         ib(k)=1+int(x*real(nbin(k),dblprec))
         ib(k)=min(nbin(k),max(1,ib(k)))
      end do
   end subroutine fractional_bin


   !----------------------------------------------------------------------------
   ! FUNCTION: valid_or_wrap_bin
   !> Wrap periodic bin indices in-place; reject out-of-range open indices.
   !----------------------------------------------------------------------------
   logical function valid_or_wrap_bin(ib,nbin,periodic)
      implicit none
      integer, dimension(3), intent(inout) :: ib
      integer, dimension(3), intent(in) :: nbin
      logical, dimension(3), intent(in) :: periodic
      integer :: k

      valid_or_wrap_bin=.true.
      do k=1,3
         if (periodic(k)) then
            ib(k)=modulo(ib(k)-1,nbin(k))+1
         else if (ib(k)<1 .or. ib(k)>nbin(k)) then
            valid_or_wrap_bin=.false.
            return
         end if
      end do
   end function valid_or_wrap_bin


   !----------------------------------------------------------------------------
   ! FUNCTION: linear_bin
   !----------------------------------------------------------------------------
   pure integer function linear_bin(ib,nbin)
      implicit none
      integer, dimension(3), intent(in) :: ib,nbin

      linear_bin=ib(1)+(ib(2)-1)*nbin(1)+ &
                 (ib(3)-1)*nbin(1)*nbin(2)
   end function linear_bin


   !----------------------------------------------------------------------------
   ! SUBROUTINE: append_unique_bin
   !> Avoid repeated bins when a periodic direction has only one or two bins.
   !----------------------------------------------------------------------------
   pure subroutine append_unique_bin(candidate,bin_ids,n_bin_ids)
      implicit none
      integer, intent(in) :: candidate
      integer, dimension(27), intent(inout) :: bin_ids
      integer, intent(inout) :: n_bin_ids
      integer :: i

      do i=1,n_bin_ids
         if (bin_ids(i)==candidate) return
      end do

      n_bin_ids=n_bin_ids+1
      bin_ids(n_bin_ids)=candidate
   end subroutine append_unique_bin


   !----------------------------------------------------------------------------
   ! SUBROUTINE: get_symops
   !> Find possible symmetry operations depending on assumed symmetry.
   !----------------------------------------------------------------------------
   subroutine get_symops(isym,NT)
      implicit none

      integer, intent(in) :: isym
      integer, intent(in) :: NT

      integer :: i,j,x,y,z,j_s,x1,x2,y1,y2,k
      integer :: i_stat
      integer, dimension(NT) :: sym_count
      real(dblprec) :: half,roothalf

      sym_count=0

      if (isym==0) then
         sym_count=1
         allocate(sym_mats(3,3,1,NT),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats), &
                     'sym_mats','get_symops')
         if (i_stat/=0) error stop 'get_symops: allocation failure'
         sym_mats=0.0_dblprec
         do k=1,NT
            do i=1,3
               sym_mats(i,i,1,k)=1.0_dblprec
            end do
         end do

      else if (isym==1) then
         allocate(sym_mats(3,3,48,NT),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats), &
                     'sym_mats','get_symops')
         if (i_stat/=0) error stop 'get_symops: allocation failure'
         sym_mats=0.0_dblprec
         sym_count=0
         do k=1,NT
            do i=1,3
               do j=0,1
                  j_s=(-1)**j
                  do x=0,1
                     do y=0,1
                        do z=0,1
                           sym_count(k)=sym_count(k)+1
                           sym_mats(1,mod(i-j_s,3)+1,sym_count(k),k)= &
                              (-1.0_dblprec)**x
                           sym_mats(2,mod(i,3)+1,sym_count(k),k)= &
                              (-1.0_dblprec)**y
                           sym_mats(3,mod(i+j_s,3)+1,sym_count(k),k)= &
                              (-1.0_dblprec)**z
                        end do
                     end do
                  end do
               end do
            end do
         end do

      else if (isym==2) then
         allocate(sym_mats(3,3,12,NT),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats), &
                     'sym_mats','get_symops')
         if (i_stat/=0) error stop 'get_symops: allocation failure'
         sym_mats=0.0_dblprec
         sym_count=0
         do k=1,NT
            do j=0,1
               do x=0,1
                  do y=0,1
                     sym_count(k)=sym_count(k)+1
                     sym_mats(1,mod(j,2)+1,sym_count(k),k)= &
                        (-1.0_dblprec)**x
                     sym_mats(2,mod(j+1,2)+1,sym_count(k),k)= &
                        (-1.0_dblprec)**y
                     sym_mats(3,3,sym_count(k),k)=1.0_dblprec
                  end do
               end do
            end do
         end do

      else if (isym==3) then
         allocate(sym_mats(3,3,24,NT),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats), &
                     'sym_mats','get_symops')
         if (i_stat/=0) error stop 'get_symops: allocation failure'
         sym_mats=0.0_dblprec
         sym_count=0
         half=0.50_dblprec
         roothalf=sqrt(3.0_dblprec)*0.50_dblprec

         do k=1,NT
            do x=0,1
               do y=0,1
                  do z=0,1
                     sym_count(k)=sym_count(k)+1
                     sym_mats(1,1,sym_count(k),k)=(-1.0_dblprec)**x
                     sym_mats(2,2,sym_count(k),k)=(-1.0_dblprec)**y
                     sym_mats(3,3,sym_count(k),k)=(-1.0_dblprec)**z
                  end do
               end do
            end do
         end do

         do k=1,NT
            do x1=0,1
               do x2=0,1
                  do y1=0,1
                     do y2=0,1
                        if ((-1.0_dblprec)**x1*(-1.0_dblprec)**x2* &
                            (-1.0_dblprec)**y1*(-1.0_dblprec)**y2 < &
                            0.0_dblprec) then
                           do z=0,1
                              sym_count(k)=sym_count(k)+1
                              sym_mats(1,1,sym_count(k),k)= &
                                 (-1.0_dblprec)**x1*half
                              sym_mats(2,1,sym_count(k),k)= &
                                 (-1.0_dblprec)**x2*roothalf
                              sym_mats(1,2,sym_count(k),k)= &
                                 (-1.0_dblprec)**y1*roothalf
                              sym_mats(2,2,sym_count(k),k)= &
                                 (-1.0_dblprec)**y2*half
                              sym_mats(3,3,sym_count(k),k)= &
                                 (-1.0_dblprec)**z
                           end do
                        end if
                     end do
                  end do
               end do
            end do
         end do

      else if (isym==4) then
         open(ifileno,file='sym.mat')
         read(ifileno,*) sym_count(1)
         allocate(sym_mats(3,3,sym_count(1),NT),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats), &
                     'sym_mats','get_symops')
         if (i_stat/=0) error stop 'get_symops: allocation failure'
         do j=1,sym_count(1)
            do x=1,3
               read(ifileno,*) (sym_mats(x,y,j,1),y=1,3)
            end do
            do k=2,NT
               sym_mats(1:3,1:3,j,k)=sym_mats(1:3,1:3,j,1)
            end do
         end do
         do k=2,NT
            sym_count(k)=sym_count(1)
         end do
         close(ifileno)

      else if (isym==5) then
         allocate(sym_mats(3,3,48,NT),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats), &
                     'sym_mats','get_symops')
         if (i_stat/=0) error stop 'get_symops: allocation failure'
         open(ifileno,file='sym.mat')
         sym_mats=0.0_dblprec
         sym_count=0
         do k=1,NT
            read(ifileno,*) sym_count(k)
            if (sym_count(k)>48) then
               error stop 'get_symops: more than 48 type-specific operations'
            end if
            do j=1,sym_count(k)
               do x=1,3
                  read(ifileno,*) (sym_mats(y,x,j,k),y=1,3)
               end do
            end do
         end do
         close(ifileno)

      else
         error stop 'get_symops: unsupported symmetry selector'
      end if

      nsym(1:NT)=sym_count(1:NT)
   end subroutine get_symops


   !----------------------------------------------------------------------------
   ! SUBROUTINE: get_fullnnlist
   !> Create the symmetry-expanded neighbour vectors and remove duplicates.
   !----------------------------------------------------------------------------
   subroutine get_fullnnlist(NT,NN,Nelem,max_no_shells,redcoord, &
                             max_no_equiv,nncoord)
      implicit none

      integer, intent(in) :: NT
      integer, intent(in) :: Nelem
      integer, intent(in) :: max_no_equiv
      integer, intent(in) :: max_no_shells
      integer, dimension(NT), intent(in) :: NN
      real(dblprec), dimension(NT,max_no_shells,3), intent(in) :: redcoord
      real(dblprec), dimension(3,max_no_equiv,max_no_shells,NT), &
         intent(out) :: nncoord

      real(dblprec), parameter :: sym_tol2 = sym_tol_default*sym_tol_default
      real(dblprec), dimension(3) :: tvect
      integer :: counter
      logical :: unique
      integer :: i,j,k,itype,ishell,isym

      ! Nelem is retained for interface compatibility with the ordinary path.
      if (Nelem /= 1) then
         error stop 'get_fullnnlist: only pair interactions are supported'
      end if

      nncoord=0.0_dblprec

      do itype=1,NT
         do ishell=1,NN(itype)
            if (nsym(itype)==1) then
               nncoord(:,1,ishell,itype)=redcoord(itype,ishell,:)
               nmdimt(ishell,itype)=1
            else
               counter=0
               do isym=1,nsym(itype)
                  tvect=0.0_dblprec
                  do i=1,3
                     do j=1,3
                        tvect(i)=tvect(i)+ &
                           redcoord(itype,ishell,j)*sym_mats(i,j,isym,itype)
                     end do
                  end do

                  unique=.true.
                  do k=1,counter
                     if (dot_product(tvect-nncoord(:,k,ishell,itype), &
                                     tvect-nncoord(:,k,ishell,itype)) < &
                         sym_tol2) then
                        unique=.false.
                        exit
                     end if
                  end do

                  if (unique) then
                     counter=counter+1
                     if (counter>max_no_equiv) then
                        error stop 'get_fullnnlist: max_no_equiv too small'
                     end if
                     nncoord(:,counter,ishell,itype)=tvect
                  end if
               end do
               nmdimt(ishell,itype)=counter
            end if
         end do
      end do
   end subroutine get_fullnnlist

end module NeighbourMapSFC
