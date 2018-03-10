!> Common Monte Carlo routines used for both non-LSF and LSF calculations
!> @brief
!! Collection of common Monte Carlo routines used in several different modes
!> @author
!! Lars Bergqvist
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!--------------------------------------------------------------------------
module montecarlo_common

   use Parameters
   implicit none

   private

   public :: choose_random_flip, flip_a, flip_g, flip_h, calculate_energy, Ising_random_flip

contains

   !> Calculates a random spin direction
   subroutine choose_random_flip(emom,newmom,Natom,Mensemble,iflip,k,temperature,rn)

      !  Hinzke-Nowak update, random between uniform rotation, gaussian rotation and spin-flips
      !  LB 150306

      use RandomNumbers, only : rng_uniform,rng_gaussian
      use Constants,only : mub,k_bolt

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in the sample
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3), intent(out) :: newmom !< New trial moment
      integer, intent(in) :: iflip     !< Atom to flip
      integer, intent(in) :: k         !< current ensemble
      real(dblprec), intent(in):: temperature !< Temperature
      real(dblprec), intent(in):: rn !< Random number

      real(dblprec), dimension(3) :: gasmom !< Gaussian vector
      real(dblprec) :: newmoml,delta
      integer :: ftype
      ftype=floor(3*rn)

      if (ftype==0) then             ! Uniform rotation
         call rng_uniform(gasmom,3)
         newmom(:)=2.0d0*gasmom(:)-1.0d0
         do while (newmom(1)**2+newmom(2)**2+newmom(3)**2>1.00d0 .or. newmom(1)**2+newmom(2)**2+newmom(3)**2<dbl_tolerance)
            call rng_uniform(gasmom,3)
            newmom(:)=2.0d0*gasmom(:)-1.0d0
         end do
         newmoml=sqrt(newmom(1)**2+newmom(2)**2+newmom(3)**2)
         newmom(:) = newmom(:)/newmoml
      elseif (ftype==1) then                            ! Gaussian rotation
         delta=(2.0/25.0)*(k_bolt*temperature/mub)**(0.2d0)
         call rng_gaussian(gasmom,3,delta)
         newmoml=sqrt((emom(1,iflip,k)+gasmom(1))**2+(emom(2,iflip,k)+gasmom(2))**2+(emom(3,iflip,k)+gasmom(3))**2)
         newmom(:)=(emom(:,iflip,k)+gasmom(:))/newmoml
      else                                             !spin-flip
         newmom(:)=-emom(:,iflip,k)
      endif

   end subroutine choose_random_flip

   !---------------------------------------------------------------------------
   !> @brief
   !> Flip selected spin if energetically favourable using Metropolis
   !>
   !> @author
   !> Lars Bergqvist
   !> Jonathan Chico ----> Implementation of induced moments fluctuations
   !---------------------------------------------------------------------------
   subroutine flip_a(Natom, Mensemble, emom, emomM, mmom, iflip,newmom,newmmom,de,&
         temperature,temprescale, do_lsf,k,flipprob,lsf_metric,ind_nlistsize,&
         ind_nlist,ind_mom_flag,max_no_neigh,sus_ind)

      use Constants

      !.. Implicit declarations
      implicit none

      !.. Input variables
      ! System/geometry variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in):: temperature !< Temperature
      real(dblprec), intent(in):: temprescale !< Temperature rescaling from QHB
      ! Moment variables
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), intent(in) :: newmmom !< New trial magnitude of moment
      real(dblprec), dimension(3), intent(in) :: newmom !< New trial moment
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! Probabilities variables
      real(dblprec), intent(in):: de  !< Energy difference
      real(dblprec), intent(in) :: flipprob !< Probability of flipping spin
      ! LSF variables
      integer, intent(in) :: lsf_metric !< LSF metric or phase space measure
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      ! Induced moment variables
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      character(len=1), intent(in) :: ind_mom_flag

      integer :: k !< Current ensemble
      integer :: ind_neigh,curr_ind,fix_neigh
      real(dblprec) :: beta,des,lmetric
      real(dblprec), dimension(3) :: ave_mom

      beta=1d0/k_bolt/(temprescale*Temperature+1.0d-15)
      if (do_lsf=='N'.and.ind_mom_flag=='N') then
         if(de<=0.0d0 .or. flipprob<exp(-beta*de)) then
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
         endif
      else if (do_lsf=='N'.and.ind_mom_flag=='Y') then
         if(de<=0.0d0 .or. flipprob<exp(-beta*de)) then
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
            ! If the flip is accepted then the magnitude of the induced moments must change
            do ind_neigh=1,ind_nlistsize(iflip)
               ! Select the current induced moment
               curr_ind=ind_nlist(ind_neigh,iflip)
               ! Sum over the nearest neighbours that are fixed
               ave_mom=0.0d0
               do fix_neigh=1, ind_nlistsize(curr_ind)
                  ave_mom(1)=ave_mom(1)+emomM(1,ind_nlist(fix_neigh,curr_ind),k)
                  ave_mom(2)=ave_mom(2)+emomM(2,ind_nlist(fix_neigh,curr_ind),k)
                  ave_mom(3)=ave_mom(3)+emomM(3,ind_nlist(fix_neigh,curr_ind),k)
               enddo
               ! Vary the magnitude of the induced moments
               emomM(:,curr_ind,k)=ave_mom(:)*sus_ind(curr_ind)
               mmom(curr_ind,k)=sqrt(emomM(1,curr_ind,k)**2+emomM(2,curr_ind,k)**2+emomM(3,curr_ind,k)**2)
               emom(:,curr_ind,k)=emomM(:,curr_ind,k)/mmom(curr_ind,k)
            enddo
         endif
      else
         if(lsf_metric==1) then
            lmetric=1.d0
         else
            lmetric=(mmom(iflip,k)/newmmom)**2
         endif
         des=-log(lmetric)/beta
         if(de<=des .or. flipprob<exp(-beta*de)/lmetric) then
            mmom(iflip,k) = newmmom
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
         endif
      endif
   end subroutine flip_a


   !> Flip selected spin if energetically favourable using Glauber dynamics
   subroutine flip_g(Natom, Mensemble, emom, emomM, mmom, iflip,newmom,newmmom,de,&
         temperature,temprescale, do_lsf,k,flipprob,lsf_metric,ind_nlistsize,&
         ind_nlist,ind_mom_flag,max_no_neigh,sus_ind)

      use Constants

      !.. Implicit declarations
      implicit none

      !.. Input variables
      ! System/geometry variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in):: temperature !< Temperature
      real(dblprec), intent(in):: temprescale !< Temperature rescaling from QHB
      ! Moment variables
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), intent(in) :: newmmom !< New trial magnitude of moment
      real(dblprec), dimension(3), intent(in) :: newmom !< New trial moment
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! Probabilities variables
      real(dblprec), intent(in):: de  !< Energy difference
      real(dblprec), intent(in) :: flipprob !< Probability of flipping spin
      ! LSF variables
      integer, intent(in) :: lsf_metric !< LSF metric or phase space measure
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      ! Induced moment variables
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      character(len=1), intent(in) :: ind_mom_flag

      integer :: k !< Current ensemble
      integer :: ind_neigh,curr_ind,fix_neigh
      real(dblprec) :: beta,lmetric
      real(dblprec), dimension(3) :: ave_mom

      beta=1d0/k_bolt/(temprescale*Temperature+1.0d-15)
      if (do_lsf=='N'.and.ind_mom_flag=='N') then
         if(flipprob<1.d0/(1.d0+exp(beta*de))) then
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
            ! If induced moments are considered
         endif
      elseif (do_lsf=='N'.and.ind_mom_flag=='Y') then
         if(de<=0.0d0 .or. flipprob<exp(-beta*de)) then
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
            ! If the flip is accepted then the magnitude of the induced moments must change
            do ind_neigh=1,ind_nlistsize(iflip)
               ! Select the current induced moment
               curr_ind=ind_nlist(ind_neigh,iflip)
               ! Sum over the nearest neighbours that are fixed
               ave_mom=0.0d0
               do fix_neigh=1, ind_nlistsize(curr_ind)
                  ave_mom(1)=ave_mom(1)+emomM(1,ind_nlist(fix_neigh,curr_ind),k)
                  ave_mom(2)=ave_mom(2)+emomM(2,ind_nlist(fix_neigh,curr_ind),k)
                  ave_mom(3)=ave_mom(3)+emomM(3,ind_nlist(fix_neigh,curr_ind),k)
               enddo
               ! Vary the magnitude of the induced moments
               emomM(:,curr_ind,k)=ave_mom(:)*sus_ind(curr_ind)
               mmom(curr_ind,k)=sqrt(emomM(1,curr_ind,k)**2+emomM(2,curr_ind,k)**2+emomM(3,curr_ind,k)**2)
               emom(:,curr_ind,k)=emomM(:,curr_ind,k)/mmom(curr_ind,k)
            enddo
         endif
      else
         if(lsf_metric==1) then
            lmetric=1.d0
         else
            lmetric=(mmom(iflip,k)/newmmom)**2
         endif
         if(flipprob<1.d0/(1.d0+lmetric*exp(beta*de))) then
            mmom(iflip,k) = newmmom
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
         endif
      endif
   end subroutine flip_g

   !> Flip selected spin to the direction of heat bath
   subroutine flip_h(Natom, Mensemble, emom, emomM, newmmom, mmom, iflip,temperature,temprescale,k,flipprob,totfield)

      use Constants
      use RandomNumbers, only : rng_uniform

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), intent(in) :: newmmom !< New trial magnitude of moment
      real(dblprec), intent(out) :: mmom !< Magnitude of magnetic moments
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec),intent(in) :: flipprob !< Probability of flipping spin
      real(dblprec), intent(in):: temperature !< Temperature
      real(dblprec), intent(in):: temprescale !< Temperature rescaling from QHB
      real(dblprec), dimension(3),intent(in) :: totfield  !<Total effective field

      integer :: k !< Current ensemble
      real(dblprec) :: beta,zarg,zctheta,zstheta,zcphi,zsphi,ctheta,stheta,phi,q(1)
      real(dblprec),dimension(3) :: stemp,zfc

      mmom=newmmom
      beta=1.d0/k_bolt/(temprescale*Temperature)
      zfc=beta*totfield*mub*mmom
      zarg=sqrt(sum(zfc(:)*zfc(:)))
      zctheta=zfc(3)/zarg
      zstheta=sqrt(1.d0-zctheta*zctheta)
      if (zstheta < 1.d-2) then
         zctheta=0.99995d0
         zstheta=1.d-2
         zcphi=1.d0
         zsphi=0.d0
      else
         zcphi=zfc(1)/(zarg*zstheta)
         zsphi=zfc(2)/(zarg*zstheta)
      endif
      ctheta=1.5d0
      do while (abs(ctheta)>1.d0)
         call rng_uniform(q,1)
         ctheta=1.d0+1.d0/zarg*log(q(1))   ! Modified Direct Heat Bath
      enddo
      stheta=sqrt(1.d0-ctheta*ctheta)
      phi=pi*(2.d0*flipprob-1.d0)
      stemp(1)=stheta*cos(phi)              ! New spin in direction of the effective field
      stemp(2)=stheta*sin(phi)
      stemp(3)=ctheta
      !--Euler rotation of spin from the Heff frame of ref. to cryst. axis-----
      emom(1,iflip,k)=zcphi*zctheta*stemp(1)-zsphi*stemp(2)+zcphi*zstheta*stemp(3)
      emom(2,iflip,k)=zsphi*zctheta*stemp(1)+zcphi*stemp(2)+zsphi*zstheta*stemp(3)
      emom(3,iflip,k)=     -zstheta*stemp(1)                     +zctheta*stemp(3)
      emomM(:,iflip,k)=mmom*emom(:,iflip,k)
   end subroutine flip_h

   !> Calculate total energy of single spin
   !> @todo Remove dipolar or activate it, right now dead piece of code
   subroutine calculate_energy(Natom, Mensemble, max_no_neigh, conf_num, ncoup, ncoupD, nlist, nlistsize, &
         do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
         do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
         taniso, eaniso, kaniso,sb,emomM, emom, mmom, iflip, newmom, extfield, de, k, &
         mult_axis,taniso_diff, eaniso_diff, kaniso_diff,sb_diff,&
         do_dip, Qdip,exc_inter)

      use Constants, only : mub

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
      integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
      integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
      integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
      integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange coupling
      integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
      integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
      real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
      integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
      integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      integer, dimension(Natom),intent(in) :: taniso_diff !< Type of anisotropy (0-2)
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
      real(dblprec), dimension(Natom), intent(in) :: sb!< Ratio between the Anisotropy constants
      real(dblprec), dimension(Natom), intent(in) :: sb_diff!< Ratio between the Anisotropy constants
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), dimension(3), intent(in) :: newmom !< New trial moment
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      real(dblprec), intent(out):: de  !< Energy difference
      integer, intent(in) :: k !< Current ensemble
      real(dblprec) :: aw1,aw2
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction
      character(len=1), intent(in) :: mult_axis
      character(len=1), intent(in) :: exc_inter !> Interpolation of Jij between FM/DLM (Y/N)

      !.. Local scalars
      integer :: j
      real(dblprec) :: tt, tta, ttb, e_c, e_t
      real(dblprec) :: bqmdot, excscale

      !.. Local arrays
      real(dblprec), dimension(3) :: beff_t, trialmom

      !.. Executable statements

      ! First calculate effective field
      beff_t(1) = 0d0
      beff_t(2) = 0d0
      beff_t(3) = 0d0
      tt=0.d0

      e_c=0.d0
      e_t=0.d0
      trialmom(:)=newmom(:)*mmom(iflip,k)

      if (exc_inter=='N') then
         ! Exchange
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
         !$omp simd reduction(+:e_c,e_t)
#endif
         do j=1,nlistsize(iflip)
            e_c=e_c-ncoup(j,iflip,1)*sum(emomM(:,iflip,k)*emomM(:,nlist(j,iflip),k))
            e_t=e_t-ncoup(j,iflip,1)*sum(trialmom(:)*emomM(:,nlist(j,iflip),k))
         end do
      else
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
         !       !$omp simd reduction(+:beff_t,tt)
#endif
         do j=1,nlistsize(iflip)
            beff_t=beff_t+emomM(:,nlist(j,iflip),k)
            tt=tt+mmom(nlist(j,iflip),k)
         enddo
         excscale=sqrt(beff_t(1)**2+beff_t(2)**2+beff_t(3)**2)/tt
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
         !$omp simd reduction(+:e_c,e_t)
#endif
         do j=1,nlistsize(iflip)
            e_c=e_c-(excscale*ncoup(j,iflip,1)+(1.d0-excscale)*ncoupD(j,iflip,1))*sum(emomM(:,iflip,k)*emomM(:,nlist(j,iflip),k))
            e_t=e_t-(excscale*ncoup(j,iflip,1)+(1.d0-excscale)*ncoupD(j,iflip,1))*sum(trialmom(:)*emomM(:,nlist(j,iflip),k))
         end do
      endif

      ! Anisotropy

      ! Uniaxial anisotropy
      if (taniso(iflip)==1) then
         tta=sum(emomM(:,iflip,k)*eaniso(:,iflip))
         ttb=sum(trialmom(:)*eaniso(:,iflip))
         e_c=e_c+kaniso(1,iflip)*(tta**2)+kaniso(2,iflip)*(tta**4)
         e_t=e_t+kaniso(1,iflip)*(ttb**2)+kaniso(2,iflip)*(ttb**4)
         ! Cubic anisotropy
      elseif (taniso(iflip)==2) then
         e_c=e_c-kaniso(1,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
            emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)- &
            kaniso(2,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
         e_t=e_t-kaniso(1,iflip)*(trialmom(1)**2*trialmom(2)**2+ &
            trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)- &
            kaniso(2,iflip)*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
      endif
      ! When both Cubic and Uniaxial are switched on
      if (taniso(iflip)==7) then

         ! Uniaxial anisotropy
         tta=(emomM(1,iflip,k)*eaniso(1,iflip)+emomM(2,iflip,k)*eaniso(2,iflip)+emomM(3,iflip,k)*eaniso(3,iflip))
         ttb=(trialmom(1)*eaniso(1,iflip)+trialmom(2)*eaniso(2,iflip)+trialmom(3)*eaniso(3,iflip))
         e_c=e_c+kaniso(1,iflip)*(tta**2)+kaniso(2,iflip)*(tta**4)
         e_t=e_t+kaniso(1,iflip)*(ttb**2)+kaniso(2,iflip)*(ttb**4)

         ! Cubic anisotropy
         aw1=kaniso(1,iflip)*sb(iflip)
         aw2=kaniso(2,iflip)*sb(iflip)

         e_c=e_c+aw1*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
            emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+ &
            aw2*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
         e_t=e_t+aw1*(trialmom(1)**2*trialmom(2)**2+ &
            trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
            aw2*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
      endif

      ! This includes a secondary anisotropy energy
      if (mult_axis=='Y') then
         ! Uniaxial anisotropy
         if (taniso_diff(iflip)==1) then
            tta=(emomM(1,iflip,k)*eaniso_diff(1,iflip)+emomM(2,iflip,k)*eaniso_diff(2,iflip)+emomM(3,iflip,k)*eaniso_diff(3,iflip))
            ttb=(trialmom(1)*eaniso_diff(1,iflip)+trialmom(2)*eaniso_diff(2,iflip)+trialmom(3)*eaniso_diff(3,iflip))
            e_c=e_c+kaniso_diff(1,iflip)*(tta**2)+kaniso_diff(2,iflip)*(tta**4)
            e_t=e_t+kaniso_diff(1,iflip)*(ttb**2)+kaniso_diff(2,iflip)*(ttb**4)

            ! Cubic anisotropy
         elseif (taniso_diff(iflip)==2) then
            e_c=e_c+kaniso_diff(1,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
               emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+ &
               kaniso_diff(2,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
            e_t=e_t+kaniso_diff(1,iflip)*(trialmom(1)**2*trialmom(2)**2+ &
               trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
               kaniso_diff(2,iflip)*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)

         endif
         ! When both Cubic and Uniaxial are switched on
         if (taniso_diff(iflip)==7) then

            ! Uniaxial anisotropy
            tta=(emomM(1,iflip,k)*eaniso_diff(1,iflip)+emomM(2,iflip,k)*eaniso_diff(2,iflip)+emomM(3,iflip,k)*eaniso_diff(3,iflip))
            ttb=(trialmom(1)*eaniso_diff(1,iflip)+trialmom(2)*eaniso_diff(2,iflip)+trialmom(3)*eaniso_diff(3,iflip))
            e_c=e_c+kaniso_diff(1,iflip)*(tta**2)+kaniso_diff(2,iflip)*(tta**4)
            e_t=e_t+kaniso_diff(1,iflip)*(ttb**2)+kaniso_diff(2,iflip)*(ttb**4)

            ! Cubic anisotropy
            aw1=kaniso_diff(1,iflip)*sb_diff(iflip)
            aw2=kaniso_diff(2,iflip)*sb_diff(iflip)

            e_c=e_c+aw1*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
               emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+ &
               aw2*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
            e_t=e_t+aw1*(trialmom(1)**2*trialmom(2)**2+ &
               trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
               aw2*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
         endif
      endif


      ! DM interaction
      if (do_dm==1) then
         do j=1,dmlistsize(iflip)
            e_c=e_c+dm_vect(1,j,iflip)*(emomM(2,iflip,k)*emomM(3,dmlist(j,iflip),k)- &
               emom(3,iflip,k)*emomM(2,dmlist(j,iflip),k))+ &
               dm_vect(2,j,iflip)*(emomM(3,iflip,k)*emomM(1,dmlist(j,iflip),k)- &
               emomM(1,iflip,k)*emomM(3,dmlist(j,iflip),k))+ &
               dm_vect(3,j,iflip)*(emom(1,iflip,k)*emomM(2,dmlist(j,iflip),k)- &
               emomM(2,iflip,k)*emomM(1,dmlist(j,iflip),k))
            e_t=e_t+dm_vect(1,j,iflip)*(trialmom(2)*emomM(3,dmlist(j,iflip),k)- &
               trialmom(3)*emomM(2,dmlist(j,iflip),k))+ &
               dm_vect(2,j,iflip)*(trialmom(3)*emomM(1,dmlist(j,iflip),k)- &
               trialmom(1)*emomM(3,dmlist(j,iflip),k))+ &
               dm_vect(3,j,iflip)*(trialmom(1)*emomM(2,dmlist(j,iflip),k)- &
               trialmom(2)*emomM(1,dmlist(j,iflip),k))

         end do
      end if

      ! PD interaction
      if(do_pd==1) then
         do j=1,pdlistsize(iflip)
            e_c=e_c-pd_vect(1,j,iflip)*emomM(1,iflip,k)*emomM(1,pdlist(j,iflip),k)- &
               pd_vect(4,j,iflip)*emomM(1,iflip,k)*emomM(2,pdlist(j,iflip),k)- &
               pd_vect(5,j,iflip)*emomM(1,iflip,k)*emomM(3,pdlist(j,iflip),k)- &
               pd_vect(4,j,iflip)*emomM(2,iflip,k)*emomM(1,pdlist(j,iflip),k)- &
               pd_vect(2,j,iflip)*emomM(2,iflip,k)*emomM(2,pdlist(j,iflip),k)- &
               pd_vect(6,j,iflip)*emomM(2,iflip,k)*emomM(3,pdlist(j,iflip),k)- &
               pd_vect(5,j,iflip)*emomM(3,iflip,k)*emomM(1,pdlist(j,iflip),k)- &
               pd_vect(6,j,iflip)*emomM(3,iflip,k)*emomM(2,pdlist(j,iflip),k)- &
               pd_vect(3,j,iflip)*emomM(3,iflip,k)*emomM(3,pdlist(j,iflip),k)

            e_t=e_t-pd_vect(1,j,iflip)*trialmom(1)*emomM(1,pdlist(j,iflip),k)- &
               pd_vect(4,j,iflip)*trialmom(1)*emomM(2,pdlist(j,iflip),k)- &
               pd_vect(5,j,iflip)*trialmom(1)*emomM(3,pdlist(j,iflip),k)- &
               pd_vect(4,j,iflip)*trialmom(2)*emomM(1,pdlist(j,iflip),k)- &
               pd_vect(2,j,iflip)*trialmom(2)*emomM(2,pdlist(j,iflip),k)- &
               pd_vect(6,j,iflip)*trialmom(2)*emomM(3,pdlist(j,iflip),k)- &
               pd_vect(5,j,iflip)*trialmom(3)*emomM(1,pdlist(j,iflip),k)- &
               pd_vect(6,j,iflip)*trialmom(3)*emomM(2,pdlist(j,iflip),k)- &
               pd_vect(3,j,iflip)*trialmom(3)*emomM(3,pdlist(j,iflip),k)
         end do
      end if

      ! BIQDM interaction
      if(do_biqdm==1) then
         do j=1,biqdmlistsize(iflip)
            e_c=e_c-biqdm_vect(1,j,iflip)*(emomM(2,iflip,k)*emomM(3,biqdmlist(j,iflip),k)- &
               emomM(3,iflip,k)*emomM(2,biqdmlist(j,iflip),k))**2- &
               biqdm_vect(1,j,iflip)*(emom(3,iflip,k)*emomM(1,biqdmlist(j,iflip),k)- &
               emomM(1,iflip,k)*emomM(3,biqdmlist(j,iflip),k))**2- &
               biqdm_vect(1,j,iflip)*(emom(1,iflip,k)*emomM(2,biqdmlist(j,iflip),k)- &
               emomM(2,iflip,k)*emomM(1,biqdmlist(j,iflip),k))**2

            e_t=e_t-biqdm_vect(1,j,iflip)*(trialmom(2)*emomM(3,biqdmlist(j,iflip),k)- &
               trialmom(3)*emomM(2,biqdmlist(j,iflip),k))**2- &
               biqdm_vect(1,j,iflip)*(trialmom(3)*emomM(1,biqdmlist(j,iflip),k)- &
               trialmom(1)*emomM(3,biqdmlist(j,iflip),k))**2- &
               biqdm_vect(1,j,iflip)*(trialmom(1)*emomM(2,biqdmlist(j,iflip),k)- &
               trialmom(2)*emomM(1,biqdmlist(j,iflip),k))**2
         end do
      end if

      ! BQ interaction
      if(do_bq==1) then

         do j=1,bqlistsize(iflip)
            ! current spin
            bqmdot=emomM(1,bqlist(j,iflip),k)*emomM(1,iflip,k)+&
               emomM(2,bqlist(j,iflip),k)*emomM(2,iflip,k)+&
               emomM(3,bqlist(j,iflip),k)*emomM(3,iflip,k)
            e_c=e_c-j_bq(j,iflip)*bqmdot**2

            !trial spin
            bqmdot=emomM(1,bqlist(j,iflip),k)*trialmom(1) + &
               emomM(2,bqlist(j,iflip),k)*trialmom(2) + &
               emomM(3,bqlist(j,iflip),k)*trialmom(3)
            e_t=e_t-j_bq(j,iflip)*bqmdot**2

         end do
      end if

      ! Add Zeeman term
      e_c=e_c-extfield(1)*emomM(1,iflip,k)-extfield(2)*emomM(2,iflip,k)-extfield(3)*emomM(3,iflip,k)
      e_t=e_t-extfield(1)*trialmom(1)-extfield(2)*trialmom(2)-extfield(3)*trialmom(3)

      !Energy difference
      tt=e_t-e_c
      de=mub*tt
      return
   end subroutine calculate_energy


   !> Flip Ising spin
   subroutine Ising_random_flip(emom,newmom,i,en,Natom,Mensemble)

      implicit none

      integer, intent(in) :: i  ! Atomic site which is going to be flipped
      integer, intent(in) :: en ! Enseble flag
      integer, intent(in) :: Natom !< Number of atoms in the sample
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3), intent(out) :: newmom !< New moment from the flipping

      newmom(1) = -1.0d0*emom(1,i,en)
      newmom(2) = -1.0d0*emom(2,i,en)
      newmom(3) = -1.0d0*emom(3,i,en)

   end subroutine Ising_random_flip

end module montecarlo_common
