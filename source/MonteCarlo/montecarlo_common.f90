!--------------------------------------------------------------------------
! MODULE: montecarlo_common
! DESCRIPTION:
!> Common Monte Carlo routines used for both non-LSF and LSF calculations
!> @brief
!> Collection of common Monte Carlo routines used in several different modes
!> @author
!> Lars Bergqvist
!> @copyright
!> GNU Public License.
!--------------------------------------------------------------------------
module montecarlo_common

   use Parameters
   implicit none

   private

   public :: choose_random_flip, flip_a, flip_g, flip_h, calculate_energy, Ising_random_flip
   public :: randomize_spins

contains

   !> Calculates a random spin direction
   subroutine choose_random_flip(emom,newmom,Natom,Mensemble,iflip,k,delta,rn,gn)

      !  Hinzke-Nowak update, random between uniform rotation, gaussian rotation and spin-flips
      !  LB 150306

      use RandomNumbers, only : rng_uniform,rng_gaussian
!      use Constants,only : mub,k_bolt,pi
      use Constants,only : pi

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in the sample
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3), intent(out) :: newmom !< New trial moment
      integer, intent(in) :: iflip     !< Atom to flip
      integer, intent(in) :: k         !< current ensemble
      real(dblprec), intent(in):: delta !< Step size in gaussian rotation
      real(dblprec), dimension(3),intent(in):: rn !< Random number
      real(dblprec), dimension(3), intent(in):: gn !< gaussian Random number

      real(dblprec), dimension(3) :: gasmom !< Gaussian vector
      real(dblprec) :: newmoml,theta,phi
      integer :: ftype

      ftype=floor(3*rn(1))

      if (ftype==0) then             ! Uniform rotation
!          call rng_uniform(gasmom(1:2),2)
          phi=rn(2)*2*pi
          theta=acos(1-2*rn(3))
          newmom(1)=sin(theta)*cos(phi)
          newmom(2)=sin(theta)*sin(phi)
          newmom(3)=cos(theta)

!         call rng_uniform(gasmom,3)
!         newmom(:)=2.0_dblprec*gasmom(:)-1.0_dblprec
!         do while (newmom(1)**2+newmom(2)**2+newmom(3)**2>1.0_dblprec .or. newmom(1)**2+newmom(2)**2+newmom(3)**2<dbl_tolerance)
!            call rng_uniform(gasmom,3)
!            newmom(:)=2.0_dblprec*gasmom(:)-1.0_dblprec
!         end do
!         newmoml=sqrt(newmom(1)**2+newmom(2)**2+newmom(3)**2)
!         newmom(:) = newmom(:)/newmoml
      elseif (ftype==1) then                            ! Gaussian rotation
!         delta=(2.0/25.0)*(k_bolt*temperature/mub)**(0.20_dblprec)
!         call rng_gaussian(gasmom,3,delta)
         gasmom=gn*delta
         newmoml=sqrt((emom(1,iflip,k)+gasmom(1))**2+(emom(2,iflip,k)+gasmom(2))**2+(emom(3,iflip,k)+gasmom(3))**2)
         newmom(:)=(emom(:,iflip,k)+gasmom(:))/newmoml
      else                                             !spin-flip
         newmom(:)=-emom(:,iflip,k)
      endif

   end subroutine choose_random_flip

   !> Calculates a random spin direction
   subroutine randomize_spins(Natom,Mensemble,emom,emomM,mmom)

      use RandomNumbers, only : rng_uniform
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in the sample
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      !
      real(dblprec), dimension(3) :: newmom !< New trial moment
      real(dblprec), dimension(3) :: urot !< Gaussian random number
      real(dblprec), dimension(3) :: zero !< zero vector
      integer :: i         !< Atom index
      integer :: k         !< current ensemble


      zero=0.0_dblprec
      urot(1)=0.0_dblprec
!      if(use_vsl) then
!#ifdef VSL
!         !$omp parallel do default(shared),private(i,k,newmom),schedule(auto),collapse(2)
!#endif
         do i=1,Natom
            do k=1,mensemble
               call rng_uniform(urot(2:3),2)
               call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,0.0_dblprec,urot,zero)
               emom(1:3,i,k)=newmom(1:3)
               emomM(1:3,i,k)=emom(1:3,i,k)*mmom(i,k)
            enddo
         enddo
!#ifdef VSL
!         !$omp end parallel do
!#endif
!      else
!         do i=1,Natom
!            do k=1,mensemble
!               call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,0.0_dblprec,0.0_dblprec,grot)
!               emom(1:3,i,k)=newmom(1:3)
!               emomM(1:3,i,k)=emom(1:3,i,k)*mmom(i,k)
!            enddo
!         enddo
!      end if


   end subroutine randomize_spins

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
         ind_nlist,ind_mom_flag,max_no_neigh_ind,sus_ind,ind_list_full,&
         do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro)

      use Constants
      use InputData, only: ind_mom_type
      use QHB, only : mix_beta, do_qhb_mix

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
      integer, intent(in) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments
      integer, dimension(Natom), intent(in) :: ind_list_full
      integer, dimension(Natom), intent(in) :: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      character(len=1), intent(in) :: ind_mom_flag
      integer, intent(in) :: do_dip
      integer, intent(in) :: Num_macro
      integer, intent(in) :: icell
      real(dblprec), intent(in) :: macro_mag_trial
      real(dblprec), dimension(3), intent(in) :: macro_trial
      real(dblprec), dimension(Num_macro,Mensemble), intent(inout) :: mmom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro

      integer :: k !< Current ensemble
      integer :: ind_neigh,curr_ind,fix_neigh
      real(dblprec) :: beta,des,lmetric
      real(dblprec), dimension(3) :: ave_mom

      beta=1.0_dblprec/k_bolt/(temprescale*Temperature+1.0d-15)
      
      !< 2022-11 New feature: Mixing scheme
      if(do_qhb_mix=='Y'.and.de>0.0_dblprec) then
         call mix_beta(Temperature,temprescale,de,beta)
      endif
      

      if (do_lsf=='N'.and.ind_mom_flag=='N') then
         if(de<=0.0_dblprec .or. flipprob<exp(-beta*de)) then
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
            if (do_dip==2) then
               emomM_macro(:,icell,k)=macro_trial(:)
               emom_macro(:,icell,k)=macro_trial(:)/macro_mag_trial
               mmom_macro(icell,k)=macro_mag_trial
            endif
         endif
         !Induced treatment v1 (Bergman): Flip both fixed and induced moments, renormalize induced magnitudes from susceptibility
         else if (do_lsf=='N'.and.ind_mom_flag=='Y'.and.ind_mom_type==1) then
            if(de<=0.0_dblprec .or. flipprob<exp(-beta*de)) then
   !           ! If the flip is accepted then the magnitude of the induced moments must change
               emom(:,iflip,k)=newmom(:)
               emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
               if (ind_list_full(iflip).eq.0) then
                  do ind_neigh=1,ind_nlistsize(iflip)
                     ! Select the current induced moment
                     curr_ind=ind_nlist(ind_neigh,iflip)
                       ! Sum over the nearest neighbours that are fixed
                       ave_mom=0.0_dblprec
                       do fix_neigh=1, ind_nlistsize(curr_ind)
                          ave_mom(1)=ave_mom(1)+emomM(1,ind_nlist(fix_neigh,curr_ind),k)*sus_ind(curr_ind)
                          ave_mom(2)=ave_mom(2)+emomM(2,ind_nlist(fix_neigh,curr_ind),k)*sus_ind(curr_ind)
                          ave_mom(3)=ave_mom(3)+emomM(3,ind_nlist(fix_neigh,curr_ind),k)*sus_ind(curr_ind)
                       enddo
                       ! Vary the magnitude of the induced moments
                       mmom(curr_ind,k)=sqrt(ave_mom(1)*ave_mom(1)+ave_mom(2)*ave_mom(2)+ave_mom(3)*ave_mom(3))+1.0e-12_dblprec
                       emomM(:,curr_ind,k)=emom(:,curr_ind,k)*mmom(curr_ind,k)
                       !emomM(:,curr_ind,k)=ave_mom(:)*sus_ind(curr_ind)
                       !mmom(curr_ind,k)=sqrt(emomM(1,curr_ind,k)**2+emomM(2,curr_ind,k)**2+emomM(3,curr_ind,k)**2)+1.0e-12_dblprec
                       !emom(:,curr_ind,k)=emomM(:,curr_ind,k)/mmom(curr_ind,k)
                  enddo
               end if
            endif
            !Induced treatment v2 (Ebert): Only flip fixed moments, recalculate induced moments from susceptibility
      else if (do_lsf=='N'.and.ind_mom_flag=='Y'.and.ind_mom_type==2) then
      !else if (do_lsf=='N'.and.ind_mom_flag=='Y') then
          if (ind_list_full(iflip).eq.0) then
            if(de<=0.0_dblprec .or. flipprob<exp(-beta*de)) then
               !           ! If the flip is accepted then the magnitude of the induced moments must change
               emom(:,iflip,k)=newmom(:)
               emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
               do ind_neigh=1,ind_nlistsize(iflip)
                  ! Select the current induced moment
                  curr_ind=ind_nlist(ind_neigh,iflip)
                  ! Sum over the nearest neighbours that are fixed
                  ave_mom=0.0_dblprec
                  do fix_neigh=1, ind_nlistsize(curr_ind)
                     ave_mom(1)=ave_mom(1)+emomM(1,ind_nlist(fix_neigh,curr_ind),k)!*sus_ind(curr_ind)
                     ave_mom(2)=ave_mom(2)+emomM(2,ind_nlist(fix_neigh,curr_ind),k)!*sus_ind(curr_ind)
                     ave_mom(3)=ave_mom(3)+emomM(3,ind_nlist(fix_neigh,curr_ind),k)!*sus_ind(curr_ind)
                  enddo
                  ! Align induced moments according to susceptibility
                  emomM(:,curr_ind,k)=ave_mom(:)*sus_ind(curr_ind)
                  mmom(curr_ind,k)=sqrt(emomM(1,curr_ind,k)**2+emomM(2,curr_ind,k)**2+emomM(3,curr_ind,k)**2)+1.0e-12_dblprec
                  emom(:,curr_ind,k)=emomM(:,curr_ind,k)/mmom(curr_ind,k)
               enddo
            end if
          endif
      else
         if(lsf_metric==1) then
            lmetric=1.0_dblprec
         elseif(lsf_metric==2) then
            lmetric=(mmom(iflip,k)/newmmom)**2
         else
            lmetric=(mmom(iflip,k)/newmmom)
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
         ind_nlist,ind_mom_flag,max_no_neigh,sus_ind,&
         do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro)

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
      integer, intent(in) :: do_dip
      integer, intent(in) :: Num_macro
      integer, intent(in) :: icell
      real(dblprec), intent(in) :: macro_mag_trial
      real(dblprec), dimension(3), intent(in) :: macro_trial
      real(dblprec), dimension(Num_macro,Mensemble), intent(inout) :: mmom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro


      integer :: k !< Current ensemble
      integer :: ind_neigh,curr_ind,fix_neigh
      real(dblprec) :: beta,lmetric
      real(dblprec), dimension(3) :: ave_mom

      beta=1.0_dblprec/k_bolt/(temprescale*Temperature+1.0d-15)
      if (do_lsf=='N'.and.ind_mom_flag=='N') then
         if(flipprob<1.0_dblprec/(1.0_dblprec+exp(beta*de))) then
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
            if (do_dip==2) then
               emomM_macro(:,icell,k)=macro_trial(:)
               emom_macro(:,icell,k)=macro_trial(:)/macro_mag_trial
               mmom_macro(icell,k)=macro_mag_trial
            endif
            ! If induced moments are considered
         endif
      elseif (do_lsf=='N'.and.ind_mom_flag=='Y') then
         if(de<=0.0_dblprec .or. flipprob<exp(-beta*de)) then
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
            ! If the flip is accepted then the magnitude of the induced moments must change
            do ind_neigh=1,ind_nlistsize(iflip)
               ! Select the current induced moment
               curr_ind=ind_nlist(ind_neigh,iflip)
               ! Sum over the nearest neighbours that are fixed
               ave_mom=0.0_dblprec
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
            lmetric=1.0_dblprec
         elseif(lsf_metric==2) then
            lmetric=(mmom(iflip,k)/newmmom)**2
         else
            lmetric=(mmom(iflip,k)/newmmom)
         endif
         if(flipprob<1.0_dblprec/(1.0_dblprec+lmetric*exp(beta*de))) then
            mmom(iflip,k) = newmmom
            emom(:,iflip,k)=newmom(:)
            emomM(:,iflip,k)=mmom(iflip,k)*newmom(:)
         endif
      endif
   end subroutine flip_g

   !> Flip selected spin to the direction of heat bath
   subroutine flip_h(Natom, Mensemble, emom, emomM, newmmom, mmom, iflip,temperature,&
         temprescale,k,flipprob,totfield,q)

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
      real(dblprec), intent(in):: q !< Random number for moment update

      integer :: k !< Current ensemble
      real(dblprec) :: beta,zarg,zctheta,zstheta,zcphi,zsphi,ctheta,stheta,phi!,q(1)
      real(dblprec),dimension(3) :: stemp,zfc

      mmom=newmmom
      beta=1.0_dblprec/k_bolt/(temprescale*Temperature)
      zfc=beta*totfield*mub*mmom
      zarg=sqrt(sum(zfc(:)*zfc(:)))
      zctheta=zfc(3)/zarg
      zstheta=sqrt(1.0_dblprec-zctheta*zctheta)
      !if (zstheta < 1.d-2) then
      if (zstheta < 1.d-6) then
         !zctheta=0.999950_dblprec
         !zctheta=0.999999995_dblprec
         zctheta=0.9999999999995_dblprec
         zstheta=1.d-6
         !zcphi=1.0_dblprec
         !zsphi=0.0_dblprec
         zcphi=zfc(1)/(zarg*zstheta)
         zsphi=zfc(2)/(zarg*zstheta)
      else
         zcphi=zfc(1)/(zarg*zstheta)
         zsphi=zfc(2)/(zarg*zstheta)
      endif
      !ctheta=1.50_dblprec
      !do while (abs(ctheta)>1.0_dblprec)
      !   call rng_uniform(q,1)
      !   ctheta=1.0_dblprec+1.0_dblprec/zarg*log(q(1))   ! Modified Direct Heat Bath
      !enddo
      ctheta=1.0_dblprec+(1.0_dblprec/zarg)*log((1.0_dblprec-exp(-2.0_dblprec*zarg))*q+exp(-2.0_dblprec*zarg)+dbl_tolerance)
      stheta=sqrt(1.0_dblprec-ctheta*ctheta)
      phi=pi*(2.0_dblprec*flipprob-1.0_dblprec)
      stemp(1)=stheta*cos(phi)              ! New spin in direction of the effective field
      stemp(2)=stheta*sin(phi)
      stemp(3)=ctheta
      !--Euler rotation of spin from the Heff frame of ref. to cryst. axis-----
      emom(1,iflip,k)=zcphi*zctheta*stemp(1)-zsphi*stemp(2)+zcphi*zstheta*stemp(3)
      emom(2,iflip,k)=zsphi*zctheta*stemp(1)+zcphi*stemp(2)+zsphi*zstheta*stemp(3)
      emom(3,iflip,k)=     -zstheta*stemp(1)                     +zctheta*stemp(3)
      emomM(:,iflip,k)=mmom*emom(:,iflip,k)
   end subroutine flip_h

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calculate_energy
   !> @brief
   !> Calculate total energy of single spin
   !> @author Lars Bergqvist
   !> Jonathan Chico --------> Adding macrocell dipole-dipole energy calculation
   !-----------------------------------------------------------------------------
   subroutine calculate_energy(Natom, Mensemble, nHam, conf_num, do_dm , do_pd, do_biqdm, do_bq, do_ring, do_chir, do_sa,&
         emomM, emom, mmom, iflip, newmom, extfield, de, k, &
         mult_axis, do_dip,Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,&
         macro_atom_nlist,emomM_macro,icell,macro_mag_trial,macro_trial, exc_inter,do_anisotropy,do_jtensor)

      use Constants, only : mub
      use macrocells, only : calc_trial_macro_mom
      use HamiltonianData, only : ham

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_ring !< Add four-spin ring (4SR) term to Hamiltonian (0/1)      
      integer, intent(in) :: do_chir !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: do_sa   !< Add symmetric anisotropic exchange (SA) term to Hamiltonian (0/1)
      integer, intent(in) :: do_anisotropy
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
      character(len=1), intent(in) :: exc_inter !< Interpolation of Jij between FM/DLM (Y/N)
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro !< The full vector of the macrocell magnetic moment
      integer, intent(out) :: icell
      real(dblprec), intent(out) :: macro_mag_trial
      real(dblprec), dimension(3), intent(out) :: macro_trial
      integer, intent(in) ::  do_jtensor  !<  Use SKKR style exchange tensor (0=off, 1=on)

      !.. Local scalars
      integer :: j,mu,nu, iflip_h
      real(dblprec) :: tt, tta, ttb, e_c, e_t
      real(dblprec) :: bqmdot, excscale
      real(dblprec) :: ringmdotij,ringmdotkl,ringmdotil,ringmdotkj,ringmdotik,ringmdotjl   
      real(dblprec), dimension(3) :: field
      integer :: im1,ip1,im2,ip2

      !.. Local arrays
      real(dblprec), dimension(3) :: beff_t, trialmom

      !.. Executable statements

      ! First calculate effective field
      !beff_t(1) = 0.0_dblprec
      !beff_t(2) = 0.0_dblprec
      !beff_t(3) = 0.0_dblprec
      tt=0.0_dblprec

      e_c=0.0_dblprec
      e_t=0.0_dblprec
      trialmom(:)=newmom(:)*mmom(iflip,k)
      iflip_h=ham%aham(iflip)

      if (do_jtensor==1) then
         do j=1,ham%nlistsize(iflip_h)
            beff_t = beff_t +  ham%j_tens(:,1,j,iflip_h)*emomM(1,ham%nlist(j,iflip),k)  &
                          +  ham%j_tens(:,2,j,iflip_h)*emomM(2,ham%nlist(j,iflip),k)  &
                          +  ham%j_tens(:,3,j,iflip_h)*emomM(3,ham%nlist(j,iflip),k)
         end do
         e_c = e_c - sum(emomM(:,iflip,k)*beff_t)
         e_t = e_t - sum(trialmom(:)*beff_t)
      else
         if (exc_inter=='N') then
            ! Exchange
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
            !$omp simd reduction(+:e_c,e_t)
#endif
            do j=1,ham%nlistsize(iflip_h)
               e_c=e_c-ham%ncoup(j,iflip_h,1)*sum(emomM(:,iflip,k)*emomM(:,ham%nlist(j,iflip),k))
               e_t=e_t-ham%ncoup(j,iflip_h,1)*sum(trialmom(:)*emomM(:,ham%nlist(j,iflip),k))
            end do
         else
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
            !$omp simd private(excscale) reduction(+:e_c,e_t)
#endif
            do j=1,ham%nlistsize(iflip_h)
               excscale=abs(sum(emom(:,ham%nlist(j,iflip),k)*emom(:,iflip,k)))
               e_c=e_c-(excscale*ham%ncoup(j,iflip_h,1)+(1.0_dblprec-excscale)*ham%ncoupD(j,iflip_h,1))*sum(emomM(:,iflip,k)*emomM(:,ham%nlist(j,iflip),k))
               e_t=e_t-(excscale*ham%ncoup(j,iflip_h,1)+(1.0_dblprec-excscale)*ham%ncoupD(j,iflip_h,1))*sum(trialmom(:)*emomM(:,ham%nlist(j,iflip),k))
            end do
         endif
      end if

      ! Anisotropy
      if (do_anisotropy==1) then
         ! Uniaxial anisotropy
         if (ham%taniso(iflip)==1) then
            tta=sum(emomM(:,iflip,k)*ham%eaniso(:,iflip))
            ttb=sum(trialmom(:)*ham%eaniso(:,iflip))
            e_c=e_c+ham%kaniso(1,iflip)*(tta**2)+ham%kaniso(2,iflip)*(tta**4)
            e_t=e_t+ham%kaniso(1,iflip)*(ttb**2)+ham%kaniso(2,iflip)*(ttb**4)
            ! Cubic anisotropy
         elseif (ham%taniso(iflip)==2) then
            e_c=e_c-ham%kaniso(1,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
               emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)- &
               ham%kaniso(2,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
               e_t=e_t-ham%kaniso(1,iflip)*(trialmom(1)**2*trialmom(2)**2+ &
               trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)- &
               ham%kaniso(2,iflip)*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
         endif
         ! When both Cubic and Uniaxial are switched on
         if (ham%taniso(iflip)==7) then
            ! Uniaxial anisotropy
            tta=(emomM(1,iflip,k)*ham%eaniso(1,iflip)+emomM(2,iflip,k)*ham%eaniso(2,iflip)+emomM(3,iflip,k)*ham%eaniso(3,iflip))
            ttb=(trialmom(1)*ham%eaniso(1,iflip)+trialmom(2)*ham%eaniso(2,iflip)+trialmom(3)*ham%eaniso(3,iflip))
            e_c=e_c+ham%kaniso(1,iflip)*(tta**2)+ham%kaniso(2,iflip)*(tta**4)
            e_t=e_t+ham%kaniso(1,iflip)*(ttb**2)+ham%kaniso(2,iflip)*(ttb**4)
            ! Cubic anisotropy
            aw1=ham%kaniso(1,iflip)*ham%sb(iflip)
            aw2=ham%kaniso(2,iflip)*ham%sb(iflip)
            ! Adding up the contributions
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
            if (ham%taniso_diff(iflip)==1) then
               tta=(emomM(1,iflip,k)*ham%eaniso_diff(1,iflip)+emomM(2,iflip,k)*ham%eaniso_diff(2,iflip)+emomM(3,iflip,k)*ham%eaniso_diff(3,iflip))
               ttb=(trialmom(1)*ham%eaniso_diff(1,iflip)+trialmom(2)*ham%eaniso_diff(2,iflip)+trialmom(3)*ham%eaniso_diff(3,iflip))
               e_c=e_c+ham%kaniso_diff(1,iflip)*(tta**2)+ham%kaniso_diff(2,iflip)*(tta**4)
               e_t=e_t+ham%kaniso_diff(1,iflip)*(ttb**2)+ham%kaniso_diff(2,iflip)*(ttb**4)

               ! Cubic anisotropy
            elseif (ham%taniso_diff(iflip)==2) then
               e_c=e_c+ham%kaniso_diff(1,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
               emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+ &
               ham%kaniso_diff(2,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
               e_t=e_t+ham%kaniso_diff(1,iflip)*(trialmom(1)**2*trialmom(2)**2+ &
               trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
               ham%kaniso_diff(2,iflip)*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
            endif
            ! When both Cubic and Uniaxial are switched on
            if (ham%taniso_diff(iflip)==7) then
               ! Uniaxial anisotropy
               tta=(emomM(1,iflip,k)*ham%eaniso_diff(1,iflip)+emomM(2,iflip,k)*ham%eaniso_diff(2,iflip)+emomM(3,iflip,k)*ham%eaniso_diff(3,iflip))
               ttb=(trialmom(1)*ham%eaniso_diff(1,iflip)+trialmom(2)*ham%eaniso_diff(2,iflip)+trialmom(3)*ham%eaniso_diff(3,iflip))
               e_c=e_c+ham%kaniso_diff(1,iflip)*(tta**2)+ham%kaniso_diff(2,iflip)*(tta**4)
               e_t=e_t+ham%kaniso_diff(1,iflip)*(ttb**2)+ham%kaniso_diff(2,iflip)*(ttb**4)

               ! Cubic anisotropy
               aw1=ham%kaniso_diff(1,iflip)*ham%sb_diff(iflip)
               aw2=ham%kaniso_diff(2,iflip)*ham%sb_diff(iflip)

               e_c=e_c+aw1*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
               emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+ &
               aw2*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
               e_t=e_t+aw1*(trialmom(1)**2*trialmom(2)**2+ &
               trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
               aw2*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
            endif
         endif
      endif
      ! DM interaction
      if (do_dm==1) then
         do j=1,ham%dmlistsize(iflip_h)
            e_c=e_c+ham%dm_vect(1,j,iflip_h)*(emomM(2,iflip,k)*emomM(3,ham%dmlist(j,iflip),k)- &
               emom(3,iflip,k)*emomM(2,ham%dmlist(j,iflip),k))+ &
               ham%dm_vect(2,j,iflip_h)*(emomM(3,iflip,k)*emomM(1,ham%dmlist(j,iflip),k)- &
               emomM(1,iflip,k)*emomM(3,ham%dmlist(j,iflip),k))+ &
               ham%dm_vect(3,j,iflip_h)*(emom(1,iflip,k)*emomM(2,ham%dmlist(j,iflip),k)- &
               emomM(2,iflip,k)*emomM(1,ham%dmlist(j,iflip),k))
            e_t=e_t+ham%dm_vect(1,j,iflip_h)*(trialmom(2)*emomM(3,ham%dmlist(j,iflip),k)- &
               trialmom(3)*emomM(2,ham%dmlist(j,iflip),k))+ &
               ham%dm_vect(2,j,iflip_h)*(trialmom(3)*emomM(1,ham%dmlist(j,iflip),k)- &
               trialmom(1)*emomM(3,ham%dmlist(j,iflip),k))+ &
               ham%dm_vect(3,j,iflip_h)*(trialmom(1)*emomM(2,ham%dmlist(j,iflip),k)- &
               trialmom(2)*emomM(1,ham%dmlist(j,iflip),k))

         end do
      end if

      ! SA interaction
      if (do_sa==1) then
         do j=1,ham%salistsize(iflip_h)
            e_c=e_c+ham%sa_vect(1,j,iflip_h)*(emomM(2,iflip,k)*emomM(3,ham%salist(j,iflip),k)+ &
               emom(3,iflip,k)*emomM(2,ham%salist(j,iflip),k))+ &
               ham%sa_vect(2,j,iflip_h)*(emomM(3,iflip,k)*emomM(1,ham%salist(j,iflip),k)+ &
               emomM(1,iflip,k)*emomM(3,ham%salist(j,iflip),k))+ &
               ham%sa_vect(3,j,iflip_h)*(emom(1,iflip,k)*emomM(2,ham%salist(j,iflip),k)+ &
               emomM(2,iflip,k)*emomM(1,ham%salist(j,iflip),k))
            e_c = -e_c
            e_t=e_t+ham%sa_vect(1,j,iflip_h)*(trialmom(2)*emomM(3,ham%salist(j,iflip),k)+ &
               trialmom(3)*emomM(2,ham%salist(j,iflip),k))+ &
               ham%sa_vect(2,j,iflip_h)*(trialmom(3)*emomM(1,ham%salist(j,iflip),k)+ &
               trialmom(1)*emomM(3,ham%salist(j,iflip),k))+ &
               ham%sa_vect(3,j,iflip_h)*(trialmom(1)*emomM(2,ham%salist(j,iflip),k)+ &
               trialmom(2)*emomM(1,ham%salist(j,iflip),k))
            e_t = -e_t

         end do
      end if

      ! PD interaction
      if(do_pd==1) then
         do j=1,ham%pdlistsize(iflip_h)
          !  e_c=e_c-ham%pd_vect(1,j,iflip_h)*emomM(1,iflip,k)*emomM(1,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(4,j,iflip_h)*emomM(1,iflip,k)*emomM(2,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(5,j,iflip_h)*emomM(1,iflip,k)*emomM(3,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(4,j,iflip_h)*emomM(2,iflip,k)*emomM(1,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(2,j,iflip_h)*emomM(2,iflip,k)*emomM(2,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(6,j,iflip_h)*emomM(2,iflip,k)*emomM(3,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(5,j,iflip_h)*emomM(3,iflip,k)*emomM(1,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(6,j,iflip_h)*emomM(3,iflip,k)*emomM(2,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(3,j,iflip_h)*emomM(3,iflip,k)*emomM(3,ham%pdlist(j,iflip),k)

          !  e_t=e_t-ham%pd_vect(1,j,iflip_h)*trialmom(1)*emomM(1,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(4,j,iflip_h)*trialmom(1)*emomM(2,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(5,j,iflip_h)*trialmom(1)*emomM(3,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(4,j,iflip_h)*trialmom(2)*emomM(1,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(2,j,iflip_h)*trialmom(2)*emomM(2,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(6,j,iflip_h)*trialmom(2)*emomM(3,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(5,j,iflip_h)*trialmom(3)*emomM(1,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(6,j,iflip_h)*trialmom(3)*emomM(2,ham%pdlist(j,iflip),k)- &
          !     ham%pd_vect(3,j,iflip_h)*trialmom(3)*emomM(3,ham%pdlist(j,iflip),k)

                field(1)=ham%pd_vect(1,j,iflip_h)*emomM(1,ham%pdlist(j,iflip),k) +&
                         ham%pd_vect(2,j,iflip_h)*emomM(2,ham%pdlist(j,iflip),k) +&
                         ham%pd_vect(3,j,iflip_h)*emomM(3,ham%pdlist(j,iflip),k) 
                field(2)=ham%pd_vect(4,j,iflip_h)*emomM(1,ham%pdlist(j,iflip),k) +&
                         ham%pd_vect(5,j,iflip_h)*emomM(2,ham%pdlist(j,iflip),k) +&
                         ham%pd_vect(6,j,iflip_h)*emomM(3,ham%pdlist(j,iflip),k) 
                field(3)=ham%pd_vect(7,j,iflip_h)*emomM(1,ham%pdlist(j,iflip),k) +&
                         ham%pd_vect(8,j,iflip_h)*emomM(2,ham%pdlist(j,iflip),k) +&
                         ham%pd_vect(9,j,iflip_h)*emomM(3,ham%pdlist(j,iflip),k) 
                           
             e_c=e_c-( emomM(1,iflip,k)*field(1)+emomM(2,iflip,k)*field(2)+emomM(3,iflip,k)*field(3) )
             e_t=e_t-( trialmom(1)*field(1)+trialmom(2)*field(2)+trialmom(3)*field(3) )                   
               
                    
         end do
      end if

      ! BIQDM interaction
      if(do_biqdm==1) then
         do j=1,ham%biqdmlistsize(iflip_h)
            e_c=e_c-ham%biqdm_vect(1,j,iflip_h)*(emomM(2,iflip,k)*emomM(3,ham%biqdmlist(j,iflip),k)- &
               emomM(3,iflip,k)*emomM(2,ham%biqdmlist(j,iflip),k))**2- &
               ham%biqdm_vect(1,j,iflip_h)*(emom(3,iflip,k)*emomM(1,ham%biqdmlist(j,iflip),k)- &
               emomM(1,iflip,k)*emomM(3,ham%biqdmlist(j,iflip),k))**2- &
               ham%biqdm_vect(1,j,iflip_h)*(emom(1,iflip,k)*emomM(2,ham%biqdmlist(j,iflip),k)- &
               emomM(2,iflip,k)*emomM(1,ham%biqdmlist(j,iflip),k))**2

            e_t=e_t-ham%biqdm_vect(1,j,iflip_h)*(trialmom(2)*emomM(3,ham%biqdmlist(j,iflip),k)- &
               trialmom(3)*emomM(2,ham%biqdmlist(j,iflip),k))**2- &
               ham%biqdm_vect(1,j,iflip_h)*(trialmom(3)*emomM(1,ham%biqdmlist(j,iflip),k)- &
               trialmom(1)*emomM(3,ham%biqdmlist(j,iflip),k))**2- &
               ham%biqdm_vect(1,j,iflip_h)*(trialmom(1)*emomM(2,ham%biqdmlist(j,iflip),k)- &
               trialmom(2)*emomM(1,ham%biqdmlist(j,iflip),k))**2
         end do
      end if

      ! BQ interaction
      if(do_bq==1) then

         do j=1,ham%bqlistsize(iflip_h)
            ! current spin
            bqmdot=emomM(1,ham%bqlist(j,iflip),k)*emomM(1,iflip,k)+&
               emomM(2,ham%bqlist(j,iflip),k)*emomM(2,iflip,k)+&
               emomM(3,ham%bqlist(j,iflip),k)*emomM(3,iflip,k)
            e_c=e_c-ham%j_bq(j,iflip_h)*bqmdot**2

            !trial spin
            bqmdot=emomM(1,ham%bqlist(j,iflip),k)*trialmom(1) + &
               emomM(2,ham%bqlist(j,iflip),k)*trialmom(2) + &
               emomM(3,ham%bqlist(j,iflip),k)*trialmom(3)
            e_t=e_t-ham%j_bq(j,iflip_h)*bqmdot**2

         end do
      end if

    ! Four-spin ring interaction
	  if(do_ring==1) then

       do j=1,ham%ringlistsize(iflip_h)
          ! current spin
             ringmdotij=emomM(1,ham%ringlist(iflip,j,1),k)*emomM(1,iflip,k)+&
                       emomM(2,ham%ringlist(iflip,j,1),k)*emomM(2,iflip,k)+&
                       emomM(3,ham%ringlist(iflip,j,1),k)*emomM(3,iflip,k)

             ringmdotkl=emomM(1,ham%ringlist(iflip,j,2),k)*emomM(1,ham%ringlist(iflip,j,3),k)+&
                       emomM(2,ham%ringlist(iflip,j,2),k)*emomM(2,ham%ringlist(iflip,j,3),k)+&
                       emomM(3,ham%ringlist(iflip,j,2),k)*emomM(3,ham%ringlist(iflip,j,3),k)

             ringmdotil=emomM(1,ham%ringlist(iflip,j,3),k)*emomM(1,iflip,k)+&
                       emomM(2,ham%ringlist(iflip,j,3),k)*emomM(2,iflip,k)+&
                       emomM(3,ham%ringlist(iflip,j,3),k)*emomM(3,iflip,k)

             ringmdotkj=emomM(1,ham%ringlist(iflip,j,2),k)*emomM(1,ham%ringlist(iflip,j,1),k)+&
                       emomM(2,ham%ringlist(iflip,j,2),k)*emomM(2,ham%ringlist(iflip,j,1),k)+&
                       emomM(3,ham%ringlist(iflip,j,2),k)*emomM(3,ham%ringlist(iflip,j,1),k)

             ringmdotik=emomM(1,ham%ringlist(iflip,j,2),k)*emomM(1,iflip,k)+&
                       emomM(2,ham%ringlist(iflip,j,2),k)*emomM(2,iflip,k)+&
                       emomM(3,ham%ringlist(iflip,j,2),k)*emomM(3,iflip,k)

             ringmdotjl=emomM(1,ham%ringlist(iflip,j,1),k)*emomM(1,ham%ringlist(iflip,j,3),k)+&
                       emomM(2,ham%ringlist(iflip,j,1),k)*emomM(2,ham%ringlist(iflip,j,3),k)+&
                       emomM(3,ham%ringlist(iflip,j,1),k)*emomM(3,ham%ringlist(iflip,j,3),k)

          e_c=e_c+ham%j_ring(iflip_h,j)*ringmdotij*ringmdotkl+&
             ham%j_ring(iflip_h,j)*ringmdotil*ringmdotkj-ham%j_ring(iflip_h,j)*ringmdotik*ringmdotjl

          !trial spin

             ringmdotij=emomM(1,ham%ringlist(iflip,j,1),k)*trialmom(1)+&
                       emomM(2,ham%ringlist(iflip,j,1),k)*trialmom(2)+&
                       emomM(3,ham%ringlist(iflip,j,1),k)*trialmom(3)

             ringmdotil=emomM(1,ham%ringlist(iflip,j,3),k)*trialmom(1)+&
                       emomM(2,ham%ringlist(iflip,j,3),k)*trialmom(2)+&
                       emomM(3,ham%ringlist(iflip,j,3),k)*trialmom(3)

             ringmdotik=emomM(1,ham%ringlist(iflip,j,2),k)*trialmom(1)+&
                       emomM(2,ham%ringlist(iflip,j,2),k)*trialmom(2)+&
                       emomM(3,ham%ringlist(iflip,j,2),k)*trialmom(3)

          e_t=e_t+ham%j_ring(iflip_h,j)*ringmdotij*ringmdotkl+&
             ham%j_ring(iflip_h,j)*ringmdotil*ringmdotkj-ham%j_ring(iflip_h,j)*ringmdotik*ringmdotjl

       end do
    end if

      ! CHIR interaction
      if (do_chir==1) then
         do j=1,ham%chirlistsize(iflip_h)

            im1=ham%chirlist(2,j,iflip)
            ip1=ham%chirlist(1,j,iflip)
            im2=ham%chirlist(2,j,im1)
            if(im2==0) im2=im1
            ip2=ham%chirlist(1,j,ip1)
            if(ip2==0) ip2=ip1
            !print '(2x,a,2i4,5x,5i4)','->  ', i,j,im2,im1,i,ip1,ip2
            field(1) = field(1)  &
               - ham%chir_coup(j,iflip_h)*emomM(2,ip1,k)*emomM(3,im1,k) + ham%chir_coup(j,iflip_h)*emomM(3,ip1,k)*emomM(2,im1,k) &
               - ham%chir_coup(j,iflip_h)*emomM(2,im2,k)*emomM(3,im1,k) + ham%chir_coup(j,iflip_h)*emomM(3,im2,k)*emomM(2,im1,k) &
               - ham%chir_coup(j,iflip_h)*emomM(2,ip1,k)*emomM(3,ip2,k) + ham%chir_coup(j,iflip_h)*emomM(3,ip1,k)*emomM(2,ip2,k) 

            field(2) = field(2)  &
               - ham%chir_coup(j,iflip_h)*emomM(3,ip1,k)*emomM(1,im1,k) + ham%chir_coup(j,iflip_h)*emomM(1,ip1,k)*emomM(3,im1,k) &
               - ham%chir_coup(j,iflip_h)*emomM(3,im2,k)*emomM(1,im1,k) + ham%chir_coup(j,iflip_h)*emomM(1,im2,k)*emomM(3,im1,k) &
               - ham%chir_coup(j,iflip_h)*emomM(3,ip1,k)*emomM(1,ip2,k) + ham%chir_coup(j,iflip_h)*emomM(1,ip1,k)*emomM(3,ip2,k) 

            field(3) = field(3)  &
               - ham%chir_coup(j,iflip_h)*emomM(1,ip1,k)*emomM(2,im1,k) + ham%chir_coup(j,iflip_h)*emomM(2,ip1,k)*emomM(1,im1,k) &
               - ham%chir_coup(j,iflip_h)*emomM(1,im2,k)*emomM(2,im1,k) + ham%chir_coup(j,iflip_h)*emomM(2,im2,k)*emomM(1,im1,k) &
               - ham%chir_coup(j,iflip_h)*emomM(1,ip1,k)*emomM(2,ip2,k) + ham%chir_coup(j,iflip_h)*emomM(2,ip1,k)*emomM(1,ip2,k) 

            e_c = e_c - emomM(1,iflip,k)*field(1) - emomM(2,iflip,k)*field(2) -emomM(3,iflip,k)*field(3) 
            e_t = e_t - trialmom(1)*field(1) - trialmom(2)*field(2) -trialmom(3)*field(3) 


!            c_vect(1)=emomM(2,ip1,k)*emomM(3,im1,k)-emomM(3,ip1,k)*emomM(2,im1,k)
!            c_vect(2)=emomM(2,ip1,k)*emomM(3,im1,k)-emomM(3,ip1,k)*emomM(2,im1,k)
!            c_vect(3)=emomM(2,ip1,k)*emomM(3,im1,k)-emomM(3,ip1,k)*emomM(2,im1,k)
!
!            b_vect(1)=emomM(2,im2,k)*emomM(3,im1,k)-emomM(3,im2,k)*emomM(2,im1,k)
!            b_vect(2)=emomM(2,im2,k)*emomM(3,im1,k)-emomM(3,im2,k)*emomM(2,im1,k)
!            b_vect(3)=emomM(2,im2,k)*emomM(3,im1,k)-emomM(3,im2,k)*emomM(2,im1,k)
!
!            f_vect(1)=emomM(2,ip1,k)*emomM(3,ip2,k)-emomM(3,ip1,k)*emomM(2,ip2,k)
!            f_vect(2)=emomM(2,ip1,k)*emomM(3,ip2,k)-emomM(3,ip1,k)*emomM(2,ip2,k)
!            f_vect(3)=emomM(2,ip1,k)*emomM(3,ip2,k)-emomM(3,ip1,k)*emomM(2,ip2,k)
!
!            e_c = e_c - ham%chir_coup(j,iflip_h)*emomM(1,iflip,k)*c_vect(1) - ham%chir_coup(j,iflip_h)*emomM(2,iflip,k)*c_vect(2)  - ham%chir_coup(j,iflip_h)*emomM(3,iflip,k)*c_vect(3) &
!                      - ham%chir_coup(j,iflip_h)*emomM(1,iflip,k)*b_vect(1) - ham%chir_coup(j,iflip_h)*emomM(2,iflip,k)*b_vect(2)  - ham%chir_coup(j,iflip_h)*emomM(3,iflip,k)*b_vect(3) &
!                      - ham%chir_coup(j,iflip_h)*emomM(1,iflip,k)*f_vect(1) - ham%chir_coup(j,iflip_h)*emomM(2,iflip,k)*f_vect(2)  - ham%chir_coup(j,iflip_h)*emomM(3,iflip,k)*f_vect(3) 
!
!            e_t = e_t - ham%chir_coup(j,iflip_h)*trialmom(1)*c_vect(1) - ham%chir_coup(j,iflip_h)*trialmom(2)*c_vect(2) - ham%chir_coup(j,iflip_h)*trialmom(3)*c_vect(3) &
!                      - ham%chir_coup(j,iflip_h)*trialmom(1)*b_vect(1) - ham%chir_coup(j,iflip_h)*trialmom(2)*b_vect(2) - ham%chir_coup(j,iflip_h)*trialmom(3)*b_vect(3) &
!                      - ham%chir_coup(j,iflip_h)*trialmom(1)*f_vect(1) - ham%chir_coup(j,iflip_h)*trialmom(2)*f_vect(2) - ham%chir_coup(j,iflip_h)*trialmom(3)*f_vect(3) 
!
!           print '(4x,2i4,5x,2f14.6)',iflip, j ,e_c,e_t
         end do
      end if

      ! Dipole-dipole interaction
      if (do_dip==1) then
         do j=1,Natom
            do mu=1,3
               do nu=1,3
                  e_c=e_c-emomM(mu,iflip,k)*ham%Qdip(nu,mu,j,iflip)*emomM(nu,j,k)
                  e_t=e_t-trialmom(mu)*ham%Qdip(nu,mu,j,iflip)*emomM(nu,j,k)
               enddo
            enddo
         end do
      elseif(do_dip==2) then
         ! Calculation of the trial moment in the macrocell approach
         call calc_trial_macro_mom(k,iflip,Natom,Mensemble,Num_macro,max_num_atom_macro_cell,&
            macro_nlistsize,macro_atom_nlist,trialmom,emomM,emomM_macro,macro_mag_trial,macro_trial)
         icell=cell_index(iflip)
         do j=1, Num_macro
            do mu=1,3
               do nu=1,3
                  e_c=e_c-emomM_macro(mu,icell,k)*ham%Qdip_macro(nu,mu,j,icell)*emomM_macro(nu,j,k)/macro_nlistsize(icell)
                  e_t=e_t-macro_trial(mu)*ham%Qdip_macro(nu,mu,j,icell)*emomM_macro(nu,j,k)/macro_nlistsize(icell)
               enddo
            enddo
         enddo
      endif

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

      newmom(1) = -1.0_dblprec*emom(1,i,en)
      newmom(2) = -1.0_dblprec*emom(2,i,en)
      newmom(3) = -1.0_dblprec*emom(3,i,en)

   end subroutine Ising_random_flip

end module montecarlo_common
