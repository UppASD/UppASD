!-------------------------------------------------------------------------------
! MODULE: Energy
!> @brief Routine for calculating the total energy of the system
!> @author Anders Bergman, Lars Bergqvist, Johan Hellsvik
!> @details The factor fcinv, used in the energy calculation, transforms all the
!> energy contributions to mRy. This is necessary since the parameters needed for the
!> fields are in Teslas.
!> @author
!> Anders Bergman, Johan Hellsvik, Lars Bergqvist, Jonathan Chico
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module Energy
   use Parameters
   use Profiling
   !  use Ewaldmom
   implicit none

   public

contains


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_energy
   !> @brief Calculates the total energy and the term projected contributions to the energy
   !> @todo Check consistency with the effective field calculation in heisge()
   !> @todo Replace use of unit length magnetic vectors emom to full length vectors emomM
   !-----------------------------------------------------------------------------
   subroutine calc_energy(Natom, Nchmax, Mensemble, conf_num, emom, emomM, mmom,simid, &
         plotenergy, mstep, extfield,time_external_fields,tenergy, eenergy, lsfenergy, &
         totene,max_no_neigh, nlistsize, nlist, ncoup, ncoupD,exc_inter, &
         do_dm, max_no_dmneigh, dmlistsize, dmlist, dm_vect, &
         do_pd, nn_pd_tot, pdlistsize, pdlist, pd_vect, &
         do_biqdm, nn_biqdm_tot, biqdmlistsize, biqdmlist, biqdm_vect, &
         do_bq, nn_bq_tot, bqlistsize, bqlist, j_bq, &
         do_dip, qdip, taniso, eaniso, kaniso,sb,&
         mult_axis,taniso_diff, eaniso_diff, kaniso_diff,sb_diff,&
         delta_t,real_time_measure,&
         do_lsf,fs_nlist,fs_nlistsize,nind,lsf_interpolate,lsf_field,Temp)
      !
      use Constants
      use LSF, only : totalenergy_LSF
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: mstep   !< Current simulation step
      integer, intent(in) :: Natom   !< Number of atoms in system
      integer, intent(in) :: Nchmax  !< Number of chemical type
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      integer, intent(in) :: do_biqdm     !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_dmneigh    !< Calculated number of neighbours with DM interactions
      integer, intent(in) :: nn_pd_tot    !< Calculated number of neighbours with PD interactions
      integer, intent(in) :: nn_bq_tot    !< Calculated number of neighbours with BQ interactions
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      integer, intent(in) :: plotenergy   !< Calculate and plot energy (0/1)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      integer, dimension(Natom),intent(in) :: taniso         !< Type of anisotropy (0-2)
      integer, dimension(Natom),intent(in) :: taniso_diff    !< Type of anisotropy (0-2)
      integer, dimension(Natom), intent(in) :: nlistsize     !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom), intent(in) :: dmlistsize    !< Size of neighbour list for DM
      integer, dimension(Natom), intent(in) :: pdlistsize    !< Size of neighbour list for PD
      integer, dimension(Natom), intent(in) :: bqlistsize    !< Size of neighbour list for BQ
      integer, dimension(Natom), intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist         !< List of neighbours for DM
      integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist         !< List of neighbours for PD
      integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist         !< List of neighbours for BQ
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist       !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
      real(dblprec), intent(in) :: delta_t            !< Current time step
      real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between Cubic and Uniaxial anisotropy
      real(dblprec), dimension(Natom),intent(in) :: sb_diff !< Ratio between Cubic and Uniaxial anisotropy
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso      !< Anisotropy constant
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
      real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
      real(dblprec), dimension(3,3,Natom,Natom), intent(in) :: Qdip !< Matrix for dipole-dipole interaction
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom     !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom    !< Current magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: extfield !< External magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: time_external_fields !< Time dependent external magnetic field
      real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup  !< Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
      real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange couplings
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: mult_axis
      character(len=1), intent(in) :: real_time_measure  !< Display measurements in real time
      character(len=1), intent(in)  ::  do_lsf     !< Including LSF energy
      character(len=1),intent(in) :: exc_inter !< Interpolation of Jij between FM/DLM
      integer, dimension(:,:), intent(in) :: fs_nlist !< First shell Neighbouring list for centered atom
      integer, dimension(:), intent(in) :: fs_nlistsize  !< Size of first shell neighbouring list for centered atom
      integer, dimension(:,:), intent(in) :: nind !< index of firstshell-neighbour-list corresponds to neighbour-list
      character(len=1), intent(in) :: lsf_interpolate !< Interpolate LSF or not
      character(len=1), intent(in) :: lsf_field       !< LSF field contribution (Local/Total)
      real(dblprec), intent(in) :: Temp               !< Temperature

      !.. Subroutine output
      real(dblprec), dimension(Mensemble), intent(out) :: tenergy !< Total energy
      real(dblprec), dimension(Mensemble), intent(out) :: eenergy !< Total exchange (transversal) energy
      real(dblprec), dimension(Mensemble), intent(out) :: lsfenergy !< Total LSF (longitudinal) energy
      real(dblprec), intent(out) :: totene  !< Total energy

      real(dblprec) :: xu1,xu2

      !.. Local scalars
      integer :: i, j, k, inttype

      !.. Local arrays
      real(dblprec) :: fcinv,fc
      real(dblprec), dimension(Mensemble) :: tt
      real(dblprec), dimension(Mensemble) :: mdot,dmx,dmy,dmz,bqmdot
      real(dblprec), dimension(Mensemble) :: pdx, pdy, pdz
      real(dblprec), dimension(Mensemble) :: sxy,syz,szx
      real(dblprec), dimension(Mensemble) :: tempk1, tempk2, texc

      real(dblprec), dimension(Mensemble) :: fenergy, dmenergy, pdenergy, biqdmenergy, bqenergy, dipenergy
      real(dblprec), dimension(Mensemble) :: aenergy, aeatom,aeatom2

      real(dblprec) :: tenergym, tenergys, aenergym, aenergys, eenergym, eenergys,lsfenergys,lsfenergym
      real(dblprec) :: dmenergym, dmenergys, pdenergym, pdenergys, biqdmenergym, biqdmenergys
      real(dblprec) :: bqenergym, bqenergys, dipenergym, dipenergys
      real(dblprec) :: fenergym, fenergys, maxfield
      character(len=30) :: filn

      real(dblprec),dimension(3) :: nbsumfield

      real(dblprec) :: excscale !< Interpolation parameter FM/DLM

      !.. Executable statements

      ! Initialize energy variables
      tenergy = 0.0_dblprec       ! Total Energy
      aenergy = 0.0_dblprec       ! Anisotropy Energy
      eenergy = 0.0_dblprec     ! Exchange Energy
      fenergy = 0.0_dblprec     ! Field Energy
      dmenergy = 0.0_dblprec    ! DM Energy
      pdenergy = 0.0_dblprec    ! PD Energy
      biqdmenergy = 0.0_dblprec ! BIQ Energy
      bqenergy = 0.0_dblprec    ! BQ Energy
      dipenergy = 0.0_dblprec   ! Dip Energy
      lsfenergy = 0.0_dblprec   ! LSF Energy
      ! Factor for energy scale
      fcinv = mub/mry
      fc = mry/mub
      inttype=1
      if (lsf_interpolate=='N') then
         inttype=0
      elseif(lsf_interpolate=='H') then
         inttype=2
      endif

      if (do_lsf=='N') then

         ! Loop over atoms
!#ifndef __PATHSCALE__
#if ((! defined  __PATHSCALE__) || (! defined __PGIF90__)) && (!_OPENMP < 201307)  &&  (!__INTEL_COMPILER >= 1800)
         !$omp parallel do default(shared),private(i,j,k,mdot,texc,nbsumfield,maxfield,excscale,tt,dmx,dmy,dmz,pdx,pdy,pdz, &
         !$omp sxy,syz,szx,bqmdot,aeatom,aeatom2,tempk1,tempk2,xu1,xu2), &
         !$omp reduction(+:eenergy,dmenergy,pdenergy,biqdmenergy,bqenergy,dipenergy,aenergy,fenergy),schedule(static)
#endif
         do i=1, Natom
            if (exc_inter=='N') then
               do j=1,nlistsize(i)
                  ! Calculate exchange energy
                  mdot(:)=emomM(1,nlist(j,i),:)*emomM(1,i,:)+&
                     emomM(2,nlist(j,i),:)*emomM(2,i,:)+&
                     emomM(3,nlist(j,i),:)*emomM(3,i,:)
                  eenergy(:) = eenergy(:) - ncoup(j,i,1)*mdot(:)
               end do
            else
               nbsumfield(:)=0._dblprec ; maxfield=0.0_dblprec
               do j=1,nlistsize(i)
                  nbsumfield(1) = nbsumfield(1) + sum(emomM(1,nlist(j,i),:))
                  nbsumfield(2) = nbsumfield(2) + sum(emomM(2,nlist(j,i),:))
                  nbsumfield(3) = nbsumfield(3) + sum(emomM(3,nlist(j,i),:))
                  maxfield   = maxfield   + sum(mmom(nlist(j,i),:))
               enddo
               excscale=sqrt(nbsumfield(1)**2+nbsumfield(2)**2+nbsumfield(3)**2)/maxfield
               do j=1,nlistsize(i)
                  ! Calculate exchange energy
                  texc=excscale*ncoup(j,i,1)+(1.0_dblprec-excscale)*ncoupD(j,i,1)
                  mdot(:)=emomM(1,nlist(j,i),:)*emomM(1,i,:)+&
                     emomM(2,nlist(j,i),:)*emomM(2,i,:)+&
                     emomM(3,nlist(j,i),:)*emomM(3,i,:)
                  eenergy(:) = eenergy(:) - texc*mdot(:)
               end do
            endif

            if(do_dm==1) then
               do j=1,dmlistsize(i)
                  ! Calculate DM energy
                  dmx(:) = emomM(3,dmlist(j,i),:)*emomM(2,i,:)-emomM(2,dmlist(j,i),:)*emomM(3,i,:)
                  dmy(:) = emomM(1,dmlist(j,i),:)*emomM(3,i,:)-emomM(3,dmlist(j,i),:)*emomM(1,i,:)
                  dmz(:) = emomM(2,dmlist(j,i),:)*emomM(1,i,:)-emomM(1,dmlist(j,i),:)*emomM(2,i,:)
                  dmenergy(:) = dmenergy(:) + dmx(:)*dm_vect(1,j,i) + dmy(:)*dm_vect(2,j,i) + dmz(:)*dm_vect(3,j,i)
               end do
            end if
            if(do_pd==1) then
               do j=1,pdlistsize(i)
                  ! Calculate PD energy
                  pdx(:) = pd_vect(1,j,i)*emomM(1,i,:)*emomM(1,pdlist(j,i),:)+&
                     pd_vect(4,j,i)*emomM(1,i,:)*emomM(2,pdlist(j,i),:)+&
                     pd_vect(5,j,i)*emomM(1,i,:)*emomM(3,pdlist(j,i),:)
                  pdy(:) = pd_vect(4,j,i)*emomM(2,i,:)*emomM(1,pdlist(j,i),:)+&
                     pd_vect(2,j,i)*emomM(2,i,:)*emomM(2,pdlist(j,i),:)+&
                     pd_vect(6,j,i)*emomM(2,i,:)*emomM(3,pdlist(j,i),:)
                  pdz(:) = pd_vect(5,j,i)*emomM(3,i,:)*emomM(1,pdlist(j,i),:)+&
                     pd_vect(6,j,i)*emomM(3,i,:)*emomM(2,pdlist(j,i),:)+&
                     pd_vect(3,j,i)*emomM(3,i,:)*emomM(3,pdlist(j,i),:)
                  pdenergy(:) = pdenergy(:) - pdx(:) - pdy(:) - pdz(:)
               end do
            end if
            if(do_biqdm==1) then
               do j=1,biqdmlistsize(i)
                  !Calculate biqDM energy
                  sxy(:) = emomM(1,i,:)*emomM(2,biqdmlist(j,i),:)-&
                     emomM(2,i,:)*emomM(1,biqdmlist(j,i),:)
                  syz(:) = emomM(2,i,:)*emomM(3,biqdmlist(j,i),:)-&
                     emomM(3,i,:)*emomM(2,biqdmlist(j,i),:)
                  szx(:) = emomM(3,i,:)*emomM(1,biqdmlist(j,i),:)-&
                     emomM(1,i,:)*emomM(3,biqdmlist(j,i),:)
                  biqdmenergy(:) = biqdmenergy(:) + biqdm_vect(1,j,i)*&
                     (sxy(:)**2+syz(:)**2+szx(:)**2)
               end do
            end if
            if(do_bq==1) then
               do j=1,bqlistsize(i)
                  ! Calculate BQ energy
                  bqmdot(:)=emomM(1,bqlist(j,i),:)*emomM(1,i,:)+&
                     emomM(2,bqlist(j,i),:)*emomM(2,i,:)+&
                     emomM(3,bqlist(j,i),:)*emomM(3,i,:)
                  bqenergy(:) = bqenergy(:) - j_bq(j,i) * bqmdot(:) * bqmdot(:)
               end do
            end if
            if(do_dip==1) then
               ! Dipolar contribution to the energy
               do j=1,Natom
                  dipenergy(:) = dipenergy(:) - Qdip(1,1,j,i)*emomM(1,i,:)*emomM(1,j,:)
                  dipenergy(:) = dipenergy(:) - Qdip(2,1,j,i)*emomM(1,i,:)*emomM(2,j,:)
                  dipenergy(:) = dipenergy(:) - Qdip(3,1,j,i)*emomM(1,i,:)*emomM(3,j,:)
                  dipenergy(:) = dipenergy(:) - Qdip(1,2,j,i)*emomM(2,i,:)*emomM(1,j,:)
                  dipenergy(:) = dipenergy(:) - Qdip(2,2,j,i)*emomM(2,i,:)*emomM(2,j,:)
                  dipenergy(:) = dipenergy(:) - Qdip(3,2,j,i)*emomM(2,i,:)*emomM(3,j,:)
                  dipenergy(:) = dipenergy(:) - Qdip(1,3,j,i)*emomM(3,i,:)*emomM(1,j,:)
                  dipenergy(:) = dipenergy(:) - Qdip(2,3,j,i)*emomM(3,i,:)*emomM(2,j,:)
                  dipenergy(:) = dipenergy(:) - Qdip(3,3,j,i)*emomM(3,i,:)*emomM(3,j,:)
               end do
            end if

            if (taniso(i)==1) then

               ! Calculate uniaxial anisotropy energy
               tt(:) = eaniso(1,i)*emomM(1,i,:)+eaniso(2,i)*emomM(2,i,:)+eaniso(3,i)*emomM(3,i,:)
               aeatom(:) = (kaniso(1,i)*tt(:)**2) + kaniso(2,i)*(tt(:)**2)**2

               aenergy(:) = aenergy(:)+aeatom(:)

            elseif (taniso(i)==2) then

               ! Calculate cubic anisotropy energy
               tempk1(:) = emomM(1,i,:)**2*emomM(2,i,:)**2 + emomM(2,i,:)**2*emomM(3,i,:)**2 +&
                  emomM(3,i,:)**2*emomM(1,i,:)**2
               tempk2(:) = emomM(1,i,:)**2 * emomM(2,i,:)**2 * emomM(3,i,:)**2

               aeatom(:) = -(kaniso(1,i)*tempk1(:) + kaniso(2,i)*tempk2(:))

               aenergy(:) = aenergy(:)+aeatom(:)

               !!! ! ------------------------------------ For Unixial+Cubic -------------------------------------------------
            elseif (taniso(i)==7) then

               ! Calculate uniaxial anisotropy energy
               tt(:) = eaniso(1,i)*emomM(1,i,:)+eaniso(2,i)*emomM(2,i,:)+eaniso(3,i)*emomM(3,i,:)
               aeatom(:) = -(kaniso(1,i)*(1-tt(:)**2) - kaniso(2,i)*(1-tt(:)**2)**2)
               aenergy(:) = aenergy(:)+aeatom(:)

               ! Calculate cubic anisotropy energy
               xu1=kaniso(1,i)*sb(i)
               xu2=kaniso(2,i)*sb(i)
               tempk1(:) = emomM(1,i,:)**2*emomM(2,i,:)**2 + emomM(2,i,:)**2*emomM(3,i,:)**2 +&
                  emomM(3,i,:)**2*emomM(1,i,:)**2
               tempk2(:) = emomM(1,i,:)**2 * emomM(2,i,:)**2 * emomM(3,i,:)**2
               aeatom(:) = -(xu1*tempk1(:) + xu2*tempk2(:))

            endif

            ! ------------------------------------ Energy Unixial+Cubic -------------------------------------------------$
            if (mult_axis=='Y') then
               if (taniso_diff(i)==1) then

                  ! Calculate uniaxial anisotropy energy
                  tt(:) = eaniso_diff(1,i)*emomM(1,i,:)+eaniso_diff(2,i)*emomM(2,i,:)+eaniso_diff(3,i)*emomM(3,i,:)
                  aeatom2(:)  = (kaniso_diff(1,i)*tt(:)**2)+kaniso_diff(2,i)*(tt(:)**2)**2
                  aenergy(:) = aenergy(:)+aeatom2(:)

               elseif (taniso_diff(i)==2) then

                  ! Calculate cubic anisotropy energy
                  tempk1(:) = emomM(1,i,:)**2*emomM(2,i,:)**2+emomM(2,i,:)**2*emomM(3,i,:)**2+&
                     emomM(3,i,:)**2*emomM(1,i,:)**2
                  tempk2(:) = emomM(1,i,:)**2*emomM(2,i,:)**2*emomM(3,i,:)**2

                  aeatom2(:) = -(kaniso_diff(1,i)*tempk1(:)+kaniso_diff(2,i)*tempk2(:))

                  !!!! ------------------------------------ For Unixial+Cubic -------------------------------------------------
               elseif (taniso_diff(i)==7) then
                  ! Calculate uniaxial anisotropy energy
                  tt(:) = eaniso_diff(1,i)*emomM(1,i,:)+eaniso_diff(2,i)*emomM(2,i,:)+eaniso_diff(3,i)*emomM(3,i,:)
                  aeatom2(:)  = -(kaniso_diff(1,i)*(1-tt(:)**2)-kaniso_diff(2,i)*(1-tt(:)**2)**2)
                  aenergy(:) = aenergy(:)+aeatom2(:)

                  ! Calculate cubic anisotropy energy
                  xu1=kaniso_diff(1,i)*sb_diff(i)
                  xu2=kaniso_diff(2,i)*sb_diff(i)
                  tempk1(:) = emomM(1,i,:)**2*emomM(2,i,:)**2+emomM(2,i,:)**2*emomM(3,i,:)**2+&
                     emomM(3,i,:)**2*emomM(1,i,:)**2
                  tempk2(:) = emomM(1,i,:)**2*emomM(2,i,:)**2*emomM(3,i,:)**2
                  aeatom(:) = -(xu1*tempk1(:) + xu2*tempk2(:))
                  aenergy(:) = aenergy(:)+aeatom(:)
               endif

            endif

            ! External field energy
            do k=1,Mensemble
               fenergy(k)=fenergy(k)-(extfield(1,i,k)+time_external_fields(1,i,k))*emomM(1,i,k)&
                  -(extfield(2,i,k)+time_external_fields(2,i,k))*emomM(2,i,k)&
                  -(extfield(3,i,k)+time_external_fields(3,i,k))*emomM(3,i,k)
            end do
         end do    ! do atom
!#ifndef __PATHSCALE__
#if ((! defined  __PATHSCALE__) || (! defined __PGIF90__)) && (!_OPENMP < 201307)  &&  (!__INTEL_COMPILER >= 1800)
         !$omp end parallel do
#endif
      else   !do_lsf
         call totalenergy_LSF(Natom, Nchmax, Mensemble, conf_num, emom, emomM, mmom,simid, &
            plotenergy, mstep, extfield,eenergy, aenergy, fenergy, lsfenergy, &
            max_no_neigh, nlistsize, nlist, ncoup, ncoupD,exc_inter, &
            taniso, eaniso, kaniso,sb,do_lsf,fs_nlist,fs_nlistsize,nind,inttype,lsf_field,Temp)
      endif

      ! Divide to get energy per atom
      eenergy(:)=0.5_dblprec*eenergy(:)/Natom
      lsfenergy(:)=0.5_dblprec*lsfenergy(:)/Natom
      fenergy(:)=fenergy(:)/Natom
      aenergy(:)=aenergy(:)/Natom
      dmenergy(:)=0.5_dblprec*dmenergy(:)/Natom
      pdenergy(:)=0.5_dblprec*pdenergy(:)/Natom
      biqdmenergy(:)=0.5_dblprec*biqdmenergy(:)/Natom
      bqenergy(:)=0.5_dblprec*bqenergy(:)/Natom
      dipenergy(:)=0.5_dblprec*dipenergy(:)/Natom
      tenergy(:)=eenergy(:)+fenergy(:)+aenergy(:)+dmenergy(:)+pdenergy(:)+biqdmenergy(:)+bqenergy(:)+dipenergy(:)+lsfenergy(:)

      ! Mean and std.dev. of  energies
      call calculate_mean_and_deviation(tenergy,Mensemble,tenergym,tenergys,fcinv)
      call calculate_mean_and_deviation(eenergy,Mensemble,eenergym,eenergys,fcinv)
      call calculate_mean_and_deviation(aenergy,Mensemble,aenergym,aenergys,fcinv)
      call calculate_mean_and_deviation(dmenergy,Mensemble,dmenergym,dmenergys,fcinv)
      call calculate_mean_and_deviation(pdenergy,Mensemble,pdenergym,pdenergys,fcinv)
      call calculate_mean_and_deviation(biqdmenergy,Mensemble,biqdmenergym,biqdmenergys,fcinv)
      call calculate_mean_and_deviation(bqenergy,Mensemble,bqenergym,bqenergys,fcinv)
      call calculate_mean_and_deviation(dipenergy,Mensemble,dipenergym,dipenergys,fcinv)
      call calculate_mean_and_deviation(fenergy,Mensemble,fenergym,fenergys,fcinv)
      call calculate_mean_and_deviation(lsfenergy,Mensemble,lsfenergym,lsfenergys,fcinv)

      ! Rescale energies for other use later (Cv)
      tenergy=tenergy*fcinv
      eenergy=eenergy*fcinv
      lsfenergy=lsfenergy*fcinv
      totene=tenergym


      ! Print to files
      write (filn,'(''totenergy.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")

      if (real_time_measure=='Y') then
         if (mstep-1==0) then
            write(ofileno,'(a)') ' # Time   Total Energy    Stdv. Tot Ene   Exchange       &
               &    Stdv. Exchange  Ani. Energy     Stdv. Ani       DM Ene         &
               &    Stdv DM Ene     PD. Ene         Stdv. PD Ene    BiqDM Ene      &
               &    Stdv BiqDM      BQ Ene          Stdv BQ Ene     Dipolar Ene    &
               &    Stdv Dipolar    Ext. Field Ene  Stdv Ext. Field  LSF Ene   Stdv LSF   '
         endif
         write(ofileno,10005) (mstep-1)*delta_t, tenergym, tenergys, eenergym, eenergys, aenergym, aenergys, &
            dmenergym, dmenergys, pdenergym, pdenergys, biqdmenergym, biqdmenergys, &
            bqenergym, bqenergys, dipenergym, dipenergys, &
            fenergym, fenergys,lsfenergym,lsfenergys
      else
         if (mstep-1==0) then
            write(ofileno,'(a)') ' # Iter.  Total Energy    Stdv. Tot Ene   Exchange       &
               &  Stdv. Exchange  Ani. Energy     Stdv. Ani       DM Ene         &
               &  Stdv DM Ene     PD. Ene         Stdv. PD Ene    BiqDM Ene      &
               &  Stdv BiqDM      BQ Ene          Stdv BQ Ene     Dipolar Ene    &
               &  Stdv Dipolar    Ext. Field Ene  Stdv Ext. Field  LSF Ene   Stdv LSF '
         endif
         write(ofileno,10004) mstep-1, tenergym, tenergys, eenergym, eenergys, aenergym, aenergys, &
            dmenergym, dmenergys, pdenergym, pdenergys, biqdmenergym, biqdmenergys, &
            bqenergym, bqenergys, dipenergym, dipenergys, &
            fenergym, fenergys,lsfenergym,lsfenergys
      endif
      close(ofileno)

      10004 format (i8,31es16.8)
      10005 format (es12.4,23es16.8)

   end subroutine calc_energy

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_en
   !> @brief Simplified version of calc_energy subroutine
   !> @todo Check consistency with the effective field calculation in heisge()
   !> @todo Replace use of unit length magnetic vectors emom to full length vectors emomM
   !-----------------------------------------------------------------------------
   subroutine calc_en(Natom, Mensemble, emomM, &
         extfield,tenergy, &
         max_no_neigh, nlistsize, nlist, ncoup, &
         do_dm, max_no_dmneigh, dmlistsize, dmlist, dm_vect, &
         do_pd, nn_pd_tot, pdlistsize, pdlist, pd_vect, &
         do_biqdm, nn_biqdm_tot, biqdmlistsize, biqdmlist, biqdm_vect, &
         do_bq, nn_bq_tot, bqlistsize, bqlist, j_bq, &
         do_dip, qdip, taniso, eaniso, kaniso,sb)
      !
      use Constants
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)

      integer, intent(in) :: Natom   !< Number of atoms in system
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      integer, intent(in) :: do_biqdm     !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_dmneigh    !< Calculated number of neighbours with DM interactions
      integer, intent(in) :: nn_pd_tot    !< Calculated number of neighbours with PD interactions
      integer, intent(in) :: nn_bq_tot    !< Calculated number of neighbours with BQ interactions
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      integer, dimension(Natom),intent(in) :: taniso         !< Type of anisotropy (0-2)
      integer, dimension(Natom), intent(in) :: nlistsize     !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom), intent(in) :: dmlistsize    !< Size of neighbour list for DM
      integer, dimension(Natom), intent(in) :: pdlistsize    !< Size of neighbour list for PD
      integer, dimension(Natom), intent(in) :: bqlistsize    !< Size of neighbour list for BQ
      integer, dimension(Natom), intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist         !< List of neighbours for DM
      integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist         !< List of neighbours for PD
      integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist         !< List of neighbours for BQ
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist       !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
      real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between Cubic and Uniaxial anisotropy
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
      real(dblprec), dimension(3,3,Natom,Natom), intent(in) :: Qdip !< Matrix for dipole-dipole interaction
      real(dblprec), dimension(3,Natom,Mensemble) :: emom     !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: extfield !< External magnetic field

      real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
      real(dblprec), dimension(max_no_neigh,Natom), intent(in) :: ncoup  !< Heisenberg exchange couplings
      real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange couplings

      !.. Subroutine output
      real(dblprec), dimension(Mensemble), intent(out) :: tenergy !< Total energy


      real(dblprec) :: xu1,xu2

      !.. Local scalars
      integer :: i, j, k

      !.. Local arrays
      real(dblprec) :: fcinv,tmp
      real(dblprec), dimension(Mensemble) :: tt,tt2,tt3,ttx,tty,ttz
      real(dblprec), dimension(Mensemble) :: mdot,dmx,dmy,dmz,bqmdot
      real(dblprec), dimension(Mensemble) :: pdx, pdy, pdz
      real(dblprec), dimension(Mensemble) :: sxy,syz,szx
      real(dblprec), dimension(Mensemble) :: tempk1, tempk2


      real(dblprec), dimension(Mensemble) :: eenergy, fenergy, dmenergy, pdenergy, biqdmenergy, bqenergy, dipenergy
      real(dblprec), dimension(Mensemble) :: aenergy, aeatom
      real(dblprec) :: tenergym

      do i=1,Mensemble
         do j=1,Natom
            tmp = 0_dblprec
            do k=1,3
               tmp = tmp+emomM(k,j,i)*emomM(k,j,i)
            end do
            tmp = sqrt(tmp)
            emom(1,j,i) = emomM(1,j,i)/tmp
            emom(2,j,i) = emomM(2,j,i)/tmp
            emom(3,j,i) = emomM(3,j,i)/tmp
         end do
      end do


      ! Initialize energy variables
      tenergy(1:Mensemble) = 0_dblprec       ! Total Energy
      aenergy(1:Mensemble) = 0_dblprec       ! Anisotropy Energy
      eenergy(1:Mensemble) = 0.0_dblprec     ! Exchange Energy
      fenergy(1:Mensemble) = 0.0_dblprec     ! Field Energy
      dmenergy(1:Mensemble) = 0.0_dblprec    ! DM Energy
      pdenergy(1:Mensemble) = 0.0_dblprec    ! PD Energy
      biqdmenergy(1:Mensemble) = 0.0_dblprec ! BIQ Energy
      bqenergy(1:Mensemble) = 0.0_dblprec    ! BQ Energy
      dipenergy(1:Mensemble) = 0.0_dblprec   ! Dip Energy

      ! Factor for energy scale
      fcinv = mub/mry

      ! Loop over atoms
      do i=1, Natom
         ttx(:)=0_dblprec
         tty(:)=0_dblprec
         ttz(:)=0_dblprec
         do j=1,nlistsize(i)
            ! Calculate exchange energy
            mdot(:)=emomM(1,nlist(j,i),:)*emomM(1,i,:)+&
               emomM(2,nlist(j,i),:)*emomM(2,i,:)+&
               emomM(3,nlist(j,i),:)*emomM(3,i,:)
            eenergy(:) = eenergy(:) - ncoup(j,i)*mdot(:)
            ! Calculate exchange field
            ttx(:) = ttx(:) + ncoup(j,i)*emomM(1,nlist(j,i),:)
            tty(:) = tty(:) + ncoup(j,i)*emomM(2,nlist(j,i),:)
            ttz(:) = ttz(:) + ncoup(j,i)*emomM(3,nlist(j,i),:)
         end do
         if(do_dm==1) then
            do j=1,dmlistsize(i)
               ! Calculate DM energy
               dmx(:) = emomM(3,dmlist(j,i),:)*emomM(2,i,:)-emomM(2,dmlist(j,i),:)*emomM(3,i,:)
               dmy(:) = emomM(1,dmlist(j,i),:)*emomM(3,i,:)-emomM(3,dmlist(j,i),:)*emomM(1,i,:)
               dmz(:) = emomM(2,dmlist(j,i),:)*emomM(1,i,:)-emomM(1,dmlist(j,i),:)*emomM(2,i,:)
               dmenergy(:) = dmenergy(:) + dmx(:)*dm_vect(1,j,i) + dmy(:)*dm_vect(2,j,i) + dmz(:)*dm_vect(3,j,i)
               !Calculate DM field
               ttx(:) = ttx(:) - dm_vect(3,j,i)*emomM(2,dmlist(j,i),:) +&
                  dm_vect(2,j,i)*emomM(3,dmlist(j,i),:)
               tty(:) = tty(:) - dm_vect(1,j,i)*emomM(3,dmlist(j,i),:) +&
                  dm_vect(3,j,i)*emomM(1,dmlist(j,i),:)
               ttz(:) = ttz(:) - dm_vect(2,j,i)*emomM(1,dmlist(j,i),:) +&
                  dm_vect(1,j,i)*emomM(2,dmlist(j,i),:)
            end do
         end if
         if(do_pd==1) then
            do j=1,pdlistsize(i)
               ! Calculate PD energy
               pdx(:) = pd_vect(1,j,i)*emomM(1,i,:)*emomM(1,pdlist(j,i),:)+&
                  pd_vect(4,j,i)*emomM(1,i,:)*emomM(2,pdlist(j,i),:)+&
                  pd_vect(5,j,i)*emomM(1,i,:)*emomM(3,pdlist(j,i),:)
               pdy(:) = pd_vect(4,j,i)*emomM(2,i,:)*emomM(1,pdlist(j,i),:)+&
                  pd_vect(2,j,i)*emomM(2,i,:)*emomM(2,pdlist(j,i),:)+&
                  pd_vect(6,j,i)*emomM(2,i,:)*emomM(3,pdlist(j,i),:)
               pdz(:) = pd_vect(5,j,i)*emomM(3,i,:)*emomM(1,pdlist(j,i),:)+&
                  pd_vect(6,j,i)*emomM(3,i,:)*emomM(2,pdlist(j,i),:)+&
                  pd_vect(3,j,i)*emomM(3,i,:)*emomM(3,pdlist(j,i),:)
               pdenergy(:) = pdenergy(:) - pdx(:) - pdy(:) - pdz(:)
               ! Calculate PD field
               ttx(:) = ttx(:) + pd_vect(1,j,i)*emomM(1,pdlist(j,i),:) +&
                  pd_vect(4,j,i)*emomM(2,pdlist(j,i),:) +&
                  pd_vect(5,j,i)*emomM(3,pdlist(j,i),:)
               tty(:) = tty(:) + pd_vect(4,j,i)*emomM(1,pdlist(j,i),:) +&
                  pd_vect(2,j,i)*emomM(2,pdlist(j,i),:) +&
                  pd_vect(6,j,i)*emomM(3,pdlist(j,i),:)
               ttz(:) = ttz(:) + pd_vect(5,j,i)*emomM(1,pdlist(j,i),:) +&
                  pd_vect(6,j,i)*emomM(2,pdlist(j,i),:) +&
                  pd_vect(3,j,i)*emomM(3,pdlist(j,i),:)
            end do
         end if
         if(do_biqdm==1) then
            do j=1,biqdmlistsize(i)
               !Calculate biqDM energy
               sxy(:) = emomM(1,i,:)*emomM(2,biqdmlist(j,i),:)-&
                  emomM(2,i,:)*emomM(1,biqdmlist(j,i),:)
               syz(:) = emomM(2,i,:)*emomM(3,biqdmlist(j,i),:)-&
                  emomM(3,i,:)*emomM(2,biqdmlist(j,i),:)
               szx(:) = emomM(3,i,:)*emomM(1,biqdmlist(j,i),:)-&
                  emomM(1,i,:)*emomM(3,biqdmlist(j,i),:)
               biqdmenergy(:) = biqdmenergy(:) + biqdm_vect(1,j,i)*&
                  (sxy(:)**2+syz(:)**2+szx(:)**2)
               !Calculate biqDM field
               ttx(:) = ttx(:) + 2.0_dblprec*biqdm_vect(1,j,i)*(&
                  szx(:)*emomM(3,biqdmlist(j,i),:)-&
                  sxy(:)*emomM(2,biqdmlist(j,i),:))
               tty(:) = tty(:) + 2.0_dblprec*biqdm_vect(1,j,i)*(&
                  sxy(:)*emomM(1,biqdmlist(j,i),:)-&
                  syz(:)*emomM(3,biqdmlist(j,i),:))
               ttz(:) = ttz(:) + 2.0_dblprec*biqdm_vect(1,j,i)*(&
                  syz(:)*emomM(2,biqdmlist(j,i),:)-&
                  szx(:)*emomM(1,biqdmlist(j,i),:))
            end do
         end if
         if(do_bq==1) then
            do j=1,bqlistsize(i)
               ! Calculate BQ energy
               bqmdot(:)=emomM(1,bqlist(j,i),:)*emomM(1,i,:)+&
                  emomM(2,bqlist(j,i),:)*emomM(2,i,:)+&
                  emomM(3,bqlist(j,i),:)*emomM(3,i,:)
               bqenergy(:) = bqenergy(:) - j_bq(j,i) * bqmdot(:) * bqmdot(:)
               ! Calculate biquadratic exchange field
               ttx(:) = ttx(:) + 2.0_dblprec*j_bq(j,i)*bqmdot(:)*emomM(1,bqlist(j,i),:)
               tty(:) = tty(:) + 2.0_dblprec*j_bq(j,i)*bqmdot(:)*emomM(2,bqlist(j,i),:)
               ttz(:) = ttz(:) + 2.0_dblprec*j_bq(j,i)*bqmdot(:)*emomM(3,bqlist(j,i),:)
            end do
         end if
         if(do_dip==1) then
            ! Dipolar contribution to the energy
            do j=1,Natom
               dipenergy(:) = dipenergy(:) - Qdip(1,1,j,i)*emomM(1,i,:)*emomM(1,j,:)
               dipenergy(:) = dipenergy(:) - Qdip(2,1,j,i)*emomM(1,i,:)*emomM(2,j,:)
               dipenergy(:) = dipenergy(:) - Qdip(3,1,j,i)*emomM(1,i,:)*emomM(3,j,:)
               dipenergy(:) = dipenergy(:) - Qdip(1,2,j,i)*emomM(2,i,:)*emomM(1,j,:)
               dipenergy(:) = dipenergy(:) - Qdip(2,2,j,i)*emomM(2,i,:)*emomM(2,j,:)
               dipenergy(:) = dipenergy(:) - Qdip(3,2,j,i)*emomM(2,i,:)*emomM(3,j,:)
               dipenergy(:) = dipenergy(:) - Qdip(1,3,j,i)*emomM(3,i,:)*emomM(1,j,:)
               dipenergy(:) = dipenergy(:) - Qdip(2,3,j,i)*emomM(3,i,:)*emomM(2,j,:)
               dipenergy(:) = dipenergy(:) - Qdip(3,3,j,i)*emomM(3,i,:)*emomM(3,j,:)
            end do
         end if

         if (taniso(i)==1) then

            ! Calculate uniaxial anisotropy energy
            tt(:) = eaniso(1,i)*emomM(1,i,:)+eaniso(2,i)*emomM(2,i,:)+eaniso(3,i)*emomM(3,i,:)
            aeatom(:) = (kaniso(1,i)*tt(:)**2) + kaniso(2,i)*(tt(:)**2)**2

            aenergy(:) = aenergy(:)+aeatom(:)

            ! Calculate uniaxial anisotropy field
            tt(:) = emomM(1,i,:)*eaniso(1,i)+emomM(2,i,:)*eaniso(2,i)+emomM(3,i,:)*eaniso(3,i)

            tt2(:) = 2*tt(:)
            tt3(:) = kaniso(1,i)+2*kaniso(2,i)*(1-tt(:)*tt(:))
            tt3(:) = tt2(:)*tt3(:)

            ttx(:) = ttx(:) - tt3(:)*eaniso(1,i)
            tty(:) = tty(:) - tt3(:)*eaniso(2,i)
            ttz(:) = ttz(:) - tt3(:)*eaniso(3,i)

         elseif (taniso(i)==2) then

            ! Calculate cubic anisotropy energy
            tempk1(:) = emomM(1,i,:)**2*emomM(2,i,:)**2 + emomM(2,i,:)**2*emomM(3,i,:)**2 +&
               emomM(3,i,:)**2*emomM(1,i,:)**2
            tempk2(:) = emomM(1,i,:)**2 * emomM(2,i,:)**2 * emomM(3,i,:)**2

            aeatom(:) = -(kaniso(1,i)*tempk1(:) + kaniso(2,i)*tempk2(:))

            aenergy(:) = aenergy(:)+aeatom(:)

            ! Calculate cubic anisotropy field
            ttx(:) = ttx(:)  &
               + 2*kaniso(1,i)*emomM(1,i,:)*(emomM(2,i,:)**2+emomM(3,i,:)**2) &
               + 2*kaniso(2,i)*emomM(1,i,:)*emomM(2,i,:)**2*emomM(3,i,:)**2
            tty(:) = tty(:)  &
               + 2*kaniso(1,i)*emomM(2,i,:)*(emomM(3,i,:)**2+emomM(1,i,:)**2) &
               + 2*kaniso(2,i)*emomM(2,i,:)*emomM(3,i,:)**2*emomM(1,i,:)**2
            ttz(:) = ttz(:)  &
               + 2*kaniso(1,i)*emomM(3,i,:)*(emomM(1,i,:)**2+emomM(2,i,:)**2) &
               + 2*kaniso(2,i)*emomM(3,i,:)*emomM(1,i,:)**2*emomM(2,i,:)**2

            !!!    ! ------------------------------------ For Unixial+Cubic -------------------------------------------------
         elseif (taniso(i)==7) then

            ! Calculate uniaxial anisotropy energy
            tt(:) = eaniso(1,i)*emomM(1,i,:)+eaniso(2,i)*emomM(2,i,:)+eaniso(3,i)*emomM(3,i,:)
            aeatom(:) = -(kaniso(1,i)*(1-tt(:)**2) - kaniso(2,i)*(1-tt(:)**2)**2)
            aenergy(:) = aenergy(:)+aeatom(:)

            ! Calculate uniaxial anisotropy field
            tt(:) = emomM(1,i,:)*eaniso(1,i)+emomM(2,i,:)*eaniso(2,i)+emomM(3,i,:)*eaniso(3,i)
            tt2(:) = 1.0_dblprec*tt(:)
            tt3(:) = kaniso(1,i)+2*kaniso(2,i)*(tt(:)**2)
            tt3(:) = tt2(:)*tt3(:)
            ttx(:) = ttx(:) - tt3(:)*eaniso(1,i)
            tty(:) = tty(:) - tt3(:)*eaniso(2,i)
            ttz(:) = ttz(:) - tt3(:)*eaniso(3,i)

            ! Calculate cubic anisotropy energy
            xu1=kaniso(1,i)*sb(i)
            xu2=kaniso(2,i)*sb(i)
            tempk1(:) = emomM(1,i,:)**2*emomM(2,i,:)**2 + emomM(2,i,:)**2*emomM(3,i,:)**2 +&
               emomM(3,i,:)**2*emomM(1,i,:)**2
            tempk2(:) = emomM(1,i,:)**2 * emomM(2,i,:)**2 * emomM(3,i,:)**2
            aeatom(:) = -(xu1*tempk1(:) + xu2*tempk2(:))

            ! Calculate cubic anisotropy field
            ttx(:) = ttx(:)  &
               + 2*xu1*emomM(1,i,:)*(emom(2,i,:)**2+emom(3,i,:)**2) &
               + 2*xu2*emomM(1,i,:)*emom(2,i,:)**2*emom(3,i,:)**2

            tty(:) = tty(:)  &
               + 2*xu1*emomM(2,i,:)*(emomM(3,i,:)**2+emomM(1,i,:)**2) &
               + 2*xu2*emomM(2,i,:)*emomM(3,i,:)**2*emomM(1,i,:)**2

            ttz(:) = ttz(:)  &
               + 2*xu1*emomM(3,i,:)*(emomM(1,i,:)**2+emomM(2,i,:)**2) &
               + 2*xu2*emomM(3,i,:)*emomM(1,i,:)**2*emomM(2,i,:)**2
         endif
         ! ------------------------------------ Energy Unixial+Cubic -------------------------------------------------

         ! External field energy
         do k=1,Mensemble
            fenergy(k)=fenergy(k)-(extfield(1,i,k))*emomM(1,i,k)&
               -(extfield(2,i,k))*emomM(2,i,k)&
               -(extfield(3,i,k))*emomM(3,i,k)
            if (extfield(1,i,k) .ne. 0.0D0 .or. extfield(2,i,k) .ne. 0.0D0 .or. extfield(3,i,k) .ne. 0.0D0 ) then

            endif
         end do

      end do

      ! Divide to get energy per atom
      eenergy(:)=0.5_dblprec*eenergy(:)/Natom
      fenergy(:)=fenergy(:)/Natom
      aenergy(:)=aenergy(:)/Natom
      dmenergy(:)=0.5_dblprec*dmenergy(:)/Natom
      pdenergy(:)=0.5_dblprec*pdenergy(:)/Natom
      biqdmenergy(:)=0.5_dblprec*biqdmenergy(:)/Natom
      bqenergy(:)=0.5_dblprec*bqenergy(:)/Natom
      dipenergy(:)=0.5_dblprec*dipenergy(:)/Natom
      tenergy(:)=eenergy(:)+fenergy(:)+aenergy(:)+dmenergy(:)+pdenergy(:)+biqdmenergy(:)+bqenergy(:)+dipenergy(:)

      ! Mean energies
      tenergym=0_dblprec


      tenergy = tenergy*Natom

   end subroutine calc_en

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calculate_mean_and_deviation
   !> @brief Subroutine to calculate the mean and standard deviation of an array
   !-----------------------------------------------------------------------------
   subroutine calculate_mean_and_deviation(vdata,vlen,mean,deviation,factor)
      !
      implicit none
      !
      integer, intent(in) :: vlen
      real(dblprec), dimension(vlen), intent(in) :: vdata
      real(dblprec), intent(out) :: mean
      real(dblprec), intent(out) :: deviation
      real(dblprec), intent(in) :: factor
      !
      !
      if(vlen>1) then
         mean=sum(vdata)
         mean=mean/(vlen*1.0_dblprec)
         deviation=sum((vdata-mean)**2)
         deviation=deviation/((vlen-1.0_dblprec)*1.0_dblprec)
      else
         mean=vdata(1)
         deviation=0.0_dblprec
      end if
      !
      mean=mean*factor
      deviation=deviation*factor
      !
      return
      !
      !
   end subroutine calculate_mean_and_deviation

   !
end module Energy
