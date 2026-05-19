!-------------------------------------------------------------------------------
!> @brief
!> Intended to use for the treatment of induced magnetic moments
!> @details Routines needed to calculate the effect of induced moments via linear response
!> treatment as outlined by S. Polesya et al. in Phys. Rev. B 82, 214409.
!> In this approach the magnitude and direction of the induced moments is determined
!> by the nearest neighbours fixed moments
!> @author
!> Jonathan Chico
!> @note The "dynamics" of the induced moments are slaved to the ones of the fixed
!> moments, as the former can only move when the later ones do. This is expected to improve
!> the \f$ T_c \f$ of certain systems.
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
!@todo Create a full list of neighbours containing all the fixed moments for all the induced
! moments
! Create a number of induced moments and a pruned list containing only them.
! Of that way the list to be evaluated are much smaller
module InducedMoments

   use Parameters
   use Profiling
   use FieldData, only : beff

   implicit none

   integer, dimension(:,:), allocatable :: nmdim_ind     !< Dimension of neighbour map for the nearest neighbours
   integer, dimension(:,:,:), allocatable :: nm_ind      !< Neighbour map for nearest neighbours
   real(dblprec), dimension(:,:), allocatable :: redcoord_mod  !< Modulus of the coordinates for Heisenberg exchange couplings for the induced moments
   real(dblprec), dimension(:,:,:), allocatable :: ind_redcoord   !< Coordinates for Heisenberg exchange couplings for the induced moments
   real(dblprec), dimension(:,:), allocatable :: emomM_trial !< Array for trial moments
   !
   real(dblprec), dimension(:,:,:), allocatable :: rng_work_u  !< work array for RNG
   real(dblprec), dimension(:,:,:), allocatable :: rng_work_g  !< work array for RNG
   real(dblprec), dimension(:,:,:), allocatable :: mom_work_i  !< work array for moments

   private
   public :: induced_mapping, mc_update_ind_mom, calculate_init_ind_sus
   public :: setup_induced_information
   !!public :: renorm_ncoup_ind,setup_induced_information

contains

   !---------------------------------------------------------------------------
   !> @brief
   !> Wrapper routine for the creation of the needed lists and arrays for the
   !> treatment of induced moments
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------
   subroutine induced_mapping(Natom,NT,NA,N1,N2,N3,sym,max_no_shells,nn,atype,      &
         ind_nlistsize,ind_nlist,do_sortcoup,Nchmax,        &
         do_ralloy,Natom_full,atype_ch,  &
         acellnumb,C1,C2,C3,Bas,BC1,BC2,BC3,ind_tol,redcoord,ind_mom,block_size,    &
         ind_list_full,max_no_neigh_ind)

      use NeighbourMap, only : setup_nm
      use HamiltonianData, only: allocate_hamiltoniandata_ind
      use InputData, only : Mensemble

      implicit none
      ! .. Input variables
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: sym !< Symmetry of system (0-3)
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, dimension(:), intent(in):: nn !< Number of neighbour shells
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), optional ,intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced (1) or fixed (0)
      character, intent(in) :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
      character(len=1), intent(in) :: BC1  !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2  !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3  !< Boundary conditions in z-direction
      real(dblprec), intent(in) :: ind_tol !< Value for the tolerance between neighbouring shells
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(3,NA), intent(in) :: Bas !< Coordinates for basis atoms
      real(dblprec), dimension(NT,max_no_shells), intent(in) :: redcoord   !< Coordinates for Heisenberg exchange couplings
      ! .. In/out variables
      integer, dimension(:), allocatable, intent(inout) :: ind_list_full   !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
      integer, dimension(:), allocatable, intent(inout) :: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(:,:), allocatable, intent(inout) :: ind_nlist !< Neighbour list for iduced moments
      !.. Output variables
      integer, intent(inout) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments

      ! .. Local variables
      integer :: i_stat, i_all
      integer :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange

      ! Allocate working variables
      allocate(ind_redcoord(NT,max_no_shells,3),stat=i_stat)
      call memocc(i_stat,product(shape(ind_redcoord))*kind(ind_redcoord),'ind_redcoord','induced_mapping')
      ind_redcoord=0.0_dblprec

      allocate(ind_list_full(Natom),stat=i_stat)
      call memocc(i_stat,product(shape(ind_list_full))*kind(ind_list_full),'ind_list_full','induced_mapping')
      ind_list_full=0

      allocate(emomM_trial(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(emomM_trial))*kind(emomM_trial),'emomM_trial','induced_mapping')
      ind_list_full=0

      allocate(rng_work_g(3,Natom,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(rng_work_g))*kind(rng_work_g),'rng_work_g','induced_mapping')

      allocate(rng_work_u(3,Natom,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(rng_work_u))*kind(rng_work_u),'rng_work_u','induced_mapping')

      allocate(mom_work_i(3,Natom,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(mom_work_i))*kind(mom_work_i),'mom_work_i','induced_mapping')


      ! Create the bonding vector for the nearest neighbours
      call setup_induced_nearest(NT,max_no_shells,ind_tol,redcoord)
      ! Using the redcoord calculated from the first nearest neighbours calculate the
      ! full neighbour list of the system using only nearest neighbours
      call setup_nm(Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,block_size,atype,Bas, &
         max_no_neigh_ind,max_no_shells,max_no_equiv,sym,nn,ind_redcoord,nm_ind,    &
         nmdim_ind,do_ralloy,Natom_full,acellnumb,atype_ch)

      call allocate_hamiltoniandata_ind(1,Natom,max_no_neigh_ind)

      ! After this the 'true' list which discriminates between induced neighbours
      ! and fixed moments must be created
      call setup_induced_list(Natom,NT,NA,atype,max_no_neigh_ind,max_no_equiv,      &
         max_no_shells,do_sortcoup,Natom_full,atype_ch,Nchmax,do_ralloy,nm_ind,nn,  &
         nmdim_ind,ind_nlistsize,ind_nlist,ind_mom)


      ! Deallocate working arrays
      i_all=-product(shape(ind_redcoord))*kind(ind_redcoord)
      deallocate(ind_redcoord,stat=i_stat)
      call memocc(i_stat,i_all,'ind_redcoord','induced_mapping')
      i_all=-product(shape(nm_ind))*kind(nm_ind)
      deallocate(nm_ind,stat=i_stat)
      call memocc(i_stat,i_all,'nm_ind','induced_mapping')
      i_all=-product(shape(nmdim_ind))*kind(nmdim_ind)
      deallocate(nmdim_ind,stat=i_stat)
      call memocc(i_stat,i_all,'nmdim_ind','induced_mapping')

   end subroutine induced_mapping

   !---------------------------------------------------------------------------
   !> @brief
   !> Creation of the bonding vector between neighbours for ony the nearest neighbours
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------
   subroutine setup_induced_nearest(NT,max_no_shells,ind_tol,redcoord)

      implicit none

      ! .. Input variables
      integer, intent(in) :: NT !< Number of types of atoms
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      real(dblprec), intent(in) :: ind_tol !< Value for the tolerance between neighbouring shells
      real(dblprec), dimension(NT,max_no_shells,3), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings

      ! .. Local variables
      integer :: IT, shells, i_stat,i_all
      integer :: temp_counter, temp_max_counter
      real(dblprec) :: temp_min,diff
      real(dblprec), dimension(NT) :: redcoord_dist
      real(dblprec), dimension(NT,3) :: redcoord_tmp

      ! The idea is to be able to use the neighbour map routine that is used for the Jij´s
      ! Hence the first step would be to grab the bonding vectors for each atom and then
      ! select for the ones that are induced create a list with only the fixed ones and
      ! then select the nearest neighbours
      ! Allocation of working arrays
      allocate(redcoord_mod(NT,max_no_shells),stat=i_stat)
      call memocc(i_stat,product(shape(redcoord_mod))*kind(redcoord_mod),'redcoord_mod','setup_induced_nearest')
      redcoord_mod=0.0_dblprec
      ! First the modulus of the bonding vectors for the different shells is calculated
      do IT=1,NT
         do shells=1,max_no_shells
            redcoord_mod(IT,shells)=sqrt(redcoord(IT,shells,1)**2+redcoord(IT,shells,2)**2+&
               redcoord(IT,shells,3)**2)
         enddo
      enddo
      ! Now the modulus per atom type and per shell has been computed
      ! What must be done now is to calculate the minimum, distance per atom,
      ! this distance would determine the nearest neighbours

      do IT=1,NT
         temp_min=0.0_dblprec
         do shells=1, max_no_shells
            if ( (redcoord_mod(IT,shells).lt.temp_min).and.(redcoord_mod(IT,shells).gt.0.0_dblprec)&
               .and.shells.gt.1) then
               temp_min=redcoord_mod(IT,shells)
            else if (shells.eq.1) then
               temp_min=redcoord_mod(IT,shells)
               redcoord_tmp(IT,1:3)=redcoord(IT,shells,1:3)
            endif
         enddo
         redcoord_dist(IT)=temp_min
         ! Ater the minimum has been found one can store the nearest neighbour list
         ! It is important to notice that until now there has been no discrimination
         ! on whether the list includes induced moments or not
         ! Need to also find which is the maximum number of neighbours which are fixed
         temp_counter=0
         temp_max_counter=0
         do shells=1,max_no_shells
            diff=redcoord_mod(IT,shells)-redcoord_dist(IT)
            if (abs(diff).lt.ind_tol) then
               ind_redcoord(IT,shells,1:3)=redcoord(IT,shells,1:3)
               temp_counter=temp_counter+1
            endif
         enddo
         if (temp_counter.gt.temp_max_counter) temp_max_counter=temp_counter
      enddo

      ! deallocate any working arrays
      i_all=-product(shape(redcoord_mod))*kind(redcoord_mod)
      deallocate(redcoord_mod,stat=i_stat)
      call memocc(i_stat,i_all,'redcoord_mod','setup_induced_nearest')

   end subroutine  setup_induced_nearest

   !----------------------------------------------------------------------------
   !> @brief
   !> Setup auxiliarly neighbour list to treat induced moments
   !> it creates and the list for the induced atom neighbours
   !> @author
   !> Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine setup_induced_list(Natom,NT,NA,atype,max_no_neigh_ind,max_no_equiv,   &
         max_no_shells,do_sortcoup,Natom_full,atype_ch,Nchmax,do_ralloy,nm_ind,nn,     &
         nmdim_ind,ind_nlistsize,ind_nlist,ind_mom)

      use Sorting, only : MergeSortIR
      !
      implicit none
      ! .. Input variables
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments
      integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, dimension(NT), intent(in):: nn !< Number of neighbour shells
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
      integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced/fixed (1/0) for the atoms in the unit cell
      integer, dimension(max_no_shells,Natom), intent(in) :: nmdim_ind !< Dimension of neighbour map
      integer, dimension(Natom,max_no_shells,max_no_equiv), intent(in) :: nm_ind !< Neighbour map
      character :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
      ! .. Output variables
      integer, dimension(Natom), intent(out):: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(max_no_neigh_ind,Natom), intent(out) :: ind_nlist !< Neighbour list for induced moments
      ! .. Local variables
      integer :: i, j, k, l
      integer :: jchem, ichem
      integer :: tempn
      integer :: ncount,ncounter
      logical :: exis

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Routine to create a neighbour list which for the induced moments returns fixed moments
      ! and that for the fixed moments, returns induced ones
      ! based on the routine setup_neighbour_hamiltonian in hamiltonianinit.f90
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Loop over atoms
      !$omp parallel do default(shared) private(i,k,j,ncount,ncounter,exis,l,jchem,ichem)
      do i=1, Natom
         ncount = 1
         ncounter = 1
         ! Shells
         do k=1, Nn(atype(i))
            ! Sites in shell
            do j=1, nmdim_ind(k,i)
               ! Existing coupling
               if (nm_ind(i,k,j)>0) then
                  ! Check if the system is a random alloy to set the correct chemical type
                  if (do_ralloy==0) then
                     ichem=1
                     jchem=1
                     ! Check if the neighbours of the iduced atom iatom are fixed
                     if (ind_mom(atype(i),ichem).eq.1 .and.ind_mom(atype(nm_ind(i,k,j)),jchem).eq.0) then
                        ! If the current atom is induced and the neighbour is fixed
                        ! add to list
                        exis=.false.
                        do l=1,ncount-1
                           if(ind_nlist(l,i)==nm_ind(i,k,j)) exis=.true.
                        end do
                        if(.not.exis) then
                           ind_nlist(ncount,i) = nm_ind(i,k,j)
                           ncount = ncount + 1
                        end if
                        ! For the fixed moments save the induced moments
                     else if (ind_mom(atype(i),ichem).eq.0.and.ind_mom(atype(nm_ind(i,k,j)),jchem).eq.1) then

                        exis=.false.
                        do l=1,ncount-1
                           if(ind_nlist(l,i)==nm_ind(i,k,j)) exis=.true.
                        end do
                        if(.not.exis) then
                           ind_nlist(ncount,i) = nm_ind(i,k,j)
                           ncount = ncount + 1
                        end if
                     endif
                  else
                     ! If one considers a random alloy
                     ichem=atype_ch(i)
                     jchem=atype_ch(nm_ind(i,j,k))
                     ! Check if the neighbours of the iduced atom iatom are fixed
                     if (ind_mom(atype(i),ichem).eq.1 .and.ind_mom(atype(nm_ind(i,j,k)),jchem).eq.0) then

                        exis=.false.
                        do l=1,ncount-1
                           if(ind_nlist(l,i)==nm_ind(i,k,j)) exis=.true.
                        end do
                        if(.not.exis) then
                           ind_nlist(ncount,i) = nm_ind(i,k,j)
                           ncount = ncount + 1
                        end if
                        ! For the fixed moments only select the induced moments
                     else if(ind_mom(atype(i),ichem).eq.0.and.ind_mom(atype(nm_ind(i,j,k)),jchem).eq.1) then

                        exis=.false.
                        do l=1,ncount-1
                           if(ind_nlist(l,i)==nm_ind(i,k,j)) exis=.true.
                        end do
                        if(.not.exis) then
                           ind_nlist(ncount,i) = nm_ind(i,k,j)
                           ncount = ncount + 1
                        end if
                     endif
                  endif
               end if
            end do
         end do
         ind_nlistsize(i) = ncount-1
      end do
      !$omp end parallel do

      if(do_sortcoup == 'Y') then
         !$omp parallel do default(shared) private(i,j,k,tempn)
         do i=1,Natom
            ! sort neighbour list - bubble sort...
            do j=1,ind_nlistsize(i)
               do k=1,ind_nlistsize(i)-j
                  if ( ind_nlist(k,i) .gt. ind_nlist(k+1,i) ) then
                     tempn        = ind_nlist(k,i)
                     ind_nlist(k,i)   = ind_nlist(k+1,i)
                     ind_nlist(k+1,i) = tempn
                  end if
               end do
            end do
         end do
         !$omp end parallel do
      end if

   end subroutine setup_induced_list

   !---------------------------------------------------------------------------
   ! SUBROUTINE: calculate_init_ind_sus
   !> @brief
   !> Calculation of the weighting/renormalization factor for the moments from
   !> linear response, this array will be then used to update the magnitude of the
   !> induced moments.
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------
   subroutine calculate_init_ind_sus(Natom,Mensemble,   &
         max_no_neigh_ind,ind_list_full,           &
         ind_nlistsize,ind_nlist,mmom,emom,emomM,sus_ind)

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh_ind   !< Calculated maximum of neighbours for induced moments
      integer, dimension(Natom), intent(in) :: ind_list_full   !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
      integer, dimension(Natom), intent(in) :: ind_nlistsize   !< Size of neighbour list for induced moments
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      ! .. In/out variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! .. Output variables
      real(dblprec), dimension(Natom), intent(inout) :: sus_ind !< Scaling factor for the magneitc moment of the induced moments
      ! .. Local variables
      integer :: fix_neigh, fix, iatom
      real(dblprec) :: ave_norm, rescale_fac, mmom_ave
      real(dblprec), dimension(3) :: ave_mom, ave_emom

      ! Calculate the factor that will weight the exchange interactions
      !     if (renorm_coll=='Y') then
      !$omp parallel do default(shared), private(fix,ave_mom,ave_norm,ave_emom,mmom_ave)
      do iatom=1,Natom
         ave_mom=0.0_dblprec
         ave_norm=0.0_dblprec
         ave_emom=0.0_dblprec
         mmom_ave=0.0_dblprec
         rescale_fac=0.0_dblprec
         ! Check if the moments are induced or fixed
         if (ind_list_full(iatom)==1) then
            ! Sum over the magnitude of the fixed magnetic moments around the induced ones
            do fix_neigh=1,ind_nlistsize(iatom)
               fix=ind_nlist(fix_neigh,iatom)
               mmom_ave      = mmom_ave      + mmom(fix,1)
               ave_mom(1:3)  = ave_mom(1:3)  + emomM(1:3,fix,1)
               ave_emom(1:3) = ave_emom(1:3) + emom(1:3,fix,1)
            enddo
            ave_norm      = sqrt(ave_mom(1)**2+ave_mom(2)**2+ave_mom(3)**2)
            ave_emom(1:3) = ave_emom(1:3)/sqrt(ave_emom(1)**2+ave_emom(2)**2+ave_emom(3)**2)
            ! Calculate the weight factor
            sus_ind(iatom)     = mmom(iatom,1)/mmom_ave
            !sus_ind(iatom)     = ave_norm/mmom_ave
            !              print *,iatom,sus_ind(iatom)
         else
            ! If the current atom is not induced just set the weight factor to 1
            sus_ind(iatom)=1.0_dblprec
         endif
      enddo
      !$omp end parallel do
      !     else
      !        sus_ind(:)=1.0_dblprec
      !     endif
      !!!! Loop over all the atoms
      !!!!$omp parallel do default(shared), private(iatom,jatom,ineigh,jneigh,kk,ave_mom,ave_norm,rescale_fac,temp_Jij,temp_Jji,temp_Dij,temp_Dji)
      !!!do ii=1,fix_num
      !!!   iatom=fix_list(ii)
      !!!   ! Loop over all the neighbours
      !!!   do ineigh=1,nlistsize(iatom)
      !!!      ! Find the neighbouring atom
      !!!      jatom=nlist(ineigh,iatom)
      !!!      ! Find the position in the list of the neighbour atom
      !!!      do kk=1,nlistsize(jatom)
      !!!         ! Find the position in the neighbour list
      !!!         if(nlist(kk,jatom)==iatom) jneigh=kk
      !!!      enddo
      !!!      ! If the current atom is a fixed moment and the neighbour is an induced moment
      !!!      if (ind_list_full(jatom)==1) then
      !!!         rescale_fac=sus_ind(iatom)*sus_ind(jatom)
      !!!      else
      !!!         rescale_fac=1.0_dblprec
      !!!      endif
      !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!      ! Renormalize the exchange interaction
      !!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!      temp_Jij=ncoup(ineigh,iatom,1)*rescale_fac
      !!!      temp_Jji=ncoup(jneigh,jatom,1)*rescale_fac
      !!!      ncoup(ineigh,iatom,1)=(temp_Jij+temp_Jji)*0.5_dblprec
      !!!      ncoup(jneigh,jatom,1)=(temp_Jij+temp_Jji)*0.5_dblprec
      !!!   enddo
      !!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!   ! If the DMI is present
      !!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!   if (do_dm==1) then
      !!!      do ineigh=1,dmlistsize(iatom)
      !!!         jatom=dmlist(ineigh,iatom)
      !!!         ! Find the position in the list of the neighbour atom
      !!!         do kk=1,dmlistsize(jatom)
      !!!            if(dmlist(kk,jatom)==iatom) jneigh=kk
      !!!         enddo
      !!!         ! If the current atom is a fixed moment and the neighbour is an induced moment
      !!!         if (ind_list_full(jatom)==1) then
      !!!            rescale_fac=sus_ind(iatom)*sus_ind(jatom)
      !!!         else
      !!!            rescale_fac=1.0_dblprec
      !!!         endif
      !!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!         ! Renormalize the Dzyaloshinskii-Moriya vectors
      !!!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!         temp_Dij(1:3)=dm_vect(1:3,ineigh,iatom)*rescale_fac
      !!!         temp_Dji(1:3)=dm_vect(1:3,jneigh,jatom)*rescale_fac
      !!!         dm_vect(1:3,ineigh,iatom)=(temp_Dij(1:3)-temp_Dji(1:3))*0.5_dblprec
      !!!         dm_vect(1:3,jneigh,jatom)=(temp_Dji(1:3)-temp_Dij(1:3))*0.5_dblprec
      !!!      enddo
      !!!   endif
      !!!enddo
      !!!$omp end parallel do

   end subroutine calculate_init_ind_sus

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_induced_information
   !> @brief Creation of an auxilary list that handles whether an atom is induced
   !> or not incuding its chemical nature
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine setup_induced_information(NA,Natom,Nchmax,do_ralloy,Natom_full,       &
         anumb,achtype,ind_mom,ind_list_full,fix_list,fix_num,restartfile,rstep,mmom,  &
         emom,emomM,initmag,Mensemble,do_mom_legacy)

      use Restart, only: read_mag_conf

      implicit none

      ! .. Input variables
      integer, intent(in) :: NA           !< Number of atoms in one cell
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Nchmax       !< Max number of chemical components on each site in cell
      integer, intent(in) :: initmag      !< Mode of initialization of magnetic moments (1-4)
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      character(len=1), intent(in) :: do_mom_legacy
      integer, dimension(Natom), intent(in) :: anumb !< Type of atom
      integer, dimension(Natom_full), intent(in) :: achtype !< Actual site of atom for dilute system
      integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced/fixed (1/0) for the unit cell
      ! .. Output variables
      integer, intent(out) :: rstep !< Starting simulation step
      integer, dimension(Natom), intent(out) :: ind_list_full !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
      ! .. In/Out variables
      integer, intent(inout) :: fix_num   !< Number of "fixed" moments
      character(len=35), intent(inout) :: restartfile !< File containing restart information
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom    !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM !< Current magnetic moment vector
      integer, dimension(:), allocatable, intent(inout) ::fix_list   !< List containing the "fixed" moments

      integer :: ii,kk,i_stat,aux_sum
      fix_num=0

      if (initmag.eq.4) then
         call read_mag_conf(Natom,Mensemble,do_mom_legacy,rstep,restartfile,mmom,   &
            emom,emomM,ind_list_full,aux_sum)
         !$omp parallel do schedule(static) reduction(+:fix_num)
         do ii=1, Natom
            if (ind_list_full(ii)==0) then
               fix_num=fix_num+1
            endif
         enddo
         !$omp end parallel do
      else
         ! Fill up the full list of whether there an atom is induced or fixed
         if (do_ralloy==0) then
            !$omp parallel do schedule(static) reduction(+:fix_num)
            do ii=1, Natom
               if (ind_mom(anumb(ii),1)==0) then
                  ind_list_full(ii)=0
                  fix_num=fix_num+1
               else
                  ind_list_full(ii)=1
               endif
            enddo
            !$omp end parallel do
         else
            !$omp parallel do schedule(static) reduction(+:fix_num)
            do ii=1, Natom
               if (ind_mom(anumb(ii),achtype(ii))==0) then
                  ind_list_full(ii)=0
                  fix_num=fix_num+1
               else
                  ind_list_full(ii)=1
               endif
            enddo
            !$omp end parallel do
         endif

      endif

      allocate(fix_list(fix_num),stat=i_stat)
      call memocc(i_stat,product(shape(fix_list))*kind(fix_list),'fix_list','setup_induced_information')
      fix_list=0

      kk=0

      ! Create a list containing only the "fixed" moments, i.e. non-induced moments
      do ii=1,Natom
         if (ind_list_full(ii)==0) then
            kk=kk+1
            fix_list(kk)=ii
         endif
      enddo

   end subroutine setup_induced_information

   !!!    !----------------------------------------------------------------------------
   !!!    ! SUBROUTINE: renorm_ncoup_ind
   !!!    !> @brief Renormalization of the exchange interaction due to the change of
   !!!    !> direction of the fixed moments.
   !!!    !> @author Jonathan Chico
   !!!    !----------------------------------------------------------------------------
   !!!    subroutine renorm_ncoup_ind(do_dm,Natom,conf_num,Mensemble,max_no_neigh,         &
   !!!       max_no_dmneigh,max_no_neigh_ind,nlistsize,dmlistsize,ind_list_full,           &
   !!!       ind_nlistsize,nlist,dmlist,ind_nlist,sus_ind,mmom,emom,emomM,ncoup,dm_vect,   &
   !!!       fix_list,fix_num)
   !!! 
   !!!       implicit none
   !!! 
   !!!       integer, intent(in) :: do_dm     !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
   !!!       integer, intent(in) :: Natom     !< Number of atoms in system
   !!!       integer, intent(in) :: fix_num   !< Number of "fixed" moments
   !!!       integer, intent(in) :: conf_num  !< number of configurations for LSF
   !!!       integer, intent(in) :: Mensemble !< Number of ensembles
   !!!       integer, intent(in) :: max_no_neigh       !< Calculated maximum of neighbours for exchange
   !!!       integer, intent(in) :: max_no_dmneigh     !< Calculated number of neighbours with DM interactions
   !!!       integer, intent(in) :: max_no_neigh_ind   !< Calculated maximum of neighbours for induced moments
   !!!       integer, dimension(fix_num), intent(in) :: fix_list      !< List containing the "fixed" moments
   !!!       integer, dimension(Natom), intent(in) :: nlistsize       !< Size of neighbour list for Heisenberg exchange couplings
   !!!       integer, dimension(Natom), intent(in) :: dmlistsize      !< Size of neighbour list for DM
   !!!       integer, dimension(Natom), intent(in) :: ind_list_full   !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
   !!!       integer, dimension(Natom), intent(in) :: ind_nlistsize   !< Size of neighbour list for induced moments
   !!!       integer, dimension(max_no_neigh,Natom), intent(in) :: nlist    !< Neighbour list for Heisenberg exchange couplings
   !!!       integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist !< List of neighbours for DM
   !!!       integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
   !!! 
   !!!       ! .. In/out variables
   !!!       real(dblprec), dimension(Natom), intent(inout) :: sus_ind !< Scaling factor for the magneitc moment of the induced moments
   !!!       real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
   !!!       real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current unit moment vector
   !!!       real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
   !!!       real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(inout) :: ncoup !< Heisenberg exchange couplings
   !!!       real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(inout) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
   !!! 
   !!!       ! .. Local variables
   !!!       integer :: iatom,jatom,ineigh,jneigh,kk,ii
   !!!       integer :: fix, fix_neigh
   !!!       real(dblprec) :: ave_norm,temp_Jij,temp_Jji,rescale_fac,mmom_ave
   !!!       real(dblprec), dimension(3) :: ave_mom,ave_emom,temp_Dij,temp_Dji
   !!! 
   !!!       ! Loop over all the atoms
   !!!       !$omp parallel do default(shared), private(iatom,jatom,jneigh,ineigh,ave_mom,ave_norm,mmom_ave,ave_emom,fix,rescale_fac,temp_Jij,temp_Jji,temp_Dij,temp_Dji)
   !!!       do ii=1,fix_num
   !!!          iatom=fix_list(ii)
   !!!          ! Loop over all the neighbours
   !!!          do ineigh=1,nlistsize(iatom)
   !!!             ! Find the neighbouring atom
   !!!             jatom=nlist(ineigh,iatom)
   !!!             ! Find the position in the list of the neighbour atom
   !!!             do kk=1,nlistsize(jatom)
   !!!                ! Find the position in the neighbour list
   !!!                if(nlist(kk,jatom)==iatom) jneigh=kk
   !!!             enddo
   !!!             ! If the current atom is a fixed moment and the neighbour is an induced moment
   !!!             if (ind_list_full(jatom)==1) then
   !!!                ave_mom=0.0_dblprec
   !!!                ave_norm=0.0_dblprec
   !!!                ave_emom=0.0_dblprec
   !!!                mmom_ave=0.0_dblprec
   !!!                rescale_fac=0.0_dblprec
   !!!                ! Calculate the average magnetic moment
   !!!                do fix_neigh=1,ind_nlistsize(jatom)
   !!!                   fix=ind_nlist(fix_neigh,jatom)
   !!!                   mmom_ave      = mmom_ave      + mmom(fix,1)
   !!!                   ave_mom(1:3)  = ave_mom(1:3)  + emomM(1:3,fix,1)
   !!!                   ave_emom(1:3) = ave_emom(1:3) + emom(1:3,fix,1)
   !!!                enddo
   !!!                ave_norm    = sqrt(ave_mom(1)**2+ave_mom(2)**2+ave_mom(3)**2)
   !!!                rescale_fac = 1.0_dblprec/(sus_ind(iatom)*sus_ind(jatom))
   !!!                rescale_fac = rescale_fac*(ave_norm/mmom_ave)
   !!!                if(do_dm/=1) then
   !!!                   sus_ind(jatom) = ave_norm/mmom_ave
   !!!                endif
   !!!             else
   !!!                rescale_fac=1.0_dblprec
   !!!             endif
   !!!             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!             ! Renormalize the exchange interaction
   !!!             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!             temp_Jij=ncoup(ineigh,iatom,1)*rescale_fac
   !!!             temp_Jji=ncoup(jneigh,jatom,1)*rescale_fac
   !!!             ncoup(ineigh,iatom,1)=(temp_Jij+temp_Jji)*0.5_dblprec
   !!!             ncoup(jneigh,jatom,1)=(temp_Jij+temp_Jji)*0.5_dblprec
   !!! 
   !!!          enddo
   !!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!          ! If the DMI is present
   !!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!          if (do_dm==1) then
   !!!             do ineigh=1,dmlistsize(iatom)
   !!!                jatom=dmlist(ineigh,iatom)
   !!!                ! Find the position in the list of the neighbour atom
   !!!                do kk=1,dmlistsize(jatom)
   !!!                   if(dmlist(kk,jatom)==iatom) jneigh=kk
   !!!                enddo
   !!!                ! If the current atom is a fixed moment and the neighbour is an induced moment
   !!!                if (ind_list_full(jatom)==1) then
   !!!                   ave_mom=0.0_dblprec
   !!!                   ave_norm=0.0_dblprec
   !!!                   ave_emom=0.0_dblprec
   !!!                   mmom_ave=0.0_dblprec
   !!!                   rescale_fac=0.0_dblprec
   !!!                   ! Calculate the average magnetic moment
   !!!                   do fix_neigh=1,ind_nlistsize(jatom)
   !!!                      fix=ind_nlist(fix_neigh,jatom)
   !!!                      mmom_ave      = mmom_ave     + mmom(fix,1)
   !!!                      ave_mom(1:3)  = ave_mom(1:3) + emomM(1:3,fix,1)
   !!!                      ave_emom(1:3) = ave_mom(1:3) + emom(1:3,fix,1)
   !!!                   enddo
   !!! 
   !!!                   ave_norm           = sqrt(ave_mom(1)**2+ave_mom(2)**2+ave_mom(3)**2)
   !!!                   rescale_fac        = 1.0_dblprec/(sus_ind(iatom)*sus_ind(jatom))
   !!!                   rescale_fac        = rescale_fac*(ave_norm/mmom_ave)
   !!!                   sus_ind(jatom)     = (ave_norm/mmom_ave)
   !!!                else
   !!!                   rescale_fac=1.0_dblprec
   !!!                endif
   !!!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!                ! Renormalize the Dzyaloshinskii-Moriya vectors
   !!!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!                temp_Dij(1:3)=dm_vect(1:3,ineigh,iatom)*rescale_fac
   !!!                temp_Dji(1:3)=dm_vect(1:3,jneigh,jatom)*rescale_fac
   !!!                dm_vect(1:3,ineigh,iatom)=(temp_Dij(1:3)-temp_Dji(1:3))*0.5_dblprec
   !!!                dm_vect(1:3,jneigh,jatom)=(temp_Dji(1:3)-temp_Dij(1:3))*0.5_dblprec
   !!!             enddo
   !!!          endif
   !!!       enddo
   !!!       !$omp end parallel do
   !!! 
   !!!    end subroutine renorm_ncoup_ind

!!!    !-----------------------------------------------------------------------------
!!!    !  SUBROUTINE: induced_loadrestart
!!!    !> @brief
!!!    !> Read magnetic moments from file for the case of induced moments
!!!    !
!!!    !> @author
!!!    !> Jonathan Chico, based on the previously existent loadrestart routine
!!!    !-----------------------------------------------------------------------------
!!!    subroutine induced_loadrestart(Natom,Mensemble,restartfile,rstep,mmom,emom,emomM,&
!!!          ind_list_full)
!!!       !
!!!       !.. Implicit declarations
!!!       implicit none
!!! 
!!!       integer, intent(in) :: Natom !< Number of atoms in system
!!!       integer, intent(in) :: Mensemble !< Number of ensembles
!!!       integer, intent(out) :: rstep !< Starting simulation step
!!!       integer, dimension(Natom), intent(inout) :: ind_list_full !< Indication of whether a given moment is induced/fixed 1/0
!!!       real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
!!!       real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
!!!       real(dblprec), dimension(Natom,Mensemble), intent(out) :: mmom !< Magnitude of magnetic moments
!!!       character(len=35), intent(inout) :: restartfile !< File containing restart information
!!! 
!!!       integer :: i, j, k, l, ios
!!!       logical :: exists
!!! 
!!!       !.. Executable statements
!!!       inquire(file=restartfile,exist=exists)
!!!       if(exists) then
!!!          open(ifileno,iostat=ios, file=restartfile, status="old")
!!!          read (ifileno,*) rstep
!!!          do i=1,Mensemble
!!!             do j=1, Natom
!!!                read (ifileno,*) k, l, mmom(j,i), emom(1,j,i), emom(2,j,i), emom(3,j,i),ind_list_full(j)
!!!                emomM(:,j,i)=emom(:,j,i)*mmom(j,i)
!!!             end do
!!!          end do
!!!          close(ifileno)
!!!       else
!!!          write(*,*) 'ERROR: Restartfile ',trim(adjustl(restartfile)), ' does not exist.'
!!!          stop
!!!       end if
!!! 
!!!    end subroutine induced_loadrestart


   !---------------------------------------------------------------------------
   !> @brief
   !> Metropolis algorithm MonteCarlo update for a system with both fixed and induced
   !> magnetic moments.
   !> @details The approach is based in the work presented by Polesya et al.
   !> in Phys. Rev. B 82, 214409.
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------
   subroutine mc_update_ind_mom(Natom,Mensemble,iflip_a,temperature,temprescale,    &
         mode,max_no_neigh,nlistsize,nlist,ncoup,conf_num,      &
         mmom,emomM,emom,extfield,do_dip, Num_macro,&
         emomM_macro,emom_macro,mmom_macro,ind_nlistsize,ind_nlist,         &
         ind_list_full,sus_ind,do_lsf,lsf_metric,ind_mom_flag,max_no_neigh_ind)
      !
      use RandomNumbers, only: rng_uniform,rng_gaussian,rng_uniformP,rng_gaussianP
      use montecarlo_common
      use Constants,only: mub,k_bolt
      use InputData, only : ind_mom_type

      !.. Input variables
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, dimension(Natom),intent(in) :: iflip_a !< Flipping pattern
      real(dblprec), intent(in) :: temperature !< Temperature
      real(dblprec), intent(in) :: temprescale !< Temperature rescaling according to QHB
      character(len=1), intent(in) :: mode !< Simulation mode (M=MC, H=MC Heat Bath)
      ! LSF variables
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      integer,intent(in) :: lsf_metric !< LSF metric or phase space measure
      character(len=1), intent(in) :: do_lsf     !< Including LSF energy
      ! Heisenberg exchange variables
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      ! Moments variables
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      ! External magnetic fields
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      ! Dipolar interaction variables
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      ! Induced moment variables
      integer, intent(in) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments
      integer, dimension(Natom), intent(in) :: ind_list_full !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
      integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      character(len=1), intent(in) :: ind_mom_flag
      ! .. Macrocell variables
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(Num_macro,Mensemble), intent(inout)   :: mmom_macro !< Magnitude of the macrocell magnetic moments
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro !< Unit vector of the macrocell magnetic moment

      !.. Local scalars
      !
      integer :: iatom, k, icell, i
      real(dblprec) :: de !< New trial magnitude of moment
      real(dblprec) :: macro_mag_trial,delta
      !
      !.. Local arrays
      !
      real(dblprec), dimension(3) :: newmom,macro_trial
      real(dblprec), dimension(natom) :: flipprob_a 
      real(dblprec) :: newmmom_a
      !     real(dblprec), dimension(3,natom,mensemble) :: newmom_a

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! The MC step for this approach consists on two operations 1) find a fixed moment
      ! rotate that moment 2) find the nearest neighbour induced moments and move them
      ! as to follow the nearest neighbour fixed moments and change their magnitude
      ! after that calculate the total change in the energy
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      delta=(2.0/25.0)*(k_bolt*temperature/mub)**(0.20_dblprec)
      call rng_uniformP(rng_work_u(:,:,:),3*natom*mensemble)
      call rng_gaussianP(rng_work_g(:,:,:),3*natom*mensemble,1.0_dblprec)
      !              if(use_vsl) then
#ifdef VSL
      !     !$omp parallel do default(shared),private(i,k,newmom),schedule(auto),collapse(2)
#endif
      do i=1,Natom
         do k=1,mensemble
            call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,delta,rng_work_u(:,i,k),rng_work_g(:,i,k))
            mom_work_i(1:3,i,k)=newmom(1:3)
         enddo
      enddo
#ifdef VSL
      !     !$omp end parallel do
#endif
      !               else
      !                  do i=1,Natom
      !                     do k=1,mensemble
      !                        call choose_random_flip(emom,newmom,Natom,Mensemble,i,k,delta,flipprob_m(:,i,k),flipprob_g(:,i,k))
      !                        mom_work_i(1:3,i,k)=newmom(1:3)
      !                     enddo
      !                  enddo
      !               end if
      !!! #ifdef VSL
      !!!          !$omp parallel do default(shared),private(iatom,newmom),schedule(static,1) collapse(2)
      !!! #endif
      !!!       do k=1,Mensemble
      !!!          do iatom=1,Natom
      !!! !           ! Only flip atoms which have a fixed magnetic moment
      !!! !           if (ind_list_full(iatom).eq.0) then
      !!!             !if (ind_list_full(iflip_a(iatom)).eq.0) then
      !!!                call rng_uniform(fran,3)
      !!!                call rng_gaussian(gran,3,1.0_dblprec)
      !!!                call choose_random_flip(emom,newmom,Natom,Mensemble,iatom,k,delta,fran,gran)
      !!!                mom_work_i(1:3,iatom)=newmom(1:3)
      !!!                !           else
      !!!                !              mom_work_i(1:3,iatom)=emom(1:3,iatom,k)
      !!!                !           endif
      !!!                !           emomM_trial(:,iatom)=emomM(:,iatom,k)
      !!!             enddo
      !!!          enddo
      !!! #ifdef VSL
      !!!          !$omp end parallel do
      !!! #endif

      call rng_uniformP(flipprob_a,natom)
      ! Calculate energy and flip spin if preferred including the contribution from induced moments
      !$omp parallel do default(shared) private(k,iatom,de,newmom) collapse(2) schedule(static,100)
      do k=1,Mensemble
         do iatom=1, Natom
            if (ind_list_full(iflip_a(iatom)).eq.0.or.ind_mom_type==1) then
               !if (ind_list_full(iflip_a(iatom)).eq.0) then
               ! Only directly update the fixed moments, induced moments are updated in an indirect fashion
               !if (ind_list_full(iatom).eq.0) then
               !print *,'-> ',iatom,iflip_a(iatom),ind_list_full(iflip_a(iatom))
               !print '(a,3f12.6)','-> ',mom_work_i(1:3,iflip_a(iatom))
               call calculate_energy_wIND_v3(k,Natom,Mensemble,max_no_neigh,nlistsize,nlist, &
                  ncoup,conf_num,iflip_a(iatom),mom_work_i(1:3,iflip_a(iatom),k),&
                  mmom,emomM,extfield,        &
                  ind_nlistsize,ind_nlist,ind_list_full,sus_ind,max_no_neigh_ind, de)
               !        if (ind_list_full(iflip_a(iatom)).eq.0) then
               !           ! Calculate the energy containing the induced moments
               !           !call calculate_energy_wIND_v2(k,Natom,Mensemble,max_no_neigh,nlistsize, &
               !           call calculate_energy_wIND(k,Natom,Mensemble,max_no_neigh,nlistsize, &
               !              nlist,ncoup,ncoupD,conf_num,exc_inter,do_dm,max_no_dmneigh,       &
               !              dmlistsize,dmlist,dm_vect,do_pd,nn_pd_tot,pdlistsize,pdlist,      &
               !              pd_vect,do_biqdm,nn_biqdm_tot,biqdmlistsize,biqdmlist,biqdm_vect, &
               !              do_bq,nn_bq_tot,bqlistsize,bqlist,j_bq,taniso,taniso_diff,eaniso, &
               !              eaniso_diff,kaniso,kaniso_diff,sb,sb_diff,mult_axis,              &
               !              iflip_a(iatom),mom_work_i(1:3,iflip_a(iatom)),mmom,emomM,emom,      &
               !              extfield,do_dip,Qdip,ind_nlistsize,ind_nlist,&
               !              ind_list_full,       &
               !              sus_ind,max_no_neigh_ind,Num_macro,max_num_atom_macro_cell,       &
               !              cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,          &
               !              Qdip_macro,icell,macro_mag_trial,macro_trial,de,do_anisotropy)
               !        else 
               !             ! Metropolis algorithm, either in Ising or Loop Algorithm form
               !             call calculate_energy(Natom, Mensemble, Natom, conf_num, do_dm , do_pd, do_biqdm, do_bq, 0,&
               !                emomM, emom, mmom, iflip_a(iatom), mom_work_i(1:3,iflip_a(iatom)), extfield, de, k, &
               !                mult_axis, do_dip,Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,&
               !                macro_atom_nlist,emomM_macro,icell,macro_mag_trial,macro_trial,exc_inter,do_anisotropy)
               !        endif

               if(mode=='D') then
                  call flip_g(Natom,Mensemble,emom,emomM,mmom,iflip_a(iatom),       &
                     mom_work_i(1:3,iflip_a(iatom),k),newmmom_a,de,     &
                     temperature,temprescale,do_lsf,k,flipprob_a(iatom),lsf_metric, &
                     ind_nlistsize,ind_nlist,ind_mom_flag,max_no_neigh_ind,sus_ind, &
                     do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro, &
                     emom_macro,emomM_macro)
               else
                  !print *,'flip_a',Natom,iatom,iflip_a(iatom)
                  call flip_a(Natom,Mensemble,emom,emomM,mmom,iflip_a(iatom),       &
                     mom_work_i(1:3,iflip_a(iatom),k),newmmom_a,de,     &
                     temperature,temprescale,do_lsf,k,flipprob_a(iatom),lsf_metric, &
                     ind_nlistsize,ind_nlist,ind_mom_flag,max_no_neigh_ind,sus_ind, &
                     ind_list_full,do_dip,Num_macro,icell,macro_mag_trial,          &
                     macro_trial,mmom_macro,emom_macro,emomM_macro)
               endif
            end if
            !endif
         enddo
      enddo
      !$omp end parallel do

   end subroutine mc_update_ind_mom

   !!!    !---------------------------------------------------------------------------
   !!!    !> @brief
   !!!    !> Calculate the total energy of a single fixed magnetic moment plus including the
   !!!    !> indirect energy contributions of the induced magnetic moments
   !!!    !>
   !!!    !> @author
   !!!    !> Jonathan Chico
   !!!    !---------------------------------------------------------------------------
   !!!    subroutine calculate_energy_wIND(k,Natom,Mensemble,max_no_neigh,nlistsize,nlist, &
   !!!       ncoup,ncoupD,conf_num,exc_inter,do_dm,max_no_dmneigh,dmlistsize,dmlist,       &
   !!!       dm_vect,do_pd,nn_pd_tot,pdlistsize,pdlist,pd_vect,do_biqdm,nn_biqdm_tot,      &
   !!!       biqdmlistsize,biqdmlist,biqdm_vect,do_bq,nn_bq_tot,bqlistsize,bqlist,j_bq,    &
   !!!       taniso,taniso_diff,eaniso,eaniso_diff,kaniso,kaniso_diff,sb,sb_diff,mult_axis,&
   !!!       iflip,newmom,mmom,emomM,emom,extfield,do_dip,Qdip,ind_nlistsize,ind_nlist,    &
   !!!       ind_list_full,sus_ind,max_no_neigh_ind,Num_macro,max_num_atom_macro_cell,     &
   !!!       cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,Qdip_macro,icell,     &
   !!!       macro_mag_trial,macro_trial,de,do_anisotropy)
   !!! 
   !!!       use Constants, only : mub
   !!!       use macrocells, only : calc_trial_macro_mom
   !!! 
   !!!       !.. Implicit declarations
   !!!       implicit none
   !!! 
   !!!       ! System variables
   !!!       integer, intent(in) :: k !< Current ensemble
   !!!       integer, intent(in) :: Natom !< Number of atoms in system
   !!!       integer, intent(in) :: Mensemble !< Number of ensembles
   !!!       ! LSF variables
   !!!       integer, intent(in) :: conf_num   !< Number of configurations for LSF
   !!!       character(len=1), intent(in) :: exc_inter !> Interpolation of Jij between FM/DLM (Y/N)
   !!!       ! Heisenberg exchange variables
   !!!       integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
   !!!       integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
   !!!       integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
   !!!       real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
   !!!       real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
   !!!       ! DMI  variables
   !!!       integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
   !!!       integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
   !!!       integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
   !!!       integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
   !!!       real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
   !!!       ! PD interactions variables
   !!!       integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
   !!!       integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
   !!!       integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
   !!!       integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
   !!!       real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
   !!!       ! BIQDM variables
   !!!       integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
   !!!       integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
   !!!       integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
   !!!       integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
   !!!       real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange coupling
   !!!       ! BQ variables
   !!!       integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
   !!!       integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
   !!!       integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
   !!!       integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
   !!!       real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
   !!!       ! Anisotropy variables
   !!!       integer, intent(in) :: do_anisotropy
   !!!       integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
   !!!       integer, dimension(Natom),intent(in) :: taniso_diff !< Type of anisotropy (0-2)
   !!!       real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
   !!!       real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
   !!!       real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
   !!!       real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
   !!!       real(dblprec), dimension(Natom), intent(in) :: sb!< Ratio between the Anisotropy constants
   !!!       real(dblprec), dimension(Natom), intent(in) :: sb_diff!< Ratio between the Anisotropy constants
   !!!       character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
   !!!       ! Moments variables
   !!!       integer, intent(in) :: iflip !< Atom to flip spin for
   !!!       real(dblprec), dimension(3), intent(in) :: newmom !< New trial moment
   !!!       real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
   !!!       real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
   !!!       real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
   !!!       ! External magnetic fields
   !!!       real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
   !!!       ! Dipolar interaction variables
   !!!       integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
   !!!       real(dblprec), dimension(3,3,Natom,Natom), intent(in) :: Qdip !< Matrix for dipole-dipole interaction
   !!!       ! Induced moment variables
   !!!       integer, intent(in) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments
   !!!       integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for induced moments
   !!!       integer, dimension(Natom), intent(in) :: ind_list_full !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
   !!!       integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
   !!!       real(dblprec), dimension(Natom), intent(in) :: sus_ind
   !!!       ! Macrocell variables
   !!!       integer, intent(in) :: Num_macro !< Number of macrocells in the system
   !!!       integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
   !!!       integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
   !!!       integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
   !!!       integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
   !!!       real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
   !!!       real(dblprec), dimension(3,3,Num_macro,Num_macro), intent(in) :: Qdip_macro !< Matrix for macro spin dipole-dipole
   !!!       ! Output variables
   !!!       real(dblprec), intent(out):: de  !< Energy difference
   !!!       integer, intent(in) :: icell
   !!!       real(dblprec), intent(in) :: macro_mag_trial
   !!!       real(dblprec), dimension(3), intent(in) :: macro_trial
   !!! 
   !!!       !.. Local scalars
   !!!       integer :: j,fix,curr_fix,neigh_test
   !!!       real(dblprec) :: tt, e_c, e_t
   !!!       real(dblprec) :: excscale
   !!! 
   !!!       !.. Local arrays
   !!!       real(dblprec), dimension(3) :: beff_t, trialmom,ave_mom
   !!! 
   !!!       neigh_test=0
   !!!       !.. Executable statements
   !!! 
   !!!       ! First calculate effective field
   !!!       beff_t(1) = 0_dblprec
   !!!       beff_t(2) = 0_dblprec
   !!!       beff_t(3) = 0_dblprec
   !!!       tt=0.0_dblprec
   !!! 
   !!!       e_c=0.0_dblprec
   !!!       e_t=0.0_dblprec
   !!!       trialmom(:)=newmom(:)*mmom(iflip,k)
   !!! 
   !!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!       ! The iflip magnetic moment should be only fixed magnetic moments, if the
   !!!       ! neigbour is fixed everything procedes as usual, if the neighbour is induced
   !!!       ! then one must rotate the induced moment according to the direction of its neighbouring
   !!!       ! fixed moments, and then the energy of that new configuration must be calculated
   !!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!       if (exc_inter=='N') then
   !!!          ! Exchange
   !!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!          ! Calculation of the exchange term
   !!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!          do j=1,nlistsize(iflip)
   !!!             neigh_test=nlist(j,iflip)
   !!!             ! Save the name of the neighbour and its chemical type
   !!!             ! If the neighbour is fixed everything goes as normal
   !!!             if (ind_list_full(neigh_test).eq.0) then
   !!!                e_c=e_c-ncoup(j,iflip,1)*sum(emomM(:,iflip,k)*emomM(:,nlist(j,iflip),k))
   !!!                e_t=e_t-ncoup(j,iflip,1)*sum(trialmom(:)*emomM(:,nlist(j,iflip),k))
   !!!                ! If the neighbour is induced rotation of the neighbour must be taken into account
   !!!             else if (ind_list_full(neigh_test).eq.1) then
   !!!                ave_mom=0.0_dblprec
   !!!                ! For each induced neighbour one must sum the corresponding fixed neighbours
   !!!                ! The difficulty now lies on selecting only nearest neighbours that are induced
   !!!                ! This command searches if the neighbour neigh_test is present in the list of induced neighbours
   !!!                if (ANY(ind_nlist(:,iflip)==neigh_test) ) then
   !!!                   ! If the neighbour is induced and it is a nearest neighbour then find its fixed neighbours
   !!!                   ! Sum over the fixed neighbours of the current induced neighbour
   !!!                   do fix=1,ind_nlistsize(neigh_test)
   !!!                      curr_fix=ind_nlist(fix,neigh_test)
   !!!                      ave_mom(1:3)=ave_mom(1:3)+emomM(1:3,curr_fix,k)
   !!!                   enddo
   !!!                   ! This calculates the energy of the previous configuration
   !!!                   e_c=e_c-ncoup(j,iflip,1)*sus_ind(neigh_test)*sum(emomM(:,iflip,k)*ave_mom(:))
   !!!                   ! This calculates the change in energy as the induced moment addapts to the changing moments
   !!!                   e_t=e_t-ncoup(j,iflip,1)*sus_ind(neigh_test)*sum(trialmom(:)*ave_mom(:))
   !!!                endif
   !!!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!                ! End calculation of the influence of the induced neighbours
   !!!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!             endif
   !!!          end do
   !!!       else
   !!!          do j=1,nlistsize(iflip)
   !!!             beff_t=beff_t+emomM(:,nlist(j,iflip),k)
   !!!             tt=tt+mmom(nlist(j,iflip),k)
   !!!          enddo
   !!!          excscale=sqrt(beff_t(1)**2+beff_t(2)**2+beff_t(3)**2)/tt
   !!!          do j=1,nlistsize(iflip)
   !!!             e_c=e_c-(excscale*ncoup(j,iflip,1)+ &
   !!!                (1._dblprec-excscale)*ncoupD(j,iflip,1))*sum(emomM(:,iflip,k)*emomM(:,nlist(j,iflip),k))
   !!!             e_t=e_t-(excscale*ncoup(j,iflip,1)+ &
   !!!                (1._dblprec-excscale)*ncoupD(j,iflip,1))*sum(trialmom(:)*emomM(:,nlist(j,iflip),k))
   !!!          end do
   !!!       endif
   !!! 
   !!! 
   !!! !!!       ! Anisotropy
   !!! !!!       if (do_anisotropy==1) then
   !!! !!!          ! Uniaxial anisotropy
   !!! !!!          if (taniso(iflip)==1) then
   !!! !!!             tta=sum(emomM(:,iflip,k)*eaniso(:,iflip))
   !!! !!!             ttb=sum(trialmom(:)*eaniso(:,iflip))
   !!! !!!             e_c=e_c+kaniso(1,iflip)*(tta**2)+kaniso(2,iflip)*(tta**4)
   !!! !!!             e_t=e_t+kaniso(1,iflip)*(ttb**2)+kaniso(2,iflip)*(ttb**4)
   !!! !!!             ! Cubic anisotropy
   !!! !!!          elseif (taniso(iflip)==2) then
   !!! !!!             e_c=e_c-kaniso(1,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
   !!! !!!                emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)- &
   !!! !!!                kaniso(2,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
   !!! !!!             e_t=e_t-kaniso(1,iflip)*(trialmom(1)**2*trialmom(2)**2+ &
   !!! !!!                trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)- &
   !!! !!!                kaniso(2,iflip)*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
   !!! !!!          endif
   !!! !!!          ! When both Cubic and Uniaxial are switched on
   !!! !!!          if (taniso(iflip)==7) then
   !!! !!!             ! Uniaxial anisotropy
   !!! !!!             tta=(emomM(1,iflip,k)*eaniso(1,iflip)+emomM(2,iflip,k)*eaniso(2,iflip)+emomM(3,iflip,k)*eaniso(3,iflip))
   !!! !!!             ttb=(trialmom(1)*eaniso(1,iflip)+trialmom(2)*eaniso(2,iflip)+trialmom(3)*eaniso(3,iflip))
   !!! !!!             e_c=e_c+kaniso(1,iflip)*(tta**2)+kaniso(2,iflip)*(tta**4)
   !!! !!!             e_t=e_t+kaniso(1,iflip)*(ttb**2)+kaniso(2,iflip)*(ttb**4)
   !!! !!!             ! Cubic anisotropy
   !!! !!!             aw1=kaniso(1,iflip)*sb(iflip)
   !!! !!!             aw2=kaniso(2,iflip)*sb(iflip)
   !!! !!!             e_c=e_c+aw1*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
   !!! !!!                emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+ &
   !!! !!!                aw2*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
   !!! !!!             e_t=e_t+aw1*(trialmom(1)**2*trialmom(2)**2+ &
   !!! !!!                trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
   !!! !!!                aw2*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
   !!! !!!          endif
   !!! !!! 
   !!! !!!          ! This includes a secondary anisotropy energy
   !!! !!!          if (mult_axis=='Y') then
   !!! !!!             ! Uniaxial anisotropy
   !!! !!!             if (taniso_diff(iflip)==1) then
   !!! !!!                tta=(emomM(1,iflip,k)*eaniso_diff(1,iflip)+emomM(2,iflip,k)*eaniso_diff(2,iflip)+emomM(3,iflip,k)*eaniso_diff(3,iflip))
   !!! !!!                ttb=(trialmom(1)*eaniso_diff(1,iflip)+trialmom(2)*eaniso_diff(2,iflip)+trialmom(3)*eaniso_diff(3,iflip))
   !!! !!!                e_c=e_c+kaniso_diff(1,iflip)*(tta**2)+kaniso_diff(2,iflip)*(tta**4)
   !!! !!!                e_t=e_t+kaniso_diff(1,iflip)*(ttb**2)+kaniso_diff(2,iflip)*(ttb**4)
   !!! !!! 
   !!! !!!             ! Cubic anisotropy
   !!! !!!             elseif (taniso_diff(iflip)==2) then
   !!! !!!                e_c=e_c+kaniso_diff(1,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
   !!! !!!                   emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+ &
   !!! !!!                   kaniso_diff(2,iflip)*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
   !!! !!!                e_t=e_t+kaniso_diff(1,iflip)*(trialmom(1)**2*trialmom(2)**2+ &
   !!! !!!                   trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
   !!! !!!                   kaniso_diff(2,iflip)*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
   !!! !!! 
   !!! !!!             endif
   !!! !!!             ! When both Cubic and Uniaxial are switched on
   !!! !!!             if (taniso_diff(iflip)==7) then
   !!! !!!                ! Uniaxial anisotropy
   !!! !!!                tta=(emomM(1,iflip,k)*eaniso_diff(1,iflip)+emomM(2,iflip,k)*eaniso_diff(2,iflip)+emomM(3,iflip,k)*eaniso_diff(3,iflip))
   !!! !!!                ttb=(trialmom(1)*eaniso_diff(1,iflip)+trialmom(2)*eaniso_diff(2,iflip)+trialmom(3)*eaniso_diff(3,iflip))
   !!! !!!                e_c=e_c+kaniso_diff(1,iflip)*(tta**2)+kaniso_diff(2,iflip)*(tta**4)
   !!! !!!                e_t=e_t+kaniso_diff(1,iflip)*(ttb**2)+kaniso_diff(2,iflip)*(ttb**4)
   !!! !!!                ! Cubic anisotropy
   !!! !!!                aw1=kaniso_diff(1,iflip)*sb_diff(iflip)
   !!! !!!                aw2=kaniso_diff(2,iflip)*sb_diff(iflip)
   !!! !!!                e_c=e_c+aw1*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2+ &
   !!! !!!                   emomM(2,iflip,k)**2*emomM(3,iflip,k)**2+emomM(3,iflip,k)**2*emomM(1,iflip,k)**2)+ &
   !!! !!!                   aw2*(emomM(1,iflip,k)**2*emomM(2,iflip,k)**2*emomM(3,iflip,k)**2)
   !!! !!!                e_t=e_t+aw1*(trialmom(1)**2*trialmom(2)**2+ &
   !!! !!!                   trialmom(2)**2*trialmom(3)**2+trialmom(3)**2*trialmom(1)**2)+ &
   !!! !!!                   aw2*(trialmom(1)**2*trialmom(2)**2*trialmom(3)**2)
   !!! !!!             endif
   !!! !!!          endif
   !!! !!!       endif
   !!! !!! 
   !!! !!!       ! DM interaction
   !!! !!!       if (do_dm==1) then
   !!! !!!          do j=1,dmlistsize(iflip)
   !!! !!!             e_c=e_c+dm_vect(1,j,iflip)*(emomM(2,iflip,k)*emomM(3,dmlist(j,iflip),k)- &
   !!! !!!                emom(3,iflip,k)*emomM(2,dmlist(j,iflip),k))+ &
   !!! !!!                dm_vect(2,j,iflip)*(emomM(3,iflip,k)*emomM(1,dmlist(j,iflip),k)- &
   !!! !!!                emomM(1,iflip,k)*emomM(3,dmlist(j,iflip),k))+ &
   !!! !!!                dm_vect(3,j,iflip)*(emom(1,iflip,k)*emomM(2,dmlist(j,iflip),k)- &
   !!! !!!                emomM(2,iflip,k)*emomM(1,dmlist(j,iflip),k))
   !!! !!!             e_t=e_t+dm_vect(1,j,iflip)*(trialmom(2)*emomM(3,dmlist(j,iflip),k)- &
   !!! !!!                trialmom(3)*emomM(2,dmlist(j,iflip),k))+ &
   !!! !!!                dm_vect(2,j,iflip)*(trialmom(3)*emomM(1,dmlist(j,iflip),k)- &
   !!! !!!                trialmom(1)*emomM(3,dmlist(j,iflip),k))+ &
   !!! !!!                dm_vect(3,j,iflip)*(trialmom(1)*emomM(2,dmlist(j,iflip),k)- &
   !!! !!!                trialmom(2)*emomM(1,dmlist(j,iflip),k))
   !!! !!!          end do
   !!! !!!       end if
   !!! !!! 
   !!! !!!       ! PD interaction
   !!! !!!       if(do_pd==1) then
   !!! !!!          do j=1,pdlistsize(iflip)
   !!! !!!             e_c=e_c-pd_vect(1,j,iflip)*emomM(1,iflip,k)*emomM(1,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(4,j,iflip)*emomM(1,iflip,k)*emomM(2,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(5,j,iflip)*emomM(1,iflip,k)*emomM(3,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(4,j,iflip)*emomM(2,iflip,k)*emomM(1,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(2,j,iflip)*emomM(2,iflip,k)*emomM(2,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(6,j,iflip)*emomM(2,iflip,k)*emomM(3,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(5,j,iflip)*emomM(3,iflip,k)*emomM(1,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(6,j,iflip)*emomM(3,iflip,k)*emomM(2,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(3,j,iflip)*emomM(3,iflip,k)*emomM(3,pdlist(j,iflip),k)
   !!! !!!             e_t=e_t-pd_vect(1,j,iflip)*trialmom(1)*emomM(1,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(4,j,iflip)*trialmom(1)*emomM(2,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(5,j,iflip)*trialmom(1)*emomM(3,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(4,j,iflip)*trialmom(2)*emomM(1,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(2,j,iflip)*trialmom(2)*emomM(2,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(6,j,iflip)*trialmom(2)*emomM(3,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(5,j,iflip)*trialmom(3)*emomM(1,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(6,j,iflip)*trialmom(3)*emomM(2,pdlist(j,iflip),k)- &
   !!! !!!                pd_vect(3,j,iflip)*trialmom(3)*emomM(3,pdlist(j,iflip),k)
   !!! !!!          end do
   !!! !!!       end if
   !!! !!! 
   !!! !!!       ! BIQDM interaction
   !!! !!!       if(do_biqdm==1) then
   !!! !!!          do j=1,biqdmlistsize(iflip)
   !!! !!!             e_c=e_c-biqdm_vect(1,j,iflip)*(emomM(2,iflip,k)*emomM(3,biqdmlist(j,iflip),k)- &
   !!! !!!                emomM(3,iflip,k)*emomM(2,biqdmlist(j,iflip),k))**2- &
   !!! !!!                biqdm_vect(1,j,iflip)*(emom(3,iflip,k)*emomM(1,biqdmlist(j,iflip),k)- &
   !!! !!!                emomM(1,iflip,k)*emomM(3,biqdmlist(j,iflip),k))**2- &
   !!! !!!                biqdm_vect(1,j,iflip)*(emom(1,iflip,k)*emomM(2,biqdmlist(j,iflip),k)- &
   !!! !!!                emomM(2,iflip,k)*emomM(1,biqdmlist(j,iflip),k))**2
   !!! !!!             e_t=e_t-biqdm_vect(1,j,iflip)*(trialmom(2)*emomM(3,biqdmlist(j,iflip),k)- &
   !!! !!!                trialmom(3)*emomM(2,biqdmlist(j,iflip),k))**2- &
   !!! !!!                biqdm_vect(1,j,iflip)*(trialmom(3)*emomM(1,biqdmlist(j,iflip),k)- &
   !!! !!!                trialmom(1)*emomM(3,biqdmlist(j,iflip),k))**2- &
   !!! !!!                biqdm_vect(1,j,iflip)*(trialmom(1)*emomM(2,biqdmlist(j,iflip),k)- &
   !!! !!!                trialmom(2)*emomM(1,biqdmlist(j,iflip),k))**2
   !!! !!!          end do
   !!! !!!       end if
   !!! !!! 
   !!! !!!       ! BQ interaction
   !!! !!!       if(do_bq==1) then
   !!! !!!          do j=1,bqlistsize(iflip)
   !!! !!!             ! current spin
   !!! !!!             bqmdot=emomM(1,bqlist(j,iflip),k)*emomM(1,iflip,k)+&
   !!! !!!                emomM(2,bqlist(j,iflip),k)*emomM(2,iflip,k)+&
   !!! !!!                emomM(3,bqlist(j,iflip),k)*emomM(3,iflip,k)
   !!! !!!             e_c=e_c-j_bq(j,iflip)*bqmdot**2
   !!! !!!             !trial spin
   !!! !!!             bqmdot=emomM(1,bqlist(j,iflip),k)*trialmom(1) + &
   !!! !!!                emomM(2,bqlist(j,iflip),k)*trialmom(2) + &
   !!! !!!                emomM(3,bqlist(j,iflip),k)*trialmom(3)
   !!! !!!             e_t=e_t-j_bq(j,iflip)*bqmdot**2
   !!! !!!          end do
   !!! !!!       end if
   !!! !!! 
   !!! !!!       ! Dipole-dipole interaction
   !!! !!!       if (do_dip==1) then
   !!! !!!          do j=1,Natom
   !!! !!!             do mu=1,3
   !!! !!!                do nu=1,3
   !!! !!!                   e_c=e_c-emomM(mu,iflip,k)*Qdip(nu,mu,j,iflip)*emomM(nu,j,k)
   !!! !!!                   e_t=e_t-trialmom(mu)*Qdip(nu,mu,j,iflip)*emomM(nu,j,k)
   !!! !!!                enddo
   !!! !!!             enddo
   !!! !!!          end do
   !!! !!!       elseif(do_dip==2) then
   !!! !!!          ! Calculation of the trial moment in the macrocell approach
   !!! !!!          call calc_trial_macro_mom(k,iflip,Natom,Mensemble,Num_macro,max_num_atom_macro_cell,&
   !!! !!!             macro_nlistsize,macro_atom_nlist,trialmom,emomM,emomM_macro,macro_mag_trial,macro_trial)
   !!! !!!          icell=cell_index(iflip)
   !!! !!!          do j=1, Num_macro
   !!! !!!             do mu=1,3
   !!! !!!                do nu=1,3
   !!! !!!                   e_c=e_c-emomM_macro(mu,icell,k)*Qdip_macro(nu,mu,j,icell)*emomM_macro(nu,j,k)/macro_nlistsize(icell)
   !!! !!!                   e_t=e_t-macro_trial(mu)*Qdip_macro(nu,mu,j,icell)*emomM_macro(nu,j,k)/macro_nlistsize(icell)
   !!! !!!                enddo
   !!! !!!             enddo
   !!! !!!          enddo
   !!! !!!       endif
   !!! 
   !!!       ! Add Zeeman term
   !!!       e_c=e_c-extfield(1)*emomM(1,iflip,k)-extfield(2)*emomM(2,iflip,k)-extfield(3)*emomM(3,iflip,k)
   !!!       e_t=e_t-extfield(1)*trialmom(1)-extfield(2)*trialmom(2)-extfield(3)*trialmom(3)
   !!! 
   !!!       !Energy difference
   !!!       tt=e_t-e_c
   !!!       de=mub*tt
   !!!       return
   !!!    end subroutine calculate_energy_wIND

   !!!    !---------------------------------------------------------------------------
   !!!    !> @brief
   !!!    !> Calculate the total energy of a single fixed magnetic moment plus including the
   !!!    !> indirect energy contributions of the induced magnetic moments
   !!!    !>
   !!!    !> @author
   !!!    !> Anders Bergman, Jonathan Chico
   !!!    !---------------------------------------------------------------------------
   !!!    subroutine calculate_energy_wIND_v2(k,Natom,Mensemble,max_no_neigh,nlistsize,nlist, &
   !!!       ncoup,ncoupD,conf_num,exc_inter,do_dm,max_no_dmneigh,dmlistsize,dmlist,       &
   !!!       dm_vect,do_pd,nn_pd_tot,pdlistsize,pdlist,pd_vect,do_biqdm,nn_biqdm_tot,      &
   !!!       biqdmlistsize,biqdmlist,biqdm_vect,do_bq,nn_bq_tot,bqlistsize,bqlist,j_bq,    &
   !!!       taniso,taniso_diff,eaniso,eaniso_diff,kaniso,kaniso_diff,sb,sb_diff,mult_axis,&
   !!!       iflip,newmom,mmom,emomM,emom,extfield,do_dip,Qdip,ind_nlistsize,ind_nlist,    &
   !!!       ind_list_full,sus_ind,max_no_neigh_ind,Num_macro,max_num_atom_macro_cell,     &
   !!!       cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,Qdip_macro,icell,     &
   !!!       macro_mag_trial,macro_trial,de,do_anisotropy)
   !!! 
   !!!       use Constants, only : mub
   !!!       use macrocells, only : calc_trial_macro_mom
   !!!       use HamiltonianActions_lite
   !!! 
   !!!       !.. Implicit declarations
   !!!       implicit none
   !!! 
   !!!       ! System variables
   !!!       integer, intent(in) :: k !< Current ensemble
   !!!       integer, intent(in) :: Natom !< Number of atoms in system
   !!!       integer, intent(in) :: Mensemble !< Number of ensembles
   !!!       ! LSF variables
   !!!       integer, intent(in) :: conf_num   !< Number of configurations for LSF
   !!!       character(len=1), intent(in) :: exc_inter !> Interpolation of Jij between FM/DLM (Y/N)
   !!!       ! Heisenberg exchange variables
   !!!       integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
   !!!       integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
   !!!       integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
   !!!       real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
   !!!       real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
   !!!       ! DMI  variables
   !!!       integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
   !!!       integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
   !!!       integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
   !!!       integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
   !!!       real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
   !!!       ! PD interactions variables
   !!!       integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
   !!!       integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
   !!!       integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
   !!!       integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
   !!!       real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
   !!!       ! BIQDM variables
   !!!       integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
   !!!       integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
   !!!       integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
   !!!       integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
   !!!       real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange coupling
   !!!       ! BQ variables
   !!!       integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
   !!!       integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
   !!!       integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
   !!!       integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
   !!!       real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
   !!!       ! Anisotropy variables
   !!!       integer, intent(in) :: do_anisotropy
   !!!       integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
   !!!       integer, dimension(Natom),intent(in) :: taniso_diff !< Type of anisotropy (0-2)
   !!!       real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
   !!!       real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
   !!!       real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
   !!!       real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
   !!!       real(dblprec), dimension(Natom), intent(in) :: sb!< Ratio between the Anisotropy constants
   !!!       real(dblprec), dimension(Natom), intent(in) :: sb_diff!< Ratio between the Anisotropy constants
   !!!       character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
   !!!       ! Moments variables
   !!!       integer, intent(in) :: iflip !< Atom to flip spin for
   !!!       real(dblprec), dimension(3), intent(in) :: newmom !< New trial moment
   !!!       real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
   !!!       real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
   !!!       real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
   !!!       ! External magnetic fields
   !!!       real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
   !!!       ! Dipolar interaction variables
   !!!       integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
   !!!       real(dblprec), dimension(3,3,Natom,Natom), intent(in) :: Qdip !< Matrix for dipole-dipole interaction
   !!!       ! Induced moment variables
   !!!       integer, intent(in) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments
   !!!       integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for induced moments
   !!!       integer, dimension(Natom), intent(in) :: ind_list_full !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
   !!!       integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
   !!!       real(dblprec), dimension(Natom), intent(in) :: sus_ind
   !!!       ! Macrocell variables
   !!!       integer, intent(in) :: Num_macro !< Number of macrocells in the system
   !!!       integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
   !!!       integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
   !!!       integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
   !!!       integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
   !!!       real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
   !!!       real(dblprec), dimension(3,3,Num_macro,Num_macro), intent(in) :: Qdip_macro !< Matrix for macro spin dipole-dipole
   !!!       ! Output variables
   !!!       real(dblprec), intent(out):: de  !< Energy difference
   !!!       integer, intent(in) :: icell
   !!!       real(dblprec), intent(in) :: macro_mag_trial
   !!!       real(dblprec), dimension(3), intent(in) :: macro_trial
   !!! 
   !!!       !.. Local scalars
   !!!       integer :: j
   !!!       integer :: jflip
   !!!       real(dblprec) :: tt, e_c, e_t
   !!!       real(dblprec) :: temp_ene
   !!! 
   !!!       !.. Local arrays
   !!! 
   !!!       !.. Executable statements
   !!! 
   !!!       ! Start with setting up the trial moments ( and copy the "active region" for reference)
   !!!       temp_ene=0.0_dblprec
   !!!       emomM_trial(:,iflip)=newmom(:)*mmom(iflip,k)
   !!!       !emomM(:,iflip,k)=newmom(:)*mmom(iflip,k)
   !!!       do j=1,ind_nlistsize(iflip)
   !!!          jflip=ind_nlist(j,iflip)
   !!!          emomM_trial(:,jflip)=emomM(:,jflip,k)-sus_ind(jflip)*emomM_trial(:,iflip)+sus_ind(jflip)*emomM(:,iflip,k)
   !!!          !emomM_trial(:,jflip)=emomM(:,jflip,k)
   !!!       end do
   !!! 
   !!! 
   !!!       e_t=0.0_dblprec
   !!!       !print '(2x,2i8,4f12.6)',iflip,iflip,newmom
   !!!       !print '(2x,2i8,4f12.6)',iflip,iflip,emomM(:,iflip,1),temp_ene
   !!!       call effective_field_extralite(Natom,Mensemble,iflip,iflip,emomM_trial,temp_ene,beff)
   !!!       e_t=e_t+temp_ene
   !!!       !print '(2x,2i8,4f12.6)',iflip,iflip,beff(:,iflip,1),e_t
   !!!       do j=1,ind_nlistsize(iflip)
   !!!          jflip=ind_nlist(j,iflip)
   !!!       !   print '(2x,2i8,4f12.6)',iflip,jflip,emomM(:,jflip,1),e_t
   !!!          call effective_field_extralite(Natom,Mensemble,jflip,jflip,emomM_trial,temp_ene,beff)
   !!!          e_t=e_t+temp_ene
   !!!       !   print '(2x,2i8,4f12.6)',iflip,jflip,beff(:,jflip,1),e_t
   !!!       end do
   !!!                   
   !!!       e_c=0.0_dblprec
   !!!       !print '(5x,2i8,4f12.6)',iflip,iflip,emomM_trial(:,iflip),temp_ene
   !!!       call effective_field_extralite(Natom,Mensemble,iflip,iflip,emomM,temp_ene,beff)
   !!!       e_c=e_c+temp_ene
   !!!       !print '(5x,2i8,4f12.6)',iflip,iflip,beff(:,iflip,1),e_c
   !!!       do j=1,ind_nlistsize(iflip)
   !!!          jflip=ind_nlist(j,iflip)
   !!!       !print '(5x,2i8,4f12.6)',iflip,jflip,emomM_trial(:,jflip),e_c
   !!!          call effective_field_extralite(Natom,Mensemble,jflip,jflip,emomM,temp_ene,beff)
   !!!          e_c=e_c+temp_ene
   !!!       !print '(5x,2i8,4f12.6)',iflip,jflip,beff(:,jflip,1),e_c
   !!!       end do
   !!! 
   !!!       !Energy difference
   !!!       tt=e_t-e_c
   !!!       de=mub*tt
   !!! 
   !!!       ! Reset trial momenta
   !!!       emomM_trial(:,iflip)=emomM(:,iflip,k)
   !!!       !emomM(:,iflip,k)=newmom(:)*mmom(iflip,k)
   !!!       do j=1,ind_nlistsize(iflip)
   !!!          jflip=ind_nlist(j,iflip)
   !!!          emomM_trial(:,jflip)=emomM(:,jflip,k)
   !!!       end do
   !!!       return
   !!!    end subroutine calculate_energy_wIND_v2


   !---------------------------------------------------------------------------
   !> @brief
   !> Calculate the total energy of a single fixed magnetic moment plus including the
   !> indirect energy contributions of the induced magnetic moments
   !>
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine calculate_energy_wIND_v3(k,Natom,Mensemble,max_no_neigh,nlistsize,nlist, &
         ncoup,conf_num,iflip,newmom,mmom,emomM,extfield,        &
         ind_nlistsize,ind_nlist,ind_list_full,sus_ind,max_no_neigh_ind, de)

      use Constants, only : mub
      use macrocells, only : calc_trial_macro_mom

      !.. Implicit declarations
      implicit none

      ! System variables
      integer, intent(in) :: k !< Current ensemble
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      ! LSF variables
      integer, intent(in) :: conf_num   !< Number of configurations for LSF
      ! Heisenberg exchange variables
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      ! Moments variables
      integer, intent(in) :: iflip !< Atom to flip spin for
      real(dblprec), dimension(3), intent(in) :: newmom !< New trial moment
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! External magnetic fields
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      !
      integer, intent(in) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments
      integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(Natom), intent(in) :: ind_list_full !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      ! Output variables
      real(dblprec), intent(out):: de  !< Energy difference

      !.. Local scalars
      integer :: j, neigh_test, ind_neigh
      integer :: ja, l, la
      real(dblprec) :: tt, e_c, e_t
      real(dblprec) :: diff_e
      real(dblprec) :: dot, mmom_diff

      !.. Local arrays
      real(dblprec), dimension(3) :: beff_t, trialmom, ave_mom, diff_mom

      neigh_test=0
      !.. Executable statements

      ! First calculate effective field
      beff_t(1) = 0_dblprec
      beff_t(2) = 0_dblprec
      beff_t(3) = 0_dblprec
      tt=0.0_dblprec

      e_c=0.0_dblprec
      e_t=0.0_dblprec
      trialmom(:)=newmom(:)*mmom(iflip,k)
      diff_mom(:)=trialmom(:)-emomM(:,iflip,k)

      diff_e=0.0_dblprec

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! The iflip magnetic moment should be only fixed magnetic moments, if the
      ! neigbour is fixed everything procedes as usual, if the neighbour is induced
      ! then one must rotate the induced moment according to the direction of its neighbouring
      ! fixed moments, and then the energy of that new configuration must be calculated
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !if (exc_inter=='N') then
      ! Exchange
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculation of the exchange term
      ! Regular interaction (compensate for induced later) 
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
      !    !$omp simd reduction(+:diff_e) private(dot)
#endif
      do j=1,nlistsize(iflip)
         ja=nlist(j,iflip)
         dot=diff_mom(1)*emomM(1,ja,k)+diff_mom(2)*emomM(2,ja,k)+diff_mom(3)*emomM(3,ja,k)
         diff_e=diff_e-ncoup(j,iflip,1)*dot
         !dot=trialmom(1)*emomM(1,ja,k)+trialmom(2)*emomM(2,ja,k)+trialmom(3)*emomM(3,ja,k)
         !e_c=e_c+ncoup(j,iflip,1)*dot
         !dot=emomM(1,iflip,k)*emomM(1,ja,k)+emomM(2,iflip,k)*emomM(2,ja,k)+emomM(3,iflip,k)*emomM(3,ja,k)
         !e_c=e_t+ncoup(j,iflip,1)*dot
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ind_list_full(iflip).eq.0) then
         ! Loop also over all close induced moments
         do j=1,ind_nlistsize(iflip)
            ! Select the current induced moment
            ja=ind_nlist(j,iflip)
            !              if (ind_list_full(ja).eq.1) then
            ! Sum over the nearest neighbours that are fixed
            ave_mom=0.0_dblprec
            do ind_neigh=1, ind_nlistsize(ja)
               ave_mom(1)=ave_mom(1)+emomM(1,ind_nlist(ind_neigh,ja),k)*sus_ind(ja)
               ave_mom(2)=ave_mom(2)+emomM(2,ind_nlist(ind_neigh,ja),k)*sus_ind(ja)
               ave_mom(3)=ave_mom(3)+emomM(3,ind_nlist(ind_neigh,ja),k)*sus_ind(ja)
            enddo
            ! Vary the magnitude of the induced moments
            mmom_diff=sqrt(ave_mom(1)*ave_mom(1)+ave_mom(2)*ave_mom(2)+ave_mom(3)*ave_mom(3))-mmom(ja,k)
            mmom_diff=(mmom_diff/mmom(ja,k)+1.0e-12_dblprec)
            ! First correct for the M->m  coupling from iflip
            dot=emomM(1,ja,k)*emomM(1,iflip,k)+emomM(2,ja,k)*emomM(2,iflip,k)+emomM(3,ja,k)*emomM(3,iflip,k)
            diff_e=diff_e-ncoup(j,iflip,1)*dot*mmom_diff
            ! Then loop over all m->M couplings
#if _OPENMP >= 201307 && ( ! defined __INTEL_COMPILER_BUILD_DATE || __INTEL_COMPILER_BUILD_DATE > 20140422)
            !           !$omp simd reduction(+:diff_e) private(dot)
#endif
            do l=1,nlistsize(ja)
               la=nlist(l,ja)
               dot=emomM(1,ja,k)*emomM(1,la,k)+emomM(2,ja,k)*emomM(2,la,k)+emomM(3,ja,k)*emomM(3,la,k)
               diff_e=diff_e-ncoup(l,ja,1)*dot*mmom_diff
            end do
            !              end if
         end do
      end if

      ! Add Zeeman term
      !e_c=e_c-extfield(1)*emomM(1,iflip,k)-extfield(2)*emomM(2,iflip,k)-extfield(3)*emomM(3,iflip,k)
      !e_t=e_t-extfield(1)*trialmom(1)-extfield(2)*trialmom(2)-extfield(3)*trialmom(3)
      diff_e=diff_e-extfield(1)*diff_mom(1)-extfield(2)*diff_mom(2)-extfield(3)*diff_mom(3)

      !Energy difference
      !tt=e_t-e_c
      !de=mub*tt
      de=mub*diff_e
      return
   end subroutine calculate_energy_wIND_v3

end module InducedMoments
