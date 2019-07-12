!-------------------------------------------------------------------------------
! MODULE: SpinIce
!> @brief Wrapper module of the Loop algorithm to treat spin ice systems
!
!> @author Anders Bergman and Jonathan Chico
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module SpinIce

   use Parameters
   use Constants
   use SpinIceData
   use RandomNumbers, only: rng_uniform,rng_gaussian, use_vsl

   implicit none

   integer, dimension(:), allocatable :: nmdimt_ver_ver ! Temporary storage of dimensions for vertex-vertex neighbour map
   integer, dimension(:), allocatable :: nmdimt_ver_atom ! Temporary storage of dimension for atom-vertex neighbour map
   real(dblprec), dimension(:,:,:), allocatable :: sym_mats !< Symmetry operation matrices
   integer :: nsym !< Number of symmetry operations

   private :: nmdimt_ver_ver,nmdimt_ver_atom

   public :: ice_rule, select_cluster, setup_neighbour_list
   public :: setup_ver_nm, setup_ver_atom_nm, mc_update_spinice

contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: mc_update_spinice
   !> @brief Driver for spin ice Monte Carlo
   !> @author Jonathan Chico and Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine mc_update_spinice(Natom,Mensemble,nHam,max_no_neigh,conf_num,ncoup,ncoupD,&
      nlist,nlistsize,aham,do_dm,max_no_dmneigh,dm_vect,dmlist,dmlistsize,&
      do_pd,nn_pd_tot,pd_vect,pdlist,pdlistsize,&
      do_biqdm,nn_biqdm_tot,biqdm_vect,biqdmlist,biqdmlistsize,&
      do_bq,nn_bq_tot,j_bq,bqlist,bqlistsize,&
      do_chir,nn_chir_tot,chir_coup,chirlist,chirlistsize,&
      taniso,eaniso,kaniso,sb,emomM,emom,mmom,iflip_a,extfield,&
      mult_axis,taniso_diff,eaniso_diff,kaniso_diff,sb_diff,&
      do_dip,Qdip,exc_inter,temperature,temprescale,&
      ind_nlistsize,ind_nlist,sus_ind,ind_mom_flag,ind_list_full,max_no_neigh_ind,&
      Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,macro_atom_nlist,&
      mmom_macro,emom_macro,emomM_macro,Qdip_macro,do_anisotropy)

      use montecarlo_common
      !
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: nHam !< Number of atoms in Hamiltonian
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: conf_num  !< number of configurations for LSF
      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(nHam),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: aham !< Hamiltonian look-up table
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
      integer, intent(in) :: do_chir   !< Add Pseudo-Dipolar (chir) term to Hamiltonian (0/1)
      integer, intent(in) :: nn_chir_tot !< Calculated number of neighbours with chir interactions
      real(dblprec), dimension(nn_chir_tot,Natom), intent(in) :: chir_coup !< Pseudo-Dipolar exchange vector
      integer, dimension(nn_chir_tot,Natom), intent(in) :: chirlist   !< List of neighbours forchir 
      integer, dimension(Natom),intent(in) :: chirlistsize !< Size of neighbour list forchir 
      integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
      integer, dimension(Natom),intent(in) :: taniso_diff !< Type of anisotropy (0-2)
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
      real(dblprec), dimension(Natom), intent(in) :: sb !< Ratio between anisotropy constants
      real(dblprec), dimension(Natom), intent(in) :: sb_diff !< Ratio between anisotropy constants
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: temperature !< Temperature
      real(dblprec), intent(in) :: temprescale !< Temperature rescaling according to QHB
      real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
      integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
      real(dblprec), dimension(3,3,Natom,Natom), intent(in) :: Qdip !< Matrix for dipole-dipole interaction
      character(len=1), intent(in) :: mult_axis
      character(len=1), intent(in)  ::  exc_inter !< Exchange interpolation between FM/DLM (Y/N)
      integer,dimension(Natom),intent(in) :: iflip_a !< Flipping pattern
      integer, intent(in) :: max_no_neigh_ind !< Calculated maximum of neighbours for induced moments
      integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(Natom),intent(in) :: ind_list_full
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      character(len=1), intent(in) :: ind_mom_flag

      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
      real(dblprec), dimension(Num_macro,Mensemble), intent(inout) :: mmom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emom_macro
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(inout) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,3,Num_macro,Num_macro), intent(in) :: Qdip_macro !< Matrix for macro spin dipole-dipole
      integer, intent(in) :: do_anisotropy

      !.. Local scalars
      integer :: i, k, loop_len, ii, ip
      integer :: iatom,icell
      real(dblprec) :: de, delt  !< Energy difference
      real(dblprec) :: macro_mag_trial

      real(dblprec) :: ve, ice_count,ice_count_up,ice_count_down !< Trial vertex energy

      !.. Local arrays
      real(dblprec), dimension(3) :: newmom !< New trial moment
      real(dblprec), dimension(Natom) :: flipprob_a ,newmmom_a
      real(dblprec), dimension(3,Natom) :: newmom_a
      real(dblprec), dimension(3) :: macro_trial

      real(dblprec), dimension(3,Natom,Mensemble) :: emomM_copy  !< Copy of current magnetic moment vector

      integer, dimension(Nvertex) :: Ncount_ver
      integer, dimension(Natom) :: Ncount_atom

      !AB
      logical, dimension(Nvertex) :: vertex_map
      logical :: flipped

      Ncount_atom=0
      Ncount_ver=0
      ver_no=0
      loop_ave_len=0.0_dblprec

      do k=1,Mensemble
         ! First Ising random atom
         !$omp parallel do default(shared),private(ip,newmom),schedule(static,1)
         do ip=1,Natom
            call Ising_random_flip(emom,newmom,iflip_a(ip),k,Natom,Mensemble)
            newmom_a(:,iflip_a(ip))=newmom(:)
         enddo
         !$omp end parallel do
         call rng_uniform(flipprob_a,natom)
         ! Calculate energy and flip spin if preferred
         !$omp parallel do default(shared), private(i,de), schedule(static,1)
         do i=1, Natom
            call  calculate_energy(Natom, Mensemble, nHam, conf_num, do_dm , do_pd, do_biqdm, do_bq, do_chir,&
            emomM, emom, mmom, iflip_a(i),newmom_a(1:3,iflip_a(i)), extfield, de, k, &
            mult_axis, do_dip,Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,&
            macro_atom_nlist,emomM_macro,icell,macro_mag_trial,macro_trial, exc_inter,do_anisotropy)

            call flip_a(Natom, Mensemble, emom, emomM, mmom, iflip_a(i),&
            newmom_a(1:3,iflip_a(i)),newmmom_a(iflip_a(i)),&
            de,temperature,temprescale,'N',k,flipprob_a(i),1,ind_nlistsize,&
            ind_nlist,ind_mom_flag,max_no_neigh_ind,sus_ind,ind_list_full,&
            do_dip,Num_macro,icell,macro_mag_trial,macro_trial,mmom_macro,emom_macro,emomM_macro)
         enddo
         !$omp end parallel do


         do i=1,Nvertex
            vertex_map(i)=ice_rule(Natom,Mensemble,k,emom,i,ve,ice_count,ice_count_up,ice_count_down)
         end do
         !print '(32L1)',vertex_map

         do iatom=1, Natom
            !call select_cluster(Natom,Mensemble,emom, Ncount_atom, loop_len, k)
            call select_loop(Natom,Mensemble,emom, Ncount_atom, loop_len, k,Nvertex,vertex_map)

            ! Calculate the energy difference that is obtained when one flips the loop
            emomM_copy=emomM

            if (loop_len.gt.0) then
               call calculate_spinice_energy(Natom,Mensemble,nHam,max_no_neigh,conf_num,ncoup,nlist,nlistsize,aham,&
               do_dm,max_no_dmneigh,dm_vect,dmlist,dmlistsize, &
               do_pd,nn_pd_tot,pd_vect,pdlist,pdlistsize,&
               do_biqdm, nn_biqdm_tot,biqdm_vect,biqdmlist,biqdmlistsize,&
               do_bq,nn_bq_tot,j_bq,bqlist,bqlistsize,&
               taniso,eaniso,kaniso,sb,emomM_copy,emom,&
               Ncount_atom,loop_len,extfield,delt,k,do_anisotropy,do_dip,Qdip)

               ! If the change is favorable one flips the spins

               call Spin_Ice_flip(emom, emomM, mmom, Ncount_atom, Natom, Mensemble, Temperature, delt, k,loop_len,flipped)
            endif

            ! Do loop over Natom to measure the average loop length
            if (mod(loop_len,2)>dbl_tolerance) write(*,*) 'WARNING'
            loop_ave_len=loop_ave_len+loop_len
         end do
         !do ii=1,Nvertex
         !  if (ice_rule(Natom,Mensemble,k,emom,ii,ve,ice_count,ice_count_up,ice_count_down)) ver_no=ver_no+1
         !end do
         if (flipped) then
            do ii=1,Nvertex
               do iatom=1,nlistsize_ver_atom(ii)
                  do i=1,loop_len
                     if(Ncount_atom(i)==iatom) then
                        vertex_map(ii)=ice_rule(Natom,Mensemble,k,emom,i,ve,ice_count,ice_count_up,ice_count_down)
                     end if
                  end do
               end do
            end do
         end if
      enddo
      loop_ave_len=loop_ave_len/(Natom*Mensemble)
      !   print *,mchits_spin_ice

   end subroutine mc_update_spinice


   ! This routine creates a list of spins that will be flipped in acordance with the spin-ice rules
   ! The list must be designed in a smart way such that one does not have to go trough all the spins everytime one selects a random one if not the process is of the order of N2 which is not efficient
   subroutine select_loop(Natom, Mensemble,emom, Ncount_atom, loop_len, k,Nvertex, vertex_map)

      implicit none

      integer, intent(in) :: k ! Current ensemble
      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      integer, intent(in) :: Nvertex
      integer, dimension(Natom), intent(inout) :: Ncount_atom

      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom

      integer, intent(out) :: loop_len ! length of spin ice loop
      logical, dimension(Nvertex), intent(in) :: vertex_map

      !AB
      logical :: loop_closed, atom_selected, failure, total_failure, vertex_selected

      integer :: start_vertex, current_vertex, next_vertex, prev_vertex, flip_atom
      integer :: flip_atom_list_num, prev_flip_atom,prev_flip_atom_num_list
      !
      integer, dimension(:), allocatable :: visited_vertices, snake, atomsnake ! This are the vertices that are visited and form the cluster
      logical, dimension(:), allocatable :: occupied
      integer :: snaketail, atomlen, i, j,snakehead  ! This are the indexes used for the cluster

      real(dblprec) :: dflip(1)
      integer :: ntrial, npossible
      integer, dimension(:), allocatable :: shortlist

      ! Logical variables for failure of the algorithm
      total_failure=.false.
      failure=.false.
      atom_selected=.false.
      loop_closed=.false.

      ! Allocating the vertex snake and the atom snake (clusters for atoms and vertices)
      !allocate(snake(Nvertex))
      allocate(atomsnake(Nvertex+1))
      !allocate(visited_vertices(Nvertex))
      allocate(occupied(Nvertex))
      allocate(shortlist(maxval(nlistsize_ver)))
      visited_vertices=0
      start_vertex=1

      ! Do select the starting vertex and afterwards try to close the vertex until this is done
      total_failure=all(.not.vertex_map)
      ! An array for keeping track of valid tries
      do while (.not.loop_closed .and. .not. total_failure)
         snakehead=1
         atomlen=0
         occupied(:)=.false.
         ! If no vertex fulfills the ice rule, declare immediate failure
         failure=all(.not.vertex_map)
         !failure=.false.

         ! Do loop until the staring vertex is found (the failure test should be moot here)
         ntrial=0
         do while (.not. vertex_selected .and. .not. failure)

            ntrial=ntrial+1
            call rng_uniform(dflip,1)
            start_vertex=ceiling(Nvertex*dflip(1))
            vertex_selected=vertex_map(start_vertex)
            if(vertex_selected) then
               current_vertex=start_vertex
               prev_vertex=start_vertex
            end if
            ! Flag for failure in too many tries have been performed (could be looser)
            failure=ntrial>Nvertex
            occupied(start_vertex)=.true.
         end do

         if (vertex_selected) then
            do while (.not. loop_closed .and. .not.failure)
               ! Find next possible vertex
               npossible=0
               shortlist=0
               loop_closed=.false.
               ! If not closed select a valid neighbour vertex
               npossible=0
               shortlist=0
               do i=1,nlistsize_ver(current_vertex)
                  ! Currently not caring about ice rules when making the loop. Uncomment below for ice-rule compliance.
                  !if(vertex_map(nlist_ver(i,current_vertex)).and..not.nlist_ver(i,current_vertex)==prev_vertex) then
                  if(.not.nlist_ver(i,current_vertex)==prev_vertex) then
                     npossible=npossible+1
                     shortlist(npossible)=nlist_ver(i,current_vertex)
                  end if
               end do
               call rng_uniform(dflip,1)
               next_vertex=shortlist(ceiling(dflip(1)*npossible))

               ! Flag for failure if no vertices are available
               failure=(npossible==0)

               if(.not. failure) then
                  ! Find atom inbetween vertices
                  do i=1, nlistsize_ver_atom(current_vertex)
                     do j=1, nlistsize_ver_atom(next_vertex)
                        ! Flip atom is the number of the spin that will be flipped and stops the count
                        if(nlist_ver_atom(i,current_vertex).eq.nlist_ver_atom(j,next_vertex).and..not.atom_selected ) then
                           flip_atom=nlist_ver_atom(j,next_vertex)  ! Atom number to be flipped
                           occupied(current_vertex)=.true.          ! Mark first vertex as occupied (for short list)
                           prev_vertex=current_vertex               ! Store previous vertex
                           current_vertex=next_vertex               ! Store the trial vertex index to the output
                           atom_selected=.true.
                           loop_closed=occupied(next_vertex)        ! Mark loop as closed if same vortex is passed twice
                        endif
                     end do
                  end do

                  ! If there was an atom selected from the select_first_atom routine one increases the atomsnake array size
                  if (atom_selected) then
                     atomlen=atomlen+1
                     atomsnake(atomlen)=flip_atom  ! Store the index of the atom that is marked in the atomsnake
                     atom_selected=.false.         ! Set the atom selection criteria to false
                     !               print '(a,2000i5)','test > ',atomsnake(1:atomlen)
                  endif

               end if

               ! Find dangling tail if loop is closed
               if(loop_closed) then
                  do i=1,atomlen-1
                     do j=1,nlistsize_ver_atom(next_vertex)
                        if(atomsnake(i)==nlist_ver_atom(j,next_vertex)) then
                           snakehead=i
                        end if
                     end do
                  end do
               end if
               ! If same atom in vertex was passed two, remove the redundant entry.
               if(atomsnake(snakehead)==atomsnake(atomlen)) snakehead=snakehead+1
            end do

         end if

         total_failure=total_failure.or.failure


      enddo

      ! If there is no total failure (it was possible to create and close a loop)
      if(.not.total_failure) then
         loop_len=0                             ! Set the length of the loop to zero
         do i=snakehead,atomlen                 ! Do a loop from the snakehead (the index of the snake that corresponds to the closing vertex)
            loop_len=loop_len+1                  ! Increase the size of the loop
            Ncount_atom(loop_len)=atomsnake(i)   ! Store the indexes of the atoms that are marked for flipping
         end do
         do i=loop_len+1, Natom
            Ncount_atom(i)=0
         enddo
      else
         Ncount_atom=0
         loop_len=0
      end if

      deallocate(atomsnake)
      deallocate(occupied)
      deallocate(shortlist)

      return
   end subroutine select_loop

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: select_cluster
   !> @brief This routine creates a list of spins that will be flipped in acordance with the spin-ice rules
   !> The list must be designed in a smart way such that one does not have to go trough all the spins everytime one selects a random one if not the process is of the order of N2 which is not efficient
   !> @author Anders Bergman and Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine select_cluster(Natom, Mensemble,emom, Ncount_atom, loop_len, k)

      implicit none

      integer, intent(in) :: k ! Current ensemble
      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble
      integer, dimension(Natom), intent(inout) :: Ncount_atom
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
      integer, intent(out) :: loop_len ! length of spin ice loop

      !AB
      logical :: loop_closed, atom_selected, failure, total_failure
      integer :: start_vertex, current_vertex, next_vertex, prev_vertex, flip_atom
      integer :: flip_atom_list_num, prev_flip_atom,prev_flip_atom_num_list
      !
      integer, dimension(:), allocatable :: visited_vertices, snake, atomsnake ! This are the vertices that are visited and form the cluster
      integer :: snaketail, atomlen, i, snakehead  ! This are the indexes used for the cluster

      ! Logical variables for failure of the algorithm
      total_failure=.false.
      failure=.false.
      atom_selected=.false.
      loop_closed=.false.

      ! Allocating the vertex snake and the atom snake (clusters for atoms and vertices)
      allocate(snake(Nvertex))
      allocate(atomsnake(Nvertex+1))
      allocate(visited_vertices(Nvertex))
      visited_vertices=0
      start_vertex=1

      ! Do select the starting vertex and afterwards try to close the vertex until this is done
      do while (.not.loop_closed .and. .not. total_failure)
         snake=0
         atomsnake=0
         snaketail=0
         snakehead=1
         atomlen=0
         ! Do loop until the staring vertex is found
         ! Try to find the first two vertices while the atom_selected is false and the failure flag is also false
         do while (.not. atom_selected .and. .not. failure)

            call select_first_atom(Natom,Nvertex,Mensemble, max_no_equiv_ver,max_no_equiv_ver_atom,nlist_ver,&
            nlist_ver_atom, nlistsize_ver, nlistsize_ver_atom, start_vertex, current_vertex,flip_atom,&
            flip_atom_list_num, atom_selected,&
            failure,k,emom,visited_vertices,vert_ice_coord)
         end do

         ! If there was an atom selected from the select_first_atom routine one increases the atomsnake array size
         if (atom_selected) then
            atomlen=atomlen+1
            atomsnake(atomlen)=flip_atom  ! Store the index of the atom that is marked in the atomsnake
            atom_selected=.false.         ! Set the atom selection criteria to false
         endif

         snaketail=snaketail+1           ! Increase the index for the tail of the snake
         snake(snaketail)=start_vertex   ! Store in the snake the tail (starting vertex)
         snaketail=snaketail+1           ! Increase the index to include the neighbouring vertex that has been found
         snake(snaketail)=current_vertex ! Store in the snake the neighbouring vertex
         prev_vertex=start_vertex        ! One renames the starting vertex to prev_vertex

         ! Try to select a new vertex while the closing and failure flags are both false
         do while (.not. loop_closed .and. .not.failure)

            prev_flip_atom=atomsnake(atomlen)
            prev_flip_atom_num_list=flip_atom_list_num
            call select_new_vertex(Natom,Nvertex,Mensemble,max_no_equiv_ver,max_no_equiv_ver_atom,nlist_ver,&
            nlist_ver_atom,nlistsize_ver, nlistsize_ver_atom,prev_vertex,current_vertex,next_vertex,&
            flip_atom,prev_flip_atom,prev_flip_atom_num_list,flip_atom_list_num,&
            atom_selected,failure,k,emom,vert_ice_coord)

            ! Failure is false then increase the size of the atomic snake (Failure is an output from select_new_vertex)
            if(.not.failure) then
               ! If there was a atom selected in the select_first_atom one starts the cluster
               if(atom_selected) then
                  atomlen=atomlen+1            ! Increase the size of the atomic snake
                  atomsnake(atomlen)=flip_atom ! Save the index of the selected atom in the atomsnake
                  flip_atom=1                  ! Set the flip atom to 1
                  atom_selected=.false.        ! Set the atom selected to false
               end if

               ! Do a loop over the snake to see if the loop should be closed
               do i=1,snaketail
                  ! If there is a vertex in the snake that is equal to the next vertex the loop is closed
                  if(snake(i)==next_vertex) then
                     loop_closed=.true.  ! If the next selected vertex is one that has been selected before the loop is closed
                     snakehead=i         ! The head of the snake is that newly selected vertex
                  end if

               end do
               ! If the loop is not closed
               if(.not.loop_closed) then
                  ! Increase the size of the vertex snake
                  snaketail=snaketail+1
                  snake(snaketail)=next_vertex ! Save in the vertex snake the most recently selected vertex
                  prev_vertex=current_vertex   ! The current vertex (not the most recently selected, the previous one) is saved as prev_vertex
                  current_vertex=next_vertex   ! The most recently selected vertex is now the current vertex
               end if

            else
               visited_vertices(next_vertex)=1 ! Flag the visited vertices
            end if

         end do
         ! Try to see if there is a total failure in the algorithm, that is if the ammount of visited vertices greater or equal to Nvertex
         failure=.false.
         total_failure=total_failure.or.(sum(visited_vertices)>=Nvertex)

      enddo

      ! If there is no total failure (it was possible to create and close a loop)
      if(.not.total_failure) then
         loop_len=0                             ! Set the length of the loop to zero
         do i=snakehead,atomlen                 ! Do a loop from the snakehead (the index of the snake that corresponds to the closing vertex)
            loop_len=loop_len+1                  ! Increase the size of the loop
            Ncount_atom(loop_len)=atomsnake(i)   ! Store the indexes of the atoms that are marked for flipping
         end do
         do i=loop_len+1, Natom
            Ncount_atom(i)=0
         enddo
      else
         Ncount_atom=0
         loop_len=0
      end if

      return
   end subroutine select_cluster

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: select_first_atom
   !> @brief This subroutine selects the first vertex by checking if it follows
   !> the ice rules and if it has neighbours that follows the ice rules
   !> @author Anders Bergman and Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine select_first_atom(Natom,Nvertex,Mensemble,max_no_equiv_ver,max_no_equiv_ver_atom,nlist_ver,&
      nlist_ver_atom,nlistsize_ver,nlistsize_ver_atom,start_vertex,current_vertex,flip_atom,&
      flip_atom_num_list,atom_selected,failure,en,emom,vertex_hit,vert_ice_coord)

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: Nvertex
      integer, intent(in) :: Mensemble
      integer, intent(in) :: en ! Current ensemble
      integer, intent(in) :: max_no_equiv_ver
      integer, intent(in) :: max_no_equiv_ver_atom
      integer, dimension(max_no_equiv_ver,Nvertex), intent(in) :: nlist_ver ! Neighbour list for vertex vertex or vertex atoms
      integer, dimension(Nvertex), intent(in):: nlistsize_ver ! Size of neighbour list for vertex vertex or vertex atoms
      integer, dimension(max_no_equiv_ver_atom,Nvertex) , intent(in) :: nlist_ver_atom
      integer, dimension(Nvertex), intent(in) :: nlistsize_ver_atom

      integer, intent(out) :: start_vertex
      integer, intent(out) :: current_vertex
      integer, intent(out) :: flip_atom
      integer, intent(out) :: flip_atom_num_list

      logical, intent(inout) :: atom_selected ! Logical variable to see if there is an atom selected
      logical, intent(out) :: failure         ! Logical variable to see if there was a failure in the algorithm
      real(dblprec), dimension(3,Nvertex,max_no_equiv_ver_atom), intent(in) :: vert_ice_coord
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
      integer, dimension(Nvertex), intent(inout) :: vertex_hit   !List to store visited vertices

      integer, dimension(:), allocatable :: nvertex_hit   ! List to store visited neighbour vertices

      real(dblprec) :: dflip(1), ve, ice_count,ice_count_up,ice_count_down ! ve is the vertex total moment

      integer :: i, l, k
      integer :: trial_vert, trial_nvert

      !AB
      ! This array was introduced by Anders to check if all the neighbours have been checked
      allocate(nvertex_hit(maxval(nlistsize_ver)))
      nvertex_hit=0

      ! Find the starting vertex that follows the ice rule and the neighbours follow the ice rule
      do while (.not.atom_selected .and. sum(vertex_hit)<Nvertex)
         call rng_uniform(dflip,1)
         ! Index of the temporally selected vertex
         trial_vert=int(dflip(1)*Nvertex)+1
         ! If the vertex has already been visited and does not fullfill the ice rules try to find a new one
         do while(vertex_hit(trial_vert)==1)
            call rng_uniform(dflip,1)
            trial_vert=int(dflip(1)*Nvertex)+1
         end do

         vertex_hit(trial_vert)=1

         if(ice_rule(Natom,Mensemble,en,emom,trial_vert,ve,ice_count,ice_count_up,ice_count_down))then
            i=0
            nvertex_hit=0
            do while(.not.atom_selected .and. sum(nvertex_hit)<nlistsize_ver(trial_vert))
               ! Select one of the neighbours of the temporal vertex to be the trial vertex
               call rng_uniform(dflip,1)
               ! Index of the temporally selected vertex
               i=int(dflip(1)*nlistsize_ver(trial_vert))+1
               do while(nvertex_hit(i)==1)
                  call rng_uniform(dflip,1)
                  i=int(dflip(1)*nlistsize_ver(trial_vert))+1
               end do
               nvertex_hit(i)=1
               trial_nvert=nlist_ver(i,trial_vert)

               if(ice_rule(Natom,Mensemble,en,emom,trial_nvert,ve,ice_count,ice_count_up,ice_count_down)&
               .and.(nlistsize_ver(trial_nvert)-1)>0)then
               ! Check the neighbours of the two selected vertices and select the shared atom
               do l=1, nlistsize_ver_atom(trial_nvert)
                  do k=1, nlistsize_ver_atom(trial_nvert)
                     ! Flip atom is the number of the spin that will be flipped and stops the count
                     if(nlist_ver_atom(l,trial_vert).eq.nlist_ver_atom(k,trial_nvert) ) then
                        flip_atom=nlist_ver_atom(k,trial_nvert)  ! Atom number to be flipped
                        current_vertex=trial_nvert               ! Store the trial vertex index to the output
                        start_vertex=trial_vert                  ! Store the starting vertex index to the output
                        atom_selected=.true.
                        flip_atom_num_list=k
                     endif
                  end do
               end do
            end if
         enddo
      end if
   enddo

   failure=failure.or..not.atom_selected

   deallocate(nvertex_hit)

   return
end subroutine select_first_atom

!-----------------------------------------------------------------------------
! SUBROUTINE: select_new_vertex
!> @brief This subroutine selects the first vertex by checking if it follows the ice
!> rules and if it has neighbours that follows the ice rules
!> @author Jonathan Chico and Anders Bergman
!-----------------------------------------------------------------------------
subroutine select_new_vertex(Natom,Nvertex,Mensemble,max_no_equiv_ver,max_no_equiv_ver_atom,nlist_ver,&
   nlist_ver_atom, nlistsize_ver, nlistsize_ver_atom, prev_vertex, current_vertex, next_vertex,&
   flip_atom, prev_flip_atom, prev_flip_atom_list_num, flip_atom_list_num,&
   atom_selected,failure,en,emom, vert_ice_coord)

   implicit none

   integer, intent(in) :: Natom
   integer, intent(in) :: Nvertex
   integer, intent(in) :: Mensemble
   integer, intent(in) :: en ! Current ensemble
   integer, intent(in) :: max_no_equiv_ver
   integer, intent(in) :: max_no_equiv_ver_atom
   integer, intent(in) :: prev_flip_atom ! previous atom marked to be flipped
   integer, intent(in) :: prev_flip_atom_list_num
   integer, dimension(max_no_equiv_ver,Nvertex), intent(in) :: nlist_ver ! Neighbour list for vertex vertex or vertex atoms
   integer, dimension(Nvertex), intent(in):: nlistsize_ver ! Size of neighbour list for vertex vertex or vertex atoms
   integer, dimension(max_no_equiv_ver_atom,Nvertex) , intent(in) :: nlist_ver_atom
   integer, dimension(Nvertex), intent(in) :: nlistsize_ver_atom
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
   real(dblprec), dimension(3,Nvertex,max_no_equiv_ver_atom), intent(in) :: vert_ice_coord
   integer, intent(in) :: prev_vertex
   integer, intent(in) :: current_vertex
   integer, intent(inout) :: next_vertex
   integer, intent(out) :: flip_atom
   integer, intent(out) :: flip_atom_list_num

   logical, intent(inout) :: atom_selected
   logical, intent(out) :: failure

   integer, dimension(:), allocatable :: nvertex_hit   !List to store visited neighbour vertices

   real(dblprec) :: dflip(1), ve, ice_count,ice_count_up,ice_count_down
   real(dblprec) :: temp_flip_sign,temp_prev_flip_sign

   integer :: i, l, k

   !AB
   allocate(nvertex_hit(maxval(nlistsize_ver)))
   nvertex_hit=0

   ! Check every neighbours of the temporal vertex
   i=0

   do while(.not.atom_selected .and. sum(nvertex_hit)<nlistsize_ver(current_vertex)-1)
      ! Select one of the neighbours of the temporal vertex to be the trial vertex
      call rng_uniform(dflip,1)
      ! Index of the temporally selected vertex
      i=int(dflip(1)*nlistsize_ver(current_vertex))+1
      ! Check if this vertex has been visited in trials or if it is the previous vertex
      ! Why should nlist_ver(i,current_vertex) should be less than 0
      do while(nvertex_hit(i)==1.or.nlist_ver(i,current_vertex)<0.or.nlist_ver(i,current_vertex)==prev_vertex)
         ! Try to find a new trial vertex
         call rng_uniform(dflip,1)
         i=int(dflip(1)*nlistsize_ver(current_vertex))+1
      end do
      ! Put the vertex hit to one to be aware if it is visited again
      nvertex_hit(i)=1
      ! Next vertex is the new trial vertex
      next_vertex=nlist_ver(i,current_vertex)
      ! Check if the ice rule is fullfilled
      if(ice_rule(Natom,Mensemble,en,emom,next_vertex,ve,ice_count,ice_count_up,ice_count_down))then
         ! Check the neighbours of the two selected vertices and select the shared atom
         do l=1, nlistsize_ver_atom(current_vertex)
            do k=1, nlistsize_ver_atom(next_vertex)
               ! Flip atom is the number of the spin that will be flipped and stops the count
               if(nlist_ver_atom(l,current_vertex).eq.nlist_ver_atom(k,next_vertex) ) then
                  temp_flip_sign=emom(1,nlist_ver_atom(l,current_vertex),en)*vert_ice_coord(1,current_vertex,l)+&
                  emom(2,nlist_ver_atom(l,current_vertex),en)*vert_ice_coord(2,current_vertex,l)+&
                  emom(3,nlist_ver_atom(l,current_vertex),en)*vert_ice_coord(3,current_vertex,l)

                  temp_prev_flip_sign=emom(1,prev_flip_atom,en)*vert_ice_coord(1,current_vertex,prev_flip_atom_list_num)+&
                  emom(2,prev_flip_atom,en)*vert_ice_coord(2,current_vertex,prev_flip_atom_list_num)+&
                  emom(3,prev_flip_atom,en)*vert_ice_coord(3,current_vertex,prev_flip_atom_list_num)
                  if ((temp_flip_sign/temp_prev_flip_sign).eq.-1.0_dblprec) then
                     flip_atom=nlist_ver_atom(k,next_vertex) ! Atom number to be flipped
                     flip_atom_list_num=k
                     atom_selected=.true.
                  endif
               endif
            end do
         end do
      end if
   enddo

   failure=failure.or..not.atom_selected

   deallocate(nvertex_hit)

   return
end subroutine select_new_vertex

!-----------------------------------------------------------------------------
! SUBROUTINE: setup_ver_nm
!> @brief Set up vertex-vertex meighbour maps based on the neighbour routines by
!> Anders Bergman
!> @author Jonathan Chico
!-----------------------------------------------------------------------------
subroutine setup_ver_nm(Nvertex,NT_ver,NA_ver,N1,N2,N3,C1_ver,C2_ver,C3_ver,&
   BC1,BC2,BC3,atype,Bas_ver,&
   max_no_equiv_ver, sym, ver_dist_coord, nm_ver_ver, nmdim_ver_ver)
   !
   !
   implicit none
   !
   ! Quatities from the repetitions of the cell that create the square lattice, they should be the same for the atom and vertex lattice
   !
   integer, intent(in) :: Nvertex ! The number of vertices in the system
   integer, intent(in) :: NT_ver ! Number of types of vertices
   integer, intent(in) :: NA_ver  ! Number of vertices in one cell
   integer, intent(in) :: N1  ! Number of cell repetitions in x direction
   integer, intent(in) :: N2  ! Number of cell repetitions in y direction
   integer, intent(in) :: N3  ! Number of cell repetitions in z direction
   real(dblprec), dimension(3), intent(in) :: C1_ver ! First lattice vector of vertex cell
   real(dblprec), dimension(3), intent(in) :: C2_ver ! Second lattice vector of vertex cell
   real(dblprec), dimension(3), intent(in) :: C3_ver ! Third lattice vector of vertex cell
   character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
   character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
   character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
   integer, dimension(Nvertex), intent(in) :: atype !< Type of atom
   integer, intent(in) :: sym !< Symmetry of system (0-3)
   !
   ! Quantites that refer to the vertices only
   !
   real(dblprec), dimension(3,NA_ver), intent(in) :: Bas_ver ! Coordinates for basis vertex elements
   real(dblprec), dimension(NT_ver,3), intent(in) :: ver_dist_coord
   ! This are the ouput variables
   integer, intent(out) :: max_no_equiv_ver ! Calculated maximum of neighbours in one shell
   integer, dimension(:,:), allocatable, intent(out) :: nm_ver_ver ! Vertex-vertex neighbour map
   integer, dimension(:), allocatable, intent(out) :: nmdim_ver_ver ! Dimension of the vertex-vertex neighbour map
   !
   integer :: i0
   integer :: i_stat,i_all
   real(dblprec) :: tol
   real(dblprec) :: detmatrix
   real(dblprec), dimension(3,3) :: invmatrix
   !
   integer :: ix,iy,iz,xc_hop,yc_hop,zc_hop,ia,counter,itype,inei
   integer :: iat,jat,j0,jx,jy,jz
   real(dblprec), dimension(3) :: bsf,cvec,icvec,rvec
   real(dblprec), dimension(:,:,:), allocatable:: nncoord_ver !< Full list of neighbours for each type
   integer, dimension(:,:), allocatable:: nm_cell
   integer, dimension(:,:,:), allocatable :: nm_trunk
   integer, dimension(:), allocatable :: nnm_cell
   logical :: is_periodic,is_dilute
   !

   tol=0.005_dblprec
   ! Calculate inverse of basis matrix
   detmatrix=C1_ver(1)*C2_ver(2)*C3_ver(3)-C1_ver(1)*C2_ver(3)*C3_ver(2)+&
   C1_ver(2)*C2_ver(3)*C3_ver(1)-C1_ver(2)*C2_ver(1)*C3_ver(3)+&
   C1_ver(3)*C2_ver(1)*C3_ver(2)-C1_ver(3)*C2_ver(2)*C3_ver(1)
   invmatrix=0.0_dblprec
   if(detmatrix/=0.0_dblprec) then
      invmatrix(1,1)=(C2_ver(2)*C3_ver(3)-C3_ver(2)*C2_ver(3))/detmatrix
      invmatrix(1,2)=(C1_ver(3)*C3_ver(2)-C3_ver(3)*C1_ver(2))/detmatrix
      invmatrix(1,3)=(C1_ver(2)*C2_ver(3)-C2_ver(2)*C1_ver(3))/detmatrix
      invmatrix(2,1)=(C2_ver(3)*C3_ver(1)-C3_ver(3)*C2_ver(1))/detmatrix
      invmatrix(2,2)=(C1_ver(1)*C3_ver(3)-C3_ver(1)*C1_ver(3))/detmatrix
      invmatrix(2,3)=(C1_ver(3)*C2_ver(1)-C2_ver(3)*C1_ver(1))/detmatrix
      invmatrix(3,1)=(C2_ver(1)*C3_ver(2)-C3_ver(1)*C2_ver(2))/detmatrix
      invmatrix(3,2)=(C1_ver(2)*C3_ver(1)-C3_ver(2)*C1_ver(1))/detmatrix
      invmatrix(3,3)=(C1_ver(1)*C2_ver(2)-C2_ver(1)*C1_ver(2))/detmatrix
   end if

   ! create all symmetry matrices wrt symmetry type
   call get_symops(sym)

   ! Allocate arrays
   ! Maximum no. of neighbours in each shell is determined by the symmetry (max=48)
   max_no_equiv_ver=nsym

   ! neighbour coordinates
   allocate(nncoord_ver(3,max_no_equiv_ver,NT_ver),stat=i_stat)
   call memocc(i_stat,product(shape(nncoord_ver))*kind(nncoord_ver),'nncoord_ver','setup_ver_nm')

   ! neighbour map information 1: atom no in first supercell
   allocate(nm_cell(max_no_equiv_ver,NA_ver),stat=i_stat)
   call memocc(i_stat,product(shape(nm_cell))*kind(nm_cell),'nm_cell','setup_ver_nm')

   ! neighbour map information 2: translation in basis vectors
   ! so that full neighbour vector is nm_cell+nm_trunk
   is_periodic=(BC1=='P'.or.BC2=='P'.or.BC3=='P')
   is_dilute=(N1*N2*N3>1.and.is_periodic)
   if(N1*N2*N3>1) then
      allocate(nm_trunk(3,max_no_equiv_ver,NA_ver),stat=i_stat)
      call memocc(i_stat,product(shape(nm_trunk))*kind(nm_trunk),'nm_trunk','setup_ver_nm')
   end if

   ! max no. neighbours in each shell
   allocate(nnm_cell(NA_ver),stat=i_stat)
   call memocc(i_stat,product(shape(nnm_cell))*kind(nnm_cell),'nnm_cell','setup_ver_nm')

   ! Allocate arrays for neighbour map dimensions
   allocate(nmdim_ver_ver(Nvertex),stat=i_stat)
   call memocc(i_stat,product(shape(nmdim_ver_ver))*kind(nmdim_ver_ver),'nmdim_ver_ver','setup_ver_nm')
   allocate(nmdimt_ver_ver(NA_ver),stat=i_stat)
   call memocc(i_stat,product(shape(nmdimt_ver_ver))*kind(nmdimt_ver_ver),'nmdimt_ver_ver','setup_ver_nm')
   nmdim_ver_ver=0
   nmdimt_ver_ver=0

   ! Create full vertex vertex neighbour list according to symmetry
   call get_full_ver_ver_list(NT_ver, ver_dist_coord,max_no_equiv_ver, nncoord_ver)

   ! Start looking in "first cell"
   !$omp parallel do default(shared) private(i0,itype,counter,inei,cvec,icvec,bsf,rvec,ia)
   do I0=1,NA_ver
      itype=atype(I0)
      ! Shell
      counter=0
      ! Symmetry equivalent sites in shell
      do inei=1,nmdimt_ver_ver(itype)
         ! Coordinate vector in cartesian coordinates
         cvec(1)=nncoord_ver(1,inei,itype)+Bas_ver(1,i0)
         cvec(2)=nncoord_ver(2,inei,itype)+Bas_ver(2,i0)
         cvec(3)=nncoord_ver(3,inei,itype)+Bas_ver(3,i0)
         ! Find coordinate vector in basis coordinates
         icvec(1)=cvec(1)*invmatrix(1,1)+cvec(2)*invmatrix(2,1)+cvec(3)*invmatrix(3,1)
         icvec(2)=cvec(1)*invmatrix(1,2)+cvec(2)*invmatrix(2,2)+cvec(3)*invmatrix(3,2)
         icvec(3)=cvec(1)*invmatrix(1,3)+cvec(2)*invmatrix(2,3)+cvec(3)*invmatrix(3,3)
         ! Fold back to original cell
         bsf(1)=floor(icvec(1)+1d-6)
         bsf(2)=floor(icvec(2)+1d-6)
         bsf(3)=floor(icvec(3)+1d-6)
         ! Corresponding position of atom in cell
         rvec(1)=cvec(1)-bsf(1)*C1_ver(1)-bsf(2)*C2_ver(1)-bsf(3)*C3_ver(1)
         rvec(2)=cvec(2)-bsf(1)*C1_ver(2)-bsf(2)*C2_ver(2)-bsf(3)*C3_ver(2)
         rvec(3)=cvec(3)-bsf(1)*C1_ver(3)-bsf(2)*C2_ver(3)-bsf(3)*C3_ver(3)
         ! loop through atoms in cell to find match
         do ia=1,NA_ver
            ! This if counter seems to be the one responsible to see if there is a neighbor, therefore if one changes the way this is done one should be able to get the neighboring spins for each vertex
            if((rvec(1)-Bas_ver(1,ia))**2+(rvec(2)-Bas_ver(2,ia))**2+(rvec(3)-Bas_ver(3,ia))**2<tol) then
               counter=counter+1
               nm_cell(counter,i0)=ia
               if(N1*N2*N3>1) then
                  nm_trunk(1,counter,i0)=bsf(1)
                  nm_trunk(2,counter,i0)=bsf(2)
                  nm_trunk(3,counter,i0)=bsf(3)
               end if
            end if
         end do
      end do
      nnm_cell(I0)=counter
   end do
   !$omp end parallel do
   !
   i_all=-product(shape(nncoord_ver))*kind(nncoord_ver)
   deallocate(nncoord_ver,stat=i_stat)
   call memocc(i_stat,i_all,'nncoord_ver','setup_ver_nm')

   ! Allocate nm_ver_ver : vertex-vertex neighbour map
   allocate(nm_ver_ver(Nvertex,max_no_equiv_ver),stat=i_stat)
   call memocc(i_stat,product(shape(nm_ver_ver))*kind(nm_ver_ver),'nm_ver_ver','setup_nm_ver')
   nm_ver_ver=0

   ! With info in nm_cell and nm_trunk for first NA atoms, make full nm list
   do iz=0, N3-1
      do iy=0, N2-1
         do ix=0, N1-1
            do I0=1, NA_ver
               itype=atype(I0)
               iat=I0+ix*NA_ver+iy*N1*NA_ver+iz*N2*N1*NA_ver
               if (iat/=0) then
                  nmdim_ver_ver(iat)=nmdimt_ver_ver(atype(iat))
                  do inei=1,nnm_cell(i0)
                     ! Designation of cell (This whould be the hopping term that one uses to find neighbours in cells that are not the primitive)
                     if(N1*N2*N3>1) then
                        xc_hop=nm_trunk(1,inei,I0)
                        yc_hop=nm_trunk(2,inei,I0)
                        zc_hop=nm_trunk(3,inei,I0)
                     else
                        xc_hop=0
                        yc_hop=0
                        zc_hop=0
                     end if
                     ! Position in cell
                     j0=nm_cell(inei,I0)
                     ! Wrap around if periodic boundaries
                     jx=xc_hop+ix
                     if(BC1=='P') then
                        jx=mod(jx+1000*N1,N1)
                     else if (N1*N2*N3<=1) then
                        jx=0
                     end if
                     jy=yc_hop+iy
                     if(BC2=='P') then
                        jy=mod(jy+1000*N2,N2)
                     else if (N1*N2*N3<=1) then
                        jy=0
                     end if
                     jz=zc_hop+iz
                     if(BC3=='P') then
                        jz=mod(jz+1000*N3,N3)
                     else if (N1*N2*N3<=1) then
                        jz=0
                     end if
                     ! See if atom exists, then add entry to neighbour map
                     ! This is an important piece for the construction of vertex spin neighbor, for a given vertex one can see the number of neighbors there are and then select them.
                     if(jx>=0.and.jx<N1.and.jy>=0.and.jy<N2.and.jz>=0.and.jz<N3) then
                        jat=j0+jx*NA_ver+jy*N1*NA_ver+jz*N2*N1*NA_ver
                        nm_ver_ver(iat,inei)=jat
                     end if
                  end do
               end if
            end do
         end do
      end do
   end do
   ! Deallocate
   if(N1*N2*N3>1) then
      i_all=-product(shape(nm_trunk))*kind(nm_trunk)
      deallocate(nm_trunk,stat=i_stat)
      call memocc(i_stat,i_all,'nm_trunk','setup_ver_nm')
   end if

   i_all=-product(shape(nm_cell))*kind(nm_cell)
   deallocate(nm_cell,stat=i_stat)
   call memocc(i_stat,i_all,'nm_cell','setup_ver_nm')
   !
   i_all=-product(shape(sym_mats))*kind(sym_mats)
   deallocate(sym_mats,stat=i_stat)
   call memocc(i_stat,i_all,'sym_mats','setup_ver_nm')
   !
   i_all=-product(shape(nnm_cell))*kind(nnm_cell)
   deallocate(nnm_cell,stat=i_stat)
   call memocc(i_stat,i_all,'nnm_cell','setup_ver_nm')
   i_all=-product(shape(nmdimt_ver_ver))*kind(nmdimt_ver_ver)
   deallocate(nmdimt_ver_ver,stat=i_stat)
   call memocc(i_stat,i_all,'nmdimt_ver_ver','setup_ver_nm')
   !
end subroutine setup_ver_nm

!-----------------------------------------------------------------------------
! SUBROUTINE: setup_ver_atom_nm
!> @brief Set up vertex-atom meighbour maps based on the neighbour routines by
!> Anders Bergman
!> @author Jonathan Chico
!-----------------------------------------------------------------------------
subroutine setup_ver_atom_nm(Nvertex,NT_ver,NA_ver,NA_atom,N1,N2,N3,&
   C1_ver,C2_ver,C3_ver,BC1,BC2,BC3,atype,Bas_ver,Bas_atom,&
   max_no_equiv_ver_atom,sym,ver_atom_dist_coord,nm_ver_atom,nmdim_ver_atom)
   !
   !
   implicit none
   !
   ! Quatities from the repetitions of the cell that create the square lattice, they should be the same for the atom and vertex lattice
   !
   integer, intent(in) :: NA_atom  ! Number of atoms in one cell
   integer, intent(in) :: N1  !< Number of cell repetitions in x direction
   integer, intent(in) :: N2  !< Number of cell repetitions in y direction
   integer, intent(in) :: N3  !< Number of cell repetitions in z direction
   real(dblprec), dimension(3,NA_atom), intent(in) :: Bas_atom ! Coordinates of the atoms in the unit cell
   character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
   character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
   character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
   integer, dimension(Nvertex), intent(in) :: atype !< Type of atom
   integer, intent(in) :: sym !< Symmetry of system (0-3)
   !
   ! Quantites that refer to the vertices only
   !
   integer, intent(in) :: Nvertex ! Number of vertices in the system
   integer, intent(in) :: NA_ver ! Number of vertices in the unit cell
   integer, intent(in) :: NT_ver ! Number of vertex types
   real(dblprec), dimension(3), intent(in) :: C1_ver ! First lattice vector of the vertex lattice
   real(dblprec), dimension(3), intent(in) :: C2_ver ! Second lattice vector of the vertex lattice
   real(dblprec), dimension(3), intent(in) :: C3_ver ! Third lattice vector of the vertex lattice
   real(dblprec), dimension(3,NA_ver), intent(in) :: Bas_ver ! Coordinates for basis vertex elements
   real(dblprec), dimension(NT_ver,3), intent(in) :: ver_atom_dist_coord ! Distance between neighbouring atoms and vertex
   ! This are the ouput variables
   integer, intent(out) :: max_no_equiv_ver_atom ! Calculated maximum of neighbours in one shell
   integer, dimension(:,:), allocatable, intent(out) :: nm_ver_atom ! Vertex-vertex neighbour map
   integer, dimension(:), allocatable, intent(out) :: nmdim_ver_atom ! Dimension of the vertex-vertex neighbour map
   !
   integer :: I0
   integer :: i_stat,i_all
   real(dblprec) :: tol
   real(dblprec) :: detmatrix
   real(dblprec), dimension(3,3) :: invmatrix
   !
   integer :: ix,iy,iz,xc_hop,yc_hop,zc_hop,ia,counter,itype,inei
   integer :: iat,jat,j0,jx,jy,jz
   real(dblprec), dimension(3) :: bsf,cvec,icvec,rvec
   real(dblprec), dimension(:,:,:), allocatable:: nncoord_ver_atom ! Full list of atoms neighbours for each vertex
   integer, dimension(:,:), allocatable:: nm_cell
   integer, dimension(:,:,:), allocatable :: nm_trunk
   integer, dimension(:), allocatable :: nnm_cell
   logical :: is_periodic,is_dilute
   !

   tol=0.005_dblprec

   ! Calculate inverse of basis matrix
   detmatrix=C1_ver(1)*C2_ver(2)*C3_ver(3)-C1_ver(1)*C2_ver(3)*C3_ver(2)+&
   C1_ver(2)*C2_ver(3)*C3_ver(1)-C1_ver(2)*C2_ver(1)*C3_ver(3)+&
   C1_ver(3)*C2_ver(1)*C3_ver(2)-C1_ver(3)*C2_ver(2)*C3_ver(1)
   invmatrix=0.0_dblprec
   if(detmatrix/=0.0_dblprec) then
      invmatrix(1,1)=(C2_ver(2)*C3_ver(3)-C3_ver(2)*C2_ver(3))/detmatrix
      invmatrix(1,2)=(C1_ver(3)*C3_ver(2)-C3_ver(3)*C1_ver(2))/detmatrix
      invmatrix(1,3)=(C1_ver(2)*C2_ver(3)-C2_ver(2)*C1_ver(3))/detmatrix
      invmatrix(2,1)=(C2_ver(3)*C3_ver(1)-C3_ver(3)*C2_ver(1))/detmatrix
      invmatrix(2,2)=(C1_ver(1)*C3_ver(3)-C3_ver(1)*C1_ver(3))/detmatrix
      invmatrix(2,3)=(C1_ver(3)*C2_ver(1)-C2_ver(3)*C1_ver(1))/detmatrix
      invmatrix(3,1)=(C2_ver(1)*C3_ver(2)-C3_ver(1)*C2_ver(2))/detmatrix
      invmatrix(3,2)=(C1_ver(2)*C3_ver(1)-C3_ver(2)*C1_ver(1))/detmatrix
      invmatrix(3,3)=(C1_ver(1)*C2_ver(2)-C2_ver(1)*C1_ver(2))/detmatrix
   end if

   ! create all symmetry matrices wrt symmetry type
   call get_symops(sym)

   ! Allocate arrays
   ! Maximum no. of neighbours in each shell is determined by the symmetry (max=48)
   max_no_equiv_ver_atom=nsym

   ! neighbour coordinates
   allocate(nncoord_ver_atom(3,max_no_equiv_ver_atom,NT_ver),stat=i_stat)
   call memocc(i_stat,product(shape(nncoord_ver_atom))*kind(nncoord_ver_atom),'nncoord_ver_atom','setup_ver_atom_nm')

   ! neighbour map information 1: vertex no in first supercell
   allocate(nm_cell(max_no_equiv_ver_atom,NA_ver),stat=i_stat)
   call memocc(i_stat,product(shape(nm_cell))*kind(nm_cell),'nm_cell','setup_ver_atom_nm')

   ! neighbour map information 2: translation in basis vectors
   ! so that full neighbour vector is nm_cell+nm_trunk
   is_periodic=(BC1=='P'.or.BC2=='P'.or.BC3=='P')
   is_dilute=(N1*N2*N3>1.and.is_periodic)
   if(N1*N2*N3>1) then
      allocate(nm_trunk(3,max_no_equiv_ver_atom,NA_ver),stat=i_stat)
      call memocc(i_stat,product(shape(nm_trunk))*kind(nm_trunk),'nm_trunk','setup_ver_atom_nm')
   end if

   ! max no. neighbours in each shell
   allocate(nnm_cell(NA_ver),stat=i_stat)
   call memocc(i_stat,product(shape(nnm_cell))*kind(nnm_cell),'nnm_cell','setup_ver_atom_nm')

   ! Allocate arrays for neighbour map dimensions
   allocate(nmdim_ver_atom(Nvertex),stat=i_stat)
   call memocc(i_stat,product(shape(nmdim_ver_atom))*kind(nmdim_ver_atom),'nmdim_ver_atom','setup_ver_atom_nm')
   allocate(nmdimt_ver_atom(NA_ver),stat=i_stat)
   call memocc(i_stat,product(shape(nmdimt_ver_atom))*kind(nmdimt_ver_atom),'nmdimt_ver_atom','setup_ver_atom_nm')
   nmdim_ver_atom=0
   nmdimt_ver_atom=0

   ! Create full vertex vertex neighbour list according to symmetry
   call get_full_ver_atom_list(NT_ver, ver_atom_dist_coord, max_no_equiv_ver_atom, nncoord_ver_atom)

   ! Start looking in "first cell"
   !$omp parallel do default(shared) private(i0,itype,counter,inei,cvec,icvec,bsf,rvec,ia)
   do I0=1,NA_ver
      itype=atype(I0)
      counter=0
      ! Symmetry equivalent sites in shell
      do inei=1,nmdimt_ver_atom(itype)
         ! Coordinate vector in cartesian coordinates of the neighbour
         cvec(1)=nncoord_ver_atom(1,inei,itype)+Bas_ver(1,I0)
         cvec(2)=nncoord_ver_atom(2,inei,itype)+Bas_ver(2,I0)
         cvec(3)=nncoord_ver_atom(3,inei,itype)+Bas_ver(3,I0)
         ! Find coordinate vector in basis coordinates
         icvec(1)=cvec(1)*invmatrix(1,1)+cvec(2)*invmatrix(2,1)+cvec(3)*invmatrix(3,1)
         icvec(2)=cvec(1)*invmatrix(1,2)+cvec(2)*invmatrix(2,2)+cvec(3)*invmatrix(3,2)
         icvec(3)=cvec(1)*invmatrix(1,3)+cvec(2)*invmatrix(2,3)+cvec(3)*invmatrix(3,3)
         ! Fold back to original cell
         bsf(1)=floor(icvec(1)+1d-6)
         bsf(2)=floor(icvec(2)+1d-6)
         bsf(3)=floor(icvec(3)+1d-6)
         ! Corresponding position of atom in cell
         rvec(1)=cvec(1)-bsf(1)*C1_ver(1)-bsf(2)*C2_ver(1)-bsf(3)*C3_ver(1)
         rvec(2)=cvec(2)-bsf(1)*C1_ver(2)-bsf(2)*C2_ver(2)-bsf(3)*C3_ver(2)
         rvec(3)=cvec(3)-bsf(1)*C1_ver(3)-bsf(2)*C2_ver(3)-bsf(3)*C3_ver(3)
         ! loop through atoms in cell to find match
         do ia=1,NA_atom
            ! If the operation finds that such position as the one created by rvec would correspond to a traslation of the basis positions then it is accepted as a neighbour
            ! If instead of using Bas_ver we use Bas_atom that would mean that we are folding to the atoms in the unit cell, of that way we would obtain the number from the atom
            if((rvec(1)-Bas_atom(1,ia))**2+(rvec(2)-Bas_atom(2,ia))**2+(rvec(3)-Bas_atom(3,ia))**2<tol) then
               counter=counter+1
               nm_cell(counter,I0)=ia
               if(N1*N2*N3>1) then
                  nm_trunk(1,counter,I0)=bsf(1)
                  nm_trunk(2,counter,I0)=bsf(2)
                  nm_trunk(3,counter,I0)=bsf(3)
               end if
            end if
         end do
      end do
      ! Number of neighbours per atom in the unit cell
      nnm_cell(I0)=counter
   end do
   !$omp end parallel do
   !
   i_all=-product(shape(nncoord_ver_atom))*kind(nncoord_ver_atom)
   deallocate(nncoord_ver_atom,stat=i_stat)
   call memocc(i_stat,i_all,'nncoord_ver_atom','setup_ver_atom_nm')

   ! Allocate nm : neighbour map
   allocate(nm_ver_atom(Nvertex,max_no_equiv_ver_atom),stat=i_stat)
   call memocc(i_stat,product(shape(nm_ver_atom))*kind(nm_ver_atom),'nm_ver_atom','setup_nm_ver_atom')
   nm_ver_atom=0

   ! With info in nm_cell and nm_trunk for first NA atoms, make full nm list
   do iz=0, N3-1
      do iy=0, N2-1
         do ix=0, N1-1
            do I0=1, NA_ver
               itype=atype(I0)
               ! This should talk about the vertex numbering as iat runs from 1 to Nvertex
               iat=I0+ix*NA_ver+iy*N1*NA_ver+iz*N2*N1*NA_ver
               if (iat/=0) then
                  nmdim_ver_atom(iat)=nmdimt_ver_atom(atype(iat)) ! The dimension of the neighbour list is the one obtained from the full list for each type
                  do inei=1,nnm_cell(i0) ! Loop over the number of neighbours per atoms in the unit cell
                     ! Designation of cell (This whould be the hopping term that one uses to find neighbours in cells that are not the primitive)
                     if(N1*N2*N3>1) then
                        ! This determies the vector in basis coordinates over which one must fold over to get the neighbour from one neighbour to the unit cell
                        xc_hop=nm_trunk(1,inei,I0)
                        yc_hop=nm_trunk(2,inei,I0)
                        zc_hop=nm_trunk(3,inei,I0)
                     else
                        xc_hop=0
                        yc_hop=0
                        zc_hop=0
                     end if
                     j0=nm_cell(inei,I0) ! This gives the atomic number in the unit cell ia
                     ! Wrap around if periodic boundaries
                     jx=xc_hop+ix ! Translation to new cell plus the number of repetitions done in the atomic cell
                     if(BC1=='P') then
                        jx=mod(jx+1000*N1,N1)
                     else if (N1*N2*N3<=1) then
                        jx=0
                     end if
                     jy=yc_hop+iy
                     if(BC2=='P') then
                        jy=mod(jy+1000*N2,N2)
                     else if (N1*N2*N3<=1) then
                        jy=0
                     end if
                     jz=zc_hop+iz
                     if(BC3=='P') then
                        jz=mod(jz+1000*N3,N3)
                     else if (N1*N2*N3<=1) then
                        jz=0
                     end if
                     ! See if atom exists, then add entry to neighbour map
                     ! This is an important piece for the construction of vertex spin neighbor, for a given vertex one can see the number of neighbors there are and then select them.
                     if(jx>=0.and.jx<N1.and.jy>=0.and.jy<N2.and.jz>=0.and.jz<N3) then
                        jat=j0+jx*NA_atom+jy*N1*NA_atom+jz*N2*N1*NA_atom
                        nm_ver_atom(iat,inei)=jat
                     end if
                  end do
               end if
            end do
         end do
      end do
   end do
   ! Deallocate
   if(N1*N2*N3>1) then
      i_all=-product(shape(nm_trunk))*kind(nm_trunk)
      deallocate(nm_trunk,stat=i_stat)
      call memocc(i_stat,i_all,'nm_trunk','setup_ver_atom_nm')
   end if

   i_all=-product(shape(nm_cell))*kind(nm_cell)
   deallocate(nm_cell,stat=i_stat)
   call memocc(i_stat,i_all,'nm_cell','setup_ver_atom_nm')
   !
   i_all=-product(shape(sym_mats))*kind(sym_mats)
   deallocate(sym_mats,stat=i_stat)
   call memocc(i_stat,i_all,'sym_mats','setup_ver_atom_nm')
   !
   i_all=-product(shape(nnm_cell))*kind(nnm_cell)
   deallocate(nnm_cell,stat=i_stat)
   call memocc(i_stat,i_all,'nnm_cell','setup_ver_atom_nm')
   i_all=-product(shape(nmdimt_ver_atom))*kind(nmdimt_ver_atom)
   deallocate(nmdimt_ver_atom,stat=i_stat)
   call memocc(i_stat,i_all,'nmdimt_ver_atom','setup_ver_atom_nm')
   !
end subroutine setup_ver_atom_nm

!> Find possible symmetry operations depending on assumed symmetry
subroutine get_symops(isym)
   !
   implicit none
   !
   integer, intent(in) :: isym !< Type of assumed symmetry (0-3)
   !
   integer :: i,j,x,y,z,j_s,x1,x2,y1,y2
   integer :: i_stat
   integer :: sym_count
   real(dblprec) :: half,roothalf
   !

   sym_count=0
   if (isym==0) then
      ! No symmetry
      sym_count=1
      allocate(sym_mats(3,3,1),stat=i_stat)
      call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
      sym_mats=0.0_dblprec
      do i=1,3
         sym_mats(i,i,1)=1.0_dblprec
      end do
   else if(isym==1) then
      ! Cubic symmetry
      allocate(sym_mats(3,3,48),stat=i_stat)
      call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
      sym_mats=0.0_dblprec
      sym_count=0
      do i=1,3
         do j=0,1
            j_s=(-1.0_dblprec)**j
            do x=0,1
               do y=0,1
                  do z=0,1
                     sym_count=sym_count+1
                     sym_mats(1,mod(i-j_s,3)+1,sym_count)=(-1.0_dblprec)**x
                     sym_mats(2,mod(i,3)+1,sym_count)=(-1.0_dblprec)**y
                     sym_mats(3,mod(i+j_s,3)+1,sym_count)=(-1.0_dblprec)**z
                  end do
               end do
            end do
         end do
      end do
   else if(isym==2) then
      ! Cubic symmetry in xy-plane
      allocate(sym_mats(3,3,12),stat=i_stat)
      call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
      sym_mats=0.0_dblprec
      sym_count=0
      do j=0,1
         do x=0,1
            do y=0,1
               sym_count=sym_count+1
               sym_mats(1,mod(j,2)+1,sym_count)=(-1.0_dblprec)**x
               sym_mats(2,mod(j+1,2)+1,sym_count)=(-1.0_dblprec)**y
               sym_mats(3,3,sym_count)=1.0_dblprec
            end do
         end do
      end do
   else if(isym==3) then
      ! Hexagonal symmetry
      allocate(sym_mats(3,3,24),stat=i_stat)
      call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
      sym_mats=0.0_dblprec
      sym_count=0
      half=0.5_dblprec
      roothalf=sqrt(3.0_dblprec)*0.5_dblprec
      ! 8 ops due to 'cartesian' inversion
      do x=0,1
         do y=0,1
            do z=0,1
               sym_count=sym_count+1
               sym_mats(1,1,sym_count)=(-1.0_dblprec)**x
               sym_mats(2,2,sym_count)=(-1.0_dblprec)**y
               sym_mats(3,3,sym_count)=(-1.0_dblprec)**z
            end do
         end do
      end do
      ! 16 ops due to 'cartesian' inversion
      do x1=0,1
         do x2=0,1
            do y1=0,1
               do y2=0,1
                  if((-1.0_dblprec)**x1*(-1.0_dblprec)**x2*(-1.0_dblprec)**y1*(-1.0_dblprec)**y2<0.0_dblprec) then
                     do z=0,1
                        sym_count=sym_count+1
                        sym_mats(1,1,sym_count)=(-1.0_dblprec)**x1*half
                        sym_mats(2,1,sym_count)=(-1.0_dblprec)**x2*roothalf
                        sym_mats(1,2,sym_count)=(-1.0_dblprec)**y1*roothalf
                        sym_mats(2,2,sym_count)=(-1.0_dblprec)**y2*half
                        sym_mats(3,3,sym_count)=(-1.0_dblprec)**z
                     end do
                  end if
               end do
            end do
         end do
      end do
   else if(isym==4) then
      open(ifileno,file='sym.mat')
      read(ifileno,*) sym_count
      allocate(sym_mats(3,3,sym_count),stat=i_stat)
      call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
      do j=1,sym_count
         do x=1,3
            read(ifileno,*) (sym_mats(y,x,j),y=1,3)
         end do
      end do
      close(ifileno)
   end if
   nsym=sym_count
   !
end subroutine get_symops

!> Create full neighbour vertex vertex list according to symmetry
subroutine get_full_ver_ver_list(NT_ver, ver_dist_coord, max_no_equiv_ver, nncoord_ver)
   !
   !
   implicit none
   !
   integer, intent(in) :: NT_ver ! Number of types of vertices
   real(dblprec), dimension(NT_ver,3), intent(in) :: ver_dist_coord
   integer, intent(in) :: max_no_equiv_ver ! Calculated maximum of neighbours in one shell
   real(dblprec), dimension(3,max_no_equiv_ver,NT_ver), intent(out) :: nncoord_ver

   real(dblprec) :: tol
   real(dblprec), dimension(3) :: tvect
   integer :: counter
   logical :: unique
   integer :: i, j, k ,itype, isym

   tol=0.0050000_dblprec

   do itype=1,NT_ver
      if (nsym==1) then
         do k=1,3
            ! If the number of symmetry operations is 1 (symmeteries 0) the only neighbour vectors are those explicitly set
            nncoord_ver(k,1,itype)=ver_dist_coord(itype,k)
         end do
         ! Then the number of neghbours by type of atom is set to 1
         nmdimt_ver_ver(itype)=1
      else
         counter=0
         ! Loop over symmetries
         do isym=1,nsym
            tvect=0.0_dblprec
            unique=.true.
            do i=1,3
               do j=1,3
                  ! An auxiliary vector is created by applying symmetry operations to the neighbour vectors
                  tvect(i)=tvect(i)+ver_dist_coord(itype,j)*sym_mats(i,j,isym)
               end do
            end do
            do k=1,counter
               ! If the norm of the difference between the created vector and the ones explicitly set are less than a threshold then the neighbour already exists
               if((tvect(1)-nncoord_ver(1,k,itype))**2+(tvect(2)-nncoord_ver(2,k,itype))**2+(tvect(3)-nncoord_ver(3,k,itype))**2 <tol) unique=.false.
            end do
            if (unique) then
               ! If the neighbour does not exist it is counted and a new vector is set
               counter=counter+1
               do i=1,3
                  nncoord_ver(i,counter,itype)=tvect(i)
               end do
            end if
         end do
         ! The total number of vectors per type is counted
         nmdimt_ver_ver(itype)=counter
      end if
   end do
   !
end subroutine get_full_ver_ver_list

! Create full neighbour vertex atom list according to symmetry
subroutine get_full_ver_atom_list(NT_ver, ver_atom_dist_coord, max_no_equiv_ver_atom, nncoord_ver_atom)
   !
   !
   implicit none
   !
   integer, intent(in) :: NT_ver ! Number of types of vertices
   real(dblprec), dimension(NT_ver,3), intent(in) :: ver_atom_dist_coord
   integer, intent(in) :: max_no_equiv_ver_atom ! Calculated maximum of neighbours in one shell for exchange
   real(dblprec), dimension(3,max_no_equiv_ver_atom,NT_ver), intent(out) :: nncoord_ver_atom


   real(dblprec) :: tol
   real(dblprec), dimension(3) :: tvect
   integer :: counter
   logical :: unique
   integer :: i, j, k ,itype, isym

   tol=0.0050000_dblprec

   do itype=1,NT_ver
      if (nsym==1) then
         do k=1,3
            ! If the number of symmetry operations is 1 (no symmetries) then the nerater neighbour coordinates is only the one specified
            nncoord_ver_atom(k,1,itype)=ver_atom_dist_coord(itype,k)
         end do
         ! If nsym is one the dymension of the neighbour list is one for each pair
         nmdimt_ver_atom(itype)=1
      else
         counter=0
         ! Loop over symmetries
         do isym=1,nsym
            tvect=0.0_dblprec
            unique=.true.
            do i=1,3
               do j=1,3
                  ! Creating new vectors from the symmetry operations
                  tvect(i)=tvect(i)+ver_atom_dist_coord(itype,j)*sym_mats(i,j,isym)
               end do
            end do
            do k=1,counter
               ! If the norm of the distance between the new vector and the position of the neighbour given is less than the treshold the neighbour already exists
               if((tvect(1)-nncoord_ver_atom(1,k,itype))**2+(tvect(2)-nncoord_ver_atom(2,k,itype))**2+(tvect(3)-nncoord_ver_atom(3,k,itype))**2 <tol) unique=.false.
            end do
            ! If the neighbour does not exists then it is created and counted (allowing to know how many neighbours per atom type there are)
            if (unique) then
               counter=counter+1
               do i=1,3
                  nncoord_ver_atom(i,counter,itype)=tvect(i)
               end do
            end if
         end do
         nmdimt_ver_atom(itype)=counter
      end if
   end do
   !
end subroutine get_full_ver_atom_list

! Mounts the actual Neighbour list
subroutine setup_neighbour_list(Nvertex, max_no_equiv, nlistsize, nlist, nm, nmdim)
   !
   implicit none
   !
   integer, intent(in) :: Nvertex ! Number of vertices in system
   integer, dimension(Nvertex), intent(out):: nlistsize ! Size of neighbour list for vertex vertex or vertex atoms
   integer, intent(in) :: max_no_equiv ! Calculated maximum of neighbours in one shell for vertex vertex of vertex atoms
   integer, dimension(max_no_equiv,Nvertex), intent(out) :: nlist ! Neighbour list for vertex vertex or vertex atoms
   integer, dimension(Nvertex,max_no_equiv), intent(in) :: nm ! Vertex-Vertex or Vertex-atom neighbour map
   integer, dimension(Nvertex), intent(in) :: nmdim ! Dimension of the vertex-vertex (vertex-atom) neighbour map
   !
   integer :: i, j, l
   integer :: iver, jver
   integer :: ncount
   logical :: exis

   iver=0

   !$omp parallel do default(shared) private(i,ncount,j,exis,l,iver,jver)
   ! Loop over atoms
   do i=1, Nvertex
      ncount = 1
      ! Sites
      do j=1, nmdim(i)
         ! Existing neighbours
         if (nm(i,j)>0) then
            exis=.false.
            do l=1,ncount-1
               if(nlist(l,i)==nm(i,j)) exis=.true.
            end do
            nlist(ncount,i) = nm(i,j)
            ncount = ncount + 1
         end if
      end do
      nlistsize(i) = ncount-1
   end do
   !omp end parallel do

end subroutine setup_neighbour_list

logical function ice_rule(Natom,Mensemble, en,emom,trial_vert,total_trial,ice_count,ice_count_up, ice_count_down)

   implicit none

   integer, intent(in) :: Natom
   integer, intent(in) :: Mensemble
   integer, intent(in) :: en ! Current ensemble
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
   integer, intent(in) :: trial_vert
   real(dblprec), intent(out) :: total_trial

   real(dblprec), intent(out) :: ice_count_up, ice_count_down
   real(dblprec), intent(out) :: ice_count

   integer :: max_ver_size
   real(dblprec) :: trial_ene
   real(dblprec) :: temp_trial
   real(dblprec), dimension(3) :: trial_ene_vec

   integer :: atom

   ! Calculate the total moment of a given vertex to check if the ice rule is fullfilled
   ice_count=0.0_dblprec
   ice_count_down=0.0_dblprec
   ice_count_up=0.0_dblprec
   trial_ene_vec=0.0_dblprec
   trial_ene=0.0_dblprec
   temp_trial=0.0_dblprec
   total_trial=0.0_dblprec

   max_ver_size=maxval(nlistsize_ver_atom)
   do atom=1, nlistsize_ver_atom(trial_vert)
      trial_ene=0.0_dblprec
      temp_trial=0.0_dblprec
      trial_ene_vec=0.0_dblprec
      trial_ene_vec(1) = emom(1,nlist_ver_atom(atom,trial_vert),en)*vert_ice_coord(1,trial_vert,atom)
      trial_ene_vec(2) = emom(2,nlist_ver_atom(atom,trial_vert),en)*vert_ice_coord(2,trial_vert,atom)
      trial_ene_vec(3) = emom(3,nlist_ver_atom(atom,trial_vert),en)*vert_ice_coord(3,trial_vert,atom)
      trial_ene=trial_ene_vec(1)+trial_ene_vec(2)+trial_ene_vec(3)

      if (trial_ene.eq.0.0_dblprec) then
         write(*,'(a)') 'WARNING'
         write(*,'(2i8,4f16.8)') trial_vert, nlist_ver_atom(atom,trial_vert),emom(1:2,nlist_ver_atom(atom,trial_vert),en),vert_ice_coord(1:2,trial_vert,atom)
         stop
      elseif (trial_ene.ne.0.0_dblprec) then
         temp_trial = trial_ene
      endif
      if (temp_trial.gt.0.0_dblprec) then
         ice_count_up = ice_count_up+1.0_dblprec
      else if (temp_trial.lt.0.0_dblprec) then
         ice_count_down = ice_count_down-1.0_dblprec
      end if
      total_trial = total_trial+temp_trial

   end do

   ice_count=ice_count_up+ice_count_down
   ! The first condition is to make sure that the edges are not set as fullfiling the ice rule

   if ( (abs(ice_count_up)+abs(ice_count_down)).lt.max_ver_size ) then
      ice_rule=.false.
      ice_count=(abs(ice_count_up)+abs(ice_count_down))
      ! Once the edges are taken out of the picture if one sums up the counters one can get one of the configurations which fullfills the ice rule
   else if ( abs(ice_count)==0.0_dblprec ) then
      ice_rule=.true.
   else
      ice_rule=.false.
   end if

   return
end function ice_rule

!> Calculate total energy for a loop in the spin Ice
subroutine calculate_spinice_energy(Natom, Mensemble, nHam, max_no_neigh, conf_num, ncoup, nlist, nlistsize, aham, &
   do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, &
   do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
   do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, &
   do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
   taniso, eaniso, kaniso,sb,emomM, emom, &
   Ncount_atom, loop_len, extfield, de, k, &
   do_anisotropy,do_dip, Qdip)

   use Constants, only : mub

   !.. Implicit declarations
   implicit none

   integer, intent(in) :: Natom !< Number of atoms in system
   integer, intent(in) :: Mensemble !< Number of ensembles
   integer, intent(in) :: nHam  !< Number of atoms in Hamiltonian
   integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
   integer, intent(in) :: conf_num   !< Number of configurations for LSF
   real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
   integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
   integer, dimension(nHam),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
   integer, dimension(Natom),intent(in) :: aham !< Hamiltonian look-up table
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
   real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
   real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
   real(dblprec), dimension(Natom), intent(in) :: sb!< Ratio between the Anisotropy constants
   real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
   integer, intent(in) :: loop_len !< Lenght of cluster loop
   integer, dimension(Natom), intent(in) :: Ncount_atom !< List of atoms in cluster loop
   real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
   real(dblprec), intent(out):: de  !< Energy difference
   integer, intent(in) :: k !< Current ensemble
   real(dblprec) :: aw1,aw2
   integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
   integer, intent(in) :: do_anisotropy
   real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction

   !.. Local scalars
   integer :: j,iloop,cloop, catom, catom_h
   real(dblprec) :: tta,ttb,e_c,e_t
   real(dblprec) :: bqmdot

   !.. Local arrays

   real(dblprec), dimension(2) :: tt_loop
   real(dblprec), dimension(:,:), allocatable :: befftemp_loop

   !.. Executable statements

   ! Energy is calculated as m(i)*Beff(i)
   allocate(befftemp_loop(3,loop_len))
   tt_loop=0.0_dblprec
   ! Loop over flipped and non-flipped cluster loop.
   !
   e_c=0.0_dblprec
   e_t=0.0_dblprec

   do iloop=1,2
      !
      if(iloop==2) then
         do cloop=1,loop_len
            catom=Ncount_atom(cloop)
            emomM(:,catom,k)=-1.0_dblprec*emomM(:,catom,k)
         end do
      end if
      befftemp_loop=0.0_dblprec

      do cloop=1,loop_len
         catom=Ncount_atom(cloop)
         catom_h=aham(catom)
         do j=1,nlistsize(catom_h)
            if (iloop==1) then
               e_c=e_c-ncoup(j,catom_h,1)*(emomM(1,catom,k)*emomM(1,nlist(j,catom),k)+ &
               emomM(2,catom,k)*emomM(2,nlist(j,catom),k)+ &
               emomM(3,catom,k)*emomM(3,nlist(j,catom),k))
            else
               e_t=e_t-ncoup(j,catom_h,1)*(emomM(1,catom,k)*emomM(1,nlist(j,catom),k)+ &
               emomM(2,catom,k)*emomM(2,nlist(j,catom),k)+ &
               emomM(3,catom,k)*emomM(3,nlist(j,catom),k))
            endif
         end do

         ! Anisotropy
         if (do_anisotropy==1) then
            ! Uniaxial anisotropy
            if (taniso(catom)==1) then
               tta=(emomM(1,catom,k)*eaniso(1,catom)+emomM(2,catom,k)*eaniso(2,catom)+emomM(3,catom,k)*eaniso(3,catom))
               ttb=(emomM(1,catom,k)*eaniso(1,catom)+emomM(2,catom,k)*eaniso(2,catom)+emomM(3,catom,k)*eaniso(3,catom))
               if (iloop==1) then
                  e_c=e_c+kaniso(1,catom)*(tta**2)+kaniso(2,catom)*(tta**4)
               else
                  e_t=e_t+kaniso(1,catom)*(ttb**2)+kaniso(2,catom)*(ttb**4)
               endif

            ! Cubic anisotropy
            elseif (taniso(catom)==2) then
               if (iloop==1) then
                  e_c=e_c+kaniso(1,catom)*(emomM(1,catom,k)**2*emomM(2,catom,k)**2+ &
                  emomM(2,catom,k)**2*emomM(3,catom,k)**2+emomM(3,catom,k)**2*emomM(1,catom,k)**2)+ &
                  kaniso(2,catom)*(emomM(1,catom,k)**2*emomM(2,catom,k)**2*emomM(3,catom,k)**2)
               else
                  e_t=e_t+kaniso(1,catom)*(emomM(1,catom,k)**2*emomM(2,catom,k)**2+ &
                  emomM(2,catom,k)**2*emomM(3,catom,k)**2+emomM(3,catom,k)**2*emomM(1,catom,k)**2)+ &
                  kaniso(2,catom)*(emomM(1,catom,k)**2*emomM(2,catom,k)**2*emomM(3,catom,k)**2)
               endif
            endif
            ! When both Cubic and Uniaxial are switched on
            if (taniso(catom)==7) then

               ! Uniaxial anisotropy
               tta=(emomM(1,catom,k)*eaniso(1,catom)+emomM(2,catom,k)*eaniso(2,catom)+emomM(3,catom,k)*eaniso(3,catom))

               if (iloop==1) then
                  e_c=e_c+kaniso(1,catom)*(tta**2)+kaniso(2,catom)*(tta**4)
               else
                  e_t=e_t+kaniso(1,catom)*(tta**2)+kaniso(2,catom)*(tta**4)
               endif

               ! Cubic anisotropy
               aw1=kaniso(1,catom)*sb(catom)
               aw2=kaniso(2,catom)*sb(catom)

               if (iloop==1) then
                  e_c=e_c+aw1*(emomM(1,catom,k)**2*emomM(2,catom,k)**2+ &
                  emomM(2,catom,k)**2*emomM(3,catom,k)**2+emomM(3,catom,k)**2*emomM(1,catom,k)**2)+ &
                  aw2*(emomM(1,catom,k)**2*emomM(2,catom,k)**2*emomM(3,catom,k)**2)
               else
                  e_t=e_t+aw1*(emomM(1,catom,k)**2*emomM(2,catom,k)**2+ &
                  emomM(2,catom,k)**2*emomM(3,catom,k)**2+emomM(3,catom,k)**2*emomM(1,catom,k)**2)+ &
                  aw2*(emomM(1,catom,k)**2*emomM(2,catom,k)**2*emomM(3,catom,k)**2)
               endif
            endif
         endif
         ! DM interaction
         if (do_dm==1) then
            do j=1,dmlistsize(catom)
               if (iloop==1) then
                  e_c=e_c+dm_vect(1,j,catom)*(emomM(2,catom,k)*emomM(3,dmlist(j,catom),k)- &
                  emom(3,catom,k)*emomM(2,dmlist(j,catom),k))+ &
                  dm_vect(2,j,catom)*(emomM(3,catom,k)*emomM(1,dmlist(j,catom),k)- &
                  emomM(1,catom,k)*emomM(3,dmlist(j,catom),k))+ &
                  dm_vect(3,j,catom)*(emom(1,catom,k)*emomM(2,dmlist(j,catom),k)- &
                  emomM(2,catom,k)*emomM(1,dmlist(j,catom),k))
               else
                  e_t=e_t+dm_vect(1,j,catom)*(emomM(2,catom,k)*emomM(3,dmlist(j,catom),k)- &
                  emomM(3,catom,k)*emomM(2,dmlist(j,catom),k))+ &
                  dm_vect(2,j,catom)*(emomM(3,catom,k)*emomM(1,dmlist(j,catom),k)- &
                  emomM(1,catom,k)*emomM(3,dmlist(j,catom),k))+ &
                  dm_vect(3,j,catom)*(emomM(1,catom,k)*emomM(2,dmlist(j,catom),k)- &
                  emomM(2,catom,k)*emomM(1,dmlist(j,catom),k))
               endif
            end do
         end if
         ! PD interaction
         if(do_pd==1) then
            do j=1,pdlistsize(catom)
               if (iloop==1) then
                  e_c=e_c-pd_vect(1,j,catom)*emomM(1,catom,k)*emomM(1,pdlist(j,catom),k)- &
                  pd_vect(4,j,catom)*emomM(1,catom,k)*emomM(2,pdlist(j,catom),k)- &
                  pd_vect(5,j,catom)*emomM(1,catom,k)*emomM(3,pdlist(j,catom),k)- &
                  pd_vect(4,j,catom)*emomM(2,catom,k)*emomM(1,pdlist(j,catom),k)- &
                  pd_vect(2,j,catom)*emomM(2,catom,k)*emomM(2,pdlist(j,catom),k)- &
                  pd_vect(6,j,catom)*emomM(2,catom,k)*emomM(3,pdlist(j,catom),k)- &
                  pd_vect(5,j,catom)*emomM(3,catom,k)*emomM(1,pdlist(j,catom),k)- &
                  pd_vect(6,j,catom)*emomM(3,catom,k)*emomM(2,pdlist(j,catom),k)- &
                  pd_vect(3,j,catom)*emomM(3,catom,k)*emomM(3,pdlist(j,catom),k)
               else
                  e_t=e_t-pd_vect(1,j,catom)*emomM(1,catom,k)*emomM(1,pdlist(j,catom),k)- &
                  pd_vect(4,j,catom)*emomM(1,catom,k)*emomM(2,pdlist(j,catom),k)- &
                  pd_vect(5,j,catom)*emomM(1,catom,k)*emomM(3,pdlist(j,catom),k)- &
                  pd_vect(4,j,catom)*emomM(2,catom,k)*emomM(1,pdlist(j,catom),k)- &
                  pd_vect(2,j,catom)*emomM(2,catom,k)*emomM(2,pdlist(j,catom),k)- &
                  pd_vect(6,j,catom)*emomM(2,catom,k)*emomM(3,pdlist(j,catom),k)- &
                  pd_vect(5,j,catom)*emomM(3,catom,k)*emomM(1,pdlist(j,catom),k)- &
                  pd_vect(6,j,catom)*emomM(3,catom,k)*emomM(2,pdlist(j,catom),k)- &
                  pd_vect(3,j,catom)*emomM(3,catom,k)*emomM(3,pdlist(j,catom),k)
               endif
            end do
         end if
         ! BIQDM interaction
         if(do_biqdm==1) then
            do j=1,biqdmlistsize(catom)
               if (iloop==1) then
                  e_c=e_c-biqdm_vect(1,j,catom)*(emomM(2,catom,k)*emomM(3,biqdmlist(j,catom),k)- &
                  emomM(3,catom,k)*emomM(2,biqdmlist(j,catom),k))**2- &
                  biqdm_vect(1,j,catom)*(emom(3,catom,k)*emomM(1,biqdmlist(j,catom),k)- &
                  emomM(1,catom,k)*emomM(3,biqdmlist(j,catom),k))**2- &
                  biqdm_vect(1,j,catom)*(emom(1,catom,k)*emomM(2,biqdmlist(j,catom),k)- &
                  emomM(2,catom,k)*emomM(1,biqdmlist(j,catom),k))**2
               else
                  e_t=e_t-biqdm_vect(1,j,catom)*(emomM(2,catom,k)*emomM(3,biqdmlist(j,catom),k)- &
                  emomM(3,catom,k)*emomM(2,biqdmlist(j,catom),k))**2- &
                  biqdm_vect(1,j,catom)*(emomM(3,catom,k)*emomM(1,biqdmlist(j,catom),k)- &
                  emomM(1,catom,k)*emomM(3,biqdmlist(j,catom),k))**2- &
                  biqdm_vect(1,j,catom)*(emomM(1,catom,k)*emomM(2,biqdmlist(j,catom),k)- &
                  emomM(2,catom,k)*emomM(1,biqdmlist(j,catom),k))**2
               endif
            end do
         end if
         ! BQ interaction
         if(do_bq==1) then
            do j=1,bqlistsize(catom)
               bqmdot=emomM(1,bqlist(j,catom),k)*emomM(1,catom,k)+&
               emomM(2,bqlist(j,catom),k)*emomM(2,catom,k)+&
               emomM(3,bqlist(j,catom),k)*emomM(3,catom,k)
               if (iloop==1) then
                  e_c=e_c-j_bq(j,catom)*bqmdot**2
               else
                  e_t=e_t-j_bq(j,catom)*bqmdot**2
               endif
            end do
         end if

         ! Dipolar Interaction Jonathan 19-07-2012
         if(present(Qdip)) then
            if(do_dip==1) then
               do j=1,Natom
                  ! Current spin
                  if (iloop==1) then
                     e_c = e_c-( Qdip(1,1,j,catom)*emomM(1,j,k)*emomM(1,catom,k) + Qdip(2,1,j,catom)*emomM(2,j,k)*emomM(1,catom,k)+ &
                     Qdip(3,1,j,catom)*emomM(3,j,k)*emomM(1,catom,k) + Qdip(1,2,j,catom)*emomM(1,j,k)*emomM(2,catom,k)+ &
                     Qdip(2,2,j,catom)*emomM(2,j,k)*emomM(2,catom,k) + Qdip(3,2,j,catom)*emomM(3,j,k)*emomM(2,catom,k)+ &
                     Qdip(1,3,j,catom)*emomM(1,j,k)*emomM(3,catom,k) + Qdip(2,3,j,catom)*emomM(2,j,k)*emomM(3,catom,k)+ &
                     Qdip(3,3,j,catom)*emomM(3,j,k)*emomM(3,catom,k))
                  else
                     ! Trial spin
                     e_t = e_t-( Qdip(1,1,j,catom)*emomM(1,j,k)*emomM(1,catom,k) + Qdip(2,1,j,catom)*emomM(2,j,k)*emomM(1,catom,k)+ &
                     Qdip(3,1,j,catom)*emomM(3,j,k)*emomM(1,catom,k) + Qdip(1,2,j,catom)*emomM(1,j,k)*emomM(2,catom,k)+ &
                     Qdip(2,2,j,catom)*emomM(2,j,k)*emomM(2,catom,k) + Qdip(3,2,j,catom)*emomM(3,j,k)*emomM(2,catom,k)+ &
                     Qdip(1,3,j,catom)*emomM(1,j,k)*emomM(3,catom,k) + Qdip(2,3,j,catom)*emomM(2,j,k)*emomM(3,catom,k)+ &
                     Qdip(3,3,j,catom)*emomM(3,j,k)*emomM(3,catom,k))
                  endif
               end do
            end if
         end if

         ! Add Zeeman term
         if (iloop==1) then
            e_c=e_c-extfield(1)*emomM(1,catom,k)-extfield(2)*emomM(2,catom,k)-extfield(3)*emomM(3,catom,k)
         else
            e_t=e_t-extfield(1)*emomM(1,catom,k)-extfield(2)*emomM(2,catom,k)-extfield(3)*emomM(3,catom,k)
         endif

         ! Add static field
         if (iloop==1) then
            tt_loop(iloop)=tt_loop(iloop)+e_c
         else
            tt_loop(iloop)=tt_loop(iloop)+e_t
         endif
      end do
   end do

   ! Energy difference = mu_b*[Eflip-Epre]
   de=-mub*(tt_loop(1)-tt_loop(2))

   deallocate(befftemp_loop)

   return

end subroutine calculate_spinice_energy


!> Flip spin in spin-ice
subroutine Spin_Ice_flip(emom, emomM, mmom, Ncount_atom, Natom, Mensemble,Temperature, de, k,loop_len, flipped)

   use Constants
   use RandomNumbers, only : rng_uniform

   !.. Implicit declarations
   implicit none

   integer, intent(in) :: Natom !< Number of atoms in system
   integer, intent(in) :: loop_len
   integer, intent(in) :: Mensemble !< Number of ensembles
   integer, dimension(Natom), intent(in) :: Ncount_atom
   real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
   real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
   real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
   real(dblprec), intent(in):: de  !< Energy difference
   real(dblprec), intent(in):: Temperature !< Temperature

   integer, intent(in) :: k !< Current ensemble
   logical, intent(out) :: flipped
   integer :: i
   real(dblprec) :: beta , flipprob(1)

   call rng_uniform(flipprob,1)
   beta=1_dblprec/k_bolt/Temperature

   if(de<=0.0_dblprec .or. flipprob(1)<=exp(-beta*de)) then
      do i=1, loop_len
         if (Ncount_atom(i).ne.0) then
            emom(:,Ncount_atom(i),k)=-1.0_dblprec*emom(:,Ncount_atom(i),k)
            !emom(:,i,k)=-1.0_dblprec*emom(:,i,k)
            emomM(:,Ncount_atom(i),k)=mmom(Ncount_atom(i),k)*emom(:,Ncount_atom(i),k)
         endif
         !emomM(:,i,k)=mmom(i,k)*emom(:,i,k)
      enddo
      mchits_spin_ice=mchits_spin_ice+1
      flipped=.true.
   else
      flipped=.false.
   endif

end subroutine Spin_Ice_flip

end module SpinIce
