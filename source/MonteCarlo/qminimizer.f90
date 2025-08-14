!-------------------------------------------------------------------------------
! MODULE: qminimizer
!> @brief
!> Spin-spiral energy minimization
!> @author
!> Anders Bergman
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module qminimizer
   use Parameters
   use Profiling
   use FieldData
   use HamiltonianActions
   !
   !
   implicit none
   !
   real(dblprec) :: fac_2d
   !
   !
   ! Input parameters
   !
   integer, parameter :: nq=1
   !
   !
   ! Data for minimization and ss-vectors
   real(dblprec), dimension(3,3) :: q  !< q-vector for spin spiral
   real(dblprec), dimension(3,3) :: s  !< 
   real(dblprec), dimension(3) :: n_vec !< unit-vector perpendicular to spins
   real(dblprec), dimension(3) :: theta, phi
   real(dblprec) :: theta_glob
   !
   character*1 :: qm_cellrot   !< Rotate each cell instead of each atom 
   character*1 :: qm_rot   !< Rotate magnetic texture instead of pure spin spirals
   character*1 :: qm_oaxis !< Ensure that the rotational axis is perpendicular to m
   character*1 :: qm_type  !< Which type of qm_axis to use (C)ycloidal, (H)elical, or (G)eneral
   integer :: qm_no_excluded !< How many types of atoms to exclude from minimization (default 0)
   character*1 :: qm_relax   !< Perform relaxation of the magnetic spiral textures (Y/N)
   character*2 :: qm_relax_mode   !< Select relaxation mode (M)etropolis or (H)eat-bath
   integer :: qm_relax_steps   !< Number of relaxation steps
   real(dblprec) :: qm_relax_temp   !< Relaxation temperature

   !
   ! Control parameters for line sweep
   real(dblpreC) :: q_min
   real(dblpreC) :: q_max
   integer :: nstep
   integer, dimension(:), allocatable :: qm_excluded_types
   logical, dimension(:), allocatable :: qm_excluded_atoms
   !
   private
   !
   !public :: mini_q, plot_q, qmc, sweep_q2, sweep_q3, plot_q3
   !public :: plot_q, qmc, sweep_q2, sweep_q3, plot_q3
   !public :: sweep_cube, plot_cube
   public :: read_parameters_qminimizer,qminimizer_init, qminimizer_wrapper
   !public :: sweep_q,read_parameters_qminimizer,qminimizer_init
   !
contains

   subroutine qminimizer_wrapper(qmode)
      !
      use InputData!, only : Natom, Mensemble, NA
      use optimizationRoutines, only : OPT_flag,max_no_constellations,maxNoConstl,unitCellType,constlNCoup, constellations,constellationsNeighType
      use macrocells, only : Num_macro,cell_index, emomM_macro,macro_nlistsize
      use MomentData, only : emom, emomM, mmom
      use SystemData, only : coord
      use Qvectors,   only : q, nq

      implicit none
      !
      character*1, intent(in) :: qmode

      !

      if (qmode=='Q') then
         ! Spin spiral minimization initial phase
         !call mini_q(Natom,Mensemble,NA,coord,do_jtensor,exc_inter,do_dm,do_pd,          &
         call sweep_q2(Natom,Mensemble,NA,coord,emomM,mmom,iphfield,    &
            OPT_flag,max_no_constellations,maxNoConstl,unitCellType,constlNCoup,    &
            constellations,constellationsNeighType,Num_macro,cell_index,  &
            emomM_macro,macro_nlistsize,simid,q,nq)
         call plot_q(Natom,Mensemble,coord,emom,emomM,mmom,simid)
      elseif (qmode=='Z') then
         ! Spin spiral minimization initial phase
         !call mini_q(Natom,Mensemble,NA,coord,do_jtensor,exc_inter,do_dm,do_pd,          &
         !call sweep_q3(Natom,Mensemble,NA,coord,do_jtensor,exc_inter,do_dm,do_pd,    &
         call sweep_cube(Natom,Mensemble,NA,coord,emomM,mmom,iphfield,    &
            OPT_flag,max_no_constellations,maxNoConstl,unitCellType,constlNCoup,    &
            constellations,constellationsNeighType,Num_macro,cell_index,  &
            emomM_macro,macro_nlistsize,simid,q,nq)
         call plot_cube(Natom,Mensemble,coord,emom,emomM,mmom,simid)
      elseif (qmode=='Y') then
         ! Spin spiral minimization initial phase
         call sweep_q3(Natom,Mensemble,NA,coord,emomM,mmom,iphfield,    &
            OPT_flag,max_no_constellations,maxNoConstl,unitCellType,constlNCoup,    &
            constellations,constellationsNeighType,Num_macro,cell_index,  &
            emomM_macro,macro_nlistsize,simid,q,nq)
         call plot_q3(Natom,Mensemble,coord,emom,emomM,mmom,simid)
      elseif (mode=='S') then
         ! Spin spiral minimization measurement phase
         call qmc(Natom,Mensemble,NA,N1,N2,N3,coord, emomM,mmom,hfield,&
            OPT_flag,max_no_constellations,maxNoConstl,unitCellType,constlNCoup,    &
            constellations,constellationsNeighType,Num_macro,cell_index,  &
            emomM_macro,macro_nlistsize)
         call plot_q(Natom, Mensemble, coord, emom, emomM, mmom,simid)

      end if
   
   end subroutine qminimizer_wrapper

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: sweep_q2
   !> @brief Stupid line search minimization of spin spirals (clone of sweep_q
   !  but for external q-point set.
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine sweep_q2(Natom,Mensemble,NA,coord,emomM,mmom,hfield,OPT_flag,           &
      max_no_constellations,maxnoconstl,unitcelltype,constlncoup,constellations,    &
      constellationsneightype,num_macro,cell_index,emomm_macro,           &
      macro_nlistsize,simid,qpts,nq)

      use randomnumbers, only: rng_uniform,rng_gaussian
      use inputdata, only : n1,n2,n3
      use math_functions, only : f_wrap_coord_diff, f_cross_product
      use depondt, only : rodmat
      use diamag, only : diamag_qvect
      use momentdata, only : emom
      use mc_driver, only : mc_minimal
      use montecarlo
      !
      !.. implicit declarations
      implicit none

      integer, intent(in) :: natom !< number of atoms in system
      integer, intent(in) :: mensemble !< number of ensembles
      integer, intent(in) :: na  !< number of atoms in one cell
      real(dblprec), dimension(3,natom), intent(in) :: coord !< coordinates of atoms
      real(dblprec), dimension(3,natom,mensemble), intent(inout) :: emomm  !< current magnetic moment vector
      real(dblprec), dimension(natom,mensemble), intent(inout) :: mmom !< magnitude of magnetic moments
      real(dblprec), dimension(3), intent(in) :: hfield !< constant effective field
      !! +++ new variables due to optimization routines +++ !!
      integer, intent(in) :: max_no_constellations ! the maximum (global) length of the constellation matrix
      ! number of entries (number of unit cells*number of atoms per unit cell) in the constellation matrix per ensemble
      integer, dimension(mensemble), intent(in) :: maxnoconstl
      ! see optimizationroutines.f90 for details on classification
      integer, dimension(natom, mensemble), intent(in) :: unitcelltype ! array of constellation id and classification (core, boundary, or noise) per atom
      ! matrix relating the interatomic exchanges for each atom in the constellation matrix
      real(dblprec), dimension(ham%max_no_neigh, max_no_constellations,mensemble), intent(in) :: constlncoup
      ! matrix storing all unit cells belonging to any constellation
      real(dblprec), dimension(3,max_no_constellations, mensemble), intent(in) :: constellations
      ! optimization flag (1 = optimization on; 0 = optimization off)
      logical, intent(in) :: opt_flag
      ! matrix storing the type of the neighbours within a given neighbourhood of a constellation; default is 1 outside the neighbourhood region
      ! the default is to achieve correct indexing. note here also that constlncoup will result in a net zero contribution to the heissenberg exchange term
      integer, dimension(ham%max_no_neigh,max_no_constellations,mensemble), intent(in) :: constellationsneightype
      ! internal effective field arising from the optimization of the heissenberg exchange term
      integer, intent(in) :: num_macro !< number of macrocells in the system
      integer, dimension(natom), intent(in) :: cell_index !< macrocell index for each atom
      integer, dimension(num_macro), intent(in) :: macro_nlistsize !< number of atoms per macrocell
      real(dblprec), dimension(3,num_macro,mensemble), intent(in) :: emomm_macro !< the full vector of the macrocell magnetic moment
      character(len=8), intent(in) :: simid  !< name of simulation
      integer, intent(in) :: nq  !< number of qpoints
      real(dblprec), dimension(3,nq), intent(in) :: qpts !< array of q-points
      !
      integer :: iq
      !
      real(dblprec), dimension(3) :: m_j
      real(dblprec) :: pi, qr
      integer :: i, k, ia, lhit, nhits, countstart
      real(dblprec) :: energy, min_energy
      character(len=30) :: filn
      !
      real(dblprec), dimension(3,3) :: q_best
      real(dblprec), dimension(3,3) :: s_save
      real(dblprec), dimension(3) :: theta_best
      real(dblprec) :: theta_glob_best
      real(dblprec), dimension(3) :: phi_best
      integer :: iter, iscale, i1 ,i2 ,i3
      real(dblprec), dimension(3) :: srvec 
      real(dblprec), dimension(3,3) :: r_mat
      real(dblprec), dimension(:,:,:), allocatable :: emomm_start

      real(dblprec), dimension(3) :: mavg
      !
      integer :: i_stat, i_all
      !
      integer :: ia_cell
      !
      pi=4._dblprec*atan(1._dblprec)
      theta_glob_best=0.0_dblprec
      theta_best=0.0_dblprec
      phi_best=0.0_dblprec
      phi=0.0_dblprec
      !
      ! normal vector
      ! read from file or default
      !n_vec(1)=0.0_dblprec;n_vec(2)=0.0_dblprec;n_vec(3)=1.0_dblprec;
      !
      ! starting atom
      i1 = n1/2
      i2 = n2/2
      i3 = n3/2

      countstart = 0+i1*na+i2*n1*na+i3*n2*n1*na
      ! build lookup table for excluding atoms in minimization steps
      if (.not. allocated(qm_excluded_atoms)) call qminimizer_build_exclude_list(natom)
      if (qm_relax=='Y') call allocate_mcdata(natom,1)
      !
      write(filn,'(''qm_sweep.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position="append")
      write(ofileno,'(a)') "#    iter                          q-vector                                 energy(mev)  "

      write(filn,'(''qm_minima.'',a,''.out'')') trim(simid)
      open(ofileno2,file=filn, position="append")
      write(ofileno2,'(a)') "#    iter                          q-vector                                 energy(mry)  "
     
    
      if (qm_rot=='y'.or.qm_cellrot=='y') then
         !print '(3f12.5)', emomm
         allocate(emomm_start(3,natom,mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(emomm_start))*kind(emomm_start),'emomm_start','sweep_q2')
         emomm_start=emomm

         if (qm_oaxis=='y') then
            mavg(1)=sum(emomm(1,:,:))/natom/mensemble
            mavg(2)=sum(emomm(2,:,:))/natom/mensemble
            mavg(3)=sum(emomm(3,:,:))/natom/mensemble
            mavg=mavg/(norm2(mavg)+1.0e-12_dblprec)
            if(abs(mavg(3)-1.0_dblprec)<1.0e-6_dblprec) then
               n_vec(1)=1.0_dblprec;n_vec(2)=0.0_dblprec;n_vec(3)=0.0_dblprec
            else
               n_vec(1)=mavg(1); n_vec(2)=-mavg(2); n_vec(3)=0.0_dblprec
            end if
            n_vec=n_vec/norm2(n_vec)
         end if

      end if
      
      !
      theta = 0.0_dblprec
      theta_glob = 0.0_dblprec
      phi = 0.0_dblprec
      min_energy=1.0d4

      !
      lhit=0
      nhits=0
      iscale=1

      ! switch rotation direction

      ! calculate total external field (not included yet)
      do k=1,mensemble
         do i=1,natom
            external_field(1:3,i,k)= hfield
            beff(1:3,i,k)=0.0_dblprec
            beff1(1:3,i,k)=0.0_dblprec
            beff2(1:3,i,k)=0.0_dblprec
         end do
      end do

      ! only use first ensemble
      k=1
      do iq=1,nq
         iter=iq
         !         !
         ! set up spin-spiral magnetization (only first cell+ neighbours)
         energy=0.0_dblprec
         ! check if diamag_qvect and then only use that spin spiral vector
         if (norm2(diamag_qvect)>0.0_dblprec) then
            q(:,1)=diamag_qvect
         else
            q(:,1)=qpts(:,iq)
         end if
         !!! s(1,1)=q(2,1)*n_vec(3)-q(3,1)*n_vec(2)
         !!! s(2,1)=q(3,1)*n_vec(1)-q(1,1)*n_vec(3)
         !!! s(3,1)=q(1,1)*n_vec(2)-q(2,1)*n_vec(1)
         !!! if(norm2(s(:,1))<1.0e-12_dblprec) then
         !!!    !print '(a,3f12.6)','s before:',s
         !!!    s(1,1)=1.0_dblprec
         !!!    s(2,1)=0.0_dblprec
         !!!    s(3,1)=0.0_dblprec
         !!!    !print '(a,3f12.6)','s after:',s
         !!! end if
         !!! s(:,1)=s(:,1)/norm2(s(:,1))

         !!! print *,'----q-and-s-vectors----',1
         !!! print '(3f10.4)', q(:,1)
         !!! print '(3f10.4)', s(:,1)
         !!! print '(3f10.4)', n_vec(:)


         !!!!stop
         ! set up magnetic order 
         ! if qm_rot=y, take the original magnetic order and rotate each spin to 
         ! create spin-spirals in an disordered background
         ! otherwise rotate all spins to create pure spin-spirals
         if (qm_rot=='y') then
            do ia=1,natom
               !
               call f_wrap_coord_diff(natom,coord,ia,countstart+1,srvec)
               !
               qr=q(1,1)*srvec(1)+q(2,1)*srvec(2)+q(3,1)*srvec(3)
               !qr=qpts(1,iq)*srvec(1)+qpts(2,iq)*srvec(2)+qpts(3,iq)*srvec(3)
               qr=2.0_dblprec*pi*qr

               call rodmat(n_vec,qr,r_mat)

               if (.not. qm_excluded_atoms(ia)) then
                  emomm(1:3,ia,k)=matmul(r_mat,emomm_start(:,ia,k))
               else
                  emomm(1:3,ia,k)=emomm_start(:,ia,k)
               end if
            end do
         elseif (qm_cellrot=='y') then
            do ia=1,natom
               !
               !ia_cell=ia/na+1
               ia_cell=((ia-1)/na)*na+1

               call f_wrap_coord_diff(natom,coord,ia_cell,countstart+1,srvec)
               !
               qr=q(1,1)*srvec(1)+q(2,1)*srvec(2)+q(3,1)*srvec(3)
               !qr=qpts(1,iq)*srvec(1)+qpts(2,iq)*srvec(2)+qpts(3,iq)*srvec(3)
               qr=2.0_dblprec*pi*qr

               call rodmat(n_vec,qr,r_mat)

               !print '(a,2i4,3x,3f12.6)','---------------',ia,ia_cell,q(:,1)
               !print '(3f10.4)', r_mat
               !print *,'---------------'

               if (.not. qm_excluded_atoms(ia)) then
                  emomm(1:3,ia,k)=matmul(r_mat,emomm_start(:,ia,k))
               else
                  emomm(1:3,ia,k)=emomm_start(:,ia,k)
               end if
            end do
            !print '(3f12.5)', emomm
         else
            call set_nsvec(qm_type,q(:,1),s(:,1),n_vec)
            !call set_nsvec(qm_type,qpts(:,iq),s(:,1),n_vec)
            !!! print *,'iq:',iq
            !!! print '(a,3f12.6)' , '   iq:', q(:,1)
            !!! !print '(a,3f12.6)' , '   iq:', qpts(:,iq)
            !!! print '(a,3f12.6)' , 's_vec:', s(:,1)
            !!! print '(a,3f12.6)' , 'n_vec:', n_vec
            !!! print '(a,3f12.6)' , 'cross:', f_cross_product(n_vec,s(:,1))

            do ia=1,natom
               !
               !srvec=coord(:,ia)-coord(:,countstart+1)
               ! possible use wrap_coord_diff() here.
               call f_wrap_coord_diff(natom,coord,ia,countstart+1,srvec)
               !
               m_j=0.0_dblprec
               qr=q(1,1)*srvec(1)+q(2,1)*srvec(2)+q(3,1)*srvec(3)
               !qr=qpts(1,iq)*srvec(1)+qpts(2,iq)*srvec(2)+qpts(3,iq)*srvec(3)
               m_j=n_vec*cos(2*pi*qr+phi(1))+s(:,1)*sin(2*pi*qr+phi(1))
               !print '(a,4f12.6)' , 'r_i  :',qr,m_j
               !print '(a,4f12.6)' , 'qr,mj:',qr,m_j
               !call normalize(m_j)
               !emom(1:3,ia,k)=m_j
               if (.not. qm_excluded_atoms(ia)) then
                emomm(1:3,ia,k)=m_j*mmom(ia,k)
                emom(1:3,ia,k)=m_j
               end if
               !print '(a,3f12.6,i8)' , 'emom :', emomm(1:3,ia,k), ia
               !print '(i7,3f12.6)', ia, emomm(1:3,ia,k)
            end do
         end if
         !!! if (qm_rot=='y') then
         !!! else
         !!!    do ia=1,natom
         !!!       !
         !!!       !srvec=coord(:,ia)-coord(:,countstart+1)
         !!!       ! possible use wrap_coord_diff() here.
         !!!       call f_wrap_coord_diff(natom,coord,ia,countstart+1,srvec)
         !!!       !
         !!!       m_j=0.0_dblprec
         !!!       qr=qpts(1,iq)*srvec(1)+qpts(2,iq)*srvec(2)+qpts(3,iq)*srvec(3)
         !!!       m_j=n_vec*cos(2*pi*qr+phi(1))+s(:,1)*sin(2*pi*qr+phi(1))
         !!!       !call normalize(m_j)
         !!!       !emom(1:3,ia,k)=m_j
         !!!       emomm(1:3,ia,k)=m_j*mmom(ia,k)
         !!!       !print '(i7,3f12.6)', ia, emomm(1:3,ia,k)
         !!!    end do
         !!! end if

         if (qm_relax=='Y') call mc_minimal(emomm,emom,mmom,qm_relax_steps,qm_relax_mode,qm_relax_temp)

         ! calculate energy for given q,s,theta combination
         ! anisotropy + external field to be added
         energy=0.0_dblprec
         !call effective_field(natom,mensemble,countstart+1,countstart+na,         &
         call effective_field(natom,mensemble,1,natom, &
            emomm,mmom,external_field,time_external_field,beff,beff1,      &
            beff2,opt_flag,max_no_constellations,maxnoconstl,unitcelltype,        &
            constlncoup,constellations,constellationsneightype,         &
            energy,num_macro,cell_index,emomm_macro,macro_nlistsize,na,n1,n2,n3)

         energy=energy/natom !/mub*mry !/mry*mub/na
         !  print '(a,3f12.6)' , 'ene  :', energy

         call f_wrap_coord_diff(natom,coord,countstart+1,countstart+1,srvec)
         !!! print '(3f10.4,5x,1f10.5,5x,3f10.4,5x,3f10.4,10x,f12.6)',qpts(:,iq),&
         !!!    srvec,&
         !!!    qpts(1,iq)*srvec(1)+qpts(2,iq)*srvec(2)+qpts(3,iq)*srvec(3),&
         !!!    emomm(:,countstart+1,1),energy

         !write(ofileno,'(i8,3g20.8,g20.8)') iq,q(:,1),energy*13.605_dblprec !/mry*mub*13.605_dblprec
         write(ofileno,'(i8,3g20.8,g20.8)') iq,q(:,1),energy !/mry*mub*13.605_dblprec
         !write(ofileno,'(i8,3g20.8,g20.8)') iq,qpts(:,iq),energy*13.605_dblprec !/mry*mub*13.605_dblprec

         ! store best energy configuration if trial energy is lower than minimum
         if(energy<min_energy) then

            min_energy=energy
            !q_best(:,1)=qpts(:,iq)
            q_best(:,1)=q(:,1)
            s_save(:,1)=s(:,1)
            lhit=iter
            nhits=nhits+1
            write(ofileno2,'(i8,3g20.8,g20.8)') iq,q(:,1),energy
            !write(ofileno2,'(i8,3g20.8,g20.8)') iq,qpts(:,iq),energy

         end if
         ! do not loop if diamag_qvect is set
         if (norm2(diamag_qvect)>0.0_dblprec) exit
      end do
      !
      !
      print '(1x,a,i6,a)','line search minimization done with ',nhits,' hits.'
      print '(1x,a)', '|-----minimum energy----|----------------q-vector-----------------|------------------s-vector----------------|'
      do iq=1,1
         print '(2x,f18.10,2x,3f14.6,2x,3f14.6)',min_energy,q_best(:,iq),s(:,iq)
      end do
      ! important: save the lowest energy q-vector
      q=q_best
      ! save best q_vector for nc-ams
      if (norm2(diamag_qvect)==0.0_dblprec) diamag_qvect = q_best(:,1)
      s=s_save
      print '(1x,a)','|-----------------------|-----------------------------------------|------------------------------------------|'

      !
      close(ofileno)
      close(ofileno2)
      !
      if (qm_rot=='y'.or.qm_cellrot=='y') then
            do ia=1,natom
               !
               ia_cell=((ia-1)/na)*na+1

               call f_wrap_coord_diff(natom,coord,ia_cell,countstart+1,srvec)
               !
               qr=q(1,1)*srvec(1)+q(2,1)*srvec(2)+q(3,1)*srvec(3)
               !qr=qpts(1,iq)*srvec(1)+qpts(2,iq)*srvec(2)+qpts(3,iq)*srvec(3)
               qr=2.0_dblprec*pi*qr

               call rodmat(n_vec,qr,r_mat)
               !print *,'---------------',ia,ia_cell
               !print '(3f10.4)', n_vec
               !print *,'---------------'
               !print '(3f10.4)', r_mat
               !print *,'---------------'

               if (.not. qm_excluded_atoms(ia)) then
                  emomm(1:3,ia,k)=matmul(r_mat,emomm_start(:,ia,k))
               else
                  emomm(1:3,ia,k)=emomm_start(:,ia,k)
               end if
            end do
            !print '(3f12.5)', emomm
         i_all=-product(shape(emomm_start))*kind(emomm_start)
         deallocate(emomm_start,stat=i_stat)
         call memocc(i_stat,i_all,'emomm_start','emomm_start')
      end if
      !
      if (qm_relax=='Y') call allocate_mcdata(natom,-1)

      return
      !
   end subroutine sweep_q2


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: sweep_q3
   !> @brief Stupid line search minimization of spin spirals (clone of sweep_q
   !  but for external q-point set.
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine sweep_q3(Natom,Mensemble,NA,coord,emomM,mmom,hfield,OPT_flag,           &
      max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations,    &
      constellationsNeighType,Num_macro,cell_index,emomM_macro,           &
      macro_nlistsize,simid,qpts,nq)

      use RandomNumbers, only: rng_uniform,rng_gaussian
      use InputData, only : N1,N2,N3
      use Math_functions, only : f_wrap_coord_diff
      use Depondt, only : rodmat
      use MomentData, only : emom
      use mc_driver, only : mc_minimal
      use MonteCarlo
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
      !! +++ New variables due to optimization routines +++ !!
      integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
      ! Number of entries (number of unit cells*number of atoms per unit cell) in the constellation matrix per ensemble
      integer, dimension(Mensemble), intent(in) :: maxNoConstl
      ! See OptimizationRoutines.f90 for details on classification
      integer, dimension(Natom, Mensemble), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      ! Matrix relating the interatomic exchanges for each atom in the constellation matrix
      real(dblprec), dimension(ham%max_no_neigh, max_no_constellations,Mensemble), intent(in) :: constlNCoup
      ! Matrix storing all unit cells belonging to any constellation
      real(dblprec), dimension(3,max_no_constellations, Mensemble), intent(in) :: constellations
      ! Optimization flag (1 = optimization on; 0 = optimization off)
      logical, intent(in) :: OPT_flag
      ! Matrix storing the type of the neighbours within a given neighbourhood of a constellation; default is 1 outside the neighbourhood region
      ! The default is to achieve correct indexing. Note here also that constlNCoup will result in a net zero contribution to the Heissenberg exchange term
      integer, dimension(ham%max_no_neigh,max_no_constellations,Mensemble), intent(in) :: constellationsNeighType
      ! Internal effective field arising from the optimization of the Heissenberg exchange term
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      character(len=8), intent(in) :: simid  !< Name of simulation
      integer, intent(in) :: nq  !< number of qpoints
      real(dblprec), dimension(3,nq), intent(in) :: qpts !< Array of q-points
      !
      integer :: iq
      !
      real(dblprec), dimension(3) :: m_j
      real(dblprec) :: pi, qr
      integer :: i, k, ia, lhit, nhits, countstart
      real(dblprec) :: energy, min_energy
      character(len=30) :: filn
      !
      real(dblprec), dimension(3,3) :: q_best
      real(dblprec), dimension(3,3) :: R_mat
      real(dblprec), dimension(3) :: theta_best
      real(dblprec) :: theta_glob_best
      real(dblprec), dimension(3) :: phi_best
      integer :: iter, iscale, i1 ,i2 ,i3, qq
      real(dblprec), dimension(3) :: srvec 
      !
      !
      !
      pi=4._dblprec*ATAN(1._dblprec)
      theta_glob_best=0.0_dblprec
      theta_best=0.0_dblprec
      phi_best=0.0_dblprec
      phi=0.0_dblprec
      !
      ! Normal vector
      ! Read from file or default
      !n_vec(1)=0.0_dblprec;n_vec(2)=0.0_dblprec;n_vec(3)=1.0_dblprec;
      !
      ! Starting atom
      I1 = N1/2
      I2 = N2/2
      I3 = N3/2

      countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
      !
      write(filn,'(''qm_sweep.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position="append")
      write(ofileno,'(a)') "#    Iter                          Q-vector                                 Energy(meV)  "

      write(filn,'(''qm_minima.'',a,''.out'')') trim(simid)
      open(ofileno2,file=filn, position="append")
      write(ofileno2,'(a)') "#    Iter                          Q-vector                                 Energy(mRy)  "
     
    
      if (.not. allocated(qm_excluded_atoms)) call qminimizer_build_exclude_list(Natom)
     
      if (qm_relax=='Y') call allocate_mcdata(Natom,1)
      !
      theta = 0.0_dblprec
      theta_glob = 0.0_dblprec
      phi = 0.0_dblprec
      min_energy=1.0d4

      !
      lhit=0
      nhits=0
      iscale=1

      ! Switch rotation direction

      ! Calculate total external field (not included yet)
      do k=1,Mensemble
         do i=1,Natom
            external_field(1:3,i,k)= hfield
            beff(1:3,i,k)=0.0_dblprec
            beff1(1:3,i,k)=0.0_dblprec
            beff2(1:3,i,k)=0.0_dblprec
         end do
      end do

      ! Only use first ensemble
      k=1
      do iq=1,nq
         iter=iq
         !         !
         ! Set up spin-spiral magnetization (only first cell+ neighbours)
         energy=0.0_dblprec
         !!!print *,'----Q-and-S-vectors----',1
         !!!print '(3f10.4)', q(:,1)
         !!!print '(3f10.4)', s(:,1)
         !!!print '(3f10.4)', n_vec(:)
         q(:,1)=qpts(:,iq)
         ! Rotate 360/iq
         theta(1)=0.0d0
         theta(2)=2.0_dblprec*pi/3.0d0*(1.0_dblprec)
         theta(3)=2.0_dblprec*pi/3.0d0*(2.0_dblprec)
         ! Replace with Rodrigues

         do qq=2,3
            call rodmat(n_vec,theta(qq),R_mat)
            q(:,qq)=matmul(R_mat,q(:,1))
            s(:,qq)=matmul(R_mat,s(:,1))
            !s(:,qq)=s(:,1)
         end do
         
         do ia=1,Natom
            !
            ! Possible use wrap_coord_diff() here.
            call f_wrap_coord_diff(Natom,coord,ia,countstart+1,srvec)
            !
            m_j=0.0_dblprec
            do qq=1,3
               qr=q(1,qq)*srvec(1)+q(2,qq)*srvec(2)+q(3,qq)*srvec(3)
               !m_j=m_j+s(:,qq)*cos(2*pi*qr)+n_vec*sin(2*pi*qr)
               m_j=m_j+n_vec*cos(2*pi*qr)+s(:,qq)*sin(2*pi*qr)
               !m_j=m_j+n_vec*cos(2*pi*qr+phi(1))+s(:,qq)*sin(2*pi*qr+phi(1))
            end do
            call normalize(m_j)
            !emom(1:3,ia,k)=m_j
            if (.not. qm_excluded_atoms(ia)) then
               emomM(1:3,ia,k)=m_j*mmom(ia,k)
               emom(1:3,ia,k)=m_j
            end if
            !print '(i7,3f12.6)', ia, emomM(1:3,ia,k)
         end do

         if (qm_relax=='Y') call mc_minimal(emomm,emom,mmom,qm_relax_steps,qm_relax_mode,qm_relax_temp)

         ! Calculate energy for given q,s,theta combination
         ! Anisotropy + external field to be added
         energy=0.0_dblprec
         !call effective_field(Natom,Mensemble,countstart+1,countstart+na,         &
         time_external_field  = 0.0
         call effective_field(Natom,Mensemble,1,Natom, &
            emomM,mmom,external_field,time_external_field,beff,beff1,      &
            beff2,OPT_flag,max_no_constellations,maxNoConstl,unitCellType,        &
            constlNCoup,constellations,constellationsNeighType,         &
            energy,Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

         energy=energy/Natom !/mub*mry !/mry*mub/NA

         call f_wrap_coord_diff(Natom,coord,countstart+1,countstart+1,srvec)
         !!! print '(3f10.4,5x,1f10.5,5x,3f10.4,5x,3f10.4,10x,f12.6)',qpts(:,iq),&
         !!!    srvec,&
         !!!    qpts(1,iq)*srvec(1)+qpts(2,iq)*srvec(2)+qpts(3,iq)*srvec(3),&
         !!!    emomM(:,countstart+1,1),energy

         write(ofileno,'(i8,3g20.8,g20.8)') iq,qpts(:,iq),energy !*13.605_dblprec !/mry*mub*13.605_dblprec

         ! Store best energy configuration if trial energy is lower than minimum
         if(energy<min_energy) then

            min_energy=energy
            !q_best(:,1)=qpts(:,iq)
            q_best=q
            lhit=iter
            nhits=nhits+1
            write(ofileno2,'(i8,3g20.8,g20.8)') iq,qpts(:,iq),energy

         end if
      end do
      !
      !
      print '(1x,a,i6,a)','Line search minimization done with ',nhits,' hits.'
      print '(1x,a)', '|-----Minimum energy----|----------------Q-vector-----------------|------------------S-vector----------------|'
      do iq=1,3
         print '(2x,f18.10,2x,3f14.6,2x,3f14.6)',min_energy,q_best(:,iq),s(:,iq)
      end do
      ! Important: Save the lowest energy q-vector
      q=q_best
      !s=s_save
      print '(1x,a)','|-----------------------|-----------------------------------------|------------------------------------------|'

      !
      close(ofileno)
      close(ofileno2)
      !
      if (qm_relax=='Y') call allocate_mcdata(Natom,-1)
      !
      return
      !
   end subroutine sweep_q3

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: sweep_cube
   !> @brief Stupid line search minimization of spin spirals (clone of sweep_q
   !  but for external q-point set.
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine sweep_cube(Natom,Mensemble,NA,coord,emomM,mmom,hfield,OPT_flag,           &
      max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations,    &
      constellationsNeighType,Num_macro,cell_index,emomM_macro,           &
      macro_nlistsize,simid,qpts,nq)

      use RandomNumbers, only: rng_uniform,rng_gaussian
      use InputData, only : N1,N2,N3
      use Math_functions, only : f_wrap_coord_diff
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
      !! +++ New variables due to optimization routines +++ !!
      integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
      ! Number of entries (number of unit cells*number of atoms per unit cell) in the constellation matrix per ensemble
      integer, dimension(Mensemble), intent(in) :: maxNoConstl
      ! See OptimizationRoutines.f90 for details on classification
      integer, dimension(Natom, Mensemble), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      ! Matrix relating the interatomic exchanges for each atom in the constellation matrix
      real(dblprec), dimension(ham%max_no_neigh, max_no_constellations,Mensemble), intent(in) :: constlNCoup
      ! Matrix storing all unit cells belonging to any constellation
      real(dblprec), dimension(3,max_no_constellations, Mensemble), intent(in) :: constellations
      ! Optimization flag (1 = optimization on; 0 = optimization off)
      logical, intent(in) :: OPT_flag
      ! Matrix storing the type of the neighbours within a given neighbourhood of a constellation; default is 1 outside the neighbourhood region
      ! The default is to achieve correct indexing. Note here also that constlNCoup will result in a net zero contribution to the Heissenberg exchange term
      integer, dimension(ham%max_no_neigh,max_no_constellations,Mensemble), intent(in) :: constellationsNeighType
      ! Internal effective field arising from the optimization of the Heissenberg exchange term
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      character(len=8), intent(in) :: simid  !< Name of simulation
      integer, intent(in) :: nq  !< number of qpoints
      real(dblprec), dimension(3,nq), intent(in) :: qpts !< Array of q-points
      !
      integer :: iq, jq
      !
      real(dblprec), dimension(3) :: m_j
      real(dblprec) :: pi, qr
      integer :: i, k, ia, lhit, nhits, countstart
      real(dblprec) :: energy, min_energy
      character(len=30) :: filn
      !
      real(dblprec), dimension(3,3) :: q_best
      real(dblprec), dimension(3,3) :: s_save
      real(dblprec), dimension(3) :: theta_best
      real(dblprec) :: theta_glob_best
      real(dblprec), dimension(3) :: phi_best
      integer :: iter, iscale, i1 ,i2 ,i3, qq
      real(dblprec), dimension(3) :: srvec 
      !
      !
      !
      pi=4._dblprec*ATAN(1._dblprec)
      theta_glob_best=0.0_dblprec
      theta_best=0.0_dblprec
      phi_best=0.0_dblprec
      phi=0.0_dblprec
      !
      ! Normal vector
      ! Read from file or default
      !n_vec(1)=0.0_dblprec;n_vec(2)=0.0_dblprec;n_vec(3)=1.0_dblprec;
      !
      ! Starting atom
      I1 = N1/2
      I2 = N2/2
      I3 = N3/2

      countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
      !
      write(filn,'(''qm_sweep.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position="append")
      write(ofileno,'(a)') "#    Iter                          Q-vector                                 Energy(meV)  "

      write(filn,'(''qm_minima.'',a,''.out'')') trim(simid)
      open(ofileno2,file=filn, position="append")
      write(ofileno2,'(a)') "#    Iter                          Q-vector                                 Energy(mRy)  "
     
    
      if (.not. allocated(qm_excluded_atoms)) call qminimizer_build_exclude_list(Natom)
      
      !
      theta = 0.0_dblprec
      theta_glob = 0.0_dblprec
      phi = 0.0_dblprec
      min_energy=1.0d4

      !
      lhit=0
      nhits=0
      iscale=1

      ! Switch rotation direction

      ! Calculate total external field (not included yet)
      do k=1,Mensemble
         do i=1,Natom
            external_field(1:3,i,k)= hfield
            beff(1:3,i,k)=0.0_dblprec
            beff1(1:3,i,k)=0.0_dblprec
            beff2(1:3,i,k)=0.0_dblprec
         end do
      end do

      ! Only use first ensemble
      k=1
      ! Hand hacked pitch/normal vectors
      q=0.0d0
      n_vec(1)=0.0d0
      n_vec(2)=0.0d0
      n_vec(3)=1.0d0
      do iq=1,1
         iter=iq
         !         !
         do jq=2, 2
            ! Set up spin-spiral magnetization (only first cell+ neighbours)
            energy=0.0_dblprec
            !!!print *,'----Q-and-S-vectors----',1
            !!!print '(3f10.4)', q(:,1)
            !!!print '(3f10.4)', s(:,1)
            !!!print '(3f10.4)', n_vec(:)
            q(:,1)=qpts(:,iq)
            q(:,2)=qpts(:,jq)

            s(1,1)=q(2,1)*n_vec(3)-q(3,1)*n_vec(2)
            s(2,1)=q(3,1)*n_vec(1)-q(1,1)*n_vec(3)
            s(3,1)=q(1,1)*n_vec(2)-q(2,1)*n_vec(1)
            s(:,1)=s(:,1)/norm2(s(:,1))

            s(1,2)=q(2,2)*n_vec(3)-q(3,2)*n_vec(2)
            s(2,2)=q(3,2)*n_vec(1)-q(1,2)*n_vec(3)
            s(3,2)=q(1,2)*n_vec(2)-q(2,2)*n_vec(1)
            s(:,2)=s(:,2)/norm2(s(:,2))
            !!! ! Rotate 360/iq
            !!! theta(1)=0.0d0
            !!! theta(2)=2.0_dblprec*pi/3.0d0*(1.0_dblprec)
            !!! theta(3)=2.0_dblprec*pi/3.0d0*(2.0_dblprec)
            !!! do qq=2,3
            !!!    rx= q(1,1)*cos(theta(qq))+q(2,1)*sin(theta(qq))
            !!!    ry=-q(1,1)*sin(theta(qq))+q(2,1)*cos(theta(qq))
            !!!    q(1,qq)=rx
            !!!    q(2,qq)=ry
            !!!    q(3,qq)=0.0_dblprec
            !!!    rx= s(1,1)*cos(theta(qq))+s(2,1)*sin(theta(qq))
            !!!    ry=-s(1,1)*sin(theta(qq))+s(2,1)*cos(theta(qq))
            !!!    s(1,qq)=rx
            !!!    s(2,qq)=ry
            !!!    !s(:,qq)=s(:,1)
            !!! end do
            !!!!stop
            print *,'----Q-and-S-vectors----',1
            print '(3f10.4)', q(:,1)
            print '(3f10.4)', s(:,1)
            print '(3f10.4)', n_vec(:)
            print *,'----Q-and-S-vectors----',2
            print '(3f10.4)', q(:,2)
            print '(3f10.4)', s(:,2)
            print '(3f10.4)', n_vec(:)
            !!! print *,'----Q-and-S-vectors----',3
            !!! print '(3f10.4)', q(:,3)
            !!! print '(3f10.4)', s(:,3)
            !!! print '(3f10.4)', n_vec(:)
            do ia=1,Natom
               !
               !srvec=coord(:,ia)-coord(:,countstart+1)
               ! Possible use wrap_coord_diff() here.
               call f_wrap_coord_diff(Natom,coord,ia,countstart+1,srvec)
               !
               m_j=0.0_dblprec
               do qq=1,2
                  !qr=qpts(1,iq)*srvec(1)+qpts(2,iq)*srvec(2)+qpts(3,iq)*srvec(3)
                  qr=q(1,qq)*srvec(1)+q(2,qq)*srvec(2)+q(3,qq)*srvec(3)
                  m_j=m_j+n_vec*cos(2*pi*qr)+s(:,qq)*sin(2*pi*qr)
                  !m_j=m_j+n_vec*cos(2*pi*qr+phi(1))+s(:,qq)*sin(2*pi*qr+phi(1))
               end do
               call normalize(m_j)
               !emom(1:3,ia,k)=m_j
               emomM(1:3,ia,k)=m_j*mmom(ia,k)
               !print '(i7,3f12.6)', ia, emomM(1:3,ia,k)
            end do

            ! Calculate energy for given q,s,theta combination
            ! Anisotropy + external field to be added
            energy=0.0_dblprec
            !call effective_field(Natom,Mensemble,countstart+1,countstart+na,         &
            call effective_field(Natom,Mensemble,1,Natom, &
               emomM,mmom,external_field,time_external_field,beff,beff1,      &
               beff2,OPT_flag,max_no_constellations,maxNoConstl,unitCellType,        &
               constlNCoup,constellations,constellationsNeighType,         &
               energy,Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)

            energy=energy/Natom !/mub*mry !/mry*mub/NA

            call f_wrap_coord_diff(Natom,coord,countstart+1,countstart+1,srvec)
            !!! print '(3f10.4,5x,1f10.5,5x,3f10.4,5x,3f10.4,10x,f12.6)',qpts(:,iq),&
            !!!    srvec,&
            !!!    qpts(1,iq)*srvec(1)+qpts(2,iq)*srvec(2)+qpts(3,iq)*srvec(3),&
            !!!    emomM(:,countstart+1,1),energy

            write(ofileno,'(i8,3g20.8,g20.8)') iq,qpts(:,iq),energy*13.605_dblprec !/mry*mub*13.605_dblprec

            ! Store best energy configuration if trial energy is lower than minimum
            if(energy<min_energy) then

               min_energy=energy
               !q_best(:,1)=qpts(:,iq)
               q_best=q
               s_save=s
               lhit=iter
               nhits=nhits+1
               write(ofileno2,'(i8,3g20.8,g20.8)') iq,qpts(:,iq),energy
               write(ofileno2,'(i8,3g20.8,g20.8)') jq,qpts(:,jq),energy

            end if
         end do
      end do
      !
      !
      print '(1x,a,i6,a)','Line search minimization done with ',nhits,' hits.'
      print '(1x,a)', '|-----Minimum energy----|----------------Q-vector-----------------|------------------S-vector----------------|'
      do iq=1,2
         print '(2x,f18.10,2x,3f14.6,2x,3f14.6)',min_energy,q_best(:,iq),s(:,iq)
      end do
      ! Important: Save the lowest energy q-vector
      q=q_best
      s=s_save
      print '(1x,a)','|-----------------------|-----------------------------------------|------------------------------------------|'

      !
      close(ofileno)
      close(ofileno2)
      !
      !
      return
      !
   end subroutine sweep_cube

   !> Plot final configuration
   subroutine plot_q(Natom, Mensemble, coord, emom, emomM, mmom,simid)
      !
      use math_functions, only : f_wrap_coord_diff
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      character(len=8), intent(in) :: simid  !< Name of simulation
      !
      real(dblprec), dimension(3) :: srvec, m_j
      integer :: lhit, ia, k, iq
      real(dblprec) :: pi, qr
      character(len=30) :: filn
      !
      !
      pi=4._dblprec*ATAN(1._dblprec)
      !
      !nplot=12
      !q=10.0_dblprec*q
      !
      write(filn,'(''qm_restart.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn)
      !write(ofileno,*) 0
      write(ofileno,'(a)') repeat("#",80)
      write(ofileno,'(a,1x,a)') "# File type:", 'R'
      write(ofileno,'(a,1x,a)') "# Simulation type:", 'Q'
      write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
      write(ofileno,'(a,1x,i8)')"# Number of ensembles: ", Mensemble
      write(ofileno,'(a)') repeat("#",80)
      write(ofileno,'(a8,a,a8,a16,a16,a16,a16)') "#Timestep","ens","iatom","|Mom|","M_x","M_y","M_z"
      if(qm_rot == 'Y'.or.qm_cellrot == 'Y') then
         do k=1,Mensemble
            do ia=1,Natom
               emom(1:3,ia,k)=emomM(1:3,ia,k)/mmom(ia,k)
               m_j=emom(1:3,ia,k)
               write(ofileno,10003) 0,1,ia,mmom(ia,k),m_j
            end do
         end do
      else
         do k=1,Mensemble
            do ia=1,Natom
               lhit=lhit+1
               !srvec=coord(:,ia)
               call f_wrap_coord_diff(Natom,coord,ia,1,srvec)
               !
               m_j=0.0_dblprec
               do iq=1,1!nq
                  call set_nsvec(qm_type,q(:,iq),s(:,iq),n_vec)
                  qr=q(1,iq)*srvec(1)+q(2,iq)*srvec(2)+q(3,iq)*srvec(3)
                  m_j=m_j+n_vec*cos(2*pi*qr+phi(iq))+s(:,iq)*sin(2*pi*qr+phi(iq))
               end do
               call normalize(m_j)
               if (.not. qm_excluded_atoms(ia)) then
                  emom(1:3,ia,k)=m_j
                  emomM(1:3,ia,k)=m_j*mmom(ia,k)
               end if
               !write(ofileno,'(2i8,4f14.6)') 1,lhit,mmom(ia,k),m_j
               !write(ofileno,'(2i8,4f14.6)') 1,ia,mmom(ia,k),m_j
               write(ofileno,10003) 0,1,ia,mmom(ia,k),emom(:,ia,k)
            end do
         end do
      end if
      close(ofileno)
      !
      10003 format(i8,i8,i8,2x,es16.8,es16.8,es16.8,es16.8)
      !10003 format(es16.8,i8,i8,2x,es16.8,es16.8,es16.8,es16.8)
      !
      !
   end subroutine plot_q


   !> Plot final configuration
   subroutine plot_q3(Natom, Mensemble, coord, emom, emomM, mmom,simid)
      use Math_functions, only : f_wrap_coord_diff
      !
      use InputData, only : N1,N2,N3, NA
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      character(len=8), intent(in) :: simid  !< Name of simulation
      !
      real(dblprec), dimension(3) :: srvec, m_j
      integer :: lhit, ia, k, iq
      integer :: I1, I2, I3, countstart
      real(dblprec) :: pi, qr, totene
      character(len=30) :: filn
      !
      !
      pi=4._dblprec*ATAN(1._dblprec)
      !
      I1 = N1/2
      I2 = N2/2
      I3 = N3/2

      countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
      !nplot=12
      !q=10.0_dblprec*q
      !
      write(filn,'(''qm_restart.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn)
      !write(ofileno,*) 0
      write(ofileno,'(a)') repeat("#",80)
      write(ofileno,'(a,1x,a)') "# File type:", 'R'
      write(ofileno,'(a,1x,a)') "# Simulation type:", 'Q'
      write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
      write(ofileno,'(a,1x,i8)')"# Number of ensembles: ", Mensemble
      write(ofileno,'(a)') repeat("#",80)
      write(ofileno,'(a8,a,a8,a16,a16,a16,a16)') "#Timestep","ens","iatom","|Mom|","M_x","M_y","M_z"
      do k=1,Mensemble
         do ia=1,Natom
            lhit=lhit+1
            !srvec=coord(:,ia)
            call f_wrap_coord_diff(Natom,coord,ia,countstart+1,srvec)
            !call wrap_coord_diff(Natom,coord,ia,1,srvec)
            !
            m_j=0.0_dblprec
            do iq=1,3
               qr=q(1,iq)*srvec(1)+q(2,iq)*srvec(2)+q(3,iq)*srvec(3)
               m_j=m_j+n_vec*cos(2*pi*qr)+s(:,iq)*sin(2*pi*qr)
               !m_j=m_j+n_vec*cos(2*pi*qr+phi(iq))+s(:,iq)*sin(2*pi*qr+phi(iq))
            end do
            call normalize(m_j)
            if (.not. qm_excluded_atoms(ia)) then
               emom(1:3,ia,k)=m_j
               emomM(1:3,ia,k)=m_j*mmom(ia,k)
            end if
            !write(ofileno,'(2i8,4f14.6)') 1,lhit,mmom(ia,k),m_j
            !write(ofileno,'(2i8,4f14.6)') 1,ia,mmom(ia,k),m_j
            write(ofileno,10003) 0,1,ia,mmom(ia,k),emom(:,ia,k)
            !print '(i6, f12.6, 4x, 3f12.6)',ia,mmom(ia,k),emom(:,ia,k)
         end do
      end do
      close(ofileno)
      !
      !
      10003 format(i8,i8,i8,2x,es16.8,es16.8,es16.8,es16.8)
      !
   end subroutine plot_q3

   !> Plot final configuration
   subroutine plot_cube(Natom, Mensemble, coord, emom, emomM, mmom,simid)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      character(len=8), intent(in) :: simid  !< Name of simulation
      !
      real(dblprec), dimension(3) :: srvec, m_j
      integer :: lhit, ia, k, iq
      real(dblprec) :: pi, qr
      character(len=30) :: filn
      !
      !
      pi=4._dblprec*ATAN(1._dblprec)
      !
      !nplot=12
      !q=10.0_dblprec*q
      !
      write(filn,'(''qm_restart.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn)
      !write(ofileno,*) 0
      write(ofileno,'(a)') repeat("#",80)
      write(ofileno,'(a,1x,a)') "# File type:", 'R'
      write(ofileno,'(a,1x,a)') "# Simulation type:", 'Q'
      write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
      write(ofileno,'(a,1x,i8)')"# Number of ensembles: ", Mensemble
      write(ofileno,'(a)') repeat("#",80)
      write(ofileno,'(a8,a,a8,a16,a16,a16,a16)') "#Timestep","ens","iatom","|Mom|","M_x","M_y","M_z"
      do k=1,Mensemble
         do ia=1,Natom
            lhit=lhit+1
            srvec=coord(:,ia)
            !
            m_j=0.0_dblprec
            do iq=1,3
               qr=q(1,iq)*srvec(1)+q(2,iq)*srvec(2)+q(3,iq)*srvec(3)
               m_j=m_j+n_vec*cos(2*pi*qr)+s(:,iq)*sin(2*pi*qr)
               !m_j=m_j+n_vec*cos(2*pi*qr+phi(iq))+s(:,iq)*sin(2*pi*qr+phi(iq))
            end do
            call normalize(m_j)
            emom(1:3,ia,k)=m_j
            emomM(1:3,ia,k)=m_j*mmom(ia,k)
            !write(ofileno,'(2i8,4f14.6)') 1,lhit,mmom(ia,k),m_j
            !write(ofileno,'(2i8,4f14.6)') 1,ia,mmom(ia,k),m_j
            write(ofileno,10003) 0,1,ia,mmom(ia,k),m_j
         end do
      end do
      close(ofileno)
      !
      10003 format(i8,i8,i8,2x,es16.8,es16.8,es16.8,es16.8)
      !10003 format(es16.8,i8,i8,2x,es16.8,es16.8,es16.8,es16.8)
      !
      !
   end subroutine plot_cube

   !> Rotation of vectors
   subroutine rotvec(s,m)
      !
      implicit none
      !
      real(dblprec), dimension(3), intent(in) :: s
      real(dblprec), dimension(3), intent(out) :: m
      !
      real(dblprec) :: theta, dot, qnorm
      real(dblprec), dimension(3) :: u,q_new, q
      real(dblprec), dimension(3,3) :: I,ux,uplus, R
      !
      I=0.0_dblprec;I(1,1)=1.0_dblprec;I(2,2)=1.0_dblprec;I(3,3)=1.0_dblprec
      q=(/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
      qnorm=1.0_dblprec
      !
      ! Perpendicular vector
      u(1)=q(2)*s(3)-q(3)*s(2)
      u(2)=q(3)*s(1)-q(1)*s(3)
      u(3)=q(1)*s(2)-q(2)*s(1)
      !
      uplus(1,1)=u(1)*u(1)
      uplus(2,1)=u(1)*u(2)
      uplus(3,1)=u(1)*u(3)
      uplus(1,2)=u(2)*u(1)
      uplus(2,2)=u(2)*u(2)
      uplus(3,2)=u(2)*u(3)
      uplus(1,3)=u(3)*u(1)
      uplus(2,3)=u(3)*u(2)
      uplus(3,3)=u(3)*u(3)
      !
      ux=0.0_dblprec
      ux(2,1)=-u(3)
      ux(3,1)= u(2)
      ux(1,2)= u(3)
      ux(3,2)=-u(1)
      ux(1,3)=-u(2)
      ux(2,3)= u(1)
      !
      dot=q(1)*s(1)+q(2)*s(2)+q(3)*s(3)
      dot=dot/qnorm
      theta=acos(dot)
      !
      R=cos(theta)*I+sin(theta)*ux+(1.0_dblprec-cos(theta))*uplus
      !
      q_new(1)=R(1,1)*m(1)+R(2,1)*m(2)+R(3,1)*m(3)
      q_new(2)=R(1,2)*m(1)+R(2,2)*m(2)+R(3,2)*m(3)
      q_new(3)=R(1,3)*m(1)+R(2,3)*m(2)+R(3,3)*m(3)
      !
      m=q_new
      !
      return
      !
   end subroutine rotvec

   !> Normalization
   subroutine normalize(v)
      !
      implicit none
      !
      real(dblprec),dimension(3), intent(inout) :: v
      !
      real(dblprec) :: vnorm
      !
      vnorm=sqrt(sum(v**2))
      if(vnorm>0.0_dblprec)  v=v/vnorm
      !
      return
   end subroutine normalize


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: qmc
   !> @ brief Energy minimization
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine qmc(Natom,Mensemble,NA,N1,N2,N3,coord,     &
      emomM,mmom,hfield,OPT_flag,     &
      max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations,    &
      constellationsNeighType,Num_macro,cell_index,emomM_macro,           &
      macro_nlistsize)
      !
      use Constants, only: k_bolt, mub
      use InputData, only: Temp
      use RandomNumbers, only : rng_uniform
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1
      integer, intent(in) :: N2
      integer, intent(in) :: N3
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
      !! +++ New variables due to optimization routines +++ !!
      integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
      ! Number of entries (number of unit cells*number of atoms per unit cell) in the constellation matrix per ensemble
      integer, dimension(Mensemble), intent(in) :: maxNoConstl
      ! See OptimizationRoutines.f90 for details on classification
      integer, dimension(Natom, Mensemble), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      ! Matrix relating the interatomic exchanges for each atom in the constellation matrix
      real(dblprec), dimension(ham%max_no_neigh, max_no_constellations,Mensemble), intent(in) :: constlNCoup
      ! Matrix storing all unit cells belonging to any constellation
      real(dblprec), dimension(3,max_no_constellations, Mensemble), intent(in) :: constellations
      ! Optimization flag (1 = optimization on; 0 = optimization off)
      logical, intent(in) :: OPT_flag
      ! Matrix storing the type of the neighbours within a given neighbourhood of a constellation; default is 1 outside the neighbourhood region
      ! The default is to achieve correct indexing. Note here also that constlNCoup will result in a net zero contribution to the Heissenberg exchange term
      integer, dimension(ham%max_no_neigh,max_no_constellations,Mensemble), intent(in) :: constellationsNeighType
      ! Internal effective field arising from the optimization of the Heissenberg exchange term
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize

      ! .. Local variables
      integer :: iq
      !
      real(dblprec), dimension(3) :: m_i, m_j
      real(dblprec) :: pi, qr, rn(3)
      integer :: i,j, k, ia, ja, lhit, nhits
      real(dblprec) :: energy, old_energy, flipprob(1), de, beta
      real(dblprec), dimension(3) :: avgmom
      real(dblprec) :: avgM
      !
      real(dblprec), dimension(3,3) :: q_trial, s_trial
      !real(dblprec), dimension(3,3) :: s_save, q_save
      real(dblprec), dimension(3) :: theta_trial, phi_trial
      !real(dblprec), dimension(3) :: phi_save, phi_diff, phi_best
      real(dblprec), dimension(3,Natom,Mensemble) :: emomM_tmp
      !
      integer :: niter, iter, iscale

      !
      pi=4._dblprec*ATAN(1._dblprec)
      !
      niter=15000000
      do iq=1,nq
         s(1,iq)=0.0_dblprec;s(2,iq)=0.0_dblprec;s(3,iq)=1.0_dblprec
         theta(iq)=90
         theta(iq)=theta(iq)*pi/180.0_dblprec
         phi(iq)=0.0_dblprec
         !
         q(1,iq)=0.0_dblprec;q(2,iq)=0.0_dblprec;q(3,iq)=1.0_dblprec
      end do

      lhit=0
      nhits=0
      iscale=1

      fac_2d=1.0_dblprec
      if(sum(coord(3,:)**2)==0) fac_2d=0.0_dblprec

      ! Calculate total external field
      do k=1,Mensemble
         do i=1,Natom
            external_field(1:3,i,k)=hfield
            beff(1:3,i,k)=0.0_dblprec
            beff1(1:3,i,k)=0.0_dblprec
            beff2(1:3,i,k)=0.0_dblprec
         end do
      end do

      energy=0.0_dblprec
      call effective_field(Natom,Mensemble,1,na, &
         emomM,mmom,external_field,       &
         time_external_field,beff,beff1,beff2,OPT_flag,max_no_constellations,       &
         maxNoConstl,unitCellType,constlNCoup,constellations,                       &
         constellationsNeighType,energy,Num_macro,cell_index,emomM_macro, &
         macro_nlistsize,NA,N1,N2,N3)
      old_energy=energy
      print *, 'Starting energy:',energy
      do iter=1,niter
         !
         !
         !
         do iq=1,nq
            ! delta q-vector
            call rng_uniform(rn,3)
            q_trial(1,iq)=2.0_dblprec*rn(1)-1.0_dblprec
            q_trial(2,iq)=2.0_dblprec*rn(2)-1.0_dblprec
            q_trial(3,iq)=2.0_dblprec*rn(3)-1.0_dblprec
            !
            ! delta s-vector
            call rng_uniform(rn,3)
            s_trial(1,iq)=2.0_dblprec*rn(1)-1.0_dblprec
            s_trial(2,iq)=2.0_dblprec*rn(2)-1.0_dblprec
            s_trial(3,iq)=2.0_dblprec*rn(3)-1.0_dblprec
            !
            ! theta angle
            call rng_uniform(rn,3)
            theta_trial(iq)=(rn(1)-0.5_dblprec)*pi
            ! phi angle
            phi_trial(iq)=(rn(2)-0.5_dblprec)*pi
            !
         end do
         ! Set up trial vectors
         do iq=1,nq
            call normalize(s_trial(1:3,iq))
         end do
         !
         avgmom=0.0_dblprec
         do k=1,Mensemble
            do ia=1,NA
               avgmom=avgmom+emomM(1:3,ia,k)/ham%nlistsize(ia)/NA
               do j=1,ham%nlistsize(ia)
                  ja=ham%nlist(j,ia)
                  avgmom=avgmom+emomM(1:3,ja,k)/ham%nlistsize(ia)/NA
               end do
            end do
         end do
         avgm=(sum(avgmom)**2)**0.5_dblprec
         ! Set up spin-spiral magnetization (only first cell+ neighbours)
         do k=1,Mensemble
            energy=0.0_dblprec
            do ia=1,NA
               m_i=emomM(1:3,ia,k)
               do iq=1,nq
                  qr=q_trial(1,iq)*coord(1,ia)+q_trial(2,iq)*coord(2,ia)+q_trial(3,iq)*coord(3,ia)
                  call normalize(m_i)
                  m_i=m_i*mmom(ia,k)/nq
                  m_i(1)=m_i(1)+cos(2*pi*qr+phi_trial(iq))*sin(theta_trial(iq))*mmom(ia,k)/nq
                  m_i(2)=m_i(2)+sin(2*pi*qr+phi_trial(iq))*sin(theta_trial(iq))*mmom(ia,k)/nq
                  m_i(3)=m_i(3)+cos(theta_trial(iq))*mmom(ia,k)/nq
                  call rotvec(s_trial(1:3,iq),m_i)
                  call normalize(m_i)
                  m_i=m_i*mmom(ia,k)/nq
               end do
               !energy=energy+kani_cell(ia)*(1-(eani_cell(1,ia)*m_i(1)+eani_cell(2,ia)*m_i(2)+eani_cell(3,ia)*m_i(3))**2)
               emomM_tmp(1:3,ia,k)=emomM(1:3,ia,k)
               emomM(1:3,ia,k)=m_i
               do j=1,ham%nlistsize(ia)
                  ja=ham%nlist(j,ia)
                  m_j=emomM(1:3,ja,k)
                  do iq=1,nq
                     call normalize(m_j)
                     m_j=m_j*mmom(ja,k)/nq
                     qr=q_trial(1,iq)*coord(1,ja)+q_trial(2,iq)*coord(2,ja)+q_trial(3,iq)*coord(3,ja)
                     m_j(1)=m_j(1)+cos(2*pi*qr+phi_trial(iq))*sin(theta_trial(iq))*mmom(ja,k)/nq
                     m_j(2)=m_j(2)+sin(2*pi*qr+phi_trial(iq))*sin(theta_trial(iq))*mmom(ja,k)/nq
                     m_j(3)=m_j(3)+cos(theta_trial(iq))*mmom(ja,k)/nq
                     call rotvec(s_trial(1:3,iq),m_j)
                     call normalize(m_j)
                     m_j=m_j*mmom(ja,k)/nq
                  end do
                  emomM_tmp(1:3,ja,k)=emomM(1:3,ja,k)
                  emomM(1:3,ja,k)=m_j
               end do
               ! Calculate energy for given q,s,theta combination
               call effective_field(Natom,Mensemble,ia,ia,emomM,mmom,   &
                  external_field,time_external_field,beff,beff1,beff2,OPT_flag,     &
                  max_no_constellations,maxNoConstl,unitCellType,constlNCoup,       &
                  constellations,constellationsNeighType,energy,Num_macro,&
                  cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
               energy = -energy
            end do
         end do
         ! Store best energy configuration
         call rng_uniform(flipprob,1)
         beta=1_dblprec/k_bolt/Temp
         de=(energy-old_energy)*mub
         if(de<=0.0_dblprec .or. flipprob(1)<=exp(-beta*de)) then
            nhits=nhits+1
            old_energy=energy
            print '(1x,a,i8,4f24.6)', 'QMC: ',nhits,avgm, energy, old_energy, de
            print '(1x,a,3f12.6)', '     ',emomM(1:3,1,1)-emomM(1:3,2,1)
         else
            avgmom=0.0_dblprec
            do k=1,Mensemble
               do ia=1,NA
                  emomM(1:3,ia,k)=emomM_tmp(1:3,ia,k)
                  avgmom=avgmom+emomM(1:3,ia,k)/ham%nlistsize(ia)/NA
                  do j=1,ham%nlistsize(ia)
                     ja=ham%nlist(j,ia)
                     emomM(1:3,ja,k)=emomM_tmp(1:3,ja,k)
                     avgmom=avgmom+emomM(1:3,ja,k)/ham%nlistsize(ia)/NA
                  end do
               end do
            end do
            avgm=(sum(avgmom)**2)**0.5_dblprec
         end if

      end do
      !
      return
      !
   end subroutine qmc

   !---------------------------------------------------------------------------
   !> @brief
   !> Setup lookup table for exclusion of spins in minimization
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine qminimizer_build_exclude_list(Natom)
      use SystemData, only : atype

      implicit none

      integer :: i_stat, iatom
      integer, intent(in) :: Natom !< Number of atoms in system

      allocate(qm_excluded_atoms(natom),stat=i_stat)
      call memocc(i_stat,product(shape(qm_excluded_atoms))*kind(qm_excluded_atoms),'qm_excluded_atoms','qminimizer_build_exclude_list')


      if (allocated(qm_excluded_types)) then
         do iatom = 1,Natom
            qm_excluded_atoms(iatom) = any(atype(iatom) == qm_excluded_types)
         end do
      else
            qm_excluded_atoms = .false.
      end if

   end subroutine qminimizer_build_exclude_list

   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_qminimizer(ifile)

      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
      integer :: rd_len, i_err, i_errb, i_stat, ii
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
         case('qm_min')
            read(ifile,*,iostat=i_err) q_min
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_max')
            read(ifile,*,iostat=i_err) q_max
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_type')
            read(ifile,*,iostat=i_err) qm_type
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_cellrot')
            read(ifile,*,iostat=i_err) qm_cellrot
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_rot')
            read(ifile,*,iostat=i_err) qm_rot
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_oaxis')
            read(ifile,*,iostat=i_err) qm_oaxis
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_qvec')
            read(ifile,*,iostat=i_err) q(:,1:nq)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_svec')
            read(ifile,*,iostat=i_err) s(:,1:nq)
            s(:,1)=s(:,1)/norm2(s(:,1)+dbl_tolerance)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_nvec')
            read(ifile,*,iostat=i_err) n_vec
            n_vec=n_vec/norm2(n_vec+dbl_tolerance)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_relax')
            read(ifile,*,iostat=i_err) qm_relax
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_relax_mode')
            read(ifile,*,iostat=i_err) qm_relax_mode
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_relax_steps')
            read(ifile,*,iostat=i_err) qm_relax_steps
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_relax_temp')
            read(ifile,*,iostat=i_err) qm_relax_temp
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_nstep')
            read(ifile,*,iostat=i_err) nstep
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('qm_exclude')
            read(ifile,*,iostat=i_err) qm_no_excluded
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            if (qm_no_excluded>0) then
               allocate(qm_excluded_types(qm_no_excluded),stat=i_stat)
               call memocc(i_stat,product(shape(qm_excluded_types))*kind(qm_excluded_types),'qm_excluded','read_parameters_qminimizer')
               do ii=1,qm_no_excluded
                  read(ifile,*,iostat=i_err) qm_excluded_types(ii)
                  if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               end do
            end if
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
end subroutine read_parameters_qminimizer

subroutine qminimizer_init()
   !
   implicit none
   !
   q_min=-1.0_dblprec
   q_max=1.0_dblprec
   q=0.0_dblprec
   q(3,1)=1.0_dblprec
   n_vec=0.0_dblprec
   n_vec(3)=1.0_dblprec
   s=0.0_dblprec
   nstep=1000
   qm_cellrot='N'
   qm_rot='N'
   qm_oaxis='N'
   qm_type='N'
   qm_relax='N'
   qm_relax_mode='H'
   qm_relax_steps=100
   qm_relax_temp=1.0e-3_dblprec
   qm_no_excluded=0

end subroutine qminimizer_init

subroutine set_nsvec(qm_type,q_vec,s_vec,n_vec)
   use RandomNumbers, only : rng_uniform
   use math_functions, only : f_cross_product

   !
   implicit none
   !
   character*1, intent(in) :: qm_type
   real(dblprec), dimension(3), intent(in) :: q_vec
   real(dblprec), dimension(3), intent(inout) :: s_vec
   real(dblprec), dimension(3), intent(inout) :: n_vec
   !
   real(dblprec), dimension(3) :: r_vec
   real(dblprec), dimension(3) :: c_vec
   real(dblprec) :: v_norm, r_norm
   !
   if(qm_type=='C') then
      n_vec = q_vec
      v_norm=norm2(n_vec)
      if(v_norm==0.0d0) then
         n_vec(1)=0.0_dblprec;n_vec(2)=0.0_dblprec;n_vec(3)=1.0_dblprec
      else
         n_vec = n_vec/v_norm
      end if

      v_norm = 0.0_dblprec
      do while (v_norm<1e-6) 
         call rng_uniform(r_vec,3)
         r_vec = 2.0_dblprec*r_vec - 1.0_dblprec
         r_norm = norm2(r_vec)
         r_vec = r_vec/r_norm

         c_vec =  f_cross_product(n_vec,r_vec)
         v_norm = norm2(c_vec)
      end do
      s_vec = c_vec/v_norm

   else if(qm_type=='H') then
      v_norm = 0.0_dblprec
      do while (v_norm<1e-6) 
         call rng_uniform(r_vec,3)
         r_vec = 2.0_dblprec*r_vec - 1.0_dblprec
         r_norm = norm2(r_vec)
         r_vec = r_vec/r_norm

         c_vec =  f_cross_product(q_vec,r_vec)
         v_norm = norm2(c_vec)
      end do
      s_vec = c_vec/v_norm
      n_vec = f_cross_product(q_vec,s_vec)
      v_norm = norm2(n_vec)
      n_vec = n_vec/v_norm
   end if

end subroutine set_nsvec
end module qminimizer
