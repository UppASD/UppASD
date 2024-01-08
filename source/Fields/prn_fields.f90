!> Data structures for printing the effective and thermal fields
!> @author
!> Jonathan Chico, Johan Hellsvik, Anders Bergman
!> @copyright
!> GNU Public License.
module prn_fields

   use Profiling
   use Constants
   use Parameters

   implicit none

   ! Input parameters to be read
   integer :: beff_step                    !< Interval between consecutive prints of the total effective field
   integer :: beff_buff                    !< Buffer size for the total field
   integer :: binteff_step                 !< Interval between consecutive prints of the internal effective field
   integer :: binteff_buff                 !< Buffer size for the internal field
   integer :: thermfield_step              !< Interval between thermal field trajectories
   integer :: thermfield_buff              !< Buffer size for the stochastic field
   integer :: torques_step                 !< Interval between consecutive prints of the resulting torques
   integer :: torques_buff                 !< Buffer size for the resulting torques
   integer :: spin_torques_step            !< Interval between consecutive prints of the resulting torques
   integer :: spin_torques_buff            !< Buffer size for the resulting torques
   integer :: larm_step                    !< Interval between consecutive prints of the larmor frequencies
   integer :: larm_buff                    !< Buffer size for the larmor frequencies
   integer :: larm_dos_size                !< Number of windows for larmor dos histogram
   character(len=1) :: do_prn_beff         !< Flag governing file output of total effective fields (Y/N)
   character(len=1) :: do_prn_binteff      !< Flag governing file output of internal effective fields (Y/N)
   character(len=1) :: do_prn_torques      !< Flag governing file output of resulting torques (Y/N)
   character(len=1) :: do_prn_spin_torques !< Flag governing file output of spin torques (Y/N)
   character(len=1) :: do_thermfield       !< Thermal fields trajectory
   character(len=1) :: do_larmor_loc       !< Calculate local precession frequencies from local field (Y/N)
   character(len=1) :: do_larmor_dos       !< Calculate average precession frequencies from local field (Y/N)

   ! Local variables for buffering and indexing of fields
   integer :: bcount_beff         !< Counter of buffer for total field
   integer :: bcount_binteff      !< Counter of buffer for internal field
   integer :: bcount_torques      !< Counter of buffer for torques
   integer :: bcount_spin_torques !< Counter of buffer for spin torques
   integer :: bcount_therm        !< Counter of buffer for stochastic field
   integer :: bcount_larm         !< Counter of buffer for Larmor frequencies
   real(dblprec), dimension(:), allocatable :: indxb_beff           !< Step counter for total field
   real(dblprec), dimension(:), allocatable :: indxb_binteff        !< Step counter for internal field
   real(dblprec), dimension(:), allocatable :: indxb_torques        !< Step counter for resulting torques
   real(dblprec), dimension(:), allocatable :: indxb_spin_torques   !< Step counter for spin torques
   real(dblprec), dimension(:), allocatable :: indxb_larm           !< Step counter for total field
   real(dblprec), dimension(:), allocatable :: indxb_therm          !< Step counter for stochastic field
   real(dblprec), dimension(:,:,:,:), allocatable :: beffb          !< Buffer the site dependent total field
   real(dblprec), dimension(:,:,:,:), allocatable :: binteffb       !< Buffer the site resulting torques
   real(dblprec), dimension(:,:,:,:), allocatable :: torquesb       !< Buffer the site dependent internal field
   real(dblprec), dimension(:,:,:,:), allocatable :: spin_torquesb  !< Buffer the site dependent internal field
   real(dblprec), dimension(:,:,:), allocatable :: larmb            !< Buffer the site dependent larmor frequencies
   real(dblprec), dimension(:,:,:,:), allocatable :: therm_fieldb   !< Buffer the site dependent stochastic field


contains


   !> Wrapper routine to print the thermal and effective fields
   !subroutine print_fields(mstep,sstep,Natom,Mensemble,simid,real_time_measure,delta_t,beff,thermal_field,emom)
   subroutine print_fields(mstep, sstep, Natom, Mensemble, simid, real_time_measure, delta_t, &
         beff, thermal_field, beff1, beff3, emom)

      implicit none

      integer, intent(in) :: mstep                      !< Current simulation step
      integer, intent(in) :: sstep                      !< Current simulation step in logarithmic scale
      integer, intent(in) :: Natom                      !< Number of atoms in the system
      integer, intent(in) :: Mensemble                  !< Number of ensembles
      character(len=8), intent(in) :: simid             !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      real(dblprec), intent(in) :: delta_t              !< Current time step
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff          !< Current site dependent total effective field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff1         !< Current site dependent internal effective field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff3         !< Current site dependent internal field from mixed spin-lattice Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field !< Current site dependent stochastic field
      real(dblprec), dimension(:,:,:), intent(in)    :: emom   !< Current unit moment vector

      ! Stochastic thermal field
      if (do_thermfield=='Y') then

         if (mod(sstep-1,thermfield_step)==0) then
            ! Write step to buffer
            call buffer_thermfield(Natom,Mensemble,mstep-1,thermal_field,&
               bcount_therm,delta_t,real_time_measure)

            if (bcount_therm==thermfield_buff) then
               ! write buffer to file
               call prn_thermfields(Natom,Mensemble,simid,real_time_measure)
               bcount_therm=1
            else
               bcount_therm=bcount_therm+1
            endif
         endif

      endif

      ! Total site dependent field
      if (do_prn_beff=='Y') then

         if (mod(sstep-1,beff_step)==0) then
            ! Write step to buffer
            call buffer_totalfield(Natom, Mensemble, mstep-1,beff,thermal_field,&
               bcount_beff,delta_t,real_time_measure)
            if (bcount_beff==beff_buff) then
               ! write buffer to file
               call prn_totalbfields(Natom, Mensemble, simid,real_time_measure)
               bcount_beff=1
            else
               bcount_beff=bcount_beff+1
            endif

         endif

      endif


      ! Total site dependent field
      if (do_prn_binteff=='Y') then

         if (mod(sstep-1,binteff_step)==0) then
            ! Write step to buffer
            call buffer_internalfield(Natom, Mensemble, mstep-1, beff1, beff3, &
               bcount_binteff, delta_t, real_time_measure)
            if (bcount_binteff==binteff_buff) then
               ! write buffer to file
               call prn_internalbfields(Natom, Mensemble, simid,real_time_measure)
               bcount_binteff=1
            else
               bcount_binteff=bcount_binteff+1
            endif

         endif

      endif


      ! Total site dependent field
      if (do_larmor_loc=='Y'.or.do_larmor_dos=='Y') then

         if (mod(sstep-1,larm_step)==0) then
            ! Write step to buffer
            call buffer_larmorfreq(Natom, Mensemble, mstep-1,beff,thermal_field,&
               bcount_larm,delta_t,real_time_measure,emom)
            if (bcount_larm==larm_buff) then
               ! write buffer to file
               call prn_larmorfreq(Natom, Mensemble, simid,real_time_measure)
               bcount_larm=1
            else
               bcount_larm=bcount_larm+1
            endif

         endif

      endif

      ! Site dependent torques
      if (do_prn_torques=='Y') then

         if (mod(sstep-1,torques_step)==0) then
            ! Write step to buffer
            call buffer_torques(Natom, Mensemble, mstep-1,beff,thermal_field,&
               bcount_torques, delta_t, real_time_measure, emom)
            if (bcount_torques==torques_buff) then
               ! write buffer to file
               call prn_torques(Natom, Mensemble, simid,real_time_measure)
               bcount_torques=1
            else
               bcount_torques=bcount_torques+1
            endif

         endif

      endif
      ! Site dependent spin torques
      if (do_prn_spin_torques=='Y') then

         if (mod(sstep-1,spin_torques_step)==0) then
            ! Write step to buffer
            call buffer_spin_torques(Natom, Mensemble, mstep-1,&
               bcount_spin_torques, delta_t, real_time_measure, emom)
            if (bcount_spin_torques==spin_torques_buff) then
               ! write buffer to file
               call prn_spin_torques(Natom, Mensemble, simid,real_time_measure)
               bcount_spin_torques=1
            else
               bcount_spin_torques=bcount_spin_torques+1
            endif

         endif

      endif

   end subroutine print_fields


   !> Flush the field measurements, i.e. print to file the fields in the last time step
   subroutine flush_prn_fields(Natom,Mensemble,simid,real_time_measure)

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid             !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      if (do_thermfield=='Y') then
         ! Write buffer to file
         bcount_therm=bcount_therm-1
         call prn_thermfields(Natom, Mensemble, simid,real_time_measure)
      endif

      if (do_prn_beff=='Y') then
         ! Write buffer to file
         bcount_beff=bcount_beff-1
         call prn_totalbfields(Natom, Mensemble, simid,real_time_measure)
      endif

      if (do_prn_binteff=='Y') then
         ! Write buffer to file
         bcount_binteff=bcount_binteff-1
         call prn_internalbfields(Natom, Mensemble, simid,real_time_measure)
      endif

      if (do_prn_torques=='Y') then
         ! Write buffer to file
         bcount_torques=bcount_torques-1
         call prn_torques(Natom, Mensemble, simid,real_time_measure)
      endif

      if (do_prn_spin_torques=='Y') then
         ! Write buffer to file
         bcount_spin_torques=bcount_spin_torques-1
         call prn_spin_torques(Natom, Mensemble, simid,real_time_measure)
      endif

   end subroutine flush_prn_fields


   !> Initialization of the variables for field printing with their default variables
   subroutine fields_prn_init()

      implicit none

      do_thermfield          = 'N'
      do_prn_beff            = 'N'
      do_prn_binteff         = 'N'
      do_larmor_loc          = 'N'
      do_larmor_dos          = 'N'
      do_prn_torques         = 'N'
      do_prn_spin_torques    = 'N'
      thermfield_step        = 1000
      thermfield_buff        = 10
      beff_step              = 1000
      beff_buff              = 10
      binteff_step           = 1000
      binteff_buff           = 10
      torques_step           = 1000
      torques_buff           = 10
      spin_torques_step      = 1000
      spin_torques_buff      = 10
      larm_step              = 1000
      larm_buff              = 10
      larm_dos_size          = 1000

   end subroutine fields_prn_init


   !> Routine for allocation and initialization of the arrays for field printing
   subroutine allocate_field_print(Natom,Mensemble,flag)

      implicit none

      integer, intent(in) :: flag      !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         bcount_beff=1
         bcount_binteff=1
         bcount_larm=1
         bcount_therm=1
         bcount_torques=1
         bcount_spin_torques=1

         if (do_thermfield=='Y') then
            allocate(therm_fieldb(3,Natom,thermfield_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(therm_fieldb))*kind(therm_fieldb),'therm_fieldb','allocate_measurements')
            allocate(indxb_therm(thermfield_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_therm))*kind(indxb_therm),'indxb_therm','allocate_measurements')
         endif

         if (do_prn_beff=='Y') then
            allocate(beffb(3,Natom,beff_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(beffb))*kind(beffb),'beffb','allocate_measurements')
            allocate(indxb_beff(beff_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_beff))*kind(indxb_beff),'indxb_beff','allocate_measurements')
         endif

         if (do_prn_binteff=='Y') then
            allocate(binteffb(6,Natom,binteff_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(binteffb))*kind(binteffb),'binteffb','allocate_measurements')
            allocate(indxb_binteff(binteff_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_binteff))*kind(indxb_binteff),'indxb_binteff','allocate_measurements')
         endif

         if (do_prn_torques=='Y') then
            allocate(torquesb(6,Natom,torques_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(torquesb))*kind(torquesb),'torquesb','allocate_measurements')
            allocate(indxb_torques(torques_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_torques))*kind(indxb_torques),'indxb_torques','allocate_measurements')
         endif

         if (do_prn_spin_torques=='Y') then
            allocate(spin_torquesb(6,Natom,spin_torques_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(spin_torquesb))*kind(spin_torquesb),'spin_torquesb','allocate_measurements')
            allocate(indxb_spin_torques(spin_torques_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_spin_torques))*kind(indxb_spin_torques),'indxb_spin_torques','allocate_measurements')
         endif

         if (do_larmor_loc=='Y'.or.do_larmor_dos=='Y') then
            allocate(larmb(Natom,larm_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(larmb))*kind(larmb),'larmb','allocate_measurements')
            allocate(indxb_larm(larm_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_larm))*kind(indxb_larm),'indxb_larm','allocate_measurements')
         endif

      else

         if (do_thermfield=='Y') then
            i_all=-product(shape(therm_fieldb))*kind(therm_fieldb)
            deallocate(therm_fieldb,stat=i_stat)
            call memocc(i_stat,i_all,'therm_fieldb','allocate_measurements')
            i_all=-product(shape(indxb_therm))*kind(indxb_therm)
            deallocate(indxb_therm,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_therm','allocate_measurements')
         endif

         if (do_prn_beff=='Y') then
            i_all=-product(shape(beffb))*kind(beffb)
            deallocate(beffb,stat=i_stat)
            call memocc(i_stat,i_all,'beffb','allocate_measurements')
            i_all=-product(shape(indxb_beff))*kind(indxb_beff)
            deallocate(indxb_beff,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_beff','allocate_measurements')
         endif

         if (do_prn_binteff=='Y') then
            i_all=-product(shape(binteffb))*kind(binteffb)
            deallocate(binteffb,stat=i_stat)
            call memocc(i_stat,i_all,'binteffb','allocate_measurements')
            i_all=-product(shape(indxb_binteff))*kind(indxb_binteff)
            deallocate(indxb_binteff,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_binteff','allocate_measurements')
         endif

         if (do_prn_torques=='Y') then
            i_all=-product(shape(torquesb))*kind(torquesb)
            deallocate(torquesb,stat=i_stat)
            call memocc(i_stat,i_all,'torquesb','allocate_measurements')
            i_all=-product(shape(indxb_torques))*kind(indxb_torques)
            deallocate(indxb_torques,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_torques','allocate_measurements')
         endif

         if (do_larmor_loc=='Y'.or.do_larmor_dos=='Y') then
            i_all=-product(shape(larmb))*kind(larmb)
            deallocate(larmb,stat=i_stat)
            call memocc(i_stat,i_all,'larmb','allocate_measurements')
            i_all=-product(shape(indxb_larm))*kind(indxb_larm)
            deallocate(indxb_larm,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_larm','allocate_measurements')
         endif

      endif

   end subroutine allocate_field_print


   !> Buffer site dependent thermal field
   subroutine buffer_thermfield(Natom, Mensemble, mstep,&
         thermal_field,bcount_therm,delta_t,real_time_measure)
      !

      implicit none

      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, intent(in) :: bcount_therm  !< Counter of buffer for thermal fields
      real(dblprec), intent(in) :: delta_t !< Current time step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: thermal_field   !< Current thermal field
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local scalar variables
      integer :: i,k

      do k=1, Mensemble
         do i=1, Natom
            therm_fieldb(1:3,i,bcount_therm,k)=thermal_field(1:3,i,k)
         end do
      end do

      if (real_time_measure=='Y') then
         indxb_therm(bcount_therm)=mstep*delta_t
      else
         indxb_therm(bcount_therm)=mstep
      endif

   end subroutine buffer_thermfield


   !> Buffer site dependent total field
   subroutine buffer_totalfield(Natom, Mensemble, mstep,beff,&
         thermal_field,bcount_beff,delta_t,real_time_measure)
      !

      implicit none

      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: bcount_beff   !< Counter of buffer for total effective field
      real(dblprec), intent(in) :: delta_t !< Current measurement time
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff          !< Current effective field from the hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: thermal_field !< Current thermal field
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local scalar variables
      integer :: i,k

      do k=1, Mensemble
         do i=1, Natom
            beffb(1:3,i,bcount_beff,k)=thermal_field(1:3,i,k)+beff(1:3,i,k)
         end do
      end do

      if (real_time_measure=='Y') then
         indxb_beff(bcount_beff)=mstep*delta_t
      else
         indxb_beff(bcount_beff)=mstep
      endif

   end subroutine buffer_totalfield


   !> Buffer site dependent internal field
   subroutine buffer_internalfield(Natom, Mensemble, mstep, beff1, beff3, &
         bcount_binteff, delta_t, real_time_measure)
      !

      implicit none

      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: bcount_binteff   !< Counter of buffer for total effective field
      real(dblprec), intent(in) :: delta_t !< Current measurement time
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff1  !< Current effective B-field from the magnetic Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff3  !< Current effective B-field from the mixed spin-lattice Hamiltonian
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local scalar variables
      integer :: i,k

      do k=1, Mensemble
         do i=1, Natom
            !beffb(1:3,i,bcount_beff,k)=thermal_field(1:3,i,k)+beff(1:3,i,k)
            binteffb(1:3,i,bcount_binteff,k) = beff1(1:3,i,k)
            binteffb(4:6,i,bcount_binteff,k) = beff3(1:3,i,k)
         end do
      end do

      if (real_time_measure=='Y') then
         indxb_binteff(bcount_binteff)=mstep*delta_t
      else
         indxb_binteff(bcount_binteff)=mstep
      endif

   end subroutine buffer_internalfield


   !> Buffer site dependent larmor frequencies
   subroutine buffer_larmorfreq(Natom, Mensemble, mstep,beff,&
         thermal_field,bcount_larm,delta_t,real_time_measure,emom)
      !

      implicit none

      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: bcount_larm   !< Counter of buffer for total effective field
      real(dblprec), intent(in) :: delta_t !< Current measurement time
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff          !< Current effective field from the hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: thermal_field !< Current thermal field
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      real(dblprec), dimension(:,:,:), intent(in)    :: emom   !< Current unit moment vector

      !.. Local scalar variables
      integer :: i,k
      real(dblprec), dimension(3) :: tmp_fld
      real(dblprec) :: fnorm, fdot

      do k=1, Mensemble
         do i=1, Natom
            tmp_fld=thermal_field(1:3,i,k)+beff(1:3,i,k)
            !Gramm-Schmidt orthogonalization
            fdot=tmp_fld(1)*emom(1,i,k)+tmp_fld(2)*emom(2,i,k)+tmp_fld(3)*emom(3,i,k)
            tmp_fld=tmp_fld-fdot*emom(1:3,i,k)
            !
            fnorm=sqrt(tmp_fld(1)*tmp_fld(1)+tmp_fld(2)*tmp_fld(2)+tmp_fld(3)*tmp_fld(3))
            larmb(i,bcount_larm,k)=fnorm
         end do
      end do

      if (real_time_measure=='Y') then
         indxb_larm(bcount_larm)=mstep*delta_t
      else
         indxb_larm(bcount_larm)=mstep
      endif

   end subroutine buffer_larmorfreq


   !> Buffer site spin transfer torques
   subroutine buffer_spin_torques(Natom, Mensemble, mstep,&
         bcount_spin_torques,delta_t,real_time_measure,emom)
      !
      use spintorques, only : btorque
      use math_functions, only : f_cross_product

      implicit none

      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: bcount_spin_torques   !< Counter of buffer for total effective field
      real(dblprec), intent(in) :: delta_t !< Current measurement time
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      real(dblprec), dimension(:,:,:), intent(in)    :: emom   !< Current unit moment vector

      !.. Local scalar variables
      integer :: i,k
      real(dblprec), dimension(3) :: prec_torque,damp_torque, tmp_fld

      ! Print precessional and damping torques from the spin transfer field
      do k=1, Mensemble
         do i=1, Natom
            tmp_fld=btorque(1:3,i,k)
            prec_torque = f_cross_product(emom(:,i,k), tmp_fld)
            damp_torque = f_cross_product(emom(:,i,k), prec_torque)
            spin_torquesb(1:3,i,bcount_spin_torques,k)=prec_torque
            spin_torquesb(4:6,i,bcount_spin_torques,k)=damp_torque
         end do
      end do

      if (real_time_measure=='Y') then
         indxb_spin_torques(bcount_spin_torques)=mstep*delta_t
      else
         indxb_spin_torques(bcount_spin_torques)=mstep
      endif

   end subroutine buffer_spin_torques

   !> Buffer site resulting torques
   subroutine buffer_torques(Natom, Mensemble, mstep,beff,thermal_field,&
         bcount_torques,delta_t,real_time_measure,emom)
      !
      use math_functions, only : f_cross_product

      implicit none

      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: bcount_torques   !< Counter of buffer for total effective field
      real(dblprec), intent(in) :: delta_t !< Current measurement time
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff          !< Current effective field from the hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: thermal_field !< Current thermal field
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      real(dblprec), dimension(:,:,:), intent(in)    :: emom   !< Current unit moment vector

      !.. Local scalar variables
      integer :: i,k
      real(dblprec), dimension(3) :: prec_torque,damp_torque, tmp_fld

      !!! ! Approach 1: Print precessional and damping torques from the total field (beff+btherm)
      !!! do k=1, Mensemble
      !!!    do i=1, Natom
      !!!       tmp_fld=thermal_field(1:3,i,k)+beff(1:3,i,k)
      !!!       prec_torque = f_cross_product(emom(:,i,k), tmp_fld)
      !!!       damp_torque = f_cross_product(emom(:,i,k), prec_torque)
      !!!       torquesb(1:3,i,bcount_torques,k)=prec_torque
      !!!       torquesb(4:6,i,bcount_torques,k)=damp_torque
      !!!    end do
      !!! end do
      ! Approach 2: Print total torques from thermal and effective fields respectively
      do k=1, Mensemble
         do i=1, Natom
            !Effective field
            tmp_fld=beff(1:3,i,k)
            prec_torque = f_cross_product(emom(:,i,k), tmp_fld)
            damp_torque = f_cross_product(emom(:,i,k), prec_torque)
            torquesb(1:3,i,bcount_torques,k)=prec_torque+damp_torque
            !Thermal field
            tmp_fld=thermal_field(1:3,i,k)
            prec_torque = f_cross_product(emom(:,i,k), tmp_fld)
            damp_torque = f_cross_product(emom(:,i,k), prec_torque)
            torquesb(4:6,i,bcount_torques,k)=prec_torque+damp_torque
         end do
      end do

      if (real_time_measure=='Y') then
         indxb_torques(bcount_torques)=mstep*delta_t
      else
         indxb_torques(bcount_torques)=mstep
      endif

   end subroutine buffer_torques


   !> Print total effective field
   subroutine prn_larmorfreq(Natom, Mensemble, simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of the simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: i,j,k, larm_idx
      character(len=30) :: filn
      real(dblprec), dimension(:), allocatable :: larm_dos       !< Histogram array for Larmor frequencies
      real(dblprec) :: larm_max

      ! Print thermal fields to output file if specified
      ! Remember to remove old data since the write statement appends new data to the file

      if(do_larmor_loc=='Y') then
         write(filn,'(''larmor_loc.'',a,''.out'')') trim(simid)
         open(ofileno,file=filn, position = 'APPEND',form = 'formatted')
         do k=1, bcount_larm
            do j=1,Mensemble
               do i=1,Natom
                  if (real_time_measure=='Y') then
                     write (ofileno,121) indxb_larm(k), i, j, larmb(i,k,j)
                  else
                     write (ofileno,120) int(indxb_larm(k)), i, j, larmb(i,k,j)
                  endif
               end do
            end do
         enddo
         close(ofileno)
      end if
      if(do_larmor_dos=='Y') then
         write(filn,'(''larmor_dos.'',a,''.out'')') trim(simid)
         open(ofileno,file=filn, position = 'APPEND',form = 'formatted')
         allocate(larm_dos(0:larm_dos_size))
         larm_dos=0.0_dblprec
         larm_max=maxval(larmb)*1.005_dblprec
         do k=1, bcount_larm
            do j=1,Mensemble
               do i=1,Natom
                  larm_idx=int(larm_dos_size*larmb(i,k,j)/larm_max)
                  larm_dos(larm_idx)=larm_dos(larm_idx)+1.0_dblprec
               end do
            end do
         enddo
         larm_dos=larm_dos/(1.0_dblprec*k*j*i)
         do larm_idx=1,larm_dos_size
            write (ofileno,122) larm_idx,hbar_mev*gama*larm_max/(1.0_dblprec*larm_dos_size)*larm_idx,larm_dos(larm_idx)
         end do
         close(ofileno)
      end if
      return

      write(*,*) 'Error writing the Larmor frequency files'

      120 format(i8,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)
      121 format(es16.4,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)
      122 format(i8,es16.8,es16.8)

   end subroutine prn_larmorfreq


   !> Print total effective field
   subroutine prn_totalbfields(Natom, Mensemble, simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of the simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: i,j,k
      character(len=30) :: filn
      real(dblprec) :: temp_norm

      ! Print thermal fields to output file if specified
      ! Remember to remove old data since the write statement appends new data to the file

      write(filn,'(''befftot.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position = 'APPEND',form = 'formatted')

      ! Write header to output files for first iteration
      if(abs(indxb_beff (1))<=0.0e0_dblprec) then
         write (ofileno,'(a)') "# Iter.     Site     Replica      B_x             B_y             B_z             B"
      end if

      do k=1, bcount_beff
         do j=1,Mensemble
            do i=1,Natom
               temp_norm=sqrt(beffb(1,i,k,j)**2+beffb(2,i,k,j)**2+beffb(3,i,k,j)**2)
               if (real_time_measure=='Y') then
                  write (ofileno,121) indxb_beff(k), i, j, beffb(1:3,i,k,j),temp_norm
               else
                  write (ofileno,120) int(indxb_beff(k)), i, j, beffb(1:3,i,k,j),temp_norm
               endif
            end do
         end do
      enddo
      close(ofileno)

      return

      write(*,*) 'Error writing the total field file'

      120 format(i8,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)
      121 format(es16.4,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)

   end subroutine prn_totalbfields


   !> Print internal magnetic field effective field from spin Hamiltonian
   !> and from spin-lattice Hamiltonian respectively
   subroutine prn_internalbfields(Natom, Mensemble, simid, real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of the simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: i,j,k
      character(len=30) :: filn

      ! Print thermal fields to output file if specified
      ! Remember to remove old data since the write statement appends new data to the file

      write(filn,'(''bintefftot.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position = 'APPEND',form = 'formatted')

      ! Write header to output files for first iteration
      if(abs(indxb_binteff (1))<=0.0e0_dblprec) then
         write (ofileno,'(a)') "# Iter.     Site     Replica      B_SD_x      B_SD_y      B_SD_z      B_SD       &
            & B_SLD_x     B_SLD_y     B_SLD_z     B_SLD"
      end if

      do k=1, bcount_binteff
         do j=1,Mensemble
            do i=1,Natom
               if (real_time_measure=='Y') then
                  write (ofileno,121) indxb_binteff(k), i, j, binteffb(1:3,i,k,j), norm2(binteffb(1:3,i,k,j)), &
                     binteffb(4:6,i,k,j), norm2(binteffb(4:6,i,k,j))
               else
                  write (ofileno,120) int(indxb_binteff(k)), i, j, binteffb(1:3,i,k,j), norm2(binteffb(1:3,i,k,j)), &
                     binteffb(4:6,i,k,j), norm2(binteffb(4:6,i,k,j))
               endif
            end do
         end do
      enddo
      close(ofileno)

      return

      write(*,*) 'Error writing the internal field file'

      120 format(i8,i8,i8,8x,8es12.4)
      121 format(es12.4,i8,i8,8x,8es12.4)

   end subroutine prn_internalbfields

   !> Print magnetic torques
   subroutine prn_torques(Natom, Mensemble, simid, real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of the simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: i,j,k
      character(len=30) :: filn

      ! Print thermal fields to output file if specified
      ! Remember to remove old data since the write statement appends new data to the file
      write(filn,'(''torques.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position = 'APPEND',form = 'formatted')

      ! Write header to output files for first iteration
      if(abs(indxb_torques(1))<=0.0e0_dblprec) then
         write (ofileno,'(a)') "# Iter.     Site     Replica        Tau_B_x     Tau_B_y     Tau_B_z     |Tau_B|     &
          Tau_b_x     Tau_b_y     Tau_b_z     |Tau_b|"
      end if

      do k=1, bcount_torques
         do j=1,Mensemble
            do i=1,Natom
               if (real_time_measure=='Y') then
                  write (ofileno,121) indxb_torques(k), i, j, torquesb(1:3,i,k,j), norm2(torquesb(1:3,i,k,j)), &
                     torquesb(4:6,i,k,j), norm2(torquesb(4:6,i,k,j))
               else
                  write (ofileno,120) int(indxb_torques(k)), i, j, torquesb(1:3,i,k,j), norm2(torquesb(1:3,i,k,j)), &
                     torquesb(4:6,i,k,j), norm2(torquesb(4:6,i,k,j))
               endif
            end do
         end do
      enddo
      close(ofileno)

      return

      write(*,*) 'Error writing the internal field file'

      120 format(i8,i8,i8,8x,8es12.4)
      121 format(es12.4,i8,i8,8x,8es12.4)

   end subroutine prn_torques

   !> Print spin transfer torques
   subroutine prn_spin_torques(Natom, Mensemble, simid, real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of the simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: i,j,k
      character(len=30) :: filn

      ! Print thermal fields to output file if specified
      ! Remember to remove old data since the write statement appends new data to the file
      write(filn,'(''spin_torques.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position = 'APPEND',form = 'formatted')

      ! Write header to output files for first iteration
      if(abs(indxb_spin_torques(1))<=0.0e0_dblprec) then
         write (ofileno,'(a)') "# Iter.     Site     Replica        STT_p_x     STT_p_y     STT_p_z     |STT_p|     &
            STT_d_x     STT_d_y     STT_d_z     |STT_d|"
      end if

      do k=1, bcount_spin_torques
         do j=1,Mensemble
            do i=1,Natom
               if (real_time_measure=='Y') then
                  write (ofileno,121) indxb_spin_torques(k), i, j, &
                  spin_torquesb(1:3,i,k,j), norm2(spin_torquesb(1:3,i,k,j)), &
                     spin_torquesb(4:6,i,k,j), norm2(spin_torquesb(4:6,i,k,j))
               else
                  write (ofileno,120) int(indxb_spin_torques(k)), i, j, &
                  spin_torquesb(1:3,i,k,j), norm2(spin_torquesb(1:3,i,k,j)), &
                     spin_torquesb(4:6,i,k,j), norm2(spin_torquesb(4:6,i,k,j))
               endif
            end do
         end do
      enddo
      close(ofileno)

      return

      write(*,*) 'Error writing the internal field file'

      120 format(i8,i8,i8,8x,8es12.4)
      121 format(es12.4,i8,i8,8x,8es12.4)

   end subroutine prn_spin_torques



   !> Print thermal field
   subroutine prn_thermfields(Natom, Mensemble, simid,real_time_measure)
      !

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Simulation name
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      character(len=30) :: filn
      ! Local variables
      integer :: i,j,k

      ! Print thermal fields to output file if specified
      ! Remember to remove old data since the write statement appends new data to the file
      write(filn,'(''thermfields.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position = 'APPEND')
      do k=1, bcount_therm
         do j=1,Mensemble
            do i=1,Natom
               if (real_time_measure=='Y') then
                  write (ofileno,111) indxb_therm(k), i, j, therm_fieldb(:,i,k,j),norm2(therm_fieldb(:,i,k,j))
               else
                  write (ofileno,110) int(indxb_therm(k)), i, j, therm_fieldb(:,i,k,j),norm2(therm_fieldb(:,i,k,j))
               endif
            end do
         end do
      enddo
      close(ofileno)

      return

      write(*,*) 'Error writting the thermal field file'

      110 format(i8,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)
      111 format(es16.4,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)

   end subroutine prn_thermfields


end module prn_fields
