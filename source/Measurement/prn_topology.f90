!------------------------------------------------------------------------------------
!> @brief
!> Data structures needed for the skyrmion number claculation and its printing
!
!> @author
!> Anders Bergman
!> Jonathan Chico ---> Reorganized in different modules, added site and type dependence
!> @copyright
!> GNU Public License
!------------------------------------------------------------------------------------
module prn_topology

   use Profiling
   use Parameters
   use Topology

   implicit none

   ! Counters for the measurements
   integer :: bcount_skyno !< Counter of buffer for skyrmion number
   integer, dimension(:), allocatable :: indxb_skyno          !< Step counter for the skyrmion number
   real(dblprec), dimension(:), allocatable :: skynob         !< Buffer for the skyrmion number
   real(dblprec), dimension(:,:), allocatable :: proj_skynob  !< Buffer for the skyrmion number
   real(dblprec), dimension(:,:), allocatable :: dens_skynob  !< Buffer fot the skyrmion density
   real(dblprec), dimension(:,:), allocatable :: cmass_skynob  !< Buffer fot the skyrmion density
   real(dblprec), dimension(:,:,:), allocatable :: spol_buff  !< Buffer data for polarization vector
   real(dblprec), dimension(:,:,:,:), allocatable :: grad_mom !< Gradient of moments
   real(dblprec), dimension(:,:,:,:,:), allocatable :: proj_grad_mom !< Gradient of moments
   real(dblprec) :: sk_num_cum=0.0_dblprec !< Cumulative sum of the skyrmion number
   real(dblprec) :: sk2_num_cum=0.0_dblprec !< Cumulative sum of the skyrmion number squared
   real(dblprec) :: sk_avrg=0.0_dblprec !< Cumulative average of the skyrmion number
   real(dblprec) :: sk2_avrg=0.0_dblprec !< Cumulative average of the skyrmion number squared
   real(dblprec) :: sk_var=0.0_dblprec !< Variance of the skyrmion number 
   integer :: sk_num_num=0 !< Number of cumulative samples of the skyrmion number

   private

   public :: print_topology, flush_topology, allocate_topology,prn_topology_init
   public :: skyno, skyno_step, skyno_buff, do_proj_skyno, do_skyno_den, do_skyno_cmass
   public :: sk_avrg, sk_var

contains

   !---------------------------------------------------------------------------------
   !> @brief
   !> Wrapper for printing the topological information of the sample
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------------
   subroutine print_topology(NT,NA,N1,N2,sstep,mstep,Natom,Mensemble,delta_t,          &
         real_time_measure,emomM,emom,simid,atype)

      implicit none

      integer, intent(in) :: NT            !< Number of types of atoms
      integer, intent(in) :: NA            !< Number atoms in unit cell
      integer, intent(in) :: N1            !< Number of cell repetitions in x direction
      integer, intent(in) :: N2            !< Number of cell repetitions in y direction
      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: sstep         !< Current simulation step in logarithmic scale
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      real(dblprec), intent(in) :: delta_t !< Time step for real time measurement
      character(len=8), intent(in) :: simid                !< Name of simulation
      character(len=1), intent(in) :: real_time_measure    !< Measurements displayed in real time
      integer, dimension(Natom), intent(in) :: atype !< Type of ato
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom   !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector

      ! Skyrmion Number
      if (mod(sstep-1,skyno_step)==0) then
         if (skyno=='Y') then

            ! Wite the total skyrmion number to the buffer
            call buffer_skyno(Natom, Mensemble,N1,N2,mstep-1,emomM,emom,            &
               bcount_skyno,delta_t,real_time_measure)

            ! Write the type projected skyrmion number to the buffer
            if ( do_proj_skyno=='Y') then
               call buffer_proj_skyno(NT,Natom, Mensemble,N1,N2,mstep-1,emomM,emom,&
                  bcount_skyno,delta_t,real_time_measure,atype)
            endif
            ! Write the site dependent skyrmion number to the buffer
            if (do_skyno_den=='Y') then
               call buffer_dens_skyno(Natom, Mensemble,N1,N2,mstep-1,emomM,emom,&
                  bcount_skyno,delta_t,real_time_measure)
            endif

         else if (skyno=='T') then
            ! Wite the total skyrmion number to the buffer
            call buffer_skyno_tri(NA, Natom, Mensemble,N1,N2,mstep-1,emom,            &
               bcount_skyno,delta_t,real_time_measure)
            ! Write the type projected skyrmion number to the buffer
            if ( do_proj_skyno=='Y' .or. do_proj_skyno=='T') then
               call buffer_proj_skyno(NA,Natom, Mensemble,N1,N2,mstep-1,emomM,emom,&
                  bcount_skyno,delta_t,real_time_measure,atype)
            endif
            ! Write the site dependent skyrmion number to the buffer
            if (do_skyno_den=='Y') then
               call buffer_dens_skyno(Natom, Mensemble,N1,N2,mstep-1,emomM,emom,&
                  bcount_skyno,delta_t,real_time_measure)
            endif

         end if


         if (skyno=='Y'.or.skyno=='T') then
            if (bcount_skyno==skyno_buff) then
               ! Print the buffer total skyrmion number to file
               call prn_skyno(simid,real_time_measure)
               ! Print the type projected skyrmion number to file
               if ( do_proj_skyno=='Y') then
                  call prn_proj_skyno(NT,simid,real_time_measure)
               endif
               ! Print the  buffer skyrion density to file
               if (do_skyno_den=='Y') then
                  call prn_dens_skyno(simid,Natom,real_time_measure)
                  if (do_skyno_cmass=='Y') then
                     call prn_cmass_skyno(simid,Natom,real_time_measure)
                  endif
               endif

               bcount_skyno=1
            else
               bcount_skyno=bcount_skyno+1
            endif
         end if

      end if

   end subroutine print_topology

   !---------------------------------------------------------------------------------
   !> @brief
   !> Allocate topological data and initialization of variables
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------------
   subroutine allocate_topology(NT,Natom,Mensemble,flag)

      implicit none

      integer, intent(in) :: NT         !< Atom types in the system
      integer, intent(in) :: flag       !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Natom      !< Number of atoms in system
      integer, intent(in) :: Mensemble  !< Number of ensembles

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         bcount_skyno=1
         if (skyno=='Y'.or.skyno=='T') then
            allocate(skynob(skyno_buff),stat=i_stat)
            call memocc(i_stat,product(shape(skynob))*kind(skynob),'skynob','allocate_topology')
            skynob=0.0_dblprec
            allocate(indxb_skyno(skyno_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_skyno))*kind(indxb_skyno),'indxb_skyno','allocate_topology')
            indxb_skyno=0
            allocate(grad_mom(3,3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(grad_mom))*kind(grad_mom),'grad_mom','allocate_topology')
            grad_mom=0.0_dblprec

            if (do_proj_skyno=='Y') then
               allocate(proj_skynob(skyno_buff,NT),stat=i_stat)
               call memocc(i_stat,product(shape(proj_skynob))*kind(proj_skynob),'proj_skynob','allocate_topology')
               proj_skynob=0.0_dblprec
               allocate(proj_grad_mom(3,3,Natom,Mensemble,NT),stat=i_stat)
               call memocc(i_stat,product(shape(proj_grad_mom))*kind(proj_grad_mom),'proj_grad_mom','allocate_topology')
               proj_grad_mom=0.0_dblprec
            endif

            if (do_skyno_den=='Y') then
               allocate(dens_skynob(skyno_buff,Natom),stat=i_stat)
               call memocc(i_stat,product(shape(dens_skynob))*kind(dens_skynob),'dens_skynob','allocate_topology')
               dens_skynob=0.0_dblprec
               if (do_skyno_cmass=='Y') then
                  allocate(cmass_skynob(3,skyno_buff),stat=i_stat)
                  call memocc(i_stat,product(shape(cmass_skynob))*kind(cmass_skynob),'cmass_skynob','allocate_topology')
                  cmass_skynob=0.0_dblprec
               endif
            endif

         endif

      else

         if (skyno=='Y'.or.skyno=='T') then
            i_all=-product(shape(skynob))*kind(skynob)
            deallocate(skynob,stat=i_stat)
            call memocc(i_stat,i_all,'skynob','allocate_topology')
            i_all=-product(shape(indxb_skyno))*kind(indxb_skyno)
            deallocate(indxb_skyno,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_skyno','allocate_topology')
            i_all=-product(shape(grad_mom))*kind(grad_mom)
            deallocate(grad_mom,stat=i_stat)
            call memocc(i_stat,i_all,'grad_mom','allocate_topology')
            ! deallocate variables for the projected skyrmion density
            if (do_proj_skyno=='Y') then
               i_all=-product(shape(proj_skynob))*kind(proj_skynob)
               deallocate(proj_skynob,stat=i_stat)
               call memocc(i_stat,i_all,'proj_skynob','allocate_topology')
               i_all=-product(shape(proj_grad_mom))*kind(proj_grad_mom)
               deallocate(proj_grad_mom,stat=i_stat)
               call memocc(i_stat,i_all,'proj_grad_mom','allocate_topology')
            endif
            ! deallocate variables for the site dependent skyrmion density
            if (do_skyno_den=='Y') then
               i_all=-product(shape(dens_skynob))*kind(dens_skynob)
               deallocate(dens_skynob,stat=i_stat)
               call memocc(i_stat,i_all,'dens_skynob','allocate_topology')
               if (do_skyno_cmass=='Y') then
                  i_all=-product(shape(cmass_skynob))*kind(cmass_skynob)
                  deallocate(cmass_skynob,stat=i_stat)
                  call memocc(i_stat,i_all,'cmass_skynob','allocate_topology')
               endif
            endif

         endif

      end if

   end subroutine allocate_topology

   !---------------------------------------------------------------------------------
   !> @brief
   !> Initialization of the variables for topological measurements
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------------
   subroutine prn_topology_init()

      implicit none
      ! Topological functions
      skyno         = "N"
      do_proj_skyno = "N"
      do_skyno_den  = "N"
      skyno_step    = 100
      skyno_buff    = 10

   end subroutine prn_topology_init

   !---------------------------------------------------------------------------------
   !> @brief
   !> Subroutine to flush the topological information, i.e. print to file
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------------
   subroutine flush_topology(NT,Natom,simid,real_time_measure)

      implicit none

      integer, intent(in) :: NT !< Number of types of atoms
      integer, intent(in) :: Natom !< Number of atoms in the system
      character(len=8), intent(in) :: simid              !< Name of simulation
      character(len=1), intent(in) :: real_time_measure  !< Measurements displayed in real time

      ! Skyrmion number
      if (skyno=='Y'.or.skyno=='T') then
         bcount_skyno=bcount_skyno-1
         ! Write the total skyrmion number buffer to file
         call prn_skyno(simid,real_time_measure)

         if (do_proj_skyno=='Y') then
            ! Write the type dependent buffer to file
            call prn_proj_skyno(NT,simid,real_time_measure)
         endif

         ! Write the site dependent buffer to file
         if (do_skyno_den=='Y') then
            call prn_dens_skyno(simid,Natom,real_time_measure)
         endif
         bcount_skyno=1
      endif

   end subroutine flush_topology

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing the type dependent projected skyrmion number
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------------
   subroutine prn_skyno(simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Simulation name
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: k
      character(len=30) :: filn
      real(dblprec) :: sk_avrg_prev

      ! write the skyrmion number
      write(filn,'(''sknumber.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position="append")

      ! Write header to output files for first iteration
      if(indxb_skyno(1)==0) then
            if (real_time_measure=='Y') then
               write (ofileno,247) "Time[s]","Skx num", "Skx avg", "Skx std"
            else
               write (ofileno,246) "# Iter","Skx num", "Skx avg", "Skx std"
            endif
      end if

      do k=1, bcount_skyno
         sk_num_cum = sk_num_cum + skynob(k)
         sk_num_num = sk_num_num + 1
         sk_avrg_prev = sk_avrg
         sk_avrg = sk_avrg_prev + (skynob(k)-sk_avrg_prev)/sk_num_num
         sk_var = sk_var + (skynob(k) - sk_avrg_prev) * (skynob(k) - sk_avrg)

         if (real_time_measure=='Y') then
            write (ofileno,241) indxb_skyno(k), skynob(k), sk_avrg, sk_var/sk_num_num
            !write (ofileno,241) indxb_skyno(k), skynob(k), sk_num_cum/sk_num_num
         else
            write (ofileno,240) int(indxb_skyno(k)), skynob(k), sk_avrg, sk_var/sk_num_num
            !write (ofileno,240) int(indxb_skyno(k)), skynob(k), sk_num_cum/sk_num_num
         endif
      enddo
      !        mean_old = mean
      !        mean = mean_old + (w / w_sum) * (x - mean_old)


      close(ofileno)

      return

      write(*,*) "Error writing the skyrmion number file"

      240 format(i8,2x,5f16.8)
      241 format(es16.6,2x,5f16.8)
      246 format(a8,2x,a10,2x,a10,2x,a10)
      247 format(a16,2x,a10,2x,a10,2x,a10)

   end subroutine prn_skyno

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing the type dependent projected skyrmion number
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine prn_proj_skyno(NT,simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NT !< Number of types of atoms
      character(len=8), intent(in) :: simid !< Simulation name
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: k,ii
      character(len=30) :: filn

      ! write the projected skyrmion number
      do ii=1, NT

         if (ii>10) then
            write(filn,'(''proj_sknumber.'',i2,''.'',a,''.out'')') ii,trim(simid)
         else
            write(filn,'(''proj_sknumber.'',i1,''.'',a,''.out'')') ii,trim(simid)

         endif
         open(ofileno,file=filn, position="append")

         ! Write header to output files for first iteration
         if(indxb_skyno(1)==0) then
            if (real_time_measure=='Y') then
               write (ofileno,247) "Time[s]","Skx num"
            else
               write (ofileno,246) "# Iter","Skx num"
            endif
         end if

         do k=1, bcount_skyno
            if (real_time_measure=='Y') then
               write (ofileno,241) indxb_skyno(k), proj_skynob(k,ii)
            else
               write (ofileno,240) int(indxb_skyno(k)), proj_skynob(k,ii)
            endif
         enddo

         close(ofileno)

      enddo

      return

      write(*,*) "Error writing the skyrmion number file"

      240 format(i8,2x,f10.4)
      241 format(es16.4,2x,f10.4)
      246 format(a8,2x,a10)
      247 format(a16,2x,a10)
   end subroutine prn_proj_skyno

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing the center of mass of the total topological charge
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine prn_cmass_skyno(simid,Natom,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in the sample
      character(len=8), intent(in) :: simid !< Simulation name
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: k,iatom
      character(len=30) :: filn

      ! write the site dependent skyrmion number
      write(filn,'(''cmass_skynum.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position="append")

      ! Write header to output files for first iteration
      if(indxb_skyno(1)==0) then
         if (real_time_measure=='Y') then
            write (ofileno,247) "Time[s]","rSk_x","rSk_y","rSk_z"
         else
            write (ofileno,246) "# Iter", "rSk_x","rSk_y","rSk_z"
         endif
      end if

      do k=1, bcount_skyno
         if (real_time_measure=='Y') then
               write (ofileno,245) indxb_skyno(k), cmass_skynob(:,k)
         else
               write (ofileno,244) int(indxb_skyno(k)), cmass_skynob(:,k)
         endif
      enddo

      close(ofileno)

      return

      write(*,*) "Error writing the skyrmion number center of mass file"

      244 format(i8,2x,3f16.8)
      245 format(es16.4,2x,3f16.8)
      246 format(a8,2x,3a10)
      247 format(a16,2x,3a10)

   end subroutine prn_cmass_skyno

   !---------------------------------------------------------------------------------
   !> @brief
   !> Printing the site dependent skyrmion number
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine prn_dens_skyno(simid,Natom,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in the sample
      character(len=8), intent(in) :: simid !< Simulation name
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: k,iatom
      character(len=30) :: filn

      ! write the site dependent skyrmion number
      write(filn,'(''dens_skynum.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position="append")

      ! Write header to output files for first iteration
      if(indxb_skyno(1)==0) then
         if (real_time_measure=='Y') then
            write (ofileno,247) "Time[s]","Site","Skx num"
         else
            write (ofileno,246) "# Iter", "Site","Skx num"
         endif
      end if

      do k=1, bcount_skyno
         if (real_time_measure=='Y') then
            do iatom=1,Natom
               write (ofileno,245) indxb_skyno(k),iatom, dens_skynob(k,iatom)
            enddo
         else
            do iatom=1,Natom
               write (ofileno,244) int(indxb_skyno(k)),iatom, dens_skynob(k,iatom)
            enddo
         endif
      enddo

      close(ofileno)

      return

      write(*,*) "Error writing the site dependent skyrmion number file"

      244 format(i8,2x,i8,2x,e16.8)
      246 format(a8,2x,a,2x,a10)
      245 format(es16.4,2x,i8,2x,e16.8)
      247 format(a16,2x,a,2x,a10)

   end subroutine prn_dens_skyno

   !---------------------------------------------------------------------------------
        ! SUBROUTINE buffer_skyno  
        !> Buffer the skyrmion number
   !---------------------------------------------------------------------------------
        subroutine buffer_skyno(Natom, Mensemble,N1,N2,mstep,emomM,emom,bcount_skyno,           &
                delta_t,real_time_measure)
      !
      use Gradients
      use Topology

      implicit none

      integer, intent(in) :: N1    !< Number of cell repetitions in x direction
      integer, intent(in) :: N2    !< Number of cell repetitions in y direction
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, intent(in) :: bcount_skyno  !< Counter of buffer for skyrmion number
      real(dblprec), intent(in) :: delta_t !< Current time step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current moment vector
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      call grad_moments(Natom, Mensemble,emom, grad_mom)
      skynob(bcount_skyno)=pontryagin_no(Natom, Mensemble,emomM, grad_mom)

      if (real_time_measure=='Y') then
         indxb_skyno(bcount_skyno)=nint(real(mstep,dblprec)*delta_t)
      else
         indxb_skyno(bcount_skyno)=mstep
      endif

   end subroutine buffer_skyno

   !---------------------------------------------------------------------------------
   ! SUBROUTINE buffer_proj_skyno  
   !> Buffer the projected skyrmion number
   !---------------------------------------------------------------------------------
   subroutine buffer_proj_skyno(NT,Natom, Mensemble,N1,N2,mstep,emomM,emom,                     &
         bcount_skyno,delta_t,real_time_measure,atype)
      !
      use Gradients
      use Topology

      implicit none

      integer, intent(in) :: NT    !< Number of types of atoms (or NA depending on call)
      integer, intent(in) :: N1    !< Number of cell repetitions in x direction
      integer, intent(in) :: N2    !< Number of cell repetitions in y direction
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, intent(in) :: bcount_skyno  !< Counter of buffer for skyrmion number
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      real(dblprec), intent(in) :: delta_t !< Current time step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      if (skyno == 'T') then
         proj_skynob(bcount_skyno,1:NT)=pontryagin_tri_proj(NT, Natom,Mensemble,emom)
      else
         call proj_grad_moments(NT,Natom, Mensemble,atype,emomM,proj_grad_mom)
         proj_skynob(bcount_skyno,1:NT)=proj_pontryagin_no(NT,Natom,Mensemble,atype,emomM,proj_grad_mom)
      end if

      if (real_time_measure=='Y') then
         indxb_skyno(bcount_skyno)=nint(real(mstep,dblprec)*delta_t)
      else
         indxb_skyno(bcount_skyno)=mstep
      endif

   end subroutine buffer_proj_skyno

   !---------------------------------------------------------------------------------
   !> @brief
   !> Buffer the site dependent skyrion number
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine buffer_dens_skyno(Natom, Mensemble,N1,N2,mstep,emomM,emom,                                &
         bcount_skyno,delta_t,real_time_measure)
      !
      use Gradients
      use Topology
      use SystemData, only : coord

      implicit none

      integer, intent(in) :: N1    !< Number of cell repetitions in x direction
      integer, intent(in) :: N2    !< Number of cell repetitions in y direction
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, intent(in) :: bcount_skyno  !< Counter of buffer for skyrmion number
      real(dblprec), intent(in) :: delta_t !< Current time step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current moment vector
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! .. Local variables
      integer :: iatom
      !real(dblprec), dimension(3) :: skmass
      !real(dblprec) :: sktot

      ! Call the calculation of the magnetization gradient
      if (skyno=='Y') then
         ! Is this necesary? it is already called for the total skyrmion number
         call grad_moments(Natom, Mensemble,emom, grad_mom)

         !skmass=0.0_dblprec
         !sktot=0.0_dblprec
         ! Loop over the atomic sites
         do iatom=1,Natom
            dens_skynob(bcount_skyno,iatom)=pontryagin_no_density(iatom,Natom, Mensemble,emomM, grad_mom)
            !sktot=sktot+pontryagin_no_density(iatom,Natom, Mensemble,emomM, grad_mom)
            !skmass(:)=skmass(:) + pontryagin_no_density(iatom,Natom, Mensemble,emomM, grad_mom)*coord(:,iatom)
         enddo
         !print '(a,3f12.6)' , '  SkMassC: ',skmass/sktot
      else if (skyno=='T') then
         do iatom=1,Natom
            dens_skynob(bcount_skyno,iatom)=pontryagin_tri_dens(iatom,Natom, Mensemble,emom)
         enddo
      end if

      if (real_time_measure=='Y') then
         indxb_skyno(bcount_skyno)=nint(real(mstep,dblprec)*delta_t)
      else
         indxb_skyno(bcount_skyno)=mstep
      endif


      if (do_skyno_cmass=='Y') then
         call buffer_cmass(Natom, Mensemble, coord, emom, dens_skynob(bcount_skyno,:),bcount_skyno)
      end if
   end subroutine buffer_dens_skyno

   !---------------------------------------------------------------------------------
   ! SUBROUTINE buffer_skyno_tri
   !> Buffer the skyrmion number using triangular mesh
   !---------------------------------------------------------------------------------
   subroutine buffer_skyno_tri(NA, Natom, Mensemble,N1,N2,mstep,emom,bcount_skyno,           &
                delta_t,real_time_measure)
      !
      use Gradients
      use Topology

      implicit none

      integer, intent(in) :: NA    !< Number of atoms in unit cell
      integer, intent(in) :: N1    !< Number of cell repetitions in x direction
      integer, intent(in) :: N2    !< Number of cell repetitions in y direction
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, intent(in) :: bcount_skyno  !< Counter of buffer for skyrmion number
      real(dblprec), intent(in) :: delta_t !< Current time step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      skynob(bcount_skyno)=pontryagin_tri(Natom, Mensemble,emom) / NA

      if (real_time_measure=='Y') then
         indxb_skyno(bcount_skyno)=nint(real(mstep,dblprec)*delta_t)
      else
         indxb_skyno(bcount_skyno)=mstep
      endif

   end subroutine buffer_skyno_tri

   subroutine buffer_cmass(Natom, Mensemble, coord, mom_arr, topo_arr,bidx)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: mom_arr
      real(dblprec), dimension(Natom), intent(in) :: topo_arr
      integer, intent(in) :: bidx  !< Index in buffer array
      !
      integer :: ia, ir, ik
      real(dblprec), dimension(3) :: mx, my, mz, tp
      real(dblprec) :: eps=1e-6_dblprec
      !
      mx = 0.0_dblprec
      my = 0.0_dblprec
      mz = 0.0_dblprec
      tp = 0.0_dblprec
      do ik=1, Mensemble
         do ia=1, Natom
            mx = mx + coord(:,ia) * mom_arr(1,ia,ik) 
            my = my + coord(:,ia) * mom_arr(2,ia,ik) 
            mz = mz + coord(:,ia) * mom_arr(3,ia,ik) 
            tp = tp + coord(:,ia) * topo_arr(ia)
         end do
      end do
      mx = mx/(sum(mom_arr(1,:,:))+eps)
      my = my/(sum(mom_arr(2,:,:))+eps)
      mz = mz/(sum(mom_arr(3,:,:))+eps)
      !!! tp = tp/(sum(topo_arr)+eps)
      cmass_skynob(:,bidx) = tp/(sum(topo_arr)+eps)
      !
      !!! print '(a, 3f12.6,f14.8)', 'Topological center of mass:', tp,(sum(topo_arr))
      !!! print '(a, 3g20.6,i)', 'Topological center of mass:', cmass_skynob(:,bidx), bidx
   end subroutine buffer_cmass



end module prn_topology
