!====================================================================!
!> @brief
!> Data structures needed for the skyrmion number claculation and its printing
!
!> @author
!> Anders Bergman
!> Jonathan Chico ---> Reorganized in different modules, added site and type dependance
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!====================================================================!
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
   real(dblprec), dimension(:,:,:), allocatable :: spol_buff  !< Buffer data for polarization vector
   real(dblprec), dimension(:,:,:,:), allocatable :: grad_mom !< Gradient of moments
   real(dblprec), dimension(:,:,:,:,:), allocatable :: proj_grad_mom !< Gradient of moments

   private

   public :: print_topology, flush_topology, allocate_topology,prn_topology_init
   public :: skyno, skyno_step, skyno_buff, do_proj_skyno, do_skyno_den

contains

   !---------------------------------------------------------------------------
   !> @brief
   !> Wrapper for printing the topological information of the sample
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------
   subroutine print_topology(NT,N1,N2,sstep,mstep,Natom,Mensemble,delta_t,&
         real_time_measure,emomM,emom,simid,atype)

      implicit none

      integer, intent(in) :: NT            !< Number of types of atoms
      integer, intent(in) :: N1            !< Number of cell repetitions in x direction
      integer, intent(in) :: N2            !< Number of cell repetitions in y direction
      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: sstep         !< Current simulation step in logarithmic scale
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, dimension(Natom), intent(in) :: atype !< Type of ato
      real(dblprec), intent(in) :: delta_t !< Time step for real time measurement
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom   !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      character(len=8), intent(in) :: simid                !< Name of simulation
      character(len=1), intent(in) :: real_time_measure    !< Measurements displayed in real time


      ! Skyrmion Number
      if (skyno=='Y') then
         if (mod(sstep-1,skyno_step)==0) then

            ! Wite the total skyrmion number to the buffer
            call buffer_skyno(Natom, Mensemble,N1,N2,mstep-1,emomM,emom,&
               bcount_skyno,delta_t,real_time_measure)

            ! Write the type projected skyrmion number to the buffer
            if ( do_proj_skyno=='Y') then
               call buffer_proj_skyno(NT,Natom, Mensemble,N1,N2,mstep-1,&
                  emomM,emom,bcount_skyno,delta_t,real_time_measure,atype)
            endif
            ! Write the site dependent skyrmion number to the buffer
            if (do_skyno_den=='Y') then
               call buffer_dens_skyno(Natom, Mensemble,N1,N2,mstep-1,emomM,emom,&
                  bcount_skyno,delta_t,real_time_measure)
            endif

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
               endif

               bcount_skyno=1
            else
               bcount_skyno=bcount_skyno+1
            endif
         else

         endif
      end if

   end subroutine print_topology

   !---------------------------------------------------------------------------
   !> @brief
   !> Allocate topological data and initialization of variables
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------
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
         if (skyno=='Y') then
            allocate(skynob(skyno_buff),stat=i_stat)
            call memocc(i_stat,product(shape(skynob))*kind(skynob),'skynob','allocate_topology')
            skynob=0.0D0
            allocate(indxb_skyno(skyno_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_skyno))*kind(indxb_skyno),'indxb_skyno','allocate_topology')
            indxb_skyno=0
            allocate(grad_mom(3,3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(grad_mom))*kind(grad_mom),'grad_mom','allocate_topology')
            grad_mom=0.0D0

            if (do_proj_skyno=='Y') then
               allocate(proj_skynob(skyno_buff,NT),stat=i_stat)
               call memocc(i_stat,product(shape(proj_skynob))*kind(proj_skynob),'proj_skynob','allocate_topology')
               proj_skynob=0.0D0
               allocate(proj_grad_mom(3,3,Natom,Mensemble,NT),stat=i_stat)
               call memocc(i_stat,product(shape(proj_grad_mom))*kind(proj_grad_mom),'proj_grad_mom','allocate_topology')
               proj_grad_mom=0.0D0
            endif

            if (do_skyno_den=='Y') then
               allocate(dens_skynob(skyno_buff,Natom),stat=i_stat)
               call memocc(i_stat,product(shape(dens_skynob))*kind(dens_skynob),'dens_skynob','allocate_topology')
               dens_skynob=0.0D0
            endif

         endif

      else

         if (skyno=='Y') then
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
            endif

         endif

      end if

   end subroutine allocate_topology

   !---------------------------------------------------------------------------
   !> @brief
   !> Initialization of the variables for topological measurements
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------
   subroutine prn_topology_init()

      implicit none
      ! Topological functions
      skyno         = "N"
      do_proj_skyno = "N"
      do_skyno_den  = "N"
      skyno_step    = 100
      skyno_buff    = 10

   end subroutine prn_topology_init

   !---------------------------------------------------------------------------
   !> @brief
   !> Subroutine to flush the topological information, i.e. print to file
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------
   subroutine flush_topology(NT,Natom,simid,real_time_measure)

      implicit none

      integer, intent(in) :: NT !< Number of types of atoms
      integer, intent(in) :: Natom !< Number of atoms in the system
      character(len=8), intent(in) :: simid              !< Name of simulation
      character(len=1), intent(in) :: real_time_measure  !< Measurements displayed in real time

      ! Skyrmion number
      if (skyno=='Y') then
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
      endif

   end subroutine flush_topology

   !---------------------------------------------------------------------------
   !> @brief
   !> Printing the type dependent projected skyrmion number
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !---------------------------------------------------------------------------
   subroutine prn_skyno(simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Simulation name
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: k
      character(len=30) :: filn

      ! write the skyrmion number
      write(filn,'(''sknumber.'',a8,''.out'')') simid
      open(ofileno,file=filn, position="append")

      ! Write header to output files for first iteration
      if(indxb_skyno(1)==0) then
         write (ofileno,'(a)') "# Iteration  sk.number"
      end if

      do k=1, bcount_skyno
         if (real_time_measure=='Y') then
            write (ofileno,241) indxb_skyno(k), skynob(k)
         else
            write (ofileno,240) int(indxb_skyno(k)), skynob(k)
         endif
      enddo

      close(ofileno)

      return

      write(*,*) "Error writing the skyrmion number file"

      240 format(i8,2x,f10.4)
      241 format(es16.4,2x,f10.4)

   end subroutine prn_skyno

   !---------------------------------------------------------------------------
   !> @brief
   !> Printing the type dependent projected skyrmion number
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------
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
            write(filn,'(''proj_sknumber.'',i2,''.'',a8,''.out'')') ii,simid
         else
            write(filn,'(''proj_sknumber.'',i1,''.'',a8,''.out'')') ii,simid

         endif
         open(ofileno,file=filn, position="append")

         ! Write header to output files for first iteration
         if(indxb_skyno(1)==0) then
            write (ofileno,'(a)') "# Iteration  sk.number"
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

   end subroutine prn_proj_skyno

   !---------------------------------------------------------------------------
   !> @brief
   !> Printing the site dependent skyrmion number
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------
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
      write(filn,'(''dens_skynum.'',a8,''.out'')') simid
      open(ofileno,file=filn, position="append")

      ! Write header to output files for first iteration
      if(indxb_skyno(1)==0) then
         write (ofileno,'(a)') "# Iteration  Atomic number  Sk. number"
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

      244 format(i8,2x,i8,2x,f10.4)
      245 format(es16.4,2x,i8,2x,f10.4)

   end subroutine prn_dens_skyno

   !> Buffer the skyrmion number
   subroutine buffer_skyno(Natom, Mensemble,N1,N2,mstep,emomM,emom,&
         bcount_skyno,delta_t,real_time_measure)
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
      skynob(bcount_skyno)=pontryagin_no(Natom, Mensemble,emomM, grad_mom)&
         *(0.5d0*(N1+N2))**0.5d0

      if (real_time_measure=='Y') then
         indxb_skyno(bcount_skyno)=idnint(real(mstep,dblprec)*delta_t)
      else
         indxb_skyno(bcount_skyno)=mstep
      endif

   end subroutine buffer_skyno

   !> Buffer the skyrmion number
   subroutine buffer_proj_skyno(NT,Natom, Mensemble,N1,N2,mstep,emomM,emom,&
         bcount_skyno,delta_t,real_time_measure,atype)
      !
      use Gradients
      use Topology

      implicit none

      integer, intent(in) :: NT    !< Number of typed of atoms
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

      call proj_grad_moments(NT,Natom, Mensemble,atype,emomM,proj_grad_mom)
      proj_skynob(bcount_skyno,1:NT)=proj_pontryagin_no(NT,Natom,Mensemble,atype,emomM,proj_grad_mom)&
         *(0.5d0*(N1+N2))**0.5d0

      if (real_time_measure=='Y') then
         indxb_skyno(bcount_skyno)=idnint(real(mstep,dblprec)*delta_t)
      else
         indxb_skyno(bcount_skyno)=mstep
      endif

   end subroutine buffer_proj_skyno

   !---------------------------------------------------------------------------
   !> @brief
   !> Buffer the site dependent skyrion number
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------
   subroutine buffer_dens_skyno(Natom, Mensemble,N1,N2,mstep,emomM,emom,&
         bcount_skyno,delta_t,real_time_measure)
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

      ! .. Local variables
      integer :: iatom

      ! Call the calculation of the magnetization gradient
      ! Is this necesary? it is already called for the total skyrmion number
      call grad_moments(Natom, Mensemble,emom, grad_mom)

      ! Loop over the atomic sites
      do iatom=1,Natom
         dens_skynob(bcount_skyno,iatom)=pontryagin_no_density(iatom,Natom, Mensemble,emomM, grad_mom)&
            *(0.5d0*(N1+N2))**0.5d0
      enddo

      if (real_time_measure=='Y') then
         indxb_skyno(bcount_skyno)=idnint(real(mstep,dblprec)*delta_t)
      else
         indxb_skyno(bcount_skyno)=mstep
      endif

   end subroutine buffer_dens_skyno

end module prn_topology
