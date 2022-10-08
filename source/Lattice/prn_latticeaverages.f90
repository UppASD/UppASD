!> Module and data structures needed for printing lattice DOF averages, cumulants,
!> energies, and temperature
module prn_latticeaverages

   use Parameters
   use Profiling
   use LatticeHamiltonianData
   use inputdata, only : do_hoc_debug

   implicit none

   ! Printing definitions
   integer :: lavrg_step !< Interval for sampling average displacements
   integer :: lavrg_buff !< Buffer size for average displacements
   character(len=1) :: do_lavrg           !< Measure average displacments (Y/N)
   character(len=1) :: do_proj_lavrg      !< Measure projected displacements (Y/A/N)
   character(len=1) :: do_projch_lavrg    !< Measure chemically projected displacements (Y/N)
   character(len=1) :: do_lenergy         !< Measure susceptibility, and specific heat(Y/N)

   real(dblprec) :: totalmass   !< Total ionic mas
   real(dblprec), dimension(3) :: masscenter   !< Center of mass
   real(dblprec), dimension(:,:), allocatable :: dcoord !< Coordinates of atoms

   !!! move to subroutine start
   ! Local calculations for printing
   real(dblprec) :: uavrg      !< Average displacement
   real(dblprec) :: uxavrg     !< Average x-component of displacement
   real(dblprec) :: uyavrg     !< Average y-component of displacement
   real(dblprec) :: uzavrg     !< Average z-component of displacement

   real(dblprec) :: vavrg      !< Average velocity
   real(dblprec) :: vxavrg     !< Average x-component of velocity
   real(dblprec) :: vyavrg     !< Average y-component of velocity
   real(dblprec) :: vzavrg     !< Average z-component of velocity

   real(dblprec) :: pavrg     !< Average ionic momentum
   real(dblprec) :: pxavrg    !< Average x-component of ionic momentum
   real(dblprec) :: pyavrg    !< Average y-component of ionic momentum
   real(dblprec) :: pzavrg    !< Average z-component of ionic momentum

   real(dblprec) :: lavrg      !< Average ionic angular momentum
   real(dblprec) :: lxavrg     !< Average x-component of ionic angular momentum
   real(dblprec) :: lyavrg     !< Average y-component of ionic angular momentum
   real(dblprec) :: lzavrg     !< Average z-component of ionic angular momentum
   !!! move to subroutine end

   real(dblprec), dimension(:), allocatable       :: indxb_uavrg       !< Step counter for average displacement
   integer :: Nlavrgcum        !< Counter for number of cumulated averages

   real(dblprec), dimension(:,:,:), allocatable   :: uavrg_buff        !< Buffer for average displacements
   real(dblprec), dimension(:,:,:), allocatable   :: uavrg2_buff_proj  !< Buffer for squared projected displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: uavrg_buff_proj   !< Buffer for projected displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: uavrg_buff_projch !< Buffer for chemical projected displacements

   real(dblprec), dimension(:,:,:), allocatable   :: vavrg_buff        !< Buffer for average velocities
   real(dblprec), dimension(:,:,:), allocatable   :: vavrg2_buff_proj  !< Buffer for squared projected velocities
   real(dblprec), dimension(:,:,:,:), allocatable :: vavrg_buff_proj   !< Buffer for projected velocities
   real(dblprec), dimension(:,:,:,:), allocatable :: vavrg_buff_projch !< Buffer for chemical projected velocities

   real(dblprec), dimension(:,:,:), allocatable   :: pavrg_buff        !< Buffer for average ionic momentum
   real(dblprec), dimension(:,:,:), allocatable   :: langavrg_buff     !< Buffer for average ionic angular momentum

   real(dblprec), dimension(:,:), allocatable     :: ldpotenrg_buff    !< Buffer for average ionic potential energy
   real(dblprec), dimension(:,:), allocatable     :: sdpotenrg_buff    !< Buffer for average magnetic energy
   real(dblprec), dimension(:,:), allocatable     :: sldpotenrg_buff   !< Buffer for average spin-lattice potential energy
   real(dblprec), dimension(:,:), allocatable     :: totpotenrg_buff   !< Buffer for average total potential energy
   real(dblprec), dimension(:,:), allocatable     :: kinenrg_buff      !< Buffer for average ionic kinetic energy
   real(dblprec), dimension(:,:), allocatable     :: totenrg_buff      !< Buffer for average ionic total energy
   real(dblprec), dimension(:,:), allocatable     :: heatcap_buff      !< Buffer for ionic heat capacity
   real(dblprec), dimension(:,:), allocatable     :: iontemp_buff      !< Buffer for ionic temperature

   real(dblprec) :: ldpotenrg_cum      !< Cumulative average ionic potential energy
   real(dblprec) :: kinenrg_cum      !< Cumulative average ionic kinetic energy
   real(dblprec) :: totenrg_cum      !< Cumulative average ionic total energy
   real(dblprec) :: ldpotenrg2_cum      !< Cumulative average ionic potential energy
   real(dblprec) :: kinenrg2_cum      !< Cumulative average ionic kinetic energy
   real(dblprec) :: totenrg2_cum      !< Cumulative average ionic total energy

   integer :: bcount_uavrg    !< Counter of buffer for averages

   public

   ! Private variables
   private :: uavrg, uxavrg, uyavrg, uzavrg, vavrg, vxavrg, vyavrg, vzavrg, indxb_uavrg, &
      uavrg_buff, uavrg2_buff_proj, uavrg_buff_proj, uavrg_buff_projch, &
      vavrg_buff, vavrg2_buff_proj, vavrg_buff_proj, vavrg_buff_projch, &
      bcount_uavrg


contains


   !---------------------------------------------------------------------------
   !> @brief
   !> Wrapper for printing the ionic observables
   !>
   !> @author
   !> johan Hellsvik
   !---------------------------------------------------------------------------
   subroutine print_lattaverages(simid, Natom, Mensemble, NT, NA, N1, N2, N3, &
         mstep, sstep, logsamp, coord, mion, uvec, vvec, lvec, eeff, Temp, &
         Nchmax, do_ralloy, Natom_full, atype, achem_ch, asite_ch, &
         ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy)

      use LatticeRestart, only : prnlattrestart

      implicit none

      character(len=8), intent(in) :: simid             !< Simulation name
      integer, intent(in) :: Natom      !< Number of atoms in the system
      integer, intent(in) :: Mensemble  !< Number of ensembles
      integer, intent(in) :: NT         !< Number of types of atoms
      integer, intent(in) :: NA         !< Number of atoms in one cell
      integer, intent(in) :: N1         !< Number of cell repetitions in x direction
      integer, intent(in) :: N2         !< Number of cell repetitions in y direction
      integer, intent(in) :: N3         !< Number of cell repetitions in z direction

      integer, intent(in) :: mstep      !< Current simulation step
      integer, intent(in) :: sstep      !< Simulation step in logarithmic scale
      character(len=1) :: logsamp       !< Sample measurements logarithmically
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mion   !< Ionic mass
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec   !< Current ionic velocity
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: lvec   !< Current local angular momentum
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff !< Total effective field from application of Hamiltonian
      real(dblprec), intent(in) :: Temp                      !< Temperature of the system (scalar)

      integer, intent(in) :: Nchmax     !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, dimension(Natom), intent(in) :: atype         !< Type of atom
      integer, dimension(Natom), intent(in) :: achem_ch      !< Chemical type of atoms (reduced list) (achem_ch(i)=achtype(acellnumbrev(i)))
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system

      real(dblprec), dimension(Mensemble), intent(in) :: ldpot_energy     !< LD potential energy
      real(dblprec), dimension(Mensemble), intent(in) :: sdpot_energy        !< SD potential energy
      real(dblprec), dimension(Mensemble), intent(in) :: sldpot_energy    !< SLD potential energy (without pure LD or SD potential energies)
      real(dblprec), dimension(Mensemble), intent(in) :: totpot_energy    !< Total potential energy: LD + SD + SLD. No kinetic energy! 

      ! Local variables


      ! Averages
      if (do_lavrg=='Y') then

         if ( mod(sstep-1,lavrg_step)==0) then

            ! Write average step to buffer
            call buffer_uavrg(Natom, Mensemble, mstep-1, coord, mion, uvec, vvec, lvec, eeff, bcount_uavrg, &
               ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy)

            ! Write projected average step to buffer
            if (do_proj_lavrg=='Y' .or. do_proj_lavrg=='A') then
               call buffer_proj_uavrg(Natom, Mensemble, NA, uvec, vvec, bcount_uavrg, &
                  do_ralloy, Natom_full, asite_ch)
            endif
            ! Write chemically projected average step to buffer
            if (do_projch_lavrg=='Y') then
               call buffer_projch_uavrg(Natom, Mensemble, Nchmax, achem_ch, uvec, vvec, bcount_uavrg)
            endif

            if (bcount_uavrg==lavrg_buff) then

               ! Write the total averages buffer to file
               call prn_uavrg(simid, Natom, Mensemble, Temp)

               ! Write the projected averages buffer to file
               if (do_proj_lavrg=='Y' .or. do_proj_lavrg=='A') then
                  call prn_proj_uavrg(simid, Natom, Mensemble, NA, NT, N1, N2, N3, atype, do_proj_lavrg)
               endif

               ! Write the chemically projected averages buffer to file
               if (do_projch_lavrg=='Y') then
                  call prn_projch_uavrg(simid, Mensemble, N1, N2, N3, Nchmax)
               endif
               ! Create a restart file
               call prnlattrestart(simid, Natom, Mensemble, mstep, uvec, vvec)

               ! Reset statistics buffer
               bcount_uavrg=1

            else

               bcount_uavrg=bcount_uavrg+1

            endif

         endif

      endif

   end subroutine print_lattaverages


   !> Flush the averages, i.e. print them to file in the last iteration
   subroutine flush_lattaverages(simid, Natom, Mensemble, NT, NA, N1, N2, N3, &
         mstep, Nchmax ,atype, mion, uvec, vvec, Temp)

      use LatticeRestart, only : prnlattrestart

      implicit none

      character(len=8), intent(in) :: simid                !< Name of simulation
      integer, intent(in) :: Natom     !< Number of atoms in the system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NT        !< Number of types of atoms
      integer, intent(in) :: NA        !< Number of atoms in one cell
      integer, intent(in) :: N1        !< Number of cell repetitions in x direction
      integer, intent(in) :: N2        !< Number of cell repetitions in y direction
      integer, intent(in) :: N3        !< Number of cell repetitions in z direction
      integer, intent(in) :: mstep     !< Current simulation step
      integer, intent(in) :: Nchmax    !< Max number of chemical components on each site in cell
      integer, dimension(Natom), intent(in) :: atype       !< Type of atom
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mion   !< Ionic mass
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: vvec   !< Current ionic velocity
      real(dblprec), intent(in) :: Temp                      !< Temperature of the system (scalar)

      ! Averages
      if (do_lavrg=='Y') then

         ! Write buffer to file
         bcount_uavrg=bcount_uavrg-1
         call prn_uavrg(simid, Natom, Mensemble, Temp)

         if(do_proj_lavrg=='Y' .or. do_proj_lavrg=='A') &
            call prn_proj_uavrg(simid, Natom, Mensemble, NA, NT, N1, N2, N3, atype, do_proj_lavrg)
         if(do_projch_lavrg=='Y') call prn_projch_uavrg(simid, Mensemble, N1, N2, N3, Nchmax)

         ! Create a restart file
         call prnlattrestart(simid, Natom, Mensemble, mstep, uvec, vvec)

      endif

   end subroutine flush_lattaverages


   !> Initialization of the printing statements
   subroutine lattavrg_init()

      implicit none

      !Averages variables
      do_lavrg        = 'Y'
      do_proj_lavrg   = 'N'
      do_projch_lavrg = 'N'
      lavrg_step      = 100
      lavrg_buff      = 10

   end subroutine lattavrg_init


   !> Initialize cumulant counters
   subroutine zero_lattcumulant_counters()
      !
      implicit none

      Nlavrgcum = 0
      ldpotenrg_cum = 0.0_dblprec
      ldpotenrg2_cum = 0.0_dblprec
      kinenrg_cum = 0.0_dblprec
      kinenrg2_cum = 0.0_dblprec
      totenrg_cum = 0.0_dblprec
      totenrg2_cum = 0.0_dblprec

   end subroutine zero_lattcumulant_counters


   !> Calculate the center of mass
   subroutine calc_masscenter(Natom, Mensemble, coord, mion, uvec)

      implicit none

      integer, intent(in) :: Natom      !< Number of atoms in the system
      integer, intent(in) :: Mensemble  !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mion   !< Ionic mass
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: uvec   !< Current ionic displacement

      ! Scalar variables
      integer :: i

      ! Executable statements
      masscenter = 0.0_dblprec
      totalmass = 0.0_dblprec

      ! Sum over ions to calculate total mass and the center of mass
      do i=1, Natom
         totalmass = totalmass + mion(i,1)
         masscenter(1:3) = masscenter(1:3) + mion(i,1) * ( coord(1:3,i) + uvec(1:3,i,1) )
      end do
      write(*,*) 'Total mass ', totalmass
      masscenter(1:3) =  masscenter(1:3) / totalmass
      write(*,'(a,3f12.6)') 'Center of mass at coordinate ', masscenter(1:3)
      ! Calculate the ionic distance relative to the center of mass
      do i=1, Natom
         dcoord(1:3,i) = coord(1:3,i) - masscenter(1:3)
      end do
      !print '(3f12.6)',dcoord

   end subroutine calc_masscenter


   !> Allocation of the necessary arrays for the measurement of the displacements and velocities
   subroutine allocate_lattaverages(Natom, Mensemble, NA, NT, Nchmax, flag)

      implicit none

      integer, intent(in) :: Natom      !< Number of atoms in the system
      integer, intent(in) :: Mensemble  !< Number of ensembles
      integer, intent(in) :: NA         !< Number of atoms in one cell
      integer, intent(in) :: NT         !< Number of types of atoms
      integer, intent(in) :: Nchmax     !< Max number of chemical components on each site in cell
      integer, intent(in) :: flag       !< Allocate or deallocate (1/-1)

      !.. Local variables
      integer :: i_stat, i_all

      if(flag>0) then

         bcount_uavrg=1
         ! Allocations for the averages
         if (do_lavrg=='Y'.or.do_proj_lavrg=='Y' .or. do_proj_lavrg=='A') then
            allocate(uavrg_buff(4,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(uavrg_buff))*kind(uavrg_buff),'uavrg_buff','allocate_lattaverages')
            allocate(uavrg_buff_proj(3,NA,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(uavrg_buff_proj))*kind(uavrg_buff_proj),'uavrg_buff_proj','allocate_lattaverages')
            allocate(uavrg2_buff_proj(NA,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(uavrg2_buff_proj))*kind(uavrg2_buff_proj),'uavrg2_buff_proj','allocate_lattaverages')
            allocate(vavrg_buff(4,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(vavrg_buff))*kind(vavrg_buff),'vavrg_buff','allocate_lattaverages')
            allocate(vavrg_buff_proj(3,NA,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(vavrg_buff_proj))*kind(vavrg_buff_proj),'vavrg_buff_proj','allocate_lattaverages')
            allocate(vavrg2_buff_proj(NA,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(vavrg2_buff_proj))*kind(vavrg2_buff_proj),'vavrg2_buff_proj','allocate_lattaverages')
            allocate(pavrg_buff(4,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(pavrg_buff))*kind(pavrg_buff),'pavrg_buff','allocate_lattaverages')
            allocate(langavrg_buff(4,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(langavrg_buff))*kind(langavrg_buff),'langavrg_buff','allocate_lattaverages')
            allocate(dcoord(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(dcoord))*kind(dcoord),'dcoord','allocate_lattaverages')

            allocate(ldpotenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(ldpotenrg_buff))*kind(ldpotenrg_buff),'ldpotenrg_buff','allocate_lattaverages')
            allocate(sdpotenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(sdpotenrg_buff))*kind(sdpotenrg_buff),'sdpotenrg_buff','allocate_lattaverages')
            allocate(sldpotenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(sldpotenrg_buff))*kind(sldpotenrg_buff),'sldpotenrg_buff','allocate_lattaverages')
            allocate(totpotenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(totpotenrg_buff))*kind(totpotenrg_buff),'totpotenrg_buff','allocate_lattaverages')
            allocate(kinenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(kinenrg_buff))*kind(kinenrg_buff),'kinenrg_buff','allocate_lattaverages')
            allocate(totenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(totenrg_buff))*kind(totenrg_buff),'totenrg_buff','allocate_lattaverages')

            allocate(heatcap_buff(lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(heatcap_buff))*kind(heatcap_buff),'heatcap_buff','allocate_lattaverages')
            allocate(iontemp_buff(lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(iontemp_buff))*kind(iontemp_buff),'iontemp_buff','allocate_lattaverages')

         endif

         ! Index array should be allocated fpr all kind of possible measurements
         if (do_lavrg=='Y'.or.do_proj_lavrg=='Y' .or. do_proj_lavrg=='A'.or.do_projch_lavrg=='Y') then
            allocate(indxb_uavrg(lavrg_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_uavrg))*kind(indxb_uavrg),'indxb_uavrg','allocate_lattaverages')
         endif

         ! Allocations for chemically projected averages
         if (do_projch_lavrg=='Y') then
            allocate(uavrg_buff_projch(4,Nchmax,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(uavrg_buff_projch))*kind(uavrg_buff_projch),'uavrg_buff_projch','allocate_lattaverages')
            allocate(vavrg_buff_projch(4,Nchmax,lavrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(vavrg_buff_projch))*kind(vavrg_buff_projch),'vavrg_buff_projch','allocate_lattaverages')
         end if

      else

         ! Deallocations for averages
         if (do_lavrg=='Y'.or.do_proj_lavrg=='Y' .or. do_proj_lavrg=='A') then
            i_all=-product(shape(uavrg_buff))*kind(uavrg_buff)
            deallocate(uavrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'uavrg_buff','allocate_lattaverages')
            i_all=-product(shape(uavrg_buff_proj))*kind(uavrg_buff_proj)
            deallocate(uavrg_buff_proj,stat=i_stat)
            call memocc(i_stat,i_all,'uavrg_buff_proj','allocate_lattaverages')
            i_all=-product(shape(uavrg2_buff_proj))*kind(uavrg2_buff_proj)
            deallocate(uavrg2_buff_proj,stat=i_stat)
            call memocc(i_stat,i_all,'uavrg2_buff_proj','allocate_lattaverages')
            i_all=-product(shape(vavrg_buff))*kind(vavrg_buff)
            deallocate(vavrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'vavrg_buff','allocate_lattaverages')
            i_all=-product(shape(vavrg_buff_proj))*kind(vavrg_buff_proj)
            deallocate(vavrg_buff_proj,stat=i_stat)
            call memocc(i_stat,i_all,'vavrg_buff_proj','allocate_lattaverages')
            i_all=-product(shape(vavrg2_buff_proj))*kind(vavrg2_buff_proj)
            deallocate(vavrg2_buff_proj,stat=i_stat)
            call memocc(i_stat,i_all,'vavrg2_buff_proj','allocate_lattaverages')
            i_all=-product(shape(pavrg_buff))*kind(pavrg_buff)
            deallocate(pavrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'pavrg_buff','allocate_lattaverages')
            i_all=-product(shape(lavrg_buff))*kind(lavrg_buff)
            deallocate(langavrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'langavrg_buff','allocate_lattaverages')
            i_all=-product(shape(dcoord))*kind(dcoord)
            deallocate(dcoord,stat=i_stat)
            call memocc(i_stat,i_all,'dcoord','allocate_lattaverages')

            i_all=-product(shape(ldpotenrg_buff))*kind(ldpotenrg_buff)
            deallocate(ldpotenrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'ldpotenrg_buff','allocate_lattaverages')
            i_all=-product(shape(sdpotenrg_buff))*kind(sdpotenrg_buff)
            deallocate(sdpotenrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'sdpotenrg_buff','allocate_lattaverages')
            i_all=-product(shape(sldpotenrg_buff))*kind(sldpotenrg_buff)
            deallocate(sldpotenrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'sldpotenrg_buff','allocate_lattaverages')
            i_all=-product(shape(totpotenrg_buff))*kind(totpotenrg_buff)
            deallocate(totpotenrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'totpotenrg_buff','allocate_lattaverages')
            i_all=-product(shape(kinenrg_buff))*kind(kinenrg_buff)
            deallocate(kinenrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'kinenrg_buff','allocate_lattaverages')
            i_all=-product(shape(totenrg_buff))*kind(totenrg_buff)
            deallocate(totenrg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'totenrg_buff','allocate_lattaverages')

            i_all=-product(shape(heatcap_buff))*kind(heatcap_buff)
            deallocate(heatcap_buff,stat=i_stat)
            call memocc(i_stat,i_all,'heatcap_buff','allocate_lattaverages')
            i_all=-product(shape(iontemp_buff))*kind(iontemp_buff)
            deallocate(iontemp_buff,stat=i_stat)
            call memocc(i_stat,i_all,'iontemp_buff','allocate_lattaverages')

         endif
         ! Index array should be allocated for all kind of possible measurements
         if (do_lavrg=='Y'.or.do_proj_lavrg=='Y' .or. do_proj_lavrg=='A'.or.do_projch_lavrg=='Y') then
            i_all=-product(shape(indxb_uavrg))*kind(indxb_uavrg)
            deallocate(indxb_uavrg,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_uavrg','allocate_lattaverages')
         endif

         ! Deallocations for chemically projected averages
         if(do_projch_lavrg=='Y') then
            i_all=-product(shape(uavrg_buff_projch))*kind(uavrg_buff_projch)
            deallocate(uavrg_buff_projch,stat=i_stat)
            call memocc(i_stat,i_all,'uavrg_buff_projch','allocate_lattaverages')
            i_all=-product(shape(vavrg_buff_projch))*kind(vavrg_buff_projch)
            deallocate(vavrg_buff_projch,stat=i_stat)
            call memocc(i_stat,i_all,'vavrg_buff_projch','allocate_lattaverages')
         end if

      endif

   end subroutine allocate_lattaverages


   !> Buffer average magnetization
   subroutine buffer_uavrg(Natom, Mensemble, mstep, coord, mion, uvec, vvec, lvec, eeff, bcount_avrg, &
         ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy)       
      !subroutine buffer_uavrg(Natom, Mensemble, mstep, coord, mion, uvec, vvec, eeff, bcount_avrg, sld_energy)

      use Constants

      !.. Implicit declarations
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble   !< Number of ensembles
      integer, intent(in) :: mstep !< Current simulation step
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mion   !< Ionic mass
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: vvec   !< Current ionic velocity
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: lvec   !< Current local angular momentum
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff !< Total effective field from application of Hamiltonian
      integer, intent(in) :: bcount_avrg !< Counter of buffer for averages
      real(dblprec), dimension(Mensemble), intent(in) :: ldpot_energy     !< LD potential energy
      real(dblprec), dimension(Mensemble), intent(in) :: sdpot_energy        !< SD potential energy
      real(dblprec), dimension(Mensemble), intent(in) :: sldpot_energy    !< SLD potential energy (without pure LD or SD potential energies)
      real(dblprec), dimension(Mensemble), intent(in) :: totpot_energy    !< Total potential energy: LD + SD + SLD. No kinetic energy! 

      !.. Scalar variables
      integer :: i, k

      !.. Local arrays
      real(dblprec), dimension(Mensemble) ::  ux, uy, uz
      real(dblprec), dimension(Mensemble) ::  vx, vy, vz
      real(dblprec), dimension(Mensemble) ::  px, py, pz
      real(dblprec), dimension(Mensemble) ::  lx, ly, lz
      real(dblprec), dimension(Mensemble) ::  ldpotenrg, sdpotenrg, sldpotenrg
      real(dblprec), dimension(Mensemble) ::  kinenrg, totpotenrg, totenrg
      real(dblprec) ::  rpx, rpy, rpz
      real(dblprec) ::  dx, dy, dz
      real(dblprec) ::  v2

      !.. Executable statements
      ux(:)=0.0_dblprec
      uy(:)=0.0_dblprec
      uz(:)=0.0_dblprec
      vx(:)=0.0_dblprec
      vy(:)=0.0_dblprec
      vz(:)=0.0_dblprec
      px(:)=0.0_dblprec
      py(:)=0.0_dblprec
      pz(:)=0.0_dblprec
      lx(:)=0.0_dblprec
      ly(:)=0.0_dblprec
      lz(:)=0.0_dblprec
      rpx=0.0_dblprec
      rpy=0.0_dblprec
      rpz=0.0_dblprec

      ldpotenrg=0.0_dblprec
      sdpotenrg=0.0_dblprec
      sldpotenrg=0.0_dblprec
      totpotenrg=0.0_dblprec
      kinenrg=0.0_dblprec
      totenrg=0.0_dblprec

      !.. Sum over atoms
      !write(*,*) '------------------------'
      ! Transpose loops over i and k for efficiency?
      do i=1, Natom
         do k=1,Mensemble
            ux(k) = ux(k) + uvec(1,i,k)
            uy(k) = uy(k) + uvec(2,i,k)
            uz(k) = uz(k) + uvec(3,i,k)
            vx(k) = vx(k) + vvec(1,i,k)
            vy(k) = vy(k) + vvec(2,i,k)
            vz(k) = vz(k) + vvec(3,i,k)
            px(k) = px(k) + mion(i,1) * vvec(1,i,k)
            py(k) = py(k) + mion(i,1) * vvec(2,i,k)
            pz(k) = pz(k) + mion(i,1) * vvec(3,i,k)
            !write(*,*) '------------------------'
            !!! !write(*,10004) 'dcoord ', i, dcoord(1:3,i)
            !!! dx = 0*dcoord(1,i) + uvec(1,i,k)
            !!! dy = 0*dcoord(2,i) + uvec(2,i,k)
            !!! dz = 0*dcoord(3,i) + uvec(3,i,k)
            !!! rpx = dy * vvec(3,i,k) - dz * vvec(2,i,k)
            !!! rpy = dz * vvec(1,i,k) - dx * vvec(3,i,k)
            !!! rpz = dx * vvec(2,i,k) - dy * vvec(1,i,k)
            !!! !write(*,10004) 'rpx, rpy, rpz ', i, rpx, rpy, rpz
            !!! lx(k) = lx(k) + mion(i,1) * rpx
            !!! ly(k) = ly(k) + mion(i,1) * rpy
            !!! lz(k) = lz(k) + mion(i,1) * rpz
            !write(*,10004) 'lx, ly, lz    ', i, lx(k), ly(k), lz(k)
            lx(k) = lx(k) + lvec(1,i,k)
            ly(k) = ly(k) + lvec(2,i,k)
            lz(k) = lz(k) + lvec(3,i,k)
            

            ldpotenrg(k)  = ldpot_energy(k)
            sdpotenrg(k)  = sdpot_energy(k)
            sldpotenrg(k) = sldpot_energy(k)
            totpotenrg(k) = totpot_energy(k)

            v2 = vvec(1,i,k)*vvec(1,i,k) + vvec(2,i,k)*vvec(2,i,k) + vvec(3,i,k)*vvec(3,i,k)
            ! Convert the kinetic energy to mRyd
            kinenrg(k) = kinenrg(k) + 0.5_dblprec * amu * mion(i,1) * angstrom**2 * v2 / mry
            !kinenrg(k) = kinenrg(k) + 0.5_dblprec * amu * mion(i,1) * angstrom**2 * v2 / ev

            totenrg(k) = totpotenrg(k) + kinenrg(k)

         end do
      end do

      !.. Save in buffer
      do k=1,Mensemble
         uavrg_buff(1,bcount_avrg,k) = ux(k)/Natom
         uavrg_buff(2,bcount_avrg,k) = uy(k)/Natom
         uavrg_buff(3,bcount_avrg,k) = uz(k)/Natom
         vavrg_buff(1,bcount_avrg,k) = vx(k)/Natom
         vavrg_buff(2,bcount_avrg,k) = vy(k)/Natom
         vavrg_buff(3,bcount_avrg,k) = vz(k)/Natom
         pavrg_buff(1,bcount_avrg,k) = px(k)/Natom
         pavrg_buff(2,bcount_avrg,k) = py(k)/Natom
         pavrg_buff(3,bcount_avrg,k) = pz(k)/Natom
         langavrg_buff(1,bcount_avrg,k) = lx(k)/Natom
         langavrg_buff(2,bcount_avrg,k) = ly(k)/Natom
         langavrg_buff(3,bcount_avrg,k) = lz(k)/Natom
         ldpotenrg_buff(bcount_avrg,k) = ldpotenrg(k)/Natom
         sdpotenrg_buff(bcount_avrg,k) = sdpotenrg(k)/Natom
         sldpotenrg_buff(bcount_avrg,k) = sldpotenrg(k)/Natom
         totpotenrg_buff(bcount_avrg,k) = totpotenrg(k)/Natom
         kinenrg_buff(bcount_avrg,k) = kinenrg(k)/Natom
         ! Decide if total energy or only lattice contribution should be measured here
         !totenrg_buff(bcount_avrg,k) = totenrg(k)/Natom - sdpotenrg(k)/Natom
         totenrg_buff(bcount_avrg,k) = kinenrg(k)/Natom + ldpotenrg(k)/Natom + sldpotenrg(k)/Natom
         !totenrg_buff(bcount_avrg,k) = totenrg(k)/Natom
      end do

      indxb_uavrg(bcount_avrg) = mstep

      10004 format (a, i8,30es16.8)

   end subroutine buffer_uavrg



   !> Buffer site projected average displacement
   subroutine buffer_proj_uavrg(Natom, Mensemble, NA, uvec, vvec, bcount_avrg, do_ralloy, Natom_full, asite_ch)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: vvec   !< Current ionic velocity
      integer, intent(in) :: bcount_avrg !< Counter of buffer for averages
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system

      !.. Locals
      real(dblprec), dimension(3,NA,Mensemble) ::  avrg_uvec !< Projected average displacement
      real(dblprec), dimension(3,NA,Mensemble) ::  avrg_vvec !< Projected average velocity

      !.. Scalar variables
      integer :: i,i_na,k

      !.. Executable statements
      avrg_uvec=0.0_dblprec
      avrg_vvec=0.0_dblprec

      !.. Sum over moments
      do k=1,Mensemble
         if(do_ralloy==0) then
            do i=1, Natom
               i_na=mod(i-1,NA)+1
               avrg_uvec(1,i_na,k) = avrg_uvec(1,i_na,k) + uvec(1,i,k)
               avrg_uvec(2,i_na,k) = avrg_uvec(2,i_na,k) + uvec(2,i,k)
               avrg_uvec(3,i_na,k) = avrg_uvec(3,i_na,k) + uvec(3,i,k)
               avrg_vvec(1,i_na,k) = avrg_vvec(1,i_na,k) + vvec(1,i,k)
               avrg_vvec(2,i_na,k) = avrg_vvec(2,i_na,k) + vvec(2,i,k)
               avrg_vvec(3,i_na,k) = avrg_vvec(3,i_na,k) + vvec(3,i,k)
            end do
         else
            do i=1, Natom
               i_na=asite_ch(i)
               avrg_uvec(1,i_na,k) = avrg_uvec(1,i_na,k) + uvec(1,i,k)
               avrg_uvec(2,i_na,k) = avrg_uvec(2,i_na,k) + uvec(2,i,k)
               avrg_uvec(3,i_na,k) = avrg_uvec(3,i_na,k) + uvec(3,i,k)
               avrg_vvec(1,i_na,k) = avrg_vvec(1,i_na,k) + vvec(1,i,k)
               avrg_vvec(2,i_na,k) = avrg_vvec(2,i_na,k) + vvec(2,i,k)
               avrg_vvec(3,i_na,k) = avrg_vvec(3,i_na,k) + vvec(3,i,k)
            end do
         end if

         !.. Save in buffer
         do i_na=1,NA
            uavrg_buff_proj(1,i_na,bcount_avrg,k) = avrg_uvec(1,i_na,k)
            uavrg_buff_proj(2,i_na,bcount_avrg,k) = avrg_uvec(2,i_na,k)
            uavrg_buff_proj(3,i_na,bcount_avrg,k) = avrg_uvec(3,i_na,k)
            uavrg2_buff_proj(i_na,bcount_avrg,k)  = (avrg_uvec(1,i_na,k)**2+avrg_uvec(3,i_na,k)**2+avrg_uvec(3,i_na,k)**2)
            vavrg_buff_proj(1,i_na,bcount_avrg,k) = avrg_vvec(1,i_na,k)
            vavrg_buff_proj(2,i_na,bcount_avrg,k) = avrg_vvec(2,i_na,k)
            vavrg_buff_proj(3,i_na,bcount_avrg,k) = avrg_vvec(3,i_na,k)
            vavrg2_buff_proj(i_na,bcount_avrg,k)  = (avrg_vvec(1,i_na,k)**2+avrg_vvec(3,i_na,k)**2+avrg_vvec(3,i_na,k)**2)
         end do

      end do
   end subroutine buffer_proj_uavrg


   !> Buffer site and chemical projected average magnetizations
   subroutine buffer_projch_uavrg(Natom, Mensemble, Nchmax, achem_ch, uvec, vvec, bcount_avrg)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, dimension(Natom), intent(in) :: achem_ch !< Chemical type of atoms (reduced list) (achem_ch(i)=achtype(acellnumbrev(i)))
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: vvec   !< Current ionic velocity
      integer, intent(in) :: bcount_avrg !< Counter of buffer for averages

      !.. Locals
      real(dblprec), dimension(4,Nchmax,Mensemble) ::  avrg_uvec !< Chemically projected average displacement
      real(dblprec), dimension(4,Nchmax,Mensemble) ::  avrg_vvec !< Chemically projected average velocity
      integer,dimension(nchmax) :: counter
      real(dblprec),dimension(nchmax) :: conc

      !.. Scalar variables
      integer :: i,ich,k

      !.. Executable statements
      avrg_uvec=0.0_dblprec
      avrg_vvec=0.0_dblprec

      !.. Sum over moments
      do k=1,Mensemble
         counter=0 ; conc=0.0_dblprec
         do i=1, Natom
            ich=achem_ch(i)
            counter(ich)=counter(ich)+1
            avrg_uvec(1,ich,k) = avrg_uvec(1,ich,k) + uvec(1,i,k)
            avrg_uvec(2,ich,k) = avrg_uvec(2,ich,k) + uvec(2,i,k)
            avrg_uvec(3,ich,k) = avrg_uvec(3,ich,k) + uvec(3,i,k)
            avrg_uvec(4,ich,k) = avrg_uvec(4,ich,k) + sqrt(uvec(1,i,k)**2+uvec(2,i,k)**2+uvec(3,i,k)**2)
            avrg_vvec(1,ich,k) = avrg_vvec(1,ich,k) + vvec(1,i,k)
            avrg_vvec(2,ich,k) = avrg_vvec(2,ich,k) + vvec(2,i,k)
            avrg_vvec(3,ich,k) = avrg_vvec(3,ich,k) + vvec(3,i,k)
            avrg_vvec(4,ich,k) = avrg_vvec(4,ich,k) + sqrt(vvec(1,i,k)**2+vvec(2,i,k)**2+vvec(3,i,k)**2)
         end do
         conc=(1._dblprec*counter/natom)
         !.. Save in buffer
         do ich=1,Nchmax
            uavrg_buff_projch(1,ich,bcount_avrg,k) = avrg_uvec(1,ich,k)/conc(ich)
            uavrg_buff_projch(2,ich,bcount_avrg,k) = avrg_uvec(2,ich,k)/conc(ich)
            uavrg_buff_projch(3,ich,bcount_avrg,k) = avrg_uvec(3,ich,k)/conc(ich)
            uavrg_buff_projch(4,ich,bcount_avrg,k) = avrg_uvec(4,ich,k)/conc(ich)
            vavrg_buff_projch(1,ich,bcount_avrg,k) = avrg_uvec(1,ich,k)/conc(ich)
            vavrg_buff_projch(2,ich,bcount_avrg,k) = avrg_uvec(2,ich,k)/conc(ich)
            vavrg_buff_projch(3,ich,bcount_avrg,k) = avrg_uvec(3,ich,k)/conc(ich)
            vavrg_buff_projch(4,ich,bcount_avrg,k) = avrg_uvec(4,ich,k)/conc(ich)
         end do
      end do

   end subroutine buffer_projch_uavrg


   !> Print instantaneous and cumulated average displacements, velocities, and energies
   subroutine prn_uavrg(simid, Natom, Mensemble, Temp)

      use Constants

      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: Temp    !< Temperature of the system (scalar)

      !.. Scalar variables
      integer :: i,k
      character(len=30) :: filn, filn2, filn3, filn4
      real(dblprec), dimension(Mensemble) :: avrgux, avrguy, avrguz, avrgu
      real(dblprec), dimension(Mensemble) :: avrgvx, avrgvy, avrgvz, avrgv
      real(dblprec), dimension(Mensemble) :: avrgpx, avrgpy, avrgpz, avrgp
      real(dblprec), dimension(Mensemble) :: avrglx, avrgly, avrglz, avrgl
      real(dblprec), dimension(Mensemble) :: ldpotenrg, sdpotenrg, sldpotenrg
      real(dblprec), dimension(Mensemble) :: totpotenrg, kinenrg, totenrg
      !real(dblprec), dimension(Mensemble) :: potenrg, kinenrg, totenrg

      real(dblprec) :: avrguxm, avrguym, avrguzm, avrgum, avrgus
      real(dblprec) :: avrgvxm, avrgvym, avrgvzm, avrgvm, avrgvs
      real(dblprec) :: avrgpxm, avrgpym, avrgpzm, avrgpm, avrgps
      real(dblprec) :: avrglxm, avrglym, avrglzm, avrglm, avrgls

      real(dblprec) :: ldpotenrgm, sdpotenrgm, sldpotenrgm
      real(dblprec) :: totpotenrgm, kinenrgm, totenrgm
      real(dblprec) :: ldpotenrg2, kinenrg2, totenrg2

      real(dblprec) :: totenrgs, heatcap, iontemp

      real(dblprec) :: spinfac

      write (filn,'(''lattaverages.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")
      write (filn2,'(''lattmomenta.'',a,''.out'')') trim(simid)
      open(ofileno2, file=filn2, position="append")
      write (filn3,'(''lattenergy.'',a,''.out'')') trim(simid)
      open(ofileno3, file=filn3, position="append")

      ! Write header to output files for first iteration
      if(abs(indxb_uavrg(1))<=0.0e0_dblprec) then
         ! Averages
         write (ofileno,'(a)') "# Iter.       u_avg_x         u_avg_y         u_avg_z         u_avg           u_stdev         &
            & v_avg_x         v_avg_y         v_avg_z         v_avg           v_stdev"
         write (ofileno2,'(a)') "# Iter.       p_avg_x         p_avg_y         p_avg_z         p_avg           p_stdev         &
            & l_avg_x         l_avg_y         l_avg_z         l_avg           l_stdev"
         write (ofileno3,'(a)') "# Iter.     ld_potenrg      sd_potenrg      sld_potenrg     totpot_enrg     kin_enrg       &
            & tot_enrg       ionic temp   "
         !& tot_enrg       tot_enrg_stdv"
         !write (ofileno3,'(a)')"# Iter.       potenrg         kinenrg         totenrg         totenrgs"
         !write (ofileno3,'(a)')"# Iter.       potenrg         potenrg2        kinenrg         kinenrg2        totenrg         &
         !   & totenerg2       totenrgs"
      end if

      !
      spinfac=0.5_dblprec*amu*angstrom/hbar*angstrom


      !
      do i=1, bcount_uavrg
         do k=1, Mensemble
            avrgux(k) = uavrg_buff(1,i,k)
            avrguy(k) = uavrg_buff(2,i,k)
            avrguz(k) = uavrg_buff(3,i,k)
            avrgvx(k) = vavrg_buff(1,i,k)
            avrgvy(k) = vavrg_buff(2,i,k)
            avrgvz(k) = vavrg_buff(3,i,k)
            avrgpx(k) = pavrg_buff(1,i,k)
            avrgpy(k) = pavrg_buff(2,i,k)
            avrgpz(k) = pavrg_buff(3,i,k)
            avrglx(k) = langavrg_buff(1,i,k)
            avrgly(k) = langavrg_buff(2,i,k)
            avrglz(k) = langavrg_buff(3,i,k)
            ldpotenrg(k) = ldpotenrg_buff(i,k)
            sdpotenrg(k) = sdpotenrg_buff(i,k)
            sldpotenrg(k) = sldpotenrg_buff(i,k)
            totpotenrg(k) = totpotenrg_buff(i,k)
            kinenrg(k) = kinenrg_buff(i,k)
            totenrg(k) = totenrg_buff(i,k)
         end do

         do k=1,Mensemble
            avrgu(k) = avrgux(k)**2+avrguy(k)**2+avrguz(k)**2
            avrgu(k) = sqrt(avrgu(k))
            avrgv(k) = avrgvx(k)**2+avrgvy(k)**2+avrgvz(k)**2
            avrgv(k) = sqrt(avrgv(k))
            avrgp(k) = avrgpx(k)**2+avrgpy(k)**2+avrgpz(k)**2
            avrgp(k) = sqrt(avrgp(k))
            avrgl(k) = avrglx(k)**2+avrgly(k)**2+avrglz(k)**2
            avrgl(k) = sqrt(avrgl(k))
         end do

         ! ensemble mean displacement and velocity
         avrguxm=0_dblprec
         avrguym=0_dblprec
         avrguzm=0_dblprec
         avrgum=0_dblprec
         avrgvxm=0_dblprec
         avrgvym=0_dblprec
         avrgvzm=0_dblprec
         avrgvm=0_dblprec
         avrgpxm=0_dblprec
         avrgpym=0_dblprec
         avrgpzm=0_dblprec
         avrgpm=0_dblprec
         avrglxm=0_dblprec
         avrglym=0_dblprec
         avrglzm=0_dblprec
         avrglm=0_dblprec
         ! ensemble stdev displacement and velocity
         avrgus=0_dblprec
         avrgvs=0_dblprec
         avrgps=0_dblprec
         avrgls=0_dblprec
         ! ensemble mean and stdev energies
         ldpotenrgm=0_dblprec
         sdpotenrgm=0_dblprec
         sldpotenrgm=0_dblprec
         totpotenrgm=0_dblprec
         kinenrgm=0_dblprec
         totenrgm=0_dblprec
         totenrgs=0_dblprec

         do k=1,Mensemble
            avrguxm = avrguxm + avrgux(k)
            avrguym = avrguym + avrguy(k)
            avrguzm = avrguzm + avrguz(k)
            avrgum = avrgum + avrgu(k)
            avrgus = avrgus + avrgu(k)**2
            avrgvxm = avrgvxm + avrgvx(k)
            avrgvym = avrgvym + avrgvy(k)
            avrgvzm = avrgvzm + avrgvz(k)
            avrgvm = avrgvm + avrgv(k)
            avrgvs = avrgvs + avrgv(k)**2
            avrgpxm = avrgpxm + avrgpx(k)
            avrgpym = avrgpym + avrgpy(k)
            avrgpzm = avrgpzm + avrgpz(k)
            avrgpm = avrgpm + avrgp(k)
            avrgps = avrgps + avrgp(k)**2
            avrglxm = avrglxm + avrglx(k)
            avrglym = avrglym + avrgly(k)
            avrglzm = avrglzm + avrglz(k)
            avrglm = avrglm + avrgl(k)
            avrgls = avrgls + avrgl(k)**2
            ldpotenrgm = ldpotenrgm + ldpotenrg(k) 
            sdpotenrgm = sdpotenrgm + sdpotenrg(k) 
            sldpotenrgm = sldpotenrgm + sldpotenrg(k) 
            totpotenrgm = totpotenrgm + totpotenrg(k) 
            kinenrgm = kinenrgm + kinenrg(k) 
            totenrgm = totenrgm + totenrg(k) 
            totenrgs = totenrgs + totenrg(k)**2 
         end do

         avrguxm=avrguxm/Mensemble
         avrguym=avrguym/Mensemble
         avrguzm=avrguzm/Mensemble
         avrgum=avrgum/Mensemble
         avrgvxm=avrgvxm/Mensemble
         avrgvym=avrgvym/Mensemble
         avrgvzm=avrgvzm/Mensemble
         avrgvm=avrgvm/Mensemble
         avrgpxm=avrgpxm/Mensemble
         avrgpym=avrgpym/Mensemble
         avrgpzm=avrgpzm/Mensemble
         avrgpm=avrgpm/Mensemble
         avrglxm=avrglxm/Mensemble
         avrglym=avrglym/Mensemble
         avrglzm=avrglzm/Mensemble
         avrglm=avrglm/Mensemble
         ldpotenrgm = ldpotenrgm/Mensemble
         sdpotenrgm = sdpotenrgm/Mensemble
         sldpotenrgm = sldpotenrgm/Mensemble
         totpotenrgm = totpotenrgm/Mensemble
         kinenrgm = kinenrgm/Mensemble
         totenrgm = totenrgm/Mensemble

         !Do not use Mensemble-1 for numerical convenience
         avrgus=avrgus/Mensemble - avrgum**2
         avrgvs=avrgvs/Mensemble - avrgvm**2
         avrgps=avrgps/Mensemble - avrgpm**2
         avrgls=avrgls/Mensemble - avrglm**2
         totenrgs=totenrgs/Mensemble - totenrgm**2

         !Filter out negative errors
         if(avrgus<0) then
            avrgus=0
         else
            avrgus=sqrt(avrgus)
         end if
         if(avrgvs<0) then
            avrgvs=0
         else
            avrgvs=sqrt(avrgvs)
         end if
         if(avrgps<0) then
            avrgps=0
         else
            avrgps=sqrt(avrgps)
         end if
         if(avrgls<0) then
            avrgls=0
         else
            avrgls=sqrt(avrgls)
         end if

         write (ofileno,10004) int(indxb_uavrg(i)), avrguxm, avrguym, avrguzm, avrgum, avrgus, &
            &  avrgvxm, avrgvym, avrgvzm, avrgvm, avrgvs
         write (ofileno2,10004) int(indxb_uavrg(i)), avrgpxm, avrgpym, avrgpzm, avrgpm, avrgps, &
            &  avrglxm*spinfac, avrglym*spinfac, avrglzm*spinfac, avrglm*spinfac, avrgls*spinfac
         write (ofileno3,10004) int(indxb_uavrg(i)), ldpotenrgm, sdpotenrgm, sldpotenrgm, totpotenrgm, kinenrgm, totenrgm, &
            kinenrgm * mry * 2 / (3 * k_bolt)
         !write (ofileno3,10004) int(indxb_uavrg(i)), ldpotenrgm, sdpotenrgm, sldpotenrgm, totpotenrgm, kinenrgm, totenrgm, &
         !     kinenrgm * ev * 2 / (3 * k_bolt)

      end do

      close(ofileno)
      close(ofileno2)
      close(ofileno3)

      write (filn4,'(''lattcumenergy.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn4, position="append")
      ! Write header to output files for first iteration
      if(abs(indxb_uavrg(1))<=0.0e0_dblprec) then
         write (ofileno,'(a)')"# Iter.       Cum potenrg     Cum kinenrg     Cum totenrg     Stdev totenrg   Heat cap        &
            & Temperature"
      end if

      do i=1, bcount_uavrg
         do k=1,Mensemble

            ldpotenrg2 = ldpotenrg_buff(i,k)**2
            kinenrg2 = kinenrg_buff(i,k)**2
            totenrg2 = totenrg_buff(i,k)**2

            ldpotenrg_cum = ( Nlavrgcum * ldpotenrg_cum + ldpotenrg_buff(i,k) ) / ( Nlavrgcum + 1 )
            ldpotenrg2_cum = ( Nlavrgcum * ldpotenrg2_cum + ldpotenrg2 ) / ( Nlavrgcum + 1 )
            kinenrg_cum = ( Nlavrgcum * kinenrg_cum + kinenrg_buff(i,k) ) / ( Nlavrgcum + 1 )
            kinenrg2_cum = ( Nlavrgcum * kinenrg2_cum + kinenrg2 ) / ( Nlavrgcum + 1 )
            !totenrg_cum = ( Nlavrgcum * totenrg_cum + totenrg_buff(i,k) ) / ( Nlavrgcum + 1 )
            !totenrg2_cum = ( Nlavrgcum * totenrg2_cum + totenrg_buff(i,k)**2 ) / ( Nlavrgcum + 1 )
            totenrg_cum = ( totenrg_cum + totenrg_buff(i,k) ) 
            totenrg2_cum = ( totenrg2_cum + totenrg_buff(i,k)**2 ) 
            !totenrg2_cum = ( Nlavrgcum * totenrg2_cum + totenrg2 ) / ( Nlavrgcum + 1 )

            !totenrgs = totenrg2_cum - totenrg_cum**2
            totenrgs = totenrg_cum**2 - totenrg2_cum
            totenrgs=totenrgs/(Nlavrgcum+1)**2

!           print '(a,i8,3g20.10)', 'Cumu: ',Nlavrgcum, totenrg_cum**2, totenrg2_cum, mry**2 / (k_bolt**2) / Temp**2
            ! cv = ( avrgen2 -avrgen**2 ) * mry**2 * Natom  / (k_bolt) / T**2 !* Natom   ! SI
            !cv = ( avrgen2 -avrgen**2 ) * mry**2 * Natom  / (k_bolt**2) / Temp**2 !* Natom   ! units of k_B
            if(Temp>0.0_dblprec) then
               heatcap = totenrgs *mry**2 / (k_bolt**2) / Temp**2 !* Natom   ! units of k_B
            else
               heatcap = 0.0_dblprec
            end if
            iontemp = kinenrg_cum * 2.0_dblprec * mry / (3.0_dblprec * k_bolt)
            write (ofileno,10004) int(indxb_uavrg(i)), ldpotenrg_cum, kinenrg_cum, totenrg_cum/ ( Nlavrgcum + 1 ), totenrgs, &
               heatcap, iontemp
            Nlavrgcum = Nlavrgcum + 1

         end do
      end do

      close(ofileno)
      return

      write (*,*) "Error writing the averages file"

      10004 format (i8,30es16.8)
      10005 format (es12.4,30es16.8)

   end subroutine prn_uavrg


   !> Print site projected average magnetizations
   subroutine prn_proj_uavrg(simid, Natom, Mensemble, NA, NT, N1, N2, N3, atype, do_proj_avrg)

      use Constants
      !.. Implicit declarations
      implicit none
      character(len=8), intent(in) :: simid          !< Name of simulation
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      character(len=1), intent(in) :: do_proj_avrg   !< Measure projected averages (Y/A/N)

      !.. Scalar variables
      integer :: i,i_na,j,k,nproj
      character(len=30) :: filn

      real(dblprec), dimension(Mensemble,NA) :: avrgux, avrguy, avrguz, avrgu, avrgu2, avrgu2t
      real(dblprec), dimension(NA) :: avrgum, avrgus, avrguxm, avrguym, avrguzm, avrgum2
      real(dblprec), dimension(Mensemble,NA) :: avrgvx, avrgvy, avrgvz, avrgv, avrgv2, avrgv2t
      real(dblprec), dimension(NA) :: avrgvm, avrgvs, avrgvxm, avrgvym, avrgvzm, avrgvm2

      i_na=-1;

      !.. Executable statements
      write (filn,'(''lattprojavgs.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")

      ! SLDTODO This if-statement seems unnecessary
      do i=1, bcount_uavrg
         if(do_proj_avrg=='Y') then
            do j=1, Mensemble
               do k=1, NT
                  avrgux(j,k)=0_dblprec
                  avrguy(j,k)=0_dblprec
                  avrguz(j,k)=0_dblprec
                  avrgu2t(j,k)=0_dblprec
                  avrgvx(j,k)=0_dblprec
                  avrgvy(j,k)=0_dblprec
                  avrgvz(j,k)=0_dblprec
                  avrgv2t(j,k)=0_dblprec

                  do i_na=1, NA
                     if(k==atype(i_na)) then
                        avrgux(j,k) = avrgux(j,k) + uavrg_buff_proj(1,i_na,i,j)/(N1*N2*N3)
                        avrguy(j,k) = avrguy(j,k) + uavrg_buff_proj(2,i_na,i,j)/(N1*N2*N3)
                        avrguz(j,k) = avrguz(j,k) + uavrg_buff_proj(3,i_na,i,j)/(N1*N2*N3)
                        avrgu2t(j,k) = avrgu2t(j,k) + uavrg2_buff_proj(i_na,i,j)/(N1*N2*N3)**2
                        avrgvx(j,k) = avrgvx(j,k) + uavrg_buff_proj(1,i_na,i,j)/(N1*N2*N3)
                        avrgvy(j,k) = avrgvy(j,k) + uavrg_buff_proj(2,i_na,i,j)/(N1*N2*N3)
                        avrgvz(j,k) = avrgvz(j,k) + uavrg_buff_proj(3,i_na,i,j)/(N1*N2*N3)
                        avrgv2t(j,k) = avrgv2t(j,k) + uavrg2_buff_proj(i_na,i,j)/(N1*N2*N3)**2
                     end if
                  end do

                  !Type k projected mean magnetization for ensemble j
                  avrgu(j,k) = sqrt(avrgux(j,k)**2+avrguy(j,k)**2+avrguz(j,k)**2)
                  avrgu2(j,k) = avrgu2t(j,k)
                  avrgv(j,k) = sqrt(avrgvx(j,k)**2+avrgvy(j,k)**2+avrgvz(j,k)**2)
                  avrgv2(j,k) = avrgv2t(j,k)
               end do
            end do
         end if
         if(do_proj_avrg=='A') then
            do j=1, Mensemble
               do k=1, NA
                  avrgux(j,k)=uavrg_buff_proj(1,k,i,j)/(N1*N2*N3)
                  avrguy(j,k)=uavrg_buff_proj(2,k,i,j)/(N1*N2*N3)
                  avrguz(j,k)=uavrg_buff_proj(3,k,i,j)/(N1*N2*N3)
                  avrgu2t(j,k) = uavrg2_buff_proj(i_na,i,j)/(N1*N2*N3)**2
                  !Site k projected mean magnetization for ensemble j
                  avrgu(j,k) = sqrt(avrgux(j,k)**2+avrguy(j,k)**2+avrguz(j,k)**2)
                  avrgu2(j,k) = avrgu2t(j,k)
                  avrgvx(j,k)=uavrg_buff_proj(1,k,i,j)/(N1*N2*N3)
                  avrgvy(j,k)=uavrg_buff_proj(2,k,i,j)/(N1*N2*N3)
                  avrgvz(j,k)=uavrg_buff_proj(3,k,i,j)/(N1*N2*N3)
                  avrgv2t(j,k) = uavrg2_buff_proj(i_na,i,j)/(N1*N2*N3)**2
                  !Site k projected mean magnetization for ensemble j
                  avrgv(j,k) = sqrt(avrgvx(j,k)**2+avrgvy(j,k)**2+avrgvz(j,k)**2)
                  avrgv2(j,k) = avrgv2t(j,k)
               end do
            end do
         end if

         !Start with calculation of mean over the ensembles
         if(do_proj_avrg=='Y') then
            nproj=NT
         else if(do_proj_avrg=='A') then
            nproj=NA
         else
            nproj=NA
         end if

         do k=1, nproj

            !mean magnetisation over the ensembles
            avrgum(k)  =0_dblprec
            avrgum2(k) =0_dblprec
            avrgvm(k)  =0_dblprec
            avrgvm2(k) =0_dblprec

            !standard deviation of the mean magnetisation over the ensembles
            avrgus(k)=0_dblprec
            avrgvs(k)=0_dblprec

            avrguxm(k)=0_dblprec
            avrguym(k)=0_dblprec
            avrguzm(k)=0_dblprec
            avrgvxm(k)=0_dblprec
            avrgvym(k)=0_dblprec
            avrgvzm(k)=0_dblprec
            do j=1, Mensemble
               avrguxm(k) = avrguxm(k) + avrgux(j,k)
               avrguym(k) = avrguym(k) + avrguy(j,k)
               avrguzm(k) = avrguzm(k) + avrguz(j,k)
               avrgum(k) = avrgum(k) + avrgu(j,k)
               avrgus(k) = avrgus(k) + avrgu(j,k)**2
               avrgum2(k) = avrgum2(k) + avrgu2(j,k)/nproj
               avrgvxm(k) = avrgvxm(k) + avrgvx(j,k)
               avrgvym(k) = avrgvym(k) + avrgvy(j,k)
               avrgvzm(k) = avrgvzm(k) + avrgvz(j,k)
               avrgvm(k) = avrgvm(k) + avrgv(j,k)
               avrgvs(k) = avrgvs(k) + avrgv(j,k)**2
               avrgvm2(k) = avrgvm2(k) + avrgv2(j,k)/nproj
            end do
            avrguxm(k)=avrguxm(k)/Mensemble
            avrguym(k)=avrguym(k)/Mensemble
            avrguzm(k)=avrguzm(k)/Mensemble
            avrgum(k)=avrgum(k)/Mensemble
            avrgus(k)=avrgus(k)/Mensemble
            avrgus(k)=avrgus(k) - avrgum(k)**2
            avrgum2(k)=avrgum2(k)/Mensemble
            avrgvxm(k)=avrgvxm(k)/Mensemble
            avrgvym(k)=avrgvym(k)/Mensemble
            avrgvzm(k)=avrgvzm(k)/Mensemble
            avrgvm(k)=avrgvm(k)/Mensemble
            avrgvs(k)=avrgvs(k)/Mensemble
            avrgvs(k)=avrgvs(k) - avrgvm(k)**2
            avrgvm2(k)=avrgvm2(k)/Mensemble

            !Filter out negative errors
            if(avrgus(k)<0) then
               avrgus(k)=0
            else
               avrgus(k)=sqrt(avrgus(k))
            end if
            if(avrgvs(k)<0) then
               avrgvs(k)=0
            else
               avrgvs(k)=sqrt(avrgvs(k))
            end if

            !Write results to file
            write (ofileno,10004) int(indxb_uavrg(i)), k, &
               avrgum(k), avrgus(k), avrguxm(k), avrguym(k), avrguzm(k)
         end do
      end do
      close(ofileno)
      return

      write (*,*) "Error writing the projected averages file"
      10004 format (i8,i8,21es16.8)
      10005 format (es16.4,i8,21es16.8)

   end subroutine prn_proj_uavrg


   !> Print chemical projected average magnetizations
   subroutine prn_projch_uavrg(simid, Mensemble, N1, N2, N3, Nchmax)
      !
      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: Nchmax    !< Max number of chemical components on each site in cell

      !.. Scalar variables
      integer :: i,j,k
      character(len=30) :: filn
      real(dblprec), dimension(Mensemble,Nchmax) :: avrgux, avrguy, avrguz, avrgu, avrguu
      real(dblprec), dimension(Nchmax) :: avrgum, avrgus, avrguxm, avrguym, avrguzm, avrguum
      real(dblprec), dimension(Mensemble,Nchmax) :: avrgvx, avrgvy, avrgvz, avrgv, avrgvv
      real(dblprec), dimension(Nchmax) :: avrgvm, avrgvs, avrgvxm, avrgvym, avrgvzm, avrgvvm

      !.. Executable statements
      write (filn,'(''lattprojchavrgs.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")

      do i=1, bcount_uavrg
         do j=1, Mensemble
            do k=1, Nchmax
               avrgux(j,k) = uavrg_buff_projch(1,k,i,j)/(N1*N2*N3)
               avrguy(j,k) = uavrg_buff_projch(2,k,i,j)/(N1*N2*N3)
               avrguz(j,k) = uavrg_buff_projch(3,k,i,j)/(N1*N2*N3)
               avrgu(j,k) = uavrg_buff_projch(4,k,i,j)/(N1*N2*N3)
               avrgvx(j,k) = uavrg_buff_projch(1,k,i,j)/(N1*N2*N3)
               avrgvy(j,k) = uavrg_buff_projch(2,k,i,j)/(N1*N2*N3)
               avrgvz(j,k) = uavrg_buff_projch(3,k,i,j)/(N1*N2*N3)
               avrgv(j,k) = uavrg_buff_projch(4,k,i,j)/(N1*N2*N3)
               !Calculate mean sublattice displacements and velocities for each ensemble
               ! average in ensemble j of sublattice k
               avrguu(j,k) = sqrt(avrgux(j,k)**2+avrguy(j,k)**2+avrguz(j,k)**2)
               avrgvv(j,k) = sqrt(avrgvx(j,k)**2+avrgvy(j,k)**2+avrgvz(j,k)**2)
            end do
         end do

         !Start with calculation of mean over the ensembles
         do k=1, Nchmax

            !mean magnetisation over the ensembles
            avrguum(k)=0_dblprec
            avrgvvm(k)=0_dblprec

            !standard deviation of the mean magnetisation over the ensembles
            avrgus(k)=0_dblprec
            avrguxm(k)=0_dblprec
            avrguym(k)=0_dblprec
            avrguzm(k)=0_dblprec
            avrgum(k)=0_dblprec
            avrgvs(k)=0_dblprec
            avrgvxm(k)=0_dblprec
            avrgvym(k)=0_dblprec
            avrgvzm(k)=0_dblprec
            avrgvm(k)=0_dblprec

            ! CHECK THESE!
            do j=1, Mensemble
               avrguxm(k) = avrguxm(k) + avrgux(j,k)
               avrguym(k) = avrguym(k) + avrguy(j,k)
               avrguzm(k) = avrguzm(k) + avrguz(j,k)
               avrgum(k) = avrgum(k) + avrgu(j,k)
               avrguum(k) = avrguum(k) + avrguu(j,k)
               avrgus(k) = avrgus(k) + avrgu(j,k)**2
               avrgvxm(k) = avrgvxm(k) + avrgvx(j,k)
               avrgvym(k) = avrgvym(k) + avrgvy(j,k)
               avrgvzm(k) = avrgvzm(k) + avrgvz(j,k)
               avrgvm(k) = avrgvm(k) + avrgv(j,k)
               avrgvvm(k) = avrgvvm(k) + avrgvv(j,k)
               avrgvs(k) = avrgvs(k) + avrgv(j,k)**2
            end do

            avrguxm(k)=avrguxm(k)/Mensemble
            avrguym(k)=avrguym(k)/Mensemble
            avrguzm(k)=avrguzm(k)/Mensemble
            avrgum(k)=avrgum(k)/Mensemble
            avrgus(k)=avrgus(k)/Mensemble
            avrguum(k)=avrguum(k)/Mensemble
            avrgus(k)=avrgus(k) - avrgum(k)**2
            avrgvxm(k)=avrgvxm(k)/Mensemble
            avrgvym(k)=avrgvym(k)/Mensemble
            avrgvzm(k)=avrgvzm(k)/Mensemble
            avrgvm(k)=avrgvm(k)/Mensemble
            avrgvs(k)=avrgvs(k)/Mensemble
            avrgvvm(k)=avrgvvm(k)/Mensemble
            avrgvs(k)=avrgvs(k) - avrgvm(k)**2

            !Filter out negative errors
            if(avrgus(k)<0) then
               avrgus(k)=0
            else
               avrgus(k)=sqrt(avrgus(k))
            end if
            if(avrgvs(k)<0) then
               avrgvs(k)=0
            else
               avrgvs(k)=sqrt(avrgvs(k))
            end if

            !Write results to file
            write (ofileno,10004) int(indxb_uavrg(i)), k, &
               avrgum(k), avrgus(k), avrguxm(k), avrguym(k), avrguzm(k),avrguum(k), &
               avrgvm(k), avrgvs(k), avrgvxm(k), avrgvym(k), avrgvzm(k),avrgvvm(k)
         end do
      end do
      close(ofileno)
      return

      write (*,*) "Error writing the chemical projected averages file"
      10004 format (i8,i8,21es16.8)
      10005 format (es16.4,i8,21es16.8)

   end subroutine prn_projch_uavrg

   !> Calculate ion temperature directly (adapted from buffer_uavrg)
   function f_iontemp(Natom, Mensemble, mion, vvec) result(iontemp)

      use Constants

      !.. Implicit declarations
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble   !< Number of ensembles
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mion   !< Ionic mass
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: vvec   !< Current ionic velocity

      !.. Scalar variables
      integer :: i, k

      real(dblprec) ::  v2, kinenergy
      real(dblprec) :: iontemp

      !.. Executable statements

      !.. Sum over atoms
      do k=1,Mensemble
         do i=1, Natom

            v2 = vvec(1,i,k)*vvec(1,i,k) + vvec(2,i,k)*vvec(2,i,k) + vvec(3,i,k)*vvec(3,i,k)

            ! Convert the kinetic energy to mRyd
            kinenergy = kinenergy + 0.5_dblprec * amu * mion(i,1) * angstrom**2 * v2 / mry

         end do
      end do

      kinenergy = kinenergy / Mensemble / Natom
      iontemp = kinenergy * mry * 2.0_dblprec / (3.0_dblprec * k_bolt)

   end function f_iontemp


end module prn_latticeaverages
