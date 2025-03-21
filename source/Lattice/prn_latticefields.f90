!> Data structures for printing the lattice DOF effective and thermal fields
module prn_latticefields

   use Profiling 
   use Constants
   use Parameters

   implicit none

   ! Input parameters to be read
   integer :: eeff_step               !< Interval between consecutive prints of the total effective field
   integer :: eeff_buff               !< Buffer size for the total field
   integer :: einteff_step               !< Interval between consecutive prints of the internal effective field
   integer :: einteff_buff               !< Buffer size for the internal field
   integer :: ethermfield_step         !< Interval between thermal field trajectories
   integer :: ethermfield_buff         !< Buffer size for the stochastic field
   character(len=1) :: do_prn_eeff    !< Flag governing file output of total effective fields (Y/N)
   character(len=1) :: do_prn_einteff    !< Flag governing file output of internal effective fields (Y/N)
   character(len=1) :: do_ethermfield  !< Thermal fields trajectory

   ! Local variables for buffering and indexing of fields
   integer :: bcount_eeff    !< Counter of buffer for total field
   integer :: bcount_einteff    !< Counter of buffer for internal field
   integer :: bcount_etherm   !< Counter of buffer for stochastic field
   real(dblprec), dimension(:), allocatable :: indxb_eeff          !< Step counter for total field
   real(dblprec), dimension(:), allocatable :: indxb_einteff       !< Step counter for internal field
   real(dblprec), dimension(:), allocatable :: indxb_etherm        !< Step counter for stochastic field
   real(dblprec), dimension(:,:,:,:), allocatable :: eeffb         !< Buffer the site dependent total field
   real(dblprec), dimension(:,:,:,:), allocatable :: einteffb      !< Buffer the site dependent internal field
   real(dblprec), dimension(:,:,:,:), allocatable :: etherm_fieldb !< Buffer the site dependent stochastic field

   public

   private :: bcount_eeff, bcount_einteff, bcount_etherm, indxb_eeff, indxb_einteff, indxb_etherm, eeffb, einteffb, etherm_fieldb


contains


   !> Wrapper routine to print the thermal and effective fields
   subroutine print_lattfields(simid, Natom, Mensemble, mstep, sstep, eeff, eeff1, eeff3, ethermal_field)

      implicit none

      character(len=8), intent(in) :: simid             !< Name of simulation 
      integer, intent(in) :: Natom                      !< Number of atoms in the system
      integer, intent(in) :: Mensemble                  !< Number of ensembles
      integer, intent(in) :: mstep                      !< Current simulation step
      integer, intent(in) :: sstep                      !< Current simulation step in logarithmic scale
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff    !< Current site dependent total effective e-field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff1   !< Current site dependent internal effective e-field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff3   !< Current site dependent internal effective e-field 
      !from mixed spin-lattice Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: ethermal_field !< Current site dependent stochastic e-field

      ! Stochastic thermal field
      if (do_ethermfield=='Y') then

         if (mod(sstep-1,ethermfield_step)==0) then
            ! Write step to buffer
            call buffer_ethermfield(Natom, Mensemble, mstep-1, ethermal_field, bcount_etherm)

            if (bcount_etherm==ethermfield_buff) then
               ! write buffer to file
               call prn_ethermfields(simid, Natom, Mensemble)
               bcount_etherm=1
            else
               bcount_etherm=bcount_etherm+1
            endif

         endif
      endif

      ! Total site dependent field
      if (do_prn_eeff=='Y') then

         if (mod(sstep-1,eeff_step)==0) then
            !write(*,*) 'do_prn_eeff sstep ', sstep
            ! Write step to buffer
            call buffer_totalefield(Natom, Mensemble, mstep-1, eeff, ethermal_field, bcount_eeff)
            if (bcount_eeff==eeff_buff) then
               ! write buffer to file
               call prn_totalefields(simid, Natom, Mensemble)
               bcount_eeff=1
            else
               bcount_eeff=bcount_eeff+1
            endif

         endif
      endif

      ! Internal site dependent field
      if (do_prn_einteff=='Y') then

         if (mod(sstep-1,einteff_step)==0) then
            !write(*,*) 'do_prn_einteff sstep ', sstep
            ! Write step to buffer
            call buffer_internalefield(Natom, Mensemble, mstep-1, eeff1, eeff3, bcount_einteff)
            if (bcount_einteff==einteff_buff) then
               ! write buffer to file
               call prn_internalefields(simid, Natom, Mensemble)
               bcount_einteff=1
            else
               bcount_einteff=bcount_einteff+1
            endif

         endif
      endif

   end subroutine print_lattfields


   !> Flush the field measurements, i.e. print to file the fields in the last time step
   subroutine flush_lattfields(simid, Natom, Mensemble)

      implicit none

      character(len=8), intent(in) :: simid    !< Name of simulation 
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 

      if (do_ethermfield=='Y') then
         ! Write buffer to file
         bcount_etherm=bcount_etherm-1
         call prn_ethermfields(simid, Natom, Mensemble)
      endif

      if (do_prn_eeff=='Y') then
         ! Write buffer to file
         bcount_eeff=bcount_eeff-1
         call prn_totalefields(simid, Natom, Mensemble)
      endif

      if (do_prn_einteff=='Y') then
         ! Write buffer to file
         bcount_einteff=bcount_einteff-1
         call prn_internalefields(simid, Natom, Mensemble)
      endif

   end subroutine flush_lattfields


   !> Initialization of the variables for lattice field printing with their default variables
   subroutine lattfields_prn_init() 

      implicit none

      do_ethermfield    = 'N'
      do_prn_eeff       = 'N'
      do_prn_einteff    = 'N'
      ethermfield_step  = 100
      ethermfield_buff  = 10
      eeff_step         = 100
      eeff_buff         = 10
      einteff_step      = 100
      einteff_buff      = 10

   end subroutine lattfields_prn_init


   !> Routine for allocation and initialization of the arrays for lattice field printing
   subroutine allocate_lattfield_print(Natom, Mensemble, flag)

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: flag      !< Allocate or deallocate (1/-1)

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         bcount_eeff=1
         bcount_einteff=1
         bcount_etherm=1

         if (do_ethermfield=='Y') then
            allocate(etherm_fieldb(3,Natom,ethermfield_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(etherm_fieldb))*kind(etherm_fieldb),'etherm_fieldb','allocate_measurements')
            allocate(indxb_etherm(ethermfield_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_etherm))*kind(indxb_etherm),'indxb_etherm','allocate_measurements')
         endif

         if (do_prn_eeff=='Y') then
            allocate(eeffb(3,Natom,eeff_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(eeffb))*kind(eeffb),'eeffb','allocate_measurements')
            allocate(indxb_eeff(eeff_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_eeff))*kind(indxb_eeff),'indxb_eeff','allocate_measurements')
         endif

         if (do_prn_einteff=='Y') then
            allocate(einteffb(6,Natom,einteff_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(einteffb))*kind(einteffb),'einteffb','allocate_measurements')
            allocate(indxb_einteff(einteff_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_einteff))*kind(indxb_einteff),'indxb_einteff','allocate_measurements')
         endif

      else 

         if (do_ethermfield=='Y') then
            i_all=-product(shape(etherm_fieldb))*kind(etherm_fieldb)
            deallocate(etherm_fieldb,stat=i_stat)
            call memocc(i_stat,i_all,'etherm_fieldb','allocate_measurements')
            i_all=-product(shape(indxb_etherm))*kind(indxb_etherm)
            deallocate(indxb_etherm,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_etherm','allocate_measurements')
         endif

         if (do_prn_eeff=='Y') then
            i_all=-product(shape(eeffb))*kind(eeffb)
            deallocate(eeffb,stat=i_stat)
            call memocc(i_stat,i_all,'eeffb','allocate_measurements')
            i_all=-product(shape(indxb_eeff))*kind(indxb_eeff)
            deallocate(indxb_eeff,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_eeff','allocate_measurements')
         endif

         if (do_prn_einteff=='Y') then
            i_all=-product(shape(einteffb))*kind(einteffb)
            deallocate(einteffb,stat=i_stat)
            call memocc(i_stat,i_all,'einteffb','allocate_measurements')
            i_all=-product(shape(indxb_einteff))*kind(indxb_einteff)
            deallocate(indxb_einteff,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_einteff','allocate_measurements')
         endif

      endif

   end subroutine allocate_lattfield_print


   !> Buffer site dependent thermal field
   subroutine buffer_ethermfield(Natom, Mensemble, mstep, ethermal_field, bcount_etherm)
      ! 

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles 
      integer, intent(in) :: mstep !< Current simulation step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: ethermal_field   !< Current thermal field
      integer, intent(in) :: bcount_etherm  !< Counter of buffer for ethermal fields

      !.. Local scalar variables
      integer :: i,k

      do k=1, Mensemble
         do i=1, Natom
            etherm_fieldb(1:3,i,bcount_etherm,k)=ethermal_field(1:3,i,k)
         end do
      end do



   end subroutine buffer_ethermfield


   !> Buffer site dependent total field
   subroutine buffer_totalefield(Natom, Mensemble, mstep, eeff, ethermal_field, bcount_eeff)
      ! 

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: mstep !< Current simulation step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff          !< Current effective field from the hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: ethermal_field !< Current thermal field
      integer, intent(in) :: bcount_eeff   !< Counter of buffer for total effective field

      !.. Local scalar variables
      integer :: i,k

      do k=1, Mensemble
         do i=1, Natom
            eeffb(1:3,i,bcount_eeff,k)=ethermal_field(1:3,i,k)+eeff(1:3,i,k)
         end do
      end do

      indxb_eeff(bcount_eeff) = mstep

   end subroutine buffer_totalefield


   !> Buffer site dependent internal field
   subroutine buffer_internalefield(Natom, Mensemble, mstep, eeff1, eeff3, bcount_eeff)
      ! 

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: mstep !< Current simulation step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff1   !< Current site dependent internal effective e-field                                    
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff3   !< Current site dependent internal effective e-field 
      !from mixed spin-lattice Hamiltonian
      integer, intent(in) :: bcount_eeff   !< Counter of buffer for total effective field

      !.. Local scalar variables
      integer :: i,k

      do k=1, Mensemble
         do i=1, Natom
            einteffb(1:3,i,bcount_einteff,k) = eeff1(1:3,i,k)
            einteffb(4:6,i,bcount_einteff,k) = eeff3(1:3,i,k)
         end do
      end do

      indxb_einteff(bcount_einteff) = mstep

   end subroutine buffer_internalefield


   !> Print total effective field 
   subroutine prn_totalefields(simid, Natom, Mensemble)
      !
      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Name of the simulation
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 

      ! Local variables
      integer :: i,j,k
      character(len=30) :: filn

      ! Print thermal fields to output file if specified
      write(filn,'(''eefftot.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position = 'APPEND')

      ! Write header to output files for first iteration
      if(abs(indxb_eeff (1))<=0.0e0_dblprec) then
         write (ofileno,'(a)') "# Iter.     Site     Replica      E_x             E_y             E_z             E"
      end if

      do k=1, bcount_eeff
         do j=1,Mensemble
            do i=1,Natom
               write (ofileno,120) int(indxb_eeff(k)), i, j, eeffb(1:3,i,k,j),norm2(eeffb(:,i,k,j))
            end do
         end do
      enddo
      close(ofileno)

      return

      write(*,*) 'Error writing the total field file'

      120 format(i8,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)
      121 format(es16.4,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)

   end subroutine prn_totalefields


   !> Print internal effective field 
   subroutine prn_internalefields(simid, Natom, Mensemble)
      !
      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Name of the simulation
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 

      ! Local variables
      integer :: i,j,k
      character(len=30) :: filn, filn2

      real(dblprec), dimension(6,Mensemble) :: einteffavrg

      ! Print internal fields to output file if specified
      write(filn,'(''eintefftot.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position = 'APPEND')
      write(filn2,'(''einteffavrg.'',a,''.out'')') trim(simid)
      open(ofileno2,file=filn2, position = 'APPEND')

      ! Write header to output files for first iteration
      if(abs(indxb_eeff (1))<=0.0e0_dblprec) then
         write (ofileno,'(a)') "# Iter.     Site     Replica      E_LD_x      E_LD_y      E_LD_z      E_LD      "// &
            & "E_SLD_x     E_SLD_y     E_SLD_z     E_SLD"
         write (ofileno2,'(a)') "# Iter.     Site     Replica      E_LD_x      E_LD_y      E_LD_z      E_LD      "// &
            & "E_SLD_x     E_SLD_y     E_SLD_z     E_SLD"
      end if

      do k=1, bcount_einteff
         do j=1,Mensemble
            einteffavrg = 0.0_dblprec
            do i=1,Natom
               einteffavrg(1:6,j) = einteffavrg(1:6,j) + einteffb(1:6,i,k,j)
               write (ofileno,120) int(indxb_einteff(k)), i, j, einteffb(1:3,i,k,j), norm2(einteffb(1:3,i,k,j)), &
                  einteffb(4:6,i,k,j), norm2(einteffb(4:6,i,k,j))
            end do
            einteffavrg(1:6,j) = einteffavrg(1:6,j) / Natom 
            write (ofileno2,120) int(indxb_einteff(k)), 1, j, einteffavrg(1:3,j), norm2(einteffavrg(1:3,j)), &
               einteffavrg(4:6,j), norm2(einteffavrg(4:6,j))
         end do
      enddo
      close(ofileno)
      close(ofileno2)

      return

      write(*,*) 'Error writing the internal field file'

      120 format(i8,i8,i8,8x,8es12.4)
      121 format(es12.4,i8,i8,8x,8es12.4)

   end subroutine prn_internalefields


   !> Print thermal field
   subroutine prn_ethermfields(simid, Natom, Mensemble)
      !

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      character(len=8), intent(in) :: simid !< Simulation name

      character(len=30) :: filn
      ! Local variables
      integer :: i,j,k

      ! Print thermal fields to output file if specified
      write(filn,'(''ethermfields.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position = 'APPEND')
      do k=1, bcount_etherm
         do j=1,Mensemble
            do i=1,Natom
               write (ofileno,110) int(indxb_etherm(k)), i, j, etherm_fieldb(:,i,k,j),norm2(etherm_fieldb(:,i,k,j))
            end do
         end do
      enddo
      close(ofileno)

      return

      write(*,*) 'Error writting the thermal field file'

      110 format(i8,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)
      111 format(es16.4,i8,i8,8x,es16.8,es16.8,es16.8,es16.8)

   end subroutine prn_ethermfields


end module prn_latticefields
