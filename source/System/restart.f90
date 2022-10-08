!------------------------------------------------------------------------------------
! MODULE: Restart
!> @brief Routine for handling the input/output of magnetic configurations
!> @author Jonathan Chico
!> @copyright
!> GNU Public License.
!------------------------------------------------------------------------------------
module Restart
   use Parameters
   use Profiling

   implicit none

   interface prn_mag_conf
      procedure prn_mag_conf_iter
      procedure prn_mag_conf_time
   end interface prn_mag_conf


   interface read_mag_conf
      procedure read_mag_conf_std
      procedure read_mag_conf_single
      procedure read_mag_conf_double
   end interface read_mag_conf

   public

contains

#ifdef USE_OVF
   !---------------------------------------------------------------------------------
   ! > @brief Print the magnetic configurations in ovf format
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine prn_mag_conf_ovf(Natom,Mensemble,type,simid,emomM)
      
      use ovf

      implicit none

      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble  !< Number of ensembles 
      character(len=1), intent(in) :: type
      character(len=8), intent(in) :: simid
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment
      ! .. Local variables
      integer :: ii,jj
      integer :: i_stat, i_all
      character(len=40) file_name
      real(kind=8), dimension(:,:), allocatable :: emom_ovf
      ! OVF data structures
      integer :: success
      type(ovf_file)      :: file
      type(ovf_segment)   :: segment
      
      if (type=='R') then
      ! Write the name of the restartfile and the position of the writing of file
         write (file_name,'(''restart.'',a,''.out'')') trim(simid)
      else if (type=='M') then
      ! Write the name of the moment and the position of the writing of file
         write (file_name,'(''moment.'',a,''.out'')') trim(simid)
      endif
      ! Initialize segment
      call segment%initialize()
      ! Open the file
      call file%open_file(file_name)
      segment%N = Natom
      segment%ValueDim = 3
      ! Allocate temporary array
      allocate(emom_ovf(3,segment%N),stat=i_stat)
      call memocc(i_stat,product(shape(emom_ovf))*kind(emom_ovf),'emom_ovf','prn_mag_conf_ovf')
      emom_ovf=0.0

      do jj=1, Mensemble
         do ii=1,Natom
            emom_ovf(:,ii)=emomM(:,ii,jj)
         enddo
      
         success = file%append_segment(segment, emom_ovf, OVF_FORMAT_TEXT)
         if ( success == OVF_ERROR) then
            write (*,*) "test write_segment did not work. Message: ", file%latest_message
            stop 1
         endif
      enddo

      ! Deallocate the temporary array
      i_all=-product(shape(emom_ovf))*kind(emom_ovf)
      deallocate(emom_ovf,stat=i_stat)
      call memocc(i_stat,i_all,'emom_ovf','prn_mag_conf_ovf')

   end subroutine prn_mag_conf_ovf

   !---------------------------------------------------------------------------------
   ! > @brief Read the magnetic configurations in ovf format
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine read_mag_conf_ovf(Natom,Mensemble,restartfile,mmom,emom,emomM)

      use ovf
      
      implicit none

      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble  !< Number of ensembles 
      character(len=35), intent(inout) :: restartfile !< File containing restart information
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom  !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment
      ! .. Local variables
      integer :: ii,jj
      integer :: i_stat, i_all
      real(kind=8), dimension(:,:), allocatable :: emom_ovf
      ! OVF data structures
      integer :: success
      type(ovf_file)      :: file
      type(ovf_segment)   :: segment
      
      ! Allocate a temporary array where the ovf data is stored
      allocate(emom_ovf(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(emom_ovf))*kind(emom_ovf),'emom_ovf','read_mag_conf_ovf')
      emom_ovf=0.0
      ! Initialize the segment
      call segment%initialize()
      ! Open the file
      call file%open_file(restartfile)
      ! Find some general information on the file
      if (success == OVF_ERROR) then
         ! Check if the file exists
         write(*,'(2x,a,i8)') "Found      = ", file%found
         ! Check if the file has the correct format
         write(*,'(2x,a,i8)') "is_ovf     = ", file%is_ovf
         ! Check how many ensembles there are in a file
         write(*,'(2x,a,i8)') "n_segments = ", file%n_segments
         stop
      endif
      ! If the number of segments is larger or equal to the number of
      ! ensembles one can do something about it
      if (file%n_segments>=Mensemble) then
         ! If there are more segments than ensembles give a warning
         if (file%n_segments>Mensemble) then
            write(*,'(2x,a,i6,a)') 'WARNING: Reading only the first ', Mensemble,' entries'
         endif
         ! Loop over the number of ensembles
         do jj=1,Mensemble
            ! Read the header of each segment to see if everything is okay
            success = file%read_segment_header(segment)
            if ( success .ne. OVF_OK) then
               write (*,*) "read_segment_header did not work. Message: ",           &
               file%latest_message
               stop 1
            endif
            ! Read the actual data and store it in the temporary array
            success = file%read_segment_data(segment,emom_ovf,jj)
            if ( success .ne. OVF_OK) then
               write (*,*) "read_segment_data on emom_ovf did not work. Message: ", &
               file%latest_message
               stop 1
            else
               do ii=1,Natom
                  emomM(:,ii,jj)=emom_ovf(:,ii)
                  mmom(ii,jj)=norm2(emomM(:,ii,jj))
                  emom(:,ii,jj)=emomM(:,ii,jj)/mmom(ii,jj)
               enddo
            endif
         enddo
         write (*,'(a)') " done"
      else
         write(*,'(2x,a)') 'Too few segments'
         stop
      endif
      ! Deallocate the temporary array
      i_all=-product(shape(emom_ovf))*kind(emom_ovf)
      deallocate(emom_ovf,stat=i_stat)
      call memocc(i_stat,i_all,'emom_ovf','read_mag_conf_ovf')

   end  subroutine read_mag_conf_ovf
#endif
   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_mag_conf_iter
   !> @brief Prints a given magnetic configuration for either a restartfile or a momentfile
   !> @details Prints a magnetic configuration, the objective is to make all the types
   !> of printing honogeneous, of that way restartfiles, momentfiles and GNEB files
   !> would have all the same structure.
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine prn_mag_conf_iter(Natom,mstep,Mensemble,type,simid,mmom,emom,suffix,  &
      mode)

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: mstep     !< Current simulation step
      integer, intent(in) :: Mensemble !< Number of ensembles 
      character(len=1), intent(in) :: type   !< type to see whether the file is a restartfile or a momentfile
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid  !< Name of simulation
      character(len=*), intent(in) :: suffix !< Suffix to be appended to the files
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector

      !.. Local variables
      integer :: ii,jj
      character(len=40) filn,fil_pos

      if (type=='R') then
      ! Write the name of the restartfile and the position of the writing of file
         write (filn,'(''restart'',a,''.'',a,''.out'')') suffix,trim(simid)
         fil_pos="rewind"
         open(ofileno,file=filn,position=trim(fil_pos))
         write(ofileno,'(a)') repeat("#",80)
         write(ofileno,'(a,1x,a)') "# File type:", type
         write(ofileno,'(a,1x,a)') "# Simulation type:", mode
         write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
         write(ofileno,'(a,1x,i8)')"# Number of ensembles: ", Mensemble
         write(ofileno,'(a)') repeat("#",80)
         write(ofileno,'(a8,a,a8,a16,a16,a16,a16)') "#iter","ens","iatom","|Mom|","M_x","M_y","M_z"
      else if (type=='M') then
      ! Write the name of the moment and the position of the writing of file
         write (filn,'(''moment'',a,''.'',a,''.out'')') suffix,trim(simid)
         fil_pos="append"
         if (mstep==0) then
            fil_pos="rewind"
            open(ofileno,file=filn,position=trim(fil_pos),access='stream',form='formatted')
            write(ofileno,'(a)') repeat("#",80)
            write(ofileno,'(a,1x,a)') "# File type:", type
            write(ofileno,'(a,1x,a)') "# Simulation type:", mode
            write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
            write(ofileno,'(a,1x,i8)')"# Number of ensembles: ", Mensemble
            write(ofileno,'(a)') repeat("#",80)
            write(ofileno,'(a8,a,a8,a16,a16,a16,a16)') "#iter","ens","iatom","|Mom|","M_x","M_y","M_z"
         else
            fil_pos="append"
            open(ofileno,file=filn,position=trim(fil_pos),access='stream',form='formatted')
         endif
      endif

      do ii=1, Mensemble
         do jj=1, Natom
            write(ofileno,10003) mstep,ii,jj, mmom(jj,ii), emom(1,jj,ii),emom(2,jj,ii),emom(3,jj,ii)
         enddo
      enddo
      close(ofileno)
      10003 format(i8,i8,i8,2x,es16.8,es16.8,es16.8,es16.8)

   end subroutine prn_mag_conf_iter

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_mag_conf_time
   !> @brief Prints a given magnetic configuration for either a restartfile or a momentfile
   !> @details Prints a magnetic configuration, the objective is to make all the types
   !> of printing honogeneous, of that way restartfiles, momentfiles and GNEB files
   !> would have all the same structure.
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine prn_mag_conf_time(Natom,time,Mensemble,type,simid,mmom,emom,suffix,   &
      mode)

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      real(dblprec), intent(in) :: time    !< Current simulation time
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=1), intent(in) :: type   !< type to see whether the file is a restartfile or a momentfile
      character(len=8), intent(in) :: simid  !< Name of simulation
      character(len=*), intent(in) :: suffix !< Suffix to be appended to the files
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector

      !.. Local variables
      integer :: ii,jj
      character(len=40) filn,fil_pos

      if (type=='R') then
      ! Write the name of the restartfile and the position of the writing of file
         write (filn,'(''restart'',a,''.'',a,''.out'')') suffix,trim(simid)
         fil_pos="rewind"
         open(ofileno,file=filn,position=trim(fil_pos),access='stream',form='formatted')
         write(ofileno,'(a)') repeat("#",80)
         write(ofileno,'(a,1x,a)') "# File type:", type
         write(ofileno,'(a,1x,a)') "# Simulation type:", mode
         write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
         write(ofileno,'(a,1x,i8)')"# Number of ensembles: ", Mensemble
         write(ofileno,'(a)') repeat("#",80)
         write(ofileno,'(a8,a,a8,a16,a16,a16,a16)') "#Time [s]","ens","iatom","|Mom|","M_x","M_y","M_z"
      else if (type=='M') then
      ! Write the name of the moment and the position of the writing of file
         write (filn,'(''moment'',a,''.'',a,''.out'')') suffix,trim(simid)
         if (time==0) then
            fil_pos="rewind"
            open(ofileno,file=filn,position=trim(fil_pos),access='stream',form='formatted')
            write(ofileno,'(a)') repeat("#",80)
            write(ofileno,'(a,1x,a)') "# File type:", type
            write(ofileno,'(a,1x,a)') "# Simulation type:", mode
            write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
            write(ofileno,'(a,1x,i8)')"# Number of ensembles: ", Mensemble
            write(ofileno,'(a)') repeat("#",80)
            write(ofileno,'(a8,a,a8,a16,a16,a16,a16)') "#Time [s]","ens","iatom","|Mom|","M_x","M_y","M_z"
         else
            fil_pos="append"
            open(ofileno,file=filn,position=trim(fil_pos),access='stream',form='formatted')
         endif
      endif

      do ii=1, Mensemble
         do jj=1, Natom
            write(ofileno,10003) time,ii,jj,mmom(jj,ii),emom(1,jj,ii),emom(2,jj,ii),emom(3,jj,ii)
         enddo
      enddo
      close(ofileno)

      10003 format(es16.8,i8,i8,2x,es16.8,es16.8,es16.8,es16.8)

   end subroutine prn_mag_conf_time

   !---------------------------------------------------------------------------------
   !> @brief Routine handling the read-in of magnetic configurations
   !---------------------------------------------------------------------------------
   subroutine read_mag_conf_std(Natom,Mensemble,do_mom_legacy,rstep,restartfile,    &
      mmom,emom,emomM)

      use math_functions, only: f_normalize_vec

      implicit none

      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=1), intent(in) :: do_mom_legacy  !< Flag to print/read moments in legacy output
      integer, intent(out) :: rstep !< Starting simulation step
      character(len=35), intent(inout) :: restartfile !< File containing restart information
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom  !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! .. Local variables
      integer :: dummy, ii, jj,ios
      logical :: exists

      inquire(file=restartfile,exist=exists)
      if(exists) then
         open(ifileno,iostat=ios, file=restartfile, status="old")
         if (do_mom_legacy.ne.'Y') then
            !------------------------------------------------------------------------
            ! If in the new format read the first 7 lines to skip the header
            !------------------------------------------------------------------------
            do ii=1,7
               read(ifileno,*)
            enddo
            !------------------------------------------------------------------------
            ! Read the magnetic configurations
            !------------------------------------------------------------------------
            do ii=1, Mensemble
               do jj=1, Natom
                  read(ifileno,*) rstep, dummy, dummy, mmom(jj,ii), emom(1,jj,ii),  &
                     emom(2,jj,ii), emom(3,jj,ii)
                  emom(:,jj,ii)=f_normalize_vec(emom(:,jj,ii),3)
                  emomM(:,jj,ii)=emom(:,jj,ii)*mmom(jj,ii)
               enddo
            enddo
         else if (do_mom_legacy.eq.'Y') then
            !------------------------------------------------------------------------
            ! If in the legacy format read the end time of the previous simulation
            !------------------------------------------------------------------------
            read(ifileno,*) rstep
            !------------------------------------------------------------------------
            ! Read the magnetic configuration
            !------------------------------------------------------------------------
            do ii=1, Mensemble
               do jj=1, Natom
                  read(ifileno,*) dummy, dummy, mmom(jj,ii), emom(1,jj,ii),         &
                     emom(2,jj,ii), emom(3,jj,ii)
                  emomM(:,jj,ii)=emom(:,jj,ii)*mmom(jj,ii)
               enddo
            enddo
         endif
         close(ifileno)
      else
         write(*,*) 'ERROR: Restartfile ',trim(adjustl(restartfile)), ' does not exist.'
         stop
      end if
   end subroutine read_mag_conf_std

   !---------------------------------------------------------------------------------
   !> @brief Routine handling the read-in of magnetic configurations with a single extra 
   !> list
   !---------------------------------------------------------------------------------
   subroutine read_mag_conf_single(Natom,Mensemble,do_mom_legacy,rstep,restartfile, &
      mmom,emom,emomM,aux_list,aux_sum)

      use math_functions, only: f_normalize_vec
      
      implicit none

      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=1), intent(in) :: do_mom_legacy  !< Flag to print/read moments in legacy output
      integer, intent(out) :: rstep !< Starting simulation step
      integer, intent(inout) :: aux_sum !< Sumation of the entries of the auxiliary array
      character(len=35), intent(inout) :: restartfile !< File containing restart information
      integer, dimension(Natom), intent(inout) :: aux_list  !< Auxiliary array
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom  !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! .. Local variables
      integer :: dummy, ii, jj,ios
      logical :: exists

      inquire(file=restartfile,exist=exists)
      if(exists) then
         open(ifileno,iostat=ios, file=restartfile, status="old")
         if (do_mom_legacy.ne.'Y') then
            !------------------------------------------------------------------------
            ! If in the new format read the first 7 lines to skip the header
            !------------------------------------------------------------------------
            do ii=1,7
               read(ifileno,*)
            enddo
            !------------------------------------------------------------------------
            ! Read the magnetic configurations
            !------------------------------------------------------------------------
            do ii=1, Mensemble
               do jj=1, Natom
                  read(ifileno,*) rstep, dummy, dummy, mmom(jj,ii), emom(1,jj,ii),  &
                     emom(2,jj,ii), emom(3,jj,ii), aux_list(jj)
                  emom(:,jj,ii)=f_normalize_vec(emom(:,jj,ii),3)
                  emomM(:,jj,ii)=emom(:,jj,ii)*mmom(jj,ii)
               enddo
            enddo
         else if (do_mom_legacy.eq.'Y') then
            !------------------------------------------------------------------------
            ! If in the legacy format read the end time of the previous simulation
            !------------------------------------------------------------------------
            read(ifileno,*) rstep
            !------------------------------------------------------------------------
            ! Read the magnetic configuration
            !------------------------------------------------------------------------
            do ii=1, Mensemble
               do jj=1, Natom
                  read(ifileno,*) dummy, dummy, mmom(jj,ii), emom(1,jj,ii),         &
                     emom(2,jj,ii), emom(3,jj,ii), aux_list(jj) 
                  emomM(:,jj,ii)=emom(:,jj,ii)*mmom(jj,ii)
               enddo
            enddo
         endif
         close(ifileno)
      else
         write(*,*) 'ERROR: Restartfile ',trim(adjustl(restartfile)), ' does not exist.'
         stop
      end if

      ! Summing over the entries of the array
      aux_sum=sum(aux_list(:))

   end subroutine read_mag_conf_single

   !---------------------------------------------------------------------------------
   !> @brief Routine handling the read-in of magnetic configurations with a single extra 
   !> list
   !---------------------------------------------------------------------------------
   subroutine read_mag_conf_double(Natom,Mensemble,do_mom_legacy,rstep,restartfile, &
      mmom,emom,emomM,aux_list_1,aux_list_2,aux_sum)

      use math_functions, only: f_normalize_vec
      
      implicit none

      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=1), intent(in) :: do_mom_legacy  !< Flag to print/read moments in legacy output
      integer, intent(out) :: rstep !< Starting simulation step
      integer, intent(inout) :: aux_sum !< Sumation of the entries of the auxiliary array
      character(len=35), intent(inout) :: restartfile !< File containing restart information
      integer, dimension(Natom), intent(inout) :: aux_list_1  !< First auxiliary array
      integer, dimension(Natom), intent(inout) :: aux_list_2  !< Second auxiliary array
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom  !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! .. Local variables
      integer :: dummy, ii, jj,ios
      logical :: exists

      inquire(file=restartfile,exist=exists)
      if(exists) then
         open(ifileno,iostat=ios, file=restartfile, status="old")
         if (do_mom_legacy.ne.'Y') then
            !------------------------------------------------------------------------
            ! If in the new format read the first 7 lines to skip the header
            !------------------------------------------------------------------------
            do ii=1,7
               read(ifileno,*)
            enddo
            !------------------------------------------------------------------------
            ! Read the magnetic configurations
            !------------------------------------------------------------------------
            do ii=1, Mensemble
               do jj=1, Natom
                  read(ifileno,*) rstep, dummy, dummy, mmom(jj,ii), emom(1,jj,ii),  &
                     emom(2,jj,ii), emom(3,jj,ii), aux_list_1(jj), aux_list_2(jj)
                  emom(:,jj,ii)=f_normalize_vec(emom(:,jj,ii),3)
                  emomM(:,jj,ii)=emom(:,jj,ii)*mmom(jj,ii)
               enddo
            enddo
         else if (do_mom_legacy.eq.'Y') then
            !------------------------------------------------------------------------
            ! If in the legacy format read the end time of the previous simulation
            !------------------------------------------------------------------------
            read(ifileno,*) rstep
            !------------------------------------------------------------------------
            ! Read the magnetic configuration
            !------------------------------------------------------------------------
            do ii=1, Mensemble
               do jj=1, Natom
                  read(ifileno,*) dummy, dummy, mmom(jj,ii), emom(1,jj,ii),         &
                     emom(2,jj,ii), emom(3,jj,ii), aux_list_1(jj), aux_list_2(jj)
                  emomM(:,jj,ii)=emom(:,jj,ii)*mmom(jj,ii)
               enddo
            enddo
         endif
         close(ifileno)
      else
         write(*,*) 'ERROR: Restartfile ',trim(adjustl(restartfile)), ' does not exist.'
         stop
      end if
      ! Summing over the entries of the array
      aux_sum=sum(aux_list_1(:))

   end subroutine read_mag_conf_double

   !---------------------------------------------------------------------------------
   !> @brief Wrapper for reading the GNEB configurations
   !---------------------------------------------------------------------------------
   subroutine GNEB_read_wrapper(Natom,Mensemble,amp_rnd,mode,relaxed_if,            &
      do_mom_legacy,rstep,exists,restartfile,mmom,emom,emomM)

      use math_functions, only: f_normalize_vec
      use RandomNumbers, only : rng_uniform

      implicit none

      ! .. Input variables
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: amp_rnd !< Amplitude of random perturbation of the components of magnetic moments
      character(len=4), intent(in) :: mode !< Type of path that is being read
      character(len=1), intent(in) :: relaxed_if !< Use relaxed ini and fin states
      character(len=1), intent(in) :: do_mom_legacy  !< Flag to print/read moments in legacy output
      ! .. Output variables
      integer, intent(out) :: rstep !< Starting simulation step
      logical, intent(out) :: exists !< See if the file actually exists
      ! .. In/Out variables
      character(len=35), intent(inout) :: restartfile !< File containing restart information
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom  !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! .. Local variables
      integer :: ii, jj
      integer :: i_stat,i_all
      real(dblprec) :: tol
      real(dblprec), dimension(3) :: u, rn
      real(dblprec), dimension(:,:), allocatable :: mmom_tmp
      real(dblprec), dimension(:,:,:), allocatable :: emom_tmp
      real(dblprec), dimension(:,:,:), allocatable :: emomM_tmp

      u  = 0.0_dblprec
      rn = 0.0_dblprec

      tol=1e-15_dblprec
      !------------------------------------------------------------------------
      ! Read the full GNEB path
      !------------------------------------------------------------------------
      if (mode=='path') then
         !------------------------------------------------------------------------
         ! The file is in the new format
         !------------------------------------------------------------------------
         if (do_mom_legacy.ne.'Y') then 
            !------------------------------------------------------------------------
            ! If the intial and final states are fixed
            !------------------------------------------------------------------------
            if (relaxed_if=='Y') then
               !---------------------------------------------------------------------
               ! Allocate auxiliary arrays
               !---------------------------------------------------------------------
               allocate(mmom_tmp(Natom,Mensemble),stat=i_stat)
               call memocc(i_stat,product(shape(mmom_tmp))*kind(mmom_tmp),'mmom_tmp','GNEB_read_wrapper')
               mmom_tmp=0.0_dblprec
               allocate(emom_tmp(3,Natom,Mensemble),stat=i_stat)
               call memocc(i_stat,product(shape(emom_tmp))*kind(emom_tmp),'emom_tmp','GNEB_read_wrapper')
               emom_tmp=0.0_dblprec
               allocate(emomM_tmp(3,Natom,Mensemble),stat=i_stat)
               call memocc(i_stat,product(shape(emomM_tmp))*kind(emomM_tmp),'emomM_tmp','GNEB_read_wrapper')
               emomM_tmp=0.0_dblprec
               !---------------------------------------------------------------------
               ! Read the configuration from file
               !---------------------------------------------------------------------
               call read_mag_conf(Natom,Mensemble,'N',rstep,restartfile,mmom_tmp,   &
                  emom_tmp,emomM_tmp)
               !---------------------------------------------------------------------
               ! Add noise to the configuration if it is necessary
               !---------------------------------------------------------------------
               if (amp_rnd>tol) then
                  !$omp parallel do default(shared), private(ii,jj,u,rn)
                  do jj=2,Mensemble-1
                     do ii=1, Natom
                        call rng_uniform(rn,3)
                        u(:)=2.0_dblprec*(rn(:)-0.50_dblprec)
                        emom(:,ii,jj)=emom_tmp(:,ii,jj)+amp_rnd*u(:)
                        emom(:,ii,jj)=f_normalize_vec(emom(:,ii,jj),3)
                        mmom(ii,jj)=mmom_tmp(ii,jj)
                        emomM(:,ii,jj)=emom(:,ii,jj)*mmom(ii,jj)
                     enddo
                  enddo
                  !$omp end parallel do
               !---------------------------------------------------------------------
               ! If there is no noise just save the data to the proper arrays
               !---------------------------------------------------------------------
               else
                  !$omp parallel do default(shared), private(ii,jj)
                  do jj=2,Mensemble-1
                     do ii=1, Natom
                        emom(:,ii,jj)=emom_tmp(:,ii,jj)
                        mmom(ii,jj)=mmom_tmp(ii,jj)
                        emomM(:,ii,jj)=emomM(:,ii,jj)
                     enddo
                  enddo
                  !$omp end parallel do
               endif
               !---------------------------------------------------------------------
               ! Deallocate auxiliary arrays
               !---------------------------------------------------------------------
               i_all=-product(shape(mmom_tmp))*kind(mmom_tmp)
               deallocate(mmom_tmp,stat=i_stat)
               call memocc(i_stat,i_all,'mmom_tmp','GNEB_read_wrapper')
               i_all=-product(shape(emom_tmp))*kind(emom_tmp)
               deallocate(emom_tmp,stat=i_stat)
               call memocc(i_stat,i_all,'emom_tmp','GNEB_read_wrapper')
               i_all=-product(shape(emomM_tmp))*kind(emomM_tmp)
               deallocate(emomM_tmp,stat=i_stat)
               call memocc(i_stat,i_all,'emom_Mtmp','GNEB_read_wrapper')
            !------------------------------------------------------------------------
            ! Read all the states
            !------------------------------------------------------------------------
            else
               !---------------------------------------------------------------------
               ! Read the configuration from file
               !---------------------------------------------------------------------
               call read_mag_conf(Natom,Mensemble,'N',rstep,restartfile,mmom,emom,  &
                  emomM)
               !---------------------------------------------------------------------
               ! Add noise to the configuration if it is necessary
               !---------------------------------------------------------------------
               if (amp_rnd>tol) then
                  !$omp parallel do default(shared), private(ii,jj,u,rn)
                  do jj=1,Mensemble
                     do ii=1, Natom
                        call rng_uniform(rn,3)
                        u(:)=2.0_dblprec*(rn(:)-0.50_dblprec)
                        emom(:,ii,jj)=emom(:,ii,jj)+amp_rnd*u(:)
                        emom(:,ii,jj)=f_normalize_vec(emom(:,ii,jj),3)
                        mmom(ii,jj)=mmom(ii,jj)
                        emomM(:,ii,jj)=emom(:,ii,jj)*mmom(ii,jj)
                     enddo
                  enddo
                  !$omp end parallel do
               endif
            endif
         !---------------------------------------------------------------------------
         ! The file is in the legacy format
         !---------------------------------------------------------------------------
         else if (do_mom_legacy.eq.'Y') then
            call load_path_legacy(Natom,Mensemble,restartfile,amp_rnd,relaxed_if,   &
               mmom,emom,emomM,exists) 
         endif
      !------------------------------------------------------------------------------
      ! Read the initial and final GNEB states
      !------------------------------------------------------------------------------
      else if (mode=='infi') then
         if (do_mom_legacy.ne.'Y') then 
            !------------------------------------------------------------------------
            ! Read the configuration from file
            !------------------------------------------------------------------------
            call read_mag_conf(Natom,2,'N',rstep,restartfile,                       &
               mmom(:,1:Mensemble:(Mensemble-1)),                                   &
               emom(:,:,1:Mensemble:(Mensemble-1)),                                 &
               emomM(:,:,1:Mensemble:(Mensemble-1)))
            !------------------------------------------------------------------------
            ! Add noise to the configuration if it is necessary
            !------------------------------------------------------------------------
            if (amp_rnd>tol) then
               !$omp parallel do default(shared), private(ii,jj,u)
               do jj=1,Mensemble,(Mensemble-1)
                  do ii=1,Natom
                     call rng_uniform(rn,3)
                     u(:)=2.0_dblprec*(rn(:)-0.50_dblprec)
                     emom(:,ii,jj)=emom(:,ii,jj)+amp_rnd*u(:)
                     emom(:,ii,jj)=f_normalize_vec(emom(:,ii,jj),3)
                     emomM(:,ii,jj)=emom(:,ii,jj)*mmom(ii,jj)
                  enddo
               enddo
               !$omp end parallel do
            endif
            !$omp parallel do default(shared), private(ii,jj), collapse(2)
            do jj=2,Mensemble-1
               do ii=1, Natom
                  mmom(ii,jj) = mmom(ii,1)
               end do
            end do
            !$omp end parallel do
         else if (do_mom_legacy.eq.'Y') then
            call load_inifin_legacy(Natom,Mensemble,restartfile,amp_rnd,mmom,emom,  &
               emomM)
         endif

      endif

   end subroutine GNEB_read_wrapper
   !---------------------------------------------------------------------------------
   ! LEGACY DEFINITIONS
   !---------------------------------------------------------------------------------
   !> Prints magnetic moment to restart file
   subroutine prnrestart(Natom, Mensemble, simid, mstep, emom, mmom)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      character(len=8), intent(in) :: simid !< Name of simulation 
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments

      integer :: ii, jj
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''restart.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)
      write (ofileno,10002) mstep

      do ii=1, Mensemble
         do jj=1, Natom
            write (ofileno,10003) ii, jj, mmom(jj,ii), emom(1,jj,ii), emom(2,jj,ii), emom(3,jj,ii)
         end do
      end do

      close(ofileno)

      10002 format (i8)
      10003 format (i8,i8,2x,es16.8,es16.8,es16.8,es16.8)

   end subroutine prnrestart

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_tottraj_legacy
   !> @brief Prints the total trajectories in the legacy format
   !> @details Prints a magnetic configuration, in the legacy format
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine prn_tottraj_legacy(Natom,Mensemble,tottraj_buff,bcount_tottraj,simid, &
      real_time_measure,indxb,emomb)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: tottraj_buff
      integer, intent(in) :: bcount_tottraj
      character(len=8), intent(in) :: simid !< Name of simulation 
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      real(dblprec), dimension(tottraj_buff), intent(in) :: indxb
      real(dblprec), dimension(3,Natom,tottraj_buff,Mensemble), intent(in) :: emomb   !< Current unit moment vector

      integer :: ii,jj
      character(len=30) :: filn

      write (filn,'(''moment.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")
      do ii=1, bcount_tottraj
         do jj=1, Natom
            if (real_time_measure=='Y') then
               write (ofileno,10003) indxb(ii), jj, emomb(1,jj,ii,Mensemble),       &
                  emomb(2,jj,ii,Mensemble), emomb(3,jj,ii,Mensemble),               &
                  emomb(1,jj,ii,1)**2+emomb(2,jj,ii,1)**2+emomb(3,jj,ii,1)**2
            else
               write (ofileno,10002) int(indxb(ii)), jj, emomb(1,jj,ii,Mensemble),  &
                  emomb(2,jj,ii,Mensemble), emomb(3,jj,ii,Mensemble),               &
                  emomb(1,jj,ii,1)**2+emomb(2,jj,ii,1)**2+emomb(3,jj,ii,1)**2
            endif
         end do
      end do

      close(ofileno)
      return
      10002 format (i8,2x,i8,2x,2x, es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)
      10003 format (es16.4,2x,i8,2x,2x, es16.8,es16.8,es16.8,es16.8)

   end subroutine prn_tottraj_legacy
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: load_path_legacy
   !> @brief Read path from file for a GNEB calculation.
   !> @details Ensembles correspond images of the GNEB method, in this approach all
   !> the images must be provided.
   !> @author Pavel Bessarab
   !----------------------------------------------------------------------------
   subroutine load_path_legacy(Natom,Mensemble,fname,amp_rnd,relaxed_if,mmom,emom,  &
      emomM,exists)
      !
      use RandomNumbers, only : rng_uniform
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=35), intent(in) :: fname !< File containing the path
      real(dblprec), intent(in) :: amp_rnd !< Amplitude of random perturbation of the components of magnetic moments
      character(len=1), intent(in) :: relaxed_if !< Use relaxed ini and fin states
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      logical, intent(out) :: exists
      integer :: ii, jj, dummy, ios
      real(dblprec) :: mmom_tmp
      real(dblprec), dimension(3) :: u,v,rn

      inquire(file=fname,exist=exists)
      if(exists) then
         open(12, iostat=ios, file=fname, status="old")
         do ii=1,Mensemble
            do jj=1, Natom
               if ((ii==1).or.(ii==Mensemble)) then
                  u(:) = 0.0_dblprec
               else
                  call rng_uniform(rn,3)
                  u(:)=2.0_dblprec*(rn(:)-0.50_dblprec)
               end if
               if (((ii==1).or.(ii==Mensemble)).and.(relaxed_if=='Y')) then
                  read (12,*) dummy, dummy, v(1), v(2), v(3), mmom_tmp
               else
                  read (12,*) dummy, dummy, emom(1,jj,ii), emom(2,jj,ii), emom(3,jj,ii), mmom(jj,ii)
                  v(1:3) = emom(1:3,jj,ii)+amp_rnd*u(1:3)
                  mmom_tmp = norm2(v)
                  emom(1:3,jj,ii) = v(1:3)/mmom_tmp
                  emomM(:,jj,ii)=emom(:,jj,ii)*mmom(jj,ii)
               end if
            end do
         end do
         close(12)
      else
         write(*,*) 'ERROR: File ',trim(adjustl(fname)), ' does not exist. Path not loaded.'
      end if

   end subroutine load_path_legacy
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: load_inifin_legacy
   !> @brief Read initial and final state from file for GNEB.
   !> @details Read of the initial and final images for a GNEB calculation.
   !> The first ensemble corresponds to the initial state, the last ensemble correspond to the final state
   !> @author Pavel Bessarab
   !----------------------------------------------------------------------------
   subroutine load_inifin_legacy(Natom,Mensemble,fname,amp_rnd,mmom,emom,emomM)
      !
      use RandomNumbers, only : rng_uniform
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=35), intent(in) :: fname !< File containing info on magnetic configuration
      real(dblprec), intent(in) :: amp_rnd !< Amplitude of random perturbation of the components of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(out) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector

      integer :: ii, jj, dummy, ios,NS
      logical :: exists
      real(dblprec) :: mmom_tmp
      real(dblprec), dimension(3) :: u,v,rn

      inquire(file=fname,exist=exists)
      if(exists) then
         NS = number_of_strings(fname)
         print *, "number of strings:",NS
         !read first image
         open(ifileno, iostat=ios, file=fname, status="old")
         do jj=1, Natom
               call rng_uniform(rn,3)
               u(:)=2.0_dblprec*(rn(:)-0.50_dblprec)
               read (ifileno,*) dummy, dummy, emom(1,jj,1), emom(2,jj,1), emom(3,jj,1), mmom(jj,1)
               v(1:3) = emom(1:3,jj,1)+amp_rnd*u(1:3)
               mmom_tmp = norm2(v)
               emom(1:3,jj,1) = v(1:3)/mmom_tmp
               emomM(:,jj,1)=emom(:,jj,1)*mmom(jj,1)
         end do
         close(ifileno)
         !read final image
         open(ifileno, iostat=ios, file=fname, status="old")
         do jj=1,NS-Natom
            read (ifileno,*) dummy, dummy, emom(1,1,Mensemble), emom(2,1,Mensemble),emom(3,1,Mensemble), mmom(1,Mensemble)
         end do
         
         do jj=1, Natom
               call rng_uniform(rn,3)
               u(:)=2.0_dblprec*(rn(:)-0.50_dblprec)
               read (ifileno,*) dummy, dummy, emom(1,jj,Mensemble), emom(2,jj,Mensemble),emom(3,jj,Mensemble), mmom(jj,Mensemble)
               v(1:3) = emom(1:3,jj,Mensemble)+amp_rnd*u(1:3)
               mmom_tmp = norm2(v)
               emom(1:3,jj,Mensemble) = v(1:3)/mmom_tmp
               emomM(:,jj,Mensemble)=emom(:,jj,Mensemble)*mmom(jj,Mensemble)
         end do
         do ii=2,Mensemble-1
            do jj=1, Natom
               mmom(jj,ii) = mmom(jj,1)
            end do
         end do
         close(ifileno)
      else
         write(*,*) 'ERROR: File ',trim(adjustl(fname)), ' does not exist.'
         stop
      end if
      CONTAINS
      function number_of_strings(filename) result(N)
         implicit none
         CHARACTER(len=*) :: filename
         INTEGER:: N
         CHARACTER(len=200)   :: STRING
         INTEGER :: ios !Input/Output Status

         open(953,file=filename,iostat=ios, action='read', status='old')
         
         N=0
         do while (ios/=-1)
            READ(unit=953,fmt='(A)', iostat=ios) STRING
            if (ios==-1) cycle
            if (len(trim(STRING)) == 0) cycle
            N=N+1
         enddo
         close(953)
      end function number_of_strings
   end subroutine load_inifin_legacy

end module Restart
