!------------------------------------------------------------------------------------
!  MODULE: InputHandler
!> @brief
!> Reads each keyword in inpsd.dat.
!> @details If the first read character of a line is a %, #, * or = the line is considered
!> to be a comment and and is ignored.
!> Comments after keyword values does not need a %, #, * or =
!
!> @author
!> Anders Bergman
!> @copyright
!> GNU Public License.
!
!> @todo
!> Put all reads in separate modules
!------------------------------------------------------------------------------------
module InputHandler

   use QHB,                only : do_qhb, qhb_mode, tcurie
   use KMCData
   use clusters
   use FixedMom,           only : inp_fixed_mom_flag, do_fixed_mom
   use stiffness
   use Profiling
   use InputData
   use prn_fields
   use Parameters
   use macrocells,         only : do_macro_cells,prn_dip_subset,dip_file
   use temperature,        only : grad, tempfile
   use ChemicalData
   use prn_averages
   use Polarization
   use prn_topology
   use prn_currents
   use ErrorHandling
   use RandomNumbers
   use prn_induced_info,   only : do_prn_induced, ind_step,ind_buff
   use prn_trajectories

   implicit none

   logical :: sane_input = .true.

   private

   public :: sane_input
   public :: allocate_hamiltonianinput, read_exchange_getMaxNoShells, read_parameters
   public :: read_positions, read_positions_alloy, change_constants, read_moments, read_exchange
   public :: read_exchange_tensor, read_exchange_build_tensor, read_anisotropy_alloy, read_anisotropy
   public :: read_dmdata, read_pddata, read_chirdata, read_biqdmdata, read_bqdata, read_sitefield
   public :: read_ip_damping, read_ip_damping_alloy, read_damping, read_damping_alloy
   public :: read_barriers, read_fixed_moments, read_exchange_getNeighVec

contains

   !--------------------------------------------------------------------------------
   !> @brief
   !> Read input parameters
   !
   !> @author
   !> Anders Bergman
   !>
   !> @date 08/02/2017 - Jonathan Chico
   !> - Reorganized input variables in blocks by where they are used.
   !> IMPORTANT TRY TO KEEP ORDER IN THE VARIABLES
   !--------------------------------------------------------------------------------
   subroutine read_parameters(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword,cache, string
      integer :: rd_len,i_err,i,i_stat,i_errb,ii, i_all
      logical :: comment
      real(dblprec) :: tmp

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

            !------------------------------------------------------------------------
            ! START OF MISC VARIABLES
            !------------------------------------------------------------------------

            case('simid')
               read(ifile,*,iostat=i_err) simid
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('aunits')
               read(ifile,*,iostat=i_err) aunits
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('set_landeg')
               read(ifile,*,iostat=i_err) set_landeg
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mensemble')
               read(ifile,*,iostat=i_err) mensemble
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF MISC VARIABLES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR CRYSTAL STRUCTURE
            !------------------------------------------------------------------------

            case('ncell')
               read(ifile,*,iostat=i_err) N1, N2, N3
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('cell')
               read(ifile,*,iostat=i_err) C1, C2, C3
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('bc')
               read(ifile,*,iostat=i_err) BC1, BC2, BC3
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('natoms')
               read(ifile,*,iostat=i_err) na
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ntypes')
               read(ifile,*,iostat=i_err) nt
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('positions')
               if(na==0) then
                  write(*,*) 'ERROR: ','natoms not set before reading ',keyword
               else
                  allocate(atype_inp(na),stat=i_stat)
                  call memocc(i_stat,product(shape(atype_inp))*kind(atype_inp),'atype_inp','read_parameters')
                  allocate(anumb_inp(na),stat=i_stat)
                  call memocc(i_stat,product(shape(anumb_inp))*kind(anumb_inp),'anumb_inp','read_parameters')
                  allocate(bas(3,na),stat=i_stat)
                  call memocc(i_stat,product(shape(bas))*kind(bas),'bas','read_parameters')
                  do i=1,na
                     read(ifile,*,iostat=i_err) anumb_inp(i), atype_inp(i), bas(1:3,i)
                     if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
                  end do
                  nt=maxval(atype_inp)
               end if

            case('posfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               posfile=adjustl(trim(cache))

            case('posfiletype')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               posfiletype=adjustl(trim(cache))

            case('alat')
               read(ifile,*,iostat=i_err) alat
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('scalefac')
               read(ifile,*,iostat=i_err) scalefac
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_ralloy')
               read(ifile,*,iostat=i_err) do_ralloy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('block_size')
               read(ifile,*,iostat=i_err) block_size
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR CRYSTAL STRUCTURE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR MOMENTS
            !------------------------------------------------------------------------

            case('momfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               momfile=adjustl(trim(cache))

            case('initmag')
               read(ifile,*,iostat=i_err) initmag
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('restartfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               restartfile=trim(adjustl(cache))

            case('mseed')
               read(ifile,*,iostat=i_err) mseed
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('theta0')
               read(ifile,*,iostat=i_err) theta0
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('phi0')
               read(ifile,*,iostat=i_err) phi0
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('roteul')
               read(ifile,*,iostat=i_err) roteul
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('rotang')
               read(ifile,*,iostat=i_err) rotang
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('initpropvec')
               read(ifile,*,iostat=i_err) initpropvec
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('initrotvec')
               read(ifile,*,iostat=i_err) initrotvec
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('initrotang')
               read(ifile,*,iostat=i_err) initrotang
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_fixed_mom')
               read(ifile,*,iostat=i_err) do_fixed_mom
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_mom_legacy')
               read(ifile,*,iostat=i_err) do_mom_legacy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR MOMENTS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR EXCHANGE
            !------------------------------------------------------------------------

            !> - exchange
            !! Name of exchange file
            case('exchange')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               allocate(jfile(conf_num),stat=i_stat)
               call memocc(i_stat,product(shape(jfile))*kind(jfile),'jfile','inputhandler')
               if (conf_num==1) then
                  jfile=adjustl(trim(cache))
                  ! For the LSF method there are several jfiles that must be read, depending on the number of configurations
               else
                  pre_jfile=adjustl(trim(cache))
                  i=len_trim(pre_jfile)-1
                  do ii=1,conf_num
                     if (ii<10) then
                        write(string,'("(a",i0,",i1)")') i
                        write(jfile(ii),string) adjustl(trim(pre_jfile)),ii
                     else if (ii>9.and.ii<100) then
                        write(string,'("(a",i0,",i2)")') i
                        write(jfile(ii),string) adjustl(trim(pre_jfile)),ii
                     else if (ii>99.and.ii<1000) then
                        write(string,'("(a",i0,",i3)")') i
                        write(jfile(ii),string) adjustl(trim(pre_jfile)),ii
                     else if (ii>999.and.ii<10000) then
                        write(string,'("(a",i0,",i4)")') i
                        write(jfile(ii),string) adjustl(trim(pre_jfile)),ii
                     endif
                  enddo
               endif

            case('exchangedlm')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               allocate(jfileD(conf_num),stat=i_stat)
               call memocc(i_stat,product(shape(jfileD))*kind(jfileD),'jfileD','inputhandler')
               if (conf_num==1) then
                  jfileD=adjustl(trim(cache))
               else
                  pre_jfile=adjustl(trim(cache))
                  i=len_trim(pre_jfile)-1
                  do ii=1,conf_num
                     if (ii<10) then
                        write(string,'("(a",i0,",i1)")') i
                        write(jfileD(ii),string) adjustl(trim(pre_jfile)),ii
                     else if (ii>9.and.ii<100) then
                        write(string,'("(a",i0,",i2)")') i
                        write(jfileD(ii),string) adjustl(trim(pre_jfile)),ii
                     else if (ii>99.and.ii<1000) then
                        write(string,'("(a",i0,",i3)")') i
                        write(jfileD(ii),string) adjustl(trim(pre_jfile)),ii
                     else if (ii>999.and.ii<10000) then
                        write(string,'("(a",i0,",i4)")') i
                        write(jfileD(ii),string) adjustl(trim(pre_jfile)),ii
                     endif
                  enddo
               endif

            case('jij_scale')
               read(ifile,*,iostat=i_err) jij_scale
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('exc_inter')
               read(ifile,*,iostat=i_err) exc_inter
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ! Type of maptype 1 for cartesian bonding vectors, 2 for units of repetitions of the lattice vectors
            case('maptype')
               read(ifile,*,iostat=i_err) maptype
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('map_multiple')
               read(ifile,*,iostat=i_err) map_multiple
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ! Symmetry of the exchange interactions, turned off by default when DMI is considered
            case('sym')
               read(ifile,*,iostat=i_err) sym
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_sortcoup')
               read(ifile,*,iostat=i_err) do_sortcoup
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_jtensor')
               read(ifile,*,iostat=i_err) do_jtensor
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('calc_jtensor')
               read(ifile,*,iostat=i_err) calc_jtensor
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR EXCHANGE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR ANISOTROPIES
            !------------------------------------------------------------------------

            case('do_anisotropy')
               read(ifile,*,iostat=i_err) do_anisotropy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('anisotropy')
               do_anisotropy=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               kfile=trim(adjustl(cache))

            case('mult_axis')
               read(ifile,*,iostat=i_err) mult_axis
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('random_anisotropy')
               do_anisotropy=1
               read(ifile,*,iostat=i_err) random_anisotropy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('random_anisotropy_density')
               read(ifile,*,iostat=i_err) random_anisotropy_density
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR ANISOTROPIES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR DMI
            !------------------------------------------------------------------------

            case('dm')
               do_dm=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               dmfile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(dmfile, &
                  'Please specify dm <dmfile> where <dmfile> is a valid dm interaction file')

            case('dm_scale')
               read(ifile,*,iostat=i_err) dm_scale
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err


            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR DMI
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR DIPOLAR
            !------------------------------------------------------------------------

            case('do_dip')
               read(ifile,*,iostat=i_err) do_dip
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('print_dip_tensor')
               read(ifile,*,iostat=i_err) print_dip_tensor
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qdip_files')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               qdip_files=trim(adjustl(cache))

            case('read_dipole')
               read(ifile,'(a)',iostat=i_err) read_dipole
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('prn_dip_subset')
               read(ifile,'(a)',iostat=i_err) prn_dip_subset
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('dip_file')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               dip_file=trim(adjustl(cache))
            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR DIPOLAR
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            !  START OF VARIABLES FOR PSEUDO DIPOLAR
            !------------------------------------------------------------------------

            case('pd')
               do_pd=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               pdfile=trim(adjustl(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PSEUDO DIPOLAR
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR Scalar Chiral
            !------------------------------------------------------------------------

            case('chir')
               do_chir=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               chirfile=trim(adjustl(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR Scalar Chiral
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR BIQUADRATIC
            !------------------------------------------------------------------------

            case('bq')
               do_bq=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               bqfile=trim(adjustl(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR BIQUADRATIC
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR BIQUADRATIC DMI
            !------------------------------------------------------------------------

            case('biqdm')
               do_biqdm=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               biqdmfile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(biqdmfile, &
                  'Please specify biqdm <biqdmfile> where <biqdmfile> is a valid biqdm interaction file')

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR BIQUADRATIC DMI
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR RANDOM NUMBERS
            !------------------------------------------------------------------------

            case('ziggurat')
               read(ifile,*,iostat=i_err) ziggurat
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('rngpol')
               read(ifile,*,iostat=i_err) rngpol
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('tseed')
               read(ifile,*,iostat=i_err) tseed
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('para_rng')
               read(ifile,*,iostat=i_err) para_rng
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('use_vsl')
               read(ifile,*,iostat=i_err) use_vsl
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR RANDOM NUMBERS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR TEMPERATURES
            !------------------------------------------------------------------------

            case('temperature')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               tempfile=adjustl(trim(cache))

            case('gradient')
               read(ifile,*,iostat=i_err) grad
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('temp')
               read(ifile,*,iostat=i_err) temp
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR TEMPERATURES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR DAMPING
            !------------------------------------------------------------------------

            case('do_site_ip_damping')
               read(ifile,*,iostat=i_err) do_site_ip_damping
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ip_damping2')
               allocate(iplambda2(1),stat=i_stat)
               call memocc(i_stat,product(shape(iplambda2))*kind(iplambda2),'iplambda2','read_parameters')
               read(ifile,*,iostat=i_err) iplambda2(1)
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_site_damping')
               read(ifile,*,iostat=i_err) do_site_damping
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('damping')
               read(ifile,*,iostat=i_err) mplambda1
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('damping2')
               read(ifile,*,iostat=i_err) mplambda2
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ip_dampfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ip_dampfile=adjustl(trim(cache))

            case('mp_dampfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               mp_dampfile=adjustl(trim(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR DAMPING
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR SIMULATION MODE
            !------------------------------------------------------------------------

            case('ip_mode')
               read(ifile,*,iostat=i_err) ipmode
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mode')
               read(ifile,*,iostat=i_err) mode
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR SIMULATION MODE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES COMMON FOR IP_MODE
            !------------------------------------------------------------------------

            case('ip_temp')
               allocate(ipTemp(1),stat=i_stat)
               call memocc(i_stat,product(shape(ipTemp))*kind(ipTemp),'ipTemp','read_parameters')
               read(ifile,*,iostat=i_err) ipTemp(1)
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ip_hfield')
               read(ifile,*,iostat=i_err) iphfield
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR COMMON FOR IP_MODE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES COMMON FOR MODE
            !------------------------------------------------------------------------

            case('hfield')
               read(ifile,*,iostat=i_err) hfield
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR COMMON FOR MODE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR ASD MODE
            !------------------------------------------------------------------------

            case('nstep')
               read(ifile,*,iostat=i_err) nstep
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(mcnstep==0.and.(mode=='M'.or.mode=='H')) then
                  mcnstep=nstep
                  sane_input=.false.
               end if

            case('timestep')
               read(ifile,*,iostat=i_err) delta_t
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('relaxtime')
               read(ifile,*,iostat=i_err) relaxtime
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sdealgh')
               read(ifile,*,iostat=i_err) sdealgh
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ipsdealgh')
               read(ifile,*,iostat=i_err) ipsdealgh
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('real_time_measure')
               read(ifile,*,iostat=i_err) real_time_measure
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('perp')
               read(ifile,*,iostat=i_err) perp
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ip_nphase')
               read(ifile,*,iostat=i_err) ipnphase
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(ipnphase>0) then
                  if (allocated(ipTemp)) then
                     i_all=-product(shape(ipTemp))*kind(ipTemp)
                     deallocate(ipTemp,stat=i_stat)
                     call memocc(i_stat,i_all,'ipTemp','read_parameters')
                  endif
                  tmp=0.0_dblprec
                  if(allocated(iplambda2)) then
                     tmp=iplambda2(1)
                     i_all=-product(shape(iplambda2))*kind(iplambda2)
                     deallocate(iplambda2,stat=i_stat)
                     call memocc(i_stat,i_all,'iplambda2','read_parameters')
                  endif
                  allocate(ipTemp(ipnphase),stat=i_stat)
                  call memocc(i_stat,product(shape(ipTemp))*kind(iptemp),'ipTemp','read_parameters')
                  allocate(ipnstep(ipnphase),stat=i_stat)
                  call memocc(i_stat,product(shape(ipnstep))*kind(ipnstep),'ipnstep','read_parameters')
                  allocate(ipdelta_t(ipnphase),stat=i_stat)
                  call memocc(i_stat,product(shape(ipdelta_t))*kind(ipdelta_t),'ipdelta_t','read_parameters')
                  allocate(iplambda1(ipnphase),stat=i_stat)
                  call memocc(i_stat,product(shape(iplambda1))*kind(iplambda1),'iplambda1','read_parameters')
                  iplambda1=0.0_dblprec
                  allocate(iplambda2(ipnphase),stat=i_stat)
                  call memocc(i_stat,product(shape(iplambda2))*kind(iplambda2),'iplambda2','read_parameters')
                  iplambda2=tmp
                  if(ipmode=='R'.or.ipmode=='P') then
                     call ErrorHandling_missing('Spin-lattice dynamics')
                     !!! allocate(iplattdamp(ipnphase),stat=i_stat)
                     !!! call memocc(i_stat,product(shape(iplattdamp))*kind(iplattdamp),'iplattdamp','read_parameters')
                     !!! do i=1,ipnphase
                     !!!    read(ifile,*,iostat=i_err) ipnstep(i),ipTemp(i),ipdelta_t(i),iplambda1(i),iplattdamp(i)
                     !!!    !read(ifile,*,iostat=i_err) ipnstep(i),ipTemp(i),ipdelta_t(i),iplambda1(i)
                     !!!    if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
                     !!! end do
                  else if (ipmode=='M'.or.ipmode=='H') then
                     ipmcnphase=ipnphase
                     if(ipmcnphase>0) then
                        if (allocated(ipTemp)) then
                           i_all=-product(shape(ipTemp))*kind(ipTemp)
                           deallocate(ipTemp,stat=i_stat)
                           call memocc(i_stat,i_all,'ipTemp','read_parameters')
                        endif
                        if (allocated(ipmcnstep)) then
                           i_all=-product(shape(ipmcnstep))*kind(ipmcnstep)
                           deallocate(ipmcnstep,stat=i_stat)
                           call memocc(i_stat,i_all,'ipmcnstep','read_parameters')
                        endif
                        allocate(ipTemp(ipmcnphase),stat=i_stat)
                        call memocc(i_stat,product(shape(ipTemp))*kind(iptemp),'ipTemp','read_parameters')
                        allocate(ipmcnstep(ipmcnphase),stat=i_stat)
                        call memocc(i_stat,product(shape(ipmcnstep))*kind(ipmcnstep),'ipmcnstep','read_parameters')
                        do i=1,ipmcnphase
                           read(ifile,*,iostat=i_err) ipmcnstep(i),ipTemp(i)
                           if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
                        end do
                     end if
                     sane_input=.false.
                  else
                     do i=1,ipnphase
                        read(ifile,*,iostat=i_err) ipnstep(i),ipTemp(i),ipdelta_t(i),iplambda1(i)
                        if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
                     end do
                  end if
               else
                  read(ifile,*)
               end if

            case('compensate_drift')
               read(ifile,*,iostat=i_err) compensate_drift
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('llg')
               read(ifile,*,iostat=i_err) llg
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR ASD MODE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR MC MODE
            !------------------------------------------------------------------------

            case('mcnstep')
               read(ifile,*,iostat=i_err) mcnstep
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(mode=='S'.and.nstep==1) then
                  nstep=mcnstep
                  sane_input=.false.
               end if

            case('ip_mcnstep')
               allocate(ipmcnstep(1),stat=i_stat)
               call memocc(i_stat,product(shape(ipmcnstep))*kind(ipmcnstep),'ipmcnstep','read_parameters')
               read(ifile,*,iostat=i_err) ipmcnstep(1)
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ip_mcanneal')
               read(ifile,*,iostat=i_err) ipmcnphase
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(ipmcnphase>0) then
                  if (allocated(ipTemp)) then
                     i_all=-product(shape(ipTemp))*kind(ipTemp)
                     deallocate(ipTemp,stat=i_stat)
                     call memocc(i_stat,i_all,'ipTemp','read_parameters')
                  endif
                  if (allocated(ipmcnstep)) then
                     i_all=-product(shape(ipmcnstep))*kind(ipmcnstep)
                     deallocate(ipmcnstep,stat=i_stat)
                     call memocc(i_stat,i_all,'ipmcnstep','read_parameters')
                  endif
                  allocate(ipTemp(ipmcnphase),stat=i_stat)
                  call memocc(i_stat,product(shape(ipTemp))*kind(iptemp),'ipTemp','read_parameters')
                  allocate(ipmcnstep(ipmcnphase),stat=i_stat)
                  call memocc(i_stat,product(shape(ipmcnstep))*kind(ipmcnstep),'ipmcnstep','read_parameters')
                  do i=1,ipmcnphase
                     read(ifile,*,iostat=i_err) ipmcnstep(i),ipTemp(i)
                     if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
                  end do
               else
                  read(ifile,*)
               endif

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR MC MODE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR GNEB MODE
            !------------------------------------------------------------------------

            case('momfile_i')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               momfile_i=trim(adjustl(cache))

            case('momfile_f')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               momfile_f=trim(adjustl(cache))

            case('spring')
               read(ifile,*,iostat=i_err) spring
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('initpath')
               read(ifile,*,iostat=i_err) initpath
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('amp_rnd')
               read(ifile,*,iostat=i_err) amp_rnd
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('amp_rnd_path')
               read(ifile,*,iostat=i_err) amp_rnd_path
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('min_algo')
               read(ifile,*,iostat=i_err) minalgo
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('min_itrmax')
               read(ifile,*,iostat=i_err) minitrmax
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mintraj_step')
               read(ifile,*,iostat=i_err) mintraj_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('min_ftol')
               read(ifile,*,iostat=i_err) minftol
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mep_itrmax')
               read(ifile,*,iostat=i_err) mepitrmax
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('meptraj_step')
               read(ifile,*,iostat=i_err) meptraj_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mep_ftol')
               read(ifile,*,iostat=i_err) mepftol
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mep_ftol_ci')
               read(ifile,*,iostat=i_err) mepftol_ci
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_gneb')
               read(ifile,'(a)',iostat=i_err) do_gneb
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_gneb_ci')
               read(ifile,'(a)',iostat=i_err) do_gneb_ci
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_norm_rx')
               read(ifile,'(a)',iostat=i_err) do_norm_rx
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('vpo_dt')
               read(ifile,*,iostat=i_err) vpodt
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('vpo_mass')
               read(ifile,*,iostat=i_err) vpomass
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sample_num')
               read(ifile,*,iostat=i_err) sample_num
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('en_zero')
               read(ifile,'(a)',iostat=i_err) en_zero
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_hess_ini')
               read(ifile,'(a)',iostat=i_err) do_hess_ini
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_hess_fin')
               read(ifile,'(a)',iostat=i_err) do_hess_fin
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_hess_sp')
               read(ifile,'(a)',iostat=i_err) do_hess_sp
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('is_afm')
               read(ifile,'(a)',iostat=i_err) is_afm
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('eig_zero')
               read(ifile,*,iostat=i_err) eig_0
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('relaxed_if')
               read(ifile,'(a)',iostat=i_err) relaxed_if
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('fixed_if')
               read(ifile,'(a)',iostat=i_err) fixed_if
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('prn_gneb_fields')
               read(ifile,'(a)',iostat=i_err) prn_gneb_fields
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR GNEB MODE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING MOMENTS TRAJECTORIES
            !------------------------------------------------------------------------

            case('do_tottraj')
               read(ifile,*,iostat=i_err) do_tottraj
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('tottraj_step')
               read(ifile,*,iostat=i_err) tottraj_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('tottraj_buff')
               read(ifile,*,iostat=i_err) tottraj_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ntraj')
               read(ifile,*,iostat=i_err) ntraj
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(ntraj>0) then
                  allocate(traj_atom(ntraj),stat=i_stat)
                  call memocc(i_stat,product(shape(traj_atom))*kind(traj_atom),'traj_atom','read_parameters')
                  allocate(traj_step(ntraj),stat=i_stat)
                  call memocc(i_stat,product(shape(traj_step))*kind(traj_step),'traj_step','read_parameters')
                  allocate(traj_buff(ntraj),stat=i_stat)
                  call memocc(i_stat,product(shape(traj_buff))*kind(traj_buff),'traj_buff','read_parameters')
                  do i=1,ntraj
                     read(ifile,*,iostat=i_err) traj_atom(i), traj_step(i), traj_buff(i)
                  end do
               else
                  read(ifile,*)
               end if

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING MOMENTS TRAJECTORIES
            !------------------------------------------------------------------------


            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR SPIN CURRENTS
            !------------------------------------------------------------------------

            case('do_currents')
               read(ifile,*,iostat=i_err) do_currents
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('current_step')
               read(ifile,*,iostat=i_err) current_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('current_buff')
               read(ifile,*,iostat=i_err) current_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('quantization_axis')
               read(ifile,*,iostat=i_err) quant_axis
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR SPIN CURRENTS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR SPIN TEMPERATURE
            !------------------------------------------------------------------------

            case('do_spintemp')
               read(ifile,*,iostat=i_err) do_spintemp
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('spintemp_step')
               read(ifile,*,iostat=i_err) spintemp_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR SPIN TEMPERATURE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR DEBUGGING
            !------------------------------------------------------------------------

            case('evolveout')
               read(ifile,*,iostat=i_err) evolveout
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('heisout')
               read(ifile,*,iostat=i_err) heisout
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('plotenergy')
               read(ifile,*,iostat=i_err) plotenergy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_sparse')
               read(ifile,*,iostat=i_err) do_sparse
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_reduced')
               read(ifile,*,iostat=i_err) do_reduced
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_prnstruct')
               read(ifile,*,iostat=i_err) do_prnstruct
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_prn_poscar')
               read(ifile,*,iostat=i_err) do_prn_poscar
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_hoc_debug')
               read(ifile,*,iostat=i_err) do_hoc_debug
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR DEBUGGING
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR POLARIZATION
            !------------------------------------------------------------------------

            case('do_pol')
               read(ifile,*,iostat=i_err) do_pol
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_loc_pol')
               read(ifile,*,iostat=i_err) do_loc_pol
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_chiral')
               read(ifile,*,iostat=i_err) do_chiral
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('pol_step')
               read(ifile,*,iostat=i_err) pol_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('pol_buff')
               read(ifile,*,iostat=i_err) pol_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('max_pol_nn')
               read(ifile,*,iostat=i_err) max_pol_nn
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR POLARIZATION
            !------------------------------------------------------------------------

            case('logsamp')
               read(ifile,*,iostat=i_err) logsamp
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING THERMAL FIELDS
            !------------------------------------------------------------------------

            case('do_thermfield') ! Flag to print the thermal field contribution
               read(ifile,*,iostat=i_err) do_thermfield
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('thermfield_step') ! Time interval between printing the thermal field
               read(ifile,*,iostat=i_err) thermfield_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('thermfield_buff') ! Buffer size to store thermal field values
               read(ifile,*,iostat=i_err) thermfield_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING THERMAL FIELDS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING TOTAL FIELDS
            !------------------------------------------------------------------------

            case('do_prn_beff') ! Flag to print the total effective field
               read(ifile,*,iostat=i_err) do_prn_beff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('beff_step') ! Time interval between printing the field
               read(ifile,*,iostat=i_err) beff_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('beff_buff') ! Buffer to save the effective field
               read(ifile,*,iostat=i_err) beff_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING TOTAL FIELDS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING STT TORQUES
            !------------------------------------------------------------------------

            case('do_prn_torques')
               read(ifile,*,iostat=i_err) do_prn_torques
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('torques_step') ! Time interval between printing the field
               read(ifile,*,iostat=i_err) torques_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('torques_buff')
               read(ifile,*,iostat=i_err) torques_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING STT TORQUES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING HAMILTONIAN FIELDS
            !------------------------------------------------------------------------

            case('do_prn_binteff') ! Flag to print the total effective field
               read(ifile,*,iostat=i_err) do_prn_binteff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('binteff_step') ! Time interval between printing the field
               read(ifile,*,iostat=i_err) binteff_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('binteff_buff') ! Buffer to save the effective field
               read(ifile,*,iostat=i_err) binteff_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING HAMILTONIAN FIELDS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING LARMOR FREQ
            !------------------------------------------------------------------------

            case('do_larmor_loc') ! Flag to print the total effective field
               read(ifile,*,iostat=i_err) do_larmor_loc
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_larmor_dos') ! Flag to print the total effective field
               read(ifile,*,iostat=i_err) do_larmor_dos
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('larm_step') ! Time interval between printing the field
               read(ifile,*,iostat=i_err) larm_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('larm_buff') ! Buffer to save the effective field
               read(ifile,*,iostat=i_err) larm_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('larm_dos_size') ! Buffer to save the effective field
               read(ifile,*,iostat=i_err) larm_dos_size
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING LARMOR FREQ
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR SKYRMION NUMBER
            !------------------------------------------------------------------------

            case('skyno')
               read(ifile,*,iostat=i_err) skyno
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('skyno_step')
               read(ifile,*,iostat=i_err) skyno_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('skyno_buff')
               read(ifile,*,iostat=i_err) skyno_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_proj_skyno')
               read(ifile,*,iostat=i_err) do_proj_skyno
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_skyno_den')
               read(ifile,*,iostat=i_err) do_skyno_den
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR SKYRMION NUMBER
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR QHB
            !------------------------------------------------------------------------

            case('do_qhb')
               read(ifile,*,iostat=i_err) do_qhb
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qhb_mode')
               read(ifile,*,iostat=i_err) qhb_mode
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('tcurie')
               read(ifile,*,iostat=i_err) tcurie
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR QHB
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR GPU
            !------------------------------------------------------------------------

            case('gpu_mode')
               read(ifile,*,iostat=i_err) gpu_mode
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('gpu_rng')
               read(ifile,*,iostat=i_err) gpu_rng
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('gpu_rng_seed')
               read(ifile,*,iostat=i_err) gpu_rng_seed

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR GPU
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR LSF
            !------------------------------------------------------------------------

            case('do_lsf')
               read(ifile,*,iostat=i_err) do_lsf
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lsf_interpolate')
               read(ifile,*,iostat=i_err) lsf_interpolate
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lsf_field')
               read(ifile,*,iostat=i_err) lsf_field
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lsf_window')
               read(ifile,*,iostat=i_err) lsf_window
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lsf_metric')
               read(ifile,*,iostat=i_err) lsf_metric
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('conf_num')
               read(ifile,*,iostat=i_err) conf_num
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('gsconf_num')
               read(ifile,*,iostat=i_err) gsconf_num
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lsffile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               lsffile=adjustl(trim(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR LSF
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR INDUCED MOMENTS
            !------------------------------------------------------------------------

            case('ind_mom_flag') ! Flag to indicate that there are induced moments being considered
               read(ifile,*,iostat=i_err) ind_mom_flag
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ind_tol') ! Value for the tolerance between neighbouring shells
               read(ifile,*,iostat=i_err) ind_tol
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_prn_induced') ! Whether or not to print extra info for the induced moments
               read(ifile,*,iostat=i_err) do_prn_induced
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ind_buff') ! Buffer for the induced moments
               read(ifile,*,iostat=i_err) ind_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ind_step') ! Steps between measurements for the induced moments
               read(ifile,*,iostat=i_err) ind_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('renorm_coll')  ! If the setup is in some sort non-collinear this forces the first calculation of the susceptibility to be collinear
               read(ifile,*,iostat=i_err) renorm_coll
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR INDUCED MOMENTES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR LOCAL FIELDS
            !------------------------------------------------------------------------

            case('locfield')
               locfield='Y'
               do_bpulse=5
               read(ifile,'(a)',iostat=i_err) cache
               locfieldfile=trim(adjustl(cache))
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('siteatomfield')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               siteatomfile=adjustl(trim(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR LOCAL FIELDS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR TIME DEPENDENT FIELDS
            !------------------------------------------------------------------------

            case('do_bpulse')
               read(ifile,*,iostat=i_err) do_bpulse
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('bpulsefile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               bpulsefile=trim(adjustl(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR TIME DEPENDENT FIELDS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR KMC
            !------------------------------------------------------------------------

            case('do_kmc') ! Perform KMC for the magnetic polaron
               read(ifile,*,iostat=i_err) do_KMC
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('kmc_method') ! Which method for the KMC will be used
               read(ifile,*,iostat=i_err)kmc_method
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('do_efield') ! Consider an electric field for the magnetic polaron
               read(ifile,*,iostat=i_err) do_efield
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('efield') ! the efield it is considered to be constant (commign from to infinite parallel plates)
               read(ifile,*,iostat=i_err) efield(1), efield(2), efield(3)
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('rate0') ! The attempt rate for the kmc process
               read(ifile,*,iostat=i_err) rate0
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('barrfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               barrfile=adjustl(trim(cache))

            case('do_prn_kmc') ! Flag for printing extra info about the KMC process
               read(ifile,*,iostat=i_err) do_prn_kmc
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('kmc_buff') ! Buffer size for measurement steps in the KMC
               read(ifile,*,iostat=i_err) kmc_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('kmc_step') ! Number of time steps between measurement steps in the KMC
               read(ifile,*,iostat=i_err) kmc_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('time_efield') ! Number of time steps between measurement steps in the KMC
               read(ifile,*,iostat=i_err) time_efield
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('kmc_posfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               kmc_posfile=adjustl(trim(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR KMC
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR STIFFNESS
            !------------------------------------------------------------------------

            case('eta_max') ! Maximum convergency parameter for the exchange stiffness
               read(ifile,*,iostat=i_err) eta_max
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('eta_min') ! Minimum convergency parameter for the exchange stiffness
               read(ifile,*,iostat=i_err) eta_min
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('do_stiffness') ! Calculate the exchange stiffness for a ferromagnet
               read(ifile,*,iostat=i_err) do_stiffness
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('do_dm_stiffness') ! Calculate the exchange stiffness for a ferromagnet
               read(ifile,*,iostat=i_err) do_dm_stiffness
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('prn_j0_matrix') ! Print the J0 matrix NAxNA matrix
               read(ifile,*,iostat=i_err) prn_J0_matrix
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR STIFFNESS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR EWALD SUMMATION
            !------------------------------------------------------------------------

            case('do_ewald')
               read(ifile,'(a)',iostat=i_err) do_ewald
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err

            case('ewald_alpha')
               read(ifile,*,iostat=i_err) Ewald_alpha
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err

            case('rmax_ewald')
               read(ifile,*,iostat=i_err) RMAX
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err

            case('kmax_ewald')
               read(ifile,*,iostat=i_err) KMAX
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR EWALD SUMMATION
            !------------------------------------------------------------------------


            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR FREQ SPIN CORRELATION
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF I/O VARIABLES
            !------------------------------------------------------------------------

            case('read_ovf')
               read(ifile,*,iostat=i_err) read_ovf
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('prn_ovf')
               read(ifile,*,iostat=i_err) prn_ovf
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF I/O VARIABLES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF LEGACY VARIABLES
            !------------------------------------------------------------------------

            case('mcavrg_step')
               read(ifile,*,iostat=i_err) mcavrg_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mcavrg_buff')
               read(ifile,*,iostat=i_err) mcavrg_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mompar')
               read(ifile,*,iostat=i_err) mompar
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF LEGACY VARIABLES
            !------------------------------------------------------------------------

            case default
               if(len(trim(keyword))>0) then
                  read(ifile,*)
               end if

            end select
         end if

         ! End of file
         if (i_errb==20) goto 20
         ! End of row
         if (i_errb==10) goto 10
      end do

      20  continue

      return
   end subroutine read_parameters

   !---------------------------------------------------------------------------------
   !> @brief
   !> Read Positions
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine read_positions()
      !
      !
      implicit none
      !
      integer :: flines,itype,mtype,iat,isite,msite,i_stat
      real(dblprec),dimension(3) :: tmp

      ! Open input file
      open(ifileno, file=posfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0
      msite=0
      ! Pre-read file to get max no. sites, types and chemical species
      do
         read(ifileno,*,end=200)  isite,itype
         flines=flines+1
         msite=max(msite,isite)
         mtype=max(mtype,itype)
      end do
      200 continue
      rewind(ifileno)

      ! Set no. sites, types and chemical species
      if(msite/=flines) write(*,*) "WARNING: Check input file ", posfile, &
         " for inconsistent information."
      na=msite
      nt=mtype
      ! Really? I would expect nchmax=1 if no random alloy and if so it will be set later.
      !nchmax=nt
      ! Allocate input arrays
      allocate(bas(3,na),stat=i_stat)
      call memocc(i_stat,product(shape(bas))*kind(bas),'bas','read_positions')
      allocate(atype_inp(na),stat=i_stat)
      call memocc(i_stat,product(shape(atype_inp))*kind(atype_inp),'atype_inp','read_positions')
      allocate(anumb_inp(na),stat=i_stat)
      call memocc(i_stat,product(shape(anumb_inp))*kind(anumb_inp),'anumb_inp','read_positions')

      ! Read basis atoms and setup type array
      ! Site, Type, Rx, Ry, Rz
      if (posfiletype=='C') then
         do iat=1, na
            read (ifileno,*) isite, itype,  bas(1,isite), bas(2,isite),&
            bas(3,isite)
            atype_inp(isite)=itype

            ! Redundant but kept for the time beeing
            anumb_inp(isite)=isite
         enddo
      elseif (posfiletype=='D') then
         do iat=1, na
            read (ifileno,*) isite, itype,  tmp(1),tmp(2),tmp(3)
            atype_inp(isite)=itype
            bas(1,isite)=tmp(1)*C1(1)+tmp(2)*C2(1)+tmp(3)*C3(1)
            bas(2,isite)=tmp(1)*C1(2)+tmp(2)*C2(2)+tmp(3)*C3(2)
            bas(3,isite)=tmp(1)*C1(3)+tmp(2)*C2(3)+tmp(3)*C3(3)

            ! Redundant but kept for the time beeing
            anumb_inp(isite)=isite
         enddo
      endif
      if (maptype==3) then
         ! Raw values of basis coordinates stored to bas0,  needed for RSPt
         allocate(bas0(3,na),stat=i_stat)
         call memocc(i_stat,product(shape(bas0))*kind(bas0),'bas0','read_positions')
         bas0(1:3,1:na)=bas(1:3,1:na)

      endif
               
      close (ifileno)
   end subroutine read_positions

   !---------------------------------------------------------------------------------
   !> @brief
   !> Read Positions for random alloys
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine read_positions_alloy()
      !
      !
      implicit none
      !
      !
      integer :: flines,isite,itype,ichem
      integer :: msite,mtype,mchem,iat, i_stat
      real(dblprec) :: rconc
      real(dblprec),dimension(3) :: tmp

      ! Open input file
      open(ifileno, file=posfile)

      ! Check if input file is for random alloy
      flines=0
      msite=0
      mtype=0
      mchem=0

      ! Pre-read file to get max no. sites, types and chemical species
      do
         read(ifileno,*,end=200)  isite,itype,ichem
         flines=flines+1
         msite=max(msite,isite)
         mtype=max(mtype,itype)
         mchem=max(mchem,ichem)
      end do
      200 continue
      rewind(ifileno)

      ! Set no. sites, types and chemical species
      na=msite
      nt=mtype
      nchmax=mchem

      ! Allocate input arrays
      allocate(bas(3,msite),stat=i_stat)
      call memocc(i_stat,product(shape(bas))*kind(bas),'bas','read_positions_alloy')
      bas=0.0_dblprec
      allocate(atype_inp(na),stat=i_stat)
      call memocc(i_stat,product(shape(atype_inp))*kind(atype_inp),'atype_inp','read_positions_alloy')
      atype_inp=0
      allocate(anumb_inp(na),stat=i_stat)
      call memocc(i_stat,product(shape(anumb_inp))*kind(anumb_inp),'anumb_inp','read_positions_alloy')
      anumb_inp=0
      ! Chemical data
      allocate(nch(na),stat=i_stat)
      call memocc(i_stat,product(shape(nch))*kind(nch),'nch','read_positions_alloy')
      nch=0
      allocate(chconc(na,nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(chconc))*kind(chconc),'chconc','read_positions_alloy')
      chconc=0.0_dblprec

      ! Read data
      ! Site  Type   Chem_type Conc, Rx, Ry, Rz
      if (posfiletype=='C') then
         do iat=1, flines
            read (ifileno,*) isite, itype, ichem, rconc, bas(1,isite), bas(2,isite), bas(3,isite)
            nch(isite)=max(nch(isite),ichem)
            atype_inp(isite)=itype
            anumb_inp(isite)=isite
            chconc(isite,ichem)=rconc
         enddo
      elseif(posfiletype=='D') then
         do iat=1, flines
            read (ifileno,*) isite, itype,ichem,rconc,  tmp(1),tmp(2),tmp(3)
            bas(1,isite)=tmp(1)*C1(1)+tmp(2)*C2(1)+tmp(3)*C3(1)
            bas(2,isite)=tmp(1)*C1(2)+tmp(2)*C2(2)+tmp(3)*C3(2)
            bas(3,isite)=tmp(1)*C1(3)+tmp(2)*C2(3)+tmp(3)*C3(3)
            nch(isite)=max(nch(isite),ichem)
            atype_inp(isite)=itype
            anumb_inp(isite)=isite
            chconc(isite,ichem)=rconc
         enddo
      endif
      close (ifileno)
   end subroutine read_positions_alloy


   !---------------------------------------------------------------------------------
   !> @brief Change constants to atomic units
   !---------------------------------------------------------------------------------
   subroutine change_constants()
      use Parameters
      use Constants

      implicit none

      !.. Scalar parameters
      gama        = 1.0_dblprec
      k_bolt      = 1.0_dblprec
      mub         = 1.0_dblprec
      mu0         = 1.0_dblprec
      mry         = 1.0_dblprec
      hbar        = 1.0_dblprec
      hbar_mev    = 1.0_dblprec
      k_bolt_ev   = 1.0_dblprec
      ry_ev       = 1.0_dblprec
      amu         = 1.0_dblprec
      angstrom    = 1.0_dblprec
      ev          = 1.0_dblprec
   end subroutine change_constants


   !--------------------------------------------------------------------------------
   !> @brief
   !> Read Magnetic moments
   !
   !> @author
   !> Anders Bergman
   !> @date 23/02/2015 - Jonathan Chico
   !> - Introducing the capacity to read whether a moment is induced or not
   !--------------------------------------------------------------------------------
   subroutine read_moments(Landeg_global)
      use Parameters
      use Profiling
      use LSF
      implicit none

      real(dblprec), intent(in) :: Landeg_global !< Default gyromagnetic ratio
      !
      integer :: i_err,isite,ichem,i_stat,iconf
      real(dblprec)  :: aemom_tmp

      iconf = 1

      open(ifileno,file=trim(momfile))

      !Allocate arrays according to data from position input
      allocate(ammom_inp(na,nchmax,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(ammom_inp))*kind(ammom_inp),'ammom_inp','read_moments')
      ammom_inp=0.0_dblprec

      allocate(aemom_inp(3,na,nchmax,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(aemom_inp))*kind(aemom_inp),'aemom_inp','read_moments')
      aemom_inp=0.0_dblprec

      allocate(Landeg_ch(na,nchmax,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(Landeg_ch))*kind(Landeg_ch),'Landeg_ch','read_moments')
      Landeg_ch=0.0_dblprec
      if (ind_mom_flag=='Y') then
         allocate(ind_mom(na,nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ind_mom))*kind(ind_mom),'ind_mom','read_moments')
         ind_mom=0
      endif

      i_err=0

      if(set_landeg==1) then
         ! If the induced magnetic moments flag is on one must read whether a certain moment is induced or not
         if (do_lsf=='N') then
            if (ind_mom_flag=='Y') then
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite, ichem, ammom_inp(isite,ichem,1),&
                     aemom_inp(1:3,isite,ichem,1), Landeg_ch(isite,ichem,1),ind_mom(isite,ichem)
                     aemom_tmp=norm2(aemom_inp(:,isite,ichem,1))
                     aemom_inp(:,isite,ichem,1)=aemom_inp(:,isite,ichem,1)/aemom_tmp
               end do
            else
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite, ichem, ammom_inp(isite,ichem,1), &
                     aemom_inp(1:3,isite,ichem,1), Landeg_ch(isite,ichem,1)
                     aemom_tmp=norm2(aemom_inp(:,isite,ichem,1))
                     aemom_inp(:,isite,ichem,1)=aemom_inp(:,isite,ichem,1)/aemom_tmp
               end do
            endif
         else ! LSF
            ! For LSF modified momfile requires configuration number as first column
            do while(i_err==0)
               read(ifileno,*,iostat=i_err) iconf, isite, ichem, ammom_inp(isite,ichem,iconf), &
                  aemom_inp(1:3,isite,ichem,iconf), Landeg_ch(isite,ichem,iconf)
                  aemom_tmp=norm2(aemom_inp(:,isite,ichem,iconf))
                  aemom_inp(:,isite,ichem,iconf)=aemom_inp(:,isite,ichem,iconf)/aemom_tmp
            end do
         endif
      else
         Landeg_ch=Landeg_global
         if (do_lsf=='N') then
            if (ind_mom_flag=='Y') then
               ! If the induced magnetic moments flag is on one must read whether a certain moment is induced or not
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite, ichem, ammom_inp(isite,ichem,1), &
                     aemom_inp(1:3,isite,ichem,1), ind_mom(isite,ichem)
                     aemom_tmp=norm2(aemom_inp(:,isite,ichem,1))
                     aemom_inp(:,isite,ichem,1)=aemom_inp(:,isite,ichem,1)/aemom_tmp
               end do
            else
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite, ichem, ammom_inp(isite,ichem,1), &
                     aemom_inp(1:3,isite,ichem,1)
                     aemom_tmp=norm2(aemom_inp(:,isite,ichem,1))
                     aemom_inp(1:3,isite,ichem,1)=aemom_inp(1:3,isite,ichem,1)/aemom_tmp
               end do
            endif
         else   ! LSF
            do while(i_err==0)
               read(ifileno,*,iostat=i_err) iconf, isite, ichem, ammom_inp(isite,ichem,iconf), &
                  aemom_inp(1:3,isite,ichem,iconf)
                  aemom_tmp=norm2(aemom_inp(:,isite,ichem,iconf))
                  aemom_inp(:,isite,ichem,iconf)=aemom_inp(:,isite,ichem,iconf)/aemom_tmp
            end do
         endif
      end if
      close(ifileno)
      !
   end subroutine read_moments

   !----------------------------------------------------------------------------------
   !  SUBROUTINE: read_fixed_moments
   !> @brief
   !> Read Magnetic moments in the case that some of the moments are kept fixed
   !
   !> @author
   !> Jonathan Chico
   !> Based in the routine by Anders Bergman, modified to deal with fixed moments
   !----------------------------------------------------------------------------------
   subroutine read_fixed_moments(Landeg_global)

      use LSF
      use FixedMom, only : inp_fixed_mom_flag
      use Profiling
      use Parameters

      implicit none

      real(dblprec), intent(in) :: Landeg_global !< Default gyromagnetic ratio
      !
      integer :: i_err,isite,ichem,i_stat,iconf
      real(dblprec)  :: aemom_tmp

      iconf = 1

      open(ifileno,file=trim(momfile))

      !Allocate arrays according to data from position input
      allocate(ammom_inp(NA,Nchmax,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(ammom_inp))*kind(ammom_inp),'ammom_inp','read_fixed_moments')
      ammom_inp=0.0_dblprec

      allocate(aemom_inp(3,NA,Nchmax,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(aemom_inp))*kind(aemom_inp),'aemom_inp','read_fixed_moments')
      aemom_inp=0.0_dblprec

      allocate(Landeg_ch(NA,Nchmax,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(Landeg_ch))*kind(Landeg_ch),'Landeg_ch','read_fixed_moments')
      Landeg_ch=0.0_dblprec

      allocate(inp_fixed_mom_flag(NA,Nchmax,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(inp_fixed_mom_flag))*kind(inp_fixed_mom_flag),'inp_fixed_mom_flag','read_fixed_moments')
      inp_fixed_mom_flag=0

      if (ind_mom_flag=='Y') then
         allocate(ind_mom(na,nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ind_mom))*kind(ind_mom),'ind_mom','read_moments')
         ind_mom=0
      endif
      i_err=0

      if(set_landeg==1) then
         if (do_lsf=='N') then
            ! If there are induced moments
            if (ind_mom_flag=='Y') then
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite, ichem, ammom_inp(isite,ichem,1), &
                     aemom_inp(1:3,isite,ichem,1), Landeg_ch(isite,ichem,1), inp_fixed_mom_flag(isite,ichem,1),ind_mom(isite,ichem)
                     aemom_tmp=norm2(aemom_inp(:,isite,ichem,1))
                     aemom_inp(:,isite,ichem,1)=aemom_inp(:,isite,ichem,1)/aemom_tmp
               end do
            else ! No induced moment
               read(ifileno,*,iostat=i_err) isite, ichem, ammom_inp(isite,ichem,1), &
                  aemom_inp(1:3,isite,ichem,1), Landeg_ch(isite,ichem,1), inp_fixed_mom_flag(isite,ichem,1)
                  aemom_tmp=norm2(aemom_inp(:,isite,ichem,1))
                  aemom_inp(:,isite,ichem,1)=aemom_inp(:,isite,ichem,1)/aemom_tmp
            endif
         else ! LSF
            ! For LSF modified momfile requires configuration number as first column
            do while(i_err==0)
               read(ifileno,*,iostat=i_err) iconf, isite, ichem, ammom_inp(isite,ichem,iconf),&
                  aemom_inp(1:3,isite,ichem,iconf), Landeg_ch(isite,ichem,iconf), inp_fixed_mom_flag(isite,ichem,iconf)
                  aemom_tmp=norm2(aemom_inp(:,isite,ichem,iconf))
                  aemom_inp(:,isite,ichem,iconf)=aemom_inp(:,isite,ichem,iconf)/aemom_tmp
            end do
         endif
      else
         Landeg_ch=Landeg_global
         if (do_lsf=='N') then
            ! Induced moment
            if (ind_mom_flag=='Y') then
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite, ichem, ammom_inp(isite,ichem,1), &
                     aemom_inp(1:3,isite,ichem,1), inp_fixed_mom_flag(isite,ichem,1),ind_mom(isite,ichem)
                     aemom_tmp=norm2(aemom_inp(:,isite,ichem,1))
                     aemom_inp(1:3,isite,ichem,1)=aemom_inp(1:3,isite,ichem,1)/aemom_tmp
               end do
            else ! No induced moment
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite, ichem, ammom_inp(isite,ichem,1), &
                     aemom_inp(1:3,isite,ichem,1), inp_fixed_mom_flag(isite,ichem,1)
                     aemom_tmp=norm2(aemom_inp(:,isite,ichem,1))
                     aemom_inp(1:3,isite,ichem,1)=aemom_inp(1:3,isite,ichem,1)/aemom_tmp
               end do
            endif
         else   ! LSF
            do while(i_err==0)
               read(ifileno,*,iostat=i_err) iconf, isite, ichem, ammom_inp(isite,ichem,iconf), &
                  aemom_inp(1:3,isite,ichem,iconf), inp_fixed_mom_flag(isite,ichem,iconf)
                  aemom_tmp=norm2(aemom_inp(:,isite,ichem,iconf))
                  aemom_inp(:,isite,ichem,iconf)=aemom_inp(:,isite,ichem,iconf)/aemom_tmp
            end do
         endif
      end if
      close(ifileno)
      !
   end subroutine read_fixed_moments

   !--------------------------------------------------------------------------------
   !> @brief
   !> Reads exchange from input file
   !>
   !> @author
   !> Anders Bergman
   !>
   !> @date 09/16/2014 - Thomas Nystrand
   !> - Splitting into subroutines and fixed bugs
   !> @date Feb 2017
   !> - Allowing reading of two different sets of files
   !> @date May 2017
   !> -  Jonathan Chico ---> Allowing for reading jfile for the cluster
   !--------------------------------------------------------------------------------
   subroutine read_exchange()
      !! @todo Consider change so that if the atom basis is given in direct coordinates, then also
      !! @todo the coordinates for the exchange coupling shells are to be stated in direct coordinates.
      !
      use clusters
      implicit none

      integer :: itype, jtype, isite, jsite, ichem, jchem, iline, ishell
      integer :: itype_clus, jtype_clus,isite_c,jsite_c,jchem_c,ichem_c
      integer :: flines,no_shells,ii,idum,flines_clus,no_shells_clus
      real(dblprec), dimension(3) :: r_red, r_tmp,r_tmp_clus,r_red_clus
      logical :: unique
      real(dblprec):: j_tmp,j_tmpD,tmp,j_tmp_clus
      real(dblprec):: tol, norm
      integer :: ifileno2,ifileno3

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      !jc = 0.0_dblprec
      ifileno2=ifileno+1
      ifileno3=ifileno+2

      do ii=1,conf_num
         open(ifileno, file=trim(jfile(ii))) ! Number of shells same for different LSF configuration
         if (exc_inter=='Y') open(ifileno2, file=trim(jfileD(ii))) ! Number of shells same for different LSF configuration
         if (ii==1) then
            call read_exchange_getMaxNoShells(no_shells,flines,jfile(1))
            call allocate_hamiltonianinput(no_shells,1)
         endif
         redcoord = 0.0_dblprec
         nn       = 0

         ! Read exchange vectors
         ! Isite, Jsite, Ichem, Jchem
         do iline=1, flines
            ! Loop through earlier vectors to find equivalent shells

            ! Read indices and coordinates
            if(do_ralloy==0) then
               read (ifileno,*) isite, jsite, r_tmp(1:3), j_tmp
               if (exc_inter=='Y') read (ifileno2,*) idum, idum, tmp, tmp, tmp, j_tmpD
               ichem=1
               jchem=1
            else
               read (ifileno,*) isite, jsite, ichem, jchem, r_tmp(1:3), j_tmp
               if (exc_inter=='Y') read (ifileno2,*) idum, idum, idum, idum, tmp, tmp, tmp, j_tmpD
            end if

            itype=atype_inp(isite)
            jtype=1
            call read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

            ! Loop through earlier vectors to find equivalent shells
            unique=.true.
            do ishell=1,nn(itype)
               norm=(r_red(1)-redcoord(itype,ishell,1))**2+ &
                  (r_red(2)-redcoord(itype,ishell,2))**2+ &
                  (r_red(3)-redcoord(itype,ishell,3))**2
               if(norm<tol) then
                  unique=.false.
                  ! If neighbour vector already exist, replace J value (could be removed)
                  if(do_ralloy==0) then
                     jc(itype,ishell,ichem,jtype,ii)=j_tmp
                     if (exc_inter=='Y') jcD(itype,ishell,ichem,jtype,ii)=j_tmpD
                  else
                     jc(itype,ishell,ichem,jchem,ii)=j_tmp
                     if (exc_inter=='Y') jcD(itype,ishell,ichem,jchem,ii)=j_tmpD
                  end if
               end if
            end do
            if (unique) then
               ! Add entry if not found earlier
               nn(itype)=nn(itype)+1
               redcoord(itype,nn(itype),1:3)=r_red(1:3)
               if(do_ralloy==0) then
                  jc(itype,nn(itype),ichem,jtype,ii)=j_tmp
                  if (exc_inter=='Y') jcD(itype,nn(itype),ichem,jtype,ii)=j_tmpD
               else
                  jc(itype,nn(itype),ichem,jchem,ii)=j_tmp
                  if (exc_inter=='Y') jcD(itype,nn(itype),ichem,jchem,ii)=j_tmpD
               end if
            end if

         enddo
         close(ifileno)
         ! Reading the cluster files
         if (do_cluster=='Y') then

            if (do_cluster=='Y'.and.exc_inter=='Y')then
               open(ifileno3, file=trim(jfile_clus(ii))) ! File for the interactions inside the embeded cluster
            else if (do_cluster=='Y') then
               open(ifileno2, file=trim(jfile_clus(ii))) ! File for the interactions inside the embeded cluster
            endif
            call read_exchange_getMaxNoShells_clus(no_shells_clus,flines_clus,ifileno2,do_ralloy)
            call allocate_hamiltonianinput_clus(conf_num,no_shells_clus,1)

            do iline=1, flines_clus
               ! Loop through earlier vectors to find equivalent shells

               ! Read indices and coordinates
               if(do_ralloy==0) then
                  if(exc_inter=='Y') then
                     read (ifileno3,*) isite_c, jsite_c, r_tmp_clus(1:3), j_tmp_clus
                  else
                     read (ifileno2,*) isite_c, jsite_c, r_tmp_clus(1:3), j_tmp_clus
                  endif
                  ichem_c=1
                  jchem_c=1
               else
                  if (exc_inter=='Y') then
                     read (ifileno3,*) isite_c, jsite_c, ichem_c, jchem_c, r_tmp_clus(1:3), j_tmp_clus
                  else
                     read (ifileno2,*) isite_c, jsite_c, ichem_c, jchem_c, r_tmp_clus(1:3), j_tmp_clus
                  endif
               end if

               itype_clus=atype_inp_clus(isite_c)
               jtype_clus=1
               call read_exchange_getNeighVec_clus(r_red_clus,r_tmp_clus,isite_c,jsite_c,&
               maptype,posfiletype)
               ! Loop through earlier vectors to find equivalent shells
               unique=.true.
               do ishell=1,NN_clus(itype_clus)
                  norm=(r_red_clus(1)-redcoord_clus(itype_clus,ishell,1))**2+ &
                  (r_red_clus(2)-redcoord_clus(itype_clus,ishell,2))**2+ &
                  (r_red_clus(3)-redcoord_clus(itype_clus,ishell,3))**2
                  if(norm<tol) then
                     unique=.false.
                     ! If neighbour vector already exist, replace J value (could be removed)
                     if(do_ralloy==0) then
                        jc_clus(itype_clus,ishell,ichem_c,jtype_clus,ii)=j_tmp_clus
                     else
                        jc_clus(itype_clus,ishell,ichem_c,jchem_c,ii)=j_tmp_clus
                     end if
                  end if
               end do
               if (unique) then
                  ! Add entry if not found earlier
                  NN_clus(itype_clus)=NN_clus(itype_clus)+1
                  redcoord_clus(itype_clus,NN_clus(itype_clus),1:3)=r_red_clus(1:3)
                  if(do_ralloy==0) then
                     jc_clus(itype_clus,NN_clus(itype_clus),ichem_c,jtype_clus,ii)=j_tmp_clus
                  else
                     jc_clus(itype_clus,NN_clus(itype_clus),ichem_c,jchem_c,ii)=j_tmp_clus
                  end if
               end if
            enddo
         endif

         ! Reducing jc size if max_no_shells are small enough !
         max_no_shells=maxval(NN)
         if(exc_inter=='Y') close(ifileno2)
         if(do_cluster=='Y') then
            max_no_shells_clus=maxval(NN_clus)
            close(ifileno3)
         endif
      enddo

      call read_exchange_reduceRedCoordMatrixSize(redcoord,nt,max_no_shells)
      call read_exchange_reduceCouplingMatrixSize(jc,nt,max_no_shells,nchmax)
      if (exc_inter=='Y') call read_exchange_reduceCouplingMatrixSize(jcD,nt,max_no_shells,nchmax)
      if (do_cluster=='Y') then
         call read_exchange_reduceRedCoordMatrixSize(redcoord_clus,NT_clus,max_no_shells_clus)
         call read_exchange_reduceCouplingMatrixSize(jc_clus,NT_clus,max_no_shells_clus,nchmax)
      endif
   end subroutine read_exchange

   !--------------------------------------------------------------------------------
   !> @brief
   !> Reads the exchange parameters from coupling jfile
   !> Put these values in a tensor
   !>
   !> @author
   !> Thomas Nystrand
   !--------------------------------------------------------------------------------
   subroutine read_exchange_tensor()
      logical :: tensor_format = .true.
      call read_exchange_tensor_base(tensor_format)
   end subroutine read_exchange_tensor

   !--------------------------------------------------------------------------------
   !> @brief
   !> Interface for read_exchange_tensor_base
   !> given that a tensor will be built from coupling file
   !>
   !> @author
   !> Thomas Nystrand
   !--------------------------------------------------------------------------------
   subroutine read_exchange_build_tensor()
      logical :: tensor_format = .false.
      call read_exchange_tensor_base(tensor_format)
   end subroutine read_exchange_build_tensor

   !--------------------------------------------------------------------------------
   !> @brief
   !> Reads the exchange tensor parameters from file
   !> Tries to limit the coupling matrix size as well
   !>
   !> @author
   !> Anders Bergman
   !>
   !> @date 09/16/2014 - Thomas Nystrand
   !> - Splitting into subroutines and fixed bugs
   !> - Added tensor construction from exchange option
   !--------------------------------------------------------------------------------
   subroutine read_exchange_tensor_base(tensor_format)
      implicit none
      logical, intent(in) :: tensor_format

      integer       :: itype,jtype,isite,jsite,ichem,jchem,iline,ishell

      integer       :: flines,no_shells
      logical       :: unique
      real(dblprec) :: j_tmp(3,3)
      real(dblprec) :: tol, norm
      real(dblprec) :: j_tmpSingle
      real(dblprec), dimension(3) :: r_tmp, r_red

      ! Set tolerance for neighbour shells
      tol=1.0d-5

     open(ifileno, file=trim(jfile(1))) ! Number of shells same for different LSF configuration

      call read_exchange_getMaxNoShells(no_shells,flines,jfile(1))
      call allocate_hamiltonianinput(no_shells,1)
      redcoord = 0.0_dblprec
      jc_tens  = 0.0_dblprec
      nn       = 0

      ! Read exchange vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines
         if(do_ralloy==0) then
            if(tensor_format) then
               read(ifileno,*) isite,jsite,r_tmp(1:3),j_tmp
            else
               read(ifileno,*) isite,jsite,r_tmp(1:3),j_tmpSingle
            endif
            ichem = 1
            jchem = 1
         else
            call ErrorHandling_ERROR('Random alloy treatment is currently disabled for '//char(13)//char(11)//char(0)// &
               ' reading exchange tensor from file')
         end if

         ! Find type of site
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)
         call read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,nn(itype)
            norm=(r_red(1)-redcoord(itype,ishell,1))**2+ &
            (r_red(2)-redcoord(itype,ishell,2))**2+ &
            (r_red(3)-redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(tensor_format) then
                  jc_tens(1:3,1:3,itype,ishell,ichem,jtype)=j_tmp(1:3,1:3)
               else
                  jc_tens(1:3,1:3,itype,ishell,ichem,jtype)= &
                  reshape((/j_tmpSingle,0.0_dblprec,0.0_dblprec, 0.0_dblprec,j_tmpSingle,0.0_dblprec, 0.0_dblprec,0.0_dblprec,j_tmpSingle/),(/3,3/))
               endif
            end if
         end do
         if (unique) then
            nn(itype)=nn(itype)+1
            redcoord(itype,nn(itype),1:3)=r_red(1:3)
            if(tensor_format) then
               jc_tens(1:3,1:3,itype,ishell,ichem,jtype)=j_tmp(1:3,1:3)
            else
               jc_tens(1:3,1:3,itype,ishell,ichem,jtype)= &
               reshape((/j_tmpSingle,0.0_dblprec,0.0_dblprec, 0.0_dblprec,j_tmpSingle,0.0_dblprec, 0.0_dblprec,0.0_dblprec,j_tmpSingle/),(/3,3/))
            endif
         end if
      enddo
      close (ifileno)

      max_no_shells=maxval(NN)
      call read_exchange_reduceRedCoordMatrixSize(redcoord,nt,max_no_shells)
      call read_exchange_reduceCouplingMatrixSizeTensor(jc_tens,nt,max_no_shells,nchmax)

      close(ifileno)

   end subroutine read_exchange_tensor_base

   !--------------------------------------------------------------------------------
   !> @brief
   !> Trying to reduce the redcoord matrix size
   !>
   !> @author
   !> Thomas Nystrand
   !--------------------------------------------------------------------------------
   subroutine read_exchange_reduceRedCoordMatrixSize(redcoord,nt,max_no_shells)
      integer, intent(in)                                         :: nt, max_no_shells
      real(dblprec), intent(inout), dimension(:,:,:), allocatable :: redcoord

      ! locals
      integer                                       :: i_stat,i_all
      real(dblprec), dimension(:,:,:), allocatable  :: redcoord_tmp

      allocate(redcoord_tmp(nt,max_no_shells,3),stat=i_stat)
      call memocc(i_stat,product(shape(redcoord_tmp))*kind(redcoord_tmp),'redcoord_tmp','read_exchange_reduceRedCoordMatrixSize')
      redcoord_tmp = redcoord(:,1:max_no_shells,:)

      i_all=-product(shape(redcoord))*kind(redcoord)
      deallocate(redcoord,stat=i_stat)
      call memocc(i_stat,i_all,'redcoord','read_exchange_reduceRedCoordMatrixSize')

      allocate(redcoord(nt,max_no_shells,3),stat=i_stat)
      call memocc(i_stat,product(shape(redcoord))*kind(redcoord),'redcoord','read_exchange_reduceRedCoordMatrixSize')

      redcoord = redcoord_tmp

      i_all=-product(shape(redcoord_tmp))*kind(redcoord_tmp)
      deallocate(redcoord_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'redcoord_tmp','read_exchange_reduceRedCoordMatrixSize')

      close (ifileno)
   end subroutine read_exchange_reduceRedCoordMatrixSize

   !--------------------------------------------------------------------------------
   !> @brief
   !> Trying to reduce the coupling matrix size
   !>
   !> @author
   !> Thomas Nystrand
   !--------------------------------------------------------------------------------
   subroutine read_exchange_reduceCouplingMatrixSize(jc,nt,max_no_shells,nchmax)
      integer, intent(in)          :: nt, max_no_shells, nchmax
      real(dblprec), intent(inout), dimension(:,:,:,:,:), allocatable :: jc

      ! locals
      integer :: i_stat,i_all
      real(dblprec), dimension(:,:,:,:,:), allocatable     :: jc_tmp

      allocate(jc_tmp(nt,max_no_shells,nchmax,nchmax,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(jc_tmp))*kind(jc_tmp),'jc_tmp','read_exchange_reduceCouplingMatrixSize')
      jc_tmp=0.0_dblprec

      jc_tmp = jc(:,1:max_no_shells,:,:,:)

      i_all=-product(shape(jc))*kind(jc)
      deallocate(jc,stat=i_stat)
      call memocc(i_stat,i_all,'jc','read_exchange_reduceCouplingMatrixSize')

      allocate(jc(nt,max_no_shells,nchmax,nchmax,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(jc))*kind(jc),'jc','read_exchange_reduceCouplingMatrixSize')
      jc=0.0_dblprec

      jc = jc_tmp

      i_all=-product(shape(jc_tmp))*kind(jc_tmp)
      deallocate(jc_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'jc_tmp','read_exchange_reduceCouplingMatrixSize')

      close (ifileno)
   end subroutine read_exchange_reduceCouplingMatrixSize

   !--------------------------------------------------------------------------------
   !> @brief
   !> Trying to reduce the coupling matrix size
   !>
   !> @author
   !> Thomas Nystrand
   !--------------------------------------------------------------------------------
   subroutine read_dmexchange_reduceCouplingMatrixSize(jc,nt,max_no_shells,nchmax)
      integer, intent(in)          :: nt, max_no_shells, nchmax
      real(dblprec), intent(inout), dimension(:,:,:,:,:), allocatable :: jc

      ! locals
      integer :: i_stat,i_all
      real(dblprec), dimension(:,:,:,:,:), allocatable     :: jc_tmp

      allocate(jc_tmp(3,nt,max_no_shells,nchmax,nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(jc_tmp))*kind(jc_tmp),'jc_tmp','read_exchange_reduceCouplingMatrixSize')
      jc_tmp=0.0_dblprec

      jc_tmp = jc(1:3,:,1:max_no_shells,:,:)

      i_all=-product(shape(jc))*kind(jc)
      deallocate(jc,stat=i_stat)
      call memocc(i_stat,i_all,'jc','read_exchange_reduceCouplingMatrixSize')

      allocate(jc(3,nt,max_no_shells,nchmax,nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(jc))*kind(jc),'jc','read_exchange_reduceCouplingMatrixSize')
      jc=0.0_dblprec

      jc = jc_tmp

      i_all=-product(shape(jc_tmp))*kind(jc_tmp)
      deallocate(jc_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'jc_tmp','read_exchange_reduceCouplingMatrixSize')

      close (ifileno)
   end subroutine read_dmexchange_reduceCouplingMatrixSize

   !--------------------------------------------------------------------------------
   !> @brief
   !> Trying to reduce the tensor coupling matrix size
   !>
   !> @author
   !> Thomas Nystrand
   !--------------------------------------------------------------------------------
   subroutine read_exchange_reduceCouplingMatrixSizeTensor(jc,nt,max_no_shells,nchmax)
      integer, intent(in)          :: nt, max_no_shells, nchmax
      real(dblprec), intent(inout), dimension(:,:,:,:,:,:), allocatable :: jc

      ! locals
      integer :: i_stat,i_all
      real(dblprec), dimension(:,:,:,:,:,:), allocatable :: jc_tmp

      allocate(jc_tmp(3,3,nt,max_no_shells,nchmax,nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(jc_tmp))*kind(jc_tmp),'jc_tmp','read_exchange_reduceCouplingMatrixSizeTensor')
      jc_tmp=0.0_dblprec

      jc_tmp = jc(:,:,:,1:max_no_shells,:,:)

      i_all=-product(shape(jc))*kind(jc)
      deallocate(jc,stat=i_stat)
      call memocc(i_stat,i_all,'jc','read_exchange_reduceCouplingMatrixSizeTensor')

      allocate(jc(3,3,nt,max_no_shells,nchmax,nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(jc))*kind(jc),'jc','read_exchange_reduceCouplingMatrixSizeTensor')
      jc=0.0_dblprec

      jc = jc_tmp

      i_all=-product(shape(jc_tmp))*kind(jc_tmp)
      deallocate(jc_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'jc_tmp','read_exchange_reduceCouplingMatrixSizeTensor')

      close (ifileno)
   end subroutine read_exchange_reduceCouplingMatrixSizeTensor

   !--------------------------------------------------------------------------------
   !> @brief
   !> Get the max no of exchange shells and lines
   !> Helper for read exchange
   !>
   !> @author
   !> Anders Bergman
   !--------------------------------------------------------------------------------
   subroutine read_exchange_getMaxNoShells(no_shells,flines,filename)

      implicit none

      integer, intent(out)                   :: no_shells,flines
      character(len=*), intent(in)           :: filename
      integer                                :: mtype
      integer                                :: itype,jtype,isite,jsite,ichem,jchem
      integer                                :: i_stat,i_all

      integer, dimension(:,:,:), allocatable :: nn_tmp

      flines=0
      mtype=0

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_exchange')
      nn_tmp=0

      ! Pre-read file to get max no. exchange shells and no. lines
      do
         if(do_ralloy==0) then
            read(ifileno,*,end=200)  isite,jsite
            ichem=1
            jchem=1
         else
            read(ifileno,*,end=200)  isite, jsite, ichem, jchem
         end if
         flines=flines+1
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)
         if(do_ralloy==0) then
            nn_tmp(itype,jtype,jchem)=nn_tmp(itype,jtype,jchem)+1
         else
            nn_tmp(itype,ichem,jchem)=nn_tmp(itype,ichem,jchem)+1
         end if
      end do
      200 continue

      rewind(ifileno)

      no_shells=0
      if (do_ralloy==0) then
         do itype=1,nt
            no_shells=max(sum(nn_tmp(itype,:,:)),no_shells)
         end do
      else
         do itype=1,nt
            do ichem=1,Nchmax
               do jchem=1,Nchmax
                  no_shells=max(sum(nn_tmp(itype,ichem,:)),no_shells)
               end do
            end do
         end do
      endif
      i_all=-product(shape(nn_tmp))*kind(nn_tmp)
      deallocate(nn_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'NN_tmp','read_exchange')

   end subroutine read_exchange_getMaxNoShells

   !--------------------------------------------------------------------------------
   !> @brief
   !> Obtaining the neighbour vector
   !>
   !> @author
   !> Anders Bergman
   !--------------------------------------------------------------------------------
   subroutine read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

      implicit none

      real(dblprec), dimension(3), intent(out) :: r_red
      integer, intent(in)                      :: isite,jsite
      real(dblprec), dimension(3), intent(in)  :: r_tmp

      real(dblprec), dimension(3) :: r_tmp2

      if(maptype==3) then
         ! Calculate proper neighbour vector (from "RSPt")
         r_tmp2=r_tmp
         r_red(1)=Bas0(1,jsite)-Bas0(1,isite)+C1(1)*r_tmp2(1)+C2(1)*r_tmp2(2)+C3(1)*r_tmp2(3)
         r_red(2)=Bas0(2,jsite)-Bas0(2,isite)+C1(2)*r_tmp2(1)+C2(2)*r_tmp2(2)+C3(2)*r_tmp2(3)
         r_red(3)=Bas0(3,jsite)-Bas0(3,isite)+C1(3)*r_tmp2(1)+C2(3)*r_tmp2(2)+C3(3)*r_tmp2(3)
         !write(1357,'(2i8,5f12.6)') isite, jsite, r_tmp2(1:3), norm2(r_red(1:3))
      elseif(maptype==2) then
         ! Calculate proper neighbour vector (from "bgfm")
         r_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(1)+C2(1)*r_tmp(2)+C3(1)*r_tmp(3)
         r_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(1)+C2(2)*r_tmp(2)+C3(2)*r_tmp(3)
         r_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(1)+C2(3)*r_tmp(2)+C3(3)*r_tmp(3)
      else
         ! Calculates neighbour vectors from direct coordinates or Cartesian
         ! coordinates, corresponding to how the atomic positions are entered
         !if (posfiletype=='C') then
            r_red=r_tmp
         !elseif (posfiletype=='D') then
         !   r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
         !   r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
         !   r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
         !else
         !   stop 'Only posfiletype = C or D is currently supported'
         !endif
      end if
   end subroutine read_exchange_getNeighVec

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_anisotropy_alloy
   !> @brief Read the anisotropy variables for random alloys
   !---------------------------------------------------------------------------------
   subroutine read_anisotropy_alloy()
      !
      !
      implicit none
      !
      integer :: iat,i_err,ichem

      open(ifileno, file=adjustl(kfile))

      i_err=0
      do while(i_err==0)
         read (ifileno,*,iostat=i_err) iat,ichem,anisotropytype(iat,ichem),anisotropy(iat,1,ichem),anisotropy(iat,2,ichem),&
            anisotropy(iat,3,ichem),anisotropy(iat,4,ichem),anisotropy(iat,5,ichem),anisotropy(iat,6,ichem)
      enddo
      close (ifileno)
   end subroutine read_anisotropy_alloy


   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_anisotropy
   !> @brief Read the anisotropy variables for pure systems
   !---------------------------------------------------------------------------------
   subroutine read_anisotropy()
      !
      implicit none
      !
      integer :: m,iat

      open(ifileno, file=adjustl(kfile))

      do m=1, na
         read (ifileno,*) iat,anisotropytype(iat,1),anisotropy(iat,1,1),anisotropy(iat,2,1),anisotropy(iat,3,1),&
            anisotropy(iat,4,1),anisotropy(iat,5,1),anisotropy(iat,6,1)
      enddo

      if (mult_axis=='Y') then
         do m=1,na
            read(ifileno,*) iat,anisotropytype_diff(iat,1),anisotropy_diff(iat,1,1),anisotropy_diff(iat,2,1),anisotropy_diff(iat,3,1),&
               anisotropy_diff(iat,4,1),anisotropy_diff(iat,5,1),anisotropy_diff(iat,6,1)
         enddo

      endif
      close (ifileno)

   end subroutine read_anisotropy

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_dmdata
   !> @brief Read the variables for the DMI
   !---------------------------------------------------------------------------------
   subroutine read_dmdata()
      !
      implicit none
      !
      integer :: flines, isite, i_stat, jsite
      integer :: flines_clus, isite_c, jsite_c
      integer :: itype, jtype, ichem, jchem, iline, ishell
      integer :: itype_clus, jtype_clus, ichem_c, jchem_c
      logical :: unique
      real(dblprec), dimension(3) :: dm_tmp,dm_tmp_clus
      real(dblprec):: tol, norm
      real(dblprec), dimension(3) :: r_tmp, r_red, r_tmp_clus,r_red_clus
      integer :: no_shells, no_shells_clus

      integer :: ifileno2

      ! Set tolerance for neighbour shells
      tol=1.0d-5

      ifileno2=ifileno+1
      open(ifileno, file=dmfile)
      if (do_cluster=='Y') open(ifileno2, file=dmfile_clus) ! File for the interactions inside the embeded cluster

      if (do_cluster=='Y') then
         call read_exchange_getMaxNoShells(no_shells,flines,dmfile)
         max_no_dmshells = no_shells
         allocate(dm_redcoord(NT,max_no_dmshells,3),stat=i_stat)
         call memocc(i_stat,product(shape(dm_redcoord))*kind(dm_redcoord),'dm_redcoord','read_dmdata')
         dm_redcoord  = 0.0_dblprec
         allocate(dm_inpvect(3,NT,max_no_dmshells,Nchmax,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(dm_inpvect))*kind(dm_inpvect),'dm_inpvect','read_dmdata')
         dm_inpvect = 0.0_dblprec
         call read_exchange_getMaxNoShells_clus(no_shells_clus,flines_clus,ifileno2,do_ralloy)
         max_no_dmshells_clus=no_shells_clus
         allocate(dm_redcoord_clus(NT_clus,max_no_dmshells_clus,3),stat=i_stat)
         call memocc(i_stat,product(shape(dm_redcoord_clus))*kind(dm_redcoord_clus),'dm_redcoord_clus','read_dmdata')
         dm_redcoord_clus  = 0.0_dblprec
         allocate(dm_inpvect_clus(3,NT_clus,max_no_dmshells_clus,Nchmax,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(dm_inpvect_clus))*kind(dm_inpvect_clus),'dm_inpvect_clus','read_dmdata')
         dm_inpvect_clus = 0.0_dblprec
      else
         call read_exchange_getMaxNoShells(no_shells,flines,dmfile)
         max_no_dmshells = no_shells
         allocate(dm_redcoord(NT,max_no_dmshells,3),stat=i_stat)
         call memocc(i_stat,product(shape(dm_redcoord))*kind(dm_redcoord),'dm_redcoord','read_dmdata')
         dm_redcoord  = 0.0_dblprec
         allocate(dm_inpvect(3,NT,max_no_dmshells,Nchmax,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(dm_inpvect))*kind(dm_inpvect),'dm_inpvect','read_dmdata')
         dm_inpvect = 0.0_dblprec
      endif

      ! Read exchange vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp, dm_tmp
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, dm_tmp
         end if

         ! Find type of site
         itype=atype_inp(isite)
         jtype=1

         call read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,dm_nn(itype)
            norm=(r_red(1)-dm_redcoord(itype,ishell,1))**2+ &
               (r_red(2)-dm_redcoord(itype,ishell,2))**2+ &
               (r_red(3)-dm_redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  dm_inpvect(1:3,itype,ishell,ichem,jtype)=dm_tmp
               else
                  dm_inpvect(1:3,itype,ishell,ichem,jchem)=dm_tmp
               end if
            end if
         end do
         if (unique) then
            dm_nn(itype)=dm_nn(itype)+1

            dm_redcoord(itype,dm_nn(itype),1:3)=r_red(1:3)
            if(do_ralloy==0) then
               dm_inpvect(1:3,itype,dm_nn(itype),ichem,1)=dm_tmp
            else
               dm_inpvect(1:3,itype,dm_nn(itype),ichem,jchem)=dm_tmp
            end if
         end if
      enddo

      if (do_cluster=='Y') then

         do iline=1, flines_clus

            ! Read indices and coordinates
            if(do_ralloy==0) then
               read (ifileno2,*) isite_c, jsite_C, r_tmp_clus, dm_tmp_clus
               ichem_c=1
               jchem_c=1
            else
               read (ifileno2,*) isite_c, jsite_c, ichem_c, jchem_c, r_tmp_clus, dm_tmp_clus
            end if

            ! Find type of site
            itype_clus=atype_inp_clus(isite_c)
            jtype_clus=1

            call read_exchange_getNeighVec_clus(r_red_clus,r_tmp_clus,isite_c,jsite_c,&
               maptype,posfiletype)

            ! Loop through earlier vectors to find equivalent shells
            unique=.true.
            do ishell=1,dm_nn_clus(itype_clus)
               norm=(r_red_clus(1)-dm_redcoord_clus(itype_clus,ishell,1))**2+ &
               (r_red_clus(2)-dm_redcoord_clus(itype_clus,ishell,2))**2+ &
               (r_red_clus(3)-dm_redcoord_clus(itype_clus,ishell,3))**2
               if(norm<tol) then
                  unique=.false.
                  if(do_ralloy==0) then
                     dm_inpvect_clus(1:3,itype_clus,ishell,ichem_c,jtype_clus)=dm_tmp_clus
                  else
                     dm_inpvect_clus(1:3,itype_clus,ishell,ichem_c,jchem_c)=dm_tmp_clus
                  end if
               end if
            end do
            if (unique) then
               dm_nn_clus(itype_clus)=dm_nn_clus(itype_clus)+1

               dm_redcoord_clus(itype_clus,dm_nn_clus(itype_clus),1:3)=r_red_clus(1:3)
               if(do_ralloy==0) then
                  dm_inpvect_clus(1:3,itype_clus,dm_nn_clus(itype_clus),ichem_c,1)=dm_tmp_clus
               else
                  dm_inpvect_clus(1:3,itype_clus,dm_nn_clus(itype_clus),ichem_c,jchem_c)=dm_tmp_clus
               end if
            end if
         enddo
      endif
      max_no_dmshells=maxval(dm_nn)
      call read_exchange_reduceRedCoordMatrixSize(dm_redcoord,nt,max_no_dmshells)
      call read_dmexchange_reduceCouplingMatrixSize(dm_inpvect,nt,max_no_dmshells,nchmax)

      if (do_cluster=='Y') then
         max_no_dmshells_clus=maxval(dm_nn_clus)
         call read_exchange_reduceRedCoordMatrixSize(dm_redcoord_clus,NT_clus,max_no_dmshells_clus)
         call read_dmexchange_reduceCouplingMatrixSize(dm_inpvect_clus,NT_clus,max_no_dmshells_clus,nchmax)
         close(ifileno2)
      endif
      close (ifileno)

   end subroutine read_dmdata


   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_chirdata
   !> @brief Read the variables for three-site scalar chirality interaction
   !---------------------------------------------------------------------------------
   !      ! Variables for CHIR exchange
   ! integer :: nn_chir_tot                                     !< Calculated number of neighbours with PD interactions
   ! integer ::  max_no_chirneigh                               !< Calculated maximum of neighbours for PD exchange
   ! integer, dimension(:), allocatable :: chirlistsize         !< Size of neighbour list for PD
   ! integer, dimension(:,:), allocatable :: chirlist           !< List of neighbours for PD
   ! real(dblprec), dimension(:,:,:), allocatable :: chir_vect  !< Pseudo-Dipolar exchange vector

   subroutine read_chirdata()
      !
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec) :: chir_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(6) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp
      integer :: i, nskip

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=chirfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      call read_exchange_getMaxNoShells(no_shells,flines,chirfile)
      !max_no_chirshells = no_shells
      !call allocate_latthamiltonianinput(no_shells,1)

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_chirdata')
      nn_tmp=0

      ! Pre-read file to get max no. exchange shells and no. lines
      do
         if(do_ralloy==0) then
            read(ifileno,*,end=200)  isite,jsite
            ichem=1
            jchem=1
         else
            read(ifileno,*,end=200)  isite,jsite, ichem, jchem
         end if
         ! Skip 0, 3 or 9 lines
         do i=1,nskip
            read(ifileno,*)
         end do
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)
         if(do_ralloy==0) then
            nn_tmp(itype,jtype,jchem)=nn_tmp(itype,jtype,jchem)+1
         else
            nn_tmp(itype,ichem,jchem)=nn_tmp(itype,ichem,jchem)+1
         end if
      end do
      200 continue

      rewind(ifileno)

      no_shells=0
      do itype=1,nt
         no_shells=max(sum(nn_tmp(itype,:,:)),no_shells)
      end do

      i_all=-product(shape(nn_tmp))*kind(nn_tmp)
      deallocate(nn_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'NN_tmp','read_chirdata')

      max_no_chirshells=no_shells
      allocate(chir_redcoord(NT,max_no_chirshells,3,2),stat=i_stat)
      call memocc(i_stat,product(shape(chir_redcoord))*kind(chir_redcoord),'chir_redcoord','read_chirdata')
      chir_redcoord=0.0_dblprec
      allocate(chir_inpval(NT,max_no_chirshells,NT,NT),stat=i_stat)
      call memocc(i_stat,product(shape(chir_inpval))*kind(chir_inpval),'chir_inpval','read_chirdata')
      chir_inpval=0.0_dblprec
      ! Read force coupling vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         ! A block of four lines is used for each input chir-tensor
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp, chir_tmp
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, chir_tmp
         end if

         ! Find type of site
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)

         if(maptype==2) then
            ! Calculate proper neighbour vector (from "bgfm")
            r_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(1)+C2(1)*r_tmp(2)+C3(1)*r_tmp(3)
            r_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(1)+C2(2)*r_tmp(2)+C3(2)*r_tmp(3)
            r_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(1)+C2(3)*r_tmp(2)+C3(3)*r_tmp(3)
            r_red(4)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(4)+C2(1)*r_tmp(5)+C3(1)*r_tmp(6)
            r_red(5)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(4)+C2(2)*r_tmp(5)+C3(2)*r_tmp(6)
            r_red(6)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(4)+C2(3)*r_tmp(5)+C3(3)*r_tmp(6)
         else
            ! Calculates neighbour vectors from direct coordinates or Cartesian
            ! coordinates, corresponding to how the atomic positions are entered
            !if (posfiletype=='C') then
               r_red=r_tmp
            !elseif (posfiletype=='D') then
            !   r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
            !   r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
            !   r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
            !   r_red(4)=r_tmp(4)*C1(1)+r_tmp(5)*C2(1)+r_tmp(6)*C3(1)
            !   r_red(5)=r_tmp(4)*C1(2)+r_tmp(5)*C2(2)+r_tmp(6)*C3(2)
            !   r_red(6)=r_tmp(4)*C1(3)+r_tmp(5)*C2(3)+r_tmp(6)*C3(3)
            !else
            !   stop 'Only posfiletype= C or D is currently supported'
            !endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,chir_nn(itype)
            norm=(r_red(1)-chir_redcoord(itype,ishell,1,1))**2+ &
               (r_red(2)-chir_redcoord(itype,ishell,2,1))**2+ &
               (r_red(3)-chir_redcoord(itype,ishell,3,1))**2+ &
               (r_red(4)-chir_redcoord(itype,ishell,1,2))**2+ &
               (r_red(5)-chir_redcoord(itype,ishell,2,2))**2+ &
               (r_red(6)-chir_redcoord(itype,ishell,3,2))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  chir_inpval(itype,ishell,ichem,jtype)=chir_tmp
               else
                  chir_inpval(itype,ishell,ichem,jchem)=chir_tmp
               end if
            end if
         end do
         if (unique) then
            chir_nn(itype)=chir_nn(itype)+1
            chir_redcoord(itype,chir_nn(itype),1:3,1)=r_red(1:3)
            chir_redcoord(itype,chir_nn(itype),1:3,2)=r_red(4:6)
            if(do_ralloy==0) then
               chir_inpval(itype,chir_nn(itype),ichem,jtype)=chir_tmp
            else
               chir_inpval(itype,chir_nn(itype),ichem,jchem)=chir_tmp
            end if
         end if
      enddo
      close (ifileno)


   end subroutine read_chirdata

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_pddata
   !> @brief Read the variables for the pseudo-dipolar interaction
   !---------------------------------------------------------------------------------
   subroutine read_pddata()
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec), dimension(6) :: pd_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(3) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=pdfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_pddata')
      nn_tmp=0

      ! Pre-read file to get max no. exchange shells and no. lines
      do
         if(do_ralloy==0) then
            read(ifileno,*,end=200)  isite,jsite
            ichem=1
            jchem=1
         else
            read(ifileno,*,end=200)  isite,jsite, ichem, jchem
         end if
         flines=flines+1
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)
         if(do_ralloy==0) then
            nn_tmp(itype,jtype,jchem)=nn_tmp(itype,jtype,jchem)+1
         else
            nn_tmp(itype,ichem,jchem)=nn_tmp(itype,ichem,jchem)+1
         end if
      end do
      200 continue

      rewind(ifileno)

      no_shells=0
      do itype=1,nt
         no_shells=max(sum(nn_tmp(itype,:,:)),no_shells)
      end do

      i_all=-product(shape(nn_tmp))*kind(nn_tmp)
      deallocate(nn_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'NN_tmp','read_pddata')

      max_no_pdshells=no_shells
      allocate(pd_redcoord(NT,max_no_pdshells,3),stat=i_stat)
      call memocc(i_stat,product(shape(pd_redcoord))*kind(pd_redcoord),'pd_redcoord','read_pddata')
      pd_redcoord=0.0_dblprec
      allocate(pd_inpvect(6,NT,max_no_pdshells,Nchmax,Nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(pd_inpvect))*kind(pd_inpvect),'pd_inpvect','read_pddata')
      pd_inpvect=0.0_dblprec
      ! Read exchange vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp, pd_tmp
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, pd_tmp
         end if

         ! Find type of site
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)

         if(maptype==2) then
            ! Calculate proper neighbour vector (from "bgfm")
            r_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(1)+C2(1)*r_tmp(2)+C3(1)*r_tmp(3)
            r_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(1)+C2(2)*r_tmp(2)+C3(2)*r_tmp(3)
            r_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(1)+C2(3)*r_tmp(2)+C3(3)*r_tmp(3)
         else
            ! Calculates neighbour vectors from direct coordinates or Cartesian
            ! coordinates, corresponding to how the atomic positions are entered
            !if (posfiletype=='C') then
               r_red=r_tmp
            !elseif (posfiletype=='D') then
            !   r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
            !   r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
            !   r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
            !else
            !   stop 'Only posfiletype= C or D is currently supported'
            !endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,pd_nn(itype)
            norm=(r_red(1)-pd_redcoord(itype,ishell,1))**2+ &
               (r_red(2)-pd_redcoord(itype,ishell,2))**2+ &
               (r_red(3)-pd_redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  pd_inpvect(1:6,itype,ishell,ichem,jtype)=pd_tmp
               else
                  pd_inpvect(1:6,itype,ishell,ichem,jchem)=pd_tmp
               end if
            end if
         end do
         if (unique) then
            pd_nn(itype)=pd_nn(itype)+1
            pd_redcoord(itype,pd_nn(itype),1:3)=r_red(1:3)
            if(do_ralloy==0) then
               pd_inpvect(1:6,itype,pd_nn(itype),ichem,jtype)=pd_tmp
            else
               pd_inpvect(1:6,itype,pd_nn(itype),ichem,jchem)=pd_tmp
            end if
         end if
      enddo
      close (ifileno)
   end subroutine read_pddata

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_biqdmdata
   !> @brief Read the variables for the biquadratic-DM interaction
   !---------------------------------------------------------------------------------
   subroutine read_biqdmdata()
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique

      real(dblprec) :: biqdm_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(3) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=biqdmfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_biqdmdata')
      nn_tmp=0

      ! Pre-read file to get max no. exchange shells and no. lines
      do
         if(do_ralloy==0) then
            read(ifileno,*,end=200)  isite,jsite
            ichem=1
            jchem=1
         else
            read(ifileno,*,end=200)  isite,jsite, ichem, jchem
         end if
         flines=flines+1
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)
         if(do_ralloy==0) then
            nn_tmp(itype,jtype,jchem)=nn_tmp(itype,jtype,jchem)+1
         else
            nn_tmp(itype,ichem,jchem)=nn_tmp(itype,ichem,jchem)+1
         end if
      end do
      200 continue

      rewind(ifileno)

      no_shells=0
      do itype=1,nt
         no_shells=max(sum(nn_tmp(itype,:,:)),no_shells)
      end do

      i_all=-product(shape(nn_tmp))*kind(nn_tmp)
      deallocate(nn_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'NN_tmp','read_biqdmdata')

      max_no_biqdmshells=no_shells
      write(*,*) 'max_no_biqdmshells ', max_no_biqdmshells
      allocate(biqdm_redcoord(NT,max_no_biqdmshells,3),stat=i_stat)
      call memocc(i_stat,product(shape(redcoord))*kind(redcoord),'biqdm_redcoord','read_biqdmdata')
      biqdm_redcoord=0.0_dblprec
      allocate(biqdm_inpvect(1,NT,max_no_biqdmshells,Nchmax,Nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(redcoord))*kind(redcoord),'biqdm_inpvect','read_biqdmdata')
      biqdm_inpvect=0.0_dblprec
      ! Read exchange vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines
         ! Read indices and coordinates
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp, biqdm_tmp
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, biqdm_tmp
         end if

         ! Find type of site
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)

         if(maptype==2) then
            ! Calculate proper neighbour vector (from "bgfm")
            r_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(1)+C2(1)*r_tmp(2)+C3(1)*r_tmp(3)
            r_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(1)+C2(2)*r_tmp(2)+C3(2)*r_tmp(3)
            r_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(1)+C2(3)*r_tmp(2)+C3(3)*r_tmp(3)
         else
            ! Calculates neighbour vectors from direct coordinates or Cartesian
            ! coordinates, corresponding to how the atomic positions are entered
            !if (posfiletype=='C') then
               r_red=r_tmp
            !elseif (posfiletype=='D') then
            !   r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
            !   r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
            !   r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
            !else
            !   stop 'Only posfiletype= C or D is currently supported'
            !endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,biqdm_nn(itype)
            norm=(r_red(1)-biqdm_redcoord(itype,ishell,1))**2+ &
               (r_red(2)-biqdm_redcoord(itype,ishell,2))**2+ &
               (r_red(3)-biqdm_redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  biqdm_inpvect(1,itype,ishell,ichem,jtype)=biqdm_tmp
               else
                  biqdm_inpvect(1,itype,ishell,ichem,jchem)=biqdm_tmp
               end if
            end if
         end do
         if (unique) then
            biqdm_nn(itype)=biqdm_nn(itype)+1
            biqdm_redcoord(itype,biqdm_nn(itype),1:3)=r_red(1:3)
            if(do_ralloy==0) then
               biqdm_inpvect(1,itype,biqdm_nn(itype),ichem,jtype)=biqdm_tmp
            else
               biqdm_inpvect(1,itype,biqdm_nn(itype),ichem,jchem)=biqdm_tmp
            end if
         end if
      enddo

      close (ifileno)

   end subroutine read_biqdmdata

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_bqdata
   !> @brief Read the variables for the biquadratic interaction
   !---------------------------------------------------------------------------------
   subroutine read_bqdata()
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec):: jbq_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(3) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=bqfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_bqdata')
      nn_tmp=0

      ! Pre-read file to get max no. exchange shells and no. lines
      do
         if(do_ralloy==0) then
            read(ifileno,*,end=200)  isite,jsite
            ichem=1
            jchem=1
         else
            read(ifileno,*,end=200)  isite,jsite, ichem, jchem
         end if
         flines=flines+1
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)
         if(do_ralloy==0) then
            nn_tmp(itype,jtype,jchem)=nn_tmp(itype,jtype,jchem)+1
         else
            nn_tmp(itype,ichem,jchem)=nn_tmp(itype,ichem,jchem)+1
         end if
      end do
      200 continue

      rewind(ifileno)

      no_shells=0
      do itype=1,nt
         no_shells=max(sum(nn_tmp(itype,:,:)),no_shells)
      end do

      i_all=-product(shape(nn_tmp))*kind(nn_tmp)
      deallocate(nn_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'NN_tmp','read_bqdata')

      max_no_bqshells=no_shells
      allocate(bq_redcoord(NT,max_no_bqshells,3),stat=i_stat)
      call memocc(i_stat,product(shape(redcoord))*kind(redcoord),'bq_redcoord','read_bqdata')
      bq_redcoord=0.0_dblprec
      allocate(jc_bq(NT,max_no_bqshells,Nchmax,Nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(redcoord))*kind(redcoord),'bq_inpvect','read_bqdata')
      jc_bq=0.0_dblprec
      ! Read exchange vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp, jbq_tmp
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, jbq_tmp
         end if

         itype=atype_inp(isite)
         jtype=atype_inp(jsite)
         if(maptype==2) then
            ! Calculate proper neighbour vector (from "bgfm")
            r_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(1)+C2(1)*r_tmp(2)+C3(1)*r_tmp(3)
            r_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(1)+C2(2)*r_tmp(2)+C3(2)*r_tmp(3)
            r_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(1)+C2(3)*r_tmp(2)+C3(3)*r_tmp(3)
         else
            ! Calculates neighbour vectors from direct coordinates or Cartesian
            ! coordinates, corresponding to how the atomic positions are entered
            !if (posfiletype=='C') then
               r_red=r_tmp
            !elseif (posfiletype=='D') then
            !   r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
            !   r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
            !   r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
            !else
            !   stop 'Only posfiletype= C or D is currently supported'
            !endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,bq_nn(itype)
            norm=(r_red(1)-bq_redcoord(itype,ishell,1))**2+ &
               (r_red(2)-bq_redcoord(itype,ishell,2))**2+ &
               (r_red(3)-bq_redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  jc_bq(itype,ishell,ichem,jtype)=jbq_tmp
               else
                  jc_bq(itype,ishell,ichem,jchem)=jbq_tmp
               end if
            end if
         end do
         if (unique) then
            bq_nn(itype)=bq_nn(itype)+1
            bq_redcoord(itype,bq_nn(itype),1:3)=r_red(1:3)
            if(do_ralloy==0) then
               jc_bq(itype,bq_nn(itype),ichem,jtype)=jbq_tmp
            else
               jc_bq(itype,bq_nn(itype),ichem,jchem)=jbq_tmp
            end if
         end if
      enddo
      close (ifileno)

   end subroutine read_bqdata

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_sitefield
   !> @brief Read the site dependent magnetic field
   !> @note One needs option do_bpulse 6 for this to work.
   !---------------------------------------------------------------------------------
   subroutine read_sitefield(Natom,sitenatomfld)

      implicit none

      integer, intent(in) :: Natom
      real(dblprec), dimension(:,:), allocatable, intent(out) :: sitenatomfld

      integer :: i,flines, isite, i_stat

      open(ifileno, file=trim(siteatomfile))

      flines=0
      ! Pre-read file to get number of lines
      do
         read(ifileno,*,end=200)  isite
         flines=flines+1
      end do

      200 continue

      rewind(ifileno)

      write(*,'(2x,a)') 'Reading site dependent fields'

      ! Allocate the site-dependent field
      allocate(sitenatomfld(3,flines),stat=i_stat)
      call memocc(i_stat,product(shape(sitenatomfld))*kind(sitenatomfld),'sitenatomfld','read_sitefield')
      sitenatomfld=0.0_dblprec
      ! If the size of the file is NATOM then there is no problem
      if ( Natom.eq.flines ) then

         do i=1, flines
            read(ifileno,*) isite, sitenatomfld(1,isite), sitenatomfld(2,isite), sitenatomfld(3,isite)
         end do
      else
         write(*,*) 'WARNING: Size of the SITEATOMFLD is not NATOM'
         do i=1, flines
            read(ifileno,*) isite, sitenatomfld(1,isite), sitenatomfld(2,isite), sitenatomfld(3,isite)
         end do

      end if

      close(ifileno)

   end subroutine read_sitefield


   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_ip_damping
   !> @brief Read the variables for the initial phase site dependent damping parameter
   !---------------------------------------------------------------------------------
   subroutine read_ip_damping()
      ! Read damping parameter from file
      !
      implicit none
      !
      integer :: i,j,i_stat,flines
      real(dblprec) :: idamp1, idamp2

      ! Open input file
      open(ifileno, file=ip_dampfile)
      ! Check if input file is for random alloy
      flines=0
      ! Pre-read file to get the number of lines
      do
         read(ifileno,*,end=200)
         flines=flines+1
      end do
      200 continue
      rewind(ifileno)

      ! In case of error genererating the input file
      if(ipnphase*NA/=flines) write(*,*) "WARNING: Check input file ", ip_dampfile, &
         " for inconsistent information."

      ! Allocate input arrays
      allocate(ipdamping1(ipnphase,NA),stat=i_stat)
      call memocc(i_stat,product(shape(ipdamping1))*kind(ipdamping1),'ipdamping1','read_ip_damping')
      ipdamping1=0.0_dblprec
      allocate(ipdamping2(ipnphase,NA),stat=i_stat)
      call memocc(i_stat,product(shape(ipdamping2))*kind(ipdamping2),'ipdamping2','read_ip_damping')
      ipdamping2=0.0_dblprec
      ! Read Site, Type, damping parameter per atom
      do i=1, ipnphase
         do j=1,NA
            read (ifileno,*) idamp1, idamp2
            ipdamping1(i,j)=idamp1
            ipdamping2(i,j)=idamp2
         end do
      end do
      close (ifileno)
   end subroutine read_ip_damping

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_ip_damping_alloy
   !> @brief Read the variables for the initial phase site dependent damping parameter for a random alloy
   !---------------------------------------------------------------------------------
   subroutine read_ip_damping_alloy()
      ! Read damping parameter from file
      !
      implicit none
      !
      integer :: i,j,i_stat,flines,ipNlines
      real(dblprec) :: idamp1, idamp2

      ! Open input file
      open(ifileno, file=ip_dampfile)

      flines=0
      ! Pre-read file to get the number of lines
      do
         read(ifileno,*,end=200)
         flines=flines+1
      end do
      200 continue
      rewind(ifileno)

      ipNlines = flines / ipnphase

      ! In case of error genererating the input file
      if((ipNlines*ipnphase)/=flines) write(*,*) "WARNING: Check input file ", ip_dampfile, &
         " for inconsistent information."

      ! Allocate input arrays
      allocate(ipdampingalloy1(ipnphase,ipNlines,ipNlines),stat=i_stat)
      call memocc(i_stat,product(shape(ipdampingalloy1))*kind(ipdampingalloy1),'ipdampingalloy1','read_ip_damping_alloy')
      ipdampingalloy1=0.0_dblprec
      allocate(ipdampingalloy2(ipnphase,ipNlines,ipNlines),stat=i_stat)
      call memocc(i_stat,product(shape(ipdampingalloy2))*kind(ipdampingalloy2),'ipdampingalloy2','read_ip_damping_alloy')
      ipdampingalloy2=0.0_dblprec
      ! Read Site, Type, damping parameter per atom
      do i=1, ipnphase
         do j=1, ipNlines
            read (ifileno,*) idamp1, idamp2
            ipdampingalloy1(i,asite(j),acomp(j))=idamp1
            ipdampingalloy2(i,asite(j),acomp(j))=idamp2
         enddo
      enddo

      close (ifileno)
   end subroutine read_ip_damping_alloy

   !---------------------------------------------------------------------------------
   ! subroutine: read_damping
   !> @brief read measurement phase site dependent damping
   !---------------------------------------------------------------------------------
   subroutine read_damping()
      ! Read damping parameter from file
      !
      implicit none
      !
      integer :: i,i_stat,flines
      real(dblprec) :: idamp1, idamp2

      ! Open input file
      open(ifileno, file=mp_dampfile)

      flines=0
      ! Pre-read file to get the number of lines
      do
         read(ifileno,*,end=200)
         flines=flines+1
      end do
      200 continue
      rewind(ifileno)

      ! In case of error genererating the input file
      if(NA/=flines) write(*,*) "WARNING: Check input file ", mp_dampfile, &
         " for inconsistent information."

      ! Allocate input arrays
      allocate(mpdamping1(NA),stat=i_stat)
      call memocc(i_stat,product(shape(mpdamping1))*kind(mpdamping1),'mpdamping1','read_damping')
      mpdamping1=0.0_dblprec
      allocate(mpdamping2(NA),stat=i_stat)
      call memocc(i_stat,product(shape(mpdamping2))*kind(mpdamping2),'mpdamping2','read_damping')
      mpdamping2=0.0_dblprec
      write(*,*) "Check", NA
      ! Read Site, Type, damping parameter per atom
      do i=1, NA
         read (ifileno,*) idamp1, idamp2
         mpdamping1(i)=idamp1
         mpdamping2(i)=idamp2
      enddo
      close (ifileno)
   end subroutine read_damping

   !---------------------------------------------------------------------------------
   ! subroutine: read_damping_alloy
   !> @brief read measurement phase site dependent damping for a a random alloy
   !---------------------------------------------------------------------------------
   subroutine read_damping_alloy()
      ! Read damping parameter from a file for a random alloy
      !
      implicit none

      !
      integer :: i_stat, i, flines
      real(dblprec) :: idamp1, idamp2

      ! Open input file
      open(ifileno, file=mp_dampfile)

      flines=0
      ! Pre-read file to get the number of lines
      do
         read(ifileno,*,end=200)
         flines=flines+1
      end do
      200 continue
      rewind(ifileno)

      mpNlines = flines

      ! Allocate input arrays
      allocate(mpdampingalloy1(flines,flines),stat=i_stat)
      call memocc(i_stat,product(shape(mpdampingalloy1))*kind(mpdampingalloy1),'mpdampingalloy1','read_damping_alloy')
      mpdampingalloy1=0.0_dblprec
      allocate(mpdampingalloy2(flines,flines),stat=i_stat)
      call memocc(i_stat,product(shape(mpdampingalloy2))*kind(mpdampingalloy2),'mpdampingalloy2','read_damping_alloy')
      mpdampingalloy2=0.0_dblprec

      ! Read data
      ! Site,  Type, Damping parameter
      do i=1, flines
         read (ifileno,*) idamp1,idamp2
         mpdampingalloy1(asite(i),acomp(i))=idamp1
         mpdampingalloy2(asite(i),acomp(i))=idamp2
      end do
      close (ifileno)
   end subroutine read_damping_alloy

   !--------------------------------------------------------------------------------
   !> @brief
   !> Modified version of the read exchange routine by Anders Bergman so that
   !> it reads the energy barriers used for the KMC
   !> @author
   !> Jonathan Chico
   !>
   !> @date 03/02/2017 - Jonathan Chico
   !--------------------------------------------------------------------------------
   subroutine read_barriers()
      !
      implicit none

      integer :: itype,jtype,isite,jsite,ichem,jchem,iline,ishell
      integer :: flines,no_shells
      real(dblprec), dimension(3) :: r_red, r_tmp
      logical :: unique
      real(dblprec):: barr_tmp
      real(dblprec):: tol, norm

      ! Set tolerance for neighbour shells
      tol=1.0d-5

      open(ifileno, file=trim(barrfile)) ! Number of shells same for different LSF configuration

      call read_exchange_getMaxNoShells(no_shells,flines,barrfile)
      call allocate_barriers(no_shells,1)

      redcoord_barriers = 0.0_dblprec
      nn_barriers       = 0

      ! Read exchange vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines
         ! Loop through earlier vectors to find equivalent shells

         ! Read indices and coordinates
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp(1:3), barr_tmp
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp(1:3), barr_tmp
         end if

         itype=atype_inp(isite)
         jtype=1
         call read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.

         do ishell=1,nn_barriers(itype)

            norm=(r_red(1)-redcoord_barriers(itype,ishell,1))**2+ &
               (r_red(2)-redcoord_barriers(itype,ishell,2))**2+ &
               (r_red(3)-redcoord_barriers(itype,ishell,3))**2

            if(norm<tol) then
               unique=.false.
               ! If neighbour vector already exist, replace energy barrier value (could be removed)
               if(do_ralloy==0) then
                  kmc_barriers(itype,ishell,ichem,jtype,1)=barr_tmp
               else
                  kmc_barriers(itype,ishell,ichem,jchem,1)=barr_tmp
               end if
            end if
         end do

         if (unique) then
            ! Add entry if not found earlier
            nn_barriers(itype)=nn_barriers(itype)+1
            redcoord_barriers(itype,nn_barriers(itype),1:3)=r_red(1:3)

            if(do_ralloy==0) then
               kmc_barriers(itype,nn_barriers(itype),ichem,jtype,1)=barr_tmp
            else
               kmc_barriers(itype,nn_barriers(itype),ichem,jchem,1)=barr_tmp
            end if
         end if
      enddo

      ! Reducing jc size if max_no_shells are small enough !
      max_no_shells_barriers=maxval(NN_barriers)
      close(ifileno)

      call read_exchange_reduceRedCoordMatrixSize(redcoord_barriers,nt,max_no_shells_barriers)
      call read_exchange_reduceCouplingMatrixSize(kmc_barriers,nt,max_no_shells_barriers,nchmax)

   end subroutine read_barriers


   !---------------------------------------------------------------------------------
   ! subroutine: allocate_hamiltonianinput
   !> @brief Allocate arrays for input for Hamiltonian
   !---------------------------------------------------------------------------------
   subroutine allocate_hamiltonianinput(no_shells, flag) !NA, limit_no_shells, Nchmax, flag)
      use Parameters
      use Profiling
      use InputData
      implicit none

      integer, intent(in),optional :: no_shells !< Parameter limiting number of exchange coupling shells
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then

         if(do_jtensor==0) then
            allocate(jc(NT,no_shells,Nchmax,Nchmax,conf_num),stat=i_stat)
            call memocc(i_stat,product(shape(jc))*kind(jc),'jc','allocate_hamiltonianinput')
            jc=0.0_dblprec
            if(exc_inter=='Y') then
               allocate(jcD(NT,no_shells,Nchmax,Nchmax,conf_num),stat=i_stat)
               call memocc(i_stat,product(shape(jcD))*kind(jcD),'jcD','allocate_hamiltonianinput')
               jcD=0.0_dblprec
            endif
         else
            allocate(jc_tens(3,3,NT,no_shells,Nchmax,Nchmax),stat=i_stat)
            call memocc(i_stat,product(shape(jc_tens))*kind(jc_tens),'jc_tens','allocate_hamiltonianinput')
            jc_tens=0.0_dblprec
         end if
         allocate(redcoord(NT,no_shells,3),stat=i_stat)
         call memocc(i_stat,product(shape(redcoord))*kind(redcoord),'redcoord','allocate_hamiltonianinput')
         redcoord=0.0_dblprec
         allocate(NN(NT),stat=i_stat)
         call memocc(i_stat,product(shape(NN))*kind(NN),'NN','allocate_hamiltonianinput')
         NN=0
         allocate(dm_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(dm_nn))*kind(dm_nn),'dm_nn','allocate_hamiltonianinput')
         dm_nn=0
         allocate(pd_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(pd_nn))*kind(pd_nn),'pd_nn','allocate_hamiltonianinput')
         pd_nn=0
         allocate(chir_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(chir_nn))*kind(chir_nn),'chir_nn','allocate_hamiltonianinput')
         chir_nn=0
         allocate(biqdm_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(biqdm_nn))*kind(biqdm_nn),'biqdm_nn','allocate_hamiltonianinput')
         biqdm_nn=0
         allocate(bq_nn(NT),stat=i_stat)
         call memocc(i_stat,product(shape(bq_nn))*kind(bq_nn),'bq_nn','allocate_hamiltonianinput')
         bq_nn=0
         allocate(anisotropytype(NA,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(anisotropytype))*kind(anisotropytype),'anisotropytype','allocate_hamiltonianinput')
         anisotropytype=0
         allocate(anisotropy(NA,6,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(anisotropy))*kind(anisotropy),'anisotropy','allocate_hamiltonianinput')
         anisotropy=0.0_dblprec

         if (mult_axis=='Y') then
            allocate(anisotropytype_diff(NA,Nchmax),stat=i_stat)
            call memocc(i_stat,product(shape(anisotropytype_diff))*kind(anisotropytype_diff),'anisotropytype_diff','allocate_hamiltonianinput')
            allocate(anisotropy_diff(NA,6,Nchmax),stat=i_stat)
            call memocc(i_stat,product(shape(anisotropy_diff))*kind(anisotropy_diff),'anisotropy_diff','allocate_hamiltonianinput')
            anisotropy_diff=0.0_dblprec
         endif

      else

         if(do_jtensor/=1) then
            i_all=-product(shape(jc))*kind(jc)
            deallocate(jc,stat=i_stat)
            call memocc(i_stat,i_all,'jc','allocate_hamiltonianinput')
            if(exc_inter=='Y') then
               deallocate(jcD,stat=i_stat)
               call memocc(i_stat,i_all,'jcD','allocate_hamiltonianinput')
            endif
         else
            i_all=-product(shape(jc_tens))*kind(jc_tens)
            deallocate(jc_tens,stat=i_stat)
            call memocc(i_stat,i_all,'jc_tens','allocate_hamiltonianinput')
         end if
         i_all=-product(shape(dm_nn))*kind(dm_nn)
         deallocate(dm_nn,stat=i_stat)
         call memocc(i_stat,i_all,'dm_nn','allocate_hamiltonianinput')
         i_all=-product(shape(pd_nn))*kind(pd_nn)
         deallocate(pd_nn,stat=i_stat)
         call memocc(i_stat,i_all,'pd_nn','allocate_hamiltonianinput')
         i_all=-product(shape(biqdm_nn))*kind(biqdm_nn)
         deallocate(biqdm_nn,stat=i_stat)
         call memocc(i_stat,i_all,'biqdm_nn','allocate_hamiltonianinput')
         i_all=-product(shape(bq_nn))*kind(bq_nn)
         deallocate(bq_nn,stat=i_stat)
         call memocc(i_stat,i_all,'bq_nn','allocate_hamiltonianinput')
         i_all=-product(shape(NN))*kind(NN)
         deallocate(NN,stat=i_stat)
         call memocc(i_stat,i_all,'NN','allocate_hamiltonianinput')
         i_all=-product(shape(jfile))*kind(jfile)
         deallocate(jfile,stat=i_stat)
         call memocc(i_stat,i_all,'jfile','allocate_hamiltonianinput')

         i_all=-product(shape(anisotropy))*kind(anisotropy)
         deallocate(anisotropy,stat=i_stat)
         call memocc(i_stat,i_all,'anisotropy','allocate_hamiltonianinput')
         i_all=-product(shape(anisotropytype))*kind(anisotropytype)
         deallocate(anisotropytype,stat=i_stat)
         call memocc(i_stat,i_all,'anisotropytype','allocate_hamiltonianinput')
         if (mult_axis=='Y') then
            i_all=-product(shape(anisotropy_diff))*kind(anisotropy_diff)
            deallocate(anisotropy_diff,stat=i_stat)
            call memocc(i_stat,i_all,'anisotropy_diff','allocate_hamiltonianinput')
            i_all=-product(shape(anisotropytype_diff))*kind(anisotropytype_diff)
            deallocate(anisotropytype_diff,stat=i_stat)
            call memocc(i_stat,i_all,'anisotropytype_diff','allocate_hamiltonianinput')
         endif

         if (do_prn_elk /= 1) then
            i_all=-product(shape(atype_inp))*kind(atype_inp)
            deallocate(atype_inp,stat=i_stat)
            call memocc(i_stat,i_all,'atype_inp','allocate_hamiltonianinput')
         end if

         i_all=-product(shape(anumb_inp))*kind(anumb_inp)
         deallocate(anumb_inp,stat=i_stat)
         call memocc(i_stat,i_all,'anumb_inp','allocate_hamiltonianinput')
         if (allocated(dm_redcoord)) then
            i_all=-product(shape(dm_redcoord))*kind(dm_redcoord)
            deallocate(dm_redcoord,stat=i_stat)
            call memocc(i_stat,i_all,'dm_redcoord','allocate_hamiltonianinput')
         endif
         if (allocated(dm_inpvect)) then
            i_all=-product(shape(dm_inpvect))*kind(dm_inpvect)
            deallocate(dm_inpvect,stat=i_stat)
            call memocc(i_stat,i_all,'dm_inpvect','allocate_hamiltonianinput')
         endif
         if (allocated(pd_redcoord)) then
            i_all=-product(shape(pd_redcoord))*kind(pd_redcoord)
            deallocate(pd_redcoord,stat=i_stat)
            call memocc(i_stat,i_all,'pd_redcoord','allocate_hamiltonianinput')
         endif
         if (allocated(pd_inpvect)) then
            i_all=-product(shape(pd_inpvect))*kind(pd_inpvect)
            deallocate(pd_inpvect,stat=i_stat)
            call memocc(i_stat,i_all,'pd_inpvect','allocate_hamiltonianinput')
         endif
         if (allocated(chir_nn)) then
            i_all=-product(shape(chir_nn))*kind(chir_nn)
            deallocate(chir_nn,stat=i_stat)
            call memocc(i_stat,i_all,'chir_nn','allocate_hamiltonianinput')
         endif
      end if

   end subroutine allocate_hamiltonianinput

   !--------------------------------------------------------------------------------
   !> @brief
   !> Allocating the necessary parameters to be able to read the barriers
   !>
   !> @author
   !> Jonathan Chico
   !>
   !> @date 03/02/2017 - Jonathan Chico
   !--------------------------------------------------------------------------------
   subroutine allocate_barriers(no_shells_barriers, flag) !NA, limit_no_shells, Nchmax, flag)

      use KMCData
      use Profiling
      use Parameters

      implicit none

      integer, intent(in),optional :: no_shells_barriers !< Parameter limiting number of exchange coupling shells
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then

         allocate(kmc_barriers(NT,no_shells_barriers,Nchmax,Nchmax,1),stat=i_stat)
         call memocc(i_stat,product(shape(kmc_barriers))*kind(kmc_barriers),'kmc_barriers','allocate_barriers')
         kmc_barriers=0.0_dblprec
         allocate(redcoord_barriers(NT,no_shells_barriers,3),stat=i_stat)
         call memocc(i_stat,product(shape(redcoord_barriers))*kind(redcoord_barriers),'redcoord_barriers','allocate_barriers')
         redcoord_barriers=0.0_dblprec
         allocate(NN_barriers(NT),stat=i_stat)
         call memocc(i_stat,product(shape(NN_barriers))*kind(NN_barriers),'NN_barriers','allocate_barriers')
         NN_barriers=0

      else

         i_all=-product(shape(kmc_barriers))*kind(kmc_barriers)
         deallocate(kmc_barriers,stat=i_stat)
         call memocc(i_stat,i_all,'kmc_barriers','allocate_barriers')
         i_all=-product(shape(NN_barriers))*kind(NN_barriers)
         deallocate(NN_barriers,stat=i_stat)
         call memocc(i_stat,i_all,'NN_barriers','allocate_barriers')
         i_all=-product(shape(barrfile))*kind(barrfile)

      end if

   end subroutine allocate_barriers

end module InputHandler
