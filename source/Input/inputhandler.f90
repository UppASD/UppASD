!------------------------------------------------------------------------------------
!  MODULE: InputHandler
!> @brief
!> Reads each keyword in inpsd.dat.
!> @details If the first read character of a line is a %, #, * or = the line is considered
!> to be a comment and and is ignored.
!> Comments after keyword values does not need a %, #, * or =
!
!> @author
!> Anders Bergman, Nikos Ntallis
!> @copyright
!> GNU Public License.
!
!> @todo
!> Put all reads in separate modules
!------------------------------------------------------------------------------------
module InputHandler

   use Parameters
   use Profiling
   use InputData
   use ErrorHandling

   implicit none

   logical :: sane_input = .true.

   private

   public :: sane_input
   public :: read_parameters, change_constants
   !!! public :: allocate_hamiltonianinput, read_exchange_getMaxNoShells, read_parameters
   !!! public :: read_positions, read_positions_alloy, change_constants, read_moments, read_exchange
   !!! public :: read_exchange_tensor, read_exchange_build_tensor, read_anisotropy_alloy, read_anisotropy
   !!! public :: read_dmdata, read_pddata, read_chirdata, read_biqdmdata, read_bqdata, read_ringdata, read_sitefield
   !!! public :: read_ip_damping, read_ip_damping_alloy, read_damping, read_damping_alloy, read_fourxdata
   !!! public :: read_barriers, read_fixed_moments, read_exchange_getNeighVec

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
      use LatticeInputData, only : iplattdamp
      use QHB,                only : do_qhb, qhb_mode, tcurie, do_qhb_mix, qhb_mix_mode
      use KMCData
      use clusters
      use FixedMom,           only : do_fixed_mom
      use stiffness
      use prn_fields
      use macrocells,         only : prn_dip_subset,dip_file
      use temperature,        only : grad, tempfile, do_3tm
      use Polarization
      use prn_topology
      use prn_currents
      use RandomNumbers
      use prn_induced_info,   only : do_prn_induced, ind_step,ind_buff

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword,cache, string
      integer :: rd_len,i_err,i,i_stat,i_errb,ii, i_all
      logical :: comment
      real(dblprec) :: tmp
      character(len=1) :: opt
      integer :: i_dum, j_dum, nlines

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
               read(ifile,*,iostat=i_err) nlines
               if(na.ne.nlines.and.(.not.na==0)) then
                  write(*,*) 'ERROR: ','natoms not consistent with ',keyword
               else
                  if (na==0) na = nlines
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

            case('moments')
               read(ifile,*,iostat=i_err) nlines
               if(na.ne.nlines.and.(.not.na==0)) then
                  write(*,*) 'ERROR: ','natoms not consistent with ',keyword
               else
                  if (na==0) na = nlines
                  allocate(ammom_inp(na,nchmax,conf_num),stat=i_stat)
                  call memocc(i_stat,product(shape(ammom_inp))*kind(ammom_inp),'ammom_inp','read_parameters')
                  ammom_inp=0.0_dblprec

                  allocate(aemom_inp(3,na,nchmax,conf_num),stat=i_stat)
                  call memocc(i_stat,product(shape(aemom_inp))*kind(aemom_inp),'aemom_inp','read_parameters')
                  aemom_inp=0.0_dblprec

                  allocate(Landeg_ch(na,nchmax,conf_num),stat=i_stat)
                  call memocc(i_stat,product(shape(Landeg_ch))*kind(Landeg_ch),'Landeg_ch','read_parameters')
                  Landeg_ch=Landeg_glob

                  do i=1,na
                     read(ifileno,*,iostat=i_err) i_dum, j_dum, ammom_inp(i_dum,j_dum,1), &
                        aemom_inp(1:3,i_dum,j_dum,1)
                     tmp=norm2(aemom_inp(:,i_dum,j_dum,1))
                     aemom_inp(1:3,i_dum,j_dum,1)=aemom_inp(1:3,i_dum,j_dum,1)/tmp
                     if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
                  end do
               end if


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

            case('mavg0')
               read(ifile,*,iostat=i_err) mavg0
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
               read(ifile,*,iostat=i_err) ham_inp%jij_scale
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ea_model')
               read(ifile,*,iostat=i_err) ham_inp%ea_model
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ea_sigma')
               read(ifile,*,iostat=i_err) ham_inp%ea_sigma
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('exc_inter')
               read(ifile,*,iostat=i_err) ham_inp%exc_inter
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ! Type of maptype 1 for cartesian bonding vectors, 2 for units of repetitions of the lattice vectors
            case('maptype')
               read(ifile,*,iostat=i_err) maptype
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('map_multiple')
               read(ifile,*,iostat=i_err) ham_inp%map_multiple
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ! Symmetry of the exchange interactions, turned off by default when DMI is considered
            case('sym')
               read(ifile,*,iostat=i_err) sym
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_sortcoup')
               read(ifile,*,iostat=i_err) do_sortcoup
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_jtensor')
               read(ifile,*,iostat=i_err) ham_inp%do_jtensor
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('calc_jtensor')
               read(ifile,*,iostat=i_err) ham_inp%calc_jtensor
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR EXCHANGE
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR ANISOTROPIES
            !------------------------------------------------------------------------

            case('do_anisotropy')
               read(ifile,*,iostat=i_err) ham_inp%do_anisotropy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('anisotropy')
               ham_inp%do_anisotropy=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%kfile=trim(adjustl(cache))

            case('mult_axis')
               read(ifile,*,iostat=i_err) ham_inp%mult_axis
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('random_anisotropy')
               ham_inp%do_anisotropy=1
               read(ifile,*,iostat=i_err) ham_inp%random_anisotropy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('random_anisotropy_density')
               read(ifile,*,iostat=i_err) ham_inp%random_anisotropy_density
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR ANISOTROPIES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR DMI
            !------------------------------------------------------------------------

            case('dm')
               ham_inp%do_dm=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%dmfile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(ham_inp%dmfile, &
                  'Please specify dm <dmfile> where <dmfile> is a valid dm interaction file')

            case('dm_scale')
               read(ifile,*,iostat=i_err) ham_inp%dm_scale
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('sa')
               ham_inp%do_sa=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%safile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(ham_inp%safile, &
                  'Please specify sa <safile> where <safile> is a valid sa interaction file')

            case('sa_scale')
               read(ifile,*,iostat=i_err) ham_inp%sa_scale
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err


            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR DMI
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR DIPOLAR
            !------------------------------------------------------------------------

            case('do_dip')
               read(ifile,*,iostat=i_err) ham_inp%do_dip
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('print_dip_tensor')
               read(ifile,*,iostat=i_err) ham_inp%print_dip_tensor
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qdip_files')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%qdip_files=trim(adjustl(cache))

            case('read_dipole')
               read(ifile,'(a)',iostat=i_err) ham_inp%read_dipole
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
               ham_inp%do_pd=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%pdfile=trim(adjustl(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PSEUDO DIPOLAR
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR Scalar Chiral
            !------------------------------------------------------------------------

            case('chir')
               ham_inp%do_chir=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%chirfile=trim(adjustl(cache))

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR Scalar four-site interactions
            !------------------------------------------------------------------------

            case('fourx')
               ham_inp%do_fourx=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%fourxfile=trim(adjustl(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR Scalar four-site interactions
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR BIQUADRATIC
            !------------------------------------------------------------------------

            case('bq')
               ham_inp%do_bq=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%bqfile=trim(adjustl(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR BIQUADRATIC
            !------------------------------------------------------------------------
            
            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR FOUR-SPIN RING EXCHANGE
            !------------------------------------------------------------------------

            case('ring')
               ham_inp%do_ring=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%ringfile=trim(adjustl(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR FOUR-SPIN RING
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR BIQUADRATIC DMI
            !------------------------------------------------------------------------

            case('biqdm')
               ham_inp%do_biqdm=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ham_inp%biqdmfile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(ham_inp%biqdmfile, &
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

            case('do_3tm')
               read(ifile,*,iostat=i_err) do_3tm
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
                     allocate(iplattdamp(ipnphase),stat=i_stat)
                     call memocc(i_stat,product(shape(iplattdamp))*kind(iplattdamp),'iplattdamp','read_parameters')
                     do i=1,ipnphase
                        read(ifile,*,iostat=i_err) ipnstep(i),ipTemp(i),ipdelta_t(i),iplambda1(i),iplattdamp(i)
                        !read(ifile,*,iostat=i_err) ipnstep(i),ipTemp(i),ipdelta_t(i),iplambda1(i)
                        if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
                     end do
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

            case('do_storeham')
               read(ifile,*,iostat=i_err) do_storeham
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_prn_poscar')
               read(ifile,*,iostat=i_err) do_prn_poscar
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_hoc_debug')
               read(ifile,*,iostat=i_err) do_hoc_debug
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_meminfo')
               read(ifile,*,iostat=i_err) do_meminfo
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

            case('do_skyno_cmass')
               read(ifile,*,iostat=i_err) do_skyno_cmass
               do_skyno_den='Y'
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

            ! Mix scheme inside QHB module

            case('do_qhb_mix')
               read(ifile,*,iostat=i_err) do_qhb_mix
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('qhb_mix_mode')
               read(ifile,*,iostat=i_err) qhb_mix_mode
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

            case('ind_mom_type') ! Flag to indicate that there are induced moments being considered
               read(ifile,*,iostat=i_err) ind_mom_type
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
               read(ifile,'(a)',iostat=i_err) ham_inp%do_ewald
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err

            case('ewald_alpha')
               read(ifile,*,iostat=i_err) ham_inp%Ewald_alpha
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err

            case('rmax_ewald')
               read(ifile,*,iostat=i_err) ham_inp%RMAX
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err

            case('kmax_ewald')
               read(ifile,*,iostat=i_err) ham_inp%KMAX
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
             
            case('multiscale')
              read(ifile,*,iostat=i_err) multiscale_file_name
              if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err
              do_multiscale = .true.

            case('prn_multiscale')
               read(ifile,*,iostat=i_err) opt
               if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err
               do_prnmultiscale = (opt == 'T' .or. opt == 't' .or. opt == '1')
            !------------------------------------------------------------------------
            ! END OF LEGACY VARIABLES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF METATYPE VARIABLES
            !------------------------------------------------------------------------

            case('metatype')
               read(ifile,*,iostat=i_err) metatype
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('metanumb')
               read(ifile,*,iostat=i_err) metanumb
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF METATYPE VARIABLES
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



end module InputHandler
