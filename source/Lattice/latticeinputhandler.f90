!---------------------------------------------------------------------------
!> @brief
!> Augments module InputHandler with read routines for lattice dynamics
!
!> @author
!> Johan Hellsvik
!---------------------------------------------------------------------------
module LatticeInputHandler

   use Profiling
   use Parameters
   use ErrorHandling

   use InputData
   use LatticeInputData

   use InputHandler_ext, only :   read_exchange_getMaxNoShells

   implicit none

   public

contains


   !--------------------------------------------------------------------------------
   !> @brief
   !> Read SLD/LD input parameters
   !
   !> @author
   !> Anders Bergman
   !>
   !> @date 08/02/2017 - Jonathan Chico
   !> - Reorganized input variables in blocks by where they are used.
   !> IMPORTANT TRY TO KEEP ORDER IN THE VARIABLES
   !--------------------------------------------------------------------------------
   subroutine read_parameters_sld(ifile)
      use FileParser
!   use Correlation
   use ChemicalData
   !use prn_averages
   !use prn_trajectories
   use prn_latticefields
   !use LatticeCorrelation
   use prn_latticeaverages
   use prn_latticetrajectories

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword, cache
      integer :: rd_len, i_err, i, i_stat, i_errb
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

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR IONS
            !------------------------------------------------------------------------

            case('phonfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               phonfile=adjustl(trim(cache))

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR IONS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR LATTICE DYNAMICS
            !------------------------------------------------------------------------

            case('do_ld')
               read(ifile,*,iostat=i_err) do_ld
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lattdamp')
               read(ifile,*,iostat=i_err) lattdamp
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_velrsc')
               read(ifile,*,iostat=i_err) do_velrsc
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('velrsc_step')
               read(ifile,*,iostat=i_err) velrsc_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('velrsc_taut')
               read(ifile,*,iostat=i_err) velrsc_taut
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_set_avrgp0')
               read(ifile,*,iostat=i_err) do_set_avrgp0
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_set_avrgu0')
               read(ifile,*,iostat=i_err) do_set_avrgu0
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ll')
               do_ll=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               llfile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(llfile, &
                  'Please specify ll <llfile> where <llfile> is a valid LL interaction file')

            case('lll')
               do_lll=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               lllfile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(lllfile, &
                  'Please specify lll <lllfile> where <lllfile> is a valid LLL interaction file')

            case('ml')
               do_ml=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               mlfile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(mlfile, &
                  'Please specify mml <mmlfile> where <mmlfile> is a valid MML interaction file')

            case('mml')
               do_mml=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               mmlfile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(mmlfile, &
                  'Please specify mml <mmlfile> where <mmlfile> is a valid MML interaction file')

            case('mml_scale')
               read(ifile,*,iostat=i_err) mml_scale
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mmll')
               do_mmll=1
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               mmllfile=trim(adjustl(cache))
               call ErrorHandling_check_file_exists(mmllfile, &
                  'Please specify mmll <mmllfile> where <mmllfile> is a valid MMLL interaction file')

            case('mml_diag')
               read(ifile,*,iostat=i_err) mml_diag
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('mml_ene_opt')
               read(ifile,*,iostat=i_err) mml_ene_opt
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ll_phonopy')
               do_ll_phonopy=1
               read(ifile,'(a)',iostat=i_err) cache
               write(*,*) 'cache ', cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ll_phonopyfile=trim(adjustl(cache))

            case('ll_phonopycoordfile')
               read(ifile,'(a)',iostat=i_err) cache
               write(*,*) 'cache ', cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               ll_phonopycoordfile=trim(adjustl(cache))

            case('i0phonopy')
               read(ifile,*,iostat=i_err) i0phonopy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('radius_phonopy')
               read(ifile,*,iostat=i_err) radius_phonopy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('scalefac_phonopy')
               read(ifile,*,iostat=i_err) scalefac_phonopy
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('imp_epsilon')
               read(ifile,*,iostat=i_err) imp_epsilon
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('imp_max_count')
               read(ifile,*,iostat=i_err) imp_max_count
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR LATTICE DYNAMICS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR INITIAL LATTICE CONF
            !------------------------------------------------------------------------

            case('initlatt')
               read(ifile,*,iostat=i_err) initlatt
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('initexc')
               read(ifile,*,iostat=i_err) initexc
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('initconc')
               read(ifile,*,iostat=i_err) initconc
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('initneigh')
               read(ifile,*,iostat=i_err) initneigh
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('initimp')
               read(ifile,*,iostat=i_err) initimp
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lattrestartfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               lattrestartfile=trim(adjustl(cache))

            case('lattroteul')
               read(ifile,*,iostat=i_err) lattroteul
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lattrotang')
               read(ifile,*,iostat=i_err) lattrotang
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR INITAL LATTICE CONF
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR CORRECTIONS TO NEWTON
            !------------------------------------------------------------------------

            case('do_n3')
               read(ifile,*,iostat=i_err) do_n3
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR CORRECTIONS TO NEWTON
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING LATTICE AVERAGES
            !------------------------------------------------------------------------

            case('do_lavrg')
               read(ifile,*,iostat=i_err) do_lavrg
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_proj_lavrg')
               read(ifile,*,iostat=i_err) do_proj_lavrg
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_projch_lavrg')
               read(ifile,*,iostat=i_err) do_projch_lavrg
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lavrg_step')
               read(ifile,*,iostat=i_err) lavrg_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lavrg_buff')
               read(ifile,*,iostat=i_err) lavrg_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING LATTICE AVERAGES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING LATTICE TRAJECTORIES
            !------------------------------------------------------------------------

            case('do_ltottraj')
               read(ifile,*,iostat=i_err) do_ltottraj
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ltottraj_step')
               read(ifile,*,iostat=i_err) ltottraj_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ltottraj_buff')
               read(ifile,*,iostat=i_err) ltottraj_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('lntraj')
               read(ifile,*,iostat=i_err) lntraj
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(lntraj>0) then
                  allocate(ltraj_atom(lntraj),stat=i_stat)
                  call memocc(i_stat,product(shape(ltraj_atom))*kind(ltraj_atom),'ltraj_atom','read_parameters')
                  allocate(ltraj_step(lntraj),stat=i_stat)
                  call memocc(i_stat,product(shape(ltraj_step))*kind(ltraj_step),'ltraj_step','read_parameters')
                  allocate(ltraj_buff(lntraj),stat=i_stat)
                  call memocc(i_stat,product(shape(ltraj_buff))*kind(ltraj_buff),'ltraj_buff','read_parameters')
                  do i=1,lntraj
                     read(ifile,*,iostat=i_err) ltraj_atom(i), ltraj_step(i), ltraj_buff(i)
                  end do
               else
                  read(ifile,*)
               end if

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING LATTICE TRAJECTORIES
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING LATTICE THERMAL FIELDS
            !------------------------------------------------------------------------

            case('do_ethermfield') ! Flag to print the thermal field contribution
               read(ifile,*,iostat=i_err) do_ethermfield
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ethermfield_step') ! Time interval between printing the thermal field
               read(ifile,*,iostat=i_err) ethermfield_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ethermfield_buff') ! Buffer size to store thermal field values
               read(ifile,*,iostat=i_err) ethermfield_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING LATTICE THERMAL FIELDS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING LATTICE TOTAL FIELDS
            !------------------------------------------------------------------------

            case('do_prn_eeff') ! Flag to print the total effective field
               read(ifile,*,iostat=i_err) do_prn_eeff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('eeff_step') ! Time interval between printing the field
               read(ifile,*,iostat=i_err) eeff_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('eeff_buff') ! Buffer to save the effective field
               read(ifile,*,iostat=i_err) eeff_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING LATTICE TOTAL FIELDS
            !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING LATTICE HAMILTONIAN FIELDS
            !------------------------------------------------------------------------

            case('do_prn_einteff') ! Flag to print the total effective field
               read(ifile,*,iostat=i_err) do_prn_einteff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('einteff_step') ! Time interval between printing the field
               read(ifile,*,iostat=i_err) einteff_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('einteff_buff') ! Buffer to save the effective field
               read(ifile,*,iostat=i_err) einteff_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING LATTICE HAMILTONIAN FIELDS
            !------------------------------------------------------------------------
!!!             !------------------------------------------------------------------------
!!!             ! START OF VARIABLES FOR UQW
!!!             !------------------------------------------------------------------------
!!! 
!!!             ! This is the flags for the U(q,w)
!!!             case('do_uc')
!!!                read(ifile,*,iostat=i_err) do_uc
!!!                if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
!!! 
!!!             case('do_ur')
!!!                read(ifile,*,iostat=i_err) do_ur
!!!                if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
!!! 
!!!             case('do_uc_proj')
!!!                read(ifile,*,iostat=i_err) do_uc_proj
!!!                if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
!!! 
!!!             case('do_uqt_traj')
!!!                read(ifile,*,iostat=i_err) do_uqt_traj
!!!                if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
!!! 
!!!             case('uc_mode')
!!!                read(ifile,*,iostat=i_err) uc_mode
!!!                if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
!!! 
!!!             case('uc_step')
!!!                read(ifile,*,iostat=i_err) uc_step
!!!                if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
!!! 
!!!             case('uc_nstep')
!!!                read(ifile,*,iostat=i_err) uc_nstep
!!!                if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
!!! 
!!!             case('uc_sep')
!!!                read(ifile,*,iostat=i_err) uc_sep
!!!                if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
!!! 
!!!             !------------------------------------------------------------------------
!!!             ! END OF VARIABLES FOR UQW
!!!             !------------------------------------------------------------------------

            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PHONON SPECTRA
            !------------------------------------------------------------------------

            case('do_phonspec')
               read(ifile,*,iostat=i_err) do_phonspec
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_phondos')
               read(ifile,*,iostat=i_err) do_phondos
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('phondosfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               phondosfile=adjustl(trim(cache))

            case('phondos_sigma')
               read(ifile,*,iostat=i_err) phondos_sigma
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('phondos_freq')
               read(ifile,*,iostat=i_err) phondos_freq
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PHONON SPECTRA
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
   end subroutine read_parameters_sld

   !
   !> @author
   !> Johan Hellsvik
   !---------------------------------------------------------------------------
   subroutine read_phonons()
      use Parameters
      use Profiling
      implicit none

      integer :: i_err,isite, ichem, i_stat


      open(801,file=trim(phonfile))

      !Allocate arrays according to data from phonon input
      allocate(mion_inp(na,nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(mion_inp))*kind(mion_inp),'mion_inp','read_phonons')
      allocate(uvec_inp(3,na,nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(uvec_inp))*kind(uvec_inp),'uvec_inp','read_phonons')
      allocate(vvec_inp(3,na,nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(vvec_inp))*kind(vvec_inp),'vvec_inp','read_phonons')

      i_err=0

      do while(i_err==0)
         read(801,*,iostat=i_err) isite, ichem, mion_inp(isite,ichem), uvec_inp(1:3,isite,ichem), vvec_inp(1:3,isite,ichem)
      end do

      close(801)

   end subroutine read_phonons


   ! The input lattice potential coupling parameters are in units mRyd / (Ang^n),
   ! where n is the order of the coupling. For the harmonic coupling
   ! ll_inptens is in units mRyd / Ang^2
   subroutine read_lldata()
      !
      use Parameters

      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec), dimension(3,3) :: ll_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(3) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      write(*,*) 'ifileno ', ifileno
      open(ifileno, file=llfile)
      !open(801, file=llfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      call read_exchange_getMaxNoShells(no_shells,flines)
      max_no_llshells = no_shells
      call allocate_latthamiltonianinput(no_shells,1)

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_lldata')
      nn_tmp=0

      ! SLDTODO This is probably redundant given the above call to read_exchange_getMaxNoShells
      ! SLDTODO Verify that no_shells is correctly set by read_exchange_getMaxNoShells[_hoc] and remove unnecessary code
      ! Pre-read file to get max no. exchange shells and no. lines
      do
         if(do_ralloy==0) then
            read(ifileno,*,end=200)  isite,jsite
            ichem=1
            jchem=1
         else
            read(ifileno,*,end=200)  isite,jsite, ichem, jchem
         end if
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
      call memocc(i_stat,i_all,'NN_tmp','read_lldata')

      max_no_llshells=no_shells
      allocate(ll_redcoord(NT,max_no_llshells,3,1),stat=i_stat)
      call memocc(i_stat,product(shape(ll_redcoord))*kind(ll_redcoord),'ll_redcoord','read_lldata')
      ll_redcoord = 0_dblprec

      !SLDTODO check if use of NT and Nchmax is consistent with use in Input/inputhandler
      allocate(ll_inptens(9,NT,max_no_llshells,NT,NT),stat=i_stat)
      call memocc(i_stat,product(shape(ll_inptens))*kind(ll_inptens),'ll_inptens','read_lldata')
      ll_inptens=0.0_dblprec

      ! Read force coupling tensors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp, ll_tmp
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, ll_tmp
         end if

         ! Find type of site
         itype=atype_inp(isite)
         !jtype=1
         jtype=atype_inp(jsite)

         if(maptype==2) then
            ! Calculate proper neighbour vector (from "bgfm")
            r_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(1)+C2(1)*r_tmp(2)+C3(1)*r_tmp(3)
            r_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(1)+C2(2)*r_tmp(2)+C3(2)*r_tmp(3)
            r_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(1)+C2(3)*r_tmp(2)+C3(3)*r_tmp(3)
         else
            ! Calculates neighbour vectors from direct coordinates or Cartesian
            ! coordinates, corresponding to how the atomic positions are entered
            if (posfiletype=='C') then
               r_red=r_tmp
            elseif (posfiletype=='D') then
               r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
               r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
               r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
            else
               stop 'Only posfiletype= C or D is currently supported'
            endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,ll_nn(itype)
            norm=(r_red(1)-ll_redcoord(itype,ishell,1,1))**2+ &
               (r_red(2)-ll_redcoord(itype,ishell,2,1))**2+ &
               (r_red(3)-ll_redcoord(itype,ishell,3,1))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  ll_inptens(1:9,itype,ishell,ichem,jtype)=reshape(ll_tmp,(/9/))
               else
                  ll_inptens(1:9,itype,ishell,ichem,jchem)=reshape(ll_tmp,(/9/))
               end if
            end if
         end do
         if (unique) then
            ll_nn(itype)=ll_nn(itype)+1
            ll_redcoord(itype,ll_nn(itype),1:3,1)=r_red(1:3)
            if(do_ralloy==0) then
               ll_inptens(1:9,itype,ll_nn(itype),ichem,jtype)=reshape(ll_tmp,(/9/))
            else
               ll_inptens(1:9,itype,ll_nn(itype),ichem,jchem)=reshape(ll_tmp,(/9/))
            end if
         end if
      enddo
      close (ifileno)

      !!! print '(2x,i4)',1
      !!! print '(9f10.6)',ll_inptens(:,1,:,:,:)
      !!! print '(2x,i4)',2
      !!! print '(9f10.6)',ll_inptens(:,2,:,:,:)
   end subroutine read_lldata


   ! Read harmonic lattice potential on phonopy format
   subroutine read_llphonopydata()
      !
      implicit none
      !
      integer :: i_stat
      integer :: itype, jtype, ichem, jchem, i_all
      logical :: unique
      real(dblprec) :: tol, cutoffradius, normph
      real(dblprec), dimension(3) :: r_tmp, r_red
      integer :: N1p, N2p, N3p

      integer :: i, j, i_atom, k, l
      real(dblprec), &
         dimension(:,:), allocatable :: cart_coord_phonopy !< LL phonopy force coupling Cartesian coordinates
      real(dblprec), &
         dimension(:,:,:), allocatable :: ll_inptens_phonopy0 !< LL phonopy force constants \Phi_{II}
      real(dblprec), &
         dimension(:,:,:), allocatable :: ll_inptens_phonopy1 !< LL phonopy force constants sum_I/=J \Phi_{IJ}

      real(dblprec), &
         dimension(:,:,:,:), allocatable :: tmp_redcoord    !< Temporary neighbour vectors for LL
      real(dblprec), &
         dimension(:,:,:,:,:), allocatable :: tmp_inptens   !< Temporary coupling tensor for LL

      ! Set tolerance for neighbour shells
      tol=1.0d-5

      ! Open input files
      open(801, file=ll_phonopyfile)
      open(102, file=ll_phonopycoordfile)

      !read(801,*) tmpreal, Natom_phonopy
      read(801,*) Natom_phonopy
      write(*,*) 'Natom_phonopy ', Natom_phonopy

      allocate(atomindex_phonopy(Natom_phonopy),stat=i_stat)
      call memocc(i_stat,product(shape(atomindex_phonopy))*kind(atomindex_phonopy),'atomindex_phonopy','read_llphonopydata')
      allocate(ll_coord_phonopy(3,Natom_phonopy),stat=i_stat)
      call memocc(i_stat,product(shape(ll_coord_phonopy))*kind(ll_coord_phonopy),'ll_coord_phonopy','read_llphonopydata')
      allocate(cart_coord_phonopy(3,Natom_phonopy),stat=i_stat)
      call memocc(i_stat,product(shape(cart_coord_phonopy))*kind(cart_coord_phonopy),'cart_coord_phonopy','read_llphonopydata')
      allocate(ll_inptens_phonopy(3,3,Natom_phonopy,Natom_phonopy),stat=i_stat)
      call memocc(i_stat,product(shape(ll_inptens_phonopy))*kind(ll_inptens_phonopy),'ll_inptens_phonopy','read_llphonopydata')
      allocate(ll_inptens_phonopy0(3,3,Natom_phonopy),stat=i_stat)
      call memocc(i_stat,product(shape(ll_inptens_phonopy0))*kind(ll_inptens_phonopy0),'ll_inptens_phonopy','read_llphonopydata')
      allocate(ll_inptens_phonopy1(3,3,Natom_phonopy),stat=i_stat)
      call memocc(i_stat,product(shape(ll_inptens_phonopy1))*kind(ll_inptens_phonopy1),'ll_inptens_phonopy','read_llphonopydata')

      do i=1,Na  !PHONOPYREVISE DONE
      !do i=1,Natom_phonopy
         do j=1,Natom_phonopy
            read(801,*)
            read(801,*) ll_inptens_phonopy(1:3,1,j,i)
            read(801,*) ll_inptens_phonopy(1:3,2,j,i)
            read(801,*) ll_inptens_phonopy(1:3,3,j,i)
         end do
      end do

      do i=1,Natom_phonopy
         ! SLDTODO remove the running index that right now goes into the first column of POSCAR.dat files
         read(102,*) atomindex_phonopy(i), ll_coord_phonopy(1:3,i)
         r_tmp(1:3) = scalefac_phonopy * ll_coord_phonopy(1:3,i)
         cart_coord_phonopy(1,i)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
         cart_coord_phonopy(2,i)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
         cart_coord_phonopy(3,i)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
      end do
      
      close(801)
      close(102)

      ll_inptens_phonopy0 = 0_dblprec
      ll_inptens_phonopy1 = 0_dblprec

      ! If an interaction radius is used, the calculation/correction of the self-interaction has to be performed later.

      ! Raw input of phonopy data, i.e. the UppASD supercell has same size and shape as the phonopy supercell.
      if(i0phonopy .eq. 0) then 
         do i=1,Natom_phonopy
            ll_inptens_phonopy0(1:3,1:3,i) = ll_inptens_phonopy(1:3,1:3,i,i);
            do j=1,Natom_phonopy
               if(i/=j) then
                  ll_inptens_phonopy1(1:3,1:3,i) = ll_inptens_phonopy1(1:3,1:3,i) + ll_inptens_phonopy(1:3,1:3,j,i);
               end if
            end do
         end do
         do i=1,Natom_phonopy
            ll_inptens_phonopy(1:3,1:3,i,i) = -ll_inptens_phonopy1(1:3,1:3,i)
         end do

         ! Input phonopy data for central atom i0phonopy
      else

         allocate(tmp_redcoord(NT,Natom_phonopy,3,1),stat=i_stat)
         call memocc(i_stat,product(shape(tmp_redcoord))*kind(tmp_redcoord),'tmp_redcoord','read_llphonopydata')
         tmp_redcoord = 0_dblprec

         allocate(tmp_inptens(9,NT,Natom_phonopy,NT,NT),stat=i_stat)
         call memocc(i_stat,product(shape(tmp_inptens))*kind(tmp_inptens),'tmp_inptens','read_llphonopydata')
         tmp_inptens = 0_dblprec

         !i_all=-product(shape(ll_redcoord))*kind(ll_redcoord)
         !deallocate(ll_redcoord,stat=i_stat)
         !call memocc(i_stat,i_all,'ll_redcoord','read_llphonopydata')

         !max_no_llshells=Natom_phonopy
         !allocate(ll_redcoord(NT,max_no_llshells,3,1),stat=i_stat)
         !call memocc(i_stat,product(shape(ll_redcoord))*kind(ll_redcoord),'ll_redcoord','read_llphonopydata')
         !ll_redcoord = 0_dblprec

         !i_all=-product(shape(ll_inptens))*kind(ll_inptens)
         !deallocate(ll_inptens,stat=i_stat)
         !call memocc(i_stat,i_all,'ll_inptens','read_llphonopydata')

         !allocate(ll_inptens(9,NT,max_no_llshells,NT,NT),stat=i_stat)
         !call memocc(i_stat,product(shape(ll_inptens))*kind(ll_inptens),'ll_inptens','read_llphonopydata')
         !ll_inptens = 0_dblprec
         
         ll_redcoord = 0_dblprec !!! Now the bond vectors of regular ll mapping is overwritten !!!
         normph = sqrt( C1(1)**2 + C2(1)**2 + C3(1)**2 )
         cutoffradius = radius_phonopy * scalefac_phonopy * normph

         if(do_hoc_debug==1) then
            write(*,'(a,f10.6)') 'cutoffradius ', cutoffradius
         end if

         itype = 1
         jtype = 1
         ichem = 1
         jchem = 1

         !r_tmp(1:3) = scalefac_phonopy * ll_coord_phonopy(1:3,i0phonopy)   !PHONOPYREVISE
         !r0(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)   !PHONOPYREVISE
         !r0(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)   !PHONOPYREVISE
         !r0(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)   !PHONOPYREVISE

         ll_nn = 0
         !ll_nn(itype) = 0

         N1p=int(scalefac_phonopy)
         N2p=int(scalefac_phonopy)
         N3p=int(scalefac_phonopy)
         
         do i=1,NA
            !do k=1, Natom_phonopy
            !   r_tmp(1:3) = scalefac_phonopy * ll_coord_phonopy(1:3,k)  !PHONOPYREVISE
            !   r0(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)  !PHONOPYREVISE
            !   r0(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)  !PHONOPYREVISE
            !   r0(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)  !PHONOPYREVISE
            !   r1=r0(1:3)-bas(1:3,i)  !PHONOPYREVISE
            !   r1norm=norm2( r1(1:3) )  !PHONOPYREVISE
            !   if ( r1norm .lt. tol ) then  !PHONOPYREVISE
            !      i_atom = k  !PHONOPYREVISE
            !      itype = atype_inp(i_atom)  !PHONOPYREVISE
            !   end if  !PHONOPYREVISE
            !end do  !PHONOPYREVISE
            i_atom=int(scalefac_phonopy**3)*(i-1)+1  !PHONOPYREVISE DONE
            itype = atype_inp(i)  !PHONOPYREVISE DONE
            write(*,*) 'i ', i, ' i_atom ', i_atom, ' i_type ', itype  !PHONOPYREVISE DONE
            do j=1,Natom_phonopy
               
               ! Calculates neighbour vectors from direct coordinates or Cartesian
               ! coordinates, corresponding to how the atomic positions are entered
               !if (posfiletype=='C') then
               !      r_red=r_tmp
               !   elseif (posfiletype=='D') then
               !r_tmp(1:3) = scalefac_phonopy * ( ll_coord_phonopy(1:3,j) - ll_coord_phonopy(1:3,i_atom) )
               r_tmp(1:3) = scalefac_phonopy * ( ll_coord_phonopy(1:3,j) - ll_coord_phonopy(1:3,i_atom) )  !PHONOPYREVISE DONE
               r_red(1:3) = r_tmp(1:3)
               !r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
               !r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
               !r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
               call wrap_phonopy_coord_diff(Natom_phonopy, N1p, N2p, N3p, cart_coord_phonopy, i, j, r_red)
               normph = sqrt(r_red(1)**2 + r_red(2)**2 + r_red(3)**2)
               !if(do_hoc_debug==1) then
               !   write(*,'(a,i4,a,3f10.6,a,f10.6)') 'j ', j, ' r_red ', r_red(1:3), ' normph', normph
               !end if
               !   else
               !      stop 'Only posfiletype= C or D is currently supported'
               !   endif
               !end if
               
               ! Check if shell is within the interaction radius
               if(normph < cutoffradius) then
                  unique=.true.
               else
                  unique=.false.
               end if
               ! Store neighbor vectors and coupling constants to temporary variables
               if (unique) then
                  !write(*,*) 'hit ', j
                  ll_nn(itype)=ll_nn(itype)+1
                  tmp_redcoord(itype,ll_nn(itype),1:3,1)=r_red(1:3)
                  !ll_redcoord(itype,ll_nn(itype),1:3,1)=r_red(1:3)
                  if(do_ralloy==0) then
                     tmp_inptens(1:9,itype,ll_nn(itype),ichem,jtype)=reshape(ll_inptens_phonopy(1:3,1:3,j,i),(/9/))  !PHONOPYREVISE DONE
                     !tmp_inptens(1:9,itype,ll_nn(itype),ichem,jtype)=reshape(ll_inptens_phonopy(1:3,1:3,j,i_atom),(/9/))
                     !ll_inptens(1:9,itype,ll_nn(itype),ichem,jtype)=reshape(ll_inptens_phonopy(1:3,1:3,j,i_atom),(/9/))
                  else
                     tmp_inptens(1:9,itype,ll_nn(itype),ichem,jchem)=reshape(ll_inptens_phonopy(1:3,1:3,j,i),(/9/))   !PHONOPYREVISE DONE
                     !tmp_inptens(1:9,itype,ll_nn(itype),ichem,jchem)=reshape(ll_inptens_phonopy(1:3,1:3,j,i_atom),(/9/))
                     !ll_inptens(1:9,itype,ll_nn(itype),ichem,jchem)=reshape(ll_inptens_phonopy(1:3,1:3,j,i_atom),(/9/))
                  end if
               end if
               
            end do
         end do

         max_no_llshells=maxval(ll_nn(1:NT))
         write(*,*) 'max_no_llshells ', max_no_llshells
         i_all=-product(shape(ll_redcoord))*kind(ll_redcoord)
         deallocate(ll_redcoord,stat=i_stat)
         call memocc(i_stat,i_all,'ll_redcoord','read_llphonopydata')

         allocate(ll_redcoord(NT,max_no_llshells,3,1),stat=i_stat)
         call memocc(i_stat,product(shape(ll_redcoord))*kind(ll_redcoord),'ll_redcoord','read_llphonopydata')
         ll_redcoord = 0_dblprec

         i_all=-product(shape(ll_inptens))*kind(ll_inptens)
         deallocate(ll_inptens,stat=i_stat)
         call memocc(i_stat,i_all,'ll_inptens','read_llphonopydata')

         allocate(ll_inptens(9,NT,max_no_llshells,NT,NT),stat=i_stat)
         call memocc(i_stat,product(shape(ll_inptens))*kind(ll_inptens),'ll_inptens','read_llphonopydata')
         ll_inptens = 0_dblprec

         ! Copy neighbor vectors and coupling constants to ll_inptens and ll_redcoord
         do i=1,NT
            do j=1,max_no_llshells
               do k=1,3
                  ll_redcoord(i,j,k,1) = tmp_redcoord(i,j,k,1)
               end do
               do l=1,9
                  ll_inptens(l,i,j,1,1) = tmp_inptens(l,i,j,1,1)
               end do
            end do
         end do

         ! Deallocate the temporary variables
         i_all=-product(shape(tmp_redcoord))*kind(tmp_redcoord)
         deallocate(tmp_redcoord,stat=i_stat)
         call memocc(i_stat,i_all,'tmp_redcoord','read_llphonopydata')

         i_all=-product(shape(tmp_inptens))*kind(tmp_inptens)
         deallocate(tmp_inptens,stat=i_stat)
         call memocc(i_stat,i_all,'tmp_inptens','read_llphonopydata')


         write(*,*) 'll_redcoord and ll_inptens'
         do i=1,NT
            write(501,'(a,i8,a,2x,i8)') 'Type', i, 'max_no_llshells', max_no_llshells
            do j=1,max_no_llshells
               write(501,'(i12,9f12.6)') j, ll_redcoord(i,j,1:3,1)
               write(501,'(9f12.6)') ll_inptens(1:9,i,j,1,1)
               write(502,'(2i6,13f12.6)') i, j, ll_redcoord(i,j,1:3,1), ll_inptens(1:9,i,j,1,1), norm2(ll_redcoord(i,j,1:3,1))
            end do
         end do

         !integer :: i, j, i_atom, j_atom, k
         !ll_redcoord(1:NT,max_no_llshells,1:3,1) = tmp_redcoord(1:NT,max_no_llshells,1:3,1)
         !ll_inptens(1:9,1:NT,max_no_llshells,1,1) =  tmp_inptens(1:9,1:NT,max_no_llshells,1,1)
         
         !write(*,*) 'll_nn ', ll_nn(itype)

      end if

   end subroutine read_llphonopydata


   subroutine wrap_phonopy_coord_diff(Natom, N1p, N2p, N3p, coord, i_atom, j_atom, cdiff)
     !
     !
     implicit none
     !
     integer, intent(in) :: Natom
     integer, intent(in) :: N1p, N2p, N3p
     real(dblprec), dimension(3,Natom), intent(in) :: coord
     integer, intent(in) :: i_atom
     integer, intent(in) :: j_atom
     real(dblprec), dimension(3), intent(out) :: cdiff
     !
     real(dblprec), dimension(3) :: odiff, oshift, mdiff
     integer :: x,y,z
     integer :: xmin,xmax,ymin,ymax,zmin,zmax
     !
     odiff=coord(:,i_atom) - coord(:,j_atom)
     !
     xmax=1;xmin=-1;ymax=1;ymin=-1;zmax=1;zmin=-1
     !
     mdiff=odiff
     if(do_hoc_debug==1) write(*,*) 'mdiff ', mdiff
     do z=zmin,zmax
        do y=ymin,ymax
           do x=xmin,xmax
              oshift = odiff + x*(N1p)*C1 + y*(N2p)*C2 + z*(N3p)*C3
              !write(*,*) 'oshift ', oshift
              !oshift = odiff + x*(N1)*C1 + y*(N2)*C2 + z*(N3)*C3
              if(norm2(oshift)<norm2(mdiff))  mdiff = oshift
           end do
        end do
     end do
     !
     !print '(2i6,6f10.6)',i_atom,j_atom,mdiff!, oshift
     cdiff=mdiff
     return
     !
   end subroutine wrap_phonopy_coord_diff
   
   
   subroutine read_llldata()
      !
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec), dimension(3,3,3) :: lll_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(6) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp

      integer :: i, nskip

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=lllfile)
      !open(801, file=lllfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      call read_exchange_getMaxNoShells_hoc(no_shells,flines,2,nskip)
      max_no_lllshells = no_shells
      !SLDTODO remove duplicated calls to allocate_latthamiltonianinput
      !call allocate_latthamiltonianinput(no_shells,1)

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_llldata')
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
      call memocc(i_stat,i_all,'NN_tmp','read_llldata')

      max_no_lllshells=no_shells
      allocate(lll_redcoord(NT,max_no_lllshells,3,2),stat=i_stat)
      call memocc(i_stat,product(shape(lll_redcoord))*kind(lll_redcoord),'lll_redcoord','read_llldata')

      allocate(lll_inptens(27,NT,max_no_lllshells,NT,NT),stat=i_stat)
      call memocc(i_stat,product(shape(lll_inptens))*kind(lll_inptens),'lll_inptens','read_llldata')

      ! Read force coupling vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         ! A block of four lines is used for each input lll-tensor
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp
            read (ifileno,*) lll_tmp(1:3,1:3,1)
            read (ifileno,*) lll_tmp(1:3,1:3,2)
            read (ifileno,*) lll_tmp(1:3,1:3,3)
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp
            read (ifileno,*) lll_tmp(1:3,1:3,1)
            read (ifileno,*) lll_tmp(1:3,1:3,2)
            read (ifileno,*) lll_tmp(1:3,1:3,3)
         end if

         ! Find type of site
         itype=atype_inp(isite)
         jtype=atype_inp(jsite)

         !SLDTODO Check that coordinates are calculated correctly for maptype 2!
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
            if (posfiletype=='C') then
               r_red=r_tmp
            elseif (posfiletype=='D') then
               r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
               r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
               r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
               r_red(4)=r_tmp(4)*C1(1)+r_tmp(5)*C2(1)+r_tmp(6)*C3(1)
               r_red(5)=r_tmp(4)*C1(2)+r_tmp(5)*C2(2)+r_tmp(6)*C3(2)
               r_red(6)=r_tmp(4)*C1(3)+r_tmp(5)*C2(3)+r_tmp(6)*C3(3)
            else
               stop 'Only posfiletype= C or D is currently supported'
            endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,lll_nn(itype)
            norm=(r_red(1)-lll_redcoord(itype,ishell,1,1))**2+ &
               (r_red(2)-lll_redcoord(itype,ishell,2,1))**2+ &
               (r_red(3)-lll_redcoord(itype,ishell,3,1))**2+ &
               (r_red(4)-lll_redcoord(itype,ishell,1,2))**2+ &
               (r_red(5)-lll_redcoord(itype,ishell,2,2))**2+ &
               (r_red(6)-lll_redcoord(itype,ishell,3,2))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  lll_inptens(1:27,itype,ishell,ichem,jtype)=reshape(lll_tmp,(/27/))
               else
                  lll_inptens(1:27,itype,ishell,ichem,jchem)=reshape(lll_tmp,(/27/))
               end if
            end if
         end do
         if (unique) then
            lll_nn(itype)=lll_nn(itype)+1
            lll_redcoord(itype,lll_nn(itype),1:3,1)=r_red(1:3)
            lll_redcoord(itype,lll_nn(itype),1:3,2)=r_red(4:6)
            if(do_ralloy==0) then
               Lll_inptens(1:27,itype,lll_nn(itype),ichem,jtype)=reshape(lll_tmp,(/27/))
            else
               lll_inptens(1:27,itype,lll_nn(itype),ichem,jchem)=reshape(lll_tmp,(/27/))
            end if
         end if


      enddo
      close (ifileno)

      if(do_hoc_debug==1) then
         do itype=1,NT
            do jtype=1,NT
               do ishell=1,lll_nn(itype)
                  write(*,*) '----'
                  write(*,'(a,i4,a,i4,a,i4)') 'itype ', itype, ' jtype ', jtype, ' ishell ', ishell
                  write(*,'(a,81f10.6)') 'lll_inptens ', lll_inptens(1:27,itype,ishell,1,jtype)
               end do
            end do
         end do
      end if


   end subroutine read_llldata


   subroutine read_lllldata()
      !
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec), dimension(3,3,3,3) :: llll_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(9) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp

      integer :: i, nskip

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=llllfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      call read_exchange_getMaxNoShells_hoc(no_shells,flines,3,nskip)
      max_no_llllshells = no_shells
      !call allocate_latthamiltonianinput(no_shells,1)

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_lllldata')
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
      call memocc(i_stat,i_all,'NN_tmp','read_lllldata')

      max_no_llllshells=no_shells
      allocate(llll_redcoord(NT,max_no_llllshells,3,3),stat=i_stat)
      call memocc(i_stat,product(shape(llll_redcoord))*kind(llll_redcoord),'llll_redcoord','read_lllldata')

      allocate(llll_inptens(81,NT,max_no_llllshells,NT,NT),stat=i_stat)
      call memocc(i_stat,product(shape(llll_inptens))*kind(llll_inptens),'llll_inptens','read_lllldata')

      ! Read force coupling vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         ! A block of ten lines is used for each input llll-tensor
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp
            read (ifileno,*) llll_tmp(1:3,1:3,1,1)
            read (ifileno,*) llll_tmp(1:3,1:3,2,1)
            read (ifileno,*) llll_tmp(1:3,1:3,3,1)
            read (ifileno,*) llll_tmp(1:3,1:3,1,2)
            read (ifileno,*) llll_tmp(1:3,1:3,2,2)
            read (ifileno,*) llll_tmp(1:3,1:3,3,2)
            read (ifileno,*) llll_tmp(1:3,1:3,1,3)
            read (ifileno,*) llll_tmp(1:3,1:3,2,3)
            read (ifileno,*) llll_tmp(1:3,1:3,3,3)
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, llll_tmp
            read (ifileno,*) llll_tmp(1:3,1:3,1,1)
            read (ifileno,*) llll_tmp(1:3,1:3,2,1)
            read (ifileno,*) llll_tmp(1:3,1:3,3,1)
            read (ifileno,*) llll_tmp(1:3,1:3,1,2)
            read (ifileno,*) llll_tmp(1:3,1:3,2,2)
            read (ifileno,*) llll_tmp(1:3,1:3,3,2)
            read (ifileno,*) llll_tmp(1:3,1:3,1,3)
            read (ifileno,*) llll_tmp(1:3,1:3,2,3)
            read (ifileno,*) llll_tmp(1:3,1:3,3,3)
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
            r_red(7)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(7)+C2(1)*r_tmp(8)+C3(1)*r_tmp(9)
            r_red(8)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(7)+C2(2)*r_tmp(8)+C3(2)*r_tmp(9)
            r_red(9)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(7)+C2(3)*r_tmp(8)+C3(3)*r_tmp(9)
         else
            ! Calculates neighbour vectors from direct coordinates or Cartesian
            ! coordinates, corresponding to how the atomic positions are entered
            if (posfiletype=='C') then
               r_red=r_tmp
            elseif (posfiletype=='D') then
               r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
               r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
               r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
               r_red(4)=r_tmp(4)*C1(1)+r_tmp(5)*C2(1)+r_tmp(6)*C3(1)
               r_red(5)=r_tmp(4)*C1(2)+r_tmp(5)*C2(2)+r_tmp(6)*C3(2)
               r_red(6)=r_tmp(4)*C1(3)+r_tmp(5)*C2(3)+r_tmp(6)*C3(3)
               r_red(7)=r_tmp(7)*C1(1)+r_tmp(8)*C2(1)+r_tmp(9)*C3(1)
               r_red(8)=r_tmp(7)*C1(2)+r_tmp(8)*C2(2)+r_tmp(9)*C3(2)
               r_red(9)=r_tmp(7)*C1(3)+r_tmp(8)*C2(3)+r_tmp(9)*C3(3)
            else
               stop 'Only posfiletype= C or D is currently supported'
            endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,llll_nn(itype)
            norm=(r_red(1)-llll_redcoord(itype,ishell,1,1))**2+ &
               (r_red(2)-llll_redcoord(itype,ishell,2,1))**2+ &
               (r_red(3)-llll_redcoord(itype,ishell,3,1))**2+ &
               (r_red(4)-llll_redcoord(itype,ishell,1,2))**2+ &
               (r_red(5)-llll_redcoord(itype,ishell,2,2))**2+ &
               (r_red(6)-llll_redcoord(itype,ishell,3,2))**2+ &
               (r_red(7)-llll_redcoord(itype,ishell,1,3))**2+ &
               (r_red(8)-llll_redcoord(itype,ishell,2,3))**2+ &
               (r_red(9)-llll_redcoord(itype,ishell,3,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  llll_inptens(1:81,itype,ishell,ichem,jtype)=reshape(llll_tmp,(/81/))
               else
                  llll_inptens(1:81,itype,ishell,ichem,jchem)=reshape(llll_tmp,(/81/))
               end if
            end if
         end do
         if (unique) then
            llll_nn(itype)=llll_nn(itype)+1
            llll_redcoord(itype,llll_nn(itype),1:3,1)=r_red(1:3)
            llll_redcoord(itype,llll_nn(itype),1:3,2)=r_red(4:6)
            llll_redcoord(itype,llll_nn(itype),1:3,3)=r_red(7:9)
            if(do_ralloy==0) then
               llll_inptens(1:81,itype,llll_nn(itype),ichem,jtype)=reshape(llll_tmp,(/81/))
            else
               llll_inptens(1:81,itype,llll_nn(itype),ichem,jchem)=reshape(llll_tmp,(/81/))
            end if
         end if
      enddo
      close (ifileno)

   end subroutine read_lllldata


   subroutine read_mldata()
      !
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec), dimension(3,3) :: ml_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(3) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp
      integer :: i, nskip

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=mlfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      !call read_exchange_getMaxNoShells_hoc(no_shells,1,flines,mlfile,1,nskip)
      nskip=0
      !max_no_mlshells = no_shells
      !call allocate_latthamiltonianinput(no_shells,1)

      call read_exchange_getMaxNoShells(no_shells,flines)
      max_no_llshells = no_shells
      !call allocate_latthamiltonianinput(no_shells,1)

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_mldata')
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
      call memocc(i_stat,i_all,'NN_tmp','read_mldata')

      max_no_mlshells=no_shells
      allocate(ml_redcoord(NT,max_no_mlshells,3,1),stat=i_stat)
      call memocc(i_stat,product(shape(ml_redcoord))*kind(ml_redcoord),'ml_redcoord','read_mldata')

      allocate(ml_inptens(9,NT,max_no_mlshells,NT,NT),stat=i_stat)
      call memocc(i_stat,product(shape(ml_inptens))*kind(ml_inptens),'ml_inptens','read_mldata')

      ! Read force coupling vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         ! A block of four lines is used for each input ml-tensor
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp, ml_tmp(1:3,1:3)
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, ml_tmp(1:3,1:3)
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
            if (posfiletype=='C') then
               r_red=r_tmp
            elseif (posfiletype=='D') then
               r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
               r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
               r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
            else
               stop 'Only posfiletype= C or D is currently supported'
            endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,ml_nn(itype)
            norm=(r_red(1)-ml_redcoord(itype,ishell,1,1))**2+ &
               (r_red(2)-ml_redcoord(itype,ishell,2,1))**2+ &
               (r_red(3)-ml_redcoord(itype,ishell,3,1))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  ml_inptens(1:9,itype,ishell,ichem,jtype)=reshape(ml_tmp,(/9/))
               else
                  ml_inptens(1:9,itype,ishell,ichem,jchem)=reshape(ml_tmp,(/9/))
               end if
            end if
         end do
         if (unique) then
            ml_nn(itype)=ml_nn(itype)+1
            ml_redcoord(itype,ml_nn(itype),1:3,1)=r_red(1:3)
            if(do_ralloy==0) then
               ml_inptens(1:9,itype,ml_nn(itype),ichem,jtype)=reshape(ml_tmp,(/9/))
            else
               ml_inptens(1:9,itype,ml_nn(itype),ichem,jchem)=reshape(ml_tmp,(/9/))
            end if
         end if
      enddo
      close (ifileno)

   end subroutine read_mldata


   subroutine read_mmldata()
      !
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec), dimension(3,3,3) :: mml_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(6) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp
      integer :: i, nskip
      !integer :: invsym_tmp

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=mmlfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      call read_exchange_getMaxNoShells_hoc(no_shells,flines,2,nskip)
      max_no_mmlshells = no_shells
      !call allocate_latthamiltonianinput(no_shells,1)

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_mmldata')
      nn_tmp=0

      allocate(mml_invsym(NT, max_no_mmlshells),stat=i_stat)
      call memocc(i_stat,product(shape(mml_invsym))*kind(mml_invsym),'mml_invsym','read_mmldata')
      mml_invsym=0

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
      call memocc(i_stat,i_all,'NN_tmp','read_mmldata')

      max_no_mmlshells=no_shells
      allocate(mml_redcoord(NT,max_no_mmlshells,3,2),stat=i_stat)
      call memocc(i_stat,product(shape(mml_redcoord))*kind(mml_redcoord),'mml_redcoord','read_mmldata')

      allocate(mml_inptens(27,NT,max_no_mmlshells,NT,NT),stat=i_stat)
      call memocc(i_stat,product(shape(mml_inptens))*kind(mml_inptens),'mml_inptens','read_mmldata')

      ! Read force coupling vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         ! A block of four lines is used for each input mml-tensor
         if(do_ralloy==0) then
            !read (ifileno,*) isite, jsite, r_tmp, invsym_tmp
            read (ifileno,*) isite, jsite, r_tmp
            read (ifileno,*) mml_tmp(1:3,1:3,1)
            read (ifileno,*) mml_tmp(1:3,1:3,2)
            read (ifileno,*) mml_tmp(1:3,1:3,3)
            ichem=1
            jchem=1
         else
            !read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, invsym_tmp
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp
            read (ifileno,*) mml_tmp(1:3,1:3,1)
            read (ifileno,*) mml_tmp(1:3,1:3,2)
            read (ifileno,*) mml_tmp(1:3,1:3,3)
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
            if (posfiletype=='C') then
               r_red=r_tmp
            elseif (posfiletype=='D') then
               r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
               r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
               r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
               r_red(4)=r_tmp(4)*C1(1)+r_tmp(5)*C2(1)+r_tmp(6)*C3(1)
               r_red(5)=r_tmp(4)*C1(2)+r_tmp(5)*C2(2)+r_tmp(6)*C3(2)
               r_red(6)=r_tmp(4)*C1(3)+r_tmp(5)*C2(3)+r_tmp(6)*C3(3)
            else
               stop 'Only posfiletype= C or D is currently supported'
            endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,mml_nn(itype)
            norm=(r_red(1)-mml_redcoord(itype,ishell,1,1))**2+ &
               (r_red(2)-mml_redcoord(itype,ishell,2,1))**2+ &
               (r_red(3)-mml_redcoord(itype,ishell,3,1))**2+ &
               (r_red(4)-mml_redcoord(itype,ishell,1,2))**2+ &
               (r_red(5)-mml_redcoord(itype,ishell,2,2))**2+ &
               (r_red(6)-mml_redcoord(itype,ishell,3,2))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  mml_inptens(1:27,itype,ishell,ichem,jtype)=reshape(mml_tmp,(/27/))
               else
                  mml_inptens(1:27,itype,ishell,ichem,jchem)=reshape(mml_tmp,(/27/))
               end if
            end if
         end do
         if (unique) then
            mml_nn(itype)=mml_nn(itype)+1
            mml_redcoord(itype,mml_nn(itype),1:3,1)=r_red(1:3)
            mml_redcoord(itype,mml_nn(itype),1:3,2)=r_red(4:6)
            !mml_invsym(itype,mml_nn(itype))=invsym_tmp
            if(do_ralloy==0) then
               mml_inptens(1:27,itype,mml_nn(itype),ichem,jtype)=reshape(mml_tmp,(/27/))
            else
               mml_inptens(1:27,itype,mml_nn(itype),ichem,jchem)=reshape(mml_tmp,(/27/))
            end if
         end if
      enddo
      close (ifileno)

      if(do_hoc_debug==1) then
         write(*,*) 'The entries of the mml_inptens'
         do itype=1,NT
            do jtype=1,NT
               do ishell=1,mml_nn(itype)
                  write(*,*) '----'
                  write(*,'(a,i4,a,i4,a,i4)') 'itype ', itype, ' jtype ', jtype, ' ishell ', ishell
                  write(*,'(a,81f10.6)') 'mml_inptens ', mml_inptens(1:27,itype,ishell,1,jtype)
               end do
            end do
         end do
      end if


   end subroutine read_mmldata


   subroutine read_mmlldata()
      !
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec), dimension(3,3,3,3) :: mmll_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(9) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp

      integer :: i, nskip

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=mmllfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      call read_exchange_getMaxNoShells_hoc(no_shells,flines,3,nskip)
      max_no_mmllshells = no_shells
      !call allocate_latthamiltonianinput(no_shells,1)

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_mmlldata')
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
      call memocc(i_stat,i_all,'NN_tmp','read_mmlldata')

      max_no_mmllshells=no_shells
      allocate(mmll_redcoord(NT,max_no_mmllshells,3,3),stat=i_stat)
      call memocc(i_stat,product(shape(mmll_redcoord))*kind(mmll_redcoord),'mmll_redcoord','read_mmlldata')

      allocate(mmll_inptens(81,NT,max_no_mmllshells,NT,NT),stat=i_stat)
      call memocc(i_stat,product(shape(mmll_inptens))*kind(mmll_inptens),'mmll_inptens','read_mmlldata')

      ! Read force coupling vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         ! A block of ten lines is used for each input mmll-tensor
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp
            read (ifileno,*) mmll_tmp(1:3,1:3,1,1)
            read (ifileno,*) mmll_tmp(1:3,1:3,2,1)
            read (ifileno,*) mmll_tmp(1:3,1:3,3,1)
            read (ifileno,*) mmll_tmp(1:3,1:3,1,2)
            read (ifileno,*) mmll_tmp(1:3,1:3,2,2)
            read (ifileno,*) mmll_tmp(1:3,1:3,3,2)
            read (ifileno,*) mmll_tmp(1:3,1:3,1,3)
            read (ifileno,*) mmll_tmp(1:3,1:3,2,3)
            read (ifileno,*) mmll_tmp(1:3,1:3,3,3)
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp
            read (ifileno,*) mmll_tmp(1:3,1:3,1,1)
            read (ifileno,*) mmll_tmp(1:3,1:3,2,1)
            read (ifileno,*) mmll_tmp(1:3,1:3,3,1)
            read (ifileno,*) mmll_tmp(1:3,1:3,1,2)
            read (ifileno,*) mmll_tmp(1:3,1:3,2,2)
            read (ifileno,*) mmll_tmp(1:3,1:3,3,2)
            read (ifileno,*) mmll_tmp(1:3,1:3,1,3)
            read (ifileno,*) mmll_tmp(1:3,1:3,2,3)
            read (ifileno,*) mmll_tmp(1:3,1:3,3,3)
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
            r_red(7)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(7)+C2(1)*r_tmp(8)+C3(1)*r_tmp(9)
            r_red(8)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(7)+C2(2)*r_tmp(8)+C3(2)*r_tmp(9)
            r_red(9)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(7)+C2(3)*r_tmp(8)+C3(3)*r_tmp(9)
         else
            ! Calculates neighbour vectors from direct coordinates or Cartesian
            ! coordinates, corresponding to how the atomic positions are entered
            if (posfiletype=='C') then
               r_red=r_tmp
            elseif (posfiletype=='D') then
               r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
               r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
               r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
               r_red(4)=r_tmp(4)*C1(1)+r_tmp(5)*C2(1)+r_tmp(6)*C3(1)
               r_red(5)=r_tmp(4)*C1(2)+r_tmp(5)*C2(2)+r_tmp(6)*C3(2)
               r_red(6)=r_tmp(4)*C1(3)+r_tmp(5)*C2(3)+r_tmp(6)*C3(3)
               r_red(7)=r_tmp(7)*C1(1)+r_tmp(8)*C2(1)+r_tmp(9)*C3(1)
               r_red(8)=r_tmp(7)*C1(2)+r_tmp(8)*C2(2)+r_tmp(9)*C3(2)
               r_red(9)=r_tmp(7)*C1(3)+r_tmp(8)*C2(3)+r_tmp(9)*C3(3)
            else
               stop 'Only posfiletype= C or D is currently supported'
            endif
         end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,mmll_nn(itype)
            norm=(r_red(1)-mmll_redcoord(itype,ishell,1,1))**2+ &
               (r_red(2)-mmll_redcoord(itype,ishell,2,1))**2+ &
               (r_red(3)-mmll_redcoord(itype,ishell,3,1))**2+ &
               (r_red(4)-mmll_redcoord(itype,ishell,1,2))**2+ &
               (r_red(5)-mmll_redcoord(itype,ishell,2,2))**2+ &
               (r_red(6)-mmll_redcoord(itype,ishell,3,2))**2+ &
               (r_red(7)-mmll_redcoord(itype,ishell,1,3))**2+ &
               (r_red(8)-mmll_redcoord(itype,ishell,2,3))**2+ &
               (r_red(9)-mmll_redcoord(itype,ishell,3,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  mmll_inptens(1:81,itype,ishell,ichem,jtype)=reshape(mmll_tmp,(/81/))
               else
                  mmll_inptens(1:81,itype,ishell,ichem,jchem)=reshape(mmll_tmp,(/81/))
               end if
            end if
         end do
         if (unique) then
            mmll_nn(itype)=mmll_nn(itype)+1
            mmll_redcoord(itype,mmll_nn(itype),1:3,1)=r_red(1:3)
            mmll_redcoord(itype,mmll_nn(itype),1:3,2)=r_red(4:6)
            mmll_redcoord(itype,mmll_nn(itype),1:3,3)=r_red(7:9)
            if(do_ralloy==0) then
               mmll_inptens(1:81,itype,mmll_nn(itype),ichem,jtype)=reshape(mmll_tmp,(/81/))
            else
               mmll_inptens(1:81,itype,mmll_nn(itype),ichem,jchem)=reshape(mmll_tmp,(/81/))
            end if
         end if
      enddo
      close (ifileno)

   end subroutine read_mmlldata


   !---------------------------------------------------------------------------
   !> @brief
   !> Get the max no of exchange shells and lines
   !> Helper for read exchange
   !>
   !> @author
   !> Anders Bergman
   !> Adopted to higher order coupling input files, Johan Hellsvik
   !---------------------------------------------------------------------------
   subroutine read_exchange_getMaxNoShells_hoc(no_shells,flines,nelem,nskip)
      integer, intent(out)                   :: no_shells,flines
      integer                                :: mtype
      integer                                :: itype,jtype,isite,jsite,ichem,jchem
      integer                                :: i_stat,i_all

      integer, dimension(:,:,:), allocatable :: nn_tmp
      integer :: nelem
      integer :: nskip
      integer :: i

      ! Open input file
      !    open(101, file=trim(filename(1))) ! Number of shells
      ! Check if input file is for random alloy
      flines=0
      mtype=0

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_exchange')
      nn_tmp=0

      if(nelem==1) nskip=0
      if(nelem==2) nskip=3
      if(nelem==3) nskip=9
      nskip = 3**(nelem-1)
      ! Pre-read file to get max no. exchange shells and no. lines
      do
         if(do_ralloy==0) then
            read(ifileno,*,end=200)  isite,jsite
            ichem=1
            jchem=1
         else
            read(ifileno,*,end=200)  isite, jsite, ichem, jchem
         end if
         ! Skip 0, 3 or 9 lines
         do i=1,nskip
            read(ifileno,*)
         end do
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

   end subroutine read_exchange_getMaxNoShells_hoc


end module LatticeInputHandler
