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
module InputHandler_ext

   use Parameters
   use Profiling
   use InputData
   use ErrorHandling
   use KMCData
   use clusters
   !use QHB,                only : do_qhb, qhb_mode, tcurie
   !use FixedMom,           only : inp_fixed_mom_flag, do_fixed_mom
   !use stiffness
   !use prn_fields
   !use macrocells,         only : do_macro_cells,prn_dip_subset,dip_file
   !use temperature,        only : grad, tempfile, do_3tm
   !use Polarization
   !use prn_topology
   !use prn_currents
   !use RandomNumbers
   !use prn_induced_info,   only : do_prn_induced, ind_step,ind_buff

   implicit none


   private

   public :: allocate_hamiltonianinput, read_exchange_getMaxNoShells
   public :: read_positions, read_positions_alloy, read_moments, read_exchange
   public :: read_exchange_tensor, read_exchange_build_tensor, read_anisotropy_alloy, read_anisotropy, read_sadata
   public :: read_dmdata, read_pddata, read_chirdata, read_biqdmdata, read_bqdata, read_ringdata, read_sitefield
   public :: read_ip_damping, read_ip_damping_alloy, read_damping, read_damping_alloy, read_fourxdata
   public :: read_barriers, read_fixed_moments, read_exchange_getNeighVec

contains


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
   subroutine read_exchange(ham_inp)
      !! @todo Consider change so that if the atom basis is given in direct coordinates, then also
      !! @todo the coordinates for the exchange coupling shells are to be stated in direct coordinates.
      !
      use Clusters
      use InputDataType

      implicit none

      type(ham_inp_t), intent(inout) :: ham_inp


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
         if (ham_inp%exc_inter=='Y') open(ifileno2, file=trim(jfileD(ii))) ! Number of shells same for different LSF configuration
         if (ii==1) then
            call read_exchange_getMaxNoShells(no_shells,flines)
            call allocate_hamiltonianinput(ham_inp,no_shells,1)
         endif
         ham_inp%redcoord = 0.0_dblprec
         ham_inp%nn       = 0

         ! Read exchange vectors
         ! Isite, Jsite, Ichem, Jchem
         do iline=1, flines
            ! Loop through earlier vectors to find equivalent shells

            ! Read indices and coordinates
            if(do_ralloy==0) then
               read (ifileno,*) isite, jsite, r_tmp(1:3), j_tmp
               if (ham_inp%exc_inter=='Y') read (ifileno2,*) idum, idum, tmp, tmp, tmp, j_tmpD
               ichem=1
               jchem=1
            else
               read (ifileno,*) isite, jsite, ichem, jchem, r_tmp(1:3), j_tmp
               if (ham_inp%exc_inter=='Y') read (ifileno2,*) idum, idum, idum, idum, tmp, tmp, tmp, j_tmpD
            end if

            itype=atype_inp(isite)
            jtype=1
            call read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

            ! Loop through earlier vectors to find equivalent shells
            unique=.true.
            do ishell=1,ham_inp%nn(itype)
               norm=(r_red(1)-ham_inp%redcoord(itype,ishell,1))**2+ &
                    (r_red(2)-ham_inp%redcoord(itype,ishell,2))**2+ &
                    (r_red(3)-ham_inp%redcoord(itype,ishell,3))**2
               if(norm<tol) then
                  unique=.false.
                  ! If neighbour vector already exist, replace J value (could be removed)
                  if(do_ralloy==0) then
                     ham_inp%jc(itype,ishell,ichem,jtype,ii)=j_tmp
                     if (ham_inp%exc_inter=='Y') ham_inp%jcD(itype,ishell,ichem,jtype,ii)=j_tmpD
                  else
                     ham_inp%jc(itype,ishell,ichem,jchem,ii)=j_tmp
                     if (ham_inp%exc_inter=='Y') ham_inp%jcD(itype,ishell,ichem,jchem,ii)=j_tmpD
                  end if
               end if
            end do
            if (unique) then
               ! Add entry if not found earlier
               ham_inp%nn(itype)=ham_inp%nn(itype)+1
               ham_inp%nntype(itype,ham_inp%nn(itype))=atype_inp(jsite)
               !write(2000,'(2i6, 2i4, 2i8)') isite, ham_inp%nn(itype),ham_inp%nntype(itype,ham_inp%nn(itype))
               ham_inp%redcoord(itype,ham_inp%nn(itype),1:3)=r_red(1:3)
               if(do_ralloy==0) then
                  ham_inp%jc(itype,ham_inp%nn(itype),ichem,jtype,ii)=j_tmp
                  if (ham_inp%exc_inter=='Y') ham_inp%jcD(itype,ham_inp%nn(itype),ichem,jtype,ii)=j_tmpD
               else
                  ham_inp%jc(itype,ham_inp%nn(itype),ichem,jchem,ii)=j_tmp
                  if (ham_inp%exc_inter=='Y') ham_inp%jcD(itype,ham_inp%nn(itype),ichem,jchem,ii)=j_tmpD
               end if
            end if

         enddo
         close(ifileno)

         ! Reading the cluster files
         if (do_cluster=='Y') then

            if (do_cluster=='Y'.and.ham_inp%exc_inter=='Y')then
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
                  if(ham_inp%exc_inter=='Y') then
                     read (ifileno3,*) isite_c, jsite_c, r_tmp_clus(1:3), j_tmp_clus
                  else
                     read (ifileno2,*) isite_c, jsite_c, r_tmp_clus(1:3), j_tmp_clus
                  endif
                  ichem_c=1
                  jchem_c=1
               else
                  if (ham_inp%exc_inter=='Y') then
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
         ham_inp%max_no_shells=maxval(ham_inp%NN)
         if(ham_inp%exc_inter=='Y') close(ifileno2)
         if(do_cluster=='Y') then
            max_no_shells_clus=maxval(NN_clus)
            close(ifileno3)
         endif
      enddo

      call read_exchange_reduceNNtypeMatrixSize(ham_inp%nntype,nt,ham_inp%max_no_shells)
      call read_exchange_reduceRedCoordMatrixSize(ham_inp%redcoord,nt,ham_inp%max_no_shells)
      call read_exchange_reduceCouplingMatrixSize(ham_inp%jc,nt,ham_inp%max_no_shells,nchmax)
      if (ham_inp%exc_inter=='Y') call read_exchange_reduceCouplingMatrixSize(ham_inp%jcD,nt,ham_inp%max_no_shells,nchmax)
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

      call read_exchange_getMaxNoShells(no_shells,flines)
      call allocate_hamiltonianinput(ham_inp,no_shells,1)
      ham_inp%redcoord = 0.0_dblprec
      ham_inp%jc_tens  = 0.0_dblprec
      ham_inp%nn       = 0

      ! Read exchange vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines
         if(do_ralloy==0) then
            if(tensor_format) then
               read(ifileno,*) isite,jsite,r_tmp(1:3),j_tmp
               j_tmp = transpose(j_tmp)
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
         jtype=1 !atype_inp(jsite)
         call read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,ham_inp%nn(itype)
            norm=(r_red(1)-ham_inp%redcoord(itype,ishell,1))**2+ &
                 (r_red(2)-ham_inp%redcoord(itype,ishell,2))**2+ &
                 (r_red(3)-ham_inp%redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(tensor_format) then
                  ham_inp%jc_tens(1:3,1:3,itype,ishell,ichem,jtype)=j_tmp(1:3,1:3)
               else
                  ham_inp%jc_tens(1:3,1:3,itype,ishell,ichem,jtype)= &
                  reshape((/j_tmpSingle,0.0_dblprec,0.0_dblprec, 0.0_dblprec,j_tmpSingle,0.0_dblprec,&
                  0.0_dblprec,0.0_dblprec,j_tmpSingle/),(/3,3/))
               endif
            end if
         end do
         if (unique.or.ham_inp%map_multiple) then
            ham_inp%nn(itype)=ham_inp%nn(itype)+1
            ham_inp%redcoord(itype,ham_inp%nn(itype),1:3)=r_red(1:3)
            if(tensor_format) then
               ham_inp%jc_tens(1:3,1:3,itype,ishell,ichem,jtype)=j_tmp(1:3,1:3)
            else
               ham_inp%jc_tens(1:3,1:3,itype,ishell,ichem,jtype)= &
               reshape( (/j_tmpSingle,0.0_dblprec,0.0_dblprec, 0.0_dblprec,j_tmpSingle,0.0_dblprec,   &
               0.0_dblprec,0.0_dblprec,j_tmpSingle/),(/3,3/))
            endif
         end if
      enddo
      close (ifileno)

      ham_inp%max_no_shells=maxval(ham_inp%NN)
      call read_exchange_reduceRedCoordMatrixSize(ham_inp%redcoord,nt,ham_inp%max_no_shells)
      call read_exchange_reduceCouplingMatrixSizeTensor(ham_inp%jc_tens,nt,ham_inp%max_no_shells,nchmax)

      close(ifileno)

   end subroutine read_exchange_tensor_base

   !--------------------------------------------------------------------------------
   !> @brief
   !> Trying to reduce the nntype matrix size
   !>
   !> @author
   !> Thomas Nystrand
   !--------------------------------------------------------------------------------
   subroutine read_exchange_reduceNNtypeMatrixSize(nntype,nt,max_no_shells)
      !
      integer, intent(in)                                         :: nt, max_no_shells
      integer, intent(inout), dimension(:,:), allocatable :: nntype

      ! locals
      integer                                       :: i_stat,i_all
      integer, dimension(:,:), allocatable  :: nntype_tmp

      allocate(nntype_tmp(nt,max_no_shells),stat=i_stat)
      call memocc(i_stat,product(shape(nntype_tmp))*kind(nntype_tmp),'nntype_tmp','read_exchange_reducenntypeMatrixSize')
      nntype_tmp = nntype(:,1:max_no_shells)

      i_all=-product(shape(nntype))*kind(nntype)
      deallocate(nntype,stat=i_stat)
      call memocc(i_stat,i_all,'nntype','read_exchange_reducenntypeMatrixSize')

      allocate(nntype(nt,max_no_shells),stat=i_stat)
      call memocc(i_stat,product(shape(nntype))*kind(nntype),'nntype','read_exchange_reducenntypeMatrixSize')

      nntype = nntype_tmp

      i_all=-product(shape(nntype_tmp))*kind(nntype_tmp)
      deallocate(nntype_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'nntype_tmp','read_exchange_reducenntypeMatrixSize')

      close (ifileno)
   end subroutine read_exchange_reduceNNtypeMatrixSize


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
   subroutine read_exchange_getMaxNoShells(no_shells,flines)
      use SystemData, only : atype

      implicit none

      integer, intent(out)                   :: no_shells,flines
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
         !itype=atype_inp(isite)
         !jtype=atype_inp(jsite)
         itype=atype(isite)
         jtype=atype(jsite)
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
         if (posfiletype=='C') then
            r_red=r_tmp
         elseif (posfiletype=='D') then
            r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
            r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
            r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
         else
            stop 'Only posfiletype = C or D is currently supported'
         endif
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

      open(ifileno, file=adjustl(ham_inp%kfile))

      i_err=0
      do while(i_err==0)
         read (ifileno,*,iostat=i_err) iat,ichem,ham_inp%anisotropytype(iat,ichem),ham_inp%anisotropy(iat,1:6,ichem)
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

      open(ifileno, file=adjustl(ham_inp%kfile))

      do m=1, na
         read (ifileno,*) iat,ham_inp%anisotropytype(iat,1),ham_inp%anisotropy(iat,1:6,1)
      enddo

      if (ham_inp%mult_axis=='Y') then
         do m=1,na
            read(ifileno,*) iat,ham_inp%anisotropytype_diff(iat,1),ham_inp%anisotropy_diff(iat,1:6,1)
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
      open(ifileno, file=ham_inp%dmfile)
      if (do_cluster=='Y') open(ifileno2, file=dmfile_clus) ! File for the interactions inside the embeded cluster

      if (do_cluster=='Y') then
         call read_exchange_getMaxNoShells(no_shells,flines)
         ham_inp%max_no_dmshells = no_shells
         allocate(ham_inp%dm_redcoord(NT,ham_inp%max_no_dmshells,3),stat=i_stat)
         call memocc(i_stat,product(shape(ham_inp%dm_redcoord))*kind(ham_inp%dm_redcoord),'dm_redcoord','read_dmdata')
         ham_inp%dm_redcoord  = 0.0_dblprec
         allocate(ham_inp%dm_inpvect(3,NT,ham_inp%max_no_dmshells,Nchmax,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ham_inp%dm_inpvect))*kind(ham_inp%dm_inpvect),'dm_inpvect','read_dmdata')
         ham_inp%dm_inpvect = 0.0_dblprec
         call read_exchange_getMaxNoShells_clus(no_shells_clus,flines_clus,ifileno2,do_ralloy)
         max_no_dmshells_clus=no_shells_clus
         allocate(dm_redcoord_clus(NT_clus,max_no_dmshells_clus,3),stat=i_stat)
         call memocc(i_stat,product(shape(dm_redcoord_clus))*kind(dm_redcoord_clus),'dm_redcoord_clus','read_dmdata')
         dm_redcoord_clus  = 0.0_dblprec
         allocate(dm_inpvect_clus(3,NT_clus,max_no_dmshells_clus,Nchmax,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(dm_inpvect_clus))*kind(dm_inpvect_clus),'dm_inpvect_clus','read_dmdata')
         dm_inpvect_clus = 0.0_dblprec
      else
         call read_exchange_getMaxNoShells(no_shells,flines)
         ham_inp%max_no_dmshells = no_shells
         allocate(ham_inp%dm_redcoord(NT,ham_inp%max_no_dmshells,3),stat=i_stat)
         call memocc(i_stat,product(shape(ham_inp%dm_redcoord))*kind(ham_inp%dm_redcoord),'dm_redcoord','read_dmdata')
         ham_inp%dm_redcoord  = 0.0_dblprec
         allocate(ham_inp%dm_inpvect(3,NT,ham_inp%max_no_dmshells,Nchmax,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ham_inp%dm_inpvect))*kind(ham_inp%dm_inpvect),'dm_inpvect','read_dmdata')
         ham_inp%dm_inpvect = 0.0_dblprec
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
         do ishell=1,ham_inp%dm_nn(itype)
            norm=(r_red(1)-ham_inp%dm_redcoord(itype,ishell,1))**2+ &
                 (r_red(2)-ham_inp%dm_redcoord(itype,ishell,2))**2+ &
                 (r_red(3)-ham_inp%dm_redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  ham_inp%dm_inpvect(1:3,itype,ishell,ichem,jtype)=dm_tmp
               else
                  ham_inp%dm_inpvect(1:3,itype,ishell,ichem,jchem)=dm_tmp
               end if
            end if
         end do
         if (unique) then
            ham_inp%dm_nn(itype)=ham_inp%dm_nn(itype)+1

            ham_inp%dm_redcoord(itype,ham_inp%dm_nn(itype),1:3)=r_red(1:3)
            if(do_ralloy==0) then
               ham_inp%dm_inpvect(1:3,itype,ham_inp%dm_nn(itype),ichem,1)=dm_tmp
            else
               ham_inp%dm_inpvect(1:3,itype,ham_inp%dm_nn(itype),ichem,jchem)=dm_tmp
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
      ham_inp%max_no_dmshells=maxval(ham_inp%dm_nn)
      call read_exchange_reduceRedCoordMatrixSize(ham_inp%dm_redcoord,nt,ham_inp%max_no_dmshells)
      call read_dmexchange_reduceCouplingMatrixSize(ham_inp%dm_inpvect,nt,ham_inp%max_no_dmshells,nchmax)

      if (do_cluster=='Y') then
         max_no_dmshells_clus=maxval(dm_nn_clus)
         call read_exchange_reduceRedCoordMatrixSize(dm_redcoord_clus,NT_clus,max_no_dmshells_clus)
         call read_dmexchange_reduceCouplingMatrixSize(dm_inpvect_clus,NT_clus,max_no_dmshells_clus,nchmax)
         close(ifileno2)
      endif
      close (ifileno)

   end subroutine read_dmdata

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_sadata
   !> @brief Read the variables for the symmetric anisotropic exchane
   !---------------------------------------------------------------------------------
   subroutine read_sadata()
      !
      implicit none
      !
      integer :: flines, isite, i_stat, jsite
      integer :: flines_clus, isite_c, jsite_c
      integer :: itype, jtype, ichem, jchem, iline, ishell
      integer :: itype_clus, jtype_clus, ichem_c, jchem_c
      logical :: unique
      real(dblprec), dimension(3) :: sa_tmp,sa_tmp_clus
      real(dblprec):: tol, norm
      real(dblprec), dimension(3) :: r_tmp, r_red, r_tmp_clus,r_red_clus
      integer :: no_shells, no_shells_clus

      integer :: ifileno2

      ! Set tolerance for neighbour shells
      tol=1.0d-5

      ifileno2=ifileno+1
      open(ifileno, file=ham_inp%safile)
      if (do_cluster=='Y') open(ifileno2, file=safile_clus) ! File for the interactions inside the embeded cluster

      if (do_cluster=='Y') then
         call read_exchange_getMaxNoShells(no_shells,flines)
         ham_inp%max_no_sashells = no_shells
         allocate(ham_inp%sa_redcoord(NT,ham_inp%max_no_sashells,3),stat=i_stat)
         call memocc(i_stat,product(shape(ham_inp%sa_redcoord))*kind(ham_inp%sa_redcoord),'sa_redcoord','read_sadata')
         ham_inp%sa_redcoord  = 0.0_dblprec
         allocate(ham_inp%sa_inpvect(3,NT,ham_inp%max_no_sashells,Nchmax,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ham_inp%sa_inpvect))*kind(ham_inp%sa_inpvect),'sa_inpvect','read_sadata')
         ham_inp%sa_inpvect = 0.0_dblprec
         call read_exchange_getMaxNoShells_clus(no_shells_clus,flines_clus,ifileno2,do_ralloy)
         max_no_sashells_clus=no_shells_clus
         allocate(sa_redcoord_clus(NT_clus,max_no_sashells_clus,3),stat=i_stat)
         call memocc(i_stat,product(shape(sa_redcoord_clus))*kind(sa_redcoord_clus),'sa_redcoord_clus','read_sadata')
         sa_redcoord_clus  = 0.0_dblprec
         allocate(sa_inpvect_clus(3,NT_clus,max_no_sashells_clus,Nchmax,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(sa_inpvect_clus))*kind(sa_inpvect_clus),'sa_inpvect_clus','read_sadata')
         sa_inpvect_clus = 0.0_dblprec
      else
         call read_exchange_getMaxNoShells(no_shells,flines)
         ham_inp%max_no_sashells = no_shells
         allocate(ham_inp%sa_redcoord(NT,ham_inp%max_no_sashells,3),stat=i_stat)
         call memocc(i_stat,product(shape(ham_inp%sa_redcoord))*kind(ham_inp%sa_redcoord),'sa_redcoord','read_sadata')
         ham_inp%sa_redcoord  = 0.0_dblprec
         allocate(ham_inp%sa_inpvect(3,NT,ham_inp%max_no_sashells,Nchmax,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(ham_inp%sa_inpvect))*kind(ham_inp%sa_inpvect),'sa_inpvect','read_sadata')
         ham_inp%sa_inpvect = 0.0_dblprec
      endif

      ! Read exchange vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp, sa_tmp
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, sa_tmp
         end if

         ! Find type of site
         itype=atype_inp(isite)
         jtype=1

         call read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,ham_inp%sa_nn(itype)
            norm=(r_red(1)-ham_inp%sa_redcoord(itype,ishell,1))**2+ &
                 (r_red(2)-ham_inp%sa_redcoord(itype,ishell,2))**2+ &
                 (r_red(3)-ham_inp%sa_redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  ham_inp%sa_inpvect(1:3,itype,ishell,ichem,jtype)=sa_tmp
               else
                  ham_inp%sa_inpvect(1:3,itype,ishell,ichem,jchem)=sa_tmp
               end if
            end if
         end do
         if (unique) then
            ham_inp%sa_nn(itype)=ham_inp%sa_nn(itype)+1

            ham_inp%sa_redcoord(itype,ham_inp%sa_nn(itype),1:3)=r_red(1:3)
            if(do_ralloy==0) then
               ham_inp%sa_inpvect(1:3,itype,ham_inp%sa_nn(itype),ichem,1)=sa_tmp
            else
               ham_inp%sa_inpvect(1:3,itype,ham_inp%sa_nn(itype),ichem,jchem)=sa_tmp
            end if
         end if
      enddo

      if (do_cluster=='Y') then

         do iline=1, flines_clus

            ! Read indices and coordinates
            if(do_ralloy==0) then
               read (ifileno2,*) isite_c, jsite_C, r_tmp_clus, sa_tmp_clus
               ichem_c=1
               jchem_c=1
            else
               read (ifileno2,*) isite_c, jsite_c, ichem_c, jchem_c, r_tmp_clus, sa_tmp_clus
            end if

            ! Find type of site
            itype_clus=atype_inp_clus(isite_c)
            jtype_clus=1

            call read_exchange_getNeighVec_clus(r_red_clus,r_tmp_clus,isite_c,jsite_c,&
               maptype,posfiletype)

            ! Loop through earlier vectors to find equivalent shells
            unique=.true.
            do ishell=1,sa_nn_clus(itype_clus)
               norm=(r_red_clus(1)-sa_redcoord_clus(itype_clus,ishell,1))**2+ &
               (r_red_clus(2)-sa_redcoord_clus(itype_clus,ishell,2))**2+ &
               (r_red_clus(3)-sa_redcoord_clus(itype_clus,ishell,3))**2
               if(norm<tol) then
                  unique=.false.
                  if(do_ralloy==0) then
                     sa_inpvect_clus(1:3,itype_clus,ishell,ichem_c,jtype_clus)=sa_tmp_clus
                  else
                     sa_inpvect_clus(1:3,itype_clus,ishell,ichem_c,jchem_c)=sa_tmp_clus
                  end if
               end if
            end do
            if (unique) then
               sa_nn_clus(itype_clus)=sa_nn_clus(itype_clus)+1

               sa_redcoord_clus(itype_clus,sa_nn_clus(itype_clus),1:3)=r_red_clus(1:3)
               if(do_ralloy==0) then
                  sa_inpvect_clus(1:3,itype_clus,sa_nn_clus(itype_clus),ichem_c,1)=sa_tmp_clus
               else
                  sa_inpvect_clus(1:3,itype_clus,sa_nn_clus(itype_clus),ichem_c,jchem_c)=sa_tmp_clus
               end if
            end if
         enddo
      endif
      ham_inp%max_no_sashells=maxval(ham_inp%sa_nn)
      call read_exchange_reduceRedCoordMatrixSize(ham_inp%sa_redcoord,nt,ham_inp%max_no_sashells)
      call read_dmexchange_reduceCouplingMatrixSize(ham_inp%sa_inpvect,nt,ham_inp%max_no_sashells,nchmax)

      if (do_cluster=='Y') then
         max_no_sashells_clus=maxval(sa_nn_clus)
         call read_exchange_reduceRedCoordMatrixSize(sa_redcoord_clus,NT_clus,max_no_sashells_clus)
         call read_dmexchange_reduceCouplingMatrixSize(sa_inpvect_clus,NT_clus,max_no_sashells_clus,nchmax)
         close(ifileno2)
      endif
      close (ifileno)

   end subroutine read_sadata


   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_chirdata
   !> @brief Read the variables for three-site scalar chirality interaction
   !---------------------------------------------------------------------------------
   !      ! Variables for CHIR exchange
   ! integer :: nn_chir_tot                                     !< Calculated number of neighbours with chir interactions
   ! integer ::  max_no_chirneigh                               !< Calculated maximum of neighbours for PD exchange
   ! integer, dimension(:), allocatable :: chirlistsize         !< Size of neighbour list for chir
   ! integer, dimension(:,:), allocatable :: chirlist           !< List of neighbours for chir
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
      open(ifileno, file=ham_inp%chirfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0
      nskip=0

      call read_exchange_getMaxNoShells(no_shells,flines)
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

      ham_inp%max_no_chirshells=no_shells
      allocate(ham_inp%chir_redcoord(NT,ham_inp%max_no_chirshells,3,2),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%chir_redcoord))*kind(ham_inp%chir_redcoord),'chir_redcoord','read_chirdata')
      ham_inp%chir_redcoord=0.0_dblprec
      allocate(ham_inp%chir_inpval(NT,ham_inp%max_no_chirshells,NT,NT),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%chir_inpval))*kind(ham_inp%chir_inpval),'chir_inpval','read_chirdata')
      ham_inp%chir_inpval=0.0_dblprec
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
         do ishell=1,ham_inp%chir_nn(itype)
            norm=(r_red(1)-ham_inp%chir_redcoord(itype,ishell,1,1))**2+ &
                 (r_red(2)-ham_inp%chir_redcoord(itype,ishell,2,1))**2+ &
                 (r_red(3)-ham_inp%chir_redcoord(itype,ishell,3,1))**2+ &
                 (r_red(4)-ham_inp%chir_redcoord(itype,ishell,1,2))**2+ &
                 (r_red(5)-ham_inp%chir_redcoord(itype,ishell,2,2))**2+ &
                 (r_red(6)-ham_inp%chir_redcoord(itype,ishell,3,2))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  ham_inp%chir_inpval(itype,ishell,ichem,jtype)=chir_tmp
               else
                  ham_inp%chir_inpval(itype,ishell,ichem,jchem)=chir_tmp
               end if
            end if
         end do
         if (unique) then
            ham_inp%chir_nn(itype)=ham_inp%chir_nn(itype)+1
            ham_inp%chir_redcoord(itype,ham_inp%chir_nn(itype),1:3,1)=r_red(1:3)
            ham_inp%chir_redcoord(itype,ham_inp%chir_nn(itype),1:3,2)=r_red(4:6)
            if(do_ralloy==0) then
               ham_inp%chir_inpval(itype,ham_inp%chir_nn(itype),ichem,jtype)=chir_tmp
            else
               ham_inp%chir_inpval(itype,ham_inp%chir_nn(itype),ichem,jchem)=chir_tmp
            end if
         end if
      enddo
      close (ifileno)


   end subroutine read_chirdata


   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_fourxdata
   !> @brief Read the variables for four-site scalar interaction
   !---------------------------------------------------------------------------------
   ! integer :: nn_fourx_tot                                     !< Calculated number of neighbours with PD interactions
   ! integer ::  max_no_fourxneigh                               !< Calculated maximum of neighbours for PD exchange
   ! integer, dimension(:), allocatable :: fourxlistsize         !< Size of neighbour list for PD
   ! integer, dimension(:,:), allocatable :: fourxlist           !< List of neighbours for PD
   ! real(dblprec), dimension(:,:,:), allocatable :: fourx_vect  !< Pseudo-Dipolar exchange vector

   subroutine read_fourxdata()
      !
      !
      implicit none
      !
      integer :: flines,mtype,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec) :: fourx_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(6) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp
      integer :: i, nskip

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=ham_inp%fourxfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0
      nskip=0

      call read_exchange_getMaxNoShells(no_shells,flines)
      !max_no_fourxshells = no_shells
      !call allocate_latthamiltonianinput(no_shells,1)

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_fourxdata')
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
      call memocc(i_stat,i_all,'NN_tmp','read_fourxdata')

      ham_inp%max_no_fourxshells=no_shells
      allocate(ham_inp%fourx_redcoord(NT,ham_inp%max_no_fourxshells,3,2),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%fourx_redcoord))*kind(ham_inp%fourx_redcoord),'fourx_redcoord','read_fourxdata')
      ham_inp%fourx_redcoord=0.0_dblprec
      allocate(ham_inp%fourx_inpval(NT,ham_inp%max_no_fourxshells,NT,NT),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%fourx_inpval))*kind(ham_inp%fourx_inpval),'fourx_inpval','read_fourxdata')
      ham_inp%fourx_inpval=0.0_dblprec
      ! Read force coupling vectors
      ! Isite, Jsite, Ichem, Jchem
      do iline=1, flines

         ! Read indices and coordinates
         ! A block of four lines is used for each input fourx-tensor
         if(do_ralloy==0) then
            read (ifileno,*) isite, jsite, r_tmp, fourx_tmp
            ichem=1
            jchem=1
         else
            read (ifileno,*) isite, jsite, ichem, jchem, r_tmp, fourx_tmp
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
         do ishell=1,ham_inp%fourx_nn(itype)
            norm=(r_red(1)-ham_inp%fourx_redcoord(itype,ishell,1,1))**2+ &
                 (r_red(2)-ham_inp%fourx_redcoord(itype,ishell,2,1))**2+ &
                 (r_red(3)-ham_inp%fourx_redcoord(itype,ishell,3,1))**2+ &
                 (r_red(4)-ham_inp%fourx_redcoord(itype,ishell,1,2))**2+ &
                 (r_red(5)-ham_inp%fourx_redcoord(itype,ishell,2,2))**2+ &
                 (r_red(6)-ham_inp%fourx_redcoord(itype,ishell,3,2))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  ham_inp%fourx_inpval(itype,ishell,ichem,jtype)=fourx_tmp
               else
                  ham_inp%fourx_inpval(itype,ishell,ichem,jchem)=fourx_tmp
               end if
            end if
         end do
         if (unique) then
            ham_inp%fourx_nn(itype)=ham_inp%fourx_nn(itype)+1
            ham_inp%fourx_redcoord(itype,ham_inp%fourx_nn(itype),1:3,1)=r_red(1:3)
            ham_inp%fourx_redcoord(itype,ham_inp%fourx_nn(itype),1:3,2)=r_red(4:6)
            if(do_ralloy==0) then
               ham_inp%fourx_inpval(itype,ham_inp%fourx_nn(itype),ichem,jtype)=fourx_tmp
            else
               ham_inp%fourx_inpval(itype,ham_inp%fourx_nn(itype),ichem,jchem)=fourx_tmp
            end if
         end if
      enddo
      close (ifileno)


   end subroutine read_fourxdata

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_pddata
   !> @brief Read the variables for the pseudo-dipolar interaction
   !---------------------------------------------------------------------------------
   subroutine read_pddata()
      !
      implicit none
      !
      integer :: flines,isite,i_stat,jsite
      integer :: itype,jtype,ichem,jchem,iline,ishell,i_all
      logical :: unique
      real(dblprec), dimension(9) :: pd_tmp
      real(dblprec):: tol, norm
      real(dblprec), dimension(3) :: r_tmp, r_red
      integer :: no_shells
      integer, dimension(:,:,:), allocatable :: nn_tmp

      ! Set tolerance for neighbour shells
      tol=1.0d-5
      ! Open input file
      open(ifileno, file=ham_inp%pdfile)
      ! Check if input file is for random alloy
      flines=0
     ! mtype=0

      allocate(nn_tmp(nt,max(nt,nchmax),nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp))*kind(nn_tmp),'NN_tmp','read_pddata')
      nn_tmp=0

      ! Pre-read file to get max no. exchange shells and no. lines
      do
         if(do_ralloy==0) then
            read(ifileno,*,end=200)  isite, jsite
            ichem=1
            jchem=1
          !print*, isite, jsite
         else
            read(ifileno,*,end=200)  isite,jsite,ichem, jchem
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
      
      !print*,'Pre read done'

      no_shells=0
      do itype=1,nt
         no_shells=max(sum(nn_tmp(itype,:,:)),no_shells)
      end do

      i_all=-product(shape(nn_tmp))*kind(nn_tmp)
      deallocate(nn_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'NN_tmp','read_pddata')

      ham_inp%max_no_pdshells=no_shells
      allocate(ham_inp%pd_redcoord(NT,ham_inp%max_no_pdshells,3),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%pd_redcoord))*kind(ham_inp%pd_redcoord),'pd_redcoord','read_pddata')
      ham_inp%pd_redcoord=0.0_dblprec
      allocate(ham_inp%pd_inpvect(9,NT,ham_inp%max_no_pdshells,Nchmax,Nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%pd_inpvect))*kind(ham_inp%pd_inpvect),'pd_inpvect','read_pddata')
      ham_inp%pd_inpvect=0.0_dblprec


        ! call read_exchange_getMaxNoShells(no_shells,flines,ham_inp%dmfile)
        ! ham_inp%max_no_dmshells = no_shells
         !allocate(ham_inp%dm_redcoord(NT,ham_inp%max_no_dmshells,3),stat=i_stat)
        ! call memocc(i_stat,product(shape(ham_inp%dm_redcoord))*kind(ham_inp%dm_redcoord),'dm_redcoord','read_dmdata')
        ! ham_inp%dm_redcoord  = 0.0_dblprec
        ! allocate(ham_inp%dm_inpvect(3,NT,ham_inp%max_no_dmshells,Nchmax,Nchmax),stat=i_stat)
       !  call memocc(i_stat,product(shape(ham_inp%dm_inpvect))*kind(ham_inp%dm_inpvect),'dm_inpvect','read_dmdata')
       !  ham_inp%dm_inpvect = 0.0_dblprec

     
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

           call read_exchange_getNeighVec(r_red,r_tmp,isite,jsite)

         !if(maptype==2) then
            ! Calculate proper neighbour vector (from "bgfm")
           ! r_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(1)+C2(1)*r_tmp(2)+C3(1)*r_tmp(3)
           ! r_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(1)+C2(2)*r_tmp(2)+C3(2)*r_tmp(3)
           ! r_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(1)+C2(3)*r_tmp(2)+C3(3)*r_tmp(3)
         !else
            ! Calculates neighbour vectors from direct coordinates or Cartesian
            ! coordinates, corresponding to how the atomic positions are entered
            !if (posfiletype=='C') then
            !   r_red=r_tmp
            !elseif (posfiletype=='D') then
            !   r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
            !   r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
            !   r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
            !else
            !   stop 'Only posfiletype= C or D is currently supported'
            !endif
         !end if

         ! Loop through earlier vectors to find equivalent shells
         unique=.true.
         do ishell=1,ham_inp%pd_nn(itype)
            norm=(r_red(1)-ham_inp%pd_redcoord(itype,ishell,1))**2+ &
                 (r_red(2)-ham_inp%pd_redcoord(itype,ishell,2))**2+ &
                 (r_red(3)-ham_inp%pd_redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  ham_inp%pd_inpvect(1:9,itype,ishell,ichem,1)=pd_tmp
               else
                  ham_inp%pd_inpvect(1:9,itype,ishell,ichem,jchem)=pd_tmp
               end if
            end if
         end do
         if (unique) then
            ham_inp%pd_nn(itype)=ham_inp%pd_nn(itype)+1
            ham_inp%pd_redcoord(itype,ham_inp%pd_nn(itype),1:3)=r_red(1:3)
            if(do_ralloy==0) then
               ham_inp%pd_inpvect(1:9,itype,ham_inp%pd_nn(itype),ichem,1)=pd_tmp
            else
               ham_inp%pd_inpvect(1:9,itype,ham_inp%pd_nn(itype),ichem,jchem)=pd_tmp
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
      open(ifileno, file=ham_inp%biqdmfile)
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

      ham_inp%max_no_biqdmshells=no_shells
      allocate(ham_inp%biqdm_redcoord(NT,ham_inp%max_no_biqdmshells,3),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%biqdm_redcoord))*kind(ham_inp%biqdm_redcoord),'biqdm_redcoord','read_biqdmdata')
      ham_inp%biqdm_redcoord=0.0_dblprec
      allocate(ham_inp%biqdm_inpvect(1,NT,ham_inp%max_no_biqdmshells,Nchmax,Nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%biqdm_inpvect))*kind(ham_inp%biqdm_inpvect),'biqdm_inpvect','read_biqdmdata')
      ham_inp%biqdm_inpvect=0.0_dblprec
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
         do ishell=1,ham_inp%biqdm_nn(itype)
            norm=(r_red(1)-ham_inp%biqdm_redcoord(itype,ishell,1))**2+ &
               (r_red(2)-ham_inp%biqdm_redcoord(itype,ishell,2))**2+ &
               (r_red(3)-ham_inp%biqdm_redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  ham_inp%biqdm_inpvect(1,itype,ishell,ichem,jtype)=biqdm_tmp
               else
                  ham_inp%biqdm_inpvect(1,itype,ishell,ichem,jchem)=biqdm_tmp
               end if
            end if
         end do
         if (unique) then
            ham_inp%biqdm_nn(itype)=ham_inp%biqdm_nn(itype)+1
            ham_inp%biqdm_redcoord(itype,ham_inp%biqdm_nn(itype),1:3)=r_red(1:3)
            if(do_ralloy==0) then
               ham_inp%biqdm_inpvect(1,itype,ham_inp%biqdm_nn(itype),ichem,jtype)=biqdm_tmp
            else
               ham_inp%biqdm_inpvect(1,itype,ham_inp%biqdm_nn(itype),ichem,jchem)=biqdm_tmp
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
      open(ifileno, file=ham_inp%bqfile)
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

      ham_inp%max_no_bqshells=no_shells
      allocate(ham_inp%bq_redcoord(NT,ham_inp%max_no_bqshells,3),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%bq_redcoord))*kind(ham_inp%bq_redcoord),'bq_redcoord','read_bqdata')
      ham_inp%bq_redcoord=0.0_dblprec
      allocate(ham_inp%jc_bq(NT,ham_inp%max_no_bqshells,Nchmax,Nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(ham_inp%jc_bq))*kind(ham_inp%jc_bq),'bq_inpvect','read_bqdata')
      ham_inp%jc_bq=0.0_dblprec
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
         do ishell=1,ham_inp%bq_nn(itype)
            norm=(r_red(1)-ham_inp%bq_redcoord(itype,ishell,1))**2+ &
               (r_red(2)-ham_inp%bq_redcoord(itype,ishell,2))**2+ &
               (r_red(3)-ham_inp%bq_redcoord(itype,ishell,3))**2
            if(norm<tol) then
               unique=.false.
               if(do_ralloy==0) then
                  ham_inp%jc_bq(itype,ishell,ichem,jtype)=jbq_tmp
               else
                  ham_inp%jc_bq(itype,ishell,ichem,jchem)=jbq_tmp
               end if
            end if
         end do
         if (unique) then
            ham_inp%bq_nn(itype)=ham_inp%bq_nn(itype)+1
            ham_inp%bq_redcoord(itype,ham_inp%bq_nn(itype),1:3)=r_red(1:3)
            if(do_ralloy==0) then
               ham_inp%jc_bq(itype,ham_inp%bq_nn(itype),ichem,jtype)=jbq_tmp
            else
               ham_inp%jc_bq(itype,ham_inp%bq_nn(itype),ichem,jchem)=jbq_tmp
            end if
         end if
      enddo
      close (ifileno)

   end subroutine read_bqdata
   
   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_ringdata
   !> @brief Read the variables for four-spin ring interaction
   !---------------------------------------------------------------------------------
    subroutine read_ringdata()
    !
    implicit none
    !
    integer :: flines, iline, ishell, i_stat, isite, jsite, ksite, lsite
    logical :: unique
    real(dblprec) :: jring_tmp
    real(dblprec) :: tol,normij,normik,normil
    real(dblprec), dimension(3) :: rij_tmp,rik_tmp,ril_tmp,rij_red,rik_red,ril_red

    ! Set tolerance for neighbour shells
    tol=1.0d-5
    
    !Proceed only in case if system is not random alloy
    if(do_ralloy==1) then
       write (*,*) 'Ring exchange is not supported for random alloys'
       stop
    else
       ! Open input file
       open(ifileno, file=ham_inp%ringfile)
       flines=0

       do
          ! Pre-read file to get max no. exchange shells and no. lines 
          read(ifileno,*,end=200)
          flines=flines+1
       end do 
    end if     
    200 continue

    ham_inp%max_no_ringshells=flines

    rewind(ifileno)

    ! At the moment support only for one type of magnetic atoms in system
    if (nt==1) then

       allocate(ham_inp%ring_redcoord_ij(NT,ham_inp%max_no_ringshells,3),stat=i_stat)
       call memocc(i_stat,product(shape(ham_inp%ring_redcoord_ij))*kind(ham_inp%ring_redcoord_ij),&
          'ring_redcoord_ij','read_ringdata')
       ham_inp%ring_redcoord_ij=0.0_dblprec

       allocate(ham_inp%ring_redcoord_ik(NT,ham_inp%max_no_ringshells,3),stat=i_stat)
       call memocc(i_stat,product(shape(ham_inp%ring_redcoord_ik))*kind(ham_inp%ring_redcoord_ik),&
          'ring_redcoord_ik','read_ringdata')
       ham_inp%ring_redcoord_ik=0.0_dblprec

       allocate(ham_inp%ring_redcoord_il(NT,ham_inp%max_no_ringshells,3),stat=i_stat)
       call memocc(i_stat,product(shape(ham_inp%ring_redcoord_il))*kind(ham_inp%ring_redcoord_il),&
          'ring_redcoord_il','read_ringdata')
       ham_inp%ring_redcoord_il=0.0_dblprec

    else

       write (*,*) 'Ring exchange is currently supported only for atoms of the same type'
       stop

    end if

    allocate(ham_inp%jc_ring(ham_inp%max_no_ringshells),stat=i_stat)
    call memocc(i_stat,product(shape(ham_inp%jc_ring))*kind(ham_inp%jc_ring),'ring_inpvect','read_ringdata')
    ham_inp%jc_ring=0.0_dblprec


    ! Read exchange vectors
    do iline=1, flines
       read (ifileno,*) isite,jsite,ksite,lsite,rij_tmp,rik_tmp,ril_tmp,jring_tmp

       if(maptype==2) then
          rij_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*rij_tmp(1)+C2(1)*rij_tmp(2)+C3(1)*rij_tmp(3)
          rij_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*rij_tmp(1)+C2(2)*rij_tmp(2)+C3(2)*rij_tmp(3)
          rij_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*rij_tmp(1)+C2(3)*rij_tmp(2)+C3(3)*rij_tmp(3)

          rik_red(1)=Bas(1,ksite)-Bas(1,isite)+C1(1)*rik_tmp(1)+C2(1)*rik_tmp(2)+C3(1)*rik_tmp(3)
          rik_red(2)=Bas(2,ksite)-Bas(2,isite)+C1(2)*rik_tmp(1)+C2(2)*rik_tmp(2)+C3(2)*rik_tmp(3)
          rik_red(3)=Bas(3,ksite)-Bas(3,isite)+C1(3)*rik_tmp(1)+C2(3)*rik_tmp(2)+C3(3)*rik_tmp(3)

          ril_red(1)=Bas(1,lsite)-Bas(1,isite)+C1(1)*ril_tmp(1)+C2(1)*ril_tmp(2)+C3(1)*ril_tmp(3)
          ril_red(2)=Bas(2,lsite)-Bas(2,isite)+C1(2)*ril_tmp(1)+C2(2)*ril_tmp(2)+C3(2)*ril_tmp(3)
          ril_red(3)=Bas(3,lsite)-Bas(3,isite)+C1(3)*ril_tmp(1)+C2(3)*ril_tmp(2)+C3(3)*ril_tmp(3)
       else                        
          ! Calculates neighbour vectors from direct coordinates or Cartesian
          ! coordinates, corresponding to how the atomic positions are entered
          if (posfiletype=='C') then
             rij_red=rij_tmp
             rik_red=rik_tmp
             ril_red=ril_tmp
          else if (posfiletype=='D') then
             rij_red(1)=rij_tmp(1)*C1(1)+rij_tmp(2)*C2(1)+rij_tmp(3)*C3(1)
             rij_red(2)=rij_tmp(1)*C1(2)+rij_tmp(2)*C2(2)+rij_tmp(3)*C3(2)
             rij_red(3)=rij_tmp(1)*C1(3)+rij_tmp(2)*C2(3)+rij_tmp(3)*C3(3)
             !
             rik_red(1)=rik_tmp(1)*C1(1)+rik_tmp(2)*C2(1)+rik_tmp(3)*C3(1)
             rik_red(2)=rik_tmp(1)*C1(2)+rik_tmp(2)*C2(2)+rik_tmp(3)*C3(2)
             rik_red(3)=rik_tmp(1)*C1(3)+rik_tmp(2)*C2(3)+rik_tmp(3)*C3(3)
             !
             ril_red(1)=ril_tmp(1)*C1(1)+ril_tmp(2)*C2(1)+ril_tmp(3)*C3(1)
             ril_red(2)=ril_tmp(1)*C1(2)+ril_tmp(2)*C2(2)+ril_tmp(3)*C3(2)
             ril_red(3)=ril_tmp(1)*C1(3)+ril_tmp(2)*C2(3)+ril_tmp(3)*C3(3) 
          end if
       end if
       ! Loop through earlier vectors to find equivalent shells
       unique=.true.
       do ishell=1,ham_inp%ring_nn(1)
          normij=(rij_red(1)-ham_inp%ring_redcoord_ij(1,ishell,1))**2+ &
             (rij_red(2)-ham_inp%ring_redcoord_ij(1,ishell,2))**2+ &
             (rij_red(3)-ham_inp%ring_redcoord_ij(1,ishell,3))**2
          !
          normik=(rik_red(1)-ham_inp%ring_redcoord_ik(1,ishell,1))**2+ &
             (rik_red(2)-ham_inp%ring_redcoord_ik(1,ishell,2))**2+ &
             (rik_red(3)-ham_inp%ring_redcoord_ik(1,ishell,3))**2
          !
          normil=(ril_red(1)-ham_inp%ring_redcoord_il(1,ishell,1))**2+ &
             (ril_red(2)-ham_inp%ring_redcoord_il(1,ishell,2))**2+ &
             (ril_red(3)-ham_inp%ring_redcoord_il(1,ishell,3))**2

          if((normij<tol).and.(normik<tol).and.(normil<tol)) then
             unique=.false.
             ham_inp%jc_ring(ishell)=jring_tmp    
          end if
       end do
       if (unique) then
          ham_inp%ring_nn(1)=ham_inp%ring_nn(1)+1
          ham_inp%ring_redcoord_ij(1,iline,1:3)=rij_red(1:3)
          ham_inp%ring_redcoord_ik(1,iline,1:3)=rik_red(1:3)        
          ham_inp%ring_redcoord_il(1,iline,1:3)=ril_red(1:3)
          ham_inp%jc_ring(iline)=jring_tmp
       end if
    end do
    close (ifileno)

 end subroutine read_ringdata


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

      call read_exchange_getMaxNoShells(no_shells,flines)
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

end module InputHandler_ext
