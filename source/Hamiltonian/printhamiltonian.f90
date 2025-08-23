!-------------------------------------------------------------------------------
! SUBROUTINE: PrintHamiltonian
!> @brief Routines for printing information about the Hamiltonian
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module PrintHamiltonian
   use Parameters
   use Profiling
   use Constants

   implicit none
   public

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: prn_exchange
   !> @brief Print strength of exchange couplings
   !> @details New version which replaces the older printing routine, it merges the
   !> information of the struct and struct1 files, into a more self-contained structure
   !> @note Index is wrong for ncoup in case of `do_recuce` but it should still be safe
   !> @note Jonathan Chico: modified the routine so that it writes the couplings in mRy
   !> for easier comparison with the jfile
   !----------------------------------------------------------------------------
   subroutine prn_exchange(NA,mdim,Natom,Nchmax,do_ralloy,Natom_full,max_no_neigh,  &
      simid,anumb,atype,nlistsize,asite_ch,achem_ch,nlist,coord,ammom_inp,ncoup,aham)

      use Math_functions, only : f_wrap_coord_diff
      use InputData, only : do_storeham

      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: mdim            !< dimension of the exchange coupling matrix (1=scalar or 9=3x3)
      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: Nchmax          !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy       !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full      !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh    !< Calculated maximum of neighbours for exchange
      character(len=8),intent(in) :: simid   !< Name of simulation
      integer, dimension(Natom), intent(in) :: anumb    !< Atom number in cell
      integer, dimension(Natom), intent(in) :: atype     !< Type of atom
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom_full), intent(in)   :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in)   :: achem_ch !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_neigh, Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax), intent(in) :: ammom_inp   !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(mdim,max_no_neigh, Natom), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(Natom), optional, intent(in) :: aham !< Hamiltonian look-up table

      !.. Local variables
      integer :: iatom,jatom,ineigh
      integer, dimension(:), allocatable :: alist
      character(len=20) :: filn
      real(dblprec) :: fc2_inv,tol
      real(dblprec) :: tmp_rij_norm
      real(dblprec), dimension(3) :: tmp_rij
      real(dblprec), dimension(mdim) :: tmp_coup
      tol=1e-5
      tmp_rij=0.0_dblprec
      tmp_rij_norm=0.0_dblprec
      tmp_coup=0.0_dblprec

      fc2_inv=mub/(mry*2.0_dblprec)
      !.. Executable statements
      write (filn,'(''struct.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      allocate(alist(Natom))
      if(present(aham)) then
         do iatom=1,Natom
            alist(iatom)=aham(iatom)
         end do
      else
         do iatom=1,Natom
            alist(iatom)=iatom
         end do
      end if
      write(ofileno,'(a)')"#######################################################"
      write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
      write(ofileno,'(a,1x,i8)')"# Maximum num of neighbours: ", max_no_neigh
      write(ofileno,'(a)')"#######################################################"
      if (mdim==1) then
         write(ofileno,10000) "#  iatom", "jatom",  "itype", "jtype", "r_{ij}^x",   &
         "r_{ij}^y","r_{ij}^z",  "J_{ij}", "|r_{ij}|"
      else
         write(ofileno,10010) "#iatom", "jatom",  "itype", "jtype", "r_{ij}^x",     &
         "r_{ij}^y","r_{ij}^z", "J_{ij}^{xx}","J_{ij}^{xy}","J_{ij}^{xz}",          &
         "J_{ij}^{yx}","J_{ij}^{yy}","J_{ij}^{yz}","J_{ij}^{zx}","Jij_zy","Jij_zz", &
         "|rij|"
      endif
      ! print neighbor list - after sort
      do iatom=1,Natom
         do ineigh=1,nlistsize(alist(iatom))
            jatom=nlist(ineigh,iatom)
            call f_wrap_coord_diff(Natom,coord,iatom,jatom,tmp_rij)
            tmp_rij_norm=norm2(tmp_rij)
            !if (tmp_rij_norm<tol) then
            !   write(*,'(1x,a,i6)') 'WARNING: Jii entry in neighbour map for atom',iatom
            !   write (ofileno,'(1x,a)') 'WARNING: Jii entry in neighbour map'
            !         write (ofileno,10001) iatom,jatom,atype(iatom),atype(jatom),   &
            !         tmp_rij(1:3),tmp_coup,tmp_rij_norm
            !else
               ! If scalar interaction
               if (mdim==1) then
                  if (do_ralloy==0) then
                     ! Calculate the coupling so that it has the same units than the jfile
                     tmp_coup=ncoup(1:mdim,ineigh,alist(iatom))*fc2_inv*            &
                     !(ammom_inp(anumb(iatom),1)*ammom_inp(anumb(jatom),1))
                     abs(ammom_inp(anumb(iatom),1)*ammom_inp(anumb(jatom),1))
                     ! Print the data
                     write (ofileno,10001) iatom,jatom,atype(iatom),atype(jatom),   &
                     tmp_rij(1:3),tmp_coup,tmp_rij_norm
                  else
                     ! Calculate the coupling so that it has the same units than the jfile
                     tmp_coup=ncoup(1:mdim,ineigh,alist(iatom))*fc2_inv*            &
                     !(ammom_inp(asite_ch(iatom),achem_ch(iatom))*ammom_inp(asite_ch(jatom),achem_ch(jatom)))
                     abs(ammom_inp(asite_ch(iatom),achem_ch(iatom))*ammom_inp(asite_ch(jatom),achem_ch(jatom)))
                     ! Print the data
                     write (ofileno,10002) iatom,jatom,atype(iatom),atype(jatom),   &
                     achem_ch(iatom),achem_ch(jatom),tmp_rij(1:3),tmp_coup,         &
                     tmp_rij_norm
                  endif
               ! If tensor
               else
                  if (do_ralloy==0) then
                     ! Calculate the coupling so that it has the same units than the jfile
                     tmp_coup=ncoup(1:mdim,ineigh,alist(iatom))*fc2_inv*            &
                     !(ammom_inp(anumb(iatom),1)*ammom_inp(anumb(jatom),1))
                     abs(ammom_inp(anumb(iatom),1)*ammom_inp(anumb(jatom),1))
                     ! Print the data
                     write (ofileno,10003) iatom,jatom,atype(iatom),atype(jatom),   &
                     tmp_rij(1:3),tmp_coup,tmp_rij_norm
                  else
                     ! Calculate the coupling so that it has the same units than the jfile
                     tmp_coup=ncoup(1:mdim,ineigh,alist(iatom))*fc2_inv*            &
                     !(ammom_inp(asite_ch(iatom),achem_ch(iatom))*ammom_inp(asite_ch(jatom),achem_ch(jatom)))
                     abs(ammom_inp(asite_ch(iatom),achem_ch(iatom))*ammom_inp(asite_ch(jatom),achem_ch(jatom)))
                     ! Print the data
                     write (ofileno,10004) iatom,jatom,atype(iatom),atype(jatom),   &
                     achem_ch(iatom),achem_ch(jatom),tmp_rij(1:3),tmp_coup,         &
                     tmp_rij_norm
                  endif
               endif
            !endif
         enddo
      end do
      !
      close(ofileno)
      !
      if(do_storeham==1) then
         write (filn,'(''hamiltonian.'',a,''.bin'')') trim(simid)
         open (ofileno,file=filn,form='unformatted')
         write(ofileno) mdim
         write(ofileno) ncoup
         close(ofileno)
      end if

      10000 format (a8,1x,a,1x,a6,1x,a6,5a16)
      10010 format (a8,1x,a,1x,a6,1x,a6,13a16)
      10001 format (i8,1X,i8,1X,i6,1x,i6,5es16.6)
      10002 format (i8,1X,i8,1X,i6,1x,i6,1x,i4,1x,i4,5es16.6)
      10003 format (i8,1X,i8,1X,i6,1x,i6,13es16.6)
      10004 format (i8,1X,i8,1X,i6,1x,i6,1x,i4,1x,i4,13es16.6)
   end subroutine prn_exchange

   !----------------------------------------------------------------------------
   ! SUBROUTINE: prn_exchange_sparse
   !> Print strength of exchange couplings in sparse matrix format
   !> The output is in Tesla
   !> @note Index is wrong for ncoup in case of `do_recuce` but it should still be safe
   !----------------------------------------------------------------------------
   subroutine prn_exchange_sparse(Natom,max_no_neigh,nlistsize,nlist,ncoup,simid,   &
      mdim)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: mdim   !< dimension of the exchange coupling matrix (1=scalar or 9=3x3)
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh, Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(mdim,max_no_neigh, Natom), intent(in) :: ncoup !< Heisenberg exchange couplings
      character(len=8),intent(in) :: simid !< Name of simulation

      !.. Local variables
      integer :: i,j
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''structsparse.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      ! print neighbor list - after sort
      do i=1,Natom
         do j=1,nlistsize(i)
            write (ofileno,10005) i, nlist(j,i), ncoup(1,j,i)
         end do
      end do
      close(ofileno)

      10005 format (2I8,2x,es16.8)

   end subroutine prn_exchange_sparse

   !----------------------------------------------------------------------------
   !> @brief  Prints the neighbour list for the induced moments
   !> @details It contains whether an atom is an iduced moment or it is fixed and
   !> its "induced" list, that is nearest neighbour list of fixed atoms (if induced)
   !> and induced (if fixed)
   !----------------------------------------------------------------------------
   subroutine prn_ind_exchange(NA,Natom,Nchmax,do_ralloy,Natom_full,                &
      max_no_neigh_ind,anumb,achtype,ind_nlistsize,ind_mom,ind_nlist,simid)

      implicit none

      integer, intent(in) :: NA           !< Number of atoms in one cell
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Nchmax       !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh_ind   !< Calculated maximum of neighbours for induced moments
      integer, dimension(Natom), intent(in)        :: anumb    !< Type of atom
      integer, dimension(Natom_full), intent(in)   :: achtype  !< Actual site of atom for dilute system
      integer, dimension(Natom), intent(in)        :: ind_nlistsize
      integer, dimension(NA,Nchmax), intent(in)    :: ind_mom  !< Indication of whether a given moment is induced/fixed (1/0) for the unit cell
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist  !< Neighbour list for induced moments
      character(len=8), intent(in) :: simid  !< Name of simulation

      ! Local variables
      integer :: i,chem
      character(len=30) :: filn

      !.. Printing the induced moments list
      write (filn,'(''indstruct.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      ! Printing the new induced moments neighbour lists
      write (ofileno,*) "Sorted data from heisge0 for induced moments"
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         if (do_ralloy==0) then
            chem=1
         else
            chem=achtype(i)
         endif
         write (ofileno,10001) i,ind_nlistsize(i),ind_mom(anumb(i),chem)
         write (ofileno,10002) ind_nlist(1:ind_nlistsize(i),i)
         if (do_ralloy==0) then
            write (ofileno,10002) ind_mom(anumb(ind_nlist(1:ind_nlistsize(i),i)),1)
         else
            write (ofileno,10002) ind_mom(anumb(ind_nlist(1:ind_nlistsize(i),i)),achtype(ind_nlist(1:ind_nlistsize(i),i)))
         endif
      end do
      close(ofileno)

      10001 format ("Atom=",i8,4x,"No neigh=",i7,2x,"Ind flag=",i7)
      10002 format ("            ",1X,5I6)

   end subroutine prn_ind_exchange

   !----------------------------------------------------------------------------
   ! SUBROUTINE: prn_dmcoup
   !> Print directions and strengths of DM couplings
   !> @details New version which replaces the older printing routine, it merges the
   !> information of the struct and struct1 files, into a more self-contained structure
   !> @note Jonathan Chico: modified the routine so that it writes the couplings in mRy
   !> for easier comparison with the jfile
   !----------------------------------------------------------------------------
   subroutine prn_dmcoup(NA,Natom,Nchmax,do_ralloy,Natom_full,max_no_dmneigh,anumb, &
      atype,dmlistsize,asite_ch,achem_ch,dmlist,coord,ammom_inp,dm_vect,simid)
      !
      use Math_functions, only : f_wrap_coord_diff
      use HamiltonianData

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: Nchmax          !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy       !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full      !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      integer, dimension(Natom), intent(in) :: anumb    !< Atom number in cell
      integer, dimension(Natom), intent(in) :: atype     !< Type of atom
      integer, dimension(Natom), intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, dimension(Natom_full), intent(in)   :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in)   :: achem_ch !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax), intent(in) :: ammom_inp   !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      character(len=8), intent(in) :: simid   !< Name of simulation
      !.. Local variables
      !.. Local variables
      integer :: iatom,jatom,ineigh
      character(len=20) :: filn
      real(dblprec) :: fc2_inv,tol
      real(dblprec) :: tmp_rij_norm
      real(dblprec), dimension(3) :: tmp_rij
      real(dblprec), dimension(3) :: tmp_coup
      integer  :: ih

      tol=1e-5
      tmp_rij=0.0_dblprec
      tmp_rij_norm=0.0_dblprec
      tmp_coup=0.0_dblprec

      fc2_inv=mub/(mry*2.0_dblprec)

      !.. Executable statements
      write (filn,'(''dmdata.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      write(ofileno,'(a)')"#######################################################"
      write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
      write(ofileno,'(a,1x,i8)')"# Maximum num of neighbours: ", max_no_dmneigh
      write(ofileno,'(a)')"#######################################################"
      write(ofileno,10000) "#  iatom", "jatom",  "itype", "jtype", "r_{ij}^x",      &
      "r_{ij}^y","r_{ij}^z", "D_{ij}^x","D_{ij}^y","D_{ij}^z", "|rij|"

      ! print neighbor list - after sort
      do iatom=1,Natom
         ih=ham%aHam(iatom)
         do ineigh=1,dmlistsize(ih)
            jatom=dmlist(ineigh,iatom)
            !tmp_rij=coord(:,jatom)-coord(:,iatom)
            call f_wrap_coord_diff(Natom,coord,iatom,jatom,tmp_rij)
            tmp_rij_norm=norm2(tmp_rij)
            if (tmp_rij_norm<tol) then
               write(*,'(1x,a,i6)') 'WARNING: Dii entry in neighbour map for atom',iatom
               write (ofileno,'(1x,a)') 'WARNING: Dii entry in neighbour map'
            else
               if (do_ralloy==0) then
                  ! Calculate the coupling so that it has the dmme units than the jfile
                  tmp_coup(1:3)=dm_vect(1:3,ineigh,iatom)*fc2_inv*                  &
                  !(ammom_inp(anumb(iatom),1)*ammom_inp(anumb(jatom),1))
                  abs(ammom_inp(anumb(iatom),1)*ammom_inp(anumb(jatom),1))
                  ! Print the data
                  write (ofileno,10001) iatom,jatom,atype(iatom),atype(jatom),      &
                  tmp_rij(1:3),tmp_coup(1:3),tmp_rij_norm
               else
                  ! Calculate the coupling so that it has the dmme units than the jfile
                  tmp_coup(1:3)=dm_vect(1:3,ineigh,iatom)*fc2_inv*                  &
                  !(ammom_inp(asite_ch(iatom),achem_ch(iatom))*ammom_inp(asite_ch(jatom),achem_ch(jatom)))
                  abs(ammom_inp(asite_ch(iatom),achem_ch(iatom))*ammom_inp(asite_ch(jatom),achem_ch(jatom)))
                  ! Print the data
                  write (ofileno,10002) iatom,jatom,atype(iatom),atype(jatom),      &
                  achem_ch(iatom),achem_ch(jatom),tmp_rij(1:3),tmp_coup(1:3),       &
                  tmp_rij_norm
               endif
            endif
         enddo
      end do
      close(ofileno)

      10000 format (a8,1x,a,1x,a6,1x,a6,7a16)
      10001 format (i8,1X,i8,1X,i5,1x,i5,7es16.4)
      10002 format (i8,1X,i8,1X,i5,1x,i5,1x,i4,1x,i4,7es16.4)

   end subroutine prn_dmcoup

   !----------------------------------------------------------------------------
   ! SUBROUTINE: prn_sacoup
   !> Print directions and strengths of SA couplings
   !> @details New version which replaces the older printing routine, it merges the
   !> information of the struct and struct1 files, into a more self-contained structure
   !> @note Jonathan Chico: modified the routine so that it writes the couplings in mRy
   !> for easier comparison with the jfile
   !----------------------------------------------------------------------------
   subroutine prn_sacoup(NA,Natom,Nchmax,do_ralloy,Natom_full,max_no_saneigh,anumb, &
      atype,salistsize,asite_ch,achem_ch,salist,coord,ammom_inp,sa_vect,simid)
      !
      use Math_functions, only : f_wrap_coord_diff

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA
      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: Nchmax          !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy       !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full      !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_saneigh !< Calculated number of neighbours with DM interactions
      integer, dimension(Natom), intent(in) :: anumb    !< Atom number in cell
      integer, dimension(Natom), intent(in) :: atype     !< Type of atom
      integer, dimension(Natom), intent(in) :: salistsize !< Size of neighbour list for DM
      integer, dimension(Natom_full), intent(in)   :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in)   :: achem_ch !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_saneigh,Natom), intent(in) :: salist   !< List of neighbours for DM
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax), intent(in) :: ammom_inp   !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(3,max_no_saneigh,Natom), intent(in) :: sa_vect !< Dzyaloshinskii-Moriya exchange vector
      character(len=8), intent(in) :: simid   !< Name of simulation
      !.. Local variables
      !.. Local variables
      integer :: iatom,jatom,ineigh
      character(len=20) :: filn
      real(dblprec) :: fc2_inv,tol
      real(dblprec) :: tmp_rij_norm
      real(dblprec), dimension(3) :: tmp_rij
      real(dblprec), dimension(3) :: tmp_coup
      tol=1e-5
      tmp_rij=0.0_dblprec
      tmp_rij_norm=0.0_dblprec
      tmp_coup=0.0_dblprec

      fc2_inv=mub/(mry*2.0_dblprec)

      !.. Executable statements
      write (filn,'(''sadata.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      write(ofileno,'(a)')"#######################################################"
      write(ofileno,'(a,1x,i8)')"# Number of atoms: ", Natom
      write(ofileno,'(a,1x,i8)')"# Maximum num of neighbours: ", max_no_saneigh
      write(ofileno,'(a)')"#######################################################"
      write(ofileno,10000) "#  iatom", "jatom",  "itype", "jtype", "r_{ij}^x",      &
      "r_{ij}^y","r_{ij}^z", "C_{ij}^x","C_{ij}^y","C_{ij}^z", "|rij|"

      ! print neighbor list - after sort
      do iatom=1,Natom
         do ineigh=1,salistsize(iatom)
            jatom=salist(ineigh,iatom)
            !tmp_rij=coord(:,jatom)-coord(:,iatom)
            call f_wrap_coord_diff(Natom,coord,iatom,jatom,tmp_rij)
            tmp_rij_norm=norm2(tmp_rij)
            if (tmp_rij_norm<tol) then
               write(*,'(1x,a,i6)') 'WARNING: Cii entry in neighbour map for atom',iatom
               write (ofileno,'(1x,a)') 'WARNING: Cii entry in neighbour map'
            else
               if (do_ralloy==0) then
                  ! Calculate the coupling so that it has the same units than the jfile
                  tmp_coup(1:3)=sa_vect(1:3,ineigh,iatom)*fc2_inv*                  &
                  !(ammom_inp(anumb(iatom),1)*ammom_inp(anumb(jatom),1))
                  abs(ammom_inp(anumb(iatom),1)*ammom_inp(anumb(jatom),1))
                  ! Print the data
                  write (ofileno,10001) iatom,jatom,atype(iatom),atype(jatom),      &
                  tmp_rij(1:3),tmp_coup(1:3),tmp_rij_norm
               else
                  ! Calculate the coupling so that it has the same units than the jfile
                  tmp_coup(1:3)=sa_vect(1:3,ineigh,iatom)*fc2_inv*                  &
                  !(ammom_inp(asite_ch(iatom),achem_ch(iatom))*ammom_inp(asite_ch(jatom),achem_ch(jatom)))
                  abs(ammom_inp(asite_ch(iatom),achem_ch(iatom))*ammom_inp(asite_ch(jatom),achem_ch(jatom)))
                  ! Print the data
                  write (ofileno,10002) iatom,jatom,atype(iatom),atype(jatom),      &
                  achem_ch(iatom),achem_ch(jatom),tmp_rij(1:3),tmp_coup(1:3),       &
                  tmp_rij_norm
               endif
            endif
         enddo
      end do
      close(ofileno)

      10000 format (a8,1x,a,1x,a6,1x,a6,7a16)
      10001 format (i8,1X,i8,1X,i5,1x,i5,7es16.4)
      10002 format (i8,1X,i8,1X,i5,1x,i5,1x,i4,1x,i4,7es16.4)

   end subroutine prn_sacoup

   !----------------------------------------------------------------------------
   !> Print directions and strengths of PD couplings
   !----------------------------------------------------------------------------
   subroutine prn_pdcoup(Natom, nn_pd_tot, pdlistsize, pdlist, pd_vect, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      integer,dimension(Natom), intent(in) :: pdlistsize !< Size of neighbour list for PD
      integer,dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
      real(dblprec),dimension(9,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=20) :: filn

      !.. Executable statements
      write (filn,'(''pddata.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      ! print neighbor list - after sort
      write (ofileno,*) "Sorted data from heisge "
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         write (ofileno,10001) i,pdlistsize(i)
         write (ofileno,10002) pdlist(1:pdlistsize(i),i)
         write (ofileno,10003) pd_vect(1:9,1:pdlistsize(i),i)*mub/mry
      end do
      close(ofileno)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (6es14.6)

   end subroutine prn_pdcoup

   !> Print directions and strengths of CHIR couplings
   subroutine prn_chircoup(Natom, max_no_chirneigh, chir_listsize, chir_list, chir_val , simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_chirneigh !< Calculated number of neighbours with chir interactions
      integer,dimension(Natom), intent(in) :: chir_listsize !< Size of neighbour list forchir 
      integer,dimension(2,max_no_chirneigh,Natom), intent(in) :: chir_list   !< List of neighbours forchir 
      real(dblprec),dimension(max_no_chirneigh,Natom), intent(in) :: chir_val  !< chir force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i,j
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''chirdata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted chir couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,chir_listsize(i)
         do j=1,chir_listsize(i)
            write (21,10004) i,chir_list(1:2,j,i),chir_val(j,i)
         end do
      end do
      close(21)

10001 format ("Atom=",i8,4x,"No neigh=",i7)
10002 format ("            ",1X,2I6)
!10002 format ("            ",1X,5I6)
10003 format (9es14.6)
10004 format (3x,3i6,27es14.6)

   end subroutine prn_chircoup

   !> Print directions and strengths of fourx couplings
   subroutine prn_fourxcoup(Natom, max_no_fourxneigh, fourx_listsize, fourx_list, fourx_val , simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_fourxneigh !< Calculated number of neighbours with fourx interactions
      integer,dimension(Natom), intent(in) :: fourx_listsize !< Size of neighbour list forfourx 
      integer,dimension(2,max_no_fourxneigh,Natom), intent(in) :: fourx_list   !< List of neighbours forfourx 
      real(dblprec),dimension(max_no_fourxneigh,Natom), intent(in) :: fourx_val  !< fourx force constant matrix
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i,j
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''fourxdata.'',a,''.out'')') trim(simid)
      open(21, file=filn)

      ! print neighbor list - after sort
      write (21,*) "Sorted fourx couplings "
      do i=1,Natom
         write (21,*) "----------------------------------"
         write (21,10001) i,fourx_listsize(i)
         do j=1,fourx_listsize(i)
            write (21,10004) i,fourx_list(1:2,j,i),fourx_val(j,i)
         end do
      end do
      close(21)

10001 format ("Atom=",i8,4x,"No neigh=",i7)
10002 format ("            ",1X,2I6)
!10002 format ("            ",1X,5I6)
10003 format (9es14.6)
10004 format (3x,3i6,27es14.6)

   end subroutine prn_fourxcoup

   !----------------------------------------------------------------------------
   !> Print directions and strengths of BIQDM couplings
   !----------------------------------------------------------------------------
   subroutine prn_biqdmcoup(Natom,nn_biqdm_tot,biqdmlistsize,biqdmlist,biqdm_vect,  &
      simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      integer,dimension(Natom), intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
      integer,dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
      real(dblprec),dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange vector
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''biqdmdata.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      ! print neighbor list - after sort
      write (ofileno,*) "Sorted data from heisge "
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         write (ofileno,10001) i,biqdmlistsize(i)
         write (ofileno,10002) biqdmlist(1:biqdmlistsize(i),i)
         write (ofileno,10003) biqdm_vect(1,1:biqdmlistsize(i),i)*mub/mry
      end do
      close(ofileno)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (6es14.6)

   end subroutine prn_biqdmcoup


   !> Print strength of exchange couplings
   subroutine prn_bqcoup(Natom, nn_bq_tot, bqlistsize, bqlist, j_bq, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
      integer,dimension(Natom), intent(in) :: bqlistsize !< Size of neighbour list for BQ
      integer,dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
      real(dblprec),dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
      character(len=8),intent(in) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=20) :: filn

      !.. Executable statements
      write (filn,'(''bqdata.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      ! print neighbor list - after sort
      write (ofileno,*) "Sorted data from heisge0"
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         write (ofileno,10001) i,bqlistsize(i)
         write (ofileno,10002) bqlist(1:bqlistsize(i),i)
         write (ofileno,10003) j_bq(1:bqlistsize(i),i)*mub/mry
      end do
      close(ofileno)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (5es16.8)

   end subroutine prn_bqcoup

  !> Print strength of four-spin ring couplings
  subroutine prn_ringcoup(Natom, nn_ring_tot, ringlist, ringlistsize, j_ring, simid)
    !
    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom        !< Number of atoms in system
    integer, intent(in) :: nn_ring_tot  !< Calculated number of neighbours with ring interactions
    integer, dimension(Natom,nn_ring_tot,3), intent(in) :: ringlist     !< List of neighbours for ring exchange
    integer, dimension(Natom), intent(in) :: ringlistsize               !< Size of neighbour list for 4SR coupling
    real(dblprec),dimension(Natom,nn_ring_tot), intent(in) :: j_ring    !< Ring exchange couplings
    character(len=8),intent(in) :: simid                                !< Name of simulation

    !.. Local variables
    integer           :: i, k, m, counter
    character(len=30) :: filn
    logical           :: selfmap
 
    !.. Executable statements
    write (filn,'(a,a,a)') 'ringdata.',trim(simid),'.out'
    open(ofileno, file=filn) 
  !print neighbor list
    do i=1,Natom
       counter=0
       write (ofileno,*) "----------------------------------"
       write (ofileno,10001) i, ringlistsize(i)   
         do k=1,ringlistsize(i)
            write (ofileno,*) "Plaquette", k
            selfmap=.false.
            if (i==ringlist(i,k,1).or.i==ringlist(i,k,2).or.i==ringlist(i,k,3)) then
                selfmap=.true.
            end if
            if (ringlist(i,k,1)==ringlist(i,k,2).or.ringlist(i,k,2)==ringlist(i,k,3).or. &
                ringlist(i,k,1)==ringlist(i,k,3)) then
                selfmap=.true.
            end if
            if (selfmap) then
            write (ofileno,*) "Some of the atoms in this plaquette are mapped onto itself!"
            write (ofileno,*) "The number of cells should be increased"
            end if
            write (ofileno,10002) ringlist(i,k,1), ringlist(i,k,2), ringlist(i,k,3)
         do m=1,counter
            if (ringlist(i,k,1)==ringlist(i,m,1).and.ringlist(i,k,2)==ringlist(i,m,2).and. &
               ringlist(i,k,3)==ringlist(i,m,3)) then
               write (ofileno,*) "Such a plaquette already exists for atom_i!"
               write (ofileno,*) "Check an input to avoid double counting"
            end if
         end do
         counter=counter+1
         write (ofileno,10003) j_ring(i,k)*mub/mry
         end do

    end do
    close(ofileno)


10001 format ("Atom_i:",i7,4x,"No plaquettes:",i7)
10002 format ("j: ",i7,4x,"k: ",i7,4x,"l: ",i7)
10003 format ("G(K)=",5es16.8)

  end subroutine prn_ringcoup


   !> Print single-ion anisotropies
   subroutine prn_anisotropy(Natom,NA,anisotropytype,taniso,eaniso,kaniso,simid,    &
      Nchmax)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, dimension(NA,Nchmax), intent(in) :: anisotropytype !< Type of anisotropies (0-2)
      integer, dimension(Natom), intent(in) :: taniso !< Type of anisotropy (0-2)
      real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
      real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i,j,is_ani
      character(len=20) :: filn

      !.. Executable statements

      ! check if no anisotropy is present, then return with out printing
      is_ani=1
      do i=1,NA
         do j=1,Nchmax
            is_ani=is_ani+anisotropytype(i,j)
         end do
      end do
      if(is_ani==0) return

      write (filn,'(''aniso1.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      ! print neighbor list - after sort
      write (ofileno,*) "List of all atoms and their anisotropy axes"
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         write (ofileno,10001) i,taniso(i)
         write (ofileno,*) eaniso(1:3,i)
         write (ofileno,10003) kaniso(1:2,i)*mub/mry
      end do
      close(ofileno)

      10001 format ("Atom=",i8,4x,"Anisotropy type:",i3)
      10003 format (5es16.8)

   end subroutine prn_anisotropy


   !----------------------------------------------------------------------------
   !> @brief Print structural information about the existing exchange couplings of the cluster
   !----------------------------------------------------------------------------
   subroutine clus_prnge(NA_clus,NT_clus,N1_clus,N2_clus,N3_clus,do_ralloy,         &
      Natom_clus,Nchmax_clus,Natom_full_clus,max_no_equiv_clus,max_no_shells_clus,  &
      NN_clus,atype_clus,achtype_clus,acellnumb_clus,acellnumbrev_clus,nmdim,nm,    &
      redcoord_clus,simid)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA_clus !< Number of atoms in the cluster
      integer, intent(in) :: NT_clus !< Number of types of atoms in the cluster
      integer, intent(in) :: N1_clus  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2_clus  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3_clus  !< Number of cell repetitions in z direction
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_clus !< Number of atoms in the cluster
      integer, intent(in) :: Nchmax_clus !< Max number of chemical components on each site in cell
      integer, intent(in) :: Natom_full_clus !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_equiv_clus !< Calculated maximum of neighbours in one shell for exchange
      integer, intent(in) :: max_no_shells_clus !< Calculated maximum of shells for exchange
      integer, dimension(NT_clus), intent(in) :: NN_clus !< Number of neighbour shells
      integer, dimension(Natom_clus), intent(in) :: atype_clus !< Type of atom
      integer, dimension(Natom_full_clus), intent(in) :: achtype_clus !< Chemical type of atoms (full list)
      integer, dimension(Natom_full_clus), intent(in) :: acellnumb_clus !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full_clus), intent(in) :: acellnumbrev_clus !< List for translating atom no. in actual cell to full cell
      integer, dimension(max_no_shells_clus,Natom_clus), intent(in) :: nmdim !< Dimension of neighbour map
      integer, dimension(Natom_clus,max_no_shells_clus,max_no_equiv_clus), intent(in) :: nm !< Neighbour map
      real(dblprec), dimension(NT_clus,max_no_shells_clus,3), intent(in) :: redcoord_clus !< Coordinates for Heisenberg exchange couplings
      character(len=8), intent(in) :: simid !< Name of simulation

      integer :: i, j, k, l, count, I0, I1, I2, I3
      integer :: iatom, jatom
      character(len=30) :: filn

      !.. Executable statements

      ! print neighbor map
      write (filn,'(''clus_struct.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      write (ofileno,*) "Clus data from heisgeinit"
      do I3=0, N3_clus-1
         do I2=0, N2_clus-1
            do I1=0, N1_clus-1
               do I0=1, NA_clus
                  i=I0+I1*NA_clus+I2*N1_clus*NA_clus+I3*N2_clus*N1_clus*NA_clus
                  if (do_ralloy==1.or.do_ralloy==2) i = acellnumb_clus(i)
                  if (i==0) cycle
                  write (ofileno,*) "------------------------------------------------------"
                  if (do_ralloy==0) then
                     write (ofileno,10001) i
                  else
                     iatom = acellnumbrev_clus(i)
                     write (ofileno,10011) i, atype_clus(i), achtype_clus(i)
                  end if
                  do k=1,NN_clus(atype_clus(i))
                     write (ofileno,10002) k, nmdim(k,i), redcoord_clus(atype_clus(I0),k,1),&
                        redcoord_clus(atype_clus(I0),k,2), redcoord_clus(atype_clus(I0),k,3),&
                        sqrt(sum(redcoord_clus(atype_clus(I0),k,:)**2))
                     write (ofileno,10003)   nm(i,k,1:nmdim(k,i))

                     ! new check for self-mapping
                     do l=1,nmdim(k,i)
                        if (nm(i,k,l)==i) then
                           write(*,'(1x,a,i6)') 'WARNING: Jii entry in neighbour map for atom',i
                           write (ofileno,'(1x,a)') 'WARNING: Jii entry in neighbour map'
                        endif
                     end do

                     if (do_ralloy==0) then
                        do j=1,NT_clus
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 if (atype_clus(nm(i,k,l))==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write (ofileno,10004) j, count
                        end do
                     else
                        do j=1,Nchmax_clus
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 jatom = acellnumbrev_clus(nm(i,k,l))
                                 if (achtype_clus(jatom)==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write (ofileno,10014) j, count
                        end do
                     end if
                  end do
               end do
            end do
         end do
      end do
      close(ofileno)

      10001 format ("Atom=",i8)
      10002 format ("Shell=",i4,2x,"Number of atoms=",i4,2x,"Shell coordinates:", 4f8.4)
      10003 format ("            ",1X,5I6)
      10004 format ("            Type=",i4,2x,"Number of atoms=",i4)

      10011 format ("Atom=",i8,2x,"Type=",i4,2x,"Chtype=",i4)
      10014 format ("            Chtype=",i4,2x,"Number of atoms=",i4)

   end subroutine clus_prnge

   !----------------------------------------------------------------------------
   !> Print structural information about the existing dm-exchange couplings of the cluster
   !----------------------------------------------------------------------------
   subroutine clus_dmprnge(NA_clus,N1_clus,N2_clus,N3_clus,NT_clus,do_ralloy,       &
      Natom_clus,Nchmax_clus,Natom_full_clus,max_no_equiv_clus,max_no_dmshells_clus,&
      dm_nn_clus,atype_clus,achtype_clus,acellnumb_clus,acellnumbrev_clus,nm,nmdim, &
      dm_redcoord_clus,simid)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA_clus !< Number of atoms in the cluster
      integer, intent(in) :: N1_clus  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2_clus  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3_clus  !< Number of cell repetitions in z direction
      integer, intent(in) :: NT_clus !< Number of types of atoms
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_clus !< Number of atoms in the cluster
      integer, intent(in) :: Nchmax_clus !< Max number of chemical components on each site in cell
      integer, intent(in) :: Natom_full_clus !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_equiv_clus !< Calculated maximum of neighbours in one shell for exchange
      integer, intent(in) :: max_no_dmshells_clus !< Calculated maximum of shells for exchange
      integer, dimension(NT_clus), intent(in) :: dm_nn_clus !< Number of neighbour shells
      integer, dimension(Natom_clus), intent(in) :: atype_clus !< Type of atom
      integer, dimension(Natom_full_clus), intent(in) :: achtype_clus !< Chemical type of atoms (full list)
      integer, dimension(Natom_full_clus), intent(in) :: acellnumb_clus !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full_clus), intent(in) :: acellnumbrev_clus !< List for translating atom no. in actual cell to full cell
      integer, dimension(Natom_clus,max_no_dmshells_clus,max_no_equiv_clus), intent(in) :: nm !< Neighbour map
      integer, dimension(max_no_dmshells_clus,Natom_clus), intent(in) :: nmdim !< Dimension of neighbour map
      real(dblprec), dimension(NT_clus,max_no_dmshells_clus,3), intent(in) :: dm_redcoord_clus !< Coordinates for Heisenberg exchange couplings
      character(len=8), intent(in) :: simid !< Name of simulation

      integer :: i, j, k, l, count, i0, i1, i2, i3
      integer :: iatom, jatom
      character(len=30) :: filn

      !.. Executable statements

      ! print neighbor map
      write (filn,'(''clus_dmstruct.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn)

      write (ofileno,*) "Data from dmheisgeinit"
      do I3=0, N3_clus-1
         do I2=0, N2_clus-1
            do I1=0, N1_clus-1
               do I0=1, NA_clus
                  i=I0+I1*NA_clus+I2*N1_clus*NA_clus+I3*N2_clus*N1_clus*NA_clus
                  if (do_ralloy==1.or.do_ralloy==2) i = acellnumb_clus(i)
                  if (i==0) cycle
                  write (ofileno,*) "------------------------------------------------------"
                  if (do_ralloy==0) then
                     write (ofileno,10001) i
                  else
                     iatom = acellnumbrev_clus(i)
                     write (ofileno,10011) i, atype_clus(i), achtype_clus(i)
                  end if
                  do k=1,dm_nn_clus(atype_clus(i))
                     write (ofileno,10002) k, nmdim(k,i), dm_redcoord_clus(atype_clus(I0),k,1),&
                        dm_redcoord_clus(atype_clus(I0),k,2), dm_redcoord_clus(atype_clus(I0),k,3),&
                        sqrt(sum(dm_redcoord_clus(atype_clus(I0),k,:)**2))
                     write (ofileno,10003)   nm(i,k,1:nmdim(k,i))

                     ! new check for self-mapping
                     do l=1,nmdim(k,i)
                        if (nm(i,k,l)==i) then
                           write(*,'(1x,a,i6)') 'WARNING: Dii entry in neighbour map for atom',i
                           write (ofileno,'(1x,a)') 'WARNING: Dii entry in neighbour map'
                        endif
                     end do

                     if (do_ralloy==0) then
                        do j=1,NT_clus
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 if (atype_clus(nm(i,k,l))==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write (ofileno,10004) j, count
                        end do
                     else
                        do j=1,Nchmax_clus
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 jatom = acellnumbrev_clus(nm(i,k,l))
                                 if (achtype_clus(jatom)==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write (ofileno,10014) j, count
                        end do
                     end if
                  end do
               end do
            end do
         end do
      end do
      close(ofileno)

      10001 format ("Atom=",i8)
      10002 format ("Shell=",i4,2x,"Number of atoms=",i4,2x,"Shell coordinates:", 4f8.4)
      10003 format ("            ",1X,5I6)
      10004 format ("            Type=",i4,2x,"Number of atoms=",i4)

      10011 format ("Atom=",i8,2x,"Type=",i4,2x,"Chtype=",i4)
      10014 format ("            Chtype=",i4,2x,"Number of atoms=",i4)

   end subroutine clus_dmprnge

end module PrintHamiltonian
