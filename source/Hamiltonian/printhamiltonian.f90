!> Routines for printing information about the Hamiltonian
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
module PrintHamiltonian
   use Parameters
   use Profiling
   use Constants

   implicit none
   public


contains


   !> Print structural information about the existing exchange couplings
   subroutine prnge(simid, Natom, NT, NA, N1, N2, N3,  atype, &
         max_no_shells, max_no_equiv, redcoord, nn, nm, nmdim, &
         do_ralloy, Natom_full, Nchmax, acellnumb, acellnumbrev, achtype)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, intent(in) :: NT !< Number of types of atoms
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      real(dblprec), dimension(NT,max_no_shells,3), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings
      integer, dimension(NT), intent(in) :: nn !< Number of neighbour shells
      integer, dimension(Natom,max_no_shells,max_no_equiv), intent(in) :: nm !< Neighbour map
      integer, dimension(max_no_shells,Natom), intent(in) :: nmdim !< Dimension of neighbour map
      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, dimension(Natom_full), intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full), intent(in) :: acellnumbrev !< List for translating atom no. in actual cell to full cell
      integer, dimension(Natom_full), intent(in) :: achtype !< Chemical type of atoms (full list)

      integer :: i, j, k, l, count, i0, i1, i2, i3
      integer :: iatom, jatom
      character(len=30) :: filn

      !.. Executable statements

      ! print neighbor map
      write (filn,'(''struct.'',a8,''.out'')') simid
      open(ofileno, file=filn)

      write (ofileno,*) "Data from heisgeinit"
      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA
                  i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                  if (do_ralloy==1) i = acellnumb(i)
                  if (i==0) cycle
                  write (ofileno,*) "------------------------------------------------------"
                  if (do_ralloy==0) then
                     write (ofileno,10001) i
                  else
                     iatom = acellnumbrev(i)
                     write (ofileno,10011) i, atype(i), achtype(i)
                  end if
                  do k=1,nn(atype(i))
                     write (ofileno,10002) k, nmdim(k,i), redcoord(atype(i0),k,1),&
                        redcoord(atype(i0),k,2), redcoord(atype(i0),k,3),&
                        sqrt(sum(redcoord(atype(i0),k,:)**2))
                     write (ofileno,10003)   nm(i,k,1:nmdim(k,i))

                     ! new check for self-mapping
                     do l=1,nmdim(k,i)
                        if (nm(i,k,l)==i) then
                           write(*,'(1x,a,i6)') 'WARNING: Jii entry in neighbour map for atom',i
                           write (ofileno,'(1x,a)') 'WARNING: Jii entry in neighbour map'
                        endif
                     end do

                     if (do_ralloy==0) then
                        do j=1,NT
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 if (atype(nm(i,k,l))==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write (ofileno,10004) j, count
                        end do
                     else
                        do j=1,Nchmax
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 jatom = acellnumbrev(nm(i,k,l))
                                 if (achtype(jatom)==j) then
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

   end subroutine prnge

   !> Print strength of exchange couplings
   subroutine prn_exchange(Natom, max_no_neigh, nlistsize, nlist, ncoup, simid, mdim)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh, Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, intent(in) :: mdim !< dimension of the exchange coupling matrix (1=scalar or 9=3x3)
      real(dblprec), dimension(mdim,max_no_neigh, Natom), intent(in) :: ncoup !< Heisenberg exchange couplings
      character(len=8),intent(in) :: simid !< Name of simulation

      !.. Local variables
      integer :: i,j
      character(len=20) :: filn

      !.. Executable statements
      write (filn,'(''struct1.'',a8,''.out'')') simid
      open(ofileno, file=filn)

      ! print neighbor list - after sort
      write (ofileno,*) "Sorted data from heisge0"
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         write (ofileno,10001) i,nlistsize(i)
         write (ofileno,10002) nlist(1:nlistsize(i),i)
         if(mdim==1) then
            write (ofileno,10003) ncoup(1:mdim,1:nlistsize(i),i)*mub/mry
         else
            do j=1,nlistsize(i)
               write (ofileno,10004) ncoup(1:mdim,j,i)*mub/mry
            end do
         end if
      end do
      close(ofileno)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (5es16.8)
      10004 format (9es16.8)

   end subroutine prn_exchange

   !> Print strength of exchange couplings in sparse matrix format
   !> The output is in Tesla
   subroutine prn_exchange_sparse(Natom, max_no_neigh, nlistsize, nlist, ncoup, simid, mdim)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh, Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, intent(in) :: mdim !< dimension of the exchange coupling matrix (1=scalar or 9=3x3)
      real(dblprec), dimension(mdim,max_no_neigh, Natom), intent(in) :: ncoup !< Heisenberg exchange couplings
      character(len=8),intent(in) :: simid !< Name of simulation

      !.. Local variables
      integer :: i,j
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''structsparse.'',a8,''.out'')') simid
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

   !>  Prints the neighbour list for the induced moments
   subroutine prn_ind_exchange(Natom,max_no_neigh,ind_nlistsize,ind_mom_nlist,simid)

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: max_no_neigh
      integer, dimension(Natom), intent(in) :: ind_nlistsize
      integer, dimension(max_no_neigh,Natom) :: ind_mom_nlist
      character(len=8), intent(in) :: simid

      ! Local variables
      integer :: i
      character(len=30) :: filn

      !.. Printing the induced moments list
      write (filn,'(''indstruct.'',a8,''.out'')') simid
      open(ofileno, file=filn)

      ! Printing the new induced moments neighbour lists
      write (ofileno,*) "Sorted data from heisge0 for induced moments"
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         write (ofileno,10001) i,ind_nlistsize(i)
         write (ofileno,10002) ind_mom_nlist(1:ind_nlistsize(i),i)
      end do
      close(ofileno)


      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)

   end subroutine prn_ind_exchange

   !> Print directions and strengths of DM couplings
   subroutine prn_dmcoup(Natom, max_no_dmneigh, dmlistsize, dmlist, dm_vect, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      integer,dimension(Natom), intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer,dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
      real(dblprec),dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=20) :: filn

      !.. Executable statements
      write (filn,'(''dmdata.'',a8,''.out'')') simid
      open(ofileno, file=filn)

      ! print neighbor list - after sort
      write (ofileno,*) "Sorted data from heisge "
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         write (ofileno,10001) i,dmlistsize(i)
         write (ofileno,10002) dmlist(1:dmlistsize(i),i)
         write (ofileno,10003) dm_vect(1:3,1:dmlistsize(i),i)*mub/mry
      end do
      close(ofileno)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (6es14.6)

   end subroutine prn_dmcoup


   !> Print directions and strengths of PD couplings
   subroutine prn_pdcoup(Natom, nn_pd_tot, pdlistsize, pdlist, pd_vect, simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      integer,dimension(Natom), intent(in) :: pdlistsize !< Size of neighbour list for PD
      integer,dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
      real(dblprec),dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
      character(len=8) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=20) :: filn

      !.. Executable statements
      write (filn,'(''pddata.'',a8,''.out'')') simid
      open(ofileno, file=filn)

      ! print neighbor list - after sort
      write (ofileno,*) "Sorted data from heisge "
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         write (ofileno,10001) i,pdlistsize(i)
         write (ofileno,10002) pdlist(1:pdlistsize(i),i)
         write (ofileno,10003) pd_vect(1:6,1:pdlistsize(i),i)*mub/mry
      end do
      close(ofileno)

      10001 format ("Atom=",i8,4x,"No neigh=",i7)
      10002 format ("            ",1X,5I6)
      10003 format (6es14.6)

   end subroutine prn_pdcoup


   !> Print directions and strengths of BIQDM couplings
   subroutine prn_biqdmcoup(Natom, nn_biqdm_tot, biqdmlistsize, biqdmlist, biqdm_vect, simid)
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
      write (filn,'(''biqdmdata.'',a8,''.out'')') simid
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
      write (filn,'(''bqdata.'',a8,''.out'')') simid
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


   !> Print single-ion anisotropies
   subroutine prn_anisotropy(Natom, NA, anisotropytype, taniso, eaniso, kaniso, simid,Nchmax)
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

      write (filn,'(''aniso1.'',a8,''.out'')') simid
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

   !> Print structural information about the existing dmexchange couplings
   subroutine dmprnge(simid, Natom, NT, NA, N1, N2, N3,  atype, &
         max_no_shells, max_no_equiv, redcoord, nn, nm, nmdim, &
         do_ralloy, Natom_full, Nchmax, acellnumb, acellnumbrev, achtype)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, intent(in) :: NT !< Number of types of atoms
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      real(dblprec), dimension(NT,max_no_shells,3), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings
      integer, dimension(NT), intent(in) :: nn !< Number of neighbour shells
      integer, dimension(Natom,max_no_shells,max_no_equiv), intent(in) :: nm !< Neighbour map
      integer, dimension(max_no_shells,Natom), intent(in) :: nmdim !< Dimension of neighbour map
      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, dimension(Natom_full), intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full), intent(in) :: acellnumbrev !< List for translating atom no. in actual cell to full cell
      integer, dimension(Natom_full), intent(in) :: achtype !< Chemical type of atoms (full list)

      integer :: i, j, k, l, count, i0, i1, i2, i3
      integer :: iatom, jatom
      character(len=30) :: filn

      !.. Executable statements

      ! print neighbor map
      write (filn,'(''dmstruct.'',a8,''.out'')') simid
      open(ofileno, file=filn)

      write (ofileno,*) "Data from dmheisgeinit"
      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA
                  i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                  if (do_ralloy==1) i = acellnumb(i)
                  if (i==0) cycle
                  write (ofileno,*) "------------------------------------------------------"
                  if (do_ralloy==0) then
                     write (ofileno,10001) i
                  else
                     iatom = acellnumbrev(i)
                     write (ofileno,10011) i, atype(i), achtype(i)
                  end if
                  do k=1,nn(atype(i))
                     write (ofileno,10002) k, nmdim(k,i), redcoord(atype(i0),k,1),&
                        redcoord(atype(i0),k,2), redcoord(atype(i0),k,3),&
                        sqrt(sum(redcoord(atype(i0),k,:)**2))
                     write (ofileno,10003)   nm(i,k,1:nmdim(k,i))

                     ! new check for self-mapping
                     do l=1,nmdim(k,i)
                        if (nm(i,k,l)==i) then
                           write(*,'(1x,a,i6)') 'WARNING: Dii entry in neighbour map for atom',i
                           write (ofileno,'(1x,a)') 'WARNING: Dii entry in neighbour map'
                        endif
                     end do

                     if (do_ralloy==0) then
                        do j=1,NT
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 if (atype(nm(i,k,l))==j) then
                                    count=count+1
                                 endif
                              end if
                           end do
                           write (ofileno,10004) j, count
                        end do
                     else
                        do j=1,Nchmax
                           count=0
                           do l=1,nmdim(k,i)
                              if (nm(i,k,l)/=0) then
                                 jatom = acellnumbrev(nm(i,k,l))
                                 if (achtype(jatom)==j) then
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

   end subroutine dmprnge

end module PrintHamiltonian
