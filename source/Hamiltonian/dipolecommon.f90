!-------------------------------------------------------------------------------
! MODULE: DipoleCommon
!> @brief Container module for auxiliary functions and routines common between
!> different methods to calculate the dipole-dipole interaction.
!> @author Jonathan Chico
!-------------------------------------------------------------------------------
module DipoleCommon

   use Profiling
   use Parameters

   implicit none

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: read_qdip
   !> @brief Routine to read in the tensor dipole-dipole interaction
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine read_qdip(N_size,Dip_tens,qdip_files)

      implicit none

      integer, intent(in) :: N_size ! Size of the number of particles to where the dipole dipole is written
      real(dblprec), dimension(3,3,N_size,N_size), intent(out) :: Dip_tens ! Tensorial dipole-dipole
      character(len=30), intent(in) :: qdip_files !< Input file that contains the dipole-dipole tensor
      ! .. Local variables
      integer :: ii, jj, mu, nu,k,l,m,n

      open(ifileno,file=trim(qdip_files))

      do ii=1,N_size
         do jj=1,N_size
            do mu=1,3
               do nu=1,3
                  read(ifileno,*)k,l,m,n,Dip_tens(nu,mu,jj,ii)
               enddo
            enddo
         enddo
      enddo

   end subroutine read_qdip

   !----------------------------------------------------------------------------
   ! FUNCTION: dipoleMatrix
   !> @brief Calculation of the entries for the distance dependence of the dipole
   !> dipole interaction
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   function dipoleMatrix(alat,Rij) result(DipMat)

      use Constants

      implicit none

      real(dblprec), intent(in) :: alat
      real(dblprec), dimension(3), intent(in) :: Rij
      real(dblprec), dimension(3,3) :: DipMat

      integer :: mu, nu
      real(dblprec) :: R2,Rm5
      real(dblprec) :: fac,tol

      tol=1e-9
      ! fac determines the strength of the dipole interaction
      fac= (1.0d-7/alat**3)*mub

      R2=Rij(1)**2+Rij(2)**2+Rij(3)**2
      Rm5=R2**(-2.5)*fac
      ! If the modulus is larger than a certain tolerance calculate the matrix
      if (R2>tol) then
         ! Calculate the entry for the dipole matrix
         do mu=1,3
            do nu=1,3
               DipMat(nu,mu)=3.0_dblprec*Rm5*Rij(nu)*Rij(mu)
            end do
            DipMat(mu,mu)=DipMat(mu,mu)-Rm5*R2
         end do
      ! If not set the matrix to zero
      else
         DipMat=0.0_dblprec
      endif
   end function dipoleMatrix

   !----------------------------------------------------------------------------
   ! FUNCTION: paddingShift
   !> @brief Function to get the correct order of the padded arrays
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   function paddingShift(N_max,I_curr) result(I_indx)

      implicit none

      integer, intent(in) :: N_max
      integer, intent(in) :: I_curr

      integer :: I_indx

      if (I_curr<N_max) then
         I_indx=-I_curr
      else
         I_indx=2*N_max-1-I_curr
      endif

   end function paddingShift

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_qdip
   !> Set up the dipole-dipole interaction matrix
   !! @todo Change to correct scale, needs length information
   !! @todo Consider Ewald summation
   !----------------------------------------------------------------------------
   subroutine setup_qdip(Natom,coord,alat,Qdip,simid,print_dip_tensor)
      !
      use Constants
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), intent(in) :: alat !< Lattice parameter
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: print_dip_tensor !< Flag to print the dipole tensor
      real(dblprec), dimension(3,3,Natom,Natom), intent(out):: Qdip !< Matrix for dipole-dipole interaction

      ! ...Local Variables...
      integer :: ii,jj,mu,nu
      real(dblprec),dimension(3) :: Rij
      character(len=23) :: output_file
      ! Creation of a file to writeout the geometry of the macro cell
      output_file = 'dip_tensor.'//simid//'.out'
      !
      !$omp parallel do private(ii,jj,Rij)
      do ii=1,Natom
         do jj=1,Natom
            Rij=(coord(:,ii)-coord(:,jj))
            Qdip(:,:,jj,ii)=dipoleMatrix(alat,Rij)
         end do
         Qdip(:,:,ii,ii)=0.0_dblprec
      end do
      !$omp end parallel do
      !
      if (print_dip_tensor=='Y') then
         open(ofileno,file=output_file)
         do ii=1,Natom
            do jj=1,Natom
               do mu=1,3
                  do nu=1,3
                     write(ofileno,'(i6,i6,i6,i6,E16.8)') ii, jj, mu, nu, Qdip(nu,mu,jj,ii)
                  enddo
               enddo
            enddo
         enddo
         close(ofileno)
      endif
      return
      !
   end subroutine setup_qdip

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: setup_macro_qdip
   !> @brief Set up the dipole-dipole interaction matrix as outlined in
   !> J. Phys.: Condens. Matter 28 (2016) 066001
   !> @details This field is in Teslas, this is determined by the factor, fac,
   !> the unit analsis shows
   !>
   !> \f$ fac=\frac{\mu_0}{4\pi a_{lat}^3}\mu_B\f$
   !>
   !> \f$ fac=\frac{N}{A^2}\frac{1}{m^3}\frac{J}{T}\f$
   !>
   !> \f$ fac=\frac{kg^2m^3As^2}{s^4A^2kgm^3}=\frac{Kg}{s^2A}=T\f$
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine setup_macro_qdip(Natom,Num_macro,max_num_atom_macro_cell,             &
      macro_nlistsize,macro_atom_nlist,coord,alat,Qdip_macro,simid,print_dip_tensor)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: max_num_atom_macro_cell !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist !< List containing the information of which atoms are in a given macrocell
      real(dblprec), intent(in) :: alat !< Lattice parameter
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: print_dip_tensor !< Flag to print the dipole tensor
      real(dblprec), dimension(3,3,Num_macro,Num_macro), intent(out):: Qdip_macro !< Matrix for the macro spin dipole-dipole interaction

      ! ...Local Variables...
      integer :: ii,jj,isite,jsite,iatom,jatom,mu,nu
      real(dblprec),dimension(3) :: Rij
      character(len=23) :: output_file
      !
      output_file = 'dip_tensor.'//simid//'.out'

      !$omp parallel do default(shared) private(ii,jj,iatom,jatom,isite,jsite,Rij)
      do ii=1,Num_macro
         do jj=1,Num_macro
            do isite=1, macro_nlistsize(ii)
               iatom=macro_atom_nlist(ii,isite)
               do jsite=1, macro_nlistsize(jj)
                  jatom=macro_atom_nlist(jj,jsite)
                  if (iatom.ne.jatom) then
                     Rij=(coord(:,iatom)-coord(:,jatom))
                     Qdip_macro(:,:,ii,jj)=Qdip_macro(:,:,ii,jj)+&
                        dipoleMatrix(alat,Rij)/(macro_nlistsize(ii)*macro_nlistsize(jj))
                  endif
               enddo
            enddo
         end do
      end do
      !$omp end parallel do
      !
      if (print_dip_tensor=='Y') then
         open(ofileno,file=output_file)
         do ii=1,Num_macro
            do jj=1,Num_macro
               do mu=1,3
                  do nu=1,3
                     write(ofileno,'(i6,i6,i6,i6,E16.8)') ii, jj, mu, nu, Qdip_macro(nu,mu,jj,ii)
                  enddo
               enddo
            enddo
         enddo
         close(ofileno)
      endif
      return
      !
   end subroutine setup_macro_qdip

end module DipoleCommon
