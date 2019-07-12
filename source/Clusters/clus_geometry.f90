!-------------------------------------------------------------------------------
! MODULE: clus_geometry
!> @brief Module taking care of the routines dealing with the geometrical information
!> of the generalized cluster method
!> @author Jonathan Chico
!-------------------------------------------------------------------------------
module clus_geometry

   use Parameters
   use Profiling
   use Geometry, only : setup_type_and_numb,setup_globcoord,setup_chemicaldata

   implicit none

   integer :: clus_expand  !< Number of atoms in the cluster that do not belong to the host
   real(dblprec) :: chconceff_clus !< Effective chemical concentration
   integer, dimension(:), allocatable :: index_clus   !< Mapping of cluster indices to host indices
   integer, dimension(:), allocatable :: anumb_clus   !< Atom number in cell
   integer, dimension(:), allocatable :: atype_clus   !< Type of atom
   integer, dimension(:), allocatable :: achtype_clus !< Chemical type of atoms (full list)
   integer, dimension(:), allocatable :: acellnumb_clus     !< List for translating atom no. in full cell to actual cell
   integer, dimension(:), allocatable :: acellnumbrev_clus  !< List for translating atom no. in actual cell to full cell
   integer, dimension(:), allocatable :: achem_ch_clus !< Chemical type of atoms (reduced list) (achem_ch(i)=achtype(acellnumbrev(i)))
   integer, dimension(:), allocatable :: atype_ch_clus !< Actual site of atom for dilute system
   integer, dimension(:), allocatable :: asite_ch_clus !< Actual site of atom for dilute system
   real(dblprec), dimension(:,:), allocatable :: coord_clus !< Coordinates of all atoms belonging to the cluster

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_clus_geometry
   !> @brief Setup the data structures and information of the cluster
   !> @details This routine allows to have generalized clusters embedded on the sample
   !> it allows for the repetition of clusters along lattice vectors.
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine setup_clus_geometry(NA,N1,N2,N3,N1_clus,N2_clus,N3_clus,NA_clus,      &
      do_ralloy,block_size,Natom_full,Nchmax_clus,Natom_full_clus,Nch_clus,tseed,   &
      C1,C2,C3,C1_clus,C2_clus,C3_clus,Bas,chconc_clus,clus_expand,index_clus,         &
      Natom_clus,atype_inp_clus,anumb_inp_clus,Bas_clus,coord_clus)

      implicit none

      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: tseed        !< Temperature seed
      integer, intent(in) :: N1_clus      !< Number of cell repetitions in x direction for the cluster
      integer, intent(in) :: N2_clus      !< Number of cell repetitions in y direction for the cluster
      integer, intent(in) :: N3_clus      !< Number of cell repetitions in z direction for the cluster
      integer, intent(in) :: NA_clus      !< Number of atoms in the cluster unit cell
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      integer, intent(in) :: Nchmax_clus     !< Max number of chemical components on each site in cell for the cluster
      integer, intent(in) :: Natom_full_clus !< Number of atoms for the cluster (=Natom if not dilute)
      integer, dimension(NA_clus), intent(in) :: Nch_clus !< Number of chemical components on each site in cell
      real(dblprec), dimension(3), intent(in) :: C1      !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2      !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3      !< Third lattice vector
      real(dblprec), dimension(3), intent(in) :: C1_clus !< First lattice vector for the cluster
      real(dblprec), dimension(3), intent(in) :: C2_clus !< Second lattice vector for the cluster
      real(dblprec), dimension(3), intent(in) :: C3_clus !< Third lattice vector for the cluster
      real(dblprec), dimension(NA_clus,Nchmax_clus), intent(in) :: chconc_clus !< Chemical concentration on sites
      !.. In/out variables
      integer, intent(inout) :: Natom_clus   !< Number of atoms in the cluster
      integer, intent(inout) :: Natom_full   !< Number of atoms for the full system (=Natom if not dilute)
      integer, intent(inout) :: clus_expand  !< Number of atoms in the cluster that do not belong to the host
      integer, dimension(NA_clus), intent(inout) :: atype_inp_clus  !< Type of atom from input
      integer, dimension(NA_clus), intent(inout) :: anumb_inp_clus !< Atom number in cell from input
      integer, dimension(:), allocatable, intent(inout) :: index_clus   !< Mapping of cluster indices to host indices
      real(dblprec), dimension(3,NA), intent(inout) :: Bas           !< Coordinates for basis atoms
      real(dblprec), dimension(3,NA_clus), intent(inout) :: Bas_clus !< Coordinates for basis atoms for the cluster
      real(dblprec), dimension(:,:), allocatable, intent(inout) :: coord_clus !< Coordinates of all atoms belonging to the cluster

      ! Define the number of atoms in the cluster sample
      Natom_clus=NA_clus*N1_clus*N2_clus*N3_clus
      ! Allocate the type and site data needed for the cluster
      call allocate_clusdata(Natom_clus,1)
      ! Setup the type and site data of the generalized cluster
      call setup_type_and_numb_clus(Natom_full_clus,NA_clus,N1_clus,N2_clus,N3_clus,&
         atype_clus,anumb_clus,atype_inp_clus,anumb_inp_clus,block_size)
      ! Setup the chemical information of the generalized cluster
      if (do_ralloy==1) then
         call allocate_chemicaldata_clus(Natom_clus,1)
         call setup_chemicaldata_clus(tseed,NA_clus,N1_clus,N2_clus,N3_clus,        &
            do_ralloy,Nchmax_clus,Natom_full_clus,Nch_clus,chconc_clus,achtype_clus,&
            atype_ch_clus,asite_ch_clus,achem_ch_clus,acellnumb_clus,               &
            acellnumbrev_clus,Natom_clus,atype_clus)
      else
         Natom_clus=Natom_full_clus
      endif
      ! Define the full coordinates of the generalized cluster
      call setup_glob_clus(NA,N1,N2,N3,N1_clus,N2_clus,N3_clus,NA_clus,do_ralloy,   &
         block_size,Natom_full,C1,C2,C3,C1_clus,C2_clus,C3_clus,Bas,clus_expand,    &
         index_clus,Natom_clus,Bas_clus,coord_clus)
      ! Set that the largest number of atoms is given by the Natom_full for the host
      ! plus the number of atoms that need to be expanded due to the cluster
      Natom_full=Natom_full+clus_expand

   end subroutine setup_clus_geometry

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_glob_clus
   !> @brief Calculation of the global coordinates for the cluster
   !> @details Routine based on the setup_globcoord subroutine, it also calculates if
   !> there are atoms outside the magnetic layers and sets up the pointer arrays
   !> mapping the cluster indices to the atom indices.
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine setup_glob_clus(NA,N1,N2,N3,N1_clus,N2_clus,N3_clus,NA_clus,do_ralloy,&
      block_size,Natom_full,C1,C2,C3,C1_clus,C2_clus,C3_clus,Bas,clus_expand,       &
      index_clus,Natom_clus,Bas_clus,coord_clus)

      implicit none

      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: N1_clus   !< Number of cell repetitions in x direction for the cluster
      integer, intent(in) :: N2_clus   !< Number of cell repetitions in y direction for the cluster
      integer, intent(in) :: N3_clus   !< Number of cell repetitions in z direction for the cluster
      integer, intent(in) :: NA_clus      !< Number of atoms in the cluster unit cell
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      integer, intent(in) :: Natom_full   !< Number of atoms for the full system (=Natom if not dilute)
      real(dblprec), dimension(3), intent(in) :: C1      !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2      !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3      !< Third lattice vector
      real(dblprec), dimension(3), intent(in) :: C1_clus !< First lattice vector for the cluster
      real(dblprec), dimension(3), intent(in) :: C2_clus !< Second lattice vector for the cluster
      real(dblprec), dimension(3), intent(in) :: C3_clus !< Third lattice vector for the cluster
      !.. In/out variables
      integer, intent(inout) :: Natom_clus   !< Number of atoms in the cluster
      integer, intent(inout) :: clus_expand  !< Number of atoms in the cluster that do not belong to the host
      integer, dimension(:), allocatable, intent(inout) :: index_clus   !< Mapping of cluster indices to host indices
      real(dblprec), dimension(3,NA), intent(inout) :: Bas  !< Coordinates for basis atoms
      real(dblprec), dimension(3,NA_clus), intent(inout) :: Bas_clus !< Coordinates for basis atoms for the cluster
      real(dblprec), dimension(:,:), allocatable, intent(inout) :: coord_clus !< Coordinates of all atoms belonging to the cluster
      ! .. Local variables
      integer :: I1,I2,I3,I0
      integer :: II3,II2,II1
      integer :: iatom,iclus,i
      integer :: i_stat,i_all
      real(dblprec) :: mod,tol
      real(dblprec) :: detmatrix
      real(dblprec), dimension(3) :: icvec,bsf
      real(dblprec), dimension(3,3) :: invmatrix
      real(dblprec), dimension(:,:), allocatable :: temp_coord
      logical :: not_found

      tol=0.005_dblprec

      not_found=.True.

      allocate(index_clus(Natom_clus),stat=i_stat)
      call memocc(i_stat,product(shape(index_clus))*kind(index_clus),'index_clus','setup_glob_clus')
      index_clus=0
      !
      ! Fold all atoms into first unit cell:
      ! Calculate inverse of basis matrix
      detmatrix=  C1_clus(1)*C2_clus(2)*C3_clus(3)-C1_clus(1)*C2_clus(3)*C3_clus(2)+&
                  C1_clus(2)*C2_clus(3)*C3_clus(1)-C1_clus(2)*C2_clus(1)*C3_clus(3)+&
                  C1_clus(3)*C2_clus(1)*C3_clus(2)-C1_clus(3)*C2_clus(2)*C3_clus(1)
      invmatrix=0.0_dblprec
      if (abs(detmatrix)>dbl_tolerance) then
         invmatrix(1,1)=(C2_clus(2)*C3_clus(3)-C3_clus(2)*C2_clus(3))/detmatrix
         invmatrix(1,2)=(C1_clus(3)*C3_clus(2)-C3_clus(3)*C1_clus(2))/detmatrix
         invmatrix(1,3)=(C1_clus(2)*C2_clus(3)-C2_clus(2)*C1_clus(3))/detmatrix
         invmatrix(2,1)=(C2_clus(3)*C3_clus(1)-C3_clus(3)*C2_clus(1))/detmatrix
         invmatrix(2,2)=(C1_clus(1)*C3_clus(3)-C3_clus(1)*C1_clus(3))/detmatrix
         invmatrix(2,3)=(C1_clus(3)*C2_clus(1)-C2_clus(3)*C1_clus(1))/detmatrix
         invmatrix(3,1)=(C2_clus(1)*C3_clus(2)-C3_clus(1)*C2_clus(2))/detmatrix
         invmatrix(3,2)=(C1_clus(2)*C3_clus(1)-C3_clus(2)*C1_clus(1))/detmatrix
         invmatrix(3,3)=(C1_clus(1)*C2_clus(2)-C2_clus(1)*C1_clus(2))/detmatrix
      end if

      do I0=1,NA_clus
         !find coordinate vector in basis coordinates
         icvec(1)=Bas_clus(1,I0)*invmatrix(1,1)+Bas_clus(2,I0)*invmatrix(2,1)+Bas_clus(3,I0)*invmatrix(3,1)
         icvec(2)=Bas_clus(1,I0)*invmatrix(1,2)+Bas_clus(2,I0)*invmatrix(2,2)+Bas_clus(3,I0)*invmatrix(3,2)
         icvec(3)=Bas_clus(1,I0)*invmatrix(1,3)+Bas_clus(2,I0)*invmatrix(2,3)+Bas_clus(3,I0)*invmatrix(3,3)
         ! fold back to original cell
         bsf(1)=floor(icvec(1)+1d-7)
         bsf(2)=floor(icvec(2)+1d-7)
         bsf(3)=floor(icvec(3)+1d-7)
         !
         Bas_clus(1,I0)=Bas_clus(1,I0)-bsf(1)*C1_clus(1)-bsf(2)*C2_clus(1)-bsf(3)*C3_clus(1)
         Bas_clus(2,I0)=Bas_clus(2,I0)-bsf(1)*C1_clus(2)-bsf(2)*C2_clus(2)-bsf(3)*C3_clus(2)
         Bas_clus(3,I0)=Bas_clus(3,I0)-bsf(1)*C1_clus(3)-bsf(2)*C2_clus(3)-bsf(3)*C3_clus(3)
      end do

      ! Allocate coordinate array
      allocate(coord_clus(3,Natom_clus),stat=i_stat)
      call memocc(i_stat,product(shape(coord_clus))*kind(coord_clus),'coord_clus','setup_glob_clus')
      coord_clus=0.0_dblprec
      i=0
      do II3=0, N3_clus-1,block_size
         do II2=0, N2_clus-1,block_size
            do II1=0, N1_clus-1,block_size
               do I3=II3, min(II3+block_size-1,N3_clus-1)
                  do I2=II2,min(II2+block_size-1,N2_clus-1)
                     do I1=II1, min(II1+block_size-1,N1_clus-1)
                        do I0=1, NA_clus
                           i=i+1
                           if (do_ralloy==0) then
                              coord_clus(1,i)=I1*C1_clus(1)+I2*C2_clus(1)+I3*C3_clus(1)+Bas_clus(1,I0)
                              coord_clus(2,i)=I1*C1_clus(2)+I2*C2_clus(2)+I3*C3_clus(2)+Bas_clus(2,I0)
                              coord_clus(3,i)=I1*C1_clus(3)+I2*C2_clus(3)+I3*C3_clus(3)+Bas_clus(3,I0)
                           else
                              if (achtype_clus(i) /= 0) then
                                 iatom=acellnumb_clus(i)
                                 coord_clus(1,iatom)=I1*C1_clus(1)+I2*C2_clus(1)+I3*C3_clus(1)+Bas_clus(1,I0)
                                 coord_clus(2,iatom)=I1*C1_clus(2)+I2*C2_clus(2)+I3*C3_clus(2)+Bas_clus(2,I0)
                                 coord_clus(3,iatom)=I1*C1_clus(3)+I2*C2_clus(3)+I3*C3_clus(3)+Bas_clus(3,I0)
                              end if
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! Fold all atoms into first unit cell:
      ! Calculate inverse of basis matrix
      detmatrix=  C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+&
                  C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+&
                  C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1)
      invmatrix=0.0_dblprec
      if (abs(detmatrix)>dbl_tolerance) then
         invmatrix(1,1)=(C2(2)*C3(3)-C3(2)*C2(3))/detmatrix
         invmatrix(1,2)=(C1(3)*C3(2)-C3(3)*C1(2))/detmatrix
         invmatrix(1,3)=(C1(2)*C2(3)-C2(2)*C1(3))/detmatrix
         invmatrix(2,1)=(C2(3)*C3(1)-C3(3)*C2(1))/detmatrix
         invmatrix(2,2)=(C1(1)*C3(3)-C3(1)*C1(3))/detmatrix
         invmatrix(2,3)=(C1(3)*C2(1)-C2(3)*C1(1))/detmatrix
         invmatrix(3,1)=(C2(1)*C3(2)-C3(1)*C2(2))/detmatrix
         invmatrix(3,2)=(C1(2)*C3(1)-C3(2)*C1(1))/detmatrix
         invmatrix(3,3)=(C1(1)*C2(2)-C2(1)*C1(2))/detmatrix
      end if

      do I0=1,NA
         !find coordinate vector in basis coordinates
         icvec(1)=Bas(1,I0)*invmatrix(1,1)+Bas(2,I0)*invmatrix(2,1)+Bas(3,I0)*invmatrix(3,1)
         icvec(2)=Bas(1,I0)*invmatrix(1,2)+Bas(2,I0)*invmatrix(2,2)+Bas(3,I0)*invmatrix(3,2)
         icvec(3)=Bas(1,I0)*invmatrix(1,3)+Bas(2,I0)*invmatrix(2,3)+Bas(3,I0)*invmatrix(3,3)
         ! fold back to original cell
         bsf(1)=floor(icvec(1)+1d-7)
         bsf(2)=floor(icvec(2)+1d-7)
         bsf(3)=floor(icvec(3)+1d-7)
         !
         Bas(1,I0)=Bas(1,I0)-bsf(1)*C1(1)-bsf(2)*C2(1)-bsf(3)*C3(1)
         Bas(2,I0)=Bas(2,I0)-bsf(1)*C1(2)-bsf(2)*C2(2)-bsf(3)*C3(2)
         Bas(3,I0)=Bas(3,I0)-bsf(1)*C1(3)-bsf(2)*C2(3)-bsf(3)*C3(3)
      end do

      ! Allocate coordinate array
      allocate(temp_coord(3,Natom_full),stat=i_stat)
      call memocc(i_stat,product(shape(temp_coord))*kind(temp_coord),'temp_coord','setup_glob_clus')
      i=0
      do II3=0, N3-1,block_size
         do II2=0, N2-1,block_size
            do II1=0, N1-1,block_size
               do I3=II3, min(II3+block_size-1,N3-1)
                  do I2=II2,min(II2+block_size-1,N2-1)
                     do I1=II1, min(II1+block_size-1,N1-1)
                        do I0=1, NA
                           i=i+1
                           temp_coord(1,i)=I1*C1(1)+I2*C2(1)+I3*C3(1)+Bas(1,I0)
                           temp_coord(2,i)=I1*C1(2)+I2*C2(2)+I3*C3(2)+Bas(2,I0)
                           temp_coord(3,i)=I1*C1(3)+I2*C2(3)+I3*C3(3)+Bas(3,I0)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Find if all the atoms belong to the cluster, if they do not find how many
      ! of them do not correspond
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      clus_expand=0
      do iclus=1,Natom_clus
         iatom=1
         not_found=.True.
         do while (iatom<=Natom_full .and. not_found)
            mod=sqrt((temp_coord(1,iatom)-coord_clus(1,iclus))**2+&
                     (temp_coord(2,iatom)-coord_clus(2,iclus))**2+&
                     (temp_coord(3,iatom)-coord_clus(3,iclus))**2)
            if (mod<tol) then
               index_clus(iclus)=iatom
               not_found=.False.
            endif
            iatom=iatom+1
         enddo
         if (not_found) then
            clus_expand=clus_expand+1
            index_clus(iclus)=Natom_full+clus_expand
         endif
      enddo
      i_all=-product(shape(temp_coord))*kind(temp_coord)
      deallocate(temp_coord,stat=i_stat)
      call memocc(i_stat,i_all,'temp_coord','setup_glob_clus')

   end subroutine setup_glob_clus

   !----------------------------------------------------------------------------
   ! SUBROUTINE: modify_global_coordinates
   !> @brief Modfify and setup the atomic types for the full system
   !> @details This routine setups the global arrays for the geometry of the sample
   !> and modifies the number of atoms to ensure that impurities outside of the
   !> magnetic layer are considered
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine modify_global_coordinates(Natom,NT,NA,N1,N2,N3,Bas,C1,C2,C3,coord,    &
      atype,anumb,atype_inp,anumb_inp,do_prnstruct,do_prn_poscar,tseed,simid,do_ralloy,Natom_full,&
      Natom_clus,Nchmax,Nch,acellnumb,acellnumbrev,achtype,chconc,atype_ch,asite_ch,&
      achem_ch,block_size)

      use RandomNumbers, only : rng_uniform
      use Sorting, only : mergesort
      implicit none

      integer, intent(in) :: NT  !< Number of atoms in one cell
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: tseed        !< Temperature seed
      integer, intent(in) :: Nchmax       !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      integer, intent(in) :: Natom_clus   !< Number of atoms in the cluster
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: do_prnstruct !< Print Hamiltonian information (0/1)
      integer, intent(in) :: do_prn_poscar    !< Print geometry on POSCAR format (0/1)
      character(len=8), intent(in) :: simid     !< Name of simulation
      integer, dimension(:), intent(in) :: Nch  !< Number of chemical components on each site in cell
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(:,:), intent(in) :: chconc !< Chemical concentration on sites
      ! .. In/Out variables
      integer, intent(inout) :: Natom !< Number of atoms in system
      integer, dimension(Natom_full), intent(inout)   :: atype       !< Type of atom
      integer, dimension(Natom_full), intent(inout)   :: anumb       !< Atom number in cell
      integer, dimension(NA), intent(inout)           :: atype_inp   !< Type of atom from input
      integer, dimension(NA), intent(inout)           :: anumb_inp   !< Atom number in cell from input
      real(dblprec), dimension(3,NA) , intent(inout) :: Bas !< Coordinates for basis atoms
      real(dblprec), dimension(:,:), allocatable, intent(inout) :: coord !< Coordinates of atoms
      ! .. Output variables
      integer, dimension(:), intent(out) :: achtype               !< Chemical type of atoms (full list)
      integer, dimension(:), intent(out) :: atype_ch              !< Actual type of atom for dilute system
      integer, dimension(:), intent(out) :: asite_ch              !< Actual site of atom for dilute system
      integer, dimension(:), intent(out) :: achem_ch              !< Chemical type of atoms (reduced list)
      integer, dimension(:), intent(out) :: acellnumb             !< List for translating atom no. in full cell to actual cell
      integer, dimension(:), intent(out) :: acellnumbrev          !< List for translating atom no. in actual cell to full cell
      ! .. Local variables
      integer :: iatom,iclus,i,count
      character(len=20) :: filn

      write (*,'(2x,a)',advance='no') 'Set up types of atoms'
      call setup_type_and_numb(Natom_full,NA,N1,N2,N3,atype,anumb,atype_inp,        &
         anumb_inp,block_size)
      write (*,'(a)') ' done'

      write (*,'(2x,a)',advance='no') 'Set up chemical information of alloy'
      if(do_ralloy==1) then
         call setup_chemicaldata(Natom,NA,N1,N2,N3,atype,tseed,do_ralloy,Natom_full,&
            Nchmax,Nch,achtype,acellnumb,acellnumbrev,chconc,atype_ch,asite_ch,     &
            achem_ch)
         ! Increase the number of atoms in the system by considering the ones outside
         ! the magnetic layers
         Natom=Natom+clus_expand
      else
         Natom=Natom_full
      end if

      write (*,'(a)') ' done'
      write (*,'(2x,a)',advance='no') 'Set up global coordinates'
      call setup_globcoord(Natom,NA,Bas,C1,C2,C3,coord,N1,N2,N3,atype,anumb,        &
         do_prnstruct,do_prn_poscar,simid,do_ralloy,Natom_full,acellnumb,achtype,atype_ch,        &
         asite_ch,achem_ch,block_size)
         ! Overwrite the coordinates of the atoms belonging to the cluster, this mostly
         ! takes care of the atoms outside the magnetic layer
         if (clus_expand>0) then
            ! Open file for writing the coordinates
            if(do_prnstruct==1.or.do_prnstruct==2.or.do_prnstruct==4) then
               write (filn,'(''coord.'',a8,''.out'')') simid
               open(ofileno, file=filn,position='append')
            end if
            ! Loop over the atoms in the cluster
            count=1
            do iclus=1,Natom_clus
               iatom=index_clus(iclus)
               coord(:,iatom)=coord_clus(:,iclus)
               if (iatom>N1*N2*N3*NA) then
                  atype(iatom)=1
                  anumb(iatom)=1
                  count=count+1
                  if (do_ralloy==0) then
                     if(do_prnstruct==1.or.do_prnstruct==2.or.do_prnstruct==4) then
                        write(ofileno,'(i7,3f12.6,2i6)') iatom,coord(1:3,iatom),       &
                           atype(iatom),anumb(iatom)
                     endif
                  else
                     i=acellnumb(iatom)
                     if(do_prnstruct==1.or.do_prnstruct==2.or.do_prnstruct==4) then
                        write(ofileno,'(i7,3f12.6,4i7)') iatom,coord(1:3,i),           &
                           atype_ch(i),asite_ch(i),achem_ch(i)
                     end if
                  endif
               endif
            enddo
            ! Close the geometry file
            if(do_prnstruct==1.or.do_prnstruct==4) then
               close(ofileno)
            end if
         endif

      write (*,'(a)') ' done'

   end subroutine modify_global_coordinates

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_type_and_numb
   !> @brief Sets up the type of atoms in the system
   !----------------------------------------------------------------------------
   subroutine setup_type_and_numb_clus(Natom_full_clus,NA_clus,N1_clus,N2_clus,     &
      N3_clus,atype_clus,anumb_clus,atype_inp_clus,anumb_inp_clus,block_size)
      !
      implicit none

      integer, intent(in) :: NA_clus  !< Number of atoms in the cluster unit cell
      integer, intent(in) :: N1_clus  !< Number of cell repetitions in x direction for the cluster
      integer, intent(in) :: N2_clus  !< Number of cell repetitions in y direction for the cluster
      integer, intent(in) :: N3_clus  !< Number of cell repetitions in z direction for the cluster
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      integer, intent(in) :: Natom_full_clus !< Number of atoms in system
      ! .. In/Out variables
      integer, dimension(Natom_full_clus), intent(inout) :: atype_clus !< Type of atom
      integer, dimension(Natom_full_clus), intent(inout) :: anumb_clus !< Atom number in cell
      integer, dimension(NA_clus), intent(inout)         :: atype_inp_clus  !< Type of atom from input
      integer, dimension(NA_clus), intent(inout)         :: anumb_inp_clus !< Atom number in cell from input
      !
      integer :: A1
      integer :: I0, I1, I2, I3, II1, II2, II3
      !
      A1=0
      do iI3=0, N3_clus-1, block_size
         do iI2=0, N2_clus-1, block_size
            do iI1=0, N1_clus-1, block_size
               do I3=ii3,min(ii3+block_size-1,N3_clus-1)
                  do I2=ii2,min(ii2+block_size-1,N2_clus-1)
                     do I1=ii1,min(ii1+block_size-1,N1_clus-1)
                        do I0=1, NA_clus
                           A1=A1+1
                           atype_clus(A1)=atype_inp_clus(I0)
                           anumb_clus(A1)=anumb_inp_clus(I0)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine setup_type_and_numb_clus

   !----------------------------------------------------------------------------
   !> @brief Sets up chemical information of the system in case of random alloy for the cluster
   !> @note IMPORTANT: In case of dilute system, Natom variable is reduced.
   !> @note Based on the routine for the host
   !----------------------------------------------------------------------------
   subroutine setup_chemicaldata_clus(tseed,NA_clus,N1_clus,N2_clus,N3_clus,        &
      do_ralloy,Nchmax_clus,Natom_full_clus,Nch_clus,chconc_clus,achtype_clus,      &
      atype_ch_clus,asite_ch_clus,achem_ch_clus,acellnumb_clus,acellnumbrev_clus,   &
      Natom_clus,atype_clus)
      !
      use RandomNumbers, only : rng_uniform
      use Sorting, only : mergesort
      !
      implicit none

      integer, intent(in) :: tseed    !< Temperature seed
      integer, intent(in) :: NA_clus  !< Number of atoms in the cluster unit cell
      integer, intent(in) :: N1_clus  !< Number of cell repetitions in x direction for the cluster
      integer, intent(in) :: N2_clus  !< Number of cell repetitions in y direction for the cluster
      integer, intent(in) :: N3_clus  !< Number of cell repetitions in z direction for the cluster
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: Nchmax_clus  !< Max number of chemical components on each site in cell
      integer, intent(in) :: Natom_full_clus !< Number of atoms in system
      integer, dimension(NA_clus), intent(in) :: Nch_clus !< Number of chemical components on each site in cell
      real(dblprec), dimension(NA_clus,Nchmax_clus), intent(in) :: chconc_clus !< Chemical concentration on sites
      ! .. Output variables
      integer, dimension(Natom_full_clus), intent(out) :: achtype_clus      !< Chemical type of atoms (full list)
      integer, dimension(Natom_full_clus), intent(out) :: atype_ch_clus     !< Actual site of atom for dilute system
      integer, dimension(Natom_full_clus), intent(out) :: asite_ch_clus     !< Actual site of atom for dilute system
      integer, dimension(Natom_full_clus), intent(out) :: achem_ch_clus     !< Chemical type of atoms (reduced list)
      integer, dimension(Natom_full_clus), intent(out) :: acellnumb_clus    !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full_clus), intent(out) :: acellnumbrev_clus !< List for translating atom no. in actual cell to full cell

      ! .. In/Out variables
      integer, intent(inout) :: Natom_clus   !< Number of atoms in the cluster
      integer, dimension(Natom_full_clus), intent(inout) :: atype_clus !< Type of atom for the cluster

      integer :: A1
      integer :: i0, i1, i2, i3
      integer :: Ncell, iatom
      integer :: ia, ich, i
      integer :: i_stat, i_all
      integer, dimension(1) :: irn
      integer, dimension(:,:), allocatable :: atoms
      integer, dimension(:,:), allocatable :: qch
      integer, dimension(:), allocatable :: ns, ne
      integer, dimension(:), allocatable :: atoms2, atoms2T
      real(dblprec), dimension(:), allocatable :: rn

      if (do_ralloy==1) then

         ! Same seed as for Monte-Carlo
         Ncell = N1_clus*N2_clus*N3_clus
         allocate(atoms(Ncell,Na_clus),stat=i_stat)
         call memocc(i_stat,product(shape(atoms))*kind(atoms),'atoms','setup_chemicaldata_clus')
         allocate(atoms2(Ncell),stat=i_stat)
         call memocc(i_stat,product(shape(atoms2))*kind(atoms2),'atoms2','setup_chemicaldata_clus')
         allocate(atoms2T(Ncell),stat=i_stat)
         call memocc(i_stat,product(shape(atoms2T))*kind(atoms2T),'atoms2T','setup_chemicaldata_clus')
         allocate(qch(Nchmax_clus,Na_clus),stat=i_stat)
         call memocc(i_stat,product(shape(qch))*kind(qch),'qch','setup_chemicaldata_clus')
         allocate(rn(Ncell),stat=i_stat)
         call memocc(i_stat,product(shape(rn))*kind(rn),'rn','setup_chemicaldata_clus')
         allocate(ns(Nchmax_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ns))*kind(ns),'ns','setup_chemicaldata_clus')
         allocate(ne(Nchmax_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ne))*kind(ne),'ne','setup_chemicaldata_clus')
         qch=0
         atoms=0
         achtype_clus=0

         do ia=1,Na_clus
            do ich=1,Nch_clus(ia)
               qch(ich,ia) = nint(chconc_clus(ia,ich)*Ncell)
            end do
         end do

         do ia=1,Na_clus
            call rng_uniform(rn,Ncell)
            !Partitioning of the Ncell random nr:s
            !to sets for each chemical type ich
            ns(1) = 1
            ne(1) = qch(1,ia)
            do ich=2,Nch_clus(ia)
               ns(ich) = ns(ich-1)+qch(ich-1,ia)
               ne(ich) = ne(ich-1)+qch(ich,ia)
            end do
            do i=1,Ncell
               irn = maxloc(rn)
               atoms(i,ia)=irn(1)
               rn(irn)=0.0_dblprec
            end do
            do ich=1,Nch_clus(ia)
               do i=ns(ich),ne(ich)
                  A1=(atoms(i,ia)-1)*Na_clus+ia
                  achtype_clus(A1)=ich
               end do
            end do
         end do
         acellnumb_clus=0
         acellnumbrev_clus=0
         iatom=0
         ! asite_ch and achem_ch contains information of species on each site
         do I3=0, N3_clus-1
            do I2=0, N2_clus-1
               do I1=0, N1_clus-1
                  do I0=1, NA_clus
                     i=I0+I1*NA_clus+I2*N1_clus*NA_clus+I3*N2_clus*N1_clus*NA_clus
                     if (achtype_clus(i) /= 0) then
                        iatom = iatom + 1
                        acellnumb_clus(i) = iatom
                        acellnumbrev_clus(iatom) = i
                        asite_ch_clus(iatom)=I0
                        atype_ch_clus(iatom)=atype_clus(i)
                        achem_ch_clus(iatom)=achtype_clus(i)
                     end if
                  end do
               end do
            end do
         end do

         ! Natom is reduced in case of dilute system
         Natom_clus = iatom
         ! Deallocate temporary arrays
         i_all=-product(shape(atoms))*kind(atoms)
         deallocate(atoms,stat=i_stat)
         call memocc(i_stat,i_all,'atoms','setup_chemicaldata_clus')
         i_all=-product(shape(atoms2))*kind(atoms2)
         deallocate(atoms2,stat=i_stat)
         call memocc(i_stat,i_all,'atoms2','setup_chemicaldata_clus')
         i_all=-product(shape(atoms2T))*kind(atoms2T)
         deallocate(atoms2T,stat=i_stat)
         call memocc(i_stat,i_all,'atoms2T','setup_chemicaldata_clus')
         i_all=-product(shape(qch))*kind(qch)
         deallocate(qch,stat=i_stat)
         call memocc(i_stat,i_all,'qch','setup_chemicaldata_clus')
         i_all=-product(shape(rn))*kind(rn)
         deallocate(rn,stat=i_stat)
         call memocc(i_stat,i_all,'rn','setup_chemicaldata_clus')
         i_all=-product(shape(ns))*kind(ns)
         deallocate(ns,stat=i_stat)
         call memocc(i_stat,i_all,'ns','setup_chemicaldata_clus')
         i_all=-product(shape(ne))*kind(ne)
         deallocate(ne,stat=i_stat)
         call memocc(i_stat,i_all,'ne','setup_chemicaldata_clus')
      end if

   end subroutine setup_chemicaldata_clus

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_clusdata
   !> @brief Routine to allocate the crystalline data for the cluster
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_clusdata(Natom_clus, flag)

      implicit none

      integer, intent(in),optional :: Natom_clus !< Number of atoms in the cluster
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_stat, i_all

      if(flag>0) then
         allocate(anumb_clus(Natom_clus),stat=i_stat)
         call memocc(i_stat,product(shape(anumb_clus))*kind(anumb_clus),'anumb_clus','allocate_clusdata')
         allocate(atype_clus(Natom_clus),stat=i_stat)
         call memocc(i_stat,product(shape(atype_clus))*kind(atype_clus),'atype_clus','allocate_clusdata')
      else
         i_all=-product(shape(anumb_clus))*kind(anumb_clus)
         deallocate(anumb_clus,stat=i_stat)
         call memocc(i_stat,i_all,'anumb_clus','allocate_clusdata')
         i_all=-product(shape(atype_clus))*kind(atype_clus)
         deallocate(atype_clus,stat=i_stat)
         call memocc(i_stat,i_all,'atype_clus','allocate_clusdata')
         i_all=-product(shape(index_clus))*kind(index_clus)
         deallocate(index_clus,stat=i_stat)
         call memocc(i_stat,i_all,'index_clus','allocate_clusdata')
         if (allocated(coord_clus)) then
            i_all=-product(shape(coord_clus))*kind(coord_clus)
            deallocate(coord_clus,stat=i_stat)
            call memocc(i_stat,i_all,'coord_clus','allocate_clusdata')
         endif
      end if

   end subroutine allocate_clusdata

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_chemicaldata_clus
   !> Allocate arrays for random alloys for the cluster
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_chemicaldata_clus(Natom_clus,flag)
      implicit none

      integer, optional, intent(in) :: Natom_clus !< Number of atoms in the cluster
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(achtype_clus(Natom_clus),stat=i_stat)
         call memocc(i_stat,product(shape(achtype_clus))*kind(achtype_clus),'achtype_clus','allocate_chemicaldata_clus')
         achtype_clus=0
         allocate(acellnumb_clus(Natom_clus),stat=i_stat)
         call memocc(i_stat,product(shape(acellnumb_clus))*kind(acellnumb_clus),'acellnumb_clus','allocate_chemicaldata_clus')
         acellnumb_clus=0
         allocate(acellnumbrev_clus(Natom_clus),stat=i_stat)
         call memocc(i_stat,product(shape(acellnumbrev_clus))*kind(acellnumbrev_clus),'acellnumbrev_clus','allocate_chemicaldata_clus')
         acellnumbrev_clus=0
         allocate(achem_ch_clus(Natom_clus),stat=i_stat)
         call memocc(i_stat,product(shape(achem_ch_clus))*kind(achem_ch_clus),'achem_ch_clus','allocate_chemicaldata_clus')
         achem_ch_clus=0
         allocate(asite_ch_clus(Natom_clus),stat=i_stat)
         call memocc(i_stat,product(shape(asite_ch_clus))*kind(asite_ch_clus),'asite_ch_clus','allocate_chemicaldata_clus')
         asite_ch_clus=0
         allocate(atype_ch_clus(Natom_clus),stat=i_stat)
         call memocc(i_stat,product(shape(atype_ch_clus))*kind(atype_ch_clus),'atype_ch','allocate_chemicaldata_clus')
         atype_ch_clus=0
      else
         i_all=-product(shape(achtype_clus))*kind(achtype_clus)
         deallocate(achtype_clus,stat=i_stat)
         call memocc(i_stat,i_all,'achtype_clus','allocate_chemicaldata_clus')
         i_all=-product(shape(acellnumb_clus))*kind(acellnumb_clus)
         deallocate(acellnumb_clus,stat=i_stat)
         call memocc(i_stat,i_all,'acellnumb_clus','allocate_chemicaldata_clus')
         i_all=-product(shape(acellnumbrev_clus))*kind(acellnumbrev_clus)
         deallocate(acellnumbrev_clus,stat=i_stat)
         call memocc(i_stat,i_all,'acellnumbrev_clus','allocate_chemicaldata_clus')
         i_all=-product(shape(achem_ch_clus))*kind(achem_ch_clus)
         deallocate(achem_ch_clus,stat=i_stat)
         call memocc(i_stat,i_all,'achem_ch_clus','allocate_chemicaldata_clus')
         i_all=-product(shape(asite_ch_clus))*kind(asite_ch_clus)
         deallocate(asite_ch_clus,stat=i_stat)
         call memocc(i_stat,i_all,'asite_ch_clus','allocate_chemicaldata_clus')
         i_all=-product(shape(atype_ch_clus))*kind(atype_ch_clus)
         deallocate(atype_ch_clus,stat=i_stat)
         call memocc(i_stat,i_all,'atype_ch_clus','allocate_chemicaldata_clus')
      end if

   end subroutine allocate_chemicaldata_clus

end module clus_geometry
