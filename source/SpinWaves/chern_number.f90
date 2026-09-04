!------------------------------------------------------------------------------------
!> @brief
!> Routine used to calculate the Chern number of the bands calculated using LSWT
!> @author
!> Manuel Pereiro 
!> Nastaran Salehi
!> @copyright
!> GNU Public License
!------------------------------------------------------------------------------------
module Chern_number
   ! use, intrinsic :: ieee_arithmetic,  only : ieee_is_nan
   use InputData,  only : Temp
   use Parameters
   use Constants
   use Profiling
   use Hamiltoniandata,    only : ham
   use InputData,   only : ham_inp
  use Diamag ,     only : clone_q,diagonalize_quad_hamiltonian,find_uv,setup_ektij,&
                   setup_jtens2_q,setup_jtens_q,sJs, setup_tensor_hamiltonian,&
                   nc_eval_complex,nc_evec_complex,boson_overlap
   !
   implicit none
   !
   character(len=1)                           :: do_chern    !< Calculate the Chern number of the bands (Y/N)
   character(len=1)                           :: do_magnon_oam !< Fishman reciprocal-space magnon OAM (Y/N)
   integer                                    :: Nx          !< Number of points of the grid in x direction
   integer                                    :: Ny          !< Number of points of the grid in y direction
   integer                                    :: Nz          !< Number of points of the grid in z direction
   real(dblprec), dimension(3)                :: Chern_qvect !< Spin spiral ordering vector
   !
   private
   ! public subroutines
   public :: do_chern,do_magnon_oam,Nx,Ny,Nz, Chern_qvect
   public :: read_parameters_chern_number,calculate_chern_number
   !
contains

   subroutine setup_chern_number()

      implicit none

      do_chern    = 'N'
      do_magnon_oam = 'N'
      Nx          = 100
      Ny          = 100
      Nz          = 1
      Chern_qvect = 0.0_dblprec


   end subroutine setup_chern_number

   subroutine calculate_chern_number(NA,Natom,Mensemble,simid,emomM,mmom,Nx,Ny,Nz,C1,C2,C3)
      ! Calculate the Chern number of the bands in the 1st BZ.
      !
      implicit none
      !
      character(LEN = 25) :: bphase_file
      character(LEN = 25) :: chern_file
      character(LEN = 25) :: oam_file
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment magnitude
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current magnetic moment vector
      integer, intent(in) :: Nx  !< Number of points of the grid in x direction
      integer, intent(in) :: Ny  !< Number of points of the grid in y direction
      integer, intent(in) :: Nz  !< Number of points of the grid in z direction
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      !
      integer                                         :: dimen                          !< Number of q-vectors
      complex(dblprec), dimension(:,:), allocatable   :: u1,u2,u3                       !< Link variable between eigenvectors
      complex(dblprec), dimension(:,:), allocatable   :: u1inv,u2inv,u3inv              !< Link variable between inverse eigenvectors
      complex(dblprec), dimension(:,:), allocatable   :: u1norm,u2norm,u3norm           !< Link unit vectors
      complex(dblprec), dimension(:,:), allocatable   :: u1invnorm,u2invnorm,u3invnorm  !< Link inverse unit vectors
      complex(dblprec), dimension(:,:), allocatable   :: Berry_cuv                      !< Berry curvature
      real(dblprec), dimension(:,:), allocatable      :: rho                            !< Bose-Einstein Distribution
      real(dblprec), dimension(:,:), allocatable      :: c2_func                        !< Function c2
      real(dblprec), dimension(:), allocatable        :: Ch_numq                        !< Chern number for the phason bands
      real(dblprec), dimension(:), allocatable        :: Ch_numqplus                    !< Chern number for the +q bands
      real(dblprec), dimension(:), allocatable        :: Ch_numqminus                   !< Chern number for the -q bands
      real(dblprec), dimension(:,:), allocatable      :: q_vchern                       !< q-vectors
      real(dblprec), dimension(:), allocatable        :: therm_conduc_band              !< Thermal conductivity per band
      integer, dimension(:), allocatable              :: indx                           !< index along x
      integer, dimension(:), allocatable              :: indy                           !< index along y
      integer, dimension(:), allocatable              :: indz                           !< index along z
      ! Fishman OAM, stored band- and k-resolved in units of hbar.
      real(dblprec), dimension(:,:), allocatable      :: oam_hbar
      complex(dblprec), dimension(:,:,:), allocatable  :: oam_evec
      integer, dimension(:,:), allocatable             :: oam_band

      real(dblprec)                                   :: therm_conduc              !< Thermal conductivity
      integer                                         :: iq,i,j,k,l,m,i_stat,nqred,icount,jcount,kcount,counter,kx,ky,kz, nmx
      ! silent NaN-cleanup before conductivity summation
      !
      print '(1x,a)', 'Calculating Chern numbers'
      !Defining variables kx,ky,kz
      kx=Nx
      ky=Ny
      kz=Nz
      !Definitions for 1d systems
      if (Nx .eq. 1 .and. Ny .eq. 1) then
      kx=2
      ky=2
      else if (Nx .eq. 1 .and. Nz .eq. 1) then
      kx=2
      kz=2
      else if (Ny .eq. 1 .and. Nz .eq.1 ) then
      ky=2
      kz=2
      end if
      !Definitions for 2d systems
      if (Nx .eq. 1) then
      kx=2
      else if (Ny .eq. 1) then
      ky=2
      else if (Nz .eq. 1) then
      kz=2
      end if
      !Defining  qred and size of the grid
      nqred=((kx-1)*(ky-1)*(kz-1))
      dimen=Nx*Ny*Nz
      !Allocate variables
      allocate(indx(3*dimen),stat=i_stat)
      call memocc(i_stat,product(shape(indx))*kind(indx),'indx','calculate_chern_number')
      allocate(indy(3*dimen),stat=i_stat)
      call memocc(i_stat,product(shape(indy))*kind(indy),'indy','calculate_chern_number')
      allocate(indz(3*dimen),stat=i_stat)
      call memocc(i_stat,product(shape(indz))*kind(indz),'indz','calculate_chern_number')
      allocate(u1(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u1))*kind(u1),'u1','calculate_chern_number')
      allocate(u1norm(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u1norm))*kind(u1norm),'u1norm','calculate_chern_number')
      allocate(u2(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u2))*kind(u2),'u2','calculate_chern_number')
      allocate(u2norm(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u2norm))*kind(u2norm),'u2norm','calculate_chern_number')
      allocate(u3(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u3))*kind(u3),'u3','calculate_chern_number')
      allocate(u3norm(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u3norm))*kind(u3norm),'u3norm','calculate_chern_number')
      allocate(u1inv(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u1inv))*kind(u1inv),'u1inv','calculate_chern_number')
      allocate(u2inv(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u2inv))*kind(u2inv),'u2inv','calculate_chern_number')
      allocate(u3inv(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u3inv))*kind(u3inv),'u3inv','calculate_chern_number')
      allocate(u1invnorm(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u1invnorm))*kind(u1invnorm),'u1invnorm','calculate_chern_number')
      allocate(u2invnorm(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u2invnorm))*kind(u2invnorm),'u2invnorm','calculate_chern_number')
      allocate(u3invnorm(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(u3invnorm))*kind(u3invnorm),'u3invnorm','calculate_chern_number')
      allocate(Berry_cuv(NA,3*nqred),stat=i_stat)
      call memocc(i_stat,product(shape(Berry_cuv))*kind(Berry_cuv),'Berry_cuv','calculate_chern_number')
      allocate(Ch_numq(NA),stat=i_stat)
      call memocc(i_stat,product(shape(Ch_numq))*kind(Ch_numq),'Ch_numq','calculate_chern_number')
      allocate(Ch_numqplus(NA),stat=i_stat)
      call memocc(i_stat,product(shape(Ch_numqplus))*kind(Ch_numqplus),'Ch_numqplus','calculate_chern_number')
      allocate(Ch_numqminus(NA),stat=i_stat)
      call memocc(i_stat,product(shape(Ch_numqminus))*kind(Ch_numqminus),'Ch_numqminus','calculate_chern_number')
      allocate(q_vchern(3,dimen),stat=i_stat)
      call memocc(i_stat,product(shape(q_vchern))*kind(q_vchern),'q_vchern','calculate_chern_number')
      allocate(rho(2*NA,dimen*3),stat=i_stat)
      call memocc(i_stat,product(shape(rho))*kind(rho),'rho','calculate_chern_number')
      allocate(c2_func(2*NA,dimen*3),stat=i_stat)
      call memocc(i_stat,product(shape(c2_func))*kind(c2_func),'c2_func','calculate_chern_number')
      allocate(therm_conduc_band(2*NA),stat=i_stat)
      call memocc(i_stat,product(shape(therm_conduc_band))*kind(therm_conduc_band),'therm_conduc_band','calculate_chern_number')
      !Calculate the grid in the reciprocal space
      call setup_grid(Nx,Ny,Nz,C1,C2,C3,dimen,q_vchern)
      !Calculate eigenvectors and eigenvalues
      call setup_tensor_hamiltonian(NA,Natom, Mensemble, simid, emomM, mmom, q_vchern, dimen,1)
      !Initialize variables
      j=0
      k=0
      l=0
      !Counter for phason mode
      counter=0
      do kcount = 1,Nz
        do jcount= 1,Ny
          do icount= 1,Nx
             counter=counter+1
             if (icount.eq.Nx .and. Nx.ne.1) then
             indx(counter)= 0
             else
             indx(counter)= 1
             end if
             if (jcount.eq.Ny .and. Ny.ne.1) then
             indy(counter)= 0
             else
             indy(counter)= 1
             end if
             if (kcount.eq.Nz .and. Nz.ne.1) then
             indz(counter)= 0
             else
             indz(counter)= 1
             end if
          end do
        end do
      end do
      !Counter for +q
      do kcount = 1,Nz
        do jcount= 1,Ny
          do icount= 1,Nx
             counter=counter+1
             if (icount.eq.Nx .and. Nx.ne.1) then
             indx(counter)= 0
             else
             indx(counter)= 2
             end if
             if (jcount.eq.Ny .and. Ny.ne.1) then
             indy(counter)= 0
             else
             indy(counter)= 2
             end if
             if (kcount.eq.Nz .and. Nz.ne.1) then
             indz(counter)= 0
             else
             indz(counter)= 2
             end if
          end do
        end do
      end do
      !Counter for -q
      do kcount = 1,Nz
        do jcount= 1,Ny
          do icount= 1,Nx
             counter=counter+1
             if (icount.eq.Nx .and. Nx.ne.1) then
             indx(counter)= 0
             else
             indx(counter)= 3
             end if
             if (jcount.eq.Ny .and. Ny.ne.1) then
             indy(counter)= 0
             else
             indy(counter)= 3
             end if
             if (kcount.eq.Nz .and. Nz.ne.1) then
             indz(counter)= 0
             else
             indz(counter)= 3
             end if
          end do
        end do
      end do
      !1d grid
      if (Nx .eq. 1 .and. Ny .eq. 1) then
        kx=Ny
        ky=0
        kz=0
      else if (Nx .eq. 1 .and. Nz .eq. 1) then
        kx=Ny
        kz=0
        ky=0
      else if (Ny .eq. 1 .and. Nz .eq.1 ) then
        kx=Nx
        ky=0
        kz=0
      end if
      ! 2D grid and 3D grid
      if (Nz .eq. 1) then
        kx=Nx
        ky=Ny
        kz=0
      else if (Ny .eq. 1) then
        kx=Nx
        ky=Nz
        kz=0
      else if (Nx .eq. 1) then
        kx=Ny
        ky=Nz
        kz=0
      else
        kx=Nx
        ky=Ny
        kz=Nx*Ny
      end if
      !Calculating the link variables
      do iq=1, 3*dimen
          if ( indx(iq).eq.1 .and. indy(iq).eq.1 .and. indz(iq).eq.1 ) then
            j=j+1
            do i=1,NA !band index
                u1(i,j)=boson_overlap(nc_evec_complex(:,i,iq),nc_evec_complex(:,i,iq+1),NA)
              if (abs(u1(i,j))== 0.0_dblprec) then
                u1norm(i,j)=(1.0_dblprec,0.0_dblprec)
              else
                u1norm(i,j)=u1(i,j)/abs(u1(i,j))
              end if

                u2(i,j)=boson_overlap(nc_evec_complex(:,i,iq+1),nc_evec_complex(:,i,iq+1+kx),NA)
              if (abs(u2(i,j))== 0.0_dblprec) then
                u2norm(i,j)=(1.0_dblprec,0.0_dblprec)
              else
                u2norm(i,j)=u2(i,j)/abs(u2(i,j))
              end if

                u3(i,j)=boson_overlap(nc_evec_complex(:,i,iq+1+kx),nc_evec_complex(:,i,iq+1+kx+kz),NA)
              if (abs(u3(i,j))== 0.0_dblprec) then
                u3norm(i,j)=(1.0_dblprec,0.0_dblprec)
              else
                u3norm(i,j)=u3(i,j)/abs(u3(i,j))
              end if

                u1inv(i,j)=1.0_dblprec/boson_overlap(nc_evec_complex(:,i,iq+kx+kz),nc_evec_complex(:,i,iq+1+kx+kz),NA)
              if (abs(u1inv(i,j))== 0.0_dblprec) then
                u1invnorm(i,j)=(1.0_dblprec,0.0_dblprec)
              else
                u1invnorm(i,j)=u1inv(i,j)/abs(u1inv(i,j))
              end if

                u2inv(i,j)=1.0_dblprec/boson_overlap(nc_evec_complex(:,i,iq+kz),nc_evec_complex(:,i,iq+kx+kz),NA)
              if (abs(u2inv(i,j))== 0.0_dblprec) then
                u2invnorm(i,j)=(1.0_dblprec,0.0_dblprec)
              else
                u2invnorm(i,j)=u2inv(i,j)/abs(u2inv(i,j))
              end if

                u3inv(i,j)=1.0_dblprec/boson_overlap(nc_evec_complex(:,i,iq),nc_evec_complex(:,i,iq+kz),NA)
              if (abs(u3inv(i,j))== 0.0_dblprec) then
                u3invnorm(i,j)=(1.0_dblprec,0.0_dblprec)
              else
                u3invnorm(i,j)=u3inv(i,j)/abs(u3inv(i,j))
              end if
            end do
          end if

          if ( indx(iq).eq.2 .and. indy(iq).eq.2 .and. indz(iq).eq.2 ) then
            k=k+1
            do i=1,NA !band index
              u1(i,j+k)=boson_overlap(nc_evec_complex(:,i,iq),nc_evec_complex(:,i,iq+1),NA)
              if (abs(u1(i,j+k))== 0.0_dblprec) then
              u1norm(i,j+k)=(1.0_dblprec,0.0_dblprec)
              else
              u1norm(i,j+k)=u1(i,j+k)/abs(u1(i,j+k))
              end if

              u2(i,j+k)=boson_overlap(nc_evec_complex(:,i,iq+1),nc_evec_complex(:,i,iq+1+kx),NA)
              if (abs(u2(i,j+k))== 0.0_dblprec) then
              u2norm(i,j+k)=(1.0_dblprec,0.0_dblprec)
              else
              u2norm(i,j+k)=u2(i,j+k)/abs(u2(i,j+k))
              end if

              u3(i,j+k)=boson_overlap(nc_evec_complex(:,i,iq+1+kx),nc_evec_complex(:,i,iq+1+kx+kz),NA)
              if (abs(u3(i,j+k))== 0.0_dblprec) then
                u3norm(i,j+k)=(1.0_dblprec,0.0_dblprec)
              else
                u3norm(i,j+k)=u3(i,j+k)/abs(u3(i,j+k))
              end if

              u1inv(i,j+k)=1.0_dblprec/boson_overlap(nc_evec_complex(:,i,iq+kx+kz),nc_evec_complex(:,i,iq+1+kx+kz),NA)
              if (abs(u1inv(i,j+k))== 0.0_dblprec) then
              u1invnorm(i,j+k)=(1.0_dblprec,0.0_dblprec)
              else
              u1invnorm(i,j+k)=u1inv(i,j+k)/abs(u1inv(i,j+k))
              end if

              u2inv(i,j+k)=1.0_dblprec/boson_overlap(nc_evec_complex(:,i,iq+kz),nc_evec_complex(:,i,iq+kx+kz),NA)
              if (abs(u2inv(i,j+k))== 0.0_dblprec) then
              u2invnorm(i,j+k)=(1.0_dblprec,0.0_dblprec)
              else
              u2invnorm(i,j+k)=u2inv(i,j+k)/abs(u2inv(i,j+k))
              end if

              u3inv(i,j+k)=1.0_dblprec/boson_overlap(nc_evec_complex(:,i,iq),nc_evec_complex(:,i,iq+kz),NA)
              if (abs(u3inv(i,j+k))== 0.0_dblprec) then
                u3invnorm(i,j+k)=(1.0_dblprec,0.0_dblprec)
              else
                u3invnorm(i,j+k)=u3inv(i,j+k)/abs(u3inv(i,j+k))
              end if
            end do
          end if

          if ( indx(iq).eq.3 .and. indy(iq).eq.3 .and. indz(iq).eq.3 ) then
            l=l+1
            do i=1,NA !band index
              u1(i,j+k+l)=boson_overlap(nc_evec_complex(:,i,iq),nc_evec_complex(:,i,iq+1),NA)
              if (abs(u1(i,j+k+l))== 0.0_dblprec) then
              u1norm(i,j+k+l)=(1.0_dblprec,0.0_dblprec)
              else
              u1norm(i,j+k+l)=u1(i,j+k+l)/abs(u1(i,j+k+l))
              end if

              u2(i,j+k+l)=boson_overlap(nc_evec_complex(:,i,iq+1),nc_evec_complex(:,i,iq+1+kx),NA)
              if (abs(u2(i,j+k+l))== 0.0_dblprec) then
              u2norm(i,j+k+l)=(1.0_dblprec,0.0_dblprec)
              else
              u2norm(i,j+k+l)=u2(i,j+k+l)/abs(u2(i,j+k+l))
              end if

              u3(i,j+k+l)=boson_overlap(nc_evec_complex(:,i,iq+1+kx),nc_evec_complex(:,i,iq+1+kx+kz),NA)
              if (abs(u3(i,j+k+l))== 0.0_dblprec) then
                u3norm(i,j+k+l)=(1.0_dblprec,0.0_dblprec)
              else
                u3norm(i,j+k+l)=u3(i,j+k+l)/abs(u3(i,j+k+l))
              end if

              u1inv(i,j+k+l)=1.0_dblprec/boson_overlap(nc_evec_complex(:,i,iq+kx+kz),nc_evec_complex(:,i,iq+1+kx+kz),NA)
              if (abs(u1inv(i,j+k+l))== 0.0_dblprec) then
              u1invnorm(i,j+k+l)=(1.0_dblprec,0.0_dblprec)
              else
              u1invnorm(i,j+k+l)=u1inv(i,j+k+l)/abs(u1inv(i,j+k+l))
              end if

              u2inv(i,j+k+l)=1.0_dblprec/boson_overlap(nc_evec_complex(:,i,iq+kz),nc_evec_complex(:,i,iq+kx+kz),NA)
              if (abs(u2inv(i,j+k+l))== 0.0_dblprec) then
              u2invnorm(i,j+k+l)=(1.0_dblprec,0.0_dblprec)
              else
              u2invnorm(i,j+k+l)=u2inv(i,j+k+l)/abs(u2inv(i,j+k+l))
              end if

              u3inv(i,j+k+l)=1.0_dblprec/boson_overlap(nc_evec_complex(:,i,iq),nc_evec_complex(:,i,iq+kz),NA)
              if (abs(u3inv(i,j+k+l))== 0.0_dblprec) then
                u3invnorm(i,j+k+l)=(1.0_dblprec,0.0_dblprec)
              else
                u3invnorm(i,j+k+l)=u3inv(i,j+k+l)/abs(u3inv(i,j+k+l))
              end if
            end do
          end if
      enddo

      Berry_cuv=log(u1norm*u2norm*u3norm*u1invnorm*u2invnorm*u3invnorm)

      if (do_magnon_oam == 'Y') then
         ! Chern number uses Berry flux.  Fishman OAM is a separate,
         ! gauge-fixed expectation value of the momentum-space angular
         ! momentum operator and is kept k- and band-resolved.
         allocate(oam_evec(2*NA,NA,dimen),stat=i_stat)
         call memocc(i_stat,product(shape(oam_evec))*kind(oam_evec),'oam_evec','calculate_chern_number')
         allocate(oam_hbar(NA,dimen),stat=i_stat)
         call memocc(i_stat,product(shape(oam_hbar))*kind(oam_hbar),'oam_hbar','calculate_chern_number')
         allocate(oam_band(NA,dimen),stat=i_stat)
         call memocc(i_stat,product(shape(oam_band))*kind(oam_band),'oam_band','calculate_chern_number')
         call prepare_oam_eigenvectors(NA,Nx,Ny,Nz,dimen,nc_evec_complex,oam_evec,oam_band)
         call calculate_fishman_oam(NA,Nx,Ny,Nz,dimen,q_vchern,oam_evec,oam_hbar)

         oam_file='oam.'//trim(simid)//'.out'
         open(ofileno,file=oam_file)
         write(ofileno,'(a)') '# Fishman reciprocal-space magnon OAM; primary units OAM/hbar'
         write(ofileno,'(a)') '# Pointwise values are gauge-fixed and gauge-dependent.'
         write(ofileno,'(a)') '# qx qy are Cartesian components of q_vchern; physical k=2*pi*q.'
         write(ofileno,'(a)') '# band qx qy energy(meV) OAM/hbar'
         do iq=1,dimen
            do i=1,NA
               write(ofileno,'(i6,4(1x,es23.15))') i,q_vchern(1,iq),q_vchern(2,iq), &
                  nc_eval_complex(oam_band(i,iq),iq),oam_hbar(i,iq)
            end do
         end do
         close(ofileno)
         print '(1x,a,a)', 'Fishman OAM written to ',trim(oam_file)
         print '(1x,a)', 'Pointwise OAM is gauge-fixed; angular-average TODO for polar meshes.'
      end if

      ! 2D grid or 3D grid
      if (kz == 0) then
        Ch_numq=(1.0_dblprec/(2*pi))*aimag(sum(Berry_cuv(:,1:j),dim=2))
        Ch_numqplus=(1.0_dblprec/(2*pi))*aimag(sum(Berry_cuv(:,j+1:j+k),dim=2))
        Ch_numqminus=(1.0_dblprec/(2*pi))*aimag(sum(Berry_cuv(:,j+k+1:j+k+l),dim=2))
      else
        Ch_numq=(1.0_dblprec/(2*pi*(max(Nx,Ny,Nz)-1)))*aimag(sum(Berry_cuv(:,1:j),dim=2))
        Ch_numqplus=(1.0_dblprec/(2*pi*(max(Nx,Ny,Nz)-1)))*aimag(sum(Berry_cuv(:,j+1:j+k),dim=2))
        Ch_numqminus=(1.0_dblprec/(2*pi*(max(Nx,Ny,Nz)-1)))*aimag(sum(Berry_cuv(:,j+k+1:j+k+l),dim=2))
      end if

      ! Calculate the thermal magnon conductivity Kxy in units W/K
      ! Definiton of the Bose-Einstein distribution
      rho=1/(exp(nc_eval_complex/(k_bolt_ev*Temp))-1)
      ! Definition of the c2 function
       do i=1,2*NA
         do m=1,3*dimen
           if (rho(i,m) <= 0.0_dblprec .or. rho(i,m) /= rho(i,m)) then
             c2_func(i,m)=0.0_dblprec
           else
             c2_func(i,m)=(1.0_dblprec+rho(i,m))*(log((1.0_dblprec+rho(i,m))/rho(i,m)))**2 - (log(rho(i,m)))**2 - 2.0_dblprec*dli2(-rho(i,m))
           end if
         end do
       end do
       ! Sum in k
      nmx = min(nqred, dimen)
      ! Silence any NaNs by zeroing affected entries (avoid noisy output)
      do i = 1, NA
        do m = 1, 3*nmx
          if (c2_func(i,m) /= c2_func(i,m)) c2_func(i,m) = 0.0_dblprec
          if (rho(i,m) /= rho(i,m)) rho(i,m) = 0.0_dblprec
          if (Berry_cuv(i,m) /= Berry_cuv(i,m)) Berry_cuv(i,m) = (0.0_dblprec,0.0_dblprec)
        end do
      end do

      therm_conduc_band = sum(c2_func(1:NA,1:3*nmx) * aimag(Berry_cuv(1:NA,1:3*nmx)**2), dim=2)
       ! Sum in band index
       therm_conduc=-(k_bolt**2)*Temp/((2*pi)**2*hbar)*sum(therm_conduc_band(1:NA),dim=1)
       !Avoid precision error at very low temperatures
       !Check for NaN: NaN != NaN
       !if (therm_conduc /= therm_conduc) then
       if (therm_conduc /= therm_conduc) then
        therm_conduc=0.0_dblprec
       end if
 
      !print '(1x,a,(2x,i2))', 'Band#', (i, i=1,NA)
      !write(*,'(1x,a,20(2x,i3))') 'Band number ->', (i, i=1,NA)
      !write(*,'(1x,a,20(2x,i3))') 'Ch_Number   ->',(nint(Ch_numq(i)), i=1,NA)
      !write(*,'(1x,a,20(2x,i3))') 'Ch_Number+Q ->',(nint(Ch_numqplus(i)), i=1,NA)
      !write(*,'(1x,a,20(2x,i3))') 'Ch_Number-Q ->',(nint(Ch_numqminus(i)), i=1,NA)
      !write(*,*) 'Thermal Conductivity in W/K ->', therm_conduc, rho
      !print in file Chern_number and Magnon Thermal Conductivity
      chern_file = 'chern.'//trim(simid)//'.out'
      open(ofileno,file=chern_file)
            write(ofileno,1006) "Band number ->",(i, i=1,NA)
            write(ofileno,1006) 'Ch_Number   ->',(nint(Ch_numq(i)), i=1,NA)
            write(ofileno,1006) 'Ch_Number+Q ->',(nint(Ch_numqplus(i)), i=1,NA)
            write(ofileno,1006) 'Ch_Number-Q ->',(nint(Ch_numqminus(i)), i=1,NA)
            write(ofileno,1008) 'Magnon Thermal Conductivity in W/K ->', therm_conduc
      close(ofileno)
      !print the Berry phase
      bphase_file = 'bphase.'//trim(simid)//'.out'
      open(ofileno,file=bphase_file)
            write(ofileno,1002) "Band #          qx          qy          qz",(i, i=1,NA)

        do m=1,j
            write(ofileno,1004)   ((q_vchern(i,m)), i=1,3), (( aimag(Berry_cuv(i,m))),i=1,NA)
            if ( mod(m,Nx-1) .eq. 0) then
            write(ofileno,1004)
            end if
        enddo
      close(ofileno)
      !
      bphase_file = 'bphase+q.'//trim(simid)//'.out'
      open(ofileno,file=bphase_file)
            write(ofileno,1002) "Band #          qx          qy          qz",(i, i=1,NA)
        do m=j+1,j+k
            write(ofileno,1004)   ((q_vchern(i,m-j)), i=1,3), (( aimag(Berry_cuv(i,m))),i=1,NA)
            if (mod(m,Nx-1) .eq. 0) then
            write(ofileno,1004)
            end if
        enddo
      close(ofileno)
      !
      bphase_file = 'bphase-q.'//trim(simid)//'.out'
      open(ofileno,file=bphase_file)
            write(ofileno,1002) "Band #          qx          qy          qz",(i, i=1,NA)
        do m=j+k+1,j+k+l
            write(ofileno,1004)   ((q_vchern(i,m-j-k)), i=1,3), (( aimag(Berry_cuv(i,m))),i=1,NA)
            if (mod(m,Nx-1) .eq. 0) then
            write(ofileno,1004)
            end if
        enddo
      close(ofileno)
      !
      1002 format (a,2000i12)
      1004 format (6x,2000f12.6,2000f12.6)
      1006 format (a,2000i12)
      1008 format (a,3x,ES23.15E3)

      !Deallocate variables
      deallocate(indx,stat=i_stat)
      call memocc(i_stat,product(shape(indx))*kind(indx),'indx','calculate_chern_number')
      deallocate(indy,stat=i_stat)
      call memocc(i_stat,product(shape(indy))*kind(indy),'indy','calculate_chern_number')
      deallocate(indz,stat=i_stat)
      call memocc(i_stat,product(shape(indz))*kind(indz),'indz','calculate_chern_number')
      deallocate(u1,stat=i_stat)
      call memocc(i_stat,product(shape(u1))*kind(u1),'u1','calculate_chern_number')
      deallocate(u1norm,stat=i_stat)
      call memocc(i_stat,product(shape(u1norm))*kind(u1norm),'u1norm','calculate_chern_number')
      deallocate(u2,stat=i_stat)
      call memocc(i_stat,product(shape(u2))*kind(u2),'u2','calculate_chern_number')
      deallocate(u2norm,stat=i_stat)
      call memocc(i_stat,product(shape(u2norm))*kind(u2norm),'u2norm','calculate_chern_number')
      deallocate(u3,stat=i_stat)
      call memocc(i_stat,product(shape(u3))*kind(u3),'u3','calculate_chern_number')
      deallocate(u3norm,stat=i_stat)
      call memocc(i_stat,product(shape(u3norm))*kind(u3norm),'u3norm','calculate_chern_number')
      deallocate(u1inv,stat=i_stat)
      call memocc(i_stat,product(shape(u1inv))*kind(u1inv),'u1inv','calculate_chern_number')
      deallocate(u2inv,stat=i_stat)
      call memocc(i_stat,product(shape(u2inv))*kind(u2inv),'u2inv','calculate_chern_number')
      deallocate(u3inv,stat=i_stat)
      call memocc(i_stat,product(shape(u3inv))*kind(u3inv),'u3inv','calculate_chern_number')
      deallocate(u1invnorm,stat=i_stat)
      call memocc(i_stat,product(shape(u1invnorm))*kind(u1invnorm),'u1invnorm','calculate_chern_number')
      deallocate(u2invnorm,stat=i_stat)
      call memocc(i_stat,product(shape(u2invnorm))*kind(u2invnorm),'u2invnorm','calculate_chern_number')
      deallocate(u3invnorm,stat=i_stat)
      call memocc(i_stat,product(shape(u3invnorm))*kind(u3invnorm),'u3invnorm','calculate_chern_number')
      deallocate(Berry_cuv,stat=i_stat)
      call memocc(i_stat,product(shape(Berry_cuv))*kind(Berry_cuv),'Berry_cuv','calculate_chern_number')
      deallocate(Ch_numq,stat=i_stat)
      call memocc(i_stat,product(shape(Ch_numq))*kind(Ch_numq),'Ch_numq','calculate_chern_number')
      deallocate(Ch_numqplus,stat=i_stat)
      call memocc(i_stat,product(shape(Ch_numqplus))*kind(Ch_numqplus),'Ch_numqplus','calculate_chern_number')
      deallocate(Ch_numqminus,stat=i_stat)
      call memocc(i_stat,product(shape(Ch_numqminus))*kind(Ch_numqminus),'Ch_numqminus','calculate_chern_number')
      deallocate(q_vchern,stat=i_stat)
      call memocc(i_stat,product(shape(q_vchern))*kind(q_vchern),'q_vchern','calculate_chern_number')
      deallocate(rho,stat=i_stat)
      call memocc(i_stat,product(shape(rho))*kind(rho),'rho','calculate_chern_number')
      deallocate(c2_func,stat=i_stat)
      call memocc(i_stat,product(shape(c2_func))*kind(c2_func),'c2_func','calculate_chern_number')
      deallocate(therm_conduc_band,stat=i_stat)
      call memocc(i_stat,product(shape(therm_conduc_band))*kind(therm_conduc_band),'therm_conduc_band','calculate_chern_number')
      deallocate(nc_eval_complex,stat=i_stat)
      call memocc(i_stat,product(shape(nc_eval_complex))*kind(nc_eval_complex),'nc_eval_complex','calculate_chern_number')
      deallocate(nc_evec_complex,stat=i_stat)
      call memocc(i_stat,product(shape(nc_evec_complex))*kind(nc_evec_complex),'nc_evec_complex','calculate_chern_number')
      if (allocated(oam_evec)) then
         deallocate(oam_evec,stat=i_stat)
         call memocc(i_stat,product(shape(oam_evec))*kind(oam_evec),'oam_evec','calculate_chern_number')
      end if
      if (allocated(oam_hbar)) then
         deallocate(oam_hbar,stat=i_stat)
         call memocc(i_stat,product(shape(oam_hbar))*kind(oam_hbar),'oam_hbar','calculate_chern_number')
      end if
      if (allocated(oam_band)) then
         deallocate(oam_band,stat=i_stat)
         call memocc(i_stat,product(shape(oam_band))*kind(oam_band),'oam_band','calculate_chern_number')
      end if
      !
      print '(1x,a)', 'Chern calculation done.'
   !
   return
   !
   end subroutine calculate_chern_number

   !> Continue positive-energy modes over the reciprocal grid.  The first
   !> point is an arbitrary anchor; every later point is matched to the
   !> previous point in the row, or to the previous row at a row boundary.
   !> The matching and phase rotation use the bosonic metric overlap.
   subroutine prepare_oam_eigenvectors(NA,Nx,Ny,Nz,dimen,eigenvectors,oam_evec,oam_band)
      implicit none
      integer, intent(in) :: NA,Nx,Ny,Nz,dimen
      complex(dblprec), intent(in) :: eigenvectors(2*NA,2*NA,*)
      complex(dblprec), intent(out) :: oam_evec(2*NA,NA,dimen)
      integer, intent(out) :: oam_band(NA,dimen)

      complex(dblprec) :: reference(2*NA,NA)
      integer :: icount,jcount,kcount,iq,band
      integer :: mode_order(NA),raw_order(NA)

      oam_evec=eigenvectors(:,1:NA,1:dimen)
      do iq=1,dimen
         do band=1,NA
            oam_band(band,iq)=band
         end do
      end do
      do kcount=1,Nz
         do jcount=1,Ny
            do icount=1,Nx
               iq=(kcount-1)*Nx*Ny+(jcount-1)*Nx+icount
               if (iq == 1) cycle
               if (icount > 1) then
                  reference=oam_evec(:,:,iq-1)
               else
                  reference=oam_evec(:,:,iq-Nx)
               end if
               call continue_boson_modes(reference,oam_evec(:,:,iq),NA,mode_order)
               raw_order=oam_band(:,iq)
               do band=1,NA
                  oam_band(band,iq)=raw_order(mode_order(band))
               end do
            end do
         end do
      end do
   end subroutine prepare_oam_eigenvectors

   !> Match a set of modes to a reference set and fix each phase so that the
   !> overlap with its reference mode is real and positive.  This is a
   !> deterministic spanning-tree gauge for the pointwise OAM diagnostic.
   subroutine continue_boson_modes(reference,current,NA,mode_order)
      implicit none
      integer, intent(in) :: NA
      complex(dblprec), intent(in) :: reference(2*NA,NA)
      complex(dblprec), intent(inout) :: current(2*NA,NA)
      integer, intent(out) :: mode_order(NA)

      complex(dblprec) :: raw_modes(2*NA,NA), reordered(2*NA,NA)
      complex(dblprec) :: overlap, phase_factor
      real(dblprec) :: best_overlap
      logical :: used(NA)
      integer :: i,j,best

      raw_modes=current
      reordered=(0.0_dblprec,0.0_dblprec)
      used=.false.
      do i=1,NA
         best=0
         best_overlap=-1.0_dblprec
         do j=1,NA
            if (.not.used(j)) then
               overlap=boson_overlap(reference(:,i),raw_modes(:,j),NA)
               if (abs(overlap) > best_overlap) then
                  best_overlap=abs(overlap)
                  best=j
               end if
            end if
         end do
         if (best == 0) error stop 'chern: failed boson band continuation'
         used(best)=.true.
         mode_order(i)=best
         overlap=boson_overlap(reference(:,i),raw_modes(:,best),NA)
         if (abs(overlap) > 1000.0_dblprec*epsilon(1.0_dblprec)) then
            phase_factor=conjg(overlap)/abs(overlap)
         else
            phase_factor=(1.0_dblprec,0.0_dblprec)
         end if
         reordered(:,i)=raw_modes(:,best)*phase_factor
      end do
      current=reordered
   end subroutine continue_boson_modes

   !> Evaluate Fishman's pointwise magnon OAM in units of hbar.
   !> The reciprocal grid is expressed in the same coordinates as q_vchern;
   !> the common 2*pi factor used by setup_ektij cancels between k and d/dk.
   subroutine calculate_fishman_oam(NA,Nx,Ny,Nz,dimen,q_vchern,oam_evec,oam_hbar)
      implicit none
      integer, intent(in) :: NA,Nx,Ny,Nz,dimen
      real(dblprec), intent(in) :: q_vchern(3,dimen)
      complex(dblprec), intent(in) :: oam_evec(2*NA,NA,dimen)
      real(dblprec), intent(out) :: oam_hbar(NA,dimen)

      real(dblprec) :: dq1(3),dq2(3),detq,qx,qy
      complex(dblprec) :: d1(2*NA),d2(2*NA),dtdx(2*NA),dtdy(2*NA)
      complex(dblprec) :: angular_derivative(2*NA)
      integer :: icount,jcount,iq,iqp,iqm,band

      if (Nx < 2 .or. Ny < 2 .or. Nz /= 1) then
         error stop 'chern: Fishman OAM requires a two-dimensional k mesh'
      end if

      ! The grid is q=x*b1/Nx+y*b2/Ny, so these are the two directional
      ! increments even for oblique reciprocal lattices.
      dq1=q_vchern(:,2)-q_vchern(:,1)
      dq2=q_vchern(:,Nx+1)-q_vchern(:,1)
      detq=dq1(1)*dq2(2)-dq1(2)*dq2(1)
      if (abs(detq) <= 1000.0_dblprec*epsilon(1.0_dblprec)) then
         error stop 'chern: singular reciprocal directions for Fishman OAM'
      end if

      do jcount=1,Ny
         do icount=1,Nx
            iq=(jcount-1)*Nx+icount
            qx=q_vchern(1,iq)
            qy=q_vchern(2,iq)
            do band=1,NA
               ! Central differences in the interior.  The two boundary
               ! rows use first-order differences because Fishman's OAM is
               ! not periodic in k; no Berry-flux quantity is substituted.
               if (icount == 1) then
                  iqp=iq+1
                  d1=oam_evec(:,band,iqp)-oam_evec(:,band,iq)
               else if (icount == Nx) then
                  iqm=iq-1
                  d1=oam_evec(:,band,iq)-oam_evec(:,band,iqm)
               else
                  iqp=iq+1
                  iqm=iq-1
                  d1=0.5_dblprec*(oam_evec(:,band,iqp)-oam_evec(:,band,iqm))
               end if

               if (jcount == 1) then
                  iqp=iq+Nx
                  d2=oam_evec(:,band,iqp)-oam_evec(:,band,iq)
               else if (jcount == Ny) then
                  iqm=iq-Nx
                  d2=oam_evec(:,band,iq)-oam_evec(:,band,iqm)
               else
                  iqp=iq+Nx
                  iqm=iq-Nx
                  d2=0.5_dblprec*(oam_evec(:,band,iqp)-oam_evec(:,band,iqm))
               end if

               ! D1 = dq1_x*dT/dx + dq1_y*dT/dy and likewise for D2.
               dtdx=(dq2(2)*d1-dq1(2)*d2)/detq
               dtdy=(-dq2(1)*d1+dq1(1)*d2)/detq
               angular_derivative=qx*dtdy-qy*dtdx

               ! UppASD uses exp(-i*k.R) in setup_ektij.  Thus l_z is
               ! -i(k_x d_y-k_y d_x), and with the
               ! eta-normalized bosonic vector this is the Fishman value
               ! Lz/hbar = 1/2 Im[T^dagger eta D T].
               oam_hbar(band,iq)=0.5_dblprec*aimag(&
                  boson_overlap(oam_evec(:,band,iq),angular_derivative,NA))
            end do
         end do
      end do
   end subroutine calculate_fishman_oam

   subroutine setup_grid(Nx,Ny,Nz,C1,C2,C3,dimen,q_vchern)
      ! Set up grid in reciprocal space (1st BZ)
      !
      implicit none
      !
      integer, intent(in) :: Nx  !< Number of points of the grid in x direction
      integer, intent(in) :: Ny  !< Number of points of the grid in y direction
      integer, intent(in) :: Nz  !< Number of points of the grid in z direction
      integer, intent(in) :: dimen !< Size of the grid
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      !
      real(dblprec), dimension(3,dimen), intent(out) :: q_vchern !< q-vectors
      !
      integer :: iq,xq,yq,zq
      integer :: i_stat, i_all
      real(dblprec), dimension(3) :: b1,r1
      real(dblprec), dimension(3) :: b2,r2
      real(dblprec), dimension(3) :: b3,r3
      real(dblprec) :: c1r1, c2r2, c3r3

      ! Calculate reciprocal lattice vectors
      ! r1 = C2xC3
      r1(1)=C2(2)*C3(3)-C2(3)*C3(2)
      r1(2)=C2(3)*C3(1)-C2(1)*C3(3)
      r1(3)=C2(1)*C3(2)-C2(2)*C3(1)
      ! r2 = C3xC1
      r2(1)=C3(2)*C1(3)-C3(3)*C1(2)
      r2(2)=C3(3)*C1(1)-C3(1)*C1(3)
      r2(3)=C3(1)*C1(2)-C3(2)*C1(1)
      ! r3 = C1xC2
      r3(1)=C1(2)*C2(3)-C1(3)*C2(2)
      r3(2)=C1(3)*C2(1)-C1(1)*C2(3)
      r3(3)=C1(1)*C2(2)-C1(2)*C2(1)
      ! cell volume C1*(C2xC3)
      c1r1=C1(1)*r1(1)+C1(2)*r1(2)+C1(3)*r1(3)
      c2r2=C2(1)*r2(1)+C2(2)*r2(2)+C2(3)*r2(3)
      c3r3=C3(1)*r3(1)+C3(2)*r3(2)+C3(3)*r3(3)
      ! b1=(2pi)*r1/(C1*r1)
      b1(1)=r1(1)/c1r1
      b1(2)=r1(2)/c1r1
      b1(3)=r1(3)/c1r1
      ! b2=(2pi)*r2/(C1*r1)
      b2(1)=r2(1)/c2r2
      b2(2)=r2(2)/c2r2
      b2(3)=r2(3)/c2r2
      ! b3=(2pi)*r3/(C1*r1)
      b3(1)=r3(1)/c3r3
      b3(2)=r3(2)/c3r3
      b3(3)=r3(3)/c3r3
      !Initialize variables
      iq=0
      !q-points expressed in the bases of the reciprocal lattice vectors
      do zq=-(Nz-1)/2,(Nz)/2
        do yq=-(Ny-1)/2,(Ny)/2
          do xq=-(Nx-1)/2,(Nx)/2
            iq=iq+1
            q_vchern(:,iq)=xq/(1.0_dblprec*Nx)*b1+yq/(1.0_dblprec*Ny)*b2+zq/(1.0_dblprec*Nz)*b3
          end do
        end do
      end do
      !
      return
   !
   end subroutine setup_grid

   subroutine read_parameters_chern_number(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
      integer :: rd_len,i_err,i_errb
      logical :: comment

      call setup_chern_number()

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

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! START OF VARIABLES FOR CALCULATING CHERN NUMBERS
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            case('do_chern') ! Calculate the Chern number of the bands
              read(ifile,*,iostat=i_err) do_chern
              if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('do_magnon_oam') ! Fishman reciprocal-space magnon OAM (requires do_chern Y)
              read(ifile,*,iostat=i_err) do_magnon_oam
              if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('kgrid') ! Read the size of the grid
              read(ifile,*,iostat=i_err) Nx, Ny, Nz
              if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err

            case('Chern_qvect') ! Ordering wave vector
            read(ifile,*,iostat=i_err) Chern_qvect
            if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword),' data',i_err
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! END OF LEGACY VARIABLES
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

   if (do_magnon_oam=='Y' .and. do_chern/='Y') then
      error stop 'do_magnon_oam requires do_chern Y'
   end if

   return
   end subroutine read_parameters_chern_number

! Dilogaritmic real function

real(dblprec) function dli2(x)
  implicit none
  double precision :: x, y, r, s, y2, y4, p, q, l
  double precision, parameter :: PI = 3.14159265358979324D0
  double precision, parameter :: cp(6) = (/ &
      0.9999999999999999502D+0,             &
     -2.6883926818565423430D+0,             &
      2.6477222699473109692D+0,             &
     -1.1538559607887416355D+0,             &
      2.0886077795020607837D-1,             &
     -1.0859777134152463084D-2              /)
  double precision, parameter :: cq(7) = (/ &
      1.0000000000000000000D+0,             &
     -2.9383926818565635485D+0,             &
      3.2712093293018635389D+0,             &
     -1.7076702173954289421D+0,             &
      4.1596017228400603836D-1,             &
     -3.9801343754084482956D-2,             &
      8.2743668974466659035D-4              /)

  ! transform to [0, 1/2]
   if (x .lt. -1) then
      l = log(1 - x)
      y = 1/(1 - x)
      r = -PI**2/6 + l*(0.5D0*l - log(-x))
      s = 1
   elseif (x .eq. -1) then
      dli2 = -PI**2/12
      return
   elseif (x .lt. 0) then
      y = x/(x - 1)
      r = -0.5D0*log(1 - x)**2
      s = -1
   elseif (x .eq. 0) then
      dli2 = 0
      return
   elseif (x .lt. 0.5D0) then
      y = x
      r = 0
      s = 1
   elseif (x .lt. 1) then
      y = 1 - x
      r = PI**2/6 - log(x)*log(y)
      s = -1
   elseif (x .eq. 1) then
      dli2 = PI**2/6
      return
   elseif (x .lt. 2) then
      l = log(x)
      y = 1 - 1/x
      r = PI**2/6 - l*(log(y) + 0.5D0*l)
      s = 1
   else
      y = 1/x
      r = PI**2/3 - 0.5D0*log(x)**2
      s = -1
   endif

  y2 = y*y
  y4 = y2*y2
  p = cp(1) + y * cp(2) + y2 * (cp(3) + y * cp(4)) +      &
      y4 * (cp(5) + y * cp(6))
  q = cq(1) + y * cq(2) + y2 * (cq(3) + y * cq(4)) +      &
      y4 * (cq(5) + y * cq(6) + y2 * cq(7))

  dli2 = r + s*y*p/q

end function dli2

! Dilogaritmic complex function
!double complex function cdli2(z)
!  implicit none
!  double complex :: z, rest, u, u2, u4, sum, fast_cdlog
!  double precision :: rz, iz, nz, sgn, dli2
!  double precision, parameter :: PI = 3.14159265358979324D0
!  double precision, parameter :: bf(10) = (/ &
!    - 1.0D0/4.0D0,                           &
!    + 1.0D0/36.0D0,                          &
!    - 1.0D0/3600.0D0,                        &
!    + 1.0D0/211680.0D0,                      &
!    - 1.0D0/10886400.0D0,                    &
!    + 1.0D0/526901760.0D0,                   &
!    - 4.0647616451442255D-11,                &
!    + 8.9216910204564526D-13,                &
!    - 1.9939295860721076D-14,                &
!    + 4.5189800296199182D-16                 /)
!
!  rz = real(z)
!  iz = aimag(z)
!
!  ! special cases
!  if (iz .eq. 0) then
!     if (rz .le. 1) cdli2 = dcmplx(dli2(rz), 0)
!     if (rz .gt. 1) cdli2 = dcmplx(dli2(rz), -PI*log(rz))
!     return
!  endif
!
!  nz = rz**2 + iz**2
!
!  if (nz .lt. EPSILON(1D0)) then
!     cdli2 = z*(1 + 0.25D0*z)
!     return
!  endif
!
!  ! transformation to |z| < 1, Re(z) <= 0.5
!  if (rz .le. 0.5D0) then
!     if (nz .gt. 1) then
!        u = -fast_cdlog(1 - 1/z)
!        rest = -0.5D0*fast_cdlog(-z)**2 - PI**2/6
!        sgn = -1
!     else ! nz <= 1
!        u = -fast_cdlog(1 - z)
!        rest = 0
!        sgn = 1
!     endif
!  else ! rz > 0.5D0
!     if (nz .le. 2*rz) then
!        u = -fast_cdlog(z)
!        rest = u*fast_cdlog(1 - z) + PI**2/6
!        sgn = -1
!     else ! nz > 2*rz
!        u = -fast_cdlog(1 - 1/z)
!        rest = -0.5D0*fast_cdlog(-z)**2 - PI**2/6
!        sgn = -1
!     endif
!  endif
!
!  u2 = u**2
!  u4 = u2**2
!  sum =                                                    &
!     u +                                                   &
!     u2 * (bf(1) +                                         &
!     u  * (bf(2) +                                         &
!     u2 * (                                                &
!         bf(3) +                                           &
!         u2*bf(4) +                                        &
!         u4*(bf(5) + u2*bf(6)) +                           &
!         u4*u4*(bf(7) + u2*bf(8) + u4*(bf(9) + u2*bf(10))) &
!     )))
!
!  cdli2 = sgn*sum + rest
!
!end function cdli2


end module Chern_number
