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
   use, intrinsic :: ieee_arithmetic,  only : ieee_is_nan
   use InputData,  only : Temp
   use Parameters
   use Constants
   use Profiling
   use Hamiltoniandata,    only : ham
   use InputData,   only : ham_inp
   use Diamag ,     only : clone_q,diagonalize_quad_hamiltonian,find_uv,setup_ektij,&
                            setup_jtens2_q,setup_jtens_q,sJs, setup_tensor_hamiltonian,&
                            nc_eval_qchern,nc_evec_qchern
   !
   implicit none
   !
   character(len=1)                           :: do_chern    !< Calculate the Chern number of the bands (Y/N)
   integer                                    :: Nx          !< Number of points of the grid in x direction
   integer                                    :: Ny          !< Number of points of the grid in y direction
   integer                                    :: Nz          !< Number of points of the grid in z direction
   real(dblprec), dimension(3)                :: Chern_qvect !< Spin spiral ordering vector
   !
   private
   ! public subroutines
   public :: do_chern,Nx,Ny,Nz, Chern_qvect
   public :: read_parameters_chern_number,calculate_chern_number
   !
contains

   subroutine setup_chern_number()

      implicit none

      do_chern    = 'N'
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
      !
      real(dblprec)                                   :: therm_conduc              !< Thermal conductivity
      integer                                         :: iq,i,j,k,l,m,i_stat,nqred,icount,jcount,kcount,counter,kx,ky,kz
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
                u1(i,j)=dot_product(nc_evec_qchern(1:NA,i,iq),nc_evec_qchern(1:NA,i,iq+1))
              if (abs(u1(i,j))== 0.0_dblprec) then
                u1norm(i,j)=0.0_dblprec
              else
                u1norm(i,j)=u1(i,j)/abs(u1(i,j))
              end if

                u2(i,j)=dot_product(nc_evec_qchern(1:NA,i,iq+1),nc_evec_qchern(1:NA,i,iq+1+kx))
              if (abs(u2(i,j))== 0.0_dblprec) then
                u2norm(i,j)=0.0_dblprec
              else
                u2norm(i,j)=u2(i,j)/abs(u2(i,j))
              end if

                u3(i,j)=dot_product(nc_evec_qchern(1:NA,i,iq+1+kx),nc_evec_qchern(1:NA,i,iq+1+kx+kz))
              if (abs(u3(i,j))== 0.0_dblprec) then
                u3norm(i,j)=0.0_dblprec
              else
                u3norm(i,j)=u3(i,j)/abs(u3(i,j))
              end if

                u1inv(i,j)=1.0_dblprec/dot_product(nc_evec_qchern(1:NA,i,iq+kx+kz),nc_evec_qchern(1:NA,i,iq+1+kx+kz))
              if (abs(u1inv(i,j))== 0.0_dblprec) then
                u1invnorm(i,j)=0.0_dblprec
              else
                u1invnorm(i,j)=u1inv(i,j)/abs(u1inv(i,j))
              end if

                u2inv(i,j)=1.0_dblprec/dot_product(nc_evec_qchern(1:NA,i,iq+kz),nc_evec_qchern(1:NA,i,iq+kx+kz))
              if (abs(u2inv(i,j))== 0.0_dblprec) then
                u2invnorm(i,j)=0.0_dblprec
              else
                u2invnorm(i,j)=u2inv(i,j)/abs(u2inv(i,j))
              end if

                u3inv(i,j)=1.0_dblprec/dot_product(nc_evec_qchern(1:NA,i,iq),nc_evec_qchern(1:NA,i,iq+kz))
              if (abs(u3inv(i,j))== 0.0_dblprec) then
                u3invnorm(i,j)=0.0_dblprec
              else
                u3invnorm(i,j)=u3inv(i,j)/abs(u3inv(i,j))
              end if
            end do
          end if

          if ( indx(iq).eq.2 .and. indy(iq).eq.2 .and. indz(iq).eq.2 ) then
            k=k+1
            do i=1,NA !band index
              u1(i,j+k)=dot_product(nc_evec_qchern(1:NA,i,iq),nc_evec_qchern(1:NA,i,iq+1))
              if (abs(u1(i,j+k))== 0.0_dblprec) then
              u1norm(i,j+k)=0.0_dblprec
              else
              u1norm(i,j+k)=u1(i,j+k)/abs(u1(i,j+k))
              end if

              u2(i,j+k)=dot_product(nc_evec_qchern(1:NA,i,iq+1),nc_evec_qchern(1:NA,i,iq+1+kx))
              if (abs(u2(i,j+k))== 0.0_dblprec) then
              u2norm(i,j+k)=0.0_dblprec
              else
              u2norm(i,j+k)=u2(i,j+k)/abs(u2(i,j+k))
              end if

              u3(i,j+k)=dot_product(nc_evec_qchern(1:NA,i,iq+1+kx),nc_evec_qchern(1:NA,i,iq+1+kx+kz))
              if (abs(u3(i,j+k))== 0.0_dblprec) then
                u3norm(i,j+k)=0.0_dblprec
              else
                u3norm(i,j+k)=u3(i,j+k)/abs(u3(i,j+k))
              end if

              u1inv(i,j+k)=1.0_dblprec/dot_product(nc_evec_qchern(1:NA,i,iq+kx+kz),nc_evec_qchern(1:NA,i,iq+1+kx+kz))
              if (abs(u1inv(i,j+k))== 0.0_dblprec) then
              u1invnorm(i,j+k)=0.0_dblprec
              else
              u1invnorm(i,j+k)=u1inv(i,j+k)/abs(u1inv(i,j+k))
              end if

              u2inv(i,j+k)=1.0_dblprec/dot_product(nc_evec_qchern(1:NA,i,iq+kz),nc_evec_qchern(1:NA,i,iq+kx+kz))
              if (abs(u2inv(i,j+k))== 0.0_dblprec) then
              u2invnorm(i,j+k)=0.0_dblprec
              else
              u2invnorm(i,j+k)=u2inv(i,j+k)/abs(u2inv(i,j+k))
              end if

              u3inv(i,j+k)=1.0_dblprec/dot_product(nc_evec_qchern(1:NA,i,iq),nc_evec_qchern(1:NA,i,iq+kz))
              if (abs(u3inv(i,j+k))== 0.0_dblprec) then
                u3invnorm(i,j+k)=0.0_dblprec
              else
                u3invnorm(i,j+k)=u3inv(i,j+k)/abs(u3inv(i,j+k))
              end if
            end do
          end if

          if ( indx(iq).eq.3 .and. indy(iq).eq.3 .and. indz(iq).eq.3 ) then
            l=l+1
            do i=1,NA !band index
              u1(i,j+k+l)=dot_product(nc_evec_qchern(1:NA,i,iq),nc_evec_qchern(1:NA,i,iq+1))
              if (abs(u1(i,j+k+l))== 0.0_dblprec) then
              u1norm(i,j+k+l)=0.0_dblprec
              else
              u1norm(i,j+k+l)=u1(i,j+k+l)/abs(u1(i,j+k+l))
              end if

              u2(i,j+k+l)=dot_product(nc_evec_qchern(1:NA,i,iq+1),nc_evec_qchern(1:NA,i,iq+1+kx))
              if (abs(u2(i,j+k+l))== 0.0_dblprec) then
              u2norm(i,j+k+l)=0.0_dblprec
              else
              u2norm(i,j+k+l)=u2(i,j+k+l)/abs(u2(i,j+k+l))
              end if

              u3(i,j+k+l)=dot_product(nc_evec_qchern(1:NA,i,iq+1+kx),nc_evec_qchern(1:NA,i,iq+1+kx+kz))
              if (abs(u3(i,j+k+l))== 0.0_dblprec) then
                u3norm(i,j+k+l)=0.0_dblprec
              else
                u3norm(i,j+k+l)=u3(i,j+k+l)/abs(u3(i,j+k+l))
              end if

              u1inv(i,j+k+l)=1.0_dblprec/dot_product(nc_evec_qchern(1:NA,i,iq+kx+kz),nc_evec_qchern(1:NA,i,iq+1+kx+kz))
              if (abs(u1inv(i,j+k+l))== 0.0_dblprec) then
              u1invnorm(i,j+k+l)=0.0_dblprec
              else
              u1invnorm(i,j+k+l)=u1inv(i,j+k+l)/abs(u1inv(i,j+k+l))
              end if

              u2inv(i,j+k+l)=1.0_dblprec/dot_product(nc_evec_qchern(1:NA,i,iq+kz),nc_evec_qchern(1:NA,i,iq+kx+kz))
              if (abs(u2inv(i,j+k+l))== 0.0_dblprec) then
              u2invnorm(i,j+k+l)=0.0_dblprec
              else
              u2invnorm(i,j+k+l)=u2inv(i,j+k+l)/abs(u2inv(i,j+k+l))
              end if

              u3inv(i,j+k+l)=1.0_dblprec/dot_product(nc_evec_qchern(1:NA,i,iq),nc_evec_qchern(1:NA,i,iq+kz))
              if (abs(u3inv(i,j+k+l))== 0.0_dblprec) then
                u3invnorm(i,j+k+l)=0.0_dblprec
              else
                u3invnorm(i,j+k+l)=u3inv(i,j+k+l)/abs(u3inv(i,j+k+l))
              end if
            end do
          end if
      enddo

      Berry_cuv=log(u1norm*u2norm*u3norm*u1invnorm*u2invnorm*u3invnorm)

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
       rho=1/(exp(nc_eval_qchern/(k_bolt_ev*Temp))-1)
      ! Definition of the c2 function
       do i=1,2*NA
         do j=1,3*dimen
           c2_func(i,j)=(1+rho(i,j))*(log((1+rho(i,j))/rho(i,j)))**2-(log(rho(i,j)))**2-2*dli2(-rho(i,j))
         end do
       end do
       ! Sum in k
       therm_conduc_band=sum(c2_func(1:NA,1:j)*aimag(Berry_cuv(1:NA,1:j)**2),dim=2)
       ! Sum in band index
       therm_conduc=-(k_bolt**2)*Temp/((2*pi)**2*hbar)*sum(therm_conduc_band(1:NA),dim=1)
       !Avoid precision error at very low temperatures
       if (ieee_is_nan(therm_conduc)) then
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
      deallocate(nc_eval_qchern,stat=i_stat)
      call memocc(i_stat,product(shape(nc_eval_qchern))*kind(nc_eval_qchern),'nc_eval_qchern','calculate_chern_number')
      deallocate(nc_evec_qchern,stat=i_stat)
      call memocc(i_stat,product(shape(nc_evec_qchern))*kind(nc_evec_qchern),'nc_evec_qchern','calculate_chern_number')
      !
      print '(1x,a)', 'Chern calculation done.'
   !
   return
   !
   end subroutine calculate_chern_number

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
