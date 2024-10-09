!> Library for random number generators
module RandomDrivers
    use Parameters
    use Profiling
    use RandomNumbers
    use omp_lib
 
 
    implicit none
 
    real(dblprec), dimension(:,:,:), allocatable :: ranv !< Work array for RNG
    real(dblprec), dimension(:,:,:), allocatable :: lattranv !< Work array for RNG
 
    public :: lattranv, lattrannum
    public :: ranv, rannum
    public :: allocate_randomwork
 
 
 contains
 

   !> Allocate work arrays for RNG
   subroutine allocate_randomwork(Natom,Mensemble,flag,do_ld)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)
      character(len=1) :: do_ld !< Do lattice dynamics ('Y'/'N')


      integer :: i_stat,i_all

      if (flag>0) then
         allocate(ranv(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ranv))*kind(ranv),'ranv','allocate_randomwork')
         if(do_ld == 'Y') then
            allocate(lattranv(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(lattranv))*kind(lattranv),'lattranv','allocate_randomwork')
         end if
      else
         i_all=-product(shape(ranv))*kind(ranv)
         deallocate(ranv,stat=i_stat)
         call memocc(i_stat,i_all,'ranv','allocate_systemdata')
         if(do_ld == 'Y') then
            i_all=-product(shape(lattranv))*kind(lattranv)
            deallocate(lattranv,stat=i_stat)
            call memocc(i_stat,i_all,'lattranv','allocate_systemdata')
         end if
      end if
   end subroutine allocate_randomwork


   !> Sets up an array of random numbers
   subroutine rannum(Natom, Mensemble, NA,  llg, lambda1_array, lambda2_array, &
         compensate_drift, bn, field1, field2, mmomi, Temp_array,temprescale)

      use Constants, only : k_bolt, gama, mub, hbar
      use InputData, only : para_rng, delta_t
      use Chroma, only : qc_generator, do_quant_colour, qc_single_generator

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom), intent(in) :: lambda2_array !< Additional damping parameter (not used for llg=1)
      integer, intent(in) :: compensate_drift !< Correct for drift in RNG
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.00_dblprec)
      real(dblprec), dimension(3,Mensemble), intent(in) :: field1 !< Average internal effective field
      real(dblprec), dimension(3,Mensemble), intent(in) :: field2 !< Average external effective field
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
      real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature (array)
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

      real(dblprec) :: D
      real(dblprec) :: delta_tau
      real(dblprec) :: mu, sigma
      real(dblprec), dimension(Mensemble,Natom) :: Dk
      real(dblprec), dimension(Mensemble) :: avf1, avf2
      real(dblprec) :: rx(NA),ry(NA),rz(NA)
      integer :: ity
      integer :: i, j

      do j=1,Mensemble !loop over simulations, avf1,2 different for each sim.
         avf1(j)=sqrt(field1(1,j)**2+field1(2,j)**2+field1(3,j)**2)
         avf2(j)=sqrt(field2(1,j)**2+field2(2,j)**2+field2(3,j)**2)
      end do

      !   LL equations ONE universal damping
      Dk(:,:)=0.00_dblprec
      if(llg==0) then
         do i=1,Natom
            Dk(:,i)=(lambda1_array(i)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         enddo
         !   LLG equations ONE universal damping
      else if (llg==1) then
         do i=1, Natom
            Dk(:,i)=(lambda1_array(i)/(1+lambda1_array(i)**2)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less
         enddo
         !   LL equations TWO damping parameters
      else if(llg==2) then
         do i=1, Natom
            do j=1,Mensemble !loop over simulations, Dk different in each sim.
               Dk(j,i)=(lambda1_array(i)*avf1(j)+lambda2_array(i)*avf2(j))/(avf1(j)+avf2(j)) &
                  *(k_bolt/gama/(mub))*(gama/bn) !last factor for dim. less.
            end do
         enddo
         !   LLG equations TWO damping parameters, fluctuations included in damping1
      else if(llg==3) then
         do i=1, Natom
            do j=1,Mensemble !loop over simulations, Dk different in each sim.
               Dk(j,i)=(lambda1_array(i)*avf1(j)+lambda2_array(i)*avf2(j))/(avf1(j)+avf2(j)) &
                  /(1+lambda1_array(i)**2)*(k_bolt/gama/(mub))*(gama/bn)!last factor for dim. les
            end do
         end do
         !   LL equations TWO damping parameters, but use thermal fluctuations corresponding to one universal damping
      else if(llg==4) then
         do i=1, Natom
            Dk(:,i)=(lambda1_array(i)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         enddo
         !   LLG equations TWO damping parameters, but use thermal fluctuations corresponding to one universl damping
      else if(llg==5) then
         do i=1,Natom
            Dk(:,i)=(lambda1_array(i)/(1+lambda1_array(i)**2)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         enddo
         !   LL equations ONE universal damping write energies
      else if(llg==6) then
         do i=1, Natom
            Dk(:,i)=(lambda1_array(i)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         enddo
         !   LLG equations ONE universal damping write energies
      else if (llg==7) then
         do i=1, Natom
            Dk(:,i)=(lambda1_array(i)/(1+lambda1_array(i)**2)*k_bolt/gama/(mub))*(gama/bn) !last factor for dim. less.
         enddo
      endif

      if (do_quant_colour=='N') then
      ! Current RNG is not parallel!
      if(.not.para_rng) then
#ifdef VSL
         if(use_vsl) then
            !            call rng_gaussian(ranv,3*Natom*Mensemble,1.0_dblprec)
            call rng_gaussianP(ranv,3*Natom*Mensemble,1.0_dblprec)
         else
            call fill_rngarray(ranv,3*Natom*Mensemble)
         end if
#else
         call fill_rngarray(ranv,3*Natom*Mensemble)
#endif
      else
#ifdef VSL
         if(use_vsl) then
            !            call rng_gaussian(ranv,3*Natom*Mensemble,1.0_dblprec)
            call rng_gaussianP(ranv,3*Natom*Mensemble,1.0_dblprec)
         else
            call fill_rngarray_para(ranv,3*Natom*Mensemble)
         end if
#else
         call fill_rngarray_para(ranv,3*Natom*Mensemble)
#endif
      end if
   else if (do_quant_colour=='Y') then
         !print *, 'Rand pre',ranv(1,1,1)
         call qc_generator(ranv,3*Natom*Mensemble,delta_t)
   else if (do_quant_colour=='S') then
         !print *, 'Rand pre',ranv(1,1,1)
         call qc_single_generator(ranv,3*Natom*Mensemble,delta_t)
  end if

      !!!!!++ DK mod in merge, should not be a problem. Thomas(2014/07/17)
      !$omp parallel do default(shared) private(i,j,D,sigma,mu) collapse(2)
      do j=1, Mensemble
         do i=1, Natom
            D=Dk(j,i)*mmomi(i,j)*Temp_array(i)*temprescale
            sigma=sqrt(2.00_dblprec*D)
            mu=00_dblprec
            ranv(1,i,j)=ranv(1,i,j)*sigma
            ranv(2,i,j)=ranv(2,i,j)*sigma
            ranv(3,i,j)=ranv(3,i,j)*sigma
         end do
      end do
      !$omp end parallel do

      ! Possibility to correct the random number distribution so that
      ! the true mean is zero. Not used by default
      if(compensate_drift==1) then
         do j=1, Mensemble
            rx=0.00_dblprec;ry=0.00_dblprec;rz=0.00_dblprec
            do i=1, Natom
               ity=mod(i-1,NA)+1
               rx(ity)=rx(ity)+ranv(1,i,j)
               ry(ity)=ry(ity)+ranv(2,i,j)
               rz(ity)=rz(ity)+ranv(3,i,j)
            end do
            rx=rx/Natom*NA
            ry=ry/Natom*NA
            rz=rz/Natom*NA
            do i=1, Natom
               ity=mod(i-1,NA)+1
               ranv(1,i,j)=ranv(1,i,j)-rx(ity)
               ranv(2,i,j)=ranv(2,i,j)-ry(ity)
               ranv(3,i,j)=ranv(3,i,j)-rz(ity)
            end do
         end do
      end if

   end subroutine rannum

   !> Sets up an array of random numbers for the lattice dynamics
   subroutine lattrannum(Natom, Mensemble, NA, deltat, lattdampvec, compensate_drift, Temp_array, temprescale)

      use Constants, only : k_bolt, angstrom
      use InputData, only : para_rng

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), dimension(Natom), intent(in) :: lattdampvec !< Ionic damping parameter
      integer, intent(in) :: compensate_drift !< Correct for drift in RNG
      real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature (array)
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

      real(dblprec) :: D
      real(dblprec) :: mu, sigma, angstrominv
      real(dblprec), dimension(Mensemble,Natom) :: Dk
      real(dblprec) :: rx(NA), ry(NA), rz(NA)
      integer :: ity
      integer :: i, j

      angstrominv = 1.00_dblprec/angstrom

      Dk(:,:) = 0.00_dblprec
      do i=1, Natom
         Dk(:,i)= lattdampvec(i) * k_bolt !Check units!
      enddo

      if(.not.para_rng) then
#ifdef VSL
         if(use_vsl) then
            call rng_gaussianP(lattranv,3*Natom*Mensemble,1.0_dblprec)
         else
            call fill_rngarray(lattranv,3*Natom*Mensemble)
         end if
#else
         call fill_rngarray(lattranv,3*Natom*Mensemble)
#endif
      else
#ifdef VSL
         if(use_vsl) then
            call rng_gaussianP(lattranv,3*Natom*Mensemble,1.0_dblprec)
         else
            call fill_rngarray_para(lattranv,3*Natom*Mensemble)
         end if
#else
         call fill_rngarray_para(lattranv,3*Natom*Mensemble)
#endif
      end if

      !$omp parallel do default(shared) private(i,j,D,sigma,mu) collapse(2)
      do j=1, Mensemble
         do i=1, Natom
            D=Dk(j,i)*deltat*Temp_array(i)*temprescale
            sigma=sqrt(2.00_dblprec*D)
            mu=00_dblprec
            lattranv(1,i,j)=lattranv(1,i,j)*sigma
            lattranv(2,i,j)=lattranv(2,i,j)*sigma
            lattranv(3,i,j)=lattranv(3,i,j)*sigma
         end do
      end do
      !$omp end parallel do

      ! Possibility to correct the random number distribution so that
      ! the true mean is zero. Not used by default
      if(compensate_drift==1) then
         do j=1, Mensemble
            rx=0.00_dblprec;ry=0.00_dblprec;rz=0.00_dblprec
            do i=1, Natom
               ity=mod(i-1,NA)+1
               rx(ity)=rx(ity)+lattranv(1,i,j)
               ry(ity)=ry(ity)+lattranv(2,i,j)
               rz(ity)=rz(ity)+lattranv(3,i,j)
            end do
            rx=rx/Natom*NA
            ry=ry/Natom*NA
            rz=rz/Natom*NA
            do i=1, Natom
               ity=mod(i-1,NA)+1
               lattranv(1,i,j)=lattranv(1,i,j)-rx(ity)
               lattranv(2,i,j)=lattranv(2,i,j)-ry(ity)
               lattranv(3,i,j)=lattranv(3,i,j)-rz(ity)
            end do
         end do
      end if

   end subroutine lattrannum


end module RandomDrivers