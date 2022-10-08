!-------------------------------------------------------------------------------
! MODULE: MagnetizationInit
!> @brief Data and routines for initializing the magnetic moments
!> @details Set of subroutines handling the initialization of the magnetic moments
!> currently the supported initial configurations are:
!>
!>    1 Random moments
!>    2 Create random distribution of moments in an area phi0, theta0 around the z-axis.
!>    3 Read from input file
!>    4 Read the moments from the momfile
!>    5 Random Ising moments with axis dictated from the momfile
!>    6 Moments according to values in inpsd.dat file for GNEB calculations.
!>      First ensemble correspond to the initial state, last ensemble correspond to the final state
!>    7 Read magnetic configuration from file. Ensembles correspond to states to be used in the GNEB calculations
!>    8 Moments according to values in input file combined with spin spiral modulation
!
!> @author Anders Bergman, Lars Bergqvist, Johan Hellsvik, Jonathan Chico, Pavel Bessarab
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module MagnetizationInit
   use Parameters
   use Profiling

   implicit none

   public

contains

   !----------------------------------------------------------------------------
   !> @brief Initializes the magnetic moments
   !>
   !> @date 2014-08-08 - Thomas Nystrand
   !> - Modified critical values xyz for initmag=1
   !> - to be stored in an array instead of separately as x, y and z.
   !> - Previous configuration caused strange round off errors when summing x**2+y**2+z**2
   !> @date 2017-08-21 - Jonathan Chico
   !> Added fixed moments variables
   !----------------------------------------------------------------------------
   subroutine magninit(Natom,Mensemble,NA,N1,N2,N3,initmag,Nchmax,         &
      aemom_inp,anumb,do_ralloy,Natom_full,achtype,acellnumb,emom,emom2,emomM,mmom, &
      rstep,theta0,phi0,mavg0,restartfile,initrotang,initpropvec,initrotvec,coord,C1,C2,  &
      C3,do_fixed_mom,Nred,red_atom_list,ind_list_full,ind_mom_flag,ind_mom,fix_num,&
      fix_list,read_ovf,do_mom_legacy,relaxed_if)
      !
      use RandomNumbers, only : rng_uniform
      use Constants, only : pi
      use InputData, only: momfile_i,momfile_f,amp_rnd
      use FixedMom, only : create_aux_fixed_arrays
      use InducedMoments, only : setup_induced_information
#ifdef USE_OVF
      use Restart, only: read_mag_conf_ovf,GNEB_read_wrapper,read_mag_conf
#else
      use Restart, only: GNEB_read_wrapper,read_mag_conf
#endif
     use ErrorHandling
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA                 !< Number of atoms in one cell
      integer, intent(in) :: N1                 !< Number of cell repetitions in x direction
      integer, intent(in) :: N2                 !< Number of cell repetitions in y direction
      integer, intent(in) :: N3                 !< Number of cell repetitions in z direction
      integer, intent(in) :: Natom              !< Number of atoms in system
      integer, intent(in) :: Nchmax             !< Max number of chemical components on each site in cell
      integer, intent(in) :: initmag            !< Mode of initialization of magnetic moments (1-4)
      integer, intent(in) :: Mensemble          !< Number of ensembles
      integer, intent(in) :: do_ralloy          !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full         !< Number of atoms for full system (=Natom if not dilute)
      real(dblprec), intent(inout) :: phi0         !< Cone angle phi
      real(dblprec), intent(inout) :: theta0       !< Cone angle theta
      real(dblprec), intent(inout) :: mavg0        !< Average magnetization for cone
      real(dblprec), intent(in) :: initrotang   !< Rotation angle phase for initial spin spiral
      character(len=1), intent(in) :: read_ovf  !< Read the magnetization data in the ovf format
      character(len=1), intent(in) :: relaxed_if
      character(len=1), intent(in) :: do_fixed_mom !< Do Fixed moment calculation (Y/N)
      character(len=1), intent(in) :: ind_mom_flag !< Flag to indicate that there are induced moments being considered
      character(len=1), intent(in) :: do_mom_legacy
      integer, dimension(Natom) , intent(in) :: anumb             !< Atom number in cell
      integer, dimension(NA,Nchmax), intent(in) :: ind_mom        !< Indication of whether a given moment is induced/fixed (1/0) for the unit cell
      integer, dimension(:), allocatable, intent(in) :: achtype   !< Chemical type of atoms (full list)
      integer, dimension(:), allocatable, intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      real(dblprec), dimension(3), intent(in) :: C1               !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2               !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3               !< Third lattice vector
      real(dblprec), dimension(3), intent(in) :: initrotvec       !< rotation vector for initial spin spiral
      real(dblprec), dimension(3), intent(in) :: initpropvec      !< propagation vector for initial spin spiral
      real(dblprec), dimension(3,Natom), intent(in) :: coord      !< Coordinates of atoms
      real(dblprec), dimension(:,:,:,:), allocatable, intent(in) :: aemom_inp  !< Magnetic moment directions from input (for alloys)
      !.. Output variables
      integer, intent(out) :: Nred  !< Number of moments that can be updated
      integer, intent(out) :: rstep !< Starting simulation step
      integer, dimension(:), allocatable, intent(out) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom, Mensemble), intent(out) :: mmom      !< Magnitude of magnetic moments
      real(dblprec), dimension(:,:,:), allocatable, intent(out) :: emom    !< Current unit moment vector
      real(dblprec), dimension(:,:,:), allocatable, intent(out) :: emom2   !< Final (or temporary) unit moment vector
      real(dblprec), dimension(:,:,:), allocatable, intent(out) :: emomM   !< Current magnetic moment vector
      !.. In/Out variables
      integer, intent(inout) :: fix_num
      character(len=35), intent(inout) :: restartfile !< File containing restart information
      integer, dimension(Natom), intent(inout) :: ind_list_full !< Indication of whether a given moment is induced/fixed 1/0
      integer, dimension(:), allocatable, intent(inout) ::fix_list

      !.. Local scalars
      integer :: i, j
      integer :: iatom
      integer :: i0, i1, i2, i3
      integer :: i_stat,i_err,isite,ichem,i_all

      real(dblprec) :: c1r1, c2r2, c3r3
      real(dblprec) :: rotang, cosra, sinra, rotmatdet
      real(dblprec) :: theta, phi, x, y, z, mmom_tmp

      logical :: exists

      !.. Local arrays
      real(dblprec), dimension(3) :: u,v
      real(dblprec), dimension(3) :: b1,r1
      real(dblprec), dimension(3) :: b2,r2
      real(dblprec), dimension(3) :: b3,r3,rn
      real(dblprec), dimension(3) :: propvec, rotvec
      real(dblprec), dimension(3) :: tmpemom1, tmpemom2
      real(dblprec), dimension(3,3) :: rotmat
      real(dblprec), dimension(:,:,:), allocatable :: rmom

      !  Initialize moment
      !-------------------------------------------------------------------------
      ! Allocate arrays for moment directions
      !-------------------------------------------------------------------------
      allocate(emom(3,Natom,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(emom))*kind(emom),'emom','magninit')
      allocate(emomM(3,Natom,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(emomM))*kind(emomM),'emomM','magninit')
      allocate(emom2(3,Natom,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(emom2))*kind(emom2),'emom2','magninit')
      !-------------------------------------------------------------------------
      ! Random directions of moments
      !-------------------------------------------------------------------------
      if(initmag==1) then
         write (*,'(2x,a)',advance='no') "Start random init spins"
         rstep=0
         do I3=0, N3-1
            do I2=0, N2-1
               do I1=0, N1-1
                  do I0=1, NA
                     i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                     ! Call the random number generator to generate the spin directions
                     call rng_uniform(rn,3)
                     x=2.0_dblprec*(rn(1)-0.50_dblprec)
                     y=2.0_dblprec*(rn(2)-0.50_dblprec)
                     z=2.0_dblprec*(rn(3)-0.50_dblprec)
                     do while (x**2+y**2+z**2>1)
                        call rng_uniform(rn,3)
                        x=1.0_dblprec*(rn(1)-0.50_dblprec)
                        y=1.0_dblprec*(rn(2)-0.50_dblprec)
                        z=1.0_dblprec*(rn(3)-0.50_dblprec)
                     end do
                     ! Normalize the spins directions
                     u(1) = x/sqrt(x**2+y**2+z**2)
                     u(2) = y/sqrt(x**2+y**2+z**2)
                     u(3) = z/sqrt(x**2+y**2+z**2)
                     ! Check if it is a random alloy
                     if(do_ralloy==0) then
                        do j=1, Mensemble
                           emom(1:3,i,j)  = u(1:3)
                           emomM(1:3,i,j) = u(1:3)*mmom(anumb(I0),j)
                        end do
                     else
                        if (achtype(i) /= 0) then
                           iatom = acellnumb(i)
                           do j=1, Mensemble
                              emom(1:3,iatom,j)    = u(1:3)
                              emomM(1:3,iatom,j)   = u(1:3)*mmom(iatom,j)
                           end do
                        end if
                     end if
                  end do
               end do
            end do
         end do
         write (*,*) " done"
      !-------------------------------------------------------------------------
      ! Random moments within a cone
      !-------------------------------------------------------------------------
      else if (initmag==2) then
         write (*,'(2x,a)',advance='no') "Start random init spins in cone"

         if (mavg0>=0.0_dblprec) then
            phi0 = 2.0_dblprec * pi
            theta0 = 2.0_dblprec*acos(mavg0) 
            print *," MAgnetization avg:",mavg0,theta0
         end if
         rstep=0
         do I3=0, N3-1
            do I2=0, N2-1
               do I1=0, N1-1
                  do I0=1, NA
                     i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                     ! Call a random number generator to perturbate the initial angles
                     call rng_uniform(rn,2)
                     theta = rn(1)
                     phi   = rn(2)
                     ! Initialize the spins via euler angles
                     u(1) = sin(theta0*theta)*cos(phi0*phi)
                     u(2) = sin(theta0*theta)*sin(phi0*phi)
                     u(3) = cos(theta0*theta)
                     print '(3f12.5)',u
                     ! Check if the system is a random alloy
                     if(do_ralloy==0) then
                        do j=1, Mensemble
                           emom(1:3,i,j)  = u(1:3)
                           emomM(1:3,i,j) = u(1:3)*mmom(anumb(I0),j)
                        end do
                     else
                        if (achtype(i) /= 0) then
                           iatom = acellnumb(i)
                           do j=1, Mensemble
                              emom(1:3,iatom,j)    = u(1:3)
                              emomM(1:3,iatom,j)   = u(1:3)*mmom(iatom,j)
                           end do
                        end if
                     end if

                  end do
               end do
            end do
         end do
         write (*,*) " done"
         !-------------------------------------------------------------------------
         ! Moments according to values in input file
         !-------------------------------------------------------------------------
         else if(initmag==3) then
            write (*,'(2x,a)',advance='no') "Moments from inpsd.dat"
            rstep=0
            iatom=0
            do I3=0, N3-1
               do I2=0, N2-1
                  do I1=0, N1-1
                     do I0=1, NA
                        i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                        ! Check if the system is a random alloy
                        if(do_ralloy==0) then
                           do j=1, Mensemble
                              emom(1:3,i,j)  = aemom_inp(1:3,I0,1,1)
                              emomM(1:3,i,j) = aemom_inp(1:3,I0,1,1)*mmom(anumb(I0),j)
                           end do
                        else
                           if (achtype(i) /= 0) then
                              iatom = acellnumb(i)
                              ! JohanM27Jul: use same initial condition for whole ensemble
                              do j=1, Mensemble
                                 emom(1:3,iatom,j)  = aemom_inp(1:3,I0,achtype(i),1)
                                 emomM(1:3,iatom,j) = aemom_inp(1:3,I0,achtype(i),1)*mmom(iatom,j)
                              end do
                           end if
                        end if
                     end do
                  end do
               end do
            end do
            write (*,*) " done"
      !-------------------------------------------------------------------------
      ! Read moments and magnetic configuration from restart file
      !-------------------------------------------------------------------------
      else if(initmag==4) then
         if (do_fixed_mom.ne.'Y'.or.ind_mom_flag.ne.'Y') then
            if (read_ovf.ne.'Y') then
               write (*,'(2x,a)',advance='no') "Read from restart file"
               call read_mag_conf(Natom,Mensemble,do_mom_legacy,rstep,restartfile,  &
                  mmom,emom,emomM)
               write (*,'(a)') " done"
            else if (read_ovf.eq.'Y') then
#ifdef USE_OVF
               write (*,'(2x,a)',advance='no') "Read from OVF restart file"
               call read_mag_conf_ovf(Natom,Mensemble,restartfile,mmom,emom,emomM)
#else
               call ErrorHandling_ERROR('`read_ovf is Y but code is not compiled with OVF support')
#endif
            endif
         endif
      !-------------------------------------------------------------------------
      ! Random moments for Ising Hamiltonian using preferred directions from input file
      !-------------------------------------------------------------------------
      else if(initmag==5) then
         write (*,'(2x,a)',advance='no') "Random Ising Moments from inpsd.dat"
         rstep=0
         iatom=0
         do I3=0, N3-1
            do I2=0, N2-1
               do I1=0, N1-1
                  do I0=1, NA
                     i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                     ! Calculate a random number
                     call rng_uniform(rn,1)
                     ! As the spins are going to be set as Ising spins see if
                     ! one is going to flip the direction of the spin
                     if(rn(1)<0.50_dblprec) then
                        y=1
                     else
                        y=-1
                     endif
                     ! Check if the system is a random alloy
                     if(do_ralloy==0) then
                        ! JohanM27Jul: use same initial condition for whole ensemble
                        do j=1, Mensemble
                           emom(1:3,i,j)  = y*aemom_inp(1:3,I0,1,1)
                           emomM(1:3,i,j) = y*aemom_inp(1:3,I0,1,1)*mmom(anumb(I0),j)
                        end do
                     else
                        if (achtype(i) /= 0) then
                           iatom = acellnumb(i)
                           ! JohanM27Jul: use same initial condition for whole ensemble
                           do j=1, Mensemble
                              emom(1:3,iatom,j)  = y*aemom_inp(1:3,I0,achtype(i),1)
                              emomM(1:3,iatom,j) = y*aemom_inp(1:3,I0,achtype(i),1)*mmom(iatom,j)
                           end do
                        end if
                     end if
                  end do
               end do
            end do
         end do
         write (*,*) " done"
      !-------------------------------------------------------------------------
      ! Moments according to values in inpsd.dat file. First ensemble correspond to the initial state,
      ! last ensemble correspond to the final state
      !-------------------------------------------------------------------------
      else if(initmag==6) then
         write (*,'(2x,a)',advance='no') "Moments from inpsd.dat"
         rstep=0
         iatom=0
         allocate(rmom(3,NA,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(rmom))*kind(rmom),'rmom','magninit')
         do j=1, Mensemble,(Mensemble-1)
            if (j==1) then
               inquire(file=trim(momfile_i),exist=exists)
               if (.not.exists) then
                  write(*,*) 'ERROR: File ',trim(adjustl(momfile_i)), ' does not exist.'
                  stop
               end if
               open(ifileno,file=trim(momfile_i))
            elseif (j==Mensemble) then
               inquire(file=trim(momfile_f),exist=exists)
               if (.not.exists) then
                  write(*,*) 'ERROR: File ',trim(adjustl(momfile_f)), ' does not exist.'
                  stop
               end if
               open(ifileno,file=trim(momfile_f))
            end if
            i_err=0
            do while(i_err==0)
               read(ifileno,*,iostat=i_err) isite, ichem, mmom_tmp, rmom(1:3,isite,ichem)
            end do
            close(ifileno)
            do I3=0, N3-1
               do I2=0, N2-1
                  do I1=0, N1-1
                     do I0=1, NA
                        i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                        call rng_uniform(rn,3)
                        u(1)=2.0_dblprec*(rn(1)-0.50_dblprec)
                        u(2)=2.0_dblprec*(rn(2)-0.50_dblprec)
                        u(3)=2.0_dblprec*(rn(3)-0.50_dblprec)

                        if(do_ralloy==0) then
                           v(1:3) = rmom(1:3,I0,1)+amp_rnd*u(1:3)
                           mmom_tmp = norm2(v)
                           emom(1:3,i,j)  = v(1:3)/mmom_tmp
                           emomM(1:3,i,j) = emom(1:3,i,j)*mmom(anumb(I0),j)
                        else
                           if (achtype(i) /= 0) then
                              iatom = acellnumb(i)
                              v(1:3) = rmom(1:3,I0,achtype(i))+amp_rnd*u(1:3)
                              mmom_tmp = norm2(v)
                              emom(1:3,iatom,j)    = v(1:3)/mmom_tmp
                              emomM(1:3,iatom,j)   = emom(1:3,iatom,j)*mmom(iatom,j)
                           end if
                        end if
                     end do
                  end do
               end do
            end do
         end do
         i_all=-product(shape(rmom))*kind(rmom)
         deallocate(rmom,stat=i_stat)
         call memocc(i_stat,i_all,'rmom','magninit')

         write (*,*) " done"
      !-------------------------------------------------------------------------
      ! Read magnetic configuration from file. Ensembles correspond to states to be used in the GNEB calculations
      !-------------------------------------------------------------------------
      else if (initmag==7) then
         write (*,'(2x,a)',advance='no') "Read moments from file"
         call GNEB_read_wrapper(Natom,Mensemble,amp_rnd,'infi',relaxed_if,          &
            do_mom_legacy,rstep,exists,restartfile,mmom,emom,emomM)
         write (*,*) " done"
      !-------------------------------------------------------------------------
      ! Moments according to values in input file combined with spin spiral modulation
      !-------------------------------------------------------------------------
      else if(initmag==8) then
         write (*,'(2x,a)') "Moments from inpsd.dat combined with"
         write (*,'(2x,a)',advance='no') "spin spiral modulation"

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

         ! Propagation vector
         propvec(1)=2*pi*( initpropvec(1)*b1(1)+initpropvec(2)*b2(1)+initpropvec(3)*b3(1) )
         propvec(2)=2*pi*( initpropvec(1)*b1(2)+initpropvec(2)*b2(2)+initpropvec(3)*b3(2) )
         propvec(3)=2*pi*( initpropvec(1)*b1(3)+initpropvec(2)*b2(3)+initpropvec(3)*b3(3) )

         ! Normalized rotation vector
         rotvec(1:3) = initrotvec(1:3)/norm2(initrotvec)

         write (*,'(2x,a,3f16.8)') "Init propagation vector", initpropvec(1:3)
         write (*,'(2x,a,3f16.8)') "Propagation vector     ", propvec(1:3)
         write (*,'(2x,a,3f16.8)') "Init rotation vector   ", initrotvec(1:3)
         write (*,'(2x,a,3f16.8)') "Rotation vector        ", rotvec(1:3)

         rstep=0
         iatom=0
         do I3=0, N3-1
            do I2=0, N2-1
               do I1=0, N1-1
                  do I0=1, NA
                     i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                     !----------------------------------------------------------
                     ! The rotation angle is calculated from the scalar product of the
                     ! propagation vector and the position of the atom. The phase can
                     ! be shifted by using a finite value for initrotang
                     !----------------------------------------------------------
                     rotang = propvec(1)*coord(1,i) + propvec(2)*coord(2,i) &
                        + propvec(3)*coord(3,i) + initrotang * pi/180.0_dblprec
                     !----------------------------------------------------------
                     ! The rotation matrix is calculated from the rotation axis and the
                     ! rotation angle
                     !----------------------------------------------------------
                     cosra = cos(rotang)
                     sinra = sin(rotang)
                     rotmat(1,1) = (1-cosra)*rotvec(1)*rotvec(1) + cosra
                     rotmat(1,2) = (1-cosra)*rotvec(1)*rotvec(2) - sinra*rotvec(3)
                     rotmat(1,3) = (1-cosra)*rotvec(1)*rotvec(3) + sinra*rotvec(2)
                     rotmat(2,1) = (1-cosra)*rotvec(1)*rotvec(2) + sinra*rotvec(3)
                     rotmat(2,2) = (1-cosra)*rotvec(2)*rotvec(2) + cosra
                     rotmat(2,3) = (1-cosra)*rotvec(2)*rotvec(3) - sinra*rotvec(1)
                     rotmat(3,1) = (1-cosra)*rotvec(1)*rotvec(3) - sinra*rotvec(2)
                     rotmat(3,2) = (1-cosra)*rotvec(2)*rotvec(3) + sinra*rotvec(1)
                     rotmat(3,3) = (1-cosra)*rotvec(3)*rotvec(3) + cosra
                     rotmatdet = &
                        rotmat(1,1)*(rotmat(2,2)*rotmat(3,3)-rotmat(2,3)*rotmat(3,2)) + &
                        rotmat(1,2)*(rotmat(2,3)*rotmat(3,1)-rotmat(2,1)*rotmat(3,3)) + &
                        rotmat(1,3)*(rotmat(2,1)*rotmat(3,2)-rotmat(2,2)*rotmat(3,1))
                     ! Calculate the rotated moments
                     tmpemom1(1:3) = aemom_inp(1:3,I0,1,1)
                     tmpemom2(1) = rotmat(1,1)*tmpemom1(1) + rotmat(1,2)*tmpemom1(2) + rotmat(1,3)*tmpemom1(3)
                     tmpemom2(2) = rotmat(2,1)*tmpemom1(1) + rotmat(2,2)*tmpemom1(2) + rotmat(2,3)*tmpemom1(3)
                     tmpemom2(3) = rotmat(3,1)*tmpemom1(1) + rotmat(3,2)*tmpemom1(2) + rotmat(3,3)*tmpemom1(3)

                     if(do_ralloy==0) then
                        do j=1, Mensemble
                           emom(1:3,i,j) = tmpemom2(1:3)
                           emomM(1:3,i,j) = tmpemom2(1:3)*mmom(anumb(I0),j)
                        end do
                     else
                        if (achtype(i) /= 0) then
                           iatom = acellnumb(i)
                           do j=1, Mensemble
                              emom(1:3,iatom,j) = tmpemom2(1:3)
                              emomM(1:3,iatom,j) = tmpemom2(1:3)*mmom(anumb(I0),j)
                           end do
                        end if
                     end if
                  end do
               end do
            end do
         end do
         write (*,*) " done"
      endif
      !-------------------------------------------------------------------------
      ! Set all the information needed for the fixed moments calculation
      !-------------------------------------------------------------------------
      if (do_fixed_mom.eq.'Y') then
         write (*,'(2x,a)',advance='no') "Setting up for fixed moments "
         call create_aux_fixed_arrays(Natom,Mensemble,initmag,Natom_full,do_ralloy, &
            anumb,achtype,rstep,red_atom_list,Nred,mmom,emom,emomM,restartfile,     &
            ind_mom_flag,ind_list_full,do_mom_legacy)
         write(*,'(a)') "done"
      else
         Nred=Natom
         allocate(red_atom_list(Nred),stat=i_stat)
         call memocc(i_stat,product(shape(red_atom_list))*kind(red_atom_list),'red_atom_list','magninit')
         red_atom_list=0
         ! Fill up the array where each entry in the index
         do iatom=1,Natom
            red_atom_list(iatom)=iatom
         enddo
      endif
      !-------------------------------------------------------------------------
      ! Set all the information needed for the induced moments calculation
      !-------------------------------------------------------------------------
      if (ind_mom_flag.eq.'Y'.and.do_fixed_mom.ne.'Y') then
         write (*,'(2x,a)',advance='no') "Setting up for induced moments "
         call setup_induced_information(NA,Natom,Nchmax,do_ralloy,Natom_full,anumb, &
            achtype,ind_mom,ind_list_full,fix_list,fix_num,restartfile,rstep,mmom,  &
            emom,emomM,initmag,Mensemble,do_mom_legacy)
         write(*,'(a)') "done"
      endif

      emom2=emom

   end subroutine magninit

   !----------------------------------------------------------------------------
   !> Set up the magnitude of magnetic moments, and Lande factors as well.
   !----------------------------------------------------------------------------
   subroutine setup_moment(Natom,Mensemble,NA,N1,N2,N3, ammom_inp,  &
      Landeg_ch,Landeg,mmom,mmom0,mmomi,do_ralloy,achtype, acellnumb,   &
      mconf)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      real(dblprec), intent(in), dimension(:,:,:) :: ammom_inp  !< Magnetic moment magnitudes from input (for alloys)
      real(dblprec), intent(in), dimension(:,:,:) :: Landeg_ch  !< Gyromagnetic ratio from input
      real(dblprec), dimension(:) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom, Mensemble) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(Natom, Mensemble) :: mmom0 !< Starting magnitude of magnetic moments
      real(dblprec), dimension(Natom, Mensemble) :: mmomi !< Inverse of magnitude of magnetic moments
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, dimension(:), intent(in) :: achtype !< Chemical type of atoms (full list)
      integer, dimension(:), intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      integer, intent(in) :: mconf !< Starting LSF configuration if moments
      !
      integer :: i
      integer :: iatom
      integer :: i0, i1, i2, i3
      !
      iatom=0
      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA
                  i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA

                  if (do_ralloy==0) then
                     mmom(i,:)=abs(ammom_inp(I0,1,mconf))
                     mmom0(i,:)=abs(ammom_inp(I0,1,mconf))
                     mmomi(i,:)=1.0_dblprec/abs(ammom_inp(I0,1,mconf))
                     !Division with 2 to store in units
                     !of g=2 (the default Lande g)
                     Landeg(i)=Landeg_ch(I0,1,1)*0.5_dblprec
                  else
                     if (achtype(i) /= 0) then
                        iatom = acellnumb(i)
                        mmom(iatom,:)=abs(ammom_inp(I0,achtype(i),mconf))
                        mmom0(iatom,:)=abs(mmom(iatom,1))
                        mmomi(iatom,:)=1.0_dblprec/abs(mmom(iatom,1))
                        !Division with 2 to store in units
                        !of g=2 (the default Lande g)
                        Landeg(iatom)=Landeg_ch(I0,achtype(i),1)*0.5_dblprec
                     end if
                  end if

               end do
            end do
         end do
      end do
   end subroutine setup_moment

   !----------------------------------------------------------------------------
   !> Rotates initial spin configuration using Euler angles
   !> @details For details on the Eusler angles one can look at the document
   !> mathworld.wolfram.com/EulerAngles.html
   !> ( /up/hellsvik/SD/RotationGroup/EulerAngles.html )
   !> or in Goldstein, Classical Mechanics
   !> Uses the x-convention
   !> roteul == 1: Rotates with Euler angles phi, theta psi
   !> roteul == 2: Rotates to x-axis, thereafter rotation
   !> with Euler angles phi, theta, psi
   !> @author Johann Hellsvik
   !----------------------------------------------------------------------------
   subroutine rotationeuler(Natom, Mensemble,roteul, rotang,emom, emomM, mmom, emom2)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: roteul !< Global rotation of magnetic moments (0/1)
      real(dblprec), dimension(3), intent(in) :: rotang !< Euler angles for global rotation of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble),intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble),intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble),intent(out) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(Natom,Mensemble),intent(in) :: mmom !< Magnitude of magnetic moments

      !.. Scalar variables
      real(dblprec) :: alpha,beta
      real(dblprec) :: phi,theta,psi
      real(dblprec) :: raddeg
      integer :: m
      integer :: i,j,k

      !.. Array variables
      real(dblprec), dimension(Mensemble) :: rdet
      real(dblprec), dimension(3,3,3,Mensemble) :: r
      real(dblprec),dimension(3,Natom,Mensemble) :: emom0
      real(dblprec),dimension(3,Natom,Mensemble) :: emomM0
      real(dblprec),dimension(3,Mensemble) :: emomM0av
      real(dblprec),dimension(3,Mensemble) :: emomM1av

      !Stores original initial configuration to emom0(1:3,m)
      emomM0av=0
      do m=1,Natom
         emom0(1:3,m,:) = emom(1:3,m,:)
         emomM0(1:3,m,:) = emomM(1:3,m,:)
         emomM0av(1:3,:) = emomM0av(1:3,:) + emomM0(1:3,m,:)
      end do

      emomM0av(1:3,:) = emomM0av(1:3,:)/Natom
      raddeg = 3.1415/180
      phi = rotang(1)*raddeg
      theta = rotang(2)*raddeg
      psi = rotang(3)*raddeg
      r(1,1,2,:) = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)
      r(1,2,2,:) = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)
      r(1,3,2,:) = sin(psi)*sin(theta)
      r(2,1,2,:) = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)
      r(2,2,2,:) = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)
      r(2,3,2,:) = cos(psi)*sin(theta)
      r(3,1,2,:) = sin(theta)*sin(phi)
      r(3,2,2,:) = -sin(theta)*cos(phi)
      r(3,3,2,:) = cos(theta)

      if (roteul==2) then
         r(:,:,3,:) = 0
         do m=1,Mensemble
            alpha = atan2(emomM0av(2,m),emomM0av(1,m))
            beta = atan2(emomM0av(3,m),sqrt(emomM0av(1,m)**2+emomM0av(2,m)**2))
            r(1,1,1,m) = cos(beta)*cos(alpha)
            r(1,2,1,m) = cos(beta)*sin(alpha)
            r(1,3,1,m) = sin(beta)
            r(2,1,1,m) = -sin(alpha)
            r(2,2,1,m) = cos(alpha)
            r(2,3,1,m) = 0
            r(3,1,1,m) = -sin(beta)*cos(alpha)
            r(3,2,1,m) = -sin(beta)*sin(alpha)
            r(3,3,1,m) = cos(beta)
            do i=1,3
               do j=1,3
                  do k=1,3
                     r(i,j,3,m) = r(i,j,3,m)+r(i,k,2,m)*r(k,j,1,m)
                  end do
               end do
            end do
         end do
      else
         r(:,:,3,:) = r(:,:,2,:)
      end if

      do j=1, Mensemble
         rdet(j) = &
            r(1,1,3,j)*(r(2,2,3,j)*r(3,3,3,j)-r(2,3,3,j)*r(3,2,3,j)) + &
            r(1,2,3,j)*(r(2,3,3,j)*r(3,1,3,j)-r(2,1,3,j)*r(3,3,3,j)) + &
            r(1,3,3,j)*(r(2,1,3,j)*r(3,2,3,j)-r(2,2,3,j)*r(3,1,3,j))
      end do

      emomM1av=0
      do m=1,Natom
         emom(1,m,:) = r(1,1,3,:)*emom0(1,m,:) + r(1,2,3,:)*emom0(2,m,:) + &
            r(1,3,3,:)*emom0(3,m,:)
         emom(2,m,:) = r(2,1,3,:)*emom0(1,m,:) + r(2,2,3,:)*emom0(2,m,:) + &
            r(2,3,3,:)*emom0(3,m,:)
         emom(3,m,:) = r(3,1,3,:)*emom0(1,m,:) + r(3,2,3,:)*emom0(2,m,:) + &
            r(3,3,3,:)*emom0(3,m,:)
      end do

      do m=1,Natom
         do j=1, Mensemble
            emomM(1:3,m,j) = emom(1:3,m,j)*mmom(m,j)
            emomM1av(1:3,j) = emomM1av(1:3,j) + emomM(1:3,m,j)
            emom2(1:3,m,j) = emom(1:3,m,j)
         end do
      end do

      emomM1av(1:3,:) = emomM1av(1:3,:)/Natom

      write (*,*) "Euler rotation angles (x-convention)"
      write (*,10001) phi/raddeg, theta/raddeg, psi/raddeg
      write (*,10002) phi, theta, psi
      write (*,*) "Rotation matrix elements"
      do j=1, Mensemble
         write (*,*) "Sample ", j
         write (*,10011) r(1,1,3,:),r(1,2,3,:),r(1,3,3,:)
         write (*,10012) r(2,1,3,:),r(2,2,3,:),r(2,3,3,:)
         write (*,10013) r(3,1,3,:),r(3,2,3,:),r(3,3,3,:)
      end do
      write (*,*) "Rotation matrix determinant"
      do j=1, Mensemble
         write (*,*) "Sample ", j
         write (*,10014) rdet
      end do
      do j=1, Mensemble
         write (*,*) "Sample ", j
         write (*,*) "Original initial configuration (average magnetic moment)"
         write (*,10021) emomM0av(1,j), emomM0av(2,j), emomM0av(3,j)
         write (*,*) "Rotated initial configuration (average magnetic moment)"
         write (*,10022) emomM1av(1,j), emomM1av(2,j), emomM1av(3,j)
      end do

      10001 format ("phi      ",f8.4,"    theta      ",f8.4,"    psi      ",f8.4)
      10002 format ("phi(rad) ",f8.4,"    theta(rad) ",f8.4,"    psi(rad) ",f8.4)
      10011 format ("a11",f8.4,"    a12",f8.4,"    a13",f8.4)
      10012 format ("a21",f8.4,"    a22",f8.4,"    a23",f8.4)
      10013 format ("a31",f8.4,"    a32",f8.4,"    a33",f8.4)
      10014 format ("rdet  ",f8.4)
      10021 format ("x0 ",f8.4,"    y0 ",f8.4,"    z0 ",f8.4)
      10022 format ("x1 ",f8.4,"    y1 ",f8.4,"    z1 ",f8.4)

   end subroutine rotationeuler

   !----------------------------------------------------------------------------
   !> @brief Excites initial spin configurations
   !----------------------------------------------------------------------------
   subroutine setinitexc(Natom, Mensemble, emom, emomM, mmom, mmom2, emom2, initexc,&
         initconc, initneigh, initimp, max_no_neigh, nlist, mseed)

      use stdtypes
      use mtprng

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble),intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble),intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble),intent(out) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(Natom,Mensemble),intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble),intent(inout) :: mmom2 !< Final (or temporary) magnitude of magnetic moments
      character(len=1), intent(in) :: initexc !< Mode of excitation of initial magnetic moments (I=vacancies, R=two magnon Raman, F=no)
      real(dblprec), intent(in)  :: initconc  !< Concentration of vacancies or two magnon Raman spin flips
      integer, intent(in)  :: initneigh !< Raman neighbour spin index
      real(dblprec), intent(in) :: initimp !< Size of impurity magnetic moment
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(max_no_neigh, Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer :: mseed !< Seed for initialization of magnetic moments

      !.. Scalar variables
      integer :: m
      integer :: i

      !.. Array variables
      real(dblprec),dimension(3,Natom,Mensemble) :: emom0
      real(dblprec),dimension(3,Natom,Mensemble) :: emomM0
      real(dblprec),dimension(3,Mensemble) :: emomM0av
      real(dblprec),dimension(3,Mensemble) :: emomM1av

      type(mtprng_state) :: state
      real(dblprec) :: x
      real(dblprec), dimension(3,Mensemble) :: tmpmom

      call mtprng_init(mseed, state)

      !Stores original initial configuration to emom0(1:3,m)
      emomM0av=0
      do m=1,Natom
         emom0(1:3,m,:) = emom(1:3,m,:)
         emomM0(1:3,m,:) = emomM(1:3,m,:)
         emomM0av(1:3,:) = emomM0av(1:3,:) + emomM0(1:3,m,:)
      end do
      emomM0av(1:3,:) = emomM0av(1:3,:)/Natom

      if(initexc == 'I') then
         emomM1av=0
         do i=1,Natom
            do m=1, Mensemble
               x=mtprng_rand_real3(state)
               if(x .lt. initconc) then
                  mmom(i,m) = initimp
               end if
            end do
         end do
      end if

      if(initexc == 'R') then
         emomM1av=0
         do i=1,Natom
            do m=1, Mensemble
               x=mtprng_rand_real3(state)
               if(x .lt. initconc) then
                  tmpmom(1:3,m) = emom(1:3,i,m)
                  emom(1:3,i,m) = emom(1:3,nlist(initneigh,i),m)
                  emom(1:3,nlist(initneigh,i),m) = tmpmom(1:3,m)
               end if
            end do
         end do
      end if

      do i=1,Natom
         do m=1, Mensemble
            emomM(1:3,i,m) = emom(1:3,i,m)*mmom(i,m)
            emomM1av(1:3,m) = emomM1av(1:3,m) + emomM(1:3,i,m)
            emom2(1:3,i,m) = emom(1:3,i,m)
            mmom2(i,m) = mmom(i,m)
         end do
      end do
      emomM1av(1:3,:) = emomM1av(1:3,:)/Natom

   end subroutine setinitexc

end module MagnetizationInit
