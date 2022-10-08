!> Data and routines for initializing the ionic displacements
module LatticeInit
   use Parameters
   use Profiling

   implicit none

   public


contains


   subroutine lattinit(Natom, Mensemble, NA, N1, N2, N3, initlatt, Nchmax, mion_inp, uvec_inp, vvec_inp, anumb, do_ralloy, &
         Natom_full, achtype, acellnumb, mion, mioninv, uvec, vvec, rstep, restartfile)
      !
      use RandomNumbers, only : rng_uniform
      !
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: initlatt !< Mode of initialization of ionic displacements (1-4)
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      real(dblprec), intent(in), dimension(NA,Nchmax) :: mion_inp  !< Ionic mass from input (for alloys)
      real(dblprec), intent(in), dimension(3,NA,Nchmax) :: uvec_inp  !< Ionic displacements from input (for alloys)
      real(dblprec), intent(in), dimension(3,NA,Nchmax) :: vvec_inp  !< Ionic velocities from input (for alloys)
      integer, dimension(Natom) , intent(in) :: anumb !< Atom number in cell
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, dimension(Natom), intent(in) :: achtype !< Chemical type of atoms (full list)
      integer, dimension(Natom), intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      real(dblprec), dimension(Natom,Mensemble), intent(out) :: mion   !< Ionic mass
      real(dblprec), dimension(Natom,Mensemble), intent(out) :: mioninv   !< Inverse ionic mass
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: uvec   !< Current ionic displacment
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: vvec   !< Current ionic velocity
      integer, intent(out) :: rstep !< Starting simulation step
      character(len=35) :: restartfile !< File containing restart information

      !.. Local scalars
      integer :: i, j
      integer :: iatom
      integer :: i0, i1, i2, i3

      !.. Local arrays
      real(dblprec), dimension(3) :: u, rn
      real(dblprec) :: x, y, z


      !
      !  Initialize ionic displacements and velocities
      !
      !     1       Random ionic displacements and velocities
      !
      !     3       Read ionic displacements and velocities from input file
      !
      !     4       Read ionic displacements and velocities from start file
      !

      write (*,'(2x,a)',advance='no') "Ionic masses from input file"
      iatom=0
      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA
                  i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                  if(do_ralloy==0) then
                     do j=1, Mensemble
                        mion(i,j) = mion_inp(I0,1)
                        mioninv(i,j) = 1_dblprec/mion(i,j)
                     end do
                  else
                     if (achtype(i) /= 0) then
                        iatom = acellnumb(i)
                        do j=1, Mensemble
                           mion(iatom,j) = mion_inp(achtype(i),1)
                           mioninv(iatom,j) = 1.0_dblprec/mion(iatom,j)
                        end do
                     end if
                  end if
               end do
            end do
         end do
      end do
      write (*,*) " done"

      ! Read ionic displacments and velocities from restart file
      if(initlatt==4) then

         write (*,'(2x,a)',advance='no') "Read from lattice restart file"
         call loadlattrestart(Natom, Mensemble, restartfile, rstep, uvec, vvec)
         write (*,*) " done"


         ! Read ionic displacments and velocities from lattice input file
      else if(initlatt==3) then

         write (*,'(2x,a)',advance='no') "Displacements and velocities from input file"
         rstep=0
         iatom=0
         do I3=0, N3-1
            do I2=0, N2-1
               do I1=0, N1-1
                  do I0=1, NA
                     i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                     if(do_ralloy==0) then
                        do j=1, Mensemble
                           uvec(1:3,i,j) = uvec_inp(1:3,I0,1)
                           vvec(1:3,i,j) = vvec_inp(1:3,I0,1)
                        end do
                     else
                        if (achtype(i) /= 0) then
                           iatom = acellnumb(i)
                           do j=1, Mensemble
                              uvec(1:3,iatom,j) = uvec_inp(1:3,I0,achtype(i))
                              vvec(1:3,iatom,j) = vvec_inp(1:3,I0,achtype(i))
                           end do
                        end if
                     end if
                  end do
               end do
            end do
         end do
         write (*,*) " done"


         ! Random displacements
      else if(initlatt==1) then

         write (*,'(2x,a)',advance='no') "Start random displacements"
         rstep=0
         do I3=0, N3-1
            do I2=0, N2-1
               do I1=0, N1-1
                  do I0=1, NA
                     i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA

                     ! SLDTODO Normalise random displacements in this fashion?
                     call rng_uniform(rn,3)
                     x=2_dblprec*(rn(1)-0.5)
                     y=2_dblprec*(rn(2)-0.5)
                     z=2_dblprec*(rn(3)-0.5)
                     do while (x**2+y**2+z**2>1)
                        call rng_uniform(rn,3)
                        x=rn(1)-0.5
                        y=rn(2)-0.5
                        z=rn(3)-0.5
                     end do
                     u(1) = 0.001*x/sqrt(x**2+y**2+z**2)
                     u(2) = 0.001*y/sqrt(x**2+y**2+z**2)
                     u(3) = 0.001*z/sqrt(x**2+y**2+z**2)
                     !u(1) = x/sqrt(x**2+y**2+z**2)
                     !u(2) = y/sqrt(x**2+y**2+z**2)
                     !u(3) = z/sqrt(x**2+y**2+z**2)

                     ! Let initial velocities be zero
                     !x=2_dblprec*(mt_ran_b()-0.5)
                     !y=2_dblprec*(mt_ran_b()-0.5)
                     !z=2_dblprec*(mt_ran_b()-0.5)
                     !do while (x**2+y**2+z**2>1)
                     !   x=mt_ran_b()-0.5
                     !   y=mt_ran_b()-0.5
                     !   z=mt_ran_b()-0.5
                     !end do
                     !v(1) = x/sqrt(x**2+y**2+z**2)
                     !v(2) = y/sqrt(x**2+y**2+z**2)
                     !v(3) = z/sqrt(x**2+y**2+z**2)

                     if(do_ralloy==0) then
                        do j=1, Mensemble
                           uvec(1:3,i,j) = u(1:3)
                           !vvec(1:3,i,j) = v(1:3)
                        end do
                     else
                        if (achtype(i) /= 0) then
                           iatom = acellnumb(i)
                           do j=1, Mensemble
                              uvec(1:3,iatom,j) = u(1:3)
                              !vvec(1:3,iatom,j) = v(1:3)
                           end do
                        end if
                     end if
                  end do
               end do
            end do
         end do
         ! Let initial velocities be zero
         vvec = 0.0_dblprec

         write (*,*) " done"

      endif

   end subroutine lattinit


   !> Read magnetic moments from file
   subroutine loadlattrestart(Natom, Mensemble, restartfile, rstep, uvec, vvec)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=35), intent(inout) :: restartfile !< File containing restart information
      integer, intent(out) :: rstep !< Starting simulation step
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: uvec !< Current ionic displacment
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: vvec !< Current ionic velocity

      integer :: i, j, k, l, ios
      logical :: exists

      !.. Executable statements
      !SLDTODE fix to use the filename passed as an argument
      !restartfile='lattrestart.dat'
      inquire(file=restartfile,exist=exists)
      if(exists) then
         open(12, iostat=ios, file=restartfile, status="old")
         read (12,*) rstep
         do i=1,Mensemble
            do j=1, Natom
               read (12,*) k, l, uvec(1,j,i), uvec(2,j,i), uvec(3,j,i), &
                  vvec(1,j,i), vvec(2,j,i), vvec(3,j,i)
            end do
         end do
         close(12)
      else
         write(*,*) 'ERROR: Restartfile ',trim(adjustl(restartfile)), ' does not exist.'
         stop
      end if

   end subroutine loadlattrestart


   !> Rotates initial ionic displacement and velocities configuration using
   !! Euler angles
   !  as described in the document
   !  mathworld.wolfram.com/EulerAngles.html
   !  ( /up/hellsvik/SD/RotationGroup/EulerAngles.html )
   !  or in Goldstein, Classical Mechanics
   !  Uses the x-convention
   !  lattroteul == 1: Rotates with Euler angles phi, theta psi
   !  lattroteul == 2: Rotates to x-axis, thereafter rotation
   !  with Euler angles phi, theta, psi
   subroutine lattrotationeuler(Natom, Mensemble, lattroteul, lattrotang, uvec, vvec)
      !subroutine lattrotationeuler(Natom, Mensemble, lattroteul, lattrotang, uvec, uvec2, vvec, vvec2)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: lattroteul !< Global rotation of magnetic moments (0/1/2)
      real(dblprec), dimension(3), intent(in) :: lattrotang !< Euler angles for global rotation of magnetic moments

      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: uvec   !< Current ionic displacment
      !real(dblprec), dimension(:,:,:), allocatable, intent(out) :: uvec2  !< Final (or temporary) ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble),intent(inout) :: vvec   !< Current ionic velocity
      !real(dblprec), dimension(:,:,:), allocatable, intent(out) :: vvec2  !< Final (or temporary) ionic velocity

      !.. Scalar variables
      real(dblprec) :: alpha,beta
      real(dblprec) :: phi,theta,psi
      real(dblprec) :: raddeg
      integer :: m
      integer :: i,j,k


      !.. Array variables
      real(dblprec), dimension(Mensemble) :: rdet
      real(dblprec), dimension(3,3,3,Mensemble) :: r
      real(dblprec), dimension(3,Natom,Mensemble) :: uvec0
      real(dblprec), dimension(3,Natom,Mensemble) :: vvec0
      real(dblprec), dimension(3,Mensemble) :: uvec0av
      real(dblprec), dimension(3,Mensemble) :: vvec0av
      real(dblprec), dimension(3,Mensemble) :: uvec1av
      real(dblprec), dimension(3,Mensemble) :: vvec1av

      !Stores original initial configuration to uvec0 and vvec0
      uvec0av = 0_dblprec
      vvec0av = 0_dblprec
      do m=1,Natom
         uvec0(1:3,m,:) = uvec(1:3,m,:)
         vvec0(1:3,m,:) = vvec(1:3,m,:)
         uvec0av(1:3,:) = uvec0av(1:3,:) + uvec0(1:3,m,:)
         vvec0av(1:3,:) = vvec0av(1:3,:) + vvec0(1:3,m,:)
      end do

      uvec0av(1:3,:) = uvec0av(1:3,:)/Natom
      raddeg = 3.1415/180
      phi = lattrotang(1)*raddeg
      theta = lattrotang(2)*raddeg
      psi = lattrotang(3)*raddeg
      r(1,1,2,:) = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)
      r(1,2,2,:) = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)
      r(1,3,2,:) = sin(psi)*sin(theta)
      r(2,1,2,:) = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)
      r(2,2,2,:) = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)
      r(2,3,2,:) = cos(psi)*sin(theta)
      r(3,1,2,:) = sin(theta)*sin(phi)
      r(3,2,2,:) = -sin(theta)*cos(phi)
      r(3,3,2,:) = cos(theta)

      if (lattroteul==2) then
         r(:,:,3,:) = 0
         do m=1,Mensemble
            alpha = atan2(uvec0av(2,m),uvec0av(1,m))
            beta = atan2(uvec0av(3,m),sqrt(uvec0av(1,m)**2+uvec0av(2,m)**2))
            r(1,1,1,m) = cos(beta)*cos(alpha)
            r(1,2,1,m) = cos(beta)*sin(alpha)
            r(1,3,1,m) = sin(beta)
            r(2,1,1,m) = -sin(alpha)
            r(2,2,1,m) = cos(alpha)
            r(2,3,1,m) = 0_dblprec
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

      uvec1av=0_dblprec
      vvec1av=0_dblprec
      do m=1,Natom
         uvec(1,m,:) = r(1,1,3,:)*uvec0(1,m,:) + r(1,2,3,:)*uvec0(2,m,:) + &
            r(1,3,3,:)*uvec0(3,m,:)
         uvec(2,m,:) = r(2,1,3,:)*uvec0(1,m,:) + r(2,2,3,:)*uvec0(2,m,:) + &
            r(2,3,3,:)*uvec0(3,m,:)
         uvec(3,m,:) = r(3,1,3,:)*uvec0(1,m,:) + r(3,2,3,:)*uvec0(2,m,:) + &
            r(3,3,3,:)*uvec0(3,m,:)
         vvec(1,m,:) = r(1,1,3,:)*vvec0(1,m,:) + r(1,2,3,:)*vvec0(2,m,:) + &
            r(1,3,3,:)*vvec0(3,m,:)
         vvec(2,m,:) = r(2,1,3,:)*vvec0(1,m,:) + r(2,2,3,:)*vvec0(2,m,:) + &
            r(2,3,3,:)*vvec0(3,m,:)
         vvec(3,m,:) = r(3,1,3,:)*vvec0(1,m,:) + r(3,2,3,:)*vvec0(2,m,:) + &
            r(3,3,3,:)*vvec0(3,m,:)
      end do

      uvec1av(1:3,:) = uvec1av(1:3,:)/Natom
      vvec1av(1:3,:) = vvec1av(1:3,:)/Natom

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
         write (*,*) "Original initial configuration (average ionic displacement)"
         write (*,10021) uvec0av(1,j), uvec0av(2,j), uvec0av(3,j)
         write (*,*) "Rotated initial configuration (average ionic displacement)"
         write (*,10022) uvec1av(1,j), uvec1av(2,j), uvec1av(3,j)
      end do


      10001 format ("phi      ",f8.4,"    theta      ",f8.4,"    psi      ",f8.4)
      10002 format ("phi(rad) ",f8.4,"    theta(rad) ",f8.4,"    psi(rad) ",f8.4)
      10011 format ("a11",f8.4,"    a12",f8.4,"    a13",f8.4)
      10012 format ("a21",f8.4,"    a22",f8.4,"    a23",f8.4)
      10013 format ("a31",f8.4,"    a32",f8.4,"    a33",f8.4)
      10014 format ("rdet  ",f8.4)
      10021 format ("x0 ",f8.4,"    y0 ",f8.4,"    z0 ",f8.4)
      10022 format ("x1 ",f8.4,"    y1 ",f8.4,"    z1 ",f8.4)

   end subroutine lattrotationeuler

end module LatticeInit
