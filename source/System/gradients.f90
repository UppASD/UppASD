!------------------------------------------------------------------------------------
!> @brief Routines for calculating dm/dr. Needed for generalized spin-torque term
!> \f$\left(\mathbf{m} \times \left(\mathbf{m}\times\frac{\partial\mathbf{m}}{\partial \mathbf{r}}\right)\right)\f$
!> @details For calculating the spin transfer torques, we use the standard adiabatic
!> and non-adiabatic terms as introduced in the LLG equation by Zhang & Li. PRL 93, 127204 (2004).
!> The terms are rewritten to suit the LL equations used in UppASD and here we use
!> the same formulas as Schieback et. al, Eur. Phys. J. B 59, 429, (2007)
!> Currently the torques are only calculated as their respective fields
!> (i.e. missing the preceeding \f$\mathbf{m}\times\f$) since that is taken care of in the Depondt solver
!> @authors 
!> A. Bergman
!> E. Mendez
!> N. Ntallis
!> M. Pereiro  
!> @copyright
!> GNU Public License.
!> @date August 2010 / April 2011
!> @todo Strenght of prefactors still not controlled
!------------------------------------------------------------------------------------
module Gradients

   use Parameters
   use Profiling
   use InputData, only: do_multiscale

   implicit none
   real(dblprec), dimension(:,:,:), allocatable :: dxyz_vec
   integer, dimension(:,:), allocatable :: dxyz_atom
   integer, dimension(:), allocatable :: dxyz_list


   private

   public :: differentiate_moments, setup_stencil_mesh, grad_moments, proj_grad_moments, deallocate_gradient_lists
   public :: dxyz_vec, dxyz_atom, dxyz_list

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: setup_stencil_mesh
   !> Find nearest neighbours in 3D stencil shape. Used for discrete differentiation.
   !---------------------------------------------------------------------------------
   subroutine setup_stencil_mesh(Natom, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3,      &
         max_no_neigh, nlistsize, nlist, coord)

      implicit none

      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh, Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms

      ! .. Local variables
      real(dblprec), parameter :: rcutoff=1.05_dblprec

      integer :: iatom, jneigh, jatom, ix, iy, iz, icount
      integer :: i_stat
      integer :: signx, signy, signz
      integer :: perx, pery, perz
      real(dblprec) :: jnorm, jpnorm, jtnorm
      real(dblprec), dimension(3) :: icoord, jcoord, jpcoord, jtcoord
      real(dblprec), dimension(3,-1:1,-1:1,-1:1) :: minmat_coord
      real(dblprec), dimension(-1:1,-1:1,-1:1) :: minmat_norm
      integer, dimension(-1:1,-1:1,-1:1) :: minmat_atom

      ! Allocate storage arrays
      if (.not.allocated(dxyz_list)) then
         allocate(dxyz_list(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(dxyz_list))*kind(dxyz_list),'dxyz_list','setup_stencil_mesh')
         dxyz_list=0
      endif
      if (.not.allocated(dxyz_atom)) then
         allocate(dxyz_atom(26,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(dxyz_atom))*kind(dxyz_atom),'dxyz_atom','setup_stencil_mesh')
         dxyz_atom=0
      endif
      if (.not.allocated(dxyz_vec)) then
         allocate(dxyz_vec(3,26,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(dxyz_vec))*kind(dxyz_vec),'dxyz_vec','setup_stencil_mesh')
         dxyz_vec=0.0_dblprec
      endif

      ! Prepare for periodic boundaries
      if(BC1=='P') then
         perx=1
      else
         perx=0
      end if
      if(BC2=='P') then
         pery=1
      else
         pery=0
      end if
      if(BC3=='P') then
         perz=1
      else
         perz=0
      end if

      ! Loop over atoms
      do iatom=1, Natom
         icoord=coord(1:3,iatom)

         ! Initialize minima matrix with large entries
         minmat_norm=rcutoff
         minmat_coord=0.0_dblprec
         minmat_atom=0

         ! Loop over neighbours
         do jneigh=1, nlistsize(iatom)
            jatom=nlist(jneigh,iatom)
            if(jatom>0) then
               jcoord=coord(1:3,jatom)-icoord
               jnorm=norm2(jcoord)

               ! Fix for periodicity if present
               jpcoord=jcoord
               jpnorm=jnorm
               do ix=-perx,perx
                  do iy=-pery,pery
                     do iz=-perz,perz
                        jtcoord(1)=jcoord(1)+ix*N1*C1(1)+iy*N2*C2(1)+iz*N3*C3(1)
                        jtcoord(2)=jcoord(2)+ix*N1*C1(2)+iy*N2*C2(2)+iz*N3*C3(2)
                        jtcoord(3)=jcoord(3)+ix*N1*C1(3)+iy*N2*C2(3)+iz*N3*C3(3)
                        jtnorm=norm2(jtcoord)
                        if(jtnorm<jpnorm) then
                           jpcoord=jtcoord
                           jpnorm=jtnorm
                        end if
                     end do
                  end do
               end do
               jcoord=jpcoord
               jnorm=jpnorm

               ! Loop over all quadrants and directions to find closest neighbours
               do ix=-1,1
                  signx=ix !(-1)**ix
                  do iy=-1,1
                     signy=iy !(-1)**iy
                     do iz=-1,1
                        signz=iz !(-1)**iz
                        if(jnorm<minmat_norm(ix,iy,iz)) then
                           if((signx>0.and.jcoord(1)>0.0_dblprec.or.signx==0.and.abs(jcoord(1))<dbl_tolerance &
                              .or.signx<0.and.jcoord(1)<0.0_dblprec) &
                              .and.(signy>0.and.jcoord(2)>0.0_dblprec.or.signy==0 &
                              .and.abs(jcoord(2))<dbl_tolerance.or.signy<0.and.jcoord(2)<0.0_dblprec) &
                              .and.(signz>0.and.jcoord(3)>0.0_dblprec.or.signz==0 &
                              .and.abs(jcoord(3))<dbl_tolerance.or.signz<0.and.jcoord(3)<0.0_dblprec) &
                              ) then
                              minmat_coord(1:3,ix,iy,iz)=jcoord
                              minmat_norm(ix,iy,iz)=jnorm
                              minmat_atom(ix,iy,iz)=jatom
                           end if
                        end if
                     end do
                  end do
               end do
            end if
         end do

         ! Make neighbour list
         icount=0
         do ix=-1,1
            do iy=-1,1
               do iz=-1,1
                  if(minmat_atom(ix,iy,iz)>0) then
                     icount=icount+1
                     dxyz_vec(:,icount,iatom)=minmat_coord(:,ix,iy,iz)
                     dxyz_atom(icount,iatom)=minmat_atom(ix,iy,iz)
                  end if
               end do
            end do
         end do
         dxyz_list(iatom)=icount
      end do

   end subroutine setup_stencil_mesh


   !> Calculate ((j * d/dr) m ) (which then ends up as one part of the spin transfer torque)
   subroutine differentiate_moments(Natom, Mensemble,emomM, dmomdr, sitenatomjvec)
    
    use MultiscaleGradients
    use MultiscaleInterpolation
    

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom, Mensemble), intent(out) :: dmomdr  !< Divergence of magnetic moment vector
      real(dblprec), dimension(3,natom), intent(in) :: sitenatomjvec !< Spin-current vector

      integer :: iatom, jneigh, jatom, k
      real(dblprec) :: dx, dy, dz, d_mom_x, d_mom_y, d_mom_z

      dmomdr=0.0_dblprec
  
    if(do_multiscale) then       
       call calculateMGradients(emomM,dmomdr)
       call multiscaleInterpolateArray(dmomdr)
    else
       

      !$omp parallel do default(shared) private(iatom,k,jneigh,jatom,d_mom_x,dx,d_mom_y,dy,d_mom_z,dz)
      do iatom=1, Natom
         do k=1, Mensemble
            do jneigh=1, dxyz_list(iatom)
               jatom=dxyz_atom(jneigh,iatom)
               d_mom_x=emomM(1,jatom,k)-emomM(1,iatom,k)
               dx=dxyz_vec(1,jneigh,iatom)
               d_mom_y=emomM(2,jatom,k)-emomM(2,iatom,k)
               dy=dxyz_vec(2,jneigh,iatom)
               d_mom_z=emomM(3,jatom,k)-emomM(3,iatom,k)
               dz=dxyz_vec(3,jneigh,iatom)

               if(abs(dx)>1.0d-7) then
                  dmomdr(1,iatom,k)=dmomdr(1,iatom,k)+(d_mom_x/dx)*sitenatomjvec(1,iatom)
                  dmomdr(2,iatom,k)=dmomdr(2,iatom,k)+(d_mom_y/dx)*sitenatomjvec(1,iatom)
                  dmomdr(3,iatom,k)=dmomdr(3,iatom,k)+(d_mom_z/dx)*sitenatomjvec(1,iatom)
               end if
               if(abs(dy)>1.0d-7) then
                  dmomdr(1,iatom,k)=dmomdr(1,iatom,k)+(d_mom_x/dy)*sitenatomjvec(2,iatom)
                  dmomdr(2,iatom,k)=dmomdr(2,iatom,k)+(d_mom_y/dy)*sitenatomjvec(2,iatom)
                  dmomdr(3,iatom,k)=dmomdr(3,iatom,k)+(d_mom_z/dy)*sitenatomjvec(2,iatom)
               end if
               if(abs(dz)>1.0d-7) then
                  dmomdr(1,iatom,k)=dmomdr(1,iatom,k)+(d_mom_x/dz)*sitenatomjvec(3,iatom)
                  dmomdr(2,iatom,k)=dmomdr(2,iatom,k)+(d_mom_y/dz)*sitenatomjvec(3,iatom)
                  dmomdr(3,iatom,k)=dmomdr(3,iatom,k)+(d_mom_z/dz)*sitenatomjvec(3,iatom)
               end if

            end do

            dmomdr(1,iatom,k)=dmomdr(1,iatom,k)/dxyz_list(iatom)
            dmomdr(2,iatom,k)=dmomdr(2,iatom,k)/dxyz_list(iatom)
            dmomdr(3,iatom,k)=dmomdr(3,iatom,k)/dxyz_list(iatom)

         end do
      end do
      !$omp end parallel do
    endif
   end subroutine differentiate_moments


   

   !> Calculate (dx m(r), dy m(r), and dz( m(r)  )
   subroutine grad_moments(Natom, Mensemble,emomM, grad_mom)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,3,Natom, Mensemble), intent(out) :: grad_mom  !< Gradient of magnetic moment vector


      integer :: iatom, jneigh, jatom, kk
      real(dblprec) :: dx, dy, dz, d_mom_x, d_mom_y, d_mom_z

      grad_mom=0.0_dblprec
      !$omp parallel do default(shared) private(iatom,kk,jneigh,jatom,d_mom_x,dx,d_mom_y,dy,d_mom_z,dz)
      do iatom=1, Natom
         do kk=1, Mensemble
            do jneigh=1, dxyz_list(iatom)
               jatom=dxyz_atom(jneigh,iatom)
               d_mom_x=emomM(1,jatom,kk)-emomM(1,iatom,kk)
               dx=dxyz_vec(1,jneigh,iatom)
               d_mom_y=emomM(2,jatom,kk)-emomM(2,iatom,kk)
               dy=dxyz_vec(2,jneigh,iatom)
               d_mom_z=emomM(3,jatom,kk)-emomM(3,iatom,kk)
               dz=dxyz_vec(3,jneigh,iatom)

               if(abs(dx)>1.0d-7) then
                  grad_mom(1,1,iatom,kk)=grad_mom(1,1,iatom,kk)+d_mom_x/dx
                  grad_mom(2,1,iatom,kk)=grad_mom(2,1,iatom,kk)+d_mom_y/dx
                  grad_mom(3,1,iatom,kk)=grad_mom(3,1,iatom,kk)+d_mom_z/dx
               end if
               if(abs(dy)>1.0d-7) then
                  grad_mom(1,2,iatom,kk)=grad_mom(1,2,iatom,kk)+d_mom_x/dy
                  grad_mom(2,2,iatom,kk)=grad_mom(2,2,iatom,kk)+d_mom_y/dy
                  grad_mom(3,2,iatom,kk)=grad_mom(3,2,iatom,kk)+d_mom_z/dy
               end if
               if(abs(dz)>1.0d-7) then
                  grad_mom(1,3,iatom,kk)=grad_mom(1,3,iatom,kk)+d_mom_x/dz
                  grad_mom(2,3,iatom,kk)=grad_mom(2,3,iatom,kk)+d_mom_y/dz
                  grad_mom(3,3,iatom,kk)=grad_mom(3,3,iatom,kk)+d_mom_z/dz
               end if
            end do
            grad_mom(1:3,1,iatom,kk)=grad_mom(1:3,1,iatom,kk)/dxyz_list(iatom)
            grad_mom(1:3,2,iatom,kk)=grad_mom(1:3,2,iatom,kk)/dxyz_list(iatom)
            grad_mom(1:3,3,iatom,kk)=grad_mom(1:3,3,iatom,kk)/dxyz_list(iatom)
         end do
      end do
      !$omp end parallel do

   end subroutine grad_moments


   !> Calculate (dx m(r), dy m(r), and dz( m(r)  ) for differen atomic species
   subroutine proj_grad_moments(NT,Natom, Mensemble,atype,emomM,proj_grad_mom)

      implicit none

      integer, intent(in) :: NT
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,3,Natom,Mensemble,NT), intent(out) :: proj_grad_mom


      integer :: iatom, jneigh, jatom, k, ii
      real(dblprec) :: dx, dy, dz, d_mom_x, d_mom_y, d_mom_z

      proj_grad_mom=0.0_dblprec

      !$omp parallel do default(shared) private(iatom,k,jneigh,jatom,d_mom_x,dx,d_mom_y,dy,d_mom_z,dz,ii)
      do iatom=1, Natom
         do k=1, Mensemble
            ii=atype(iatom)
            do jneigh=1, dxyz_list(iatom)
               jatom=dxyz_atom(jneigh,iatom)
               d_mom_x=emomM(1,jatom,k)-emomM(1,iatom,k)
               dx=dxyz_vec(1,jneigh,iatom)
               d_mom_y=emomM(2,jatom,k)-emomM(2,iatom,k)
               dy=dxyz_vec(2,jneigh,iatom)
               d_mom_z=emomM(3,jatom,k)-emomM(3,iatom,k)
               dz=dxyz_vec(3,jneigh,iatom)

               if(abs(dx)>1.0d-7) then
                  proj_grad_mom(1,1,iatom,k,ii)=proj_grad_mom(1,1,iatom,k,ii)+d_mom_x/dx
                  proj_grad_mom(2,1,iatom,k,ii)=proj_grad_mom(2,1,iatom,k,ii)+d_mom_y/dx
                  proj_grad_mom(3,1,iatom,k,ii)=proj_grad_mom(3,1,iatom,k,ii)+d_mom_z/dx
               end if
               if(abs(dy)>1.0d-7) then
                  proj_grad_mom(1,2,iatom,k,ii)=proj_grad_mom(1,2,iatom,k,ii)+d_mom_x/dy
                  proj_grad_mom(2,2,iatom,k,ii)=proj_grad_mom(2,2,iatom,k,ii)+d_mom_y/dy
                  proj_grad_mom(3,2,iatom,k,ii)=proj_grad_mom(3,2,iatom,k,ii)+d_mom_z/dy
               end if
               if(abs(dz)>1.0d-7) then
                  proj_grad_mom(1,3,iatom,k,ii)=proj_grad_mom(1,3,iatom,k,ii)+d_mom_x/dz
                  proj_grad_mom(2,3,iatom,k,ii)=proj_grad_mom(2,3,iatom,k,ii)+d_mom_y/dz
                  proj_grad_mom(3,3,iatom,k,ii)=proj_grad_mom(3,3,iatom,k,ii)+d_mom_z/dz
               end if
            end do
            proj_grad_mom(1:3,1,iatom,k,ii)=proj_grad_mom(1:3,1,iatom,k,ii)/dxyz_list(iatom)
            proj_grad_mom(1:3,2,iatom,k,ii)=proj_grad_mom(1:3,2,iatom,k,ii)/dxyz_list(iatom)
            proj_grad_mom(1:3,3,iatom,k,ii)=proj_grad_mom(1:3,3,iatom,k,ii)/dxyz_list(iatom)

         end do
      end do
      !$omp end parallel do

   end subroutine proj_grad_moments

   subroutine deallocate_gradient_lists()

      implicit none

      integer :: i_all, i_stat

      if (allocated(dxyz_list)) then
         i_all=-product(shape(dxyz_list))*kind(dxyz_list)
         deallocate(dxyz_list,stat=i_stat)
         call memocc(i_stat,i_all,'dxyz_list','deallocate_gradient_lists')
      endif
      if (allocated(dxyz_atom)) then
         i_all=-product(shape(dxyz_atom))*kind(dxyz_atom)
         deallocate(dxyz_atom,stat=i_stat)
         call memocc(i_stat,i_all,'dxyz_atom','deallocate_gradient_lists')
      endif
      if (allocated(dxyz_vec)) then
         i_all=-product(shape(dxyz_vec))*kind(dxyz_vec)
         deallocate(dxyz_vec,stat=i_stat)
         call memocc(i_stat,i_all,'dxyz_vec','deallocate_gradient_lists')
      endif

   end subroutine deallocate_gradient_lists

end module Gradients
