program test_mesh2d
   use Parameters, only : dblprec
   use Mesh2D
   implicit none

   integer, parameter :: n=4, na=1, nlayer=2, natom=n*n, natom_layers=n*n*nlayer
   real(dblprec), allocatable :: coord(:,:), coord_before(:,:)
   real(dblprec), allocatable :: coord_layers(:,:)
   real(dblprec) :: c1(3), c2(3), c3(3), total_area
   complex(dblprec) :: psi(natom), psi_tri(3), gx, gy
   complex(dblprec), parameter :: alpha=(0.3_dblprec,-0.2_dblprec)
   complex(dblprec), parameter :: beta=(1.25_dblprec,0.75_dblprec)
   complex(dblprec), parameter :: gamma=(-0.4_dblprec,0.6_dblprec)
   integer :: ix, iy, iz, isite, ip, itri, iv, index, index_layer
   real(dblprec) :: p(3,2), delta(2)
   complex(dblprec) :: gx_tri, gy_tri
   real(dblprec), parameter :: tol=1.0e-12_dblprec

   allocate(coord(3,natom),coord_before(3,natom))
   c1=(/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
   c2=(/0.0_dblprec,1.0_dblprec,0.0_dblprec/)
   c3=(/0.0_dblprec,0.0_dblprec,1.0_dblprec/)

   do iy=0,n-1
      do ix=0,n-1
         index=iy*n+ix+1
         coord(:,index)=(/real(ix,dblprec),real(iy,dblprec),0.0_dblprec/)
      end do
   end do
   coord_before=coord

   ! The periodic mesh must tile the cell without changing SystemData::coord.
   call mesh2d_build(natom,n,n,1,na,coord,c1,c2,c3,'P','P','0')
   if (nsimp/=2*n*n) error stop 'mesh2d periodic: wrong triangle count'
   total_area=sum(tri_area)
   if (abs(total_area-real(n*n,dblprec))>tol) error stop 'mesh2d periodic: wrong area'
   if (maxval(abs(coord-coord_before))>0.0_dblprec) error stop 'mesh2d changed coordinates'
   do itri=1,nsimp
      if (tri_area(itri)<=0.0_dblprec) error stop 'mesh2d periodic: non-positive area'

      ! Verify the linear FEM coefficients in each minimum-image triangle.
      p=0.0_dblprec
      do iv=1,3
         delta=coord(1:2,simp(iv,itri))-coord(1:2,simp(1,itri))
         if (delta(1)>0.5_dblprec*real(n,dblprec)) delta(1)=delta(1)-real(n,dblprec)
         if (delta(1)<-0.5_dblprec*real(n,dblprec)) delta(1)=delta(1)+real(n,dblprec)
         if (delta(2)>0.5_dblprec*real(n,dblprec)) delta(2)=delta(2)-real(n,dblprec)
         if (delta(2)<-0.5_dblprec*real(n,dblprec)) delta(2)=delta(2)+real(n,dblprec)
         p(iv,:) = delta
      end do
      do iv=1,3
         psi_tri(iv)=alpha+beta*cmplx(p(iv,1),0.0_dblprec,dblprec)+ &
            gamma*cmplx(p(iv,2),0.0_dblprec,dblprec)
      end do
      gx_tri=sum(grad_b(:,itri)*psi_tri)
      gy_tri=sum(grad_c(:,itri)*psi_tri)
      if (abs(gx_tri-beta)>tol .or. abs(gy_tri-gamma)>tol) then
         error stop 'mesh2d periodic: FEM coefficient regression failed'
      end if
   end do

   ! On an open mesh, an affine field is globally single-valued.  Gather the
   ! triangle gradients through the CSR adjacency and check every supported site.
   call mesh2d_build(natom,n,n,1,na,coord,c1,c2,c3,'0','0','0')
   do isite=1,natom
      psi(isite)=alpha+beta*real(coord(1,isite),dblprec)+gamma*real(coord(2,isite),dblprec)
   end do
   do isite=1,natom
      gx=(0.0_dblprec,0.0_dblprec)
      gy=(0.0_dblprec,0.0_dblprec)
      do ip=site_tri_ptr(isite),site_tri_ptr(isite+1)-1
         itri=site_tri_idx(ip)
         gx=gx+tri_area(itri)*sum(grad_b(:,itri)*psi(simp(:,itri)))
         gy=gy+tri_area(itri)*sum(grad_c(:,itri)*psi(simp(:,itri)))
      end do
      if (site_wsum(isite)>0.0_dblprec) then
         gx=gx/site_wsum(isite)
         gy=gy/site_wsum(isite)
         if (abs(gx-beta)>tol .or. abs(gy-gamma)>tol) then
            error stop 'mesh2d open: site FEM gradient regression failed'
         end if
      end if
   end do

   ! Each z layer must receive an independent xy mesh.  In particular, no
   ! triangle may connect the two layers through the layer offset.
   allocate(coord_layers(3,natom_layers))
   do iz=0,nlayer-1
      do iy=0,n-1
         do ix=0,n-1
            index_layer=iz*n*n+iy*n+ix+1
            coord_layers(:,index_layer)=(/real(ix,dblprec),real(iy,dblprec),real(iz,dblprec)/)
         end do
      end do
   end do
   call mesh2d_build(natom_layers,n,n,nlayer,na,coord_layers,c1,c2,c3,'P','P','0')
   if (nsimp/=2*n*n*nlayer) error stop 'mesh2d multilayer: wrong triangle count'
   total_area=sum(tri_area)
   if (abs(total_area-real(n*n*nlayer,dblprec))>tol) error stop 'mesh2d multilayer: wrong area'
   do itri=1,nsimp
      if (tri_area(itri)<=0.0_dblprec) error stop 'mesh2d multilayer: non-positive area'
      if (any((simp(:,itri)-1)/(n*n) /= (simp(1,itri)-1)/(n*n))) then
         error stop 'mesh2d multilayer: triangle crosses z layers'
      end if
   end do

   call mesh2d_release()
   deallocate(coord_layers)
end program test_mesh2d
