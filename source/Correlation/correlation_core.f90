!------------------------------------------------------------------------------------
!> @brief Data and routines for calculate connected spin correlation function \f$ \mathbf{S}\left(\mathbf{r},t\right)\f$
!> and Fourier transforms \f$\mathbf{S}\left(\mathbf{q},\omega\right)\f$
!> @details
!> In this routine magnon disperion relations, magnon density of states, etc. are calculated from the
!> time and spatially displaced correlation function \f$ C^k \left(\mathbf{r}-\mathbf{r'},t\right)\f$ which is defined as
!> \f$ C^k (\mathbf{r}-\mathbf{r'},t) = \langle m^k_{\mathbf{r}}(t) m^k_{\mathbf{r'}}(0) \rangle - \langle m^k_{\mathbf{r}}(t) \rangle \langle m^k_{\mathbf{r'}}(0) \rangle\f$.
!> Using this one can then calculate the dynamical structure factor \f$ \mathbf{S}\left(\mathbf{q},\omega\right)\f$, via a Fourier transforms
!> of the correlation function
!> \f$ S^k(\mathbf{q},\omega) = \frac{1}{\sqrt{2\pi}N} \sum_{\mathbf{r},\mathbf{r'}} e^{i\mathbf{q}\cdot(\mathbf{r}-\mathbf{r'})} \int_{-\infty}^{\infty} e^{i\omega t} C^k (\mathbf{r}-\mathbf{r'},t) dt\f$
!> @author
!> A. Bergman, L. Bergqvist, J. Hellsvik, J. Chico
!> @todo Automatic generation of q-points for do_sc="Y"
!------------------------------------------------------------------------------------
module Correlation_core
   use Parameters
   use Profiling
   use Correlation_type
   use Correlation_utils
   use Correlation_print
   use Correlation_kernels
   use Qvectors
   use Omegas
   !
   implicit none
   !
contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_gk2
   !> @brief Calculate \f$\mathbf{S}\left(\mathbf{q},t\right)\f$ for obtaining \f$\mathbf{S}\left(\mathbf{q},\omega\right)\f$ after FT
   !---------------------------------------------------------------------------------
   subroutine calc_gk2(Natom, Mensemble, NT,atype,Nchmax,achtype, cc, coord, simid, SA, flag)
      !
      use Constants
      use Math_functions, only : gramms
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NT           !< Number of types of atoms
      integer, dimension(Natom), intent(in) :: atype        !< Type of atom
      integer, intent(in) :: Nchmax       !< Number of chemical types
      integer, dimension(Natom), intent(in) :: achtype      !< Chemistry of atom

      type(corr_t), intent(inout) :: cc !< Derived type for correlation data
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: SA     !< First vector to correlate
      integer, intent(inout) :: flag  !< Setup, sample, or print
      !
      integer ::  l, i_stat, i_all
      complex(dblprec) :: i, iqfac
      real(dblprec) :: nainv, mavg_norm, win_fac
      real(dblprec), dimension(3) :: SA_avrg
      real(dblprec) :: qfac
      !
      !

      i=(0.0_dblprec,1.0_dblprec)

      if(cc%gk_flag==0) then

         ! First call, allocate and clear arrays
         allocate(cc%m_k(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(cc%m_k))*kind(cc%m_k),'m_k','calc_gk2')
         cc%m_k=0.0_dblprec

         if(cc%do_proj=='C'.or.cc%do_proj=='Y') then
            allocate(cc%m_k_proj(3,nt,nq),stat=i_stat)
            call memocc(i_stat,product(shape(cc%m_k_proj))*kind(cc%m_k_proj),'m_k_proj','calc_gk2')
            cc%m_k_proj=0.0_dblprec
         end if

         if(cc%do_projch=='C'.or.cc%do_projch=='Y') then
            allocate(cc%m_k_projch(3,Nchmax,nq),stat=i_stat)
            call memocc(i_stat,product(shape(cc%m_k_projch))*kind(cc%m_k_projch),'m_k_projch','calc_gk2')
            cc%m_k_projch=0.0_dblprec
         end if

         if (cc%do_sc_local_axis=='Y'.or.cc%do_sc_local_axis=='B') then

            call calc_mavrg_vec(Natom,Mensemble,SA,SA_avrg,cc%SA_axis)
            cc%mavg_local_axis=SA

            call find_local_rotmat(Natom*Mensemble,SA,cc%mavg_local_rotmat)

         else

            call calc_mavrg_vec(Natom,Mensemble,SA,SA_avrg,cc%SA_axis)

            do l=1,Mensemble
               !SA_axis(:,l)=SA_axis(:,l)/Natom
               mavg_norm=sum(cc%SA_axis(:,l)*cc%SA_axis(:,l))**0.5_dblprec
               if(mavg_norm>1.0d-2) then
                  cc%SA_axis(:,l)=cc%SA_axis(:,l)/mavg_norm
               else
                  cc%SA_axis(:,l)=(/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
               end if

            end do

            !call gramms(cc%SA_axis,cc%mort_axis,Mensemble)

         end if
         !
         call find_rmid(r_mid,coord,Natom)

         cc%gk_flag=1
         cc%sc_nsamp=0
         flag=1

      end if

      nainv=1.0_dblprec/Natom
      qfac=2.0_dblprec*pi

      ! Calculate g(k) for the current iteration and add to G(k,t)
      if (cc%gk_flag==1.and.flag<2) then

         iqfac=i*qfac

         ! Project measurable to local axis if supposed to
         call project_local_axis(Natom,Mensemble,cc,SA)

         win_fac=1.0_dblprec*nainv
         call corr_kernel_single(Natom,Mensemble,nq,coord,r_mid,iqfac,win_fac,cc%m_loc,cc%m_k(1:3,1:nq))


         if(cc%do_proj=='C'.or.cc%do_proj=='Y') then
            call corr_kernel_proj(Natom,NT,Mensemble,nq,coord,r_mid,atype,iqfac,win_fac,cc%m_loc,cc%m_k_proj(:,:,1:nq))
         end if

         if(cc%do_projch=='C'.or.cc%do_projch=='Y') then
            call corr_kernel_proj(Natom,Nchmax,Mensemble,nq,coord,r_mid,achtype,iqfac,win_fac,cc%m_loc,cc%m_k_projch(:,:,1:nq))
         end if

         cc%sc_nsamp=cc%sc_nsamp +  1

      end if

      ! Final operations, transform and print
      if (flag==2) then

         if (flag==2) cc%gk_flag=2

         call print_gk(NT, Nchmax, cc, cc, simid, cc%label)

         call print_gr(Natom, cc,  coord, simid)

         i_all=-product(shape(cc%m_k))*kind(cc%m_k)
         deallocate(cc%m_k,stat=i_stat)
         call memocc(i_stat,i_all,'m_k','calc_gk2')

         if(cc%do_proj=='C'.or.cc%do_proj=='Y') then
            i_all=-product(shape(cc%m_k_proj))*kind(cc%m_k_proj)
            deallocate(cc%m_k_proj,stat=i_stat)
            call memocc(i_stat,i_all,'m_k_proj','calc_gk2')
         end if

         if(cc%do_projch=='C'.or.cc%do_projch=='Y') then
            i_all=-product(shape(cc%m_k_projch))*kind(cc%m_k_projch)
            deallocate(cc%m_k_projch,stat=i_stat)
            call memocc(i_stat,i_all,'m_k_projch','calc_gk2')
         end if

      end if

      if (flag==3) then
         !call deallocate_gk(cc)
      end if

      return
      !
      !
   end subroutine calc_gk2

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_gkt
   !> @brief Calculate \f$\mathbf{S}\left(\mathbf{q},t\right)\f$ for obtaining \f$\mathbf{S}\left(\mathbf{q},\omega\right)\f$ after FT
   !---------------------------------------------------------------------------------
   subroutine calc_gkt(Natom, Mensemble, NT,atype,Nchmax,achtype, cc, coord, SA, flag)
      !
      use Constants
      use Math_functions, only : gramms
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NT           !< Number of types of atoms
      integer, dimension(Natom), intent(in) :: atype        !< Type of atom
      integer, intent(in) :: Nchmax       !< Number of chemical types
      integer, dimension(Natom), intent(in) :: achtype      !< Chemistry of atom

      type(corr_t), intent(inout) :: cc !< Derived type for correlation data
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: SA     !< First vector to correlate
      integer, intent(inout) :: flag  !< Setup, sample, or print
      !
      integer ::  l, i_stat
      complex(dblprec) :: i, iqfac
      real(dblprec) :: nainv, mavg_norm, win_fac
      real(dblprec), dimension(3) :: SA_avrg
      real(dblprec) :: qfac
      !
      !

      i=(0.0_dblprec,1.0_dblprec)

      if(cc%gkt_flag==0) then

         ! First call, allocate and clear arrays
         allocate(cc%m_kt(3,nq,cc%sc_max_nstep),stat=i_stat)
         call memocc(i_stat,product(shape(cc%m_kt))*kind(cc%m_kt),'cc%m_kt','calc_gkt')
         cc%m_kt=0.0_dblprec

         if(cc%do_proj=='Q'.or.cc%do_proj=='T'.or.cc%do_proj=='Y') then
            allocate(cc%m_kt_proj(3,nt,nq,cc%sc_max_nstep),stat=i_stat)
            call memocc(i_stat,product(shape(cc%m_kt_proj))*kind(cc%m_kt_proj),'cc%m_kt_proj','calc_gkt')
            cc%m_kt_proj=0.0_dblprec
         end if

         if(cc%do_projch=='Q'.or.cc%do_projch=='T'.or.cc%do_projch=='Y') then
            allocate(cc%m_kt_projch(3,Nchmax,nq,cc%sc_max_nstep),stat=i_stat)
            call memocc(i_stat,product(shape(cc%m_kt_projch))*kind(cc%m_kt_projch),'cc%m_kt_projch','calc_gkt')
            cc%m_kt_projch=0.0_dblprec
         end if


         !call allocate_corr(Natom,Mensemble, cc)

         if (cc%do_sc_local_axis=='Y'.or.cc%do_sc_local_axis=='B') then

            call calc_mavrg_vec(Natom,Mensemble,SA,SA_avrg,cc%SA_axis)
            cc%mavg_local_axis=SA

            call find_local_rotmat(Natom*Mensemble,SA,cc%mavg_local_rotmat)

         else

            call calc_mavrg_vec(Natom,Mensemble,SA,SA_avrg,cc%SA_axis)

            do l=1,Mensemble
               !SA_axis(:,l)=SA_axis(:,l)/Natom
               mavg_norm=sum(cc%SA_axis(:,l)*cc%SA_axis(:,l))**0.5_dblprec
               if(mavg_norm>1.0d-2) then
                  cc%SA_axis(:,l)=cc%SA_axis(:,l)/mavg_norm
               else
                  cc%SA_axis(:,l)=(/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
               end if

            end do


         end if
         !

         cc%gkt_flag=1
         cc%sc_samp_done=0
         flag=1

      end if

      nainv=1.0_dblprec/Natom
      qfac=2.0_dblprec*pi

      ! Calculate g(k) for the current iteration and add to G(k,t)
      if (cc%gkt_flag==1.and.flag<2) then

         iqfac=i*qfac

         ! Project measurable to local axis if supposed to
         call project_local_axis(Natom,Mensemble,cc,SA)

         win_fac=1.0_dblprec*nainv
         call corr_kernel_single(Natom,Mensemble,nq,coord,r_mid,iqfac,win_fac,cc%m_loc,cc%m_kt(1:3,1:nq,cc%sc_tidx))

         if(cc%do_proj=='Q'.or.cc%do_proj=='T'.or.cc%do_proj=='Y') then
            call corr_kernel_proj(Natom,NT,Mensemble,nq,coord,r_mid,atype,iqfac,win_fac,cc%m_loc,cc%m_kt_proj(:,:,1:nq,cc%sc_tidx))
         end if

         if(cc%do_projch=='Q'.or.cc%do_projch=='T'.or.cc%do_projch=='Y') then
            call corr_kernel_proj(Natom,Nchmax,Mensemble,nq,coord,r_mid,achtype,iqfac,win_fac, &
                     cc%m_loc,cc%m_kt_projch(:,:,1:nq,cc%sc_tidx))
         end if

      end if

      ! Final operations, transform and print
      if (flag==2) then

         if (flag==2) cc%gkt_flag=2

         call calc_gkw(NT, Nchmax, cc)

      end if

      !!!if (flag==3) then

      !!!   call print_gkw(Natom, Mensemble,NT,atype,Nchmax,achtype, cc, coord, simid, cc%label)

      !!!end if

      !!!if (flag==4) then

      !!!   call print_gkt(Natom, Mensemble,NT,atype,Nchmax,achtype, cc, coord, simid, cc%label)

      !!!end if

      !!!if (flag==-1) then

      !!!   call deallocate_gkw(cc)

      !!!end if

      return
      !
      !
   end subroutine calc_gkt

subroutine calc_gkw(NT, Nchmax, cc)

      use Constants
      use Math_functions, only : gramms
      !
      implicit none
      !
      integer, intent(in) :: NT           !< Number of types of atoms
      integer, intent(in) :: Nchmax       !< Number of chemical types

      type(corr_t), intent(inout) :: cc !< Derived type for correlation data

      !
      integer  :: i_stat, j
      complex(dblprec)  :: i
      !
      !

      i=(0.0_dblprec,1.0_dblprec)


         !print *,'-------------------------------------------------------------', label


         allocate(cc%m_kw(3,nq,cc%nw),stat=i_stat)
         call memocc(i_stat,product(shape(cc%m_kw))*kind(cc%m_kw),'cc%m_kw','calc_gkw')
         cc%m_kw=0.0_dblprec

         allocate(cc%dt(cc%sc_max_nstep),stat=i_stat)
         call memocc(i_stat,product(shape(cc%dt))*kind(cc%dt),'cc%dt','calc_gkw')

         allocate(cc%time(cc%sc_max_nstep),stat=i_stat)
         call memocc(i_stat,product(shape(cc%time))*kind(cc%time),'cc%time','calc_gkw')

         ! Fill time arrays
         do j=1,cc%sc_max_nstep
            cc%dt(j) = cc%scstep_arr(j)*cc%deltat_corr(j)
         end do
         cc%time(1)=0.0_dblprec
         do j=2,cc%sc_max_nstep
            cc%time(j)=cc%time(j-1)+cc%dt(j)
         end do

         call calc_sqt_conv(3,cc,cc%dt,cc%m_kt)

         cc%m_kw=0.0_dblprec
         call corr_kernel_time(nq,cc%nw,3,cc,cc%dt,sc_window_fun,cc%m_kt,cc%m_kw)

         ! Calculate the convolution for sqw
         call calc_sqw_conv(3,cc%nw,cc%m_kw)


         if(cc%do_proj=='Q'.or.cc%do_proj=='Y') then

            ! Allocate arrays
            allocate(cc%m_kw_proj(3,nt,nq,cc%nw),stat=i_stat)
            call memocc(i_stat,product(shape(cc%m_kw_proj))*kind(cc%m_kw_proj),'cc%m_kw_proj','calc_gkw')
            cc%m_kw_proj=0.0_dblprec


            call calc_sqt_conv(3*nt,cc,cc%dt,cc%m_kt_proj)

            call corr_kernel_time(nq,cc%nw,3*nt,cc,cc%dt,sc_window_fun,cc%m_kt_proj,cc%m_kw_proj)

            ! Calculate the convolution for sqw
            call calc_sqw_conv(3*nt,cc%nw,cc%m_kw_proj)


         end if

         if(cc%do_projch=='Q'.or.cc%do_projch=='Y') then

            ! Allocate arrays
            allocate(cc%m_kw_projch(3,Nchmax,nq,cc%nw),stat=i_stat)
            call memocc(i_stat,product(shape(cc%m_kw_projch))*kind(cc%m_kw_projch),'cc%m_kw_projch','calc_gkw')
            cc%m_kw_projch=0.0_dblprec

            call calc_sqt_conv(3*Nchmax,cc,cc%dt,cc%m_kt_projch)

            call corr_kernel_time(nq,cc%nw,3*Nchmax,cc,cc%dt,sc_window_fun,cc%m_kt_projch,cc%m_kw_projch)

            ! Calculate the convolution for sqw
            call calc_sqw_conv(3*Nchmax*Nchmax,cc%nw,cc%m_kw_projch)

         end if

         !
      !
   end subroutine calc_gkw

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_sr
   !> @brief Perform only spatial correlation in real space to be able to deal with non periodic systems
   !---------------------------------------------------------------------------------
   subroutine calc_sr(Natom, Mensemble, cc, simid, emomM, cr_flag)
      !

      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      type(corr_t) :: cc !< Correlation struct
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(inout) :: cr_flag  !< Allocate or deallocate (1/-1)
      !
      integer :: iatom,r,l,i_stat,i_all, c_idx
      character(len=30) :: filn
      real(dblprec) :: nainv

      if(cr_flag==0) then
         ! First call, allocate and clear arrays
         
         !!! Create neighbour list for truncated S(r) correlations
         !!call f_get_nnlist(Natom,coord,cc%cr_list,cc%cr_vect,cc%cr_nnmax,1,cc%cr_cut)
         !!! Identify identical r-vectors and create look-up table
         !!call f_get_nnhist(Natom,cc%cr_list,cc%cr_vect,cc%cr_nnmax,cc%cr_nhist,cc%cr_uniq_vect,cc%cr_lut,1)
        
         allocate(cc%corr_sr(3,cc%cr_nhist),stat=i_stat)
         call memocc(i_stat,product(shape(cc%corr_sr))*kind(cc%corr_sr),'cc%corr_sr','calc_sr')

         allocate(cc%cr_hist(cc%cr_nhist),stat=i_stat)
         call memocc(i_stat,product(shape(cc%cr_hist))*kind(cc%cr_hist),'cc%cr_hist','calc_sr')

         cc%corr_sr=0.0_dblprec
         cr_flag=1
         cc%sc_samp_done_sr=0
      
      end if

      nainv=1.0_dblprec/Natom

      if (cr_flag==1) then
         ! Calculate s(k) for the current iteration and add to average of G(k)

         !!! if (cc%do_connected=='Y') then
         !!!    !$omp parallel do default(shared) private(iatom,l) schedule(static)
         !!!    do l=1,Mensemble
         !!!       do iatom=1,Natom
         !!!          connected(1,iatom)=connected(1,iatom)+emomM(1,iatom,l)
         !!!          connected(2,iatom)=connected(2,iatom)+emomM(2,iatom,l)
         !!!          connected(3,iatom)=connected(3,iatom)+emomM(3,iatom,l)
         !!!       enddo
         !!!    enddo
         !!!    !$omp end parallel do
         !!! else
         !!!    connected(:,:)=0.0_dblprec
         !!! endif

         !!!$omp parallel do default(shared) private(r,iatom,l,c_idx) schedule(static)
         do l=1,Mensemble
            do iatom=1,Natom
               !do r=1,Natom
               do r=1,cc%cr_list(1,iatom)
                  c_idx = cc%cr_lut(r,iatom)
                  cc%corr_sr(:,c_idx) = cc%corr_sr(:,c_idx) + emomM(:,iatom,l) * emomM(:,cc%cr_list(r+1,iatom),l)  
                  !cc%corr_sr(1,iatom)=cc%corr_sr(1,iatom)+emomM(1,iatom,l)*emomM(1,r,l)*nainv-connected(1,iatom)*connected(1,r)/Mensemble
                  !cc%corr_sr(2,iatom)=cc%corr_sr(2,iatom)+emomM(2,iatom,l)*emomM(2,r,l)*nainv-connected(2,iatom)*connected(2,r)/Mensemble
                  !cc%corr_sr(3,iatom)=cc%corr_sr(3,iatom)+emomM(3,iatom,l)*emomM(3,r,l)*nainv-connected(3,iatom)*connected(3,r)/Mensemble
               end do
            enddo
         end do
         !!!$omp end parallel do

         cc%sc_samp_done_sr=cc%sc_samp_done_sr+1
      end if

      if (cr_flag==2) then
         ! Finish sampling and write S(q)
         !cc%corr_sr=cc%corr_sr/(cc%sc_samp_done_sr)
         cc%corr_sr=cc%corr_sr/(cc%sc_samp_done_sr*Mensemble*Natom)

         ! Write G(r)
         write (filn,'(''dir_sr.'',a,''.out'')') trim(simid)
         open(ofileno,file=filn,status='replace')
         do r=1,cc%cr_nhist
            write(ofileno,'(i10,3f10.4,5f18.8)') r,cc%cr_uniq_vect(:,r),cc%corr_sr(:,r),&
               norm2(cc%corr_sr(:,r)),sum(cc%corr_sr(:,r))
         end do
         !do r=1,Natom
         !   write(ofileno,'(i10,3f10.4,5f18.8)') r,(coord(l,r),l=1,3),(((cc%corr_sr(l,r))),l=1,3),&
         !      sqrt(cc%corr_sr(1,r)**2+cc%corr_sr(2,r)**2+cc%corr_sr(3,r)**2),cc%corr_sr(1,r)+cc%corr_sr(2,r)+cc%corr_sr(3,r)
         !end do
         close(ofileno)

         ! Write G(|r|)
         write (filn,'(''dir_sra.'',a,''.out'')') trim(simid)
         open(ofileno,file=filn,status='replace')
         do r=1,cc%cr_nhist
            write(ofileno,'(7f18.8)') norm2(cc%cr_uniq_vect(:,r)),cc%corr_sr(:,r),norm2(cc%corr_sr(:,r))
         end do
         !do r=1,Natom
         !   write(ofileno,'(7f18.8)') sqrt((coord(1,r)-coord(1,1))**2+(coord(2,r)-coord(2,1))**2+(coord(3,r)-coord(3,1))**2),&
         !      (((cc%corr_sr(l,r))),l=1,3),&
         !      sqrt(cc%corr_sr(1,r)**2+cc%corr_sr(2,r)**2+cc%corr_sr(3,r)**2),cc%corr_sr(1,r)+cc%corr_sr(2,r)+cc%corr_sr(3,r)
         !end do
         close(ofileno)

         ! Deallocate arrays
         i_all=-product(shape(cc%corr_sr))*kind(cc%corr_sr)
         deallocate(cc%corr_sr,stat=i_stat)
         call memocc(i_stat,i_all,'cc%corr_sr','calc_sr')

         i_all=-product(shape(cc%cr_hist))*kind(cc%cr_hist)
         deallocate(cc%cr_hist,stat=i_stat)
         call memocc(i_stat,i_all,'cc%cr_hist','calc_sr')
      end if
      return

   end subroutine calc_sr



   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_srt
   !> @brief Perform direct real-space correlation in real space and time
   !---------------------------------------------------------------------------------
   subroutine calc_srt(Natom, Mensemble, cc, coord, simid, emomM, cr_flag)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      type(corr_t) :: cc !< Correlation struct
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(inout) :: cr_flag  !< Allocate or deallocate (1/-1)
      !
      integer :: iatom,r,l,i_stat,i_all, jatom
      character(len=30) :: filn
      real(dblprec), dimension(3,Natom) :: connected
      real(dblprec) :: nainv
      real(dblprec) :: r_max

      if(cr_flag==0) then
         ! First call, allocate and clear arrays

         ! Start by making look-up tables for correlation measurements

         r_max = 3.0_dblprec

!        hit = 0 

         do iatom=1,Natom
!           dist_arr=0.0_dblprec
            do jatom=1,Natom
!              dist_arr(jatom) = norm2(coord(:,iatom) - coord(:,jatom))
            end do
         end do

         allocate(cc%corr_srt(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(cc%corr_sr))*kind(cc%corr_sr),'corr_sr','calc_sr')
         cc%corr_sr=0.0_dblprec
         cr_flag=1
         cc%sc_samp_done_sr=0
      end if

      nainv=1.0_dblprec/Natom

      if (cr_flag==1) then
         ! Calculate s(k) for the current iteration and add to average of G(k)

         if (cc%do_connected=='Y') then
            !$omp parallel do default(shared) private(iatom,l) schedule(static)
            do l=1,Mensemble
               do iatom=1,Natom
                  connected(1,iatom)=connected(1,iatom)+emomM(1,iatom,l)
                  connected(2,iatom)=connected(2,iatom)+emomM(2,iatom,l)
                  connected(3,iatom)=connected(3,iatom)+emomM(3,iatom,l)
               enddo
            enddo
            !$omp end parallel do
         else
            connected(:,:)=0.0_dblprec
         endif

         !$omp parallel do default(shared) private(r,iatom,l) schedule(static)
         do l=1,Mensemble
            do iatom=1,Natom
               do r=1,Natom
                  cc%corr_sr(1,iatom)= &
                     cc%corr_sr(1,iatom)+emomM(1,iatom,l)*emomM(1,r,l)*nainv-connected(1,iatom)*connected(1,r)/Mensemble
                  cc%corr_sr(2,iatom)= &
                     cc%corr_sr(2,iatom)+emomM(2,iatom,l)*emomM(2,r,l)*nainv-connected(2,iatom)*connected(2,r)/Mensemble
                  cc%corr_sr(3,iatom)= &
                     cc%corr_sr(3,iatom)+emomM(3,iatom,l)*emomM(3,r,l)*nainv-connected(3,iatom)*connected(3,r)/Mensemble
               end do
            enddo
         end do
         !$omp end parallel do

         cc%sc_samp_done_sr=cc%sc_samp_done_sr+1
      end if

      if (cr_flag==2) then
         ! Finish sampling and write S(q)
         cc%corr_sr=cc%corr_sr/(cc%sc_samp_done_sr*Mensemble)

         ! Write G(r)
         write (filn,'(''dir_sr.'',a,''.out'')') trim(simid)
         open(ofileno,file=filn,status='replace')
         do r=1,Natom
            write(ofileno,'(i10,3f10.4,5f18.8)') r,(coord(l,r),l=1,3),(((cc%corr_sr(l,r))),l=1,3),&
               sqrt(cc%corr_sr(1,r)**2+cc%corr_sr(2,r)**2+cc%corr_sr(3,r)**2),cc%corr_sr(1,r)+cc%corr_sr(2,r)+cc%corr_sr(3,r)
         end do
         close(ofileno)

         ! Write G(|r|)
         write (filn,'(''dir_sra.'',a,''.out'')') trim(simid)
         open(ofileno,file=filn,status='replace')
         do r=1,Natom
            write(ofileno,'(7f18.8)') sqrt((coord(1,r)-coord(1,1))**2+(coord(2,r)-coord(2,1))**2+(coord(3,r)-coord(3,1))**2),&
               (((cc%corr_sr(l,r))),l=1,3),&
               sqrt(cc%corr_sr(1,r)**2+cc%corr_sr(2,r)**2+cc%corr_sr(3,r)**2),cc%corr_sr(1,r)+cc%corr_sr(2,r)+cc%corr_sr(3,r)
         end do
         close(ofileno)

         ! Deallocate arrays
         i_all=-product(shape(cc%corr_sr))*kind(cc%corr_sr)
         deallocate(cc%corr_sr,stat=i_stat)
         call memocc(i_stat,i_all,'cc%corr_sr','calc_sr')
      end if
      return

   end subroutine calc_srt


end module Correlation_core
