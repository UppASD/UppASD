 module Correlation_Print
   use Parameters
   use Profiling
   use Correlation_type
   use Correlation_utils
   use Correlation_kernels
   use Qvectors
   use Omegas
   !
   implicit none

   private

   public :: print_gkw, print_gkt, print_gk, print_gr

contains

   subroutine print_gkw(NT,Nchmax, cc, dc, simid, label)

      use Constants
      !
      implicit none
      !
      integer, intent(in) :: NT           !< Number of types of atoms
      integer, intent(in) :: Nchmax       !< Number of chemical types

      type(corr_t), intent(inout) :: cc !< Derived type for correlation data
      type(corr_t), intent(inout) :: dc !< Derived type for correlation data

      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=2), intent(in) :: label

      !
      integer  :: iw, iq, i_stat, i_all
      integer :: ia,ib, it, jt
      character(len=30) :: filn
      complex(dblprec)  :: i
      real(dblprec), dimension(3,3) :: unit3
      real(dblprec) :: qnorm, qnorm2, polfac
      !
      !
      complex(dblprec), dimension(:,:,:),     allocatable :: c_kw                ! Correlation for G(k,w)
      complex(dblprec), dimension(:,:,:,:,:), allocatable :: c_kw_proj           ! Correlation for G(k,w)
      complex(dblprec), dimension(:,:,:,:,:), allocatable :: c_kw_projch         ! Correlation for G(k,w)
      complex(dblprec), dimension(:,:,:,:),   allocatable :: c_kw_tens           ! Correlation for G(k,w)

      complex(dblprec), dimension(:,:),       allocatable :: sqwintensity        ! Intensity for G(k,w)

      i=(0.0_dblprec,1.0_dblprec)


      ! Allocate arrays

      allocate(c_kw(3,nq,cc%nw),stat=i_stat)
      call memocc(i_stat,product(shape(c_kw))*kind(c_kw),'c_kw','print_gkw')
      c_kw=0.0_dblprec

      call combine_corr_scalar(nq, 3, cc%nw, cc%m_kw, dc%m_kw, c_kw)

      c_kw = c_kw / (0.5_dblprec*(cc%sc_tidx-1))  / (0.5_dblprec*(dc%sc_tidx-1))

      if(cc%do_proj=='Q'.or.cc%do_proj=='Y') then

         ! Allocate arrays
         allocate(c_kw_proj(3,nt,nt,nq,cc%nw),stat=i_stat)
         call memocc(i_stat,product(shape(c_kw_proj))*kind(c_kw_proj),'c_kw_proj','print_gkw')
         c_kw_proj=0.0_dblprec

         call combine_corr_proj_scalar(nt, nq, 3, cc%nw, cc%m_kw_proj, cc%m_kw_proj, c_kw_proj)

         ! Write S(q,w)
         do it=1,nt
            do jt=1,nt
               write (filn,'(a,''qw_proj.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
               call corr_write_abscorr(nq,cc%nw,3,filn,c_kw_proj(:,it,jt,:,:),cc%w)
               if(cc%do_sc_list=='Y') then
                  write (filn,'(a,''qw_proj_l.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                  call corr_write_abscorr_list(nq,cc%nw,3,filn,c_kw_proj(:,it,jt,:,:),cc%w)
               end if
               if(cc%do_sc_complex=='Y') then
                  write (filn,'(''c'',a,''qw_proj.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                  call corr_write_compcorr(nq,cc%nw,3,filn,c_kw_proj(:,it,jt,:,:),cc%w)
                  if(cc%do_sc_list=='Y') then
                     write (filn,'(''c'',a,''qw_proj_l.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                     call corr_write_compcorr_list(nq,cc%nw,3,filn,c_kw_proj(:,it,jt,:,:),cc%w)
                  end if
               end if
            end do
         end do

         i_all=-product(shape(c_kw_proj))*kind(c_kw_proj)
         deallocate(c_kw_proj,stat=i_stat)
         call memocc(i_stat,i_all,'c_kw_proj','print_gkw')

      end if

      if(cc%do_projch=='Q'.or.cc%do_projch=='Y') then

         ! Allocate arrays
         allocate(c_kw_projch(3,Nchmax,Nchmax,nq,cc%nw),stat=i_stat)
         call memocc(i_stat,product(shape(c_kw_projch))*kind(c_kw_projch),'c_kw_projch','print_gkw')
         c_kw_projch=0.0_dblprec

         call combine_corr_proj_scalar(Nchmax, nq, 3, cc%nw, dc%m_kw_projch, dc%m_kw_projch, c_kw_projch)

         ! Write S(q,w)
         do it=1,nchmax
            do jt=1,nchmax
               write (filn,'(a,''qw_projch.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
               call corr_write_abscorr(nq,cc%nw,3,filn,c_kw_projch(:,it,jt,:,:),cc%w)
               if(cc%do_sc_list=='Y') then
                  write (filn,'(a,''qw_projch_l.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                  call corr_write_abscorr_list(nq,cc%nw,3,filn,c_kw_projch(:,it,jt,:,:),cc%w)
               end if
               if(cc%do_sc_complex=='Y') then
                  write (filn,'(''c'',a,''qw_projch_l.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                  call corr_write_compcorr(nq,cc%nw,3,filn,c_kw_projch(:,it,jt,:,:),cc%w)
                  if(cc%do_sc_list=='Y') then
                     write (filn,'(''c'',a,''qw_projch_list.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                     call corr_write_compcorr_list(nq,cc%nw,3,filn,c_kw_projch(:,it,jt,:,:),cc%w)
                  end if
               end if
            end do
         end do

         i_all=-product(shape(c_kw_projch))*kind(c_kw_projch)
         deallocate(c_kw_projch,stat=i_stat)
         call memocc(i_stat,i_all,'c_kw_projch','print_gkw')

      end if

      if(cc%do_sc_dosonly/="Y") then

         ! Write S(q,w)
         write (filn,'(a,''qw.'',a,''.out'')') trim(label),trim(simid)
         call corr_write_abscorr(nq,cc%nw,3,filn,c_kw,cc%w)
         if(cc%do_sc_list=='Y') then
            write (filn,'(a,''qw_l.'',a,''.out'')') trim(label),trim(simid)
            call corr_write_abscorr_list(nq,cc%nw,3,filn,c_kw,cc%w)
         end if

         ! Write S(q,w) decomposed in real and imaginary parts
         if(cc%do_sc_complex=='Y') then
            write (filn,'(''c'',a,''qw.'',a,''.out'')') trim(label),trim(simid)
            call corr_write_compcorr(nq,cc%nw,3,filn,c_kw,cc%w)
            if(cc%do_sc_list=='Y') then
               write (filn,'(''c'',a,''qw_l.'',a,''.out'')') trim(label),trim(simid)
               call corr_write_compcorr_list(nq,cc%nw,3,filn,c_kw,cc%w)
            end if
         end if

         ! Write the elements of tensorial S(q,w)
         if(cc%do_sc_tens=='Y') then

            ! Allocate the tensorial S(q,w)
            allocate(c_kw_tens(3,3,nq,cc%nw),stat=i_stat)
            call memocc(i_stat,product(shape(c_kw_tens))*kind(c_kw_tens),'c_kw_tens','print_gkw')
            c_kw_tens=0.0_dblprec

            call combine_corr_tensor(nq, 3, cc%nw, cc%m_kw, dc%m_kw, c_kw_tens)


            ! Write absolute values of the elements of tensorial S(q,w)
            write (filn,'(a,''qwtensa.'',a,''.out'')') trim(label),trim(simid)
            call corr_write_abscorr(nq,cc%nw,9,filn,c_kw_tens,cc%w)

            ! Write the real values of the elements of tensorial S(q,w)
            write (filn,'(a,''qwtensr.'',a,''.out'')') trim(label),trim(simid)
            call corr_write_recorr(nq,cc%nw,9,filn,c_kw_tens,cc%w)

            ! Write the imaginary values of the elements of tensorial S(q,w)
            ! Optionally write instead on polar notation, c.f. csqw.
            write (filn,'(a,''qwtensi.'',a,''.out'')') trim(label),trim(simid)
            call corr_write_aicorr(nq,cc%nw,9,filn,c_kw_tens,cc%w)

            ! Calculate the scattering intensity from the tensorial S(q,w)
            ! Write its real, imaginary, and absolute value to file.
            allocate(sqwintensity(nq,cc%nw),stat=i_stat)
            call memocc(i_stat,product(shape(sqwintensity))*kind(sqwintensity),'sqwintensity','print_gkw')
            sqwintensity=0.0_dblprec


            ! 
            unit3 = 0.0_dblprec
            unit3(1,1) = 1.0_dblprec;  unit3(2,2) = 1.0_dblprec; unit3(3,3) = 1.0_dblprec
            do iq=1,Nq
               qnorm = norm2(q(1:3,iq))
               qnorm2 = qnorm**2+1.0e-12_dblprec
               do iw=1,cc%Nw/2
                  do ia=1,3
                     do ib=1,3
                        polfac = unit3(ia,ib) - q(ia,iq) * q(ib,iq) / qnorm2
                        sqwintensity(iq,iw) = sqwintensity(iq,iw) + polfac * c_kw_tens(ia,ib,iq,iw)
                     end do
                  end do
               end do
            end do

            write (filn,'(a,''qwintensity.'',a,''.out'')') trim(label),trim(simid)
            open(ofileno, file=filn)
            write (ofileno,'(a)') &
               "#    iq    q_x       q_y       q_x          iw       Re(Int(S))      Im(Int(S))      |Int(S)|   "
            do iq=1,Nq
               do iw=1,cc%Nw/2
                  write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, real(sqwintensity(iq,iw)), &
                     aimag(sqwintensity(iq,iw)), abs(sqwintensity(iq,iw))
               end do
            end do
            close(ofileno)
            !
            i_all=-product(shape(c_kw_tens))*kind(c_kw_tens)
            deallocate(c_kw_tens,stat=i_stat)
            call memocc(i_stat,i_all,'c_kw_tens','print_gkw')
            !

         end if

      end if

      ! Calculate and print DOS
      write (filn,'(a,''wdos.'',a,''.out'')') trim(label),trim(simid)
      call extract_dos(nq,cc%nw,3,cc%w,c_kw,filn)

      if (cc%do_sc_tens=='Y') then
         write (filn,'(a,''wdos_int.'',a,''.out'')') trim(label),trim(simid)
         call extract_dos(nq,cc%nw,1,cc%w,sqwintensity,filn)
      end if


      if (cc%do_sc_tens=='Y') then
         i_all=-product(shape(sqwintensity))*kind(sqwintensity)
         deallocate(sqwintensity,stat=i_stat)
         call memocc(i_stat,i_all,'sqwintensity','print_gkw')
      end if

      i_all=-product(shape(c_kw))*kind(c_kw)
      deallocate(c_kw,stat=i_stat)
      call memocc(i_stat,i_all,'c_kw','print_gkw')
      !

      10005 format (i7,3f10.6,2x,i7,2x,9es16.8)
      !
   end subroutine print_gkw

   subroutine print_gkt(NT, Nchmax, cc, dc, simid, label)

      use Constants
      !
      implicit none
      !
      integer, intent(in) :: NT           !< Number of types of atoms
      integer, intent(in) :: Nchmax       !< Number of chemical types

      type(corr_t), intent(inout) :: cc !< Derived type for correlation data
      type(corr_t), intent(inout) :: dc !< Derived type for correlation data

      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=2), intent(in) :: label

      !
      integer  :: i_stat, i_all
      integer :: it, jt
      character(len=30) :: filn
      complex(dblprec)  :: i
      !
      !
      complex(dblprec), dimension(:,:,:),     allocatable :: c_kt                ! Correlation for G(k,w)
      complex(dblprec), dimension(:,:,:,:,:), allocatable :: c_kt_proj           ! Correlation for G(k,w)
      complex(dblprec), dimension(:,:,:,:,:), allocatable :: c_kt_projch         ! Correlation for G(k,w)


      i=(0.0_dblprec,1.0_dblprec)



      ! Allocate arrays

      allocate(c_kt(3,nq,cc%sc_max_nstep),stat=i_stat)
      call memocc(i_stat,product(shape(c_kt))*kind(c_kt),'c_kt','print_gkt')
      c_kt=0.0_dblprec

      call combine_corr_scalar(nq, 3, cc%sc_max_nstep, cc%m_kt, dc%m_kt, c_kt)
      c_kt = c_kt / (0.5_dblprec*(cc%sc_tidx-1))  / (0.5_dblprec*(dc%sc_tidx-1))

      ! Write S(q,t)
      write (filn,'(a,''qt.'',a,''.out'')') trim(label),trim(simid)
      call corr_write_abscorr(nq,cc%sc_max_nstep,3,filn,c_kt)

      ! Write S(q,w) decomposed in real and imaginary parts
      if(cc%do_sc_complex=='Y') then
         write (filn,'(''c'',a,''qt.'',a,''.out'')') trim(label),trim(simid)
         call corr_write_compcorr(nq,cc%sc_max_nstep,3,filn,c_kt)
      end if


      i_all=-product(shape(c_kt))*kind(c_kt)
      deallocate(c_kt,stat=i_stat)
      call memocc(i_stat,i_all,'c_kt','print_gkt')
      !

      if(cc%do_proj=='T'.or.cc%do_proj=='Y') then

         ! Allocate arrays
         allocate(c_kt_proj(3,nt,nt,nq,cc%sc_max_nstep),stat=i_stat)
         call memocc(i_stat,product(shape(c_kt_proj))*kind(c_kt_proj),'c_kt_proj','print_gkt')
         c_kt_proj=0.0_dblprec

         call combine_corr_proj_scalar(nt, nq, 3, cc%sc_max_nstep, cc%m_kt_proj, cc%m_kt_proj, c_kt_proj)

         ! Write S(q,w)
         do it=1,nt
            do jt=1,nt
               write (filn,'(a,''qt_proj.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
               call corr_write_abscorr(nq,cc%sc_max_nstep,3,filn,c_kt_proj(:,it,jt,:,:))
               if(cc%do_sc_complex=='Y') then
                  write (filn,'(''c'',a,''qt_proj.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                  call corr_write_compcorr(nq,cc%sc_max_nstep,3,filn,c_kt_proj(:,it,jt,:,:))
               end if
            end do
         end do

         i_all=-product(shape(c_kt_proj))*kind(c_kt_proj)
         deallocate(c_kt_proj,stat=i_stat)
         call memocc(i_stat,i_all,'c_kt_proj','print_gkt')

      end if

      if(cc%do_projch=='T'.or.cc%do_projch=='Y') then

         ! Allocate arrays
         allocate(c_kt_projch(3,Nchmax,Nchmax,nq,cc%sc_max_nstep),stat=i_stat)
         call memocc(i_stat,product(shape(c_kt_projch))*kind(c_kt_projch),'c_kt_projch','print_gkt')
         c_kt_projch=0.0_dblprec

         call combine_corr_proj_scalar(Nchmax, nq, 3, cc%sc_max_nstep, dc%m_kt_projch, dc%m_kt_projch, c_kt_projch)

         ! Write S(q,t)
         do it=1,nchmax
            do jt=1,nchmax
               write (filn,'(a,''qt_projch.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
               call corr_write_abscorr(nq,cc%sc_max_nstep,3,filn,c_kt_projch(:,it,jt,:,:))
               if(cc%do_sc_complex=='Y') then
                  write (filn,'(''c'',a,''qt_projch.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                  call corr_write_compcorr(nq,cc%sc_max_nstep,3,filn,c_kt_projch(:,it,jt,:,:))
               end if
            end do
         end do

         i_all=-product(shape(c_kt_projch))*kind(c_kt_projch)
         deallocate(c_kt_projch,stat=i_stat)
         call memocc(i_stat,i_all,'c_kt_projch','print_gkt')

      end if


      return

      !
   end subroutine print_gkt

   subroutine print_gk(NT, Nchmax, cc, dc, simid, label)

      use Constants
      !
      implicit none
      !
      integer, intent(in) :: NT           !< Number of types of atoms
      integer, intent(in) :: Nchmax       !< Number of chemical types

      type(corr_t), intent(inout) :: cc !< Derived type for correlation data
      type(corr_t), intent(inout) :: dc !< Derived type for correlation data

      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=2), intent(in) :: label

      !
      integer  :: iq, iw, i_stat, i_all
      integer :: ia,ib, it, jt
      character(len=30) :: filn
      complex(dblprec)  :: i
      real(dblprec), dimension(3,3) :: unit3
      real(dblprec) :: qnorm, qnorm2, polfac
      !
      !
      complex(dblprec), dimension(:,:),     allocatable :: c_k                ! Correlation for G(k,w)
      complex(dblprec), dimension(:,:,:,:), allocatable :: c_k_proj           ! Correlation for G(k,w)
      complex(dblprec), dimension(:,:,:,:), allocatable :: c_k_projch         ! Correlation for G(k,w)
      complex(dblprec), dimension(:,:,:),   allocatable :: c_k_tens           ! Correlation for G(k,w)

      complex(dblprec), dimension(:),   allocatable :: sqintensity

      i=(0.0_dblprec,1.0_dblprec)

      ! Allocate arrays

      allocate(c_k(3,nq),stat=i_stat)
      call memocc(i_stat,product(shape(c_k))*kind(c_k),'c_k','print_gk')
      c_k=0.0_dblprec

      ! AB here
      call combine_corr_scalar(nq, 3, 1, cc%m_k, dc%m_k, c_k)
      c_k = c_k / cc%sc_nsamp / dc%sc_nsamp

      ! Write S(q)
      write (filn,'(a,''q.'',a,''.out'')') trim(label),trim(simid)
      call corr_write_signcorr(nq,1,3,filn,c_k)

      ! Write S(q,w) decomposed in real and imaginary parts
      if(cc%do_sc_complex=='Y') then
         write (filn,'(''c'',a,''q.'',a,''.out'')') trim(label),trim(simid)
         call corr_write_compcorr(nq,1,3,filn,c_k)
      end if

      if(cc%do_proj=='C'.or.cc%do_proj=='Y') then

         ! Allocate arrays
         allocate(c_k_proj(3,nt,nt,nq),stat=i_stat)
         call memocc(i_stat,product(shape(c_k_proj))*kind(c_k_proj),'c_k_proj','print_gk')
         c_k_proj=0.0_dblprec

         call combine_corr_proj_scalar(nt, nq, 3, 1, cc%m_k_proj, cc%m_k_proj, c_k_proj)
         c_k_proj = c_k_proj / cc%sc_nsamp / dc%sc_nsamp

         ! Write S(q)
         do it=1,nt
            do jt=1,nt
               write (filn,'(a,''q_proj.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
               call corr_write_signcorr(nq,1,3,filn,c_k_proj(:,it,jt,:))
               if(cc%do_sc_complex=='Y') then
                  write (filn,'(''c'',a,''q_proj.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                  call corr_write_compcorr(nq,1,3,filn,c_k_proj(:,it,jt,:))
               end if
            end do
         end do

         i_all=-product(shape(c_k_proj))*kind(c_k_proj)
         deallocate(c_k_proj,stat=i_stat)
         call memocc(i_stat,i_all,'c_k_proj','print_gkw')

      end if

      if(cc%do_projch=='C'.or.cc%do_projch=='Y') then

         ! Allocate arrays
         allocate(c_k_projch(3,Nchmax,Nchmax,nq),stat=i_stat)
         call memocc(i_stat,product(shape(c_k_projch))*kind(c_k_projch),'c_k_projch','print_gk')
         c_k_projch=0.0_dblprec

         call combine_corr_proj_scalar(Nchmax, nq, 3, 1, dc%m_k_projch, dc%m_k_projch, c_k_projch)
         c_k_projch = c_k_projch / cc%sc_nsamp / dc%sc_nsamp

         ! Write S(q,w)
         do it=1,nchmax
            do jt=1,nchmax
               write (filn,'(a,''q_projch.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
               call corr_write_signcorr(nq,1,3,filn,c_k_projch(:,it,jt,:))
               if(cc%do_sc_complex=='Y') then
                  write (filn,'(''c'',a,''q_projch.'',i0,''.'',i0,''.'',a,''.out'')') trim(label),it,jt,trim(simid)
                  call corr_write_compcorr(nq,1,3,filn,c_k_projch(:,it,jt,:))
               end if
            end do
         end do

         i_all=-product(shape(c_k_projch))*kind(c_k_projch)
         deallocate(c_k_projch,stat=i_stat)
         call memocc(i_stat,i_all,'c_k_projch','print_gkw')

      end if


      ! Write the elements of tensorial S(q,w)
      if(cc%do_sc_tens=='Y') then

         ! Allocate the tensorial S(q,w)
         allocate(c_k_tens(3,3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(c_k_tens))*kind(c_k_tens),'c_k_tens','print_gk')
         c_k_tens=0.0_dblprec

         call combine_corr_tensor(nq, 3, 1, cc%m_k, dc%m_k, c_k_tens)
         c_k_tens = c_k_tens / cc%sc_nsamp / dc%sc_nsamp


         ! Write absolute values of the elements of tensorial S(q,w)
         write (filn,'(a,''qtensa.'',a,''.out'')') trim(label),trim(simid)
         call corr_write_abscorr(nq,1,9,filn,c_k_tens)

         ! Write the real values of the elements of tensorial S(q,w)
         write (filn,'(a,''qtensr.'',a,''.out'')') trim(label),trim(simid)
         call corr_write_recorr(nq,1,9,filn,c_k_tens)

         ! Write the imaginary values of the elements of tensorial S(q,w)
         ! Optionally write instead on polar notation, c.f. csqw.
         write (filn,'(a,''qtensi.'',a,''.out'')') trim(label),trim(simid)
         call corr_write_aicorr(nq,1,9,filn,c_k_tens)

         ! Calculate the scattering intensity from the tensorial S(q,w)
         ! Write its real, imaginary, and absolute value to file.
         allocate(sqintensity(nq),stat=i_stat)
         call memocc(i_stat,product(shape(sqintensity))*kind(sqintensity),'sqintensity','print_gk')
         sqintensity=0.0_dblprec


         ! 
         unit3 = 0.0_dblprec
         unit3(1,1) = 1.0_dblprec;  unit3(2,2) = 1.0_dblprec; unit3(3,3) = 1.0_dblprec
         do iq=1,Nq
            qnorm = norm2(q(1:3,iq))
            qnorm2 = qnorm**2+1.0e-12_dblprec
            do ia=1,3
               do ib=1,3
                  polfac = unit3(ia,ib) - q(ia,iq) * q(ib,iq) / qnorm2
                  sqintensity(iq) = sqintensity(iq) + polfac * c_k_tens(ia,ib,iq)
               end do
            end do
         end do

         write (filn,'(a,''qintensity.'',a,''.out'')') trim(label),trim(simid)
         open(ofileno, file=filn)
         iw=1
          write (ofileno,'(a)') &
         " #   iq    q_x        q_y      q_x          qw      Re(Int(S))      Im(Int(S))        |Int(S)|"
         do iq=1,Nq
            write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, real(sqintensity(iq)), &
               aimag(sqintensity(iq)), abs(sqintensity(iq))
         end do
         close(ofileno)
         !
         i_all=-product(shape(c_k_tens))*kind(c_k_tens)
         deallocate(c_k_tens,stat=i_stat)
         call memocc(i_stat,i_all,'c_k_tens','print_gk')
         !

      end if


      if (cc%do_sc_tens=='Y') then
         i_all=-product(shape(sqintensity))*kind(sqintensity)
         deallocate(sqintensity,stat=i_stat)
         call memocc(i_stat,i_all,'sqintensity','print_gk')
      end if

      i_all=-product(shape(c_k))*kind(c_k)
      deallocate(c_k,stat=i_stat)
      call memocc(i_stat,i_all,'c_k','print_gk')
      !

      10005 format (i7,3f10.6,2x,i7,2x,9es16.8)
      !
   end subroutine print_gk

   subroutine print_gr(Natom,  cc,  coord, simid)

      use Constants
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system

      type(corr_t), intent(inout) :: cc !< Derived type for correlation data

      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      character(len=8), intent(in) :: simid !< Name of simulation

      !
      integer  :: iq, r, i_all, i_stat
      integer :: l
      character(len=30) :: filn
      complex(dblprec)  :: i, epowqr, iqfac
      real(dblprec), dimension(3) :: cl
      complex(dblprec), dimension(3) :: cl_step, s0, sp
      real(dblprec) :: qdr, k_min, qfac
      !
      !
      complex(dblprec), dimension(:,:),     allocatable :: c_r                ! Correlation for G(k,w)

      !-- local 
      qfac=2._dblprec*pi
      iqfac=-2._dblprec*pi*(0.0_dblprec,1.0_dblprec)
      i=(0.0_dblprec,1.0_dblprec)

      ! Calculate the correlation length following the Katzgraber recipe
      s0=cc%m_k(:,qmin(1))*conjg(cc%m_k(:,qmin(1)))
      sp=cc%m_k(:,qmin(2))*conjg(cc%m_k(:,qmin(2)))
      k_min=sqrt(q(1,qmin(2))**2+q(2,qmin(2))**2+q(3,qmin(2))**2)
      cl_step = (s0/sp-1.0_dblprec)/4.0_dblprec/sin((qfac*k_min/2.0_dblprec))**2
      cl=real(sqrt(cl_step))
      write(*,'(2x,a20,2x,f11.5,2x,f11.5,2x,g15.5)')'Correlation lengths:',cl(1),cl(2),cl(3)

      ! Transform G(k) to G(r)

      allocate(c_r(3,Natom),stat=i_stat)
      call memocc(i_stat,product(shape(c_r))*kind(c_r),'c_r','print_gr')
      c_r=0.0_dblprec

      !$omp parallel do default(shared) private(r,iq,qdr,epowqr) schedule(static)
      do r=1,Natom
         do iq=1,nq
            qdr=q(1,iq)*(coord(1,r)-r_mid(1))+q(2,iq)*(coord(2,r)-r_mid(2))+q(3,iq)*(coord(3,r)-r_mid(3))
            epowqr=exp(-iqfac*qdr)*sc_window_fac(sc_window_fun,iq,nq)
            c_r(:,r)=c_r(:,r)+epowqr*cc%m_k(:,iq)*conjg(cc%m_k(:,iq))
         end do
      end do
      !$omp end parallel do

      ! Write G(r)
      !r_mid=0.0_dblprec
      write (filn,'(a,''r.'',a,''.out'')') trim(cc%label),trim(simid)
      open(ofileno,file=filn,status='replace')
      write (ofileno,'(a)') &
      "#       ir     r_x        r_y      r_x         Re(S(r)_xx)      Re(S(r)_yy)       Re(S(r)_zz)        |S(r)|        Sx+Sy+Sz    "
      do r=1,Natom
         write(ofileno,'(i10,3f10.4,5f18.8)') r,(coord(l,r)-r_mid(l),l=1,3),((real(c_r(l,r))),l=1,3),&
            real(sqrt(c_r(1,r)**2+c_r(2,r)**2+c_r(3,r)**2)),real(c_r(1,r)+c_r(2,r)+c_r(3,r))
      end do
      close(ofileno)

      ! Write G(|r|)
      write (filn,'(a,''ra.'',a,''.out'')') trim(cc%label),trim(simid)
      open(ofileno,file=filn,status='replace')
      write (ofileno,'(a)') &
      "#        |r|                Re(S(r)_xx)      Re(S(r)_yy)       Re(S(r)_zz)           |S(r)|            Sx+Sy+Sz      "
      do r=1,Natom
         !write(ofileno,'(7f18.8)') sqrt( (coord(1,r)-r_mid(1))**2+(coord(2,r)-r_mid(2))**2+(coord(3,r)-r_mid(3))**2),&
         !   (((c_r(l,r))),l=1,3),&
         !   sqrt(c_r(1,r)**2+c_r(2,r)**2+c_r(3,r)**2),c_r(1,r)+c_r(2,r)+c_r(3,r)
         write(ofileno,'(f18.8,3x,3f18.8,2x,2f18.8)') norm2(coord(:,r)-r_mid), real(c_r(:,r)), & 
            sqrt(real(sum(conjg(c_r(:,r))*(c_r(:,r))))),real(sum(c_r(:,r)))
      end do
      close(ofileno)

      ! Deallocate arrays
      i_all=-product(shape(c_r))*kind(c_r)
      deallocate(c_r,stat=i_stat)
      call memocc(i_stat,i_all,'corr_r','print_gr')

      return

   end subroutine print_gr

   subroutine corr_write_abscorr(nq,nw,nelem,filn,corr_out,w_arr)
      !
      implicit none

      integer, intent(in) :: nq
      integer, intent(in) :: nw
      integer, intent(in) :: nelem
      character(len=30) :: filn
      complex(dblprec), dimension(nelem,nq,nw) :: corr_out
      real(dblprec), dimension(nw), intent(in), optional :: w_arr

      integer :: iq, iw

      open(ofileno, file=filn)
      if (print_real_w .and. present(w_arr))  then
         if (nelem==9) then
            write (ofileno,'(a,a)') &
               "#    iq    q_x        q_y      q_x             w(Hz)            |S_xx|           |S_xy|          |S_xz|",&
               "         |S_yx|            |yy|          |S_yz|         |S_zx|          |S_zy|          |S_zz|         |S|"
         else
            write (ofileno,'(a)') &
               "#    iq    q_x        q_y      q_x             w(Hz)             |S_xx|          |S_yy|          |S_zz|         |S|      "
         end if
         do iq=1,Nq
            do iw=1,max(Nw/2,1)
               write (ofileno,10015) iq,q(1,iq), q(2,iq),q(3,iq),w_arr(iw), &
                  abs(corr_out(:,iq,iw)),abs(sum(conjg(corr_out(:,iq,iw))*corr_out(:,iq,iw))**0.5_dblprec)
            end do
         end do
      else
         if (nelem==9) then
            write (ofileno,'(a,a)') &
               "#    iq    q_x        q_y      q_x         qw         |S_xx|           |S_xy|          |S_xz|",&
               "         |S_yx|            |yy|          |S_yz|         |S_zx|          |S_zy|          |S_zz|         |S|"
         else
            write (ofileno,'(a)') &
               "#    iq    q_x        q_y      q_x        qw           |S_xx|          |S_yy|          |S_zz|         |S|      "
         end if
         do iq=1,Nq
            do iw=1,max(Nw/2,1)
               write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                  abs(corr_out(:,iq,iw)),abs(sum(conjg(corr_out(:,iq,iw))*corr_out(:,iq,iw))**0.5_dblprec)
            end do
         end do
      end if
      close(ofileno)

      10005 format (i7,3f10.6,2x,i7,2x,17es16.8)
      10015 format (i7,3f10.6,2x,es16.8,2x,17es16.8)

   end subroutine corr_write_abscorr

   subroutine corr_write_signcorr(nq,nw,nelem,filn,corr_out,w_arr)
      !
      implicit none

      integer, intent(in) :: nq
      integer, intent(in) :: nw
      integer, intent(in) :: nelem
      character(len=30) :: filn
      complex(dblprec), dimension(nelem,nq,nw) :: corr_out
      real(dblprec), dimension(nw), intent(in), optional :: w_arr

      integer :: iq, iw

      open(ofileno, file=filn)
      if (print_real_w .and. present(w_arr))  then
         write (ofileno,'(a)') &
            "#    iq    q_x        q_y      q_x              w(Hz)           Re(S_xx)        Re(S_yy)        Re(S_zz)        |S|      "
         do iq=1,Nq
            do iw=1,max(Nw/2,1)
               write (ofileno,10015) iq,q(1,iq), q(2,iq),q(3,iq),w_arr(iw), &
                  real(corr_out(:,iq,iw)),abs(sum(conjg(corr_out(:,iq,iw))*corr_out(:,iq,iw))**0.5_dblprec)
            end do
         end do
      else
         write (ofileno,'(a)') &
            "#    iq    q_x        q_y      q_x          qw       Re(S_xx)        Re(S_yy)        Re(S_zz)         |S|      "
         do iq=1,Nq
            do iw=1,max(Nw/2,1)
               write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                  real(corr_out(:,iq,iw)),abs(sum(conjg(corr_out(:,iq,iw))*corr_out(:,iq,iw))**0.5_dblprec)
            end do
         end do
      end if
      close(ofileno)

      10005 format (i7,3f10.6,2x,i7,2x,17es16.8)
      10015 format (i7,3f10.6,2x,es16.8,2x,17es16.8)

   end subroutine corr_write_signcorr

   subroutine corr_write_abscorr_list(nq,nw,nelem,filn,corr_out,w_arr)
      !
      implicit none

      integer, intent(in) :: nq
      integer, intent(in) :: nw
      integer, intent(in) :: nelem
      character(len=30) :: filn
      complex(dblprec), dimension(nelem,nq,nw) :: corr_out
      real(dblprec), dimension(nw), intent(in), optional :: w_arr

      integer :: iq
      integer, dimension(nelem) :: iidx

      open(ofileno, file=filn)
      write (ofileno,'(a)') &
         "#    iq    q_x        q_y      q_x          |S(q)|"
      do iq=1,Nq
         iidx=maxloc(abs(corr_out(:,iq,1:nw/2+1)),2)
         write (ofileno,10006) iq,q(1,iq), q(2,iq),q(3,iq),iidx
      end do
      close(ofileno)

      10006 format (i7,3f10.6,2x,12i6)

   end subroutine corr_write_abscorr_list

   subroutine corr_write_recorr(nq,nw,nelem,filn,corr_out,w_arr)
      !
      implicit none

      integer, intent(in) :: nq
      integer, intent(in) :: nw
      integer, intent(in) :: nelem
      character(len=30) :: filn
      complex(dblprec), dimension(nelem,nq,nw) :: corr_out
      real(dblprec), dimension(nw), intent(in), optional :: w_arr

      integer :: iq, iw

      open(ofileno, file=filn)
      if (print_real_w .and. present(w_arr))  then
         if (nelem==9) then
            write (ofileno,'(a,a)') &
               "#    iq    q_x        q_y      q_x             w(Hz)           Re(S_xx)        Re(S_xy)        Re(S_xz)",&
               "       Re(S_yx)        Re(S_yy)        Re(S_yz)       Re(S_zx)        Re(S_zy)        Re(S_zz)         |S|"
         else
            write (ofileno,'(a)') &
               "#    iq    q_x        q_y      q_x             w(Hz)           Re(S_xx)        Re(S_yy)        Re(S_zz)         |S|      "
         end if
         do iq=1,Nq
            do iw=1,max(1,Nw/2)
               write (ofileno,10015) iq,q(1,iq), q(2,iq),q(3,iq),w_arr(iw), &
                  real(corr_out(:,iq,iw)),real(sum(conjg(corr_out(:,iq,iw))*corr_out(:,iq,iw))**0.5_dblprec)
            end do
         end do
      else
         if (nelem==9) then
            write (ofileno,'(a,a)') &
               "#    iq    q_x        q_y      q_x          qw       Re(S_xx)        Re(S_xy)        Re(S_xz)",&
               "       Re(S_yx)        Re(S_yy)        Re(S_yz)       Re(S_zx)        Re(S_zy)        Re(S_zz)         |S|"
         else
            write (ofileno,'(a)') &
               "#    iq    q_x        q_y      q_x          qw       Re(S_xx)        Re(S_yy)        Re(S_zz)         |S|      "
         end if
         do iq=1,Nq
            do iw=1,max(1,Nw/2)
               write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                  real(corr_out(:,iq,iw)),real(sum(conjg(corr_out(:,iq,iw))*corr_out(:,iq,iw))**0.5_dblprec)
            end do
         end do
      end if
      close(ofileno)

      10005 format (i7,3f10.6,2x,i7,2x,17es16.8)
      10015 format (i7,3f10.6,2x,es16.8,2x,17es16.8)

   end subroutine corr_write_recorr

   subroutine corr_write_aicorr(nq,nw,nelem,filn,corr_out,w_arr)
      !
      implicit none

      integer, intent(in) :: nq
      integer, intent(in) :: nw
      integer, intent(in) :: nelem
      character(len=30) :: filn
      complex(dblprec), dimension(nelem,nq,nw) :: corr_out
      real(dblprec), dimension(nw), intent(in), optional :: w_arr

      integer :: iq, iw

      open(ofileno, file=filn)
      if (print_real_w .and. present(w_arr))  then
         if (nelem==9) then
            write (ofileno,'(a,a)') &
               "#    iq    q_x        q_y      q_x            w(Hz)         Im(S_xx)        Im(S_xy)        Im(S_xz)",&
               "       Im(S_yx)        Im(S_yy)        Im(S_yz)       Im(S_zx)        Im(S_zy)        Im(S_zz)         |S|"
         else
            write (ofileno,'(a)') &
               "#    iq    q_x        q_y      q_x            w(Hz)         Im(S_xx)        Im(S_yy)        Im(S_zz)         |S|      "
         end if
         do iq=1,Nq
            do iw=1,max(1,Nw/2)
               write (ofileno,10015) iq,q(1,iq), q(2,iq),q(3,iq),w_arr(iw), &
                  aimag(corr_out(:,iq,iw))!,real(sum(conjg(corr_out(:,iq,iw))*corr_out(:,iq,iw))**0.5_dblprec)
            end do
         end do
      else
         if (nelem==9) then
            write (ofileno,'(a,a)') &
               "#    iq    q_x        q_y      q_x          iw       Im(S_xx)        Im(S_xy)        Im(S_xz)",&
               "       Im(S_yx)        Im(S_yy)        Im(S_yz)       Im(S_zx)        Im(S_zy)        Im(S_zz)         |S|"
         else
            write (ofileno,'(a)') &
               "#    iq    q_x        q_y      q_x          iw       Im(S_xx)        Im(S_yy)        Im(S_zz)         |S|      "
         end if
         do iq=1,Nq
            do iw=1,max(1,Nw/2)
               write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
                  aimag(corr_out(:,iq,iw))!,real(sum(conjg(corr_out(:,iq,iw))*corr_out(:,iq,iw))**0.5_dblprec)
            end do
         end do
      end if
      close(ofileno)

      10005 format (i7,3f10.6,2x,i7,2x,17es16.8)
      10015 format (i7,3f10.6,2x,es16.8,2x,17es16.8)

   end subroutine corr_write_aicorr

   subroutine corr_write_compcorr(nq,nw,nelem,filn,corr_out,w_arr)
      !
      implicit none

      integer, intent(in) :: nq
      integer, intent(in) :: nw
      integer, intent(in) :: nelem
      character(len=30) :: filn
      complex(dblprec), dimension(nelem,nq,nw) :: corr_out
      real(dblprec), dimension(nw), intent(in), optional :: w_arr

      integer :: iq, iw



      open(ofileno, file=filn)
      do iq=1,Nq
         do iw=1,max(1,Nw/2)
            write (ofileno,10005) iq,q(1,iq), q(2,iq),q(3,iq),iw, &
               abs(corr_out(1,iq,iw)),atan2(aimag(corr_out(1,iq,iw)),real(corr_out(1,iq,iw))), &
               abs(corr_out(2,iq,iw)),atan2(aimag(corr_out(2,iq,iw)),real(corr_out(2,iq,iw))), &
               abs(corr_out(3,iq,iw)),atan2(aimag(corr_out(3,iq,iw)),real(corr_out(3,iq,iw))), &
               abs(sum(conjg(corr_out(:,iq,iw))*corr_out(:,iq,iw))**0.5_dblprec)
         end do
      end do
      close(ofileno)

      10005 format (i7,3f10.6,2x,i7,2x,17es16.8)

   end subroutine corr_write_compcorr

   subroutine corr_write_compcorr_list(nq,nw,nelem,filn,corr_out,w_arr)
      !
      implicit none

      integer, intent(in) :: nq
      integer, intent(in) :: nw
      integer, intent(in) :: nelem
      character(len=30) :: filn
      complex(dblprec), dimension(nelem,nq,nw) :: corr_out
      real(dblprec), dimension(nw), intent(in), optional :: w_arr

      integer :: iq
      integer, dimension(nelem) :: iidx



      open(ofileno, file=filn)
      do iq=1,Nq
         iidx=maxloc(abs(corr_out(:,iq,1:nw/2+1)),2)
         write (ofileno,10006) iq,q(1,iq), q(2,iq),q(3,iq), &
               iidx(1), abs(corr_out(1,iq,iidx(1))),atan2(aimag(corr_out(1,iq,iidx(1))),real(corr_out(1,iq,iidx(1)))), &
               iidx(2), abs(corr_out(2,iq,iidx(2))),atan2(aimag(corr_out(2,iq,iidx(2))),real(corr_out(2,iq,iidx(2)))), &
               iidx(3), abs(corr_out(3,iq,iidx(3))),atan2(aimag(corr_out(3,iq,iidx(3))),real(corr_out(3,iq,iidx(3))))
      end do
      close(ofileno)

      10006 format (i6,3f10.6,2x,i6,2es16.8,1x,i6,2es16.8,1x,i6,2es16.8)

   end subroutine corr_write_compcorr_list



end module Correlation_Print
