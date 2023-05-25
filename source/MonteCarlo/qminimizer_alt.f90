!!!    !-----------------------------------------------------------------------------
!!!    ! SUBROUTINE: mini_q
!!!    !> @brief Main driver for the minimization of spin spirals
!!!    !> @author Anders Bergman
!!!    !-----------------------------------------------------------------------------
!!!    subroutine mini_q(Natom,Mensemble,NA,coord,do_jtensor,exc_inter,do_dm,do_pd,     &
!!!       do_biqdm,do_bq,do_chir,taniso,sb,do_dip,emomM,mmom,hfield,OPT_flag,           &
!!!       max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations,    &
!!!       constellationsNeighType,mult_axis,Num_macro,cell_index,emomM_macro,           &
!!!       macro_nlistsize,do_anisotropy)
!!! 
!!!       use RandomNumbers, only: rng_uniform,rng_gaussian, use_vsl
!!!       use Constants, only : mub, mry
!!!       use InputData, only : N1,N2,N3
!!!       !
!!!       !.. Implicit declarations
!!!       implicit none
!!! 
!!!       integer, intent(in) :: Natom !< Number of atoms in system
!!!       integer, intent(in) :: Mensemble !< Number of ensembles
!!!       integer, intent(in) :: NA  !< Number of atoms in one cell
!!!       real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
!!!       integer, intent(in) :: do_jtensor   !<  Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
!!!       character(len=1),intent(in) :: exc_inter !< Interpolate Jij (Y/N)
!!!       integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_biqdm   !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_chir  !< Add scalar chirality exchange (CHIR) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
!!!       integer, intent(in) :: do_anisotropy
!!!       integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
!!!       real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between the anisotropies
!!!       real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
!!!       real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
!!!       real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
!!!       character(len=1), intent(in) :: mult_axis
!!!       !! +++ New variables due to optimization routines +++ !!
!!!       integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
!!!       ! Number of entries (number of unit cells*number of atoms per unit cell) in the constellation matrix per ensemble
!!!       integer, dimension(Mensemble), intent(in) :: maxNoConstl
!!!       ! See OptimizationRoutines.f90 for details on classification
!!!       integer, dimension(Natom, Mensemble), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
!!!       ! Matrix relating the interatomic exchanges for each atom in the constellation matrix
!!!       real(dblprec), dimension(ham%max_no_neigh, max_no_constellations,Mensemble), intent(in) :: constlNCoup
!!!       ! Matrix storing all unit cells belonging to any constellation
!!!       real(dblprec), dimension(3,max_no_constellations, Mensemble), intent(in) :: constellations
!!!       ! Optimization flag (1 = optimization on; 0 = optimization off)
!!!       logical, intent(in) :: OPT_flag
!!!       ! Matrix storing the type of the neighbours within a given neighbourhood of a constellation; default is 1 outside the neighbourhood region
!!!       ! The default is to achieve correct indexing. Note here also that constlNCoup will result in a net zero contribution to the Heissenberg exchange term
!!!       integer, dimension(ham%max_no_neigh,max_no_constellations,Mensemble), intent(in) :: constellationsNeighType
!!!       ! Internal effective field arising from the optimization of the Heissenberg exchange term
!!!       integer, intent(in) :: Num_macro !< Number of macrocells in the system
!!!       integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
!!!       integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
!!!       real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
!!!       !
!!!       integer :: iq
!!!       !
!!!       real(dblprec), dimension(3) :: m_i, m_j, r_i, r_j, m_avg
!!!       real(dblprec) :: pi, qr, rx,ry, rz
!!!       integer :: i,j, k, iia, ia, ja, lhit, nhits, countstart
!!!       real(dblprec) :: energy, min_energy
!!!       !
!!!       real(dblprec), dimension(3,3) :: q_best, q_diff, s_diff, q_0, s_0, q_norm
!!!       real(dblprec), dimension(3,3) :: s_save, q_save
!!!       real(dblprec), dimension(3) :: theta_save, theta_diff, theta_best
!!!       real(dblprec) :: theta_glob_save, cone_ang, cone_ang_save
!!!       real(dblprec) :: theta_glob_diff, theta_glob_best
!!!       real(dblprec), dimension(3) :: phi_save, phi_diff, phi_best
!!!       real(dblprec), dimension(1) :: cone_ang_diff
!!!       real(dblprec), dimension(1) :: rng_tmp_arr
!!!       integer :: niter, iter, iscale, i1 ,i2 ,i3
!!!       real(dblprec) :: q_range, theta_range, s_range, phi_range, theta_glob_range, cone_ang_range
!!!       real(dblprec) :: q_range0, theta_range0, s_range0, phi_range0, theta_glob_range0
!!!       real(dblprec), dimension(3) :: srvec 
!!!       !
!!!       !
!!!       real(dblprec) :: q_start=1.0_dblprec
!!!       !
!!!       !
!!!       pi=4._dblprec*ATAN(1._dblprec)
!!!       theta_glob_best=0.0_dblprec
!!!       theta_best=0.0_dblprec
!!!       phi_best=0.0_dblprec
!!!       !
!!!       ! Normal vector
!!!       n_vec(1)=0.0_dblprec;n_vec(2)=0.0_dblprec;n_vec(3)=1.0_dblprec;
!!!       !
!!!       ! Starting atom
!!!       I1 = N1/2
!!!       I2 = N2/2
!!!       I3 = N3/2
!!! 
!!!       countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
!!!       print *, 'CountStart: ',countstart
!!!       !
!!!       niter=10000
!!!       ! Hard wired to planar spin spirals with 360/nq degree angle inbetween
!!!       !
!!!       theta = 0.0_dblprec
!!!       theta_glob = 0.0_dblprec
!!!       !theta_glob = pi/4.0_dblprec
!!!       phi = 0.0_dblprec
!!!       do iq=1,nq
!!!          ! Start along [100] direction
!!!          q(1,iq)=1.0_dblprec;q(2,iq)=0.0_dblprec;q(3,iq)=0.0_dblprec
!!!          !! Start along [100] direction
!!!          call normalize(q(1:3,iq))
!!!          ! Rotate 360/iq
!!!          theta(iq)=2.0_dblprec*pi/nq*(iq-1.0_dblprec)
!!!          rx= q(1,iq)*cos(theta(iq)+theta_glob)+q(2,iq)*sin(theta(iq)+theta_glob)
!!!          ry=-q(1,iq)*sin(theta(iq)+theta_glob)+q(2,iq)*cos(theta(iq)+theta_glob)
!!!          q(1,iq)=rx
!!!          q(2,iq)=ry
!!!          q(3,iq)=0.0_dblprec
!!!          !
!!!          q_norm(:,iq)=q(:,iq)/sqrt(sum(q(:,iq)*q(:,iq))+1.0e-12_dblprec)
!!!          !
!!!          ! Create pitch vectors perpendicular to q
!!!          if(norm2(s(:,iq))==0.0_dblprec) then
!!!             s(1,iq)=q_norm(2,iq)*n_vec(3)-q_norm(3,iq)*n_vec(2)
!!!             s(2,iq)=q_norm(3,iq)*n_vec(1)-q_norm(1,iq)*n_vec(3)
!!!             s(3,iq)=q_norm(1,iq)*n_vec(2)-q_norm(2,iq)*n_vec(1)
!!!          end if
!!!          call normalize(s(1:3,iq))
!!!          !
!!!          print *,'----Q-and-S-vectors----',iq
!!!          print '(3f10.4)', q(:,iq)
!!!          print '(3f10.4)', s(:,iq)
!!!          !
!!!          ! Currently hard wired starting guess
!!!          q(:,iq)=q(:,iq)/sqrt(sum(q(:,iq)*q(:,iq))) *q_start  !/30.0_dblprec
!!!          !
!!!       end do
!!!       q_0=q
!!!       ! For Neel spirals:
!!!       !s=q_norm
!!!       !do iq=1,nq
!!!       !   print *,'----Q_0 vector----',iq
!!!       !   print '(3f10.4)', q_0(:,iq)
!!!       !end do
!!!       s_0=s
!!!       min_energy=1.0d4
!!!       ! Set starting minimization ranges
!!!       q_range0=1.0_dblprec
!!!       s_range0=0
!!!       theta_range0=0
!!!       theta_glob_range0= 1.0_dblprec
!!!       phi_range0=0.0_dblprec
!!!       cone_ang_range=0
!!!       ! Real ranges (scaled during the minimization)
!!!       q_range=q_range0
!!!       s_range=s_range0
!!!       theta_range=theta_range0
!!!       theta_glob_range=theta_glob_range0
!!!       phi_range=phi_range0
!!!       !
!!!       lhit=0
!!!       nhits=0
!!!       iscale=1
!!! 
!!!       ! Legacy code for 2d-systems
!!!       fac_2d=1.0_dblprec
!!!       ! Switch rotation direction
!!! 
!!!       ! Calculate total external field (not included yet)
!!!       do k=1,Mensemble
!!!          do i=1,Natom
!!!             external_field(1:3,i,k)= hfield
!!!             beff(1:3,i,k)=0.0_dblprec
!!!             beff1(1:3,i,k)=0.0_dblprec
!!!             beff2(1:3,i,k)=0.0_dblprec
!!!          end do
!!!       end do
!!! 
!!!       do iter=0,niter
!!! !         !
!!! !         q_diff=0.0_dblprec
!!! !         call rng_uniform(q_diff(1,1),1)
!!! !         q_diff(1,1)=2.0_dblprec*q_diff(1,1)-1.0_dblprec
!!! !         q_diff(2,1)=0.0_dblprec
!!! !         q_diff(3,1)=0.0_dblprec
!!! !         q_diff=q_diff*q_range
!!! !         !
!!! !         theta_save=theta
!!! !         phi_save=phi
!!! !         do iq=1,nq
!!! !            theta_diff(iq)=0.0_dblprec
!!! !            ! phi angle
!!! !            call rng_uniform(phi_diff(iq),1)
!!! !            phi_diff(iq)=0.0_dblprec
!!! !            !
!!! !         end do
!!! !         theta=theta+theta_diff
!!! !         phi=phi+phi_diff
!!! !         !
!!! !         ! global theta angle
!!! !         theta_glob_save=theta_glob
!!! !         call rng_uniform(rng_tmp_arr,1)
!!! !         theta_glob_diff=rng_tmp_arr(1)
!!! !         theta_glob_diff=(theta_glob_diff-0.5_dblprec)*theta_glob_range*2*pi
!!! !         theta_glob=theta_glob+theta_glob_diff
!!! !         ! Set up trial vectors
!!! !         energy=0.0_dblprec
!!! !         q_save=q
!!! !         q(:,1)=q(:,1)+q_diff(:,1)
!!! !         ! Local rotations
!!! !         do iq=1,nq
!!! !            rx= q(1,1)*cos(theta(iq))+q(2,1)*sin(theta(iq))
!!! !            ry=-q(1,1)*sin(theta(iq))+q(2,1)*cos(theta(iq))
!!! !            q(1,iq)=rx
!!! !            q(2,iq)=ry
!!! !            q(3,iq)=0.0_dblprec
!!! !         end do
!!! !         ! Global rotation
!!! !         do iq=1,nq
!!! !            rx= q(1,iq)*cos(theta_glob)+q(2,iq)*sin(theta_glob)
!!! !            ry=-q(1,iq)*sin(theta_glob)+q(2,iq)*cos(theta_glob)
!!! !            q(1,iq)=2.0_dblprec*((0.1_dblprec*iter)/(1.0_dblprec*niter)-0.05_dblprec)
!!! !            q(2,iq)=0.0_dblprec
!!! !            q(3,iq)=0.0_dblprec
!!! !         end do
!!!          do iq=1,nq
!!!             q(:,iq)=q_0(:,iq)-q_0(:,iq)/sqrt(sum(q_0(:,iq)**2))*(2.0_dblprec*q_start*iter)/(1.0_dblprec*niter)
!!!          end do
!!! !         !
!!!          s_save=s
!!!          cone_ang=0.0_dblprec
!!! 
!!!          ! Fold back to wanted range
!!!          do iq=1,nq
!!!             if(theta(iq)>pi) theta(iq)=theta(iq)-pi*2.0_dblprec
!!!             if(theta(iq)<-pi) theta(iq)=theta(iq)+pi*2.0_dblprec
!!!             if(phi(iq)>pi) phi(iq)=phi(iq)-pi*2.0_dblprec
!!!             if(phi(iq)<-pi) phi(iq)=phi(iq)+pi*2.0_dblprec
!!!             if(q(1,iq)>1.0_dblprec) q(1,iq)=q(1,iq)-2.0_dblprec
!!!             if(q(2,iq)>1.0_dblprec) q(2,iq)=q(2,iq)-2.0_dblprec
!!!             if(q(3,iq)>1.0_dblprec) q(3,iq)=q(3,iq)-2.0_dblprec
!!!             if(q(1,iq)<-1.0_dblprec) q(1,iq)=q(1,iq)+2.0_dblprec
!!!             if(q(2,iq)<-1.0_dblprec) q(2,iq)=q(2,iq)+2.0_dblprec
!!!             if(q(3,iq)<-1.0_dblprec) q(3,iq)=q(3,iq)+2.0_dblprec
!!!          end do
!!!          !
!!!          ! Set up spin-spiral magnetization (only first cell+ neighbours)
!!!          do k=1,Mensemble
!!!             energy=0.0_dblprec
!!!             do ia=1,Natom
!!!             !  lhit=lhit+1
!!!                srvec=coord(:,ia)-coord(:,countstart+1)
!!!                !
!!!                m_j=0.0_dblprec
!!!                do iq=1,nq
!!!                   qr=q(1,iq)*srvec(1)+q(2,iq)*srvec(2)+q(3,iq)*srvec(3)
!!!                   m_j=m_j+n_vec*cos(2*pi*qr+phi(iq))+s(:,iq)*sin(2*pi*qr+phi(iq))
!!!                end do
!!!                call normalize(m_j)
!!!                !emom(1:3,ia,k)=m_j
!!!                emomM(1:3,ia,k)=m_j*mmom(ia,k)
!!!                !write(ofileno,'(2i8,4f14.6)') 1,lhit,mmom(ia,k),m_j
!!!             !  write(ofileno,'(2i8,4f14.6)') 1,ia,mmom(ia,k),m_j
!!!             end do
!!! !        do ia=1,NA
!!! !        !  lhit=lhit+1
!!! !           iia=ia+countstart
!!! !           srvec=coord(:,iia)-coord(:,countstart+1)
!!! !           !
!!! !           m_j=0.0_dblprec
!!! !           do iq=1,nq
!!! !              qr=q(1,iq)*srvec(1)+q(2,iq)*srvec(2)+q(3,iq)*srvec(3)
!!! !              m_j=m_j+n_vec*cos(2*pi*qr+phi(iq))+s(:,iq)*sin(2*pi*qr+phi(iq))
!!! !           end do
!!! !           call normalize(m_j)
!!! !           !emom(1:3,ia,k)=m_j
!!! !           emomM(1:3,iia,k)=m_j*mmom(iia,k)
!!! !           do j=1,ham%nlistsize(iia)
!!! !              ja=ham%nlist(j,iia)
!!! !              srvec=coord(:,ja) -coord(:,countstart+1)
!!! !              m_j=0.0_dblprec
!!! !              do iq=1,nq
!!! !                 !
!!! !                 qr=q(1,iq)*srvec(1)+q(2,iq)*srvec(2)+q(3,iq)*srvec(3)
!!! !                 m_j=m_j+n_vec*cos(2*pi*qr+phi(iq))+s(:,iq)*sin(2*pi*qr+phi(iq))
!!! !                 !
!!! !              end do
!!! !              call normalize(m_j)
!!! !              m_j=m_j*mmom(ja,k)
!!! !              emomM(1:3,ja,k)=m_j
!!! !              !
!!! !           end do
!!! !        end do
!!!          !  do ia=1,NA
!!!          !     iia=ia+countstart
!!!          !     m_i=0.0_dblprec
!!!          !     do iq=1,nq
!!!          !        qr=q(1,iq)*(coord(1,iia)-coord(1,countstart+1))+q(2,iq)*(coord(2,iia)-coord(2,countstart+1))+q(3,iq)*(coord(3,iia)-coord(3,countstart+1))
!!!          !        !
!!!          !        r_i=n_vec*cos(2.0_dblprec*pi*qr+phi(iq))+s(:,iq)*sin(2.0_dblprec*pi*qr+phi(iq))
!!!          !        m_i=m_i+r_i
!!!          !        !
!!!          !     end do
!!!          !     call normalize(m_i)
!!!          !     m_i=m_i*mmom(iia,k)
!!!          !     !
!!!          !     emomM(1:3,iia,k)=m_i
!!!          !     do j=1,ham%nlistsize(iia)
!!!          !        ja=ham%nlist(j,iia)
!!!          !        m_j=0.0_dblprec
!!!          !        do iq=1,nq
!!!          !           qr=q(1,iq)*(coord(1,ja)-coord(1,countstart+1))+q(2,iq)*(coord(2,ja)-coord(2,countstart+1))+q(3,iq)*(coord(3,ja)-coord(3,countstart+1))
!!!          !           !
!!!          !           r_j=n_vec*cos(2.0_dblprec*pi*qr+phi(iq))+s(:,iq)*sin(2.0_dblprec*pi*qr+phi(iq))
!!!          !           m_j=m_j+r_j
!!!          !           !
!!!          !        end do
!!!          !        call normalize(m_j)
!!!          !        m_j=m_j*mmom(ja,k)
!!!          !        emomM(1:3,ja,k)=m_j
!!!          !        !
!!!          !     end do
!!!                ! Calculate energy for given q,s,theta combination
!!!                call effective_field(Natom,Mensemble,countstart+1,countstart+na,         &
!!!                   do_jtensor,do_anisotropy,exc_inter,do_dm,do_pd,do_biqdm,do_bq,do_chir,&
!!!                   do_dip,emomM,mmom,external_field,time_external_field,beff,beff1,      &
!!!                   beff2,OPT_flag,max_no_constellations,maxNoConstl,unitCellType,        &
!!!                   constlNCoup,constellations,constellationsNeighType,mult_axis,         &
!!!                   energy,Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
!!!                ! Anisotropy + external field to be added
!!! 
!!! !           end do
!!!          end do
!!!          energy=energy/mry*mub/NA/4
!!!          write(2000,'(i8,g20.8,g20.8)') iter,q(1,1),energy
!!!          ! Store best energy configuration if trial energy is lower than minimum
!!!          if(energy<min_energy) then
!!!             do iq=1,nq
!!!                write(*,'(i9,a,3f12.6,a,f8.3,a,f8.3,a,f8.3,a,2g14.7)',advance="yes")   &
!!!                   iter,'  New Q: ',q(:,iq),'  Theta: ',theta(iq)/pi*180, '  Phi: ',phi(iq)/pi*180, &
!!!                   ' Cone angle: ',cone_ang,'  dE: ',energy-min_energy, energy
!!!             end do
!!! 
!!!             min_energy=energy
!!!             q_best=q
!!!             theta_best=theta
!!!             theta_glob_best=theta_glob
!!!             phi_best=phi
!!!             lhit=iter
!!!             nhits=nhits+1
!!! 
!!!             write(200,'(i8,f18.8,3f14.6,20f14.6)') iter,energy,q, theta/pi*180,theta_glob/pi*180
!!! 
!!!             ! Restore previous configuration
!!!             q=q_save
!!!             s=s_save
!!!             theta=theta_save
!!!             theta_glob=theta_glob_save
!!!             phi=phi_save
!!!             cone_ang=cone_ang_save
!!! 
!!!          end if
!!!          ! Reduce range for global search
!!!          q_range=q_range0*(1.0_dblprec/iter**0.5)
!!!          theta_glob_range=theta_glob_range0*(1.0_dblprec/iter**0.5)
!!!       end do
!!!       !
!!!       !
!!!       print '(1x,a,i6,a)','Stochastic minimization done with ',nhits,' hits.'
!!!       print '(1x,a)','--------Energy---------|------------Q-and-S-vectors------------------|-------Theta-----'
!!!       do iq=1,nq
!!!          print '(2x,f18.10,2x,3f14.6,2x,2f18.6)',min_energy,q_best(:,iq), theta_best(iq)/pi*180,phi_best(iq)/pi*180
!!!          print '(2x,f18.10,2x,3f14.6,2x,2f18.6)',min_energy,s(:,iq), cone_ang     ,phi_best(iq)/pi*180
!!!          print '(1x,a)','-----------------------|---------------------------------------------|-----------------'
!!!       end do
!!!       q=q_best
!!!       !q=q_norm*0.05000_dblprec
!!!       print '(1x,a)','-----------------------|---------------------------------------------|-----------------'
!!!       !
!!!       !
!!!       return
!!!       !
!!!    end subroutine mini_q
!!! 
!!!    !-----------------------------------------------------------------------------
!!!    ! SUBROUTINE: sweep_q
!!!    !> @brief Stupid line search minimization of spin spirals
!!!    !> @author Anders Bergman
!!!    !-----------------------------------------------------------------------------
!!!    subroutine sweep_q(Natom,Mensemble,NA,coord,do_jtensor,exc_inter,do_dm,do_pd,     &
!!!       do_biqdm,do_bq,do_chir,taniso,sb,do_dip,emomM,mmom,hfield,OPT_flag,           &
!!!       max_no_constellations,maxNoConstl,unitCellType,constlNCoup,constellations,    &
!!!       constellationsNeighType,mult_axis,Num_macro,cell_index,emomM_macro,           &
!!!       macro_nlistsize,do_anisotropy,simid)
!!! 
!!!       use RandomNumbers, only: rng_uniform,rng_gaussian, use_vsl
!!!       use Constants, only : mub, mry
!!!       use InputData, only : N1,N2,N3
!!!       use AMS, only : wrap_coord_diff
!!!       !
!!!       !.. Implicit declarations
!!!       implicit none
!!! 
!!!       integer, intent(in) :: Natom !< Number of atoms in system
!!!       integer, intent(in) :: Mensemble !< Number of ensembles
!!!       integer, intent(in) :: NA  !< Number of atoms in one cell
!!!       real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
!!!       integer, intent(in) :: do_jtensor   !<  Use SKKR style exchange tensor (0=off, 1=on, 2=with biquadratic exchange)
!!!       character(len=1),intent(in) :: exc_inter !< Interpolate Jij (Y/N)
!!!       integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_biqdm   !< Add Biquadratic DM (BIQDM) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_chir  !< Add scalar chirality exchange (CHIR) term to Hamiltonian (0/1)
!!!       integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
!!!       integer, intent(in) :: do_anisotropy
!!!       integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
!!!       real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between the anisotropies
!!!       real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
!!!       real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
!!!       real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
!!!       character(len=1), intent(in) :: mult_axis
!!!       !! +++ New variables due to optimization routines +++ !!
!!!       integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
!!!       ! Number of entries (number of unit cells*number of atoms per unit cell) in the constellation matrix per ensemble
!!!       integer, dimension(Mensemble), intent(in) :: maxNoConstl
!!!       ! See OptimizationRoutines.f90 for details on classification
!!!       integer, dimension(Natom, Mensemble), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
!!!       ! Matrix relating the interatomic exchanges for each atom in the constellation matrix
!!!       real(dblprec), dimension(ham%max_no_neigh, max_no_constellations,Mensemble), intent(in) :: constlNCoup
!!!       ! Matrix storing all unit cells belonging to any constellation
!!!       real(dblprec), dimension(3,max_no_constellations, Mensemble), intent(in) :: constellations
!!!       ! Optimization flag (1 = optimization on; 0 = optimization off)
!!!       logical, intent(in) :: OPT_flag
!!!       ! Matrix storing the type of the neighbours within a given neighbourhood of a constellation; default is 1 outside the neighbourhood region
!!!       ! The default is to achieve correct indexing. Note here also that constlNCoup will result in a net zero contribution to the Heissenberg exchange term
!!!       integer, dimension(ham%max_no_neigh,max_no_constellations,Mensemble), intent(in) :: constellationsNeighType
!!!       ! Internal effective field arising from the optimization of the Heissenberg exchange term
!!!       integer, intent(in) :: Num_macro !< Number of macrocells in the system
!!!       integer, dimension(Natom), intent(in) :: cell_index !< Macrocell index for each atom
!!!       integer, dimension(Num_macro), intent(in) :: macro_nlistsize !< Number of atoms per macrocell
!!!       real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
!!!       character(len=8), intent(in) :: simid  !< Name of simulation
!!!       !
!!!       integer :: iq
!!!       !
!!!       real(dblprec), dimension(3) :: m_i, m_j, r_i, r_j, m_avg
!!!       real(dblprec) :: pi, qr, rx,ry, rz
!!!       integer :: i,j, k, iia, ia, ja, lhit, nhits, countstart
!!!       real(dblprec) :: energy, min_energy
!!!       character(len=30) :: filn
!!!       !
!!!       real(dblprec), dimension(3,3) :: q_best, q_diff, s_diff, q_0, s_0, q_norm
!!!       real(dblprec), dimension(3,3) :: s_save, q_save
!!!       real(dblprec), dimension(3) :: theta_save, theta_diff, theta_best
!!!       real(dblprec) :: theta_glob_save, cone_ang, cone_ang_save
!!!       real(dblprec) :: theta_glob_diff, theta_glob_best
!!!       real(dblprec), dimension(3) :: phi_save, phi_diff, phi_best
!!!       real(dblprec), dimension(1) :: cone_ang_diff
!!!       real(dblprec), dimension(1) :: rng_tmp_arr
!!!       integer :: iter, iscale, i1 ,i2 ,i3
!!!       real(dblprec) :: q_range, theta_range, s_range, phi_range, theta_glob_range, cone_ang_range
!!!       real(dblprec) :: q_range0, theta_range0, s_range0, phi_range0, theta_glob_range0
!!!       real(dblprec), dimension(3) :: srvec 
!!!       !
!!!       !
!!!       real(dblprec) :: q_start=1.0_dblprec
!!!       !
!!!       !
!!!       pi=4._dblprec*ATAN(1._dblprec)
!!!       theta_glob_best=0.0_dblprec
!!!       theta_best=0.0_dblprec
!!!       phi_best=0.0_dblprec
!!!       !
!!!       ! Normal vector
!!!       ! Read from file or default
!!!       !n_vec(1)=0.0_dblprec;n_vec(2)=0.0_dblprec;n_vec(3)=1.0_dblprec;
!!!       !
!!!       ! Starting atom
!!!       I1 = N1/2
!!!       I2 = N2/2
!!!       I3 = N3/2
!!! 
!!!       countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
!!!       !
!!!       write(filn,'(''qm_sweep.'',a,''.out'')') trim(simid)
!!!       open(ofileno,file=filn, position="append")
!!!       write(ofileno,'(a)') "#    Iter                          Q-vector                                 Energy  "
!!! 
!!!       write(filn,'(''qm_minima.'',a,''.out'')') trim(simid)
!!!       open(ofileno2,file=filn, position="append")
!!!       write(ofileno2,'(a)') "#    Iter                          Q-vector                                 Energy  "
!!!       ! Read from ip_mcnstep
!!!       !niter=10000
!!!       ! Hard wired to planar spin spirals with 360/nq degree angle inbetween
!!!       !
!!!       theta = 0.0_dblprec
!!!       theta_glob = 0.0_dblprec
!!!       !theta_glob = pi/4.0_dblprec
!!!       phi = 0.0_dblprec
!!!       ! NQ=1 for now
!!!       do iq=1,nq
!!! 
!!!          ! Q read from file or default
!!!          !q(1,iq)=1.0_dblprec;q(2,iq)=0.0_dblprec;q(3,iq)=0.0_dblprec
!!!          !! Start along [100] direction
!!!          call normalize(q(1:3,iq))
!!!          !
!!!          q_norm(:,iq)=q(:,iq) !/sqrt(sum(q(:,iq)*q(:,iq))+1.0e-12_dblprec)
!!!          !
!!!          if(norm2(s(:,iq))==0.0_dblprec) then
!!!             ! Create pitch vectors perpendicular to q and n
!!!             s(1,iq)=q_norm(2,iq)*n_vec(3)-q_norm(3,iq)*n_vec(2)
!!!             s(2,iq)=q_norm(3,iq)*n_vec(1)-q_norm(1,iq)*n_vec(3)
!!!             s(3,iq)=q_norm(1,iq)*n_vec(2)-q_norm(2,iq)*n_vec(1)
!!!          end if
!!!          call normalize(s(1:3,iq))
!!!          !
!!!          print *,'----Q-and-S-vectors----',iq
!!!          print '(3f10.4)', q(:,iq)
!!!          print '(3f10.4)', s(:,iq)
!!!          print '(3f10.4)', n_vec(:)
!!!          !
!!!          ! Currently hard wired starting guess (q_start from file or default)
!!!          !q(:,iq)=q(:,iq)/sqrt(sum(q(:,iq)*q(:,iq))) *q_start  !/30.0_dblprec
!!!       end do
!!!       !
!!!       !For Neel spirals:
!!!       !s=q_norm
!!!       !do iq=1,nq
!!!       !   print *,'----Q_0 vector----',iq
!!!       !   print '(3f10.4)', q_0(:,iq)
!!!       !end do
!!!       q_0=q
!!!       s_0=s
!!!       q_scale=(q_max-q_min)/(1.0_dblprec*nstep)
!!!       min_energy=1.0d4
!!!       !
!!!       lhit=0
!!!       nhits=0
!!!       iscale=1
!!! 
!!!       ! Switch rotation direction
!!! 
!!!       ! Calculate total external field (not included yet)
!!!       do k=1,Mensemble
!!!          do i=1,Natom
!!!             external_field(1:3,i,k)= hfield
!!!             beff(1:3,i,k)=0.0_dblprec
!!!             beff1(1:3,i,k)=0.0_dblprec
!!!             beff2(1:3,i,k)=0.0_dblprec
!!!          end do
!!!       end do
!!! 
!!!       ! Only use first ensemble
!!!       k=1
!!!       do iter=0,nstep
!!!          !print *,'------------------'
!!!          ! Loop over q
!!!          do iq=1,nq
!!!             !q(:,iq)=q_0(:,iq)-q_0(:,iq)/sqrt(sum(q_0(:,iq)**2))*(2.0_dblprec*q_start*iter)/(1.0_dblprec*niter)
!!!             q(:,iq)=q_norm(:,iq)*(q_min+q_scale*iter)
!!!          end do
!!!          !  print '(a,i7,3f12.6)', '--->',iter, q(:,1)
!!!          !         !
!!!          ! Set up spin-spiral magnetization (only first cell+ neighbours)
!!!          energy=0.0_dblprec
!!!          !!!print *,'----Q-and-S-vectors----',1
!!!          !!!print '(3f10.4)', q(:,1)
!!!          !!!print '(3f10.4)', s(:,1)
!!!          !!!print '(3f10.4)', n_vec(:)
!!!          !!!!stop
!!!          do ia=1,Natom
!!!             !
!!!             !srvec=coord(:,ia)-coord(:,countstart+1)
!!!             ! Possible use wrap_coord_diff() here.
!!!             call wrap_coord_diff(Natom,coord,ia,countstart+1,srvec)
!!!             !
!!!             m_j=0.0_dblprec
!!!             do iq=1,nq
!!!                qr=q(1,iq)*srvec(1)+q(2,iq)*srvec(2)+q(3,iq)*srvec(3)
!!!                m_j=m_j+n_vec*cos(2*pi*qr+phi(iq))+s(:,iq)*sin(2*pi*qr+phi(iq))
!!!             end do
!!!             call normalize(m_j)
!!!             !emom(1:3,ia,k)=m_j
!!!             emomM(1:3,ia,k)=m_j*mmom(ia,k)
!!!             !print '(i7,3f12.6)', ia, emomM(1:3,ia,k)
!!! 
!!!          end do
!!!          ! Calculate energy for given q,s,theta combination
!!!          ! Anisotropy + external field to be added
!!!          call effective_field(Natom,Mensemble,countstart+1,countstart+na,         &
!!!             do_jtensor,do_anisotropy,exc_inter,do_dm,do_pd,do_biqdm,do_bq,do_chir,&
!!!             do_dip,emomM,mmom,external_field,time_external_field,beff,beff1,      &
!!!             beff2,OPT_flag,max_no_constellations,maxNoConstl,unitCellType,        &
!!!             constlNCoup,constellations,constellationsNeighType,mult_axis,         &
!!!             energy,Num_macro,cell_index,emomM_macro,macro_nlistsize,NA,N1,N2,N3)
!!! 
!!!          energy=energy/NA !/mry*mub/NA
!!! 
!!!          write(ofileno,'(i8,3g20.8,g20.8)') iter,q(:,1),energy
!!!          ! Store best energy configuration if trial energy is lower than minimum
!!!          if(energy<min_energy) then
!!!             !do iq=1,nq
!!!             !   write(*,'(i9,a,3f12.6,a,f8.3,a,f8.3,a,f8.3,a,2g14.7)',advance="yes")   &
!!!             !      iter,'  New Q: ',q(:,iq),'  Theta: ',theta(iq)/pi*180, '  Phi: ',phi(iq)/pi*180, &
!!!             !      ' Cone angle: ',cone_ang,'  dE: ',energy-min_energy, energy
!!!             !end do
!!! 
!!!             min_energy=energy
!!!             q_best=q
!!!             lhit=iter
!!!             nhits=nhits+1
!!! 
!!!             write(ofileno2,'(i8,3g20.8,g20.8)') iter,q(:,1),energy
!!! 
!!!          end if
!!!       end do
!!!       !
!!!       !
!!!       print '(1x,a,i6,a)','Line search minimization done with ',nhits,' hits.'
!!!       print '(1x,a)', '|-----Minimum energy----|----------------Q-vector-----------------|------------------S-vector----------------|'
!!!       do iq=1,nq
!!!          print '(2x,f18.10,2x,3f14.6,2x,3f14.6)',min_energy,q_best(:,iq),s(:,iq)
!!!       end do
!!!       ! Important: Save the lowest energy q-vector
!!!       q=q_best
!!!       print '(1x,a)','|-----------------------|-----------------------------------------|------------------------------------------|'
!!! 
!!!       !
!!!       close(ofileno)
!!!       close(ofileno2)
!!!       !
!!!       !
!!!       return
!!!       !
!!!    end subroutine sweep_q

!!!    !> Anisotropy
!!!    subroutine spinspiral_ani_field(Natom,Mensemble,NA,mmom,taniso,sb,hfield,energy)
!!!       !
!!!       implicit none
!!!       !
!!!       integer, intent(in) :: Natom !< Number of atoms in system
!!!       integer, intent(in) :: Mensemble !< Number of ensembles
!!!       integer, intent(in) :: NA  !< Number of atoms in one cell
!!!       real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
!!!       integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
!!!       real(dblprec), dimension(Natom),intent(in) :: sb !< Ratio between the anisotropies
!!!       real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
!!!       real(dblprec), intent(inout), optional :: energy !< Total energy
!!!       !
!!!       !
!!!       integer :: iq,i
!!!       real(dblprec) :: tt1, tt2, tt3, totmom, qfac
!!!       real(dblprec), dimension(3) :: field
!!!       !
!!!       totmom=0.0_dblprec
!!!       do i=1,NA
!!!          do iq=1,nq
!!!             field=0.0_dblprec
!!!             if (taniso(i)==1.or.taniso(i)==7) then
!!!                ! Uniaxial anisotropy
!!!                tt1=s(1,iq)*ham%eaniso(1,i)+s(2,iq)*ham%eaniso(2,i)+s(3,iq)*ham%eaniso(3,i)
!!!                tt1=tt1*mmom(i,1)*theta(iq)
!!! 
!!!                tt2=ham%kaniso(1,i)+2.0_dblprec*ham%kaniso(2,i)*(1-tt1*tt1)
!!!                !
!!!                tt3= 2.0_dblprec*tt1*tt2
!!! 
!!!                field(1)  = field(1) - tt3*ham%eaniso(1,i)
!!!                field(2)  = field(2) - tt3*ham%eaniso(2,i)
!!!                field(3)  = field(3) - tt3*ham%eaniso(3,i)
!!! 
!!!             end if
!!!             if (ham%taniso(i)==2.or.ham%taniso(i)==7) then
!!!                qfac=1.00
!!!                if(ham%taniso(i)==7) qfac=sb(i)
!!!                ! Cubic anisotropy
!!!                field(1) = field(1)  &
!!!                   + qfac*2.0_dblprec*ham%kaniso(1,i)*mmom(i,1)*s(1,iq)*(mmom(i,1)*s(2,iq)**2+mmom(i,1)*s(3,iq)**2)*theta(iq)**3 &
!!!                   + qfac*2.0_dblprec*ham%kaniso(2,i)*mmom(i,1)*s(1,iq)*mmom(i,1)*s(2,iq)**2*mmom(i,1)*s(3,iq)**2*theta(iq)**5
!!!                field(2) = field(2)  &
!!!                   + qfac*2.0_dblprec*ham%kaniso(1,i)*mmom(i,1)*s(2,iq)*(mmom(i,1)*s(3,iq)**2+mmom(i,1)*s(1,iq)**2) *theta(iq)**3&
!!!                   + qfac*2.0_dblprec*ham%kaniso(2,i)*mmom(i,1)*s(2,iq)*mmom(i,1)*s(3,iq)**2*mmom(i,1)*s(1,iq)**2*theta(iq)**5
!!!                field(3) = field(3)  &
!!!                   + qfac*2.0_dblprec*ham%kaniso(1,i)*mmom(i,1)*s(3,iq)*(mmom(i,1)*s(1,iq)**2+mmom(i,1)*s(2,iq)**2) *theta(iq)**3&
!!!                   + qfac*2.0_dblprec*ham%kaniso(2,i)*mmom(i,1)*s(3,iq)*mmom(i,1)*s(1,iq)**2*mmom(i,1)*s(2,iq)**2*theta(iq)**5
!!!                !
!!!             end if
!!!             energy=energy-(2.0_dblprec*field(1)*mmom(i,1)*theta(iq)*s(1,iq)+2.0_dblprec*field(1)*mmom(2,1)*theta(iq)*s(2,iq)+2.0_dblprec*field(3)*mmom(i,1)*theta(iq)*s(3,iq))
!!!             energy=energy+(hfield(1)*mmom(i,1)*theta(iq)*s(1,iq)+hfield(2)*mmom(i,1)*theta(iq)*s(2,iq)+hfield(3)*mmom(i,1)*theta(iq)*s(3,iq))/nq*0.5_dblprec
!!! 
!!!          end do
!!!       end do
!!!       return
!!!       !
!!!    end subroutine spinspiral_ani_field
