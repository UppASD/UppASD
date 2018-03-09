!====================================================================!
!> @brief
!> Intended to use for calculation of the stiffness constant for FM materials
!
!> @author
!> Jonathan Chico
!> M. Pereiro ---> Improvment of the regresion method by using Pade aproximants
!> L. Bergqvist -> Added Tc-MFA
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!====================================================================!

! Routines to calculate the ferromagnetic exchange stiffness for ordered systems and random alloys
! The formalism used here is the one used by Pajda et al. PRB 64, 174402, where the exchange interactions are summed
! and weighted by the square of the distance to obtain the spin wave stiffness D. This formalism is then exted to multi-sublattice
! cases and random alloys to be able to explore more complex situations
module stiffness

   use Constants
   use Profiling

   implicit none

   integer :: eta_max !< Number of convergence parameters for the stiffness
   integer :: eta_min !< Minimum convergence parameters for the stiffness
   character(len=1) :: do_stiffness !< Calculate the spin wave stiffness of a ferromagnet
   character(len=1) :: do_dm_stiffness !< Calculate the DMI spiralization
   character(len=1) :: prn_J0_matrix   !< Print the full site dependent J0 matrix
   real(dblprec) :: A_xc     !< Exchange stiffness constant from rational fit J/m
   real(dblprec) :: M_sat    !< Saturation magnetization muB/m^3
   real(dblprec) :: Dxc_fit  !< Spin wave stiffnes from rational fit
   real(dblprec) :: cell_vol !< Unit cell volume m^3
   real(dblprec) :: A_xc_lsq !< Exchange stiffness constant from LSQ fit
   real(dblprec) :: D_err_fit   !< Error form spin wave stiffnes rational fit
   real(dblprec) :: Dxc_fit_lsq !< Spin wave stiffness LSQ fit
   real(dblprec), dimension(:,:), allocatable :: J0_matrix !< Matrix being used to calculate the Tc-MFA
   real(dblprec), dimension(:,:), allocatable :: D_xc_stiffness_matrix !< Spin wave stiffness tensor rational fit
   real(dblprec), dimension(:,:), allocatable :: A_xc_stiffness_matrix !< Exchange stiffness constant rational fit
   real(dblprec), dimension(:,:), allocatable :: D_xc_stiffness_matrix_lsq !< Spin wave stiffness tensor LSQ fit
   real(dblprec), dimension(:,:), allocatable :: A_xc_stiffness_matrix_lsq !< Exchange stiffness constant LSQ fit
   real(dblprec), dimension(:,:), allocatable :: DM0_mat !< DMI spiralization matrix [meVA]
   real(dblprec), dimension(:,:), allocatable :: DM0_mat_lsq !< DMI spiralization matrix (LSQ fit) [meVA]

   public


contains

!---------------------------------------------------------------------------
!> @brief
!> Calculate the exchange stiffness and the spin wave stiffness for a ferromagnet
!
!> @author
!> Jonathan Chico
!>
!> @date 10/02/2017 - Jonathan Chico
!> - Added the calculation of the stiffness tensor, that is now also a matrix
!>   will be written which can be useful for non cubic systems
!---------------------------------------------------------------------------
  subroutine ferro_stiffness(NT,NA,N1,N2,N3,hdim,mconf,Natom,Nchmax,eta_max,eta_min,&
             conf_num,do_ralloy,Natom_full,max_no_neigh,Nch,anumb,atype,nlistsize,&
             atype_ch,asite_ch,achem_ch,nlist,alat,C1,C2,C3,coord,chconc,ammom_inp,ncoup,&
             A_xc,M_sat,Dxc_fit,cell_vol,A_xc_lsq,D_err_fit,Dxc_fit_lsq,J0_matrix,&
             D_xc_stiffness_matrix,A_xc_stiffness_matrix,D_xc_stiffness_matrix_lsq,&
             A_xc_stiffness_matrix_lsq)

    use Constants

    implicit none

    ! .. Input variables
    integer, intent(in) :: NT  !< Number of types of atoms
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    integer, intent(in) :: hdim   !< Number of elements in Hamiltonian element (scalar or vector)
    integer, intent(in) :: mconf  !< LSF ground state configuration
    integer, intent(in) :: Natom  !< Number of atoms in system
    integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
    integer, intent(in) :: eta_max  !< Number of convergence parameters for the stiffness
    integer, intent(in) :: eta_min  !< Minimum  convergence parameters for the stiffness
    integer, intent(in) :: conf_num !< Number of LSF configurations
    integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
    integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, dimension(NA), intent(in) :: Nch !< Number of chemical components on each site in cell
    integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist  !< Neighbour list for Heisenberg exchange couplings

    real(dblprec), intent(in) :: alat !< Lattice parameter
    real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
    real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
    real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
    real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
    real(dblprec), dimension(NA,Nchmax), intent(in) :: chconc !< Chemical concentration on sites
    real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp !< Magnetic moment directions from input (for alloys)
    real(dblprec), dimension(hdim,max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings

    ! .. Output variables
    real(dblprec), intent(out) :: A_xc     !< Exchange stiffness constant from rational fit J/m
    real(dblprec), intent(out) :: M_sat    !< Saturation magnetization muB/m^3
    real(dblprec), intent(out) :: Dxc_fit  !< Spin wave stiffnes from rational fit
    real(dblprec), intent(out) :: cell_vol !< Unit cell volume m^3
    real(dblprec), intent(out) :: A_xc_lsq !< Exchange stiffness constant from LSQ fit
    real(dblprec), intent(out) :: D_err_fit   !< Error form spin wave stiffnes rational fit
    real(dblprec), intent(out) :: Dxc_fit_lsq !< Spin wave stiffness LSQ fit
    real(dblprec), dimension(NA,NA), intent(out) :: J0_matrix !< Matrix being used to calculate the Tc-MFA
    real(dblprec), dimension(3,3), intent(out) :: D_xc_stiffness_matrix !< Spin wave stiffness tensor rational fit
    real(dblprec), dimension(3,3), intent(out) :: A_xc_stiffness_matrix !< Exchange stiffness constant rational fit
    real(dblprec), dimension(3,3), intent(out) :: D_xc_stiffness_matrix_lsq !< Spin wave stiffness tensor LSQ fit
    real(dblprec), dimension(3,3), intent(out) :: A_xc_stiffness_matrix_lsq !< Exchange stiffness constant LSQ fit

    ! .. Local variables
    integer :: ii, jj
    integer :: lwork, info, alwork
    integer :: iatom, jatom
    integer :: katom, I1, I2, I3, countstart
    integer :: i, k, eta, ich, eta_redu

    real(dblprec) :: total_mom
    real(dblprec) :: fcinv, jij, jijsign
    real(dblprec) :: rij2, rij, rfit, drfit, stiff_par

    real(dblprec), dimension(3) :: rcoord
    real(dblprec), dimension(eta_max) :: temp_x, eig_val, eig_val_temp
    real(dblprec), dimension(eta_max-(eta_min-1),3) :: lmatrix
    real(dblprec), dimension(eta_max-(eta_min-1)) :: dvector
    real(dblprec), dimension(eta_max,3,3) :: eig_val_mat

    real(dblprec), dimension(:), allocatable :: wres
    real(dblprec), dimension(:), allocatable :: awork
    real(dblprec), dimension(:,:,:), allocatable :: etemp
    real(dblprec), dimension(:,:,:), allocatable :: stiff_matrix !< Matrix being used to calculate the ferromagnetic exchange stiffness
    real(dblprec), dimension(:,:,:,:,:), allocatable :: D_matrix !< Matrix used to calculate the FM exchange stiffness tensor

    real(dblprec), dimension(:), allocatable :: work
    real(dblprec), dimension(:), allocatable :: rwork
    real(dblprec), dimension(:), allocatable :: iwres
    real(dblprec), dimension(:), allocatable :: iwres_mat
    real(dblprec), dimension(:), allocatable :: rwres
    real(dblprec), dimension(:), allocatable :: rwres_mat
    real(dblprec), dimension(:), allocatable :: cwres
    real(dblprec), dimension(:), allocatable :: cwres_mat
    real(dblprec), dimension(:,:),allocatable :: ctemp
    real(dblprec), dimension(:,:),allocatable :: A_inplace
    real(dblprec), dimension(:,:),allocatable :: A_mat_inplace

    lwork=16*na
    alwork=((eta_max-(eta_min-1))*16)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate work arrays
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(cwres(NA))
    allocate(rwres(NA))
    allocate(iwres(NA))
    allocate(work(lwork))
    allocate(cwres_mat(NA))
    allocate(rwres_mat(NA))
    allocate(iwres_mat(NA))
    allocate(ctemp(na,na))
    allocate(rwork(lwork))
    allocate(wres(eta_max))
    allocate(awork(alwork))
    allocate(A_inplace(NA,NA))
    allocate(A_mat_inplace(NA,NA))
    allocate(etemp(NA,NA,eta_max))
    allocate(stiff_matrix(NA,NA,eta_max))
    allocate(D_matrix(3,3,NA,NA,eta_max))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialization of variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    work=0.0d0
    wres=0.0d0
    etemp=0.0d0
    ctemp=0.0d0
    rwork=0.0d0
    cwres=0.0d0
    iwres=0.0d0
    rwres=0.0d0
    temp_x=0.0d0
    eig_val=0.0d0
    lmatrix=0.0d0
    dvector=0.0d0
    cwres_mat=0.0d0
    iwres_mat=0.0d0
    rwres_mat=0.0d0
    total_mom=0.0d0
    J0_matrix=0.0d0
    A_inplace=0.0d0
    eig_val_mat=0.0d0
    eig_val_temp=0.0d0
    stiff_matrix=0.0d0
    A_mat_inplace=0.0d0
    D_matrix=0.0d0
    D_xc_stiffness_matrix=0.0d0
    A_xc_stiffness_matrix=0.0d0
    D_xc_stiffness_matrix_lsq=0.0d0
    A_xc_stiffness_matrix_lsq=0.0d0
    D_err_fit=0.0d0

    fcinv=mub/mry
    I1 = N1/2
    I2 = N2/2
    I3 = N3/2
    countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate the volume of the cell in meters
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cell_vol=(C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+&
          C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+&
          C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1))*alat**3
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating the total moment of the cell
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (do_ralloy==0) then
       do i=1,NA
          total_mom=total_mom+ammom_inp(i,1,mconf)
       enddo
    ! In case of random alloy calculate the weighted average
    else if (do_ralloy==1) then
       do i=1,NA
          do ich=1,Nch(i)
             total_mom=total_mom+ammom_inp(i,ich,mconf)*chconc(i,ich)
           enddo
       enddo
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating the saturation magnetization
    M_sat=total_mom/cell_vol
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate the stiffness for checmically pure systems
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (do_ralloy==0) then

       ! Need to create a matrix which includes inter and intra sublattice interactions
       ! Now must loop over the convergency factor to make sure that the sum is well defined
       do eta=1, eta_max

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Loop over atoms in the unit cell and sum up the exchange interactions
          do i=1, NA
             iatom=i+countstart
             ! Loop over the neighbors for each atom in the unit cell
             do jatom=1, nlistsize(iatom)
                ! Neighbouring atom
                katom=nlist(jatom,iatom)

                ! Distace vector between the atoms
                rcoord(1)=coord(1,katom)-coord(1,iatom)
                rcoord(2)=coord(2,katom)-coord(2,iatom)
                rcoord(3)=coord(3,katom)-coord(3,iatom)

                ! Distance vector between neighbouring atoms
                rij2=rcoord(1)**2+rcoord(2)**2+rcoord(3)**2
                ! Distance between neighbouring atoms
                rij=sqrt(rij2)

                ! Calculating the "real" Jijs in mRyd
                jij=ncoup(1,jatom,iatom,mconf)*fcinv*ammom_inp(anumb(katom),1,mconf)*ammom_inp(anumb(iatom),1,mconf)*0.5d0
                ! Which is the sign between the magnetic moments (To take care of AFM interactions)
                jijsign=sign(ammom_inp(anumb(katom),1,mconf),ammom_inp(anumb(iatom),1,mconf))/abs(ammom_inp(anumb(katom),1,mconf))
                ! Calculate the proportionality parameter
                stiff_par=(2.0d0/3.0d0)*jij*jijsign*rij2*(alat**2)/(sqrt(abs(ammom_inp(anumb(katom),1,mconf))*abs(ammom_inp(anumb(iatom),1,mconf))))
                ! The actual stiffness matrix
                stiff_matrix(i,anumb(katom),eta)=stiff_matrix(i,anumb(katom),eta)+stiff_par*exp(-0.10d0*eta*rij)

                ! Calculate the stiffness matrix
                do ii=1,3
                   do jj=1,3
                      ! Create a parameter for the calculation of the stiffness
                      stiff_par=2.0d0*jij*jijsign*rcoord(ii)*rcoord(jj)*(alat**2)/&
                      !stiff_par=2.0d0*jij*jijsign*rcoord(ii)*rcoord(jj)/&
                      (sqrt(abs(ammom_inp(anumb(katom),1,mconf))*abs(ammom_inp(anumb(iatom),1,mconf))))
                      ! Save the parameter for the stiffness matrix
                      D_matrix(ii,jj,i,anumb(katom),eta)=D_matrix(ii,jj,i,anumb(katom),eta)+&
                      stiff_par*exp(-0.10d0*eta*sqrt(abs(rcoord(ii)*rcoord(jj))))
                   enddo
                enddo
                J0_matrix(i,anumb(katom))=J0_matrix(i,anumb(katom))+(2.0d0/3.0d0)*jij*jijsign
             enddo

          enddo

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating the eigenvalues for the scalar stiffness
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! The matrix is transformed to the more familiar units of meVÅ^2
          A_inplace(1:NA,1:NA)=stiff_matrix(1:NA,1:NA,eta)*1d20*ry_ev

          ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
          call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,eta), NA, WORK, LWORK, INFO)
          if(info.ne.0) then
              print '(2x,a,i4)', 'Problem in zgeev 1:',info
          end if

          ! Temporal x-axis for the fitting to a polynomial
          temp_x(eta)=0.10d0*eta
          ! Finding the maximum eigenvalue of the exchange stiffness matrix
          eig_val(eta)=maxval((rwres))
          !eig_val(eta)=maxval(real(cwres))
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! End of calculation of eignevalues for the scalar stiffness
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating the eigenvalues for the stiffness tensor
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           do ii=1,3
              do jj=1,3
                 ! Transform to meVA^2
                 A_mat_inplace(1:NA,1:NA)=D_matrix(ii,jj,1:NA,1:NA,eta)*1d20*ry_ev
                 ! Calculate eigenvalues for each component of the matrix
                 call dgeev('N','N',NA, A_mat_inplace, NA, rwres_mat, iwres_mat, ctemp, NA, etemp(1,1,eta), NA, WORK, LWORK, INFO)
                 if(info.ne.0) then
                     print '(2x,a,i4)', 'Problem in zgeev 2:',info
                 end if
                 eig_val_mat(eta,ii,jj)=maxval((rwres_mat))
              enddo
           enddo
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! End of calculation of eignevalues for the stiffness tensor
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       enddo

       ! Defining the reduced eta, common for scalar and tensor
       eta_redu=eta_max-(eta_min-1)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculating the scalar stiffness with rational polynomials
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! A polynomial fit is performed to obtain the spin wave stiffness
       call ratint(temp_x(eta_min:eta_max),eig_val(eta_min:eta_max),eta_redu,0.0d0,rfit,drfit)
       ! Calculate the exchange stiffness from micromagnetics with rational
       ! polynomials
       A_xc=rfit*M_sat*1e-20/(1000*Joule_ev*4.0d0)
       ! Storing the spin wave stiffness scalar
       Dxc_fit=rfit
       D_err_fit=drfit

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculating the scalar stiffness with rational polynomials
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculating the stiffness tensor with rational polynomials
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do ii=1,3
          do jj=1,3
             eig_val_temp(1:eta_max)=eig_val_mat(1:eta_max,ii,jj)
             ! A polynomial fit is performed to obtain the spin wave stiffness
              call ratint(temp_x(eta_min:eta_max),eig_val_temp(eta_min:eta_max),eta_redu,0.0d0,rfit,drfit)
             ! Calculate the exchange stiffness from micromagnetics with rational
             ! polynomials
             A_xc_stiffness_matrix(ii,jj)=rfit*1e-20*M_sat/(1000*Joule_ev*4.0d0)
             ! Saving the spin wave stiffness matrix
             D_xc_stiffness_matrix(ii,jj)=rfit
         enddo
       enddo
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculating the stiffness tensor with rational polynomials
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculating the scalar exchange stiffness with LSQ
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Creating a vector for the eta parameter
       do eta=eta_min,eta_max
           k=eta-(eta_min-1)
           lmatrix(k,1)=1.d0 ; lmatrix(k,2)=0.1d0*eta ; lmatrix(k,3)=(0.1d0*eta)**2; dvector(k)=eig_val(eta)
       enddo

       call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
       if(info.ne.0) then
           print '(2x,a,i4)', 'Problem in dgels:',info
       end if
       ! Calculate the Exchange stiffness from micromagnetics
       A_xc_lsq=dvector(1)*1e-20*M_sat/(1000*Joule_ev*4.0d0)

       Dxc_fit_lsq=dvector(1)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculation of the scalar exchange stiffness with LSQ
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculating the exchange stiffness tensor with LSQ
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do ii=1, 3
          do jj=1,3
          ! Creating a vector for the eta parameter
          do eta=eta_min,eta_max
              k=eta-(eta_min-1)
              lmatrix(k,1)=1.d0 ; lmatrix(k,2)=0.1d0*eta ; lmatrix(k,3)=(0.1d0*eta)**2; dvector(k)=eig_val_mat(eta,ii,jj)
          enddo

          call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
          if(info.ne.0) then
              print '(2x,a,i4)', 'Problem in dgels:',info
          end if
          ! Calculate the exchange stiffness from micromagnetics
          A_xc_stiffness_matrix_lsq(ii,jj)=dvector(1)*1e-20*M_sat/(1000*Joule_ev*4.0d0)

          D_xc_stiffness_matrix_lsq(ii,jj)=dvector(1)

         enddo
       enddo

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculation of the exchange stiffness tensor with LSQ
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculation of the MF Tc
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       J0_matrix=J0_matrix/eta_max
       A_inplace(1:NA,1:NA)=(J0_matrix*ry_ev*1e-3)/k_bolt_ev
       ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
       call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,1), NA, WORK, LWORK, INFO)
       if(info.ne.0) then
           print '(2x,a,i4)', 'Problem in zgeev 3:',info
       end if
       eig_val(1)=maxval((rwres))
       !eig_val(1)=maxval(real(cwres))
       write(*,1009) 'Tc-MFA from stiffness :',eig_val(1),'K'

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculation of the MF Tc
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate the stiffness for the random alloy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else if (do_ralloy==1) then
      ! Loop over the convergency factor to make sure that the sum is well defined

      do eta=1, eta_max

         ! Loop over atoms the atoms in the unit cell
         do i=1,NA
             iatom=i+countstart
             ! Loop over the neighbors for each atom in the unit cell
             do jatom=1, nlistsize(iatom)

                ! Neighbouring atom
                katom=nlist(jatom,iatom)

                ! Distace vector between the atoms
                rcoord(1)=coord(1,katom)-coord(1,iatom)
                rcoord(2)=coord(2,katom)-coord(2,iatom)
                rcoord(3)=coord(3,katom)-coord(3,iatom)

                ! Distance vector between neighbouring atoms
                rij2=rcoord(1)**2+rcoord(2)**2+rcoord(3)**2
                ! Distance between neighbouring atoms
                rij=sqrt(rij2)

                ! Calculating the "real" Jijs in mRyd
                jij=ncoup(1,jatom,iatom,mconf)*fcinv*ammom_inp(asite_ch(katom),achem_ch(katom),mconf)*&
                    ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf)*0.5d0
                ! Which is the sign between the magnetic moments (To take care of AFM interactions)
                jijsign=sign(ammom_inp(asite_ch(katom),achem_ch(katom),mconf),ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf))/&
                        abs(ammom_inp(asite_ch(katom),achem_ch(katom),mconf))

                ! Parameter for the stiffness
                stiff_par=(2.0d0/3.0d0)*jij*jijsign*rij2*(alat**2)/(sqrt(abs(ammom_inp(asite_ch(katom),achem_ch(katom),mconf))*&
                           abs(ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf))))
                ! The actual stiffness matrix
                stiff_matrix(i,asite_ch(katom),eta)=stiff_matrix(i,asite_ch(katom),eta)+&
                stiff_par*exp(-0.10d0*eta*rij)

                ! Calculate the stiffness matrix
                do ii=1,3
                   do jj=1,3
                      D_matrix(ii,jj,i,anumb(katom),eta)=D_matrix(ii,jj,i,asite_ch(katom),eta)+&
                      2.0d0*jij*jijsign*rcoord(ii)*rcoord(jj)*(alat**2)*exp(-0.10d0*eta*rij)/&
                      (sqrt(abs(ammom_inp(asite_ch(katom),achem_ch(katom),mconf))*abs(ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf))))
                   enddo
                enddo

                ! Filling the J0 matrix
                J0_matrix(i,asite_ch(katom))=J0_matrix(i,asite_ch(katom))+(2.0d0/3.0d0)*jij*jijsign
             enddo

         enddo

        ! Calculation of the eigenvalues for the stiffness matrix

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating the eigenvalues for the scalar stiffness
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! The matrix is transformed to the more familiar units of meVÅ^2
          A_inplace(1:NA,1:NA)=stiff_matrix(1:NA,1:NA,eta)*1d20*ry_ev

          ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
          call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,eta), NA, WORK, LWORK, INFO)
          if(info.ne.0) then
              print '(2x,a,i4)', 'Problem in zgeev 4:',info
          end if

          ! Temporal x-axis for the fitting to a polynomial
          temp_x(eta)=0.10d0*eta
          ! Finding the maximum eigenvalue of the exchange stiffness matrix
          eig_val(eta)=maxval((rwres))
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! End of calculation of eignevalues for the scalar stiffness
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating the eigenvalues for the stiffness tensor
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           do ii=1,3
              do jj=1,3
                 ! Transform to meVA^2
                 A_mat_inplace(1:NA,1:NA)=D_matrix(ii,jj,1:NA,1:NA,eta)*1d20*ry_ev
                 ! Calculate eigenvalues for each component of the matrix
                 call dgeev('N','N',NA, A_mat_inplace, NA, rwres_mat, iwres_mat, ctemp, NA, etemp(1,1,eta), NA, WORK, LWORK, INFO)
                 eig_val_mat(eta,ii,jj)=maxval((rwres_mat))
              enddo
           enddo
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! End of calculation of eignevalues for the stiffness tensor
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo

       ! Defining the reduced eta, common for scalar and tensor
       eta_redu=eta_max-(eta_min-1)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculating the scalar stiffness with rational polynomials
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! A polynomial fit is performed to obtain the spin wave stiffness
       call ratint(temp_x(eta_min:eta_max),eig_val(eta_min:eta_max),eta_redu,0.0d0,rfit,drfit)
       ! Calculate the exchange stiffness from micromagnetics with rational
       ! polynomials
       A_xc=rfit*M_sat*1e-20/(1000*Joule_ev*4.0d0)
       ! Storing the spin wave stiffness scalar
       Dxc_fit=rfit
       D_err_fit=drfit

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculating the scalar stiffness with rational polynomials
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculating the stiffness tensor with rational polynomials
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do ii=1,3
          do jj=1,3
             eig_val_temp(1:eta_max)=eig_val_mat(1:eta_max,ii,jj)
             ! A polynomial fit is performed to obtain the spin wave stiffness
             call ratint(temp_x(eta_min:eta_max),eig_val_temp(eta_min:eta_max),eta_redu,0.0d0,rfit,drfit)
             ! Calculate the exchange stiffness from micromagnetics with rational
             ! polynomials
             A_xc_stiffness_matrix(ii,jj)=rfit*1e-20*M_sat/(1000*Joule_ev*4.0d0)
             ! Saving the spin wave stiffness matrix
             D_xc_stiffness_matrix(ii,jj)=rfit
         enddo
       enddo
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculating the stiffness tensor with rational polynomials
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculating the scalar exchange stiffness with LSQ
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Creating a vector for the eta parameter
       do eta=eta_min,eta_max
           k=eta-(eta_min-1)
           lmatrix(k,1)=1.d0 ; lmatrix(k,2)=0.1d0*eta ; lmatrix(k,3)=(0.1d0*eta)**2; dvector(k)=eig_val(eta)
       enddo

       call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
       if(info.ne.0) then
           print '(2x,a,i4)', 'Problem in dgels:',info
       end if
       ! Calculate the Exchange stiffness from micromagnetics
       A_xc_lsq=dvector(1)*1e-20*M_sat/(1000*Joule_ev*4)

       Dxc_fit_lsq=dvector(1)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculation of the scalar exchange stiffness with LSQ
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculating the exchange stiffness tensor with LSQ
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do ii=1, 3
          do jj=1,3
          ! Creating a vector for the eta parameter
          do eta=eta_min,eta_max
              k=eta-(eta_min-1)
              lmatrix(k,1)=1.d0 ; lmatrix(k,2)=0.1d0*eta ; lmatrix(k,3)=(0.1d0*eta)**2; dvector(k)=eig_val_mat(eta,ii,jj)
          enddo

          call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
          if(info.ne.0) then
              print '(2x,a,i4)', 'Problem in dgels:',info
          end if
          ! Calculate the Exchange stiffness from micromagnetics
          A_xc_stiffness_matrix_lsq(ii,jj)=dvector(1)*1e-20*M_sat/(1000*Joule_ev*4)

          D_xc_stiffness_matrix_lsq(ii,jj)=dvector(1)

         enddo
       enddo

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculation of the exchange stiffness tensor with LSQ
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculation of the MF Tc
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       J0_matrix=J0_matrix/eta_max
       A_inplace(1:NA,1:NA)=(J0_matrix*ry_ev*1e-3)/k_bolt_ev
       ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
       call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,eta), NA, WORK, LWORK, INFO)
       if(info.ne.0) then
           print '(2x,a,i4)', 'Problem in zgeev 5:',info
       end if
       eig_val(1)=maxval((rwres))
       write(*,1009) 'Tc-MFA from stiffness :',eig_val(1),'K'

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! End of calculation of the MF Tc
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! End of stiffness calculations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! The issue for the chemical dissorder is that there might not be different chemical types in the same cell
   deallocate(work)
   deallocate(wres)
   deallocate(cwres)
   deallocate(ctemp)
   deallocate(rwork)
   deallocate(awork)
   deallocate(etemp)
   deallocate(D_matrix)
   deallocate(A_inplace)
   deallocate(cwres_mat)
   deallocate(stiff_matrix)
   deallocate(A_mat_inplace)

1009 format (2x,a,f10.1,2x,a)
  end subroutine ferro_stiffness

!---------------------------------------------------------------------------
!> @brief
!> Calculation of the tensorial DMI stiffness
!
!> @author
!> Anders Bergman
!> @date 10/02/2017 - Jonathan Chico
!> - Routine originally written by Anders Bergman, implemented in ASD by
!> Jonathan Chico
!> - Modification to make compatible with new printing routine
!---------------------------------------------------------------------------
  subroutine DMI_stiffness(NT,NA,N1,N2,N3,Natom,Nchmax,eta_max,eta_min,max_no_dmneigh,anumb,&
             dmlistsize,dmlist,alat,coord,ammom_inp,dm_vect,DM0_mat,DM0_mat_lsq)

    use Constants

    implicit none

    ! .. Input variables
    integer, intent(in) :: NT  !< Number of types of atoms
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    integer, intent(in) :: Natom  !< Number of atoms in system
    integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
    integer, intent(in) :: eta_max  !< Number of convergence parameters for the stiffness
    integer, intent(in) :: eta_min  !< Minimum  convergence parameters for the stiffness
    integer, intent(in) :: max_no_dmneigh !< Calculated maximum of neighbours for exchange
    integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
    integer, dimension(Natom), intent(in) :: dmlistsize !< Size of neighbour list for DM
    integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist !< Neighbour list for Heisenberg exchange couplings

    real(dblprec), intent(in) :: alat !< Lattice parameter
    real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
    real(dblprec), dimension(NA,Nchmax), intent(in) :: ammom_inp !< Magnetic moment directions from input (for alloys)
    real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Heisenberg exchange couplings

    ! .. Output variables
    real(dblprec), dimension(3,3), intent(out) :: DM0_mat !< DMI spiralization matrix [meVA]
    real(dblprec), dimension(3,3), intent(out) :: DM0_mat_lsq !< DMI spiralization matrix [meVA]

    ! .. Local variables
    integer :: info,ii,jj
    integer :: i, j, k, lwork, alwork, eta
    integer :: iatom, jatom
    integer :: katom, I1, I2, I3, countstart
    real(dblprec) :: rij2, rij
    real(dblprec) :: fcinv, dmsign
    real(dblprec) :: dm_mag_par
    real(dblprec), dimension(500) :: lp_work
    real(dblprec), dimension(3) :: rcoord, DM_xc, dm_stiff,dij
    real(dblprec), dimension(:),allocatable :: eigenvals
    real(dblprec), dimension(:,:,:,:,:), allocatable :: dm_mat !< Matrix being used to calculate the DM stiffness
    complex(dblprec), dimension(:,:),allocatable :: eigenvecs
    real(dblprec), dimension(eta_max-(eta_min-1),3) :: lmatrix
    real(dblprec), dimension(eta_max-(eta_min-1)) :: dvector
    real(dblprec), dimension(:), allocatable :: awork
    real(dblprec), dimension(:,:,:), allocatable :: DM0_mat_eta

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate working arrays
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    alwork=((eta_max-(eta_min-1))*16)
    allocate(eigenvals(NA))
    allocate(eigenvecs(NA,NA))
    allocate(dm_mat(3,3,NA,NA,0:eta_max))
    allocate(DM0_mat_eta(3,3,0:eta_max))
    allocate(awork(alwork))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DM_xc=0.0d0
    rcoord=0.0d0
    dm_mat=0.0d0
    DM0_mat=0.0d0
    DM0_mat_eta=0.0d0
    dm_stiff=0.0d0
    eigenvals=0.0d0
    lmatrix=0.0d0
    dvector=0.0d0

    I1 = N1/2
    I2 = N2/2
    I3 = N3/2
    fcinv=mub/mry
    countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! This routine might be redundant
    if (anumb(countstart+1)/=1) then ! could be mod(countstart,NA) ==0 aswell
       do i = 1,NA
         if (anumb(countstart+i+1)==1) then
            countstart = countstart+i
            exit
         else if (anumb(countstart-i+1)==1) then
            countstart = countstart-i
            exit
         end if
       end do
    end if


    do eta=0,eta_max
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Loop over atoms in the unit cell and sum up the DMI vectors
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do i=1, NA
          iatom=i+countstart
          ! Loop over the neighbors for each atom in the unit cell
          do jatom=1, dmlistsize(iatom)

             ! Neighbouring atom
             katom=dmlist(jatom,iatom)
             ! Distance VECTOR betwwen neighbouring atoms
             rcoord(1)=coord(1,katom)-coord(1,iatom)
             rcoord(2)=coord(2,katom)-coord(2,iatom)
             rcoord(3)=coord(3,katom)-coord(3,iatom)

             ! Distance between the atoms
             rij2=rcoord(1)**2+rcoord(2)**2+rcoord(3)**2
             rij=sqrt(rij2)

             ! Calculating the "real" DMI in mRyd
             dij(1)=dm_vect(1,jatom,iatom)*fcinv*ammom_inp(anumb(katom),1)*ammom_inp(anumb(iatom),1)*0.5d0
             dij(2)=dm_vect(2,jatom,iatom)*fcinv*ammom_inp(anumb(katom),1)*ammom_inp(anumb(iatom),1)*0.5d0
             dij(3)=dm_vect(3,jatom,iatom)*fcinv*ammom_inp(anumb(katom),1)*ammom_inp(anumb(iatom),1)*0.5d0
             ! Which is the sign between the magnetic moments (To take care of AFM interactions)
             dmsign=sign(ammom_inp(anumb(katom),1),ammom_inp(anumb(iatom),1))/abs(ammom_inp(anumb(katom),1))

             ! Variable to weight in the moments and lattice parameter
             dm_mag_par=alat/sqrt(abs(ammom_inp(anumb(katom),1))*abs(ammom_inp(anumb(iatom),1)))

             ! The actual DM stiffness matrix in here it will be looped over spin
             ! and lattice space
             do ii=1,3
                do jj=1,3
                   ! Notice that this matrix mixes both lattice and spin degrees of
                   ! freedom as the dij acts on the spins and rij is a vector in
                   ! real space
                   !dm_mat(ii,jj,i,anumb(katom))=dm_mat(ii,jj,i,anumb(katom))+dmsign*rcoord(ii)*dij(jj)*&
                   !                             dm_mag_par
                   dm_mat(ii,jj,i,anumb(katom),eta)=dm_mat(ii,jj,i,anumb(katom),eta)+dmsign*rcoord(ii)*dij(jj)*&
                      dm_mag_par*exp(-0.10d0*eta*rij)
                enddo
             enddo

          enddo

       enddo
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of calculation of DM spiralization matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! ! Loop to calculate the eigenvalues of the DM spiralization matrix
    !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DM0_mat_eta=0.0d0
     do eta=0,eta_max
        do i=1,3
           do j=1,3
              do ii=1,na
                 do jj=1,na
                    DM0_mat_eta(i,j,eta)=DM0_mat_eta(i,j,eta)+dm_mat(i,j,ii,jj,eta)*1d10*ry_ev
                 end do
              end do
           enddo
        enddo
     end do
     DM0_mat=DM0_mat_eta(:,:,0)
    !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! ! End of calculation of the DM spiralzation eigenvalues
    !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LSQ over eta
    do i=1,3
       do j=1,3
          ! Creating a vector for the eta parameter
          do eta=eta_min,eta_max
             k=eta-(eta_min-1)
             lmatrix(k,1)=1.d0 ; lmatrix(k,2)=0.1d0*eta ; lmatrix(k,3)=(0.1d0*eta)**2; dvector(k)=DM0_mat_eta(i,j,eta)
          enddo

          call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
          if(info.ne.0) then
             print '(2x,a,i4)', 'Problem in dgels:',info
          end if
          DM0_mat_lsq(i,j)=dvector(1)
       enddo
    enddo

    deallocate(eigenvals)
    deallocate(eigenvecs)
    deallocate(dm_mat)
    deallocate(DM0_mat_eta)
    deallocate(awork)

  end subroutine DMI_stiffness

!---------------------------------------------------------------------------
!> @brief
!> Fitting via a rational polynomial approach
!
!> @author
!> Manuel Pereiro
!---------------------------------------------------------------------------
  subroutine ratint(xa,ya,n,x,y,dy)
   ! Largest expected value of n, and a small number.
   ! Given arrays xa and ya, each of length n, and given a value of x,
   ! this routine returns a value of y and an accuracy estimate dy. The
   ! value returned is that of the diagonal rational function, evaluated at x,
   ! which passes through the n points (xai, yai), i = 1...n.
   implicit none

   ! .. Input variables
   integer :: n
   real(dblprec), intent(in) :: xa(:),ya(:),x

   ! .. Output variables
   real(dblprec), intent(out) :: y,dy

   ! .. Local variables
   integer :: i,m,ns
   real :: dd,h,hh,t,w
   real(dblprec) :: TINY
   parameter (TINY=1.e-25)
   real(dblprec), dimension(n) :: c,d

   ns=1
   hh=abs(x-xa(1))

   do i=1,n
      h=abs(x-xa(i))
      if (h.eq.0.) then
          y=ya(i)
          dy=0.0
          return
       else if (h.lt.hh) then
          ns=i
          hh=h
       endif
       c(i)=ya(i)
       d(i)=ya(i)+TINY !The TINY part is needed to prevent a rare zero over zero condition.
   end do

   y=ya(ns)
   ns=ns-1

   do m=1,n-1
     do i=1,n-m
        w=c(i+1)-d(i)
        h=xa(i+m)-x !h will never be zero, since this was tested in the initialization loop
        w=c(i+1)-d(i)
        h=xa(i+m)-x
        t=(xa(i)-x)*d(i)/h
        dd=t-c(i+1)
        !if(dd.eq.0.) pause 'failure in ratint' !This error condition indicates that the interpolating
        if(dd.eq.0.) stop 'failure in ratint' !This error condition indicates that the interpolating
        dd=w/dd                                !function has a pole at the requested value of x.
        d(i)=c(i+1)*dd
        c(i)=t*dd
     end do

     if (2*ns.lt.n-m) then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif

     y=y+dy
   end do

end subroutine ratint

end module stiffness
