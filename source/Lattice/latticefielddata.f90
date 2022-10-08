!> Data for storing external and effective fields
module LatticeFieldData
   use Parameters
   use Profiling
   !
   implicit none
   !
   !real(dblprec), dimension(:,:), allocatable :: sitefld !< Site dependent applied efield
   real(dblprec), dimension(:,:,:), allocatable :: eeff !< Total effective efield from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: eeff1 !< Internal effective efield from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: eeff2 !< External efield from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: eeff3 !< Internal effective efield from mixed spin-lattice Hamiltonians
   real(dblprec), dimension(:,:,:), allocatable :: e2eff !< Next step total effective efield from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: latt_external_field !< External static electric field
   real(dblprec), dimension(:,:,:), allocatable :: latt_time_external_field !< External time-dependent electric field
   real(dblprec), dimension(:,:,:), allocatable :: ethermal_field


contains


   !> Allocate/deallocate arrays for external and effective fields
   subroutine latt_allocate_fields(Natom,Mensemble,flag)

      implicit none

      integer, intent(in),optional :: Natom !< Number of atoms in system
      integer, intent(in),optional :: Mensemble !< Number of ensembles 
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      ! If flag > 0 allocate arrays
      if(flag>0) then
         allocate(eeff(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(eeff))*kind(eeff),'eeff','latt_allocate_fields')
         eeff=0.0_dblprec
         allocate(eeff1(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(eeff1))*kind(eeff1),'eeff1','latt_allocate_fields')
         eeff1=0.0_dblprec
         allocate(eeff2(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(eeff2))*kind(eeff2),'eeff2','latt_allocate_fields')
         eeff2=0.0_dblprec
         allocate(eeff3(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(eeff3))*kind(eeff3),'eeff3','latt_allocate_fields')
         eeff3=0.0_dblprec
         allocate(e2eff(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(e2eff))*kind(e2eff),'e2eff','latt_allocate_fields')
         e2eff=0.0_dblprec
         allocate(latt_external_field(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(latt_external_field))*kind(latt_external_field),'latt_external_field','latt_allocate_fields')
         latt_external_field=0.0_dblprec
         allocate(latt_time_external_field(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(latt_time_external_field))*kind(latt_time_external_field),'latt_time_external_field','latt_allocate_fields')
         latt_time_external_field=0.0_dblprec
         allocate(ethermal_field(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ethermal_field))*kind(ethermal_field),'ethermal_field','latt_allocate_fields')
         ethermal_field=0.0_dblprec
      else
         ! Otherwise deallocate arrays
         i_all=-product(shape(eeff))*kind(eeff)
         deallocate(eeff,stat=i_stat)
         call memocc(i_stat,i_all,'eeff','latt_allocate_fields')
         i_all=-product(shape(eeff1))*kind(eeff1)
         deallocate(eeff1,stat=i_stat)
         call memocc(i_stat,i_all,'eeff1','latt_allocate_fields')
         i_all=-product(shape(eeff2))*kind(eeff2)
         deallocate(eeff2,stat=i_stat)
         call memocc(i_stat,i_all,'eeff2','latt_allocate_fields')
         i_all=-product(shape(eeff3))*kind(eeff3)
         deallocate(eeff3,stat=i_stat)
         call memocc(i_stat,i_all,'eeff3','latt_allocate_fields')
         i_all=-product(shape(e2eff))*kind(e2eff)
         deallocate(e2eff,stat=i_stat)
         call memocc(i_stat,i_all,'e2eff','latt_allocate_fields')
         i_all=-product(shape(latt_external_field))*kind(latt_external_field)
         deallocate(latt_external_field,stat=i_stat)
         call memocc(i_stat,i_all,'latt_external_field','latt_allocate_fields')

         if (allocated(latt_time_external_field)) then
            i_all=-product(shape(latt_time_external_field))*kind(latt_time_external_field)
            deallocate(latt_time_external_field,stat=i_stat)
            call memocc(i_stat,i_all,'latt_time_external_field','latt_allocate_fields')
         endif

         i_all=-product(shape(ethermal_field))*kind(ethermal_field)
         deallocate(ethermal_field,stat=i_stat)
         call memocc(i_stat,i_all,'ethermal_field','latt_allocate_fields')

      end if
   end subroutine latt_allocate_fields


end module LatticeFieldData
