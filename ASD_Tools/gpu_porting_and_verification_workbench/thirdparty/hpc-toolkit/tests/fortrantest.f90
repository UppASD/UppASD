! SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
!
! SPDX-License-Identifier: Apache-2.0

program fortrantest
    use iso_c_binding, only: c_size_t, c_double
    implicit none

    interface
        subroutine reduce_sum(count, data, result) bind(C, name="reduce_sum")
            import :: c_size_t, c_double
            integer(c_size_t), value    :: count
            real(c_double), intent(in)  :: data(*)
            real(c_double), intent(out) :: result
        end subroutine reduce_sum
    end interface

    integer(c_size_t), parameter :: n = 1000
    real(c_double), allocatable  :: array(:)
    real(c_double)               :: total_c, total_fortran

    allocate(array(n))

    call random_number(array)

    call reduce_sum(n, array, total_c)
    total_fortran = sum(array)

    print *, "reduce_sum (C):", total_c
    print *, "sum (Fortran): ", total_fortran

    deallocate(array)
end program fortrantest
