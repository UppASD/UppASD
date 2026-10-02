! SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
!
! SPDX-License-Identifier: Apache-2.0

program fortrantest_device
    use iso_c_binding, only: c_size_t, c_int64_t, c_ptr, c_double, c_f_pointer
    implicit none

    type, bind(C) :: buffer_device_double
        integer(c_size_t) :: count
        type(c_ptr)        :: data
    end type buffer_device_double

    type, bind(C) :: device_rng_t
        integer(c_size_t) :: count
        type(c_ptr)        :: states
    end type device_rng_t

    type, bind(C) :: device_reduce_workspace_t
        integer(c_size_t) :: bytes
        type(c_ptr)        :: data
    end type device_reduce_workspace_t

    type, bind(C) :: buffer_pinned_double
        integer(c_size_t) :: count
        type(c_ptr)        :: data
    end type buffer_pinned_double

    interface
        function buffer_create_device_double(count) &
                bind(C, name="buffer_create_device_double") result(buffer)
            import :: c_size_t, buffer_device_double
            integer(c_size_t), value   :: count
            type(buffer_device_double) :: buffer
        end function buffer_create_device_double

        subroutine buffer_destroy_device_double(buffer) &
                bind(C, name="buffer_destroy_device_double")
            import :: buffer_device_double
            type(buffer_device_double), intent(inout) :: buffer
        end subroutine buffer_destroy_device_double

        function device_rng_create(seed, count) bind(C, name="device_rng_create") result(rng)
            import :: c_int64_t, c_size_t, device_rng_t
            integer(c_int64_t), value :: seed
            integer(c_size_t), value  :: count
            type(device_rng_t)        :: rng
        end function device_rng_create

        subroutine device_rng_destroy(rng) bind(C, name="device_rng_destroy")
            import :: device_rng_t
            type(device_rng_t), intent(inout) :: rng
        end subroutine device_rng_destroy

        subroutine device_randomize(count, data, rng) bind(C, name="device_randomize")
            import :: c_size_t, c_ptr, device_rng_t
            integer(c_size_t), value          :: count
            type(c_ptr), value                :: data
            type(device_rng_t), intent(inout) :: rng
        end subroutine device_randomize

        function device_reduce_workspace_required_bytes(count, data, res) &
                bind(C, name="device_reduce_workspace_required_bytes") result(required_bytes)
            import :: c_size_t, c_ptr
            integer(c_size_t), value :: count
            type(c_ptr), value       :: data
            type(c_ptr), value       :: res
            integer(c_size_t)        :: required_bytes
        end function device_reduce_workspace_required_bytes

        function device_reduce_workspace_create(bytes) &
                bind(C, name="device_reduce_workspace_create") result(workspace)
            import :: c_size_t, device_reduce_workspace_t
            integer(c_size_t), value        :: bytes
            type(device_reduce_workspace_t) :: workspace
        end function device_reduce_workspace_create

        subroutine device_reduce_workspace_destroy(workspace) &
                bind(C, name="device_reduce_workspace_destroy")
            import :: device_reduce_workspace_t
            type(device_reduce_workspace_t), intent(inout) :: workspace
        end subroutine device_reduce_workspace_destroy

        subroutine device_reduce_sum(count, data, res, workspace) &
                bind(C, name="device_reduce_sum")
            import :: c_size_t, c_ptr, device_reduce_workspace_t
            integer(c_size_t), value                       :: count
            type(c_ptr), value                              :: data
            type(c_ptr), value                              :: res
            type(device_reduce_workspace_t), intent(inout) :: workspace
        end subroutine device_reduce_sum

        function buffer_create_pinned_double(count) &
                bind(C, name="buffer_create_pinned_double") result(buffer)
            import :: c_size_t, buffer_pinned_double
            integer(c_size_t), value   :: count
            type(buffer_pinned_double) :: buffer
        end function buffer_create_pinned_double

        subroutine buffer_destroy_pinned_double(buffer) &
                bind(C, name="buffer_destroy_pinned_double")
            import :: buffer_pinned_double
            type(buffer_pinned_double), intent(inout) :: buffer
        end subroutine buffer_destroy_pinned_double

        subroutine buffer_copy_d2h_double(src, dst) bind(C, name="buffer_copy_d2h_double")
            import :: buffer_device_double, buffer_pinned_double
            type(buffer_device_double), intent(in)    :: src
            type(buffer_pinned_double), intent(inout) :: dst
        end subroutine buffer_copy_d2h_double
    end interface

    integer(c_size_t), parameter  :: n    = 1000
    integer(c_int64_t), parameter :: seed = 1763249876_c_int64_t

    type(buffer_device_double)      :: data_buf, result_buf
    type(device_rng_t)              :: rng
    type(device_reduce_workspace_t) :: workspace
    type(buffer_pinned_double)      :: host_result
    integer(c_size_t)               :: required_bytes
    real(c_double), pointer         :: total_device

    data_buf = buffer_create_device_double(n)
    rng = device_rng_create(seed, data_buf%count)
    call device_randomize(data_buf%count, data_buf%data, rng)

    result_buf = buffer_create_device_double(1_c_size_t)
    required_bytes = device_reduce_workspace_required_bytes(data_buf%count, data_buf%data, &
                                                            result_buf%data)
    workspace = device_reduce_workspace_create(required_bytes)
    call device_reduce_sum(data_buf%count, data_buf%data, result_buf%data, workspace)

    host_result = buffer_create_pinned_double(result_buf%count)
    call buffer_copy_d2h_double(result_buf, host_result)
    call c_f_pointer(host_result%data, total_device)

    print *, "device_reduce_sum:", total_device

    call buffer_destroy_pinned_double(host_result)
    call device_reduce_workspace_destroy(workspace)
    call buffer_destroy_device_double(result_buf)
    call device_rng_destroy(rng)
    call buffer_destroy_device_double(data_buf)
end program fortrantest_device
