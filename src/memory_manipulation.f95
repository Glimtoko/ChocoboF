module memory_management
use iso_fortran_env, only: int32, real64

interface set_data
    module procedure set_data_int1d
    module procedure set_data_int2d
    module procedure set_data_real1d
    module procedure set_data_real2d
end interface

contains


subroutine set_data_int1d(array, n, start, default)
    integer(kind=int32), dimension(:), allocatable, intent(out) :: array
    integer(kind=int32), intent(in) :: n
    integer(kind=int32), intent(in), optional :: start
    integer(kind=int32), intent(in), optional :: default

    integer(kind=int32) :: d, s

    if (present(start)) then
        s = start
    else
        s = 1
    end if

    if (present(default)) then
        d = default
    else
        d = 0_int32
    end if

    if (allocated(array)) deallocate(array)
    allocate(array(s:n))
    array = d
end subroutine


subroutine set_data_int2d(array, n1, n2, start, default)
    integer(kind=int32), dimension(:,:), allocatable, intent(out) :: array
    integer(kind=int32), intent(in) :: n1
    integer(kind=int32), intent(in) :: n2
    integer(kind=int32), intent(in), optional :: start
    integer(kind=int32), intent(in), optional :: default

    integer(kind=int32) :: d, s

    if (present(start)) then
        s = start
    else
        s = 1
    end if

    if (present(default)) then
        d = default
    else
        d = 0_int32
    end if

    if (allocated(array)) deallocate(array)
    allocate(array(s:n1, s:n2))
    array = d
end subroutine


subroutine set_data_real1d(array, n, start, default)
    real(kind=real64), dimension(:), allocatable, intent(out) :: array
    integer(kind=int32), intent(in) :: n
    integer(kind=int32), intent(in), optional :: start
    real(kind=real64), intent(in), optional :: default

    integer(kind=int32) :: d, s

    if (present(start)) then
        s = start
    else
        s = 1
    end if

    if (present(default)) then
        d = default
    else
        d = 0.0_real64
    end if

    if (allocated(array)) deallocate(array)
    allocate(array(s:n))
    array = d
end subroutine


subroutine set_data_real2d(array, n1, n2, start, default)
    real(kind=real64), dimension(:,:), allocatable, intent(out) :: array
    integer(kind=int32), intent(in) :: n1
    integer(kind=int32), intent(in) :: n2
    integer(kind=int32), intent(in), optional :: start
    real(kind=real64), intent(in), optional :: default

    integer(kind=int32) :: d, s

    if (present(start)) then
        s = start
    else
        s = 1
    end if

    if (present(default)) then
        d = default
    else
        d = 0.0_real64
    end if

    if (allocated(array)) deallocate(array)
    allocate(array(s:n1, s:n2))
    array = d
end subroutine


end module
