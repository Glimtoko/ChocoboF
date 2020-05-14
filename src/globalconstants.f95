MODULE GLOBALCONSTANTS
use iso_fortran_env, only: int32, real64
! defines double precision
IMPLICIT NONE

! INTEGER, PARAMETER        :: DP=selected_real_kind(15,300)

REAL (kind=real64), PARAMETER :: zero = 0.00000000000000000_real64,             &
                              one = 1.00000000000000000_real64,             &
                              two = 2.00000000000000000_real64,             &
                            three = 3.00000000000000000_real64,             &
                             four = 4.00000000000000000_real64,             &
                             five = 5.00000000000000000_real64,             &
                              six = 6.00000000000000000_real64,             &
                            seven = 7.00000000000000000_real64,             &
                            eight = 8.00000000000000000_real64,             &
                             nine = 9.00000000000000000_real64,             &
                              ten = 10.0000000000000000_real64,             &
                           twelve = 12.0000000000000000_real64,             &
                           thirty = 30.0000000000000000_real64,             &
                       thirtyfive = 35.0000000000000000_real64,             &
                             half = 5.00000000000000000e-1_real64,          &
                          quarter = 2.50000000000000000e-1_real64,          &
                               pi = 3.14159265358979324_real64,             &
                            twopi = two*pi
COMPLEX (kind=real64), PARAMETER :: im =(0.00000000000000000_real64,             &
                                    1.00000000000000000_real64)
END MODULE GLOBALCONSTANTS
