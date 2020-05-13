MODULE GLOBALCONSTANTS
! defines double precision
IMPLICIT NONE

INTEGER, PARAMETER        :: DP=selected_real_kind(15,300)

REAL (KIND=DP), PARAMETER :: zero = 0.00000000000000000_DP,             &
                              one = 1.00000000000000000_DP,             &
                              two = 2.00000000000000000_DP,             &
                            three = 3.00000000000000000_DP,             &
                             four = 4.00000000000000000_DP,             &
                             five = 5.00000000000000000_DP,             &
                              six = 6.00000000000000000_DP,             &
                            seven = 7.00000000000000000_DP,             &
                            eight = 8.00000000000000000_DP,             &
                             nine = 9.00000000000000000_DP,             &
                              ten = 10.0000000000000000_DP,             &
                           twelve = 12.0000000000000000_DP,             &
                           thirty = 30.0000000000000000_DP,             &
                       thirtyfive = 35.0000000000000000_DP,             &
                             half = 5.00000000000000000e-1_DP,          &
                          quarter = 2.50000000000000000e-1_DP,          &
                               pi = 3.14159265358979324_DP,             &
                            twopi = two*pi
COMPLEX (KIND=DP), PARAMETER :: im =(0.00000000000000000_DP,             &
                                    1.00000000000000000_DP)
END MODULE GLOBALCONSTANTS
