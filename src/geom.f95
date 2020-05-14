
! ===================================================================
MODULE geom_data
use iso_fortran_env, only: int32, real64
!constructing problem geometry, x,y vectors, lement-node array.


IMPLICIT NONE
INTEGER, PUBLIC, PARAMETER :: maxreg=10

REAL (kind=real64), PUBLIC, ALLOCATABLE :: xv(:) ! nodal x points node order
REAL (kind=real64), PUBLIC, ALLOCATABLE :: yv(:) ! nodal y points node order

!     .. Module Scalars ..
integer, public    :: nreg
INTEGER, PUBLIC    :: nel    !max no. elements
INTEGER, PUBLIC    :: nnod  !max no. nodes

! "logicals" - actually integers
INTEGER, PUBLIC, ALLOCATABLE :: znodbound(:)

!       ..arrays..
INTEGER, DIMENSION(:,:),PUBLIC, ALLOCATABLE :: nodelist !element-node array

INTEGER, PUBLIC    :: ireg    ! region number
INTEGER, PUBLIC    :: inod    !node number
INTEGER, PUBLIC    :: iel    !element number
INTEGER :: i,j

END MODULE geom_data
