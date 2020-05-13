module geom_data
use globalconstants
implicit none

integer, public, parameter :: maxreg=10

real (kind=dp), public, allocatable :: xv(:) ! nodal x points node order
real (kind=dp), public, allocatable :: yv(:) ! nodal y points node order

!     .. module scalars ..
integer, public    :: nreg
integer, public    :: nel    !max no. elements
integer, public    :: nnod  !max no. nodes

! "logicals" - actually integers
integer, public, allocatable :: znodbound(:)

!       ..arrays..
integer, dimension(:,:),public, allocatable :: nodelist !element-node array

integer, public    :: ireg    ! region number
integer, public    :: inod    !node number
integer, public    :: iel    !element number
integer :: i,j

end module geom_data


module mesh_data
USE GlobalConstants
!        .. vectors ..
real (kind=dp), public, allocatable ::pre(:) !element pressures
real (kind=dp), public, allocatable ::rho(:)!element density
real (kind=dp), public, pointer ::en(:)!element energy
real (kind=dp), public, allocatable ::cc(:) !element sound speed
real (kind=dp), public, allocatable ::qq(:) !element artificial viscosity
real (kind=dp), public, allocatable ::massel(:) !element mass
real (kind=dp), public, allocatable ::volel(:) !element volume
real (kind=dp), public, allocatable ::volelold(:) !element old vol
real (kind=dp), public, allocatable ::area(:) !element area
real (kind=dp),dimension(:),public, pointer ::uv !node x velocity
real (kind=dp),dimension(:),public, pointer ::vv!node y velocity

real (kind=dp), public, allocatable ::uvold(:) !old x velocity
real (kind=dp), public, allocatable ::vvold(:) !old y velocity
real (kind=dp), public, allocatable ::uvbar(:) !av x velocity
real (kind=dp), public, allocatable ::vvbar(:) !av y velocity
real (kind=dp), public, allocatable ::xv05(:) !1/2dt node x
real (kind=dp), public, allocatable ::yv05(:) !1/2dt node y
real (kind=dp), public, allocatable ::pre05(:) !1/2dt pressures
real (kind=dp), public, allocatable ::rho05(:) !1/2dt density
real (kind=dp), public, pointer ::en05(:)!1/2dt energy
real (kind=dp), public, allocatable ::volel05(:) !1/2dt volume
real (kind=dp), public, allocatable ::divint(:) !int divergence v
real (kind=dp), public, allocatable ::divvel(:) !int divergence v
real (kind=dp), public, allocatable ::dxtb(:),dxlr(:)
real (kind=dp), public,allocatable ::dudxt(:),dudxb(:),dudxl(:),dudxr(:)
real (kind=dp), public,allocatable ::phib(:), phit(:), phil(:), phir(:)

! used in momentum and hourgl
real (kind=dp), public, allocatable ::forcenodx(:) !force mass x
real (kind=dp), public, allocatable ::forcenody(:) !force mass y

!fem
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::nint(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::elwtc(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::dndx(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::dndy(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::pdndx(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::pdndy(:,:)

REAL(KIND=DP), PUBLIC :: time
INTEGER, PUBLIC       ::prout, stepno
REAL(KIND=DP), PUBLIC :: dt  !time step size
end module mesh_data


module cutoffs
use globalConstants
! cut off parameters
REAL(KIND=DP), PUBLIC, PARAMETER :: mindt=1.00000000000e-4
REAL(KIND=DP), PUBLIC, PARAMETER ::  zcut=1.00000000000e-8
REAL(KIND=DP), PUBLIC, PARAMETER :: dtminhg=0.000000008
REAL(KIND=DP), PUBLIC, PARAMETER ::  zerocut=1.0000000e-18
! stop divide by 0's
REAL(KIND=DP), PUBLIC, PARAMETER :: dencut=1.000000000e-6
end module cutoffs


module core_input
use globalConstants
use geom_data
REAL(KIND=DP), PUBLIC:: gamma
! artificial viscosity
REAL(KIND=DP), PUBLIC:: cq
REAL(KIND=DP), PUBLIC:: cl
!max time step
REAL(KIND=DP), PUBLIC:: maxallstep
! initial time step
REAL(KIND=DP), PUBLIC:: dtinit
INTEGER, PUBLIC :: dtoption ! timestep option
! time step groth factor
REAL(KIND=DP), PUBLIC :: t0  !initial time
REAL(KIND=DP), PUBLIC :: tf  !final time
REAL(KIND=DP), PUBLIC :: growth
! if 0 cartesian 1 axisymmetric
INTEGER, PUBLIC :: zaxis
INTEGER, PUBLIC :: zintdivvol   ! 0 not vol 1 indiv volume change
INTEGER, PUBLIC :: avtype ! type of artificial viscosity
INTEGER, PUBLIC :: zantihg  ! if 0 no hourglassing if 1 hourglassing
INTEGER, PUBLIC :: hgregtyp(1:maxreg)
REAL(KIND=DP), PUBLIC :: kappareg(1:maxreg)

real(kind=dp) :: dtsilo, lastsilo
INTEGER:: nadvect, stepcnt, h5type
logical :: tioonefile

NAMELIST /tinp/t0,tf,gamma,cq,cl,maxallstep,dtinit,dtoption,growth,zaxis,  &
        zintdivvol,avtype,zantihg,hgregtyp,kappareg,stepcnt,dtsilo,h5type,tioonefile

end module core_input
