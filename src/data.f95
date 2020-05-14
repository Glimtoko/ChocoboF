module geom_data
use iso_fortran_env, only: int32, real64
use globalconstants
implicit none

integer, public, parameter :: maxreg=10

end module geom_data


module mesh_mod
    use iso_fortran_env, only: int32, real64
    type :: mesh
        integer(kind=int32) :: nreg
        integer(kind=int32) :: nel    !max no. elements
        integer(kind=int32) :: nnod  !max no. nodes

        ! "logicals" - actually integers
        integer(kind=int32), allocatable :: znodbound(:)

        !       ..arrays..
        integer(kind=int32), dimension(:,:), allocatable :: nodelist !element-node array
        real(kind=real64), dimension(:), allocatable :: xv ! nodal x points node order
        real(kind=real64), dimension(:), allocatable :: yv ! nodal y points node order
        real(kind=real64), dimension(:), allocatable :: pre !element pressures
        real(kind=real64), dimension(:), allocatable :: rho !element density
        real(kind=real64), dimension(:), allocatable :: en !element energy
        real(kind=real64), dimension(:), allocatable :: cc !element sound speed
        real(kind=real64), dimension(:), allocatable :: qq !element artificial viscosity
        real(kind=real64), dimension(:), allocatable :: massel !element mass
        real(kind=real64), dimension(:), allocatable :: volel !element volume
        real(kind=real64), dimension(:), allocatable :: volelold !element old vol
        real(kind=real64), dimension(:), allocatable :: area  !element area
        real(kind=real64), dimension(:), allocatable :: uv !node x velocity
        real(kind=real64), dimension(:), allocatable :: vv!node y velocity

        real(kind=real64), dimension(:), allocatable :: uvold !old x velocity
        real(kind=real64), dimension(:), allocatable :: vvold !old y velocity
        real(kind=real64), dimension(:), allocatable :: uvbar !av x velocity
        real(kind=real64), dimension(:), allocatable :: vvbar !av y velocity
        real(kind=real64), dimension(:), allocatable :: xv05 !1/2dt node x
        real(kind=real64), dimension(:), allocatable :: yv05 !1/2dt node y
        real(kind=real64), dimension(:), allocatable :: pre05 !1/2dt pressures
        real(kind=real64), dimension(:), allocatable :: rho05 !1/2dt density
        real(kind=real64), dimension(:), allocatable :: en05 !1/2dt energy
        real(kind=real64), dimension(:), allocatable :: volel05 !1/2dt volume
        real(kind=real64), dimension(:), allocatable :: divint !int divergence v
        real(kind=real64), dimension(:), allocatable :: divvel !int divergence v
        real(kind=real64), dimension(:), allocatable :: dxtb,dxlr
        real(kind=real64), dimension(:), allocatable :: dudxt,dudxb,dudxl,dudxr
        real(kind=real64), dimension(:), allocatable :: phib, phit, phil, phir

        ! used in momentum and hourgl
        real(kind=real64), dimension(:), allocatable ::forcenodx !force mass x
        real(kind=real64), dimension(:), allocatable ::forcenody !force mass y

        !fem
        real(kind=real64), dimension(:,:), allocatable ::nint
        real(kind=real64), dimension(:,:), allocatable ::elwtc
        real(kind=real64), dimension(:,:), allocatable ::dndx
        real(kind=real64), dimension(:,:), allocatable ::dndy
        real(kind=real64), dimension(:,:), allocatable ::pdndx
        real(kind=real64), dimension(:,:), allocatable ::pdndy
    end type
end module


module mesh_data
USE GlobalConstants
!     .. module scalars ..
integer, public    :: nreg
integer, public    :: nel    !max no. elements
integer, public    :: nnod  !max no. nodes

! "logicals" - actually integers
integer, public, allocatable :: znodbound(:)

!       ..arrays..
integer, dimension(:,:), public, allocatable :: nodelist !element-node array
real(kind=real64), public, allocatable :: xv(:) ! nodal x points node order
real(kind=real64), public, allocatable :: yv(:) ! nodal y points node order
real(kind=real64), public, allocatable ::pre(:) !element pressures
real(kind=real64), public, allocatable ::rho(:)!element density
real(kind=real64), public, allocatable ::en(:)!element energy
real(kind=real64), public, allocatable ::cc(:) !element sound speed
real(kind=real64), public, allocatable ::qq(:) !element artificial viscosity
real(kind=real64), public, allocatable ::massel(:) !element mass
real(kind=real64), public, allocatable ::volel(:) !element volume
real(kind=real64), public, allocatable ::volelold(:) !element old vol
real(kind=real64), public, allocatable ::area(:) !element area
real(kind=real64),dimension(:),public, allocatable ::uv !node x velocity
real(kind=real64),dimension(:),public, allocatable ::vv!node y velocity

real(kind=real64), public, allocatable ::uvold(:) !old x velocity
real(kind=real64), public, allocatable ::vvold(:) !old y velocity
real(kind=real64), public, allocatable ::uvbar(:) !av x velocity
real(kind=real64), public, allocatable ::vvbar(:) !av y velocity
real(kind=real64), public, allocatable ::xv05(:) !1/2dt node x
real(kind=real64), public, allocatable ::yv05(:) !1/2dt node y
real(kind=real64), public, allocatable ::pre05(:) !1/2dt pressures
real(kind=real64), public, allocatable ::rho05(:) !1/2dt density
real(kind=real64), public, allocatable ::en05(:)!1/2dt energy
real(kind=real64), public, allocatable ::volel05(:) !1/2dt volume
real(kind=real64), public, allocatable ::divint(:) !int divergence v
real(kind=real64), public, allocatable ::divvel(:) !int divergence v
real(kind=real64), public, allocatable ::dxtb(:),dxlr(:)
real(kind=real64), public,allocatable ::dudxt(:),dudxb(:),dudxl(:),dudxr(:)
real(kind=real64), public,allocatable ::phib(:), phit(:), phil(:), phir(:)

! used in momentum and hourgl
real(kind=real64), public, allocatable ::forcenodx(:) !force mass x
real(kind=real64), public, allocatable ::forcenody(:) !force mass y

!fem
real(kind=real64), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::nint(:,:)
real(kind=real64), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::elwtc(:,:)
real(kind=real64), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::dndx(:,:)
real(kind=real64), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::dndy(:,:)
real(kind=real64), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::pdndx(:,:)
real(kind=real64), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::pdndy(:,:)
end module mesh_data


module cutoffs
use globalConstants
! cut off parameters
REAL(kind=real64), PUBLIC, PARAMETER :: mindt=1.00000000000e-4
REAL(kind=real64), PUBLIC, PARAMETER ::  zcut=1.00000000000e-8
REAL(kind=real64), PUBLIC, PARAMETER :: dtminhg=0.000000008
REAL(kind=real64), PUBLIC, PARAMETER ::  zerocut=1.0000000e-18
! stop divide by 0's
REAL(kind=real64), PUBLIC, PARAMETER :: dencut=1.000000000e-6
end module cutoffs


module core_input
use globalConstants
use geom_data
REAL(kind=real64), PUBLIC:: gamma
! artificial viscosity
REAL(kind=real64), PUBLIC:: cq
REAL(kind=real64), PUBLIC:: cl
!max time step
REAL(kind=real64), PUBLIC:: maxallstep
! initial time step
REAL(kind=real64), PUBLIC:: dtinit
INTEGER, PUBLIC :: dtoption ! timestep option
! time step groth factor
REAL(kind=real64), PUBLIC :: t0  !initial time
REAL(kind=real64), PUBLIC :: tf  !final time
REAL(kind=real64), PUBLIC :: growth
! if 0 cartesian 1 axisymmetric
INTEGER, PUBLIC :: zaxis
INTEGER, PUBLIC :: zintdivvol   ! 0 not vol 1 indiv volume change
INTEGER, PUBLIC :: avtype ! type of artificial viscosity
INTEGER, PUBLIC :: zantihg  ! if 0 no hourglassing if 1 hourglassing
INTEGER, PUBLIC :: hgregtyp(1:maxreg)
REAL(kind=real64), PUBLIC :: kappareg(1:maxreg)

real(kind=real64) :: dtsilo, lastsilo
INTEGER:: nadvect, stepcnt, h5type
logical :: tioonefile

NAMELIST /tinp/t0,tf,gamma,cq,cl,maxallstep,dtinit,dtoption,growth,zaxis,  &
        zintdivvol,avtype,zantihg,hgregtyp,kappareg,stepcnt,dtsilo,h5type,tioonefile

end module core_input
