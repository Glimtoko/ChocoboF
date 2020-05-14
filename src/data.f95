module mesh_mod
use iso_fortran_env, only: int32, real64
type :: MeshT
    integer(kind=int32) :: nel  ! Number of cells
    integer(kind=int32) :: nnod ! Number of nodes
    integer(kind=int32) :: nreg ! Number of regions

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


module core_input
use iso_fortran_env, only: int32, real64


! Maximum number of regions allowed
integer(kind=int32), public, parameter :: maxreg=10

real(kind=real64), public:: gamma

! artificial viscosity
real(kind=real64), public:: cq
real(kind=real64), public:: cl

!max time step
real(kind=real64), public:: maxallstep

! initial time step
real(kind=real64), public:: dtinit

! timestep option
integer, public :: dtoption

! Problem time extent
real(kind=real64), public :: t0  !initial time
real(kind=real64), public :: tf  !final time

! time step growth factor
real(kind=real64), public :: growth

! if 0 cartesian 1 axisymmetric
integer(kind=int32), public :: zaxis

! 0 not vol 1 indiv volume change
integer(kind=int32), public :: zintdivvol

! type of artificial viscosity
integer(kind=int32), public :: avtype

! if 0 no hourglassing if 1 hourglassing
integer(kind=int32), public :: zantihg

! Anti-hourglass filter type
integer(kind=int32), public :: hgregtyp(1:maxreg)

! Kappa for A-H filter
real(kind=real64), public :: kappareg(1:maxreg)

! SILO/TIO output controls
real(kind=real64) :: dtsilo
integer(kind=int32) :: h5type
logical :: tioonefile

! Maximum number of timesteps - used for debugging
integer(kind=int32) :: stepcnt

namelist /tinp/t0,tf,gamma,cq,cl,maxallstep,dtinit,dtoption,growth,zaxis,  &
        zintdivvol,avtype,zantihg,hgregtyp,kappareg,stepcnt,dtsilo,h5type,tioonefile

end module core_input
