!===================================================================
program chocobof
use iso_fortran_env, only: int32

use globalconstants
use geom_data
use sod_init
use lagrangian_hydro
use mesh_data
use core_input
use graphics
use write_silo
use write_tio
use cutoffs, only: dtminhg

implicit none
real(kind=dp) :: ctime0, ctime
real(kind=dp) :: dt05
real(kind=dp) :: totalenergy,totalke,totalie
integer(kind=int32) :: dtcontrol


nadvect=0
stepcnt = 0

! CALL CPU_TIME(ctime0)

! Get mesh size
CALL geominit(nel, nnod, nreg)

! Allocate core mesh arrays
ALLOCATE (xv(1:nnod), yv(1:nnod))
ALLOCATE (znodbound(1:nnod))
ALLOCATE (nodelist(1:4,1:nel))
!*******************************************
!initialise scalars, vectors, matrices     *
!*******************************************
xv=zero
yv=zero
nodelist=0

! set logicals to 0
znodbound=0

CALL geomcalc(nreg, nodelist, znodbound, xv, yv)

!
!============end of geom initialise and connectivity arrays

!==============initialisation for lagstep===============

! time input namelist and eos,artif visc variables

READ(212,NML=tinp)
WRITE(*,NML=tinp)

! initialise physical variables
allocate (pre(1:nel),rho(1:nel),en(1:nel),cc(1:nel),qq(1:nel))
allocate (massel(1:nel),area(1:nel),volel(1:nel),volelold(1:nel))
allocate (uv(1:nnod),vv(1:nnod))
allocate (pre05(1:nel),rho05(1:nel),en05(1:nel),volel05(1:nel))
allocate (divint(1:nel),divvel(1:nel))
allocate (uvold(1:nnod),vvold(1:nnod),uvbar(1:nnod),vvbar(1:nnod))
allocate (xv05(1:nnod),yv05(1:nnod))
allocate (nint(1:4,1:nel),dndx(1:4,1:nel),dndy(1:4,1:nel))
allocate (elwtc(1:4,1:nel))
allocate (pdndx(1:4,1:nel),pdndy(1:4,1:nel))

divint=zero
divvel=zero
en=zero
cc=zero
qq=zero
massel=zero
volel=zero
volelold=zero
area=zero

!=====initialise prdcor var
pre05=zero
rho05=zero
en05=zero
volel05=zero
xv05=zero
yv05=zero
uvold=zero
vvold=zero
uvbar=zero
vvbar=zero

!=initialise fem
elwtc=zero
nint=zero
dndx=zero
dndy=zero

CALL init(nel, nodelist, xv, yv, pre, rho, uv, vv)

! calculate element volume and mass
call calculate_volume(xv, yv, nodelist, nel, volel, area)
call calculate_mass(volel, rho, nel, massel)


Do iel=1,nel
    en(iel)=pre(iel)/((gamma-one)*rho(iel))
END DO

! CALL totalen(en,rho,uv,vv,totalenergy,totalke,totalie)
call calculate_total_energy( &
    en, rho, uv, vv, massel, elwtc, nodelist, nel, zaxis, &
    totalenergy, totalke, totalie &
)

! ====== Time Zero Dump ======
if (h5type == 1) then
    call write_silo_file(1, "gasout", nel, nnod, nodelist, t0, 0, xv, yv, rho, pre, en, uv, vv)
else
    call write_tio_file(1, "gasout", nel, nnod, nodelist, t0, 0, xv, yv, rho, pre, en, uv, vv)
end if

!===================Lagstep predictor corrector drive========
prout=0
time=t0
lastsilo = t0
dt=dtinit
stepno = 1
CALL CPU_TIME(ctime0)
do while (time <= tf)
    !calculate the FEM elements
    call calculate_finite_elements( &
        xv, yv, nodelist, nel, nint, dndx, dndy, pdndx, pdndy, elwtc &
    )

    !calc divergence of v
    call caluclate_div_v(uv, vv, pdndx, pdndy, nodelist, nel, divvel)

    !calculate element sound speed
    call calculate_soundspeed(pre, rho, gamma, nel, cc)

    ! calc element artificial viscosity
    call calculate_q(rho, cc, divvel, area, cq, cl, nel, qq)

    ! calculate stable timestep
    call get_dt( &
        rho, area, cc, qq, time, t0, dtinit,  &
        maxallstep, growth, nel, dt, dtcontrol &
    )
    time = time + dt
    WRITE(*,"(i4, 2x, f10.7, 2x, f15.12, 2x, i4)") stepno,time,dt, dtcontrol

    !calculate 1/2 time step nodal positions
    dt05 = half * dt
    call move_nodes(dt05, xv, yv, uv, vv, nnod, xv05, yv05)

    !calculate the FEM elements
    call calculate_finite_elements( &
        xv05, yv05, nodelist, nel, nint, dndx, dndy, pdndx, pdndy, elwtc &
    )

    ! calculate 1/2 time step volume
    volelold = volel
    call calculate_volume(xv05, yv05, nodelist, nel, volel05, area)

    ! calculate 1/2 time step density
    call calculate_density(massel, volel05, nel, rho05)

    !calc integral of divergence of v
    call calculate_int_divv(zintdivvol, dt05, volel, volelold, uv, vv, dndx, dndy, nodelist, nel, divint)

    !calculate 1/2 time step energy
    call calculate_energy(dt05, pre, qq, massel, en, divint, nel, en05)

    ! calculate 1/2 time step pressure using eos
    call perfect_gas(en05, rho05, gamma, nel, pre05)

    ! momentum equation end timestep
    uvold=uv     !store old velocity values
    vvold=vv

!     CALL momentum(dt,uvold,vvold,rho05,pre05,qq, nint, dndx, dndy,uv,vv)
    call momentum_calculation(dt, dtminhg, zantihg, hgregtyp, kappareg, &
        uvold, vvold, xv05, yv05, rho05, pre05, area, cc, qq,  &
        nint, dndx, dndy, nodelist, znodbound, nel, nnod, uv, vv)

    uvbar=half*(uvold+uv) !calculate an average
    vvbar=half*(vvold+vv)

    !calculate full time step nodal positions
    call move_nodes(dt, xv, yv, uvbar, vvbar, nnod, xv, yv)

    !calculate the FEM elements
    call calculate_finite_elements( &
        xv, yv, nodelist, nel, nint, dndx, dndy, pdndx, pdndy, elwtc &
    )

    ! calculate full time step volume
    call calculate_volume(xv, yv, nodelist, nel, volel, area)

    ! calculate full time step density
    call calculate_density(massel, volel, nel, rho)

    !calculate div with ubar,vbar
    call calculate_int_divv(zintdivvol, dt, volel, volelold, uvbar, vvbar, dndx, dndy, nodelist, nel, divint)

    !calculate full time step energy
    call calculate_energy(dt, pre05, qq, massel, en, divint, nel, en)

    ! calculate full time step pressure using eos
    call perfect_gas(en, rho, gamma, nel, pre)

    if ((time - lastsilo) >= dtsilo) then
        write(*,'("SILO Output file created at time = ",f7.5)') time
        lastsilo = lastsilo + dtsilo
        if (h5type == 1) then
            call write_silo_file(stepno, "gasout", nel, nnod, nodelist, time, stepno, xv, yv, rho, pre, en, uv, vv)
        else
            call write_tio_file(stepno, "gasout", nel, nnod, nodelist, t0, 0, xv, yv, rho, pre, en, uv, vv)
        end if
    end if
    !+++++++++++++++ end Lagstep

    call output(stepcnt, nel, nnod, nodelist, time, stepno, xv, yv, rho, pre, en, uv, vv, volel, prout)
    stepno = stepno + 1
    if (stepno == stepcnt .and. stepcnt > 0) stop
end do

! Time taken
CALL CPU_TIME(ctime)
WRITE (*, *) ' CPU TIME ' , ctime - ctime0, 'seconds'
! READ (*, *)



END PROGRAM ChocoboF
