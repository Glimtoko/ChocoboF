!===================================================================
program chocobof
use iso_fortran_env, only: int32, real64

use memory_management
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
real(kind=real64) :: ctime0, ctime
real(kind=real64) :: dt05
real(kind=real64) :: totalenergy,totalke,totalie
integer(kind=int32) :: dtcontrol

integer(kind=int32) :: iel

REAL(kind=real64) :: time
INTEGER(kind=int32) :: prout, stepno
REAL(kind=real64) :: dt


nadvect = 0
stepcnt = 0

! Get mesh size
CALL geominit(nel, nnod, nreg)

!*******************************************
!initialise scalars, vectors, matrices     *
!*******************************************
call set_data(xv, nnod)
call set_data(yv, nnod)
call set_data(xv05, nnod)
call set_data(yv05, nnod)

call set_data(nodelist, 4, nel)
call set_data(znodbound, nnod)

call set_data(pre, nel)
call set_data(rho, nel)
call set_data(en, nel)
call set_data(cc, nel)
call set_data(qq, nel)

call set_data(massel, nel)
call set_data(area, nel)
call set_data(volel, nel)
call set_data(volelold, nel)

call set_data(uv, nnod)
call set_data(vv, nnod)
call set_data(uvold, nnod)
call set_data(vvold, nnod)
call set_data(uvbar, nnod)
call set_data(vvbar, nnod)

call set_data(pre05, nel)
call set_data(rho05, nel)
call set_data(en05, nel)
call set_data(volel05, nel)

call set_data(divint, nel)
call set_data(divvel, nel)

call set_data(nint, 4, nel)
call set_data(dndx, 4, nel)
call set_data(dndy, 4, nel)
call set_data(elwtc, 4, nel)
call set_data(pdndx, 4, nel)
call set_data(pdndy, 4, nel)

CALL geomcalc(nreg, nodelist, znodbound, xv, yv)

!
!============end of geom initialise and connectivity arrays

!==============initialisation for lagstep===============

! time input namelist and eos,artif visc variables

READ(212,NML=tinp)
WRITE(*,NML=tinp)

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
