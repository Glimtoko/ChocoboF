!===================================================================
PROGRAM ChocoboF

!This is the program for the mesh

USE GLOBALCONSTANTS
USE geom_data
use sod_init
USE LAGSTEP
use lagrangian_hydro
use mesh_data
use core_input
USE GRAPHICS
use write_silo
use write_tio

IMPLICIT NONE
REAL     :: ctime0, ctime
! real(kind=dp) :: dtsilo, lastsilo
! INTEGER:: nadvect, stepcnt, h5type
! INTEGER:: nadvect
! NAMELIST /tinp/t0,tf,gamma,cq,cl,maxallstep,dtinit,dtoption,growth,zaxis,  &
!         zintdivvol,avtype,zantihg,hgregtyp,kappareg,stepcnt,dtsilo,h5type,tioonefile

REAL(KIND=DP)     :: dt05
REAL(KIND=DP) :: totalenergy,totalke,totalie
INTEGER :: dtcontrol
nadvect=0
stepcnt = 0

CALL CPU_TIME(ctime0)

! calc xv,yv, nodelist
CALL geominit()
CALL geomcalc()

!
!============end of geom initialise and connectivity arrays

!==============initialisation for lagstep===============

! time input namelist and eos,artif visc variables

READ(212,NML=tinp)
WRITE(*,NML=tinp)

! initialise physical variables
CALL declare()
CALL init()
! calculate element volume and mass
CALL volcalc(xv,yv,volel,area)
CALL masscalc()


Do iel=1,nel
    en(iel)=pre(iel)/((gamma-one)*rho(iel))
END DO

CALL totalen(en,rho,uv,vv,totalenergy,totalke,totalie)

! ====== Time Zero Dump ======
if (h5type == 1) then
    call write_silo_file(1, "gasout")
else
    call write_tio_file(1, "gasout")
end if

!===================Lagstep predictor corrector drive========
prout=0
time=t0
lastsilo = t0
dt=dtinit
stepno = 1
DO WHILE (time.LE.tf)
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
!     CALL energy(dt05,pre,qq,massel,en, divint,en05)
    call calculate_energy(dt05, pre, qq, massel, en, divint, nel, en05)

    ! calculate 1/2 time step pressure using eos
    CALL eos(en05,rho05,pre05)

    ! momentum equation end timestep
    uvold=uv     !store old velocity values
    vvold=vv

    CALL momentum(dt,uvold,vvold,rho05,pre05,qq, nint, dndx, dndy,uv,vv)

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
!     CALL energy(dt,pre05,qq,massel,en,divint,en)
    call calculate_energy(dt, pre05, qq, massel, en, divint, nel, en)

    ! calculate full time step pressure using eos
    CALL eos(en,rho,pre)

    if ((time - lastsilo) >= dtsilo) then
        write(*,'("SILO Output file created at time = ",f7.5)') time
        lastsilo = lastsilo + dtsilo
        if (h5type == 1) then
            call write_silo_file(stepno, "gasout")
        else
            call write_tio_file(stepno, "gasout")
        end if
    end if
    ! error check total energy
!     CALL totalen(en,rho,uv,vv,totalenergy,totalke,totalie)
!     WRITE(27,*) time, prout, totalenergy
    !+++++++++++++++ end Lagstep



    CALL output(stepcnt)
    stepno = stepno + 1
    if (stepno == stepcnt .and. stepcnt > 0) stop
END DO

! Time taken
CALL CPU_TIME(ctime)
WRITE (*, *) ' CPU TIME ' , ctime - ctime0, 'seconds'
! READ (*, *)



END PROGRAM ChocoboF
!
! SUBROUTINE stabletimestep(dtcontrol)
! USE GLOBALCONSTANTS
! USE geom_data
! USE LAGSTEP
! IMPLICIT NONE
!
! integer, intent(out) :: dtcontrol
! REAL (KIND=DP)::dtmin,deltat(1:nel),dtold
!
!
! dtmin=one
! dtold=dt
!
! ! print *, area(1)
! dtcontrol = 0
! if (dtoption == 1) then
!     do ireg=1,nreg
! !         do iel=maxel(ireg-1)+1, maxel(ireg)
!         do iel=1,nel
!             deltat(iel)=area(iel)/MAX(dencut,((cc(iel)**2)+two*(qq(iel)/rho(iel))))
!
!             deltat(iel)=(sqrt(deltat(iel)))/two
!             if (area(iel) < 0.0) print *, "Negative area in cell", iel, deltat(iel)
!             IF (deltat(iel).LT.dtmin) then
!                 dtcontrol = iel
!                 dtmin=deltat(iel)
!             END IF
!         end do
!     end do
! end if
!
! IF (time.EQ.t0) THEN
!     dt=MIN(dtmin, maxallstep, dtinit)
! ELSE
!     dt=MIN(dtmin, maxallstep, growth*dtold)
!     if (dt ==growth*dtold) dtcontrol = -1
! END IF
! IF ((time.GT.0.2500000000000_DP-dt).AND.(time.LT.0.2510050000)) THEN
!     dt= 0.25000000000000_DP-time
! END IF
! WRITE(21,*) time,dt,dtmin
!
! RETURN
! END SUBROUTINE stabletimestep
