!===================================================================
program chocobof
use iso_fortran_env, only: int32, real64

use fileunits
use memory_management
use mesh_mod

use mesh_reader
use sod_init
use lagrangian_hydro
use core_input
use text_output
use write_silo
use write_tio

implicit none
real(kind=real64) :: ctime0, ctime
real(kind=real64) :: totalenergy,totalke,totalie
integer(kind=int32) :: dtcontrol

integer(kind=int32) :: iel, mat

REAL(kind=real64) :: time
INTEGER(kind=int32) :: prout, stepno
REAL(kind=real64) :: dt, dt05

integer(kind=int32) :: lastsilo

character(len=15) :: problemname
logical :: use_spherical_sod
integer(kind=int32) :: status

type(MeshT) :: mesh

! Get the problem name from the command line
use_spherical_sod = .false.
call get_command_argument(number=1, value=problemname, status=status)

if (status > 0) then
    use_spherical_sod = .true.
elseif (status < 0) then
    print *, "Error: Invalid problem name (probably > 15 characters)"
    stop 99
end if

! Open input file
open(unit=control, file='param.dat')

! Get mesh size
if (use_spherical_sod) then
    call sod_get_mesh_size(mesh%nel, mesh%nnod, mesh%nreg)
else
    call get_mesh_size(problemname, mesh%nel, mesh%nnod, mesh%nreg)
end if

! For now, assume # materials = # regions
mesh%nmat = mesh%nreg


! Initialise scalars, vectors, matrices
call set_data(mesh%xv, mesh%nnod)
call set_data(mesh%yv, mesh%nnod)
call set_data(mesh%xv05, mesh%nnod)
call set_data(mesh%yv05, mesh%nnod)

call set_data(mesh%nodelist, 4, mesh%nel)
call set_data(mesh%region, mesh%nel)
call set_data(mesh%material, mesh%nel)
call set_data(mesh%znodbound, mesh%nnod)
call set_data(mesh%regiontocell, 2, mesh%nreg)

call set_data(mesh%gamma, mesh%nmat)

call set_data(mesh%pre, mesh%nel)
call set_data(mesh%rho, mesh%nel)
call set_data(mesh%en, mesh%nel)
call set_data(mesh%cc, mesh%nel)
call set_data(mesh%qq, mesh%nel)

call set_data(mesh%massel, mesh%nel)
call set_data(mesh%area, mesh%nel)
call set_data(mesh%volel, mesh%nel)
call set_data(mesh%volelold, mesh%nel)

call set_data(mesh%uv, mesh%nnod)
call set_data(mesh%vv, mesh%nnod)
call set_data(mesh%uvold, mesh%nnod)
call set_data(mesh%vvold, mesh%nnod)
call set_data(mesh%uvbar, mesh%nnod)
call set_data(mesh%vvbar, mesh%nnod)

call set_data(mesh%pre05, mesh%nel)
call set_data(mesh%rho05, mesh%nel)
call set_data(mesh%en05, mesh%nel)
call set_data(mesh%volel05, mesh%nel)

call set_data(mesh%divint, mesh%nel)
call set_data(mesh%divvel, mesh%nel)

call set_data(mesh%nint, 4, mesh%nel)
call set_data(mesh%dndx, 4, mesh%nel)
call set_data(mesh%dndy, 4, mesh%nel)
call set_data(mesh%elwtc, 4, mesh%nel)
call set_data(mesh%pdndx, 4, mesh%nel)
call set_data(mesh%pdndy, 4, mesh%nel)

! Generate mesh geometry
if (use_spherical_sod) then
    call sod_generate_mesh(mesh%nreg, mesh%nodelist, mesh%znodbound, mesh%xv, mesh%yv)
    mesh%region = 1
    mesh%material = 1
else
    call generate_mesh(problemname, mesh%nodelist, mesh%znodbound, mesh%region, mesh%material, mesh%regiontocell,mesh%xv, mesh%yv)
end if


! Read user input
read(control, nml=tinp)
write(*, nml=tinp)

! Initialise spherical sod problem
if (use_spherical_sod) then
    call sod_initial_conditions(mesh%nel, mesh%nodelist, mesh%xv, mesh%yv, mesh%pre, mesh%rho, mesh%uv, mesh%vv)
    mesh%gamma(1) = gamma
else
    call initial_conditions(problemname, mesh%nel, mesh%nmat, mesh%material, mesh%pre, mesh%rho, mesh%gamma)
end if

! Calculate element volume and mass
call calculate_volume(mesh%xv, mesh%yv, mesh%nodelist, mesh%nel, mesh%volel, mesh%area)
call calculate_mass(mesh%volel, mesh%rho, mesh%nel, mesh%massel)


do iel = 1, mesh%nel
    mat = mesh%material(iel)
    mesh%en(iel) = mesh%pre(iel)/((mesh%gamma(mat) - 1.0_real64)*mesh%rho(iel))
end do

call calculate_total_energy( &
    mesh%en, mesh%rho, mesh%uv, mesh%vv, mesh%massel, mesh%elwtc, mesh%nodelist, mesh%nel, zaxis, &
    totalenergy, totalke, totalie &
)

! Time zero visualisation
if (h5type == 1) then
    call write_silo_file(1, "gasout", mesh%nel, mesh%nnod, mesh%nodelist, t0, &
        0, mesh%xv, mesh%yv, mesh%rho, mesh%pre, mesh%en, mesh%uv, mesh%vv)
else
    call write_tio_file(1, "gasout", mesh%nel, mesh%nnod, mesh%nodelist, t0, &
        0, mesh%xv, mesh%yv, mesh%rho, mesh%pre, mesh%en, mesh%uv, mesh%vv)
end if


prout = 0
time = t0
lastsilo = t0
dt = dtinit
stepcnt = 0
stepno = 1
call CPU_TIME(ctime0)
do while (time <= tf)
    !calculate the FEM elements
    call calculate_finite_elements( &
        mesh%xv, mesh%yv, mesh%nodelist, mesh%nel, mesh%nint, mesh%dndx, mesh%dndy, mesh%pdndx, mesh%pdndy, mesh%elwtc &
    )

    !calc divergence of v
    call caluclate_div_v(mesh%uv, mesh%vv, mesh%pdndx, mesh%pdndy, mesh%nodelist, mesh%nel, mesh%divvel)

    !calculate element sound speed
    call calculate_soundspeed(mesh%pre, mesh%rho, mesh%material, mesh%gamma, mesh%nel, mesh%cc)

    ! calc element artificial viscosity
    call calculate_q(mesh%rho, mesh%cc, mesh%divvel, mesh%area, cq, cl, mesh%nel, mesh%qq)

    ! calculate stable timestep
    call get_dt( &
        mesh%rho, mesh%area, mesh%cc, mesh%qq, time, t0, dtinit,  &
        maxallstep, growth, mesh%nel, dt, dtcontrol &
    )
    time = time + dt
    WRITE(*,"(i4, 2x, f10.7, 2x, f15.12, 2x, i4)") stepno,time,dt, dtcontrol

    !calculate 1/2 time step nodal positions
    dt05 = 0.5_real64 * dt
    call move_nodes(dt05, mesh%xv, mesh%yv, mesh%uv, mesh%vv, mesh%nnod, mesh%xv05, mesh%yv05)

    !calculate the FEM elements
    call calculate_finite_elements( &
        mesh%xv05, mesh%yv05, mesh%nodelist, mesh%nel, mesh%nint, mesh%dndx, mesh%dndy, mesh%pdndx, mesh%pdndy, mesh%elwtc &
    )

    ! calculate 1/2 time step volume
    mesh%volelold = mesh%volel
    call calculate_volume(mesh%xv05, mesh%yv05, mesh%nodelist, mesh%nel, mesh%volel05, mesh%area)

    ! calculate 1/2 time step density
    call calculate_density(mesh%massel, mesh%volel05, mesh%nel, mesh%rho05)

    !calc integral of divergence of v
    call calculate_int_divv(zintdivvol, dt05, mesh%volel, mesh%volelold, &
        mesh%uv, mesh%vv, mesh%dndx, mesh%dndy, mesh%nodelist, mesh%nel, mesh%divint)

    !calculate 1/2 time step energy
    call calculate_energy(dt05, mesh%pre, mesh%qq, mesh%massel, mesh%en, mesh%divint, mesh%nel, mesh%en05)

    ! calculate 1/2 time step mesh%pressure using eos
    call perfect_gas(mesh%en05, mesh%rho05, mesh%material, mesh%gamma, mesh%nel, mesh%pre05)

    ! momentum equation end timestep
    mesh%uvold=mesh%uv     !store old velocity values
    mesh%vvold=mesh%vv

    call momentum_calculation(dt, zantihg, hgregtyp, kappareg, &
        mesh%uvold, mesh%vvold, mesh%xv05, mesh%yv05, mesh%rho05, mesh%pre05, mesh%area, mesh%cc, mesh%qq,  &
        mesh%nint, mesh%dndx, mesh%dndy, mesh%nodelist, mesh%znodbound, mesh%nel, mesh%nnod, mesh%uv, mesh%vv)

    mesh%uvbar=0.5_real64*(mesh%uvold+mesh%uv) !calculate an average
    mesh%vvbar=0.5_real64*(mesh%vvold+mesh%vv)

    !calculate full time step nodal positions
    call move_nodes(dt, mesh%xv, mesh%yv, mesh%uvbar, mesh%vvbar, mesh%nnod, mesh%xv, mesh%yv)

    !calculate the FEM elements
    call calculate_finite_elements( &
        mesh%xv, mesh%yv, mesh%nodelist, mesh%nel, mesh%nint, mesh%dndx, mesh%dndy, mesh%pdndx, mesh%pdndy, mesh%elwtc &
    )

    ! calculate full time step volume
    call calculate_volume(mesh%xv, mesh%yv, mesh%nodelist, mesh%nel, mesh%volel, mesh%area)

    ! calculate full time step density
    call calculate_density(mesh%massel, mesh%volel, mesh%nel, mesh%rho)

    !calculate div with ubar,vbar
    call calculate_int_divv(zintdivvol, dt, mesh%volel, mesh%volelold, mesh%uvbar, &
        mesh%vvbar, mesh%dndx, mesh%dndy, mesh%nodelist, mesh%nel, mesh%divint)

    !calculate full time step energy
    call calculate_energy(dt, mesh%pre05, mesh%qq, mesh%massel, mesh%en, mesh%divint, mesh%nel, mesh%en)

    ! calculate full time step mesh%pressure using eos
    call perfect_gas(mesh%en, mesh%rho, mesh%material, mesh%gamma, mesh%nel, mesh%pre)

    if ((time - lastsilo) >= dtsilo) then
        write(*,'("SILO Output file created at time = ",f7.5)') time
        lastsilo = lastsilo + dtsilo
        if (h5type == 1) then
            call write_silo_file(stepno, "gasout", mesh%nel, mesh%nnod, &
                mesh%nodelist, time, stepno, mesh%xv, mesh%yv, mesh%rho, mesh%pre, mesh%en, mesh%uv, mesh%vv)
        else
            call write_tio_file(stepno, "gasout", mesh%nel, mesh%nnod, &
                mesh%nodelist, t0, 0, mesh%xv, mesh%yv, mesh%rho, mesh%pre, mesh%en, mesh%uv, mesh%vv)
        end if
    end if

    if (time >= 0.20) then
        call output(problemname, stepcnt, mesh%nel, mesh%nnod, mesh%nodelist, time, stepno, &
            mesh%xv, mesh%yv, mesh%rho, mesh%pre, mesh%en, mesh%uv, mesh%vv, mesh%volel, prout)
    endif
    stepno = stepno + 1
    if (stepno == stepcnt .and. stepcnt > 0) stop
end do

! Time taken
call cpu_time(ctime)
write(*,*)
write(*,*) 'Time taken = ' , ctime - ctime0, 'seconds'

end program ChocoboF
