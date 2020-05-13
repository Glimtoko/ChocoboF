module lagrangian_hydro
! This is a copy of the lagstep module, but all routines in it have fully
! declared APIs and do not rely on global data. Also, I've given them more
! expressive names
use globalconstants
use iso_fortran_env, only: int32
contains

subroutine calculate_finite_elements(x, y, nodelist, nel, ni, dndx, dndy, pdndx, pdndy, elwtc)
    implicit none
    real(kind=dp), dimension(:), intent(in) ::  x, y
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel

    real(kind=dp), dimension(:,:), intent(out) :: ni, dndx, dndy
    real(kind=dp), dimension(:,:), intent(out) :: pdndx, pdndy
    real(kind=dp), dimension(:,:), intent(out) :: elwtc


    real(kind=dp) ::jacob

    integer :: iel, jjj

    real(kind=dp) :: a1,a2,a3,b1,b2,b3


    do iel = 1,nel
        a1 = quarter*(-x(nodelist(1,iel))+x(nodelist(2,iel))   &
                    +x(nodelist(3,iel))-x(nodelist(4,iel)))
        a2 = quarter*(x(nodelist(1,iel))-x(nodelist(2,iel))   &
                    +x(nodelist(3,iel))-x(nodelist(4,iel)))
        a3 = quarter*(-x(nodelist(1,iel))-x(nodelist(2,iel))   &
                    +x(nodelist(3,iel))+x(nodelist(4,iel)))
        b1 = quarter*(-y(nodelist(1,iel))+y(nodelist(2,iel))   &
                    +y(nodelist(3,iel))-y(nodelist(4,iel)))
        b2 = quarter*(y(nodelist(1,iel))-y(nodelist(2,iel))   &
                    +y(nodelist(3,iel))-y(nodelist(4,iel)))
        b3 = quarter*(-y(nodelist(1,iel))-y(nodelist(2,iel))   &
                    +y(nodelist(3,iel))+y(nodelist(4,iel)))
        ! for momentum
        ni(1,iel) = ((three*b3-b2)*(three*a1-a2)-(three*a3-a2)*(three*b1-b2))/nine
        ni(2,iel) = ((three*b3+b2)*(three*a1-a2)-(three*a3+a2)*(three*b1-b2))/nine
        ni(3,iel) = ((three*b3+b2)*(three*a1+a2)-(three*a3+a2)*(three*b1+b2))/nine
        ni(4,iel) = ((three*b3-b2)*(three*a1+a2)-(three*a3-a2)*(three*b1+b2))/nine

        ! integrated ddndx's for intdiv and energy
        dndx(1,iel) = -b3+b1
        dndx(2,iel) =  b3+b1
        dndx(3,iel) =  b3-b1
        dndx(4,iel) = -b3-b1

        dndy(1,iel) =  a3-a1
        dndy(2,iel) = -a3-a1
        dndy(3,iel) = -a3+a1
        dndy(4,iel) =  a3+a1

        ! jacobian and ddndx's for div and viscosity
        jacob = a1*b3-a3*b1
        do jjj = 1,4
            pdndx(jjj,iel) = quarter*dndx(jjj,iel)/jacob
            pdndy(jjj,iel) = quarter*dndy(jjj,iel)/jacob
        end do

        ! calculate elwtc for sale integral x over volume
        elwtc(1,iel) = ni(1,iel)
        elwtc(2,iel) = ni(2,iel)
        elwtc(3,iel) = ni(3,iel)
        elwtc(4,iel) = ni(4,iel)
    end do
end subroutine calculate_finite_elements


subroutine caluclate_div_v(u, v, pdndx, pdndy, nodelist, nel, divvel)
    ! Not really divergence just a 2d measure of gradient accross cell
    ! Calculates divergence for viscosity
    implicit none
    real (kind=dp), dimension(:), intent(in) ::  u, v
    real (kind=dp), dimension(:,:), intent(in) ::  pdndx, pdndy
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel

    real (kind=dp), dimension(:), intent(out):: divvel

    real (kind=dp) ::dterm
    integer(kind=int32) :: iel, j

    divvel = zero
    do iel = 1,nel
        do j = 1,4
            dterm = u(nodelist(j,iel))*pdndx(j,iel) + v(nodelist(j,iel))*pdndy(j,iel)
            divvel(iel) = divvel(iel) + dterm
        end do
    end do

    return
end subroutine caluclate_div_v


subroutine calculate_soundspeed(pressure, rho, gamma, nel, cc)
    !calculate element sound speed
    ! ideal gas equation gamma is defined as parameter at init
    implicit none
    real(kind=dp), dimension(:), intent(in) :: pressure, rho
    real(kind=dp), intent(in) :: gamma
    integer(kind=int32) :: nel

    real(kind=dp), dimension(:), intent(out) :: cc

    integer(kind=int32) :: iel

    do iel=1,nel
        cc(iel)=sqrt(gamma*pressure(iel)/rho(iel))
    end do

end subroutine calculate_soundspeed


subroutine calculate_q(rho, cc, divvel, area, cq, cl, nel, q)
    ! use core_input, only: cq, cl
    implicit none
    real(kind=dp), dimension(:), intent(in) :: divvel, rho, cc, area
    real(kind=dp), intent(in) :: cq, cl
    integer(kind=int32), intent(in) :: nel

    real (kind=dp), dimension(:), intent(out) :: q

    real(kind=dp) :: dudx
    integer(kind=int32) :: iel

    ! bulk q note our cl is andy's/2
    do iel = 1,nel
        if(divvel(iel) < zero) then
            dudx = sqrt(area(iel))*divvel(iel)
            q(iel) = (cq*rho(iel)*(dudx)**2) +    &
                (cl*rho(iel)*cc(iel)*abs(dudx))
        else
            q(iel) = zero
        end if
    end do
end subroutine calculate_q


subroutine get_dt(rho, area, cc, q, time, t0, dtinit, dtmax, growth, nel, dt, dtcontrol)
    use cutoffs, only: dencut
    implicit none
    real(kind=dp), dimension(:), intent(in) :: rho, area, cc, q
    real(kind=dp), intent(in) :: time, t0, dtinit, dtmax, growth
    integer(kind=int32), intent(in) :: nel

    real(kind=dp), intent(inout) :: dt
    integer(kind=int32), intent(out) :: dtcontrol

    real(kind=dp) :: dtmin, dtold
    real(kind=dp), dimension(nel) :: deltat
    integer(kind=int32) :: iel


    dtmin = one
    dtold = dt

    dtcontrol = 0
    do iel=1,nel
        deltat(iel) = area(iel)/max(dencut,((cc(iel)**2)+two*(q(iel)/rho(iel))))

        deltat(iel) = (sqrt(deltat(iel)))/two
        if (area(iel) < 0.0_dp) print *, "negative area in cell", iel, deltat(iel)
        if (deltat(iel) < dtmin) then
            dtcontrol = iel
            dtmin = deltat(iel)
        end if
    end do

    if (time == t0) then
        dt = min(dtmin, dtmax, dtinit)
    else
        dt = min(dtmin, dtmax, growth*dtold)
        if (dt == growth*dtold) dtcontrol = -1
    end if
    if ((time > 0.25_dp - dt) .and. (time < 0.251005_dp)) then
        dt = 0.25_dp - time
    end if
    write(21,*) time, dt, dtmin

    return
end subroutine get_dt


subroutine move_nodes(dt, x, y, u, v, nnod, xout, yout)
    implicit none
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:), intent(in) :: x, y, u, v
    integer(kind=int32), intent(in) :: nnod

    real(kind=dp), dimension(:), intent(out) :: xout, yout

    integer(kind=int32) :: inod

    do inod = 1,nnod
        xout(inod) = x(inod) + dt*u(inod)
        yout(inod) = y(inod) + dt*v(inod)
    end do
end subroutine move_nodes


subroutine calculate_volume(x, y, nodelist, nel, volume, area)
    !calculate element area of quadralateral using jacobian
    implicit none
    real(kind=dp), dimension(:), intent(in)::  x, y
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel

    real(kind=dp), dimension(:), intent(out) :: volume, area

    real(kind=dp) :: a1,a2,a3,b1,b2,b3

    integer(kind=int32) :: iel

    do iel = 1,nel
        a1 = quarter*(-x(nodelist(1,iel))+x(nodelist(2,iel))   &
                      +x(nodelist(3,iel))-x(nodelist(4,iel)))
        a2 = quarter*(x(nodelist(1,iel))-x(nodelist(2,iel))   &
                      +x(nodelist(3,iel))-x(nodelist(4,iel)))
        a3 = quarter*(-x(nodelist(1,iel))-x(nodelist(2,iel))   &
                      +x(nodelist(3,iel))+x(nodelist(4,iel)))
        b1 = quarter*(-y(nodelist(1,iel))+y(nodelist(2,iel))   &
                      +y(nodelist(3,iel))-y(nodelist(4,iel)))
        b2 = quarter*(y(nodelist(1,iel))-y(nodelist(2,iel))   &
                      +y(nodelist(3,iel))-y(nodelist(4,iel)))
        b3 = quarter*(-y(nodelist(1,iel))-y(nodelist(2,iel))   &
                      +y(nodelist(3,iel))+y(nodelist(4,iel)))

        volume(iel) = four*(a1*b3 - a3*b1)
        area(iel) = volume(iel)
    end do

    return
end subroutine calculate_volume


subroutine calculate_density(mass, volume, nel, rho)
    implicit none
    real(kind=dp), dimension(:), intent(in)::  mass, volume
    integer(kind=int32) :: nel

    real(kind=dp), dimension(:), intent(out):: rho

    integer(kind=int32) :: iel

    do iel=1,nel
        rho(iel)=mass(iel)/volume(iel)
    end do
end subroutine calculate_density


subroutine calculate_int_divv(zintdivvol, dt, vol, volold, u, v, dndx, dndy, nodelist, nel, intdiv)
    implicit none
    integer(kind=int32), intent(in) :: zintdivvol
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:), intent(in) :: u, v, vol, volold
    real(kind=dp), dimension(:,:), intent(in) :: dndx, dndy
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel

    real(kind=dp), dimension(:), intent(out) :: intdiv

    integer(kind=int32) :: iel, j
    real(kind=dp) :: eterm

    if (zintdivvol == 0) then
        intdiv=zero
        do iel=1, nel
            do j=1, 4
                eterm = u(nodelist(j,iel))*dndx(j,iel)         &
                       +v(nodelist(j,iel))*dndy(j,iel)
                intdiv(iel) = intdiv(iel)+eterm
            end do
        end do
    else
        do iel=1,nel
            intdiv(iel)=(vol(iel)-volold(iel))/(dt)
        end do
    end if
end subroutine calculate_int_divv


subroutine calculate_energy(dt, press, visc, mass, ener, intdiv, nel, enout)
    implicit none
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:), intent(in) ::  press, visc
    real(kind=dp), dimension(:), intent(in) ::  mass
    real(kind=dp), dimension(:), pointer :: ener
    real(kind=dp), dimension(:), intent(in) ::  intdiv
    integer(kind=int32), intent(in) :: nel

    real (kind=dp), dimension(:), pointer :: enout

    integer(kind=int32) :: iel

    do iel=1,nel
        enout(iel) = ener(iel) - dt*(press(iel)+visc(iel))*intdiv(iel)/mass(iel)
    end do
end subroutine calculate_energy

end module lagrangian_hydro
