module lagrangian_hydro
! This is a copy of the lagstep module, but all routines in it have fully
! declared APIs and do not rely on global data. Also, I've given them more
! expressive names

use iso_fortran_env, only: int32, real64

! Cut-off values
real(kind=real64), public, parameter :: dtminhg = 0.000000008
real(kind=real64), public, parameter :: dencut = 1.000000000e-6

contains


subroutine calculate_mass(volel, rho, nel, massel)
    !calculate element mass=element density* element volume
    implicit none
    real(kind=real64), dimension(:), intent(in) :: volel, rho
    integer(kind=int32), intent(in) :: nel

    real(kind=real64), dimension(:), intent(out) :: massel

    integer(kind=int32) :: iel

    do iel=1,nel
        massel(iel)=volel(iel)*rho(iel)
    end do
    return
end subroutine calculate_mass


subroutine calculate_total_energy( &
    energy, rho, u, v, mass, elwtc, nodelist, nel, zaxis, &
    total_energy, total_ke, total_ie &
)
    implicit none
    real (kind=real64), dimension(:), intent(in) ::  energy, rho
    real (kind=real64), dimension(:), intent(in) ::  u, v, mass
    real(kind=real64), dimension(:,:), intent(in) :: elwtc
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel, zaxis

    real (kind=real64), intent(out):: total_energy, total_ie, total_ke

    real (kind=real64) :: elenergy
    real (kind=real64) :: two_pi, tek
    integer(kind=int32) :: iel, j

    real(kind=real64), parameter ::pi = 3.14159265358979324_real64

    if (zaxis == 0) then
        two_pi = 1.0_real64
    else if (zaxis == 1) then
        two_pi = 2.0_real64 * pi
    end if

    total_energy = 0.0_real64
    total_ke = 0.0_real64
    total_ie = 0.0_real64
    do iel = 1,nel
        tek = 0.0_real64

        do j = 1,4
            tek =  tek + 0.5_real64*rho(iel)*elwtc(j,iel)*   &
                   (u(nodelist(j,iel))**2 + v(nodelist(j,iel))**2)*two_pi
        end do
        elenergy = mass(iel)*energy(iel)*two_pi + tek
        total_ke = total_ke + tek
        total_ie = total_ie + mass(iel)*energy(iel)*two_pi
        total_energy = total_energy+elenergy
    end do
end subroutine calculate_total_energy


subroutine calculate_finite_elements(x, y, nodelist, nel, ni, dndx, dndy, pdndx, pdndy, elwtc)
    implicit none
    real(kind=real64), dimension(:), intent(in) ::  x, y
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel

    real(kind=real64), dimension(:,:), intent(out) :: ni, dndx, dndy
    real(kind=real64), dimension(:,:), intent(out) :: pdndx, pdndy
    real(kind=real64), dimension(:,:), intent(out) :: elwtc


    real(kind=real64) ::jacob

    integer :: iel, jjj

    real(kind=real64) :: a1, a2, a3, b1, b2, b3


    do iel = 1,nel
        a1 = 0.25_real64*(-x(nodelist(1,iel))+x(nodelist(2,iel))   &
                    +x(nodelist(3,iel))-x(nodelist(4,iel)))
        a2 = 0.25_real64*(x(nodelist(1,iel))-x(nodelist(2,iel))   &
                    +x(nodelist(3,iel))-x(nodelist(4,iel)))
        a3 = 0.25_real64*(-x(nodelist(1,iel))-x(nodelist(2,iel))   &
                    +x(nodelist(3,iel))+x(nodelist(4,iel)))
        b1 = 0.25_real64*(-y(nodelist(1,iel))+y(nodelist(2,iel))   &
                    +y(nodelist(3,iel))-y(nodelist(4,iel)))
        b2 = 0.25_real64*(y(nodelist(1,iel))-y(nodelist(2,iel))   &
                    +y(nodelist(3,iel))-y(nodelist(4,iel)))
        b3 = 0.25_real64*(-y(nodelist(1,iel))-y(nodelist(2,iel))   &
                    +y(nodelist(3,iel))+y(nodelist(4,iel)))
        ! for momentum
        ni(1,iel) = ((3.0_real64*b3-b2)*(3.0_real64*a1-a2)-(3.0_real64*a3-a2)*(3.0_real64*b1-b2))/9.0_real64
        ni(2,iel) = ((3.0_real64*b3+b2)*(3.0_real64*a1-a2)-(3.0_real64*a3+a2)*(3.0_real64*b1-b2))/9.0_real64
        ni(3,iel) = ((3.0_real64*b3+b2)*(3.0_real64*a1+a2)-(3.0_real64*a3+a2)*(3.0_real64*b1+b2))/9.0_real64
        ni(4,iel) = ((3.0_real64*b3-b2)*(3.0_real64*a1+a2)-(3.0_real64*a3-a2)*(3.0_real64*b1+b2))/9.0_real64

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
            pdndx(jjj,iel) = 0.25_real64*dndx(jjj,iel)/jacob
            pdndy(jjj,iel) = 0.25_real64*dndy(jjj,iel)/jacob
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
    real (kind=real64), dimension(:), intent(in) ::  u, v
    real (kind=real64), dimension(:,:), intent(in) ::  pdndx, pdndy
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel

    real (kind=real64), dimension(:), intent(out):: divvel

    real (kind=real64) ::dterm
    integer(kind=int32) :: iel, j

    divvel = 0.0_real64
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
    real(kind=real64), dimension(:), intent(in) :: pressure, rho
    real(kind=real64), intent(in) :: gamma
    integer(kind=int32) :: nel

    real(kind=real64), dimension(:), intent(out) :: cc

    integer(kind=int32) :: iel

    do iel=1,nel
        cc(iel)=sqrt(gamma*pressure(iel)/rho(iel))
    end do

end subroutine calculate_soundspeed


subroutine calculate_q(rho, cc, divvel, area, cq, cl, nel, q)
    ! use core_input, only: cq, cl
    implicit none
    real(kind=real64), dimension(:), intent(in) :: divvel, rho, cc, area
    real(kind=real64), intent(in) :: cq, cl
    integer(kind=int32), intent(in) :: nel

    real (kind=real64), dimension(:), intent(out) :: q

    real(kind=real64) :: dudx
    integer(kind=int32) :: iel

    ! bulk q note our cl is andy's/2
    do iel = 1,nel
        if(divvel(iel) < 0.0_real64) then
            dudx = sqrt(area(iel))*divvel(iel)
            q(iel) = (cq*rho(iel)*(dudx)**2) +    &
                (cl*rho(iel)*cc(iel)*abs(dudx))
        else
            q(iel) = 0.0_real64
        end if
    end do
end subroutine calculate_q


subroutine get_dt(rho, area, cc, q, time, t0, dtinit, dtmax, growth, nel, dt, dtcontrol)
    implicit none
    real(kind=real64), dimension(:), intent(in) :: rho, area, cc, q
    real(kind=real64), intent(in) :: time, t0, dtinit, dtmax, growth
    integer(kind=int32), intent(in) :: nel

    real(kind=real64), intent(inout) :: dt
    integer(kind=int32), intent(out) :: dtcontrol

    real(kind=real64) :: dtmin, dtold
    real(kind=real64), dimension(nel) :: deltat
    integer(kind=int32) :: iel


    dtmin = 1.0_real64
    dtold = dt

    dtcontrol = 0
    do iel=1,nel
        deltat(iel) = area(iel)/max(dencut,((cc(iel)**2)+2.0_real64*(q(iel)/rho(iel))))

        deltat(iel) = (sqrt(deltat(iel)))/2.0_real64
        if (area(iel) < 0.0_real64) print *, "negative area in cell", iel, deltat(iel)
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
    if ((time > 0.25_real64 - dt) .and. (time < 0.251005_real64)) then
        dt = 0.25_real64 - time
    end if
    write(21,*) time, dt, dtmin

    return
end subroutine get_dt


subroutine move_nodes(dt, x, y, u, v, nnod, xout, yout)
    implicit none
    real(kind=real64), intent(in) :: dt
    real(kind=real64), dimension(:), intent(in) :: x, y, u, v
    integer(kind=int32), intent(in) :: nnod

    real(kind=real64), dimension(:), intent(out) :: xout, yout

    integer(kind=int32) :: inod

    do inod = 1,nnod
        xout(inod) = x(inod) + dt*u(inod)
        yout(inod) = y(inod) + dt*v(inod)
    end do
end subroutine move_nodes


subroutine calculate_volume(x, y, nodelist, nel, volume, area)
    !calculate element area of quadralateral using jacobian
    implicit none
    real(kind=real64), dimension(:), intent(in)::  x, y
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel

    real(kind=real64), dimension(:), intent(out) :: volume, area

    real(kind=real64) :: a1,a2,a3,b1,b2,b3

    integer(kind=int32) :: iel

    do iel = 1,nel
        a1 = 0.25_real64*(-x(nodelist(1,iel))+x(nodelist(2,iel))   &
                      +x(nodelist(3,iel))-x(nodelist(4,iel)))
        a2 = 0.25_real64*(x(nodelist(1,iel))-x(nodelist(2,iel))   &
                      +x(nodelist(3,iel))-x(nodelist(4,iel)))
        a3 = 0.25_real64*(-x(nodelist(1,iel))-x(nodelist(2,iel))   &
                      +x(nodelist(3,iel))+x(nodelist(4,iel)))
        b1 = 0.25_real64*(-y(nodelist(1,iel))+y(nodelist(2,iel))   &
                      +y(nodelist(3,iel))-y(nodelist(4,iel)))
        b2 = 0.25_real64*(y(nodelist(1,iel))-y(nodelist(2,iel))   &
                      +y(nodelist(3,iel))-y(nodelist(4,iel)))
        b3 = 0.25_real64*(-y(nodelist(1,iel))-y(nodelist(2,iel))   &
                      +y(nodelist(3,iel))+y(nodelist(4,iel)))

        volume(iel) = 4.0_real64*(a1*b3 - a3*b1)
        area(iel) = volume(iel)
    end do

    return
end subroutine calculate_volume


subroutine calculate_density(mass, volume, nel, rho)
    implicit none
    real(kind=real64), dimension(:), intent(in)::  mass, volume
    integer(kind=int32) :: nel

    real(kind=real64), dimension(:), intent(out):: rho

    integer(kind=int32) :: iel

    do iel=1,nel
        rho(iel)=mass(iel)/volume(iel)
    end do
end subroutine calculate_density


subroutine calculate_int_divv(zintdivvol, dt, vol, volold, u, v, dndx, dndy, nodelist, nel, intdiv)
    implicit none
    integer(kind=int32), intent(in) :: zintdivvol
    real(kind=real64), intent(in) :: dt
    real(kind=real64), dimension(:), intent(in) :: u, v, vol, volold
    real(kind=real64), dimension(:,:), intent(in) :: dndx, dndy
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel

    real(kind=real64), dimension(:), intent(out) :: intdiv

    integer(kind=int32) :: iel, j
    real(kind=real64) :: eterm

    if (zintdivvol == 0) then
        intdiv=0.0_real64
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
    real(kind=real64), intent(in) :: dt
    real(kind=real64), dimension(:), intent(in) ::  press, visc
    real(kind=real64), dimension(:), intent(in) ::  mass
    real(kind=real64), dimension(:), intent(in) :: ener
    real(kind=real64), dimension(:), intent(in) ::  intdiv
    integer(kind=int32), intent(in) :: nel

    real (kind=real64), dimension(:), intent(out) :: enout

    integer(kind=int32) :: iel

    do iel=1,nel
        enout(iel) = ener(iel) - dt*(press(iel)+visc(iel))*intdiv(iel)/mass(iel)
    end do
end subroutine calculate_energy


subroutine perfect_gas(energy, rho, gamma, nel, pressure)
    implicit none
    real(kind=real64), dimension(:), intent(in) ::  energy, rho
    real(kind=real64), intent(in) :: gamma
    integer(kind=int32), intent(in) :: nel

    real(kind=real64), dimension(:), intent(out) :: pressure

    integer(kind=int32) :: iel

    do iel=1,nel
        pressure(iel) = (gamma - 1.0_real64)*rho(iel)*energy(iel)
    end do
end subroutine perfect_gas


subroutine hourglass_filter( &
    u, v, x05, y05, rho, area, cc, &
    dt, dndx, dndy, hgregtyp, kappareg, nodelist, nel, &
    fx, fy &
)

    ! 26th june
    !calculate anti-hourglass filters
    !called in momentum subroutine

    implicit none
    real (kind=real64), intent(in), dimension(:) :: u, v, x05, y05, rho, area, cc
    real (kind=real64), intent(in) :: dt
    real (kind=real64), intent(in), dimension(:,:) :: dndx, dndy
    integer(kind=int32), intent(in), dimension(:) :: hgregtyp
    real (kind=real64), intent(in), dimension(:) :: kappareg
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), intent(in) :: nel

    real (kind=real64), intent(inout)::fx(:),fy(:)


    real (kind=real64)  :: ugam, vgam, temp, biibii, xdiff, ydiff
    real (kind=real64)  :: gam1, gam2, gam3, gam4, qx, qy, kap
    integer(kind=int32) :: hgtyp, iel

    real(kind=real64) :: a1,a2,a3,b1,b2,b3

    ! do ireg=1,nreg   ! loop regions

        hgtyp = hgregtyp(1)

        if (hgtyp == 3) then
        ! dyna hourglass type 2 with different l and no speed
            do iel=1, nel
            ! calculate u  gamma
            !            ik     k
            ugam=u(nodelist(1,iel))-u(nodelist(2,iel))    &
                +u(nodelist(3,iel))-u(nodelist(4,iel))
            vgam=v(nodelist(1,iel))-v(nodelist(2,iel))    &
                +v(nodelist(3,iel))-v(nodelist(4,iel))

            ! add hg restoring force to nodal components
            kap=kappareg(1)
            temp = -1.0*kap*rho(iel)*abs(area(iel))/max(dt,dtminhg)

            fx(nodelist(1,iel)) = fx(nodelist(1,iel)) + temp*ugam
            fy(nodelist(1,iel)) = fy(nodelist(1,iel)) + temp*vgam
            fx(nodelist(2,iel)) = fx(nodelist(2,iel)) - temp*ugam
            fy(nodelist(2,iel)) = fy(nodelist(2,iel)) - temp*vgam
            fx(nodelist(3,iel)) = fx(nodelist(3,iel)) + temp*ugam
            fy(nodelist(3,iel)) = fy(nodelist(3,iel)) + temp*vgam
            fx(nodelist(4,iel)) = fx(nodelist(4,iel)) - temp*ugam
            fy(nodelist(4,iel)) = fy(nodelist(4,iel)) - temp*vgam

            end do
        end if

        if (hgtyp == 2) then
        ! as used in dyna
        ! involves speed
        ! charcteristic length sqrt(area)
            do iel=1, nel
            ! calculate u  gamma
            !            ik     k
            ugam=u(nodelist(1,iel))-u(nodelist(2,iel))    &
                +u(nodelist(3,iel))-u(nodelist(4,iel))
            vgam=v(nodelist(1,iel))-v(nodelist(2,iel))    &
                +v(nodelist(3,iel))-v(nodelist(4,iel))

            ! add hg restoring force to nodal components
            kap=kappareg(1)
            temp = -1.0*kap*rho(iel)*sqrt(area(iel))*cc(iel)

            fx(nodelist(1,iel)) = fx(nodelist(1,iel)) + temp*ugam
            fy(nodelist(1,iel)) = fy(nodelist(1,iel)) + temp*vgam
            fx(nodelist(2,iel)) = fx(nodelist(2,iel)) - temp*ugam
            fy(nodelist(2,iel)) = fy(nodelist(2,iel)) - temp*vgam
            fx(nodelist(3,iel)) = fx(nodelist(3,iel)) + temp*ugam
            fy(nodelist(3,iel)) = fy(nodelist(3,iel)) + temp*vgam
            fx(nodelist(4,iel)) = fx(nodelist(4,iel)) - temp*ugam
            fy(nodelist(4,iel)) = fy(nodelist(4,iel)) - temp*vgam
            end do
        end if

        if (hgtyp == 1) then
        ! type 1 is artificial damping and stiffness as described
        ! by belystchko and flanagan
            do iel=1, nel

            a1=0.25_real64*(-x05(nodelist(1,iel))+x05(nodelist(2,iel))   &
                    +x05(nodelist(3,iel))-x05(nodelist(4,iel)))
            a2=0.25_real64*(x05(nodelist(1,iel))-x05(nodelist(2,iel))   &
                    +x05(nodelist(3,iel))-x05(nodelist(4,iel)))
            a3=0.25_real64*(-x05(nodelist(1,iel))-x05(nodelist(2,iel))   &
                    +x05(nodelist(3,iel))+x05(nodelist(4,iel)))
            b1=0.25_real64*(-y05(nodelist(1,iel))+y05(nodelist(2,iel))   &
                    +y05(nodelist(3,iel))-y05(nodelist(4,iel)))
            b2=0.25_real64*(y05(nodelist(1,iel))-y05(nodelist(2,iel))   &
                    +y05(nodelist(3,iel))-y05(nodelist(4,iel)))
            b3=0.25_real64*(-y05(nodelist(1,iel))-y05(nodelist(2,iel))   &
                    +y05(nodelist(3,iel))+y05(nodelist(4,iel)))
            biibii=4.0_real64*(b3**2+b1**2+a3**2+a1**2)
            xdiff=x05(nodelist(1,iel))-x05(nodelist(2,iel))    &
                    +x05(nodelist(3,iel))-x05(nodelist(4,iel))
            ydiff=y05(nodelist(1,iel))-y05(nodelist(2,iel))    &
                    +y05(nodelist(3,iel))-y05(nodelist(4,iel))
            gam1= 0.5_real64-0.5_real64*(dndx(1,iel)*xdiff+dndy(1,iel)*ydiff)/area(iel)
            gam2=-0.5_real64-0.5_real64*(dndx(2,iel)*xdiff+dndy(2,iel)*ydiff)/area(iel)
            gam3= 0.5_real64-0.5_real64*(dndx(3,iel)*xdiff+dndy(3,iel)*ydiff)/area(iel)
            gam4=-0.5_real64-0.5_real64*(dndx(4,iel)*xdiff+dndy(4,iel)*ydiff)/area(iel)

            ! qx qy different forms according to
            ! whether artificial damping or stiffness are required

            kap=kappareg(1)
            qx= -kap*rho(iel)*cc(iel)*sqrt(biibii)/10.0_real64  &!damping
                -0.5_real64*kap*rho(iel)*dt*(cc(iel)**2)*biibii/area(iel) !stiffness

            qy= qx*(gam1*v(nodelist(1,iel))+gam2*v(nodelist(2,iel))  &
                    +gam3*v(nodelist(3,iel))+gam4*v(nodelist(4,iel)))
            qx= qx*(gam1*u(nodelist(1,iel))+gam2*u(nodelist(2,iel))  &
                    +gam3*u(nodelist(3,iel))+gam4*u(nodelist(4,iel)))

            fx(nodelist(1,iel)) = fx(nodelist(1,iel)) + gam1*qx
            fy(nodelist(1,iel)) = fy(nodelist(1,iel)) + gam1*qy
            fx(nodelist(2,iel)) = fx(nodelist(2,iel)) + gam2*qx
            fy(nodelist(2,iel)) = fy(nodelist(2,iel)) + gam2*qy
            fx(nodelist(3,iel)) = fx(nodelist(3,iel)) + gam3*qx
            fy(nodelist(3,iel)) = fy(nodelist(3,iel)) + gam3*qy
            fx(nodelist(4,iel)) = fx(nodelist(4,iel)) + gam4*qx
            fy(nodelist(4,iel)) = fy(nodelist(4,iel)) + gam4*qy

            end do
        end if
    ! end do

end subroutine hourglass_filter


subroutine momentum_calculation( &
    dt, &
    zantihg, hgregtyp, kappareg, &
    u, v, x, y, rho, pressure, area, cc, q,  &  ! X, Y are 0.5_real64 timestep
    nint, dndx, dndy, &
    nodelist, znodbound, nel, nnod, &
    uout, vout &
)
    implicit none
    real(kind=real64), intent(in) :: dt
    integer(kind=int32), intent(in) :: zantihg
    integer(kind=int32), dimension(:), intent(in) :: hgregtyp
    real(kind=real64), dimension(:), intent(in) :: kappareg
    real(kind=real64), dimension(:), intent(in) :: u, v, x, y, rho, pressure
    real(kind=real64), dimension(:), intent(in) :: area, cc, q
    real(kind=real64), dimension(:,:), intent(in) :: nint, dndx, dndy
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    integer(kind=int32), dimension(:), intent(in) :: znodbound
    integer(kind=int32), intent(in) :: nel, nnod

    real(kind=real64), dimension(:), intent(out) :: uout, vout

    real(kind=real64), allocatable :: massnod(:) !node mass
    real(kind=real64), allocatable :: massj(:) !node mass el
    real(kind=real64), allocatable :: forcejx(:) !force mass el
    real(kind=real64), allocatable :: forcejy(:) !force mass el
    real(kind=real64), allocatable :: forcenodx(:) !force mass x
    real(kind=real64), allocatable :: forcenody(:) !force mass y

    integer(kind=int32) :: iel, inod, j

    allocate(massnod(1:nnod), massj(1:nnod), forcejx(1:nnod))
    allocate(forcejy(1:nnod))
    allocate(forcenodx(nnod))
    allocate(forcenody(nnod))

    massnod=0.0_real64
    forcenodx=0.0_real64
    forcenody=0.0_real64
    do inod=1,nnod
        do iel=1,nel
            do j=1,4
                ! this will calculate fine cells contrib to
                !disjoint nodes aswell
                if (nodelist(j,iel) == inod) then
                    massnod(inod) = massnod(inod) + rho(iel)*nint(j,iel)
                    forcenodx(inod) = forcenodx(inod) + (pressure(iel)+q(iel))*dndx(j,iel)
                    forcenody(inod) = forcenody(inod) + (pressure(iel)+q(iel))*dndy(j,iel)
                end if
            end do
        end do
    end do

    if (zantihg == 1) then
        call hourglass_filter( &
            u, v, x, y, rho, area, cc, &
            dt, dndx, dndy, hgregtyp, kappareg, nodelist, nel, &
            forcenodx, forcenody &
        )

    end if

    do inod=1,nnod
        uout(inod)=u(inod)+dt*forcenodx(inod)/massnod(inod)
        vout(inod)=v(inod)+dt*forcenody(inod)/massnod(inod)

        if ((znodbound(inod) == -1).or.(znodbound(inod) == -3)) then
            uout(inod)=u(inod)
        end if

        if ((znodbound(inod) == -2).or.(znodbound(inod) == -3)) then
            vout(inod)=v(inod)
        end if
    end do

end subroutine momentum_calculation

end module lagrangian_hydro
