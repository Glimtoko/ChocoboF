module sod_init
use iso_fortran_env, only: int32, real64
use core_input, only: maxreg
implicit none

! User inputs
real (kind=real64), private :: xo(1:maxreg) ! x points region start
real (kind=real64), private :: yo(1:maxreg) ! y points region start
real (kind=real64), private :: xf(1:maxreg) !  x points region end
real (kind=real64), private :: yf(1:maxreg) !  x points region end
integer(kind=int32), private :: meshx(1:maxreg) ! x mesh for each region
integer(kind=int32), private :: meshy(1:maxreg) ! y mesh for each region

real (kind=real64), private, allocatable :: deltax(:) !  x orig spacing/ region
real (kind=real64), private, allocatable :: deltay(:) ! y orig spacing/ region
integer(kind=int32), private, allocatable :: maxnod(:) ! mx node for each region
integer(kind=int32), private, allocatable :: maxel(:) ! mx element for each region



contains

subroutine geominit(nel, nnod, nreg)
    ! set up variables, allocate, set g coeffs, error coeffs
    implicit none
    integer(kind=int32), intent(out) :: nel, nnod, nreg

    integer(kind=int32) :: uplow(maxreg)  !no el on top bottom boundary each reg
    integer(kind=int32) :: lefrig(maxreg)   !no el on left right boundary counters
    integer(kind=int32) :: maxup, maxlef ! largest amount of boundary nodes over regions
    integer(kind=int32) :: ireg
    !**********************
    !input region pointers*
    ! and meshing         *
    !**********************
    namelist /inp/xo,xf,yo,yf,meshx,meshy



    open (unit=212, file='param.dat')

    read(212,nml=inp)

    write(*,nml=inp)

    ! set nreg
    nreg = 1

    ! **************************
    ! allocate  arrays         *
    ! **************************
    allocate (deltax(1:maxreg), deltay(1:maxreg), maxnod(0:maxreg), maxel(0:maxreg))

    ! global values
    nel = 0
    nnod = 0

    ! local values
    maxel = 0
    maxnod = 0
    uplow = 0
    lefrig = 0
    maxup = 0
    maxlef = 0
    do ireg = 1,nreg
        deltax(ireg) = (xf(ireg) - xo(ireg))/real(meshx(ireg), real64)
        deltay(ireg) = (yf(ireg) - yo(ireg))/real(meshy(ireg), real64)

        ! number nodes on each boundary
        uplow(ireg) = meshx(ireg) + 1
        lefrig(ireg) = meshy(ireg) + 1
        if (uplow(ireg) > maxup) then
            maxup = uplow(ireg)
        end if
        if (lefrig(ireg) > maxlef) then
            maxlef = lefrig(ireg)
        end if

        ! general element and node no's at end each region
        maxel(ireg) = maxel(ireg-1) + meshx(ireg)*meshy(ireg)
        maxnod(ireg) = maxnod(ireg-1) + (meshx(ireg)+1)*(meshy(ireg)+1)

        ! number elements and nodes in whole problem
        nel = nel + meshx(ireg)*meshy(ireg)   ! max el whole grid
        nnod = nnod + (meshx(ireg)+1)*(meshy(ireg)+1)!max nod whole grid
    end do

    write(*,*) 'maxlef, maxup', maxlef, maxup,nreg


    ! output files
    open (unit=21, file='results/time.txt')
    open (unit=27, file='results/toten.txt')
    return

END SUBROUTINE geominit


! ===================================================================
subroutine geomcalc(nreg, nodelist, znodbound, xv, yv)

    ! calculates x array in terms of nodes
    implicit none

    integer(kind=int32), intent(in) :: nreg
    integer(kind=int32), dimension(:,:), intent(out) :: nodelist
    integer(kind=int32), dimension(:), intent(out) :: znodbound
    real(kind=real64), dimension(:), intent(out) :: xv, yv

    integer(kind=int32) :: low      ! boundary node counters
    integer(kind=int32) :: up
    integer(kind=int32) :: left
    integer(kind=int32) :: right

    integer(kind=int32) :: i, j, inod, ireg
    integer(kind=int32) :: node1, node2, node3, node4

    inod=1

    do ireg=1,nreg
        left=1
        right=1
        up=1
        low=1

        do j=0,meshy(ireg)
            do i=0,meshx(ireg)
                xv(inod) = xo(ireg) + i*deltax(ireg)
                yv(inod) = yo(ireg) + j*deltay(ireg)

                node1 = inod-j-maxnod(ireg-1)+maxel(ireg-1)
                node2 = inod-j-1-maxnod(ireg-1)+maxel(ireg-1)
                node3 = inod-j-1-meshx(ireg)-maxnod(ireg-1)+maxel(ireg-1)
                node4 = inod-j-meshx(ireg)-maxnod(ireg-1)+maxel(ireg-1)


                ! this part calculates
                ! element node connectivity array
                ! nodelist( gen el no, local nod no)=gen nod no

                !general node
                if ((i /= 0) .and. (i /= meshx(ireg))) then
                    if ((j /= 0) .and. (j /= meshy(ireg))) then
                        nodelist(1,node1)=inod
                        nodelist(2,node2)=inod
                        nodelist(3,node3)=inod
                        nodelist(4,node4)=inod
                    end if
                end if


                !boundary and corner nodes
                if (j == 0) then
                    if ((i /= 0) .and. (i /= meshx(ireg))) then
                        nodelist(2,node2)=inod
                        nodelist(1,node1)=inod
                        low=low+1
                        znodbound(inod)=-2
                    end if
                end if

                if (j == meshy(ireg)) then
                    if((i /= 0) .and. (i /= meshx(ireg))) then
                        nodelist(3,node3)=inod
                        nodelist(4,node4)=inod
                        up=up+1
                        znodbound(inod)=-2
                    end if
                end if

                if (i == 0) then
                    if (j == 0) then  !corner node
                        nodelist(1,node1)=inod
                        low=low+1
                        left=left+1
                        znodbound(inod)=-3
                    else if (j == meshy(ireg)) then !corner node
                        nodelist(4,node4)=inod
                        left=left+1
                        up=up+1
                        znodbound(inod)=-3
                    else   ! boundary node
                        nodelist(4,node4)=inod
                        nodelist(1,node1)=inod
                        left=left+1
                        znodbound(inod)=-1
                    end if
                end if

                if (i == meshx(ireg)) then
                    if (j == 0) then     !corner node
                        nodelist(2,node2)=inod
                        right=right+1
                        low=low+1
                        znodbound(inod)=-3
                    else if (j == meshy(ireg)) then  !corner node
                        nodelist(3,node3)=inod
                        right=right+1
                        up=up+1
                        znodbound(inod)=-3
                    else   ! boundary node
                        nodelist(3,node3)=inod
                        nodelist(2,node2)=inod
                        right=right+1
                        znodbound(inod)=-1
                    end if
                end if

                inod=inod+1
            end do
        end do
    end do
end subroutine geomcalc


subroutine init(nel, nodelist, xv, yv, pre, rho, uv, vv)

    implicit none
    integer(kind=int32), intent(in) :: nel
    integer(kind=int32), dimension(:,:), intent(in) :: nodelist
    real(kind=real64), dimension(:), intent(in) :: xv, yv

    real(kind=real64), dimension(:), intent(out) :: pre, rho, uv, vv

    real (kind=real64) ::   xbubble, ybubble
    integer(kind=int32) :: iel


    ! sod circle quarter
    pre = 1.0_real64/10.0_real64
    rho = 1.0_real64/8.0_real64
!     uv=zero
!     vv=zero
    do iel=1,nel
        xbubble=xv(nodelist(1,iel))+xv(nodelist(2,iel))  &
                +xv(nodelist(3,iel))+xv(nodelist(4,iel))
        xbubble=xbubble/4.0_real64
        ybubble=yv(nodelist(1,iel))+yv(nodelist(2,iel))  &
                +yv(nodelist(3,iel))+yv(nodelist(4,iel))
        ybubble=ybubble/4.0_real64
        if ((xbubble)**2+(ybubble)**2 <= ((0.4_real64*(yf(1)-yo(1)))**2)+0.0000001_real64) then
            rho(iel) = 1.0_real64
            pre(iel) = 1.0_real64
        end if
    end do


end subroutine init


end module
