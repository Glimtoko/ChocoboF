module sod_init
use globalconstants
use geom
implicit none

!        .. vectors ..
REAL (KIND=DP), PUBLIC :: xo(1:maxreg) ! x points region start
REAL (KIND=DP), PUBLIC :: yo(1:maxreg) ! y points region start
REAL (KIND=DP), PUBLIC :: xf(1:maxreg) !  x points region end
REAL (KIND=DP), PUBLIC :: yf(1:maxreg) !  x points region end
INTEGER, PUBLIC :: meshx(1:maxreg) ! x mesh for each region
INTEGER, PUBLIC :: meshy(1:maxreg) ! y mesh for each region

REAL (KIND=DP), PUBLIC, ALLOCATABLE :: deltax(:) !  x orig spacing/ region
REAL (KIND=DP), PUBLIC, ALLOCATABLE :: deltay(:) ! y orig spacing/ region
INTEGER, PUBLIC, ALLOCATABLE :: maxnod(:) ! mx node for each region
INTEGER, PUBLIC, ALLOCATABLE :: maxel(:) ! mx element for each region
INTEGER, PUBLIC :: uplow(1:maxreg)  !no el on top bottom boundary each reg
INTEGER, PUBLIC :: lefrig(1:maxreg)   !no el on left right boundary counters
INTEGER, PUBLIC :: low      ! boundary node counters
INTEGER, PUBLIC :: up
INTEGER, PUBLIC :: left
INTEGER, PUBLIC :: right

!     .. Module Scalars ..
INTEGER, PUBLIC    :: maxup,maxlef ! largest amount of boundary nodes over regions


contains
SUBROUTINE geominit()
! set up variables, allocate, set G coeffs, error coeffs
IMPLICIT NONE

!**********************
!input region pointers*
! and meshing         *
!**********************
NAMELIST /inp/xo,xf,yo,yf,meshx,meshy

OPEN (UNIT=212, FILE='param.dat')

READ(212,NML=inp)

WRITE(*,NML=inp)

! Set nreg
nreg = maxreg

! **************************
! ALLOCATE  ARRAYS         *
! **************************
ALLOCATE (deltax(1:maxreg), deltay(1:maxreg), maxnod(0:maxreg), maxel(0:maxreg))

!calculate deltax, deltay, nel, maxnod
nel=0
nnod=0
maxel=0
maxnod=0
uplow=0
lefrig=0
maxup=0
maxlef=0
Do ireg=1,nreg
    deltax(ireg) = (xf(ireg) - xo(ireg))/REAL(meshx(ireg), DP)
    deltay(ireg) = (yf(ireg) - yo(ireg))/REAL(meshy(ireg), DP)

    ! number nodes on each boundary
    uplow(ireg) = meshx(ireg) + 1
    lefrig(ireg) = meshy(ireg) + 1
    IF (uplow(ireg).GT.maxup) THEN
        maxup=uplow(ireg)
    END IF
    IF (lefrig(ireg).GT.maxlef) THEN
        maxlef=lefrig(ireg)
    END IF

    ! GENERAL ELEMENT AND NODE NO'S AT END EACH REGION
    maxel(ireg) = maxel(ireg-1) + meshx(ireg)*meshy(ireg)
    maxnod(ireg) = maxnod(ireg-1) + (meshx(ireg)+1)*(meshy(ireg)+1)

    ! number elements and nodes in whole problem
    nel = nel + meshx(ireg)*meshy(ireg)   ! max el whole grid
    nnod = nnod + (meshx(ireg)+1)*(meshy(ireg)+1)!max nod whole grid
END DO

write(*,*) 'maxlef, maxup', maxlef, maxup,nreg
! **************************
! ALLOCATE  ARRAYS         *
! **************************
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

! Output fils
OPEN (UNIT=21, FILE='results/time.txt')
OPEN (UNIT=27, FILE='results/toten.txt')
RETURN

END SUBROUTINE geominit


! ===================================================================
SUBROUTINE geomcalc()
! calculates x array in terms of nodes
IMPLICIT NONE

integer :: node1, node2, node3, node4

inod=1

Do ireg=1,nreg
    left=1
    right=1
    up=1
    low=1

    DO j=0,meshy(ireg)
        DO i=0,meshx(ireg)
            xv(inod) = xo(ireg) + i*deltax(ireg)
            yv(inod) = yo(ireg) + j*deltay(ireg)

            node1 = inod-j-maxnod(ireg-1)+maxel(ireg-1)
            node2 = inod-j-1-maxnod(ireg-1)+maxel(ireg-1)
            node3 = inod-j-1-meshx(ireg)-maxnod(ireg-1)+maxel(ireg-1)
            node4 = inod-j-meshx(ireg)-maxnod(ireg-1)+maxel(ireg-1)


            ! This part calculates
            ! element node connectivity array
            ! nodelist( gen el no, local nod no)=gen nod no

            !general node
            IF ((i /= 0).AND.(i /= meshx(ireg))) then
                IF ((j /= 0).AND.(j /= meshy(ireg))) then
                    nodelist(1,node1)=inod
                    nodelist(2,node2)=inod
                    nodelist(3,node3)=inod
                    nodelist(4,node4)=inod
                END IF
            END IF


            !boundary and corner nodes
            IF (j == 0) then
                IF ((i /= 0).AND.(i /= meshx(ireg))) then
                    nodelist(2,node2)=inod
                    nodelist(1,node1)=inod
                    low=low+1
                    znodbound(inod)=-2
                END IF
            END IF

            IF (j == meshy(ireg)) then
                IF((i /= 0).AND.(i /= meshx(ireg))) then
                    nodelist(3,node3)=inod
                    nodelist(4,node4)=inod
                    up=up+1
                    znodbound(inod)=-2
                END IF
            END IF

            IF (i == 0) then
                IF (j == 0) then  !corner node
                    nodelist(1,node1)=inod
                    low=low+1
                    left=left+1
                    znodbound(inod)=-3
                ELSE IF (j == meshy(ireg)) then !corner node
                    nodelist(4,node4)=inod
                    left=left+1
                    up=up+1
                    znodbound(inod)=-3
                ELSE   ! boundary node
                    nodelist(4,node4)=inod
                    nodelist(1,node1)=inod
                    left=left+1
                    znodbound(inod)=-1
                END IF
            END IF

            IF (i == meshx(ireg)) then
                IF (j == 0) then     !corner node
                    nodelist(2,node2)=inod
                    right=right+1
                    low=low+1
                    znodbound(inod)=-3
                ELSE IF (j == meshy(ireg)) then  !corner node
                    nodelist(3,node3)=inod
                    right=right+1
                    up=up+1
                    znodbound(inod)=-3
                ELSE   ! boundary node
                    nodelist(3,node3)=inod
                    nodelist(2,node2)=inod
                    right=right+1
                    znodbound(inod)=-1
                END IF
            END IF

            inod=inod+1
        END DO
    END DO
END DO

!need boundary arrays
RETURN
END SUBROUTINE geomcalc

end module
