
! ===================================================================
MODULE GEOM
!constructing problem geometry, x,y vectors, lement-node array.
USE GlobalConstants

IMPLICIT NONE
INTEGER, PUBLIC, PARAMETER :: maxreg=10

REAL (KIND=DP), PUBLIC, ALLOCATABLE :: xv(:) ! nodal x points node order
REAL (KIND=DP), PUBLIC, ALLOCATABLE :: yv(:) ! nodal y points node order

!     .. Module Scalars ..
integer, public    :: nreg
INTEGER, PUBLIC    :: nel    !max no. elements
INTEGER, PUBLIC    :: nnod  !max no. nodes

! "logicals" - actually integers
INTEGER, PUBLIC, ALLOCATABLE :: znodbound(:)

!       ..arrays..
INTEGER, DIMENSION(:,:),PUBLIC, ALLOCATABLE :: nodelist !element-node array

INTEGER, PUBLIC    :: ireg    ! region number
INTEGER, PUBLIC    :: inod    !node number
INTEGER, PUBLIC    :: iel    !element number
INTEGER :: i,j


CONTAINS

!=================================================

! SUBROUTINE geominit()
! ! set up variables, allocate, set G coeffs, error coeffs
! IMPLICIT NONE
!
! !**********************
! !input region pointers*
! ! and meshing         *
! !**********************
! NAMELIST /inp/xo,xf,yo,yf,meshx,meshy
!
! OPEN (UNIT=212, FILE='param.dat')
!
! READ(212,NML=inp)
!
! WRITE(*,NML=inp)
!
! ! **************************
! ! ALLOCATE  ARRAYS         *
! ! **************************
! ALLOCATE (deltax(1:nreg), deltay(1:nreg), maxnod(0:nreg), maxel(0:nreg))
!
! !calculate deltax, deltay, nel, maxnod
! nel=0
! nnod=0
! maxel=0
! maxnod=0
! uplow=0
! lefrig=0
! maxup=0
! maxlef=0
! Do ireg=1,nreg
!     deltax(ireg) = (xf(ireg) - xo(ireg))/REAL(meshx(ireg), DP)
!     deltay(ireg) = (yf(ireg) - yo(ireg))/REAL(meshy(ireg), DP)
!
!     ! number nodes on each boundary
!     uplow(ireg) = meshx(ireg) + 1
!     lefrig(ireg) = meshy(ireg) + 1
!     IF (uplow(ireg).GT.maxup) THEN
!         maxup=uplow(ireg)
!     END IF
!     IF (lefrig(ireg).GT.maxlef) THEN
!         maxlef=lefrig(ireg)
!     END IF
!
!     ! GENERAL ELEMENT AND NODE NO'S AT END EACH REGION
!     maxel(ireg) = maxel(ireg-1) + meshx(ireg)*meshy(ireg)
!     maxnod(ireg) = maxnod(ireg-1) + (meshx(ireg)+1)*(meshy(ireg)+1)
!
!     ! number elements and nodes in whole problem
!     nel = nel + meshx(ireg)*meshy(ireg)   ! max el whole grid
!     nnod = nnod + (meshx(ireg)+1)*(meshy(ireg)+1)!max nod whole grid
! END DO
!
! write(*,*) 'maxlef, maxup', maxlef, maxup,nreg
! ! **************************
! ! ALLOCATE  ARRAYS         *
! ! **************************
! ALLOCATE (xv(1:nnod), yv(1:nnod))
! ALLOCATE (znodbound(1:nnod))
! ALLOCATE (nodelist(1:4,1:nel))
! !*******************************************
! !initialise scalars, vectors, matrices     *
! !*******************************************
! xv=zero
! yv=zero
! nodelist=0
!
! ! set logicals to 0
! znodbound=0
!
! ! Output fils
! OPEN (UNIT=21, FILE='results/time.txt')
! OPEN (UNIT=27, FILE='results/toten.txt')
! RETURN
!
! END SUBROUTINE geominit
!
!
! ! ===================================================================
! SUBROUTINE geomcalc()
! ! calculates x array in terms of nodes
! IMPLICIT NONE
!
! integer :: node1, node2, node3, node4
!
! inod=1
!
! Do ireg=1,nreg
!     left=1
!     right=1
!     up=1
!     low=1
!
!     DO j=0,meshy(ireg)
!         DO i=0,meshx(ireg)
!             xv(inod) = xo(ireg) + i*deltax(ireg)
!             yv(inod) = yo(ireg) + j*deltay(ireg)
!
!             node1 = inod-j-maxnod(ireg-1)+maxel(ireg-1)
!             node2 = inod-j-1-maxnod(ireg-1)+maxel(ireg-1)
!             node3 = inod-j-1-meshx(ireg)-maxnod(ireg-1)+maxel(ireg-1)
!             node4 = inod-j-meshx(ireg)-maxnod(ireg-1)+maxel(ireg-1)
!
!
!             ! This part calculates
!             ! element node connectivity array
!             ! nodelist( gen el no, local nod no)=gen nod no
!
!             !general node
!             IF ((i /= 0).AND.(i /= meshx(ireg))) then
!                 IF ((j /= 0).AND.(j /= meshy(ireg))) then
!                     nodelist(1,node1)=inod
!                     nodelist(2,node2)=inod
!                     nodelist(3,node3)=inod
!                     nodelist(4,node4)=inod
!                 END IF
!             END IF
!
!
!             !boundary and corner nodes
!             IF (j == 0) then
!                 IF ((i /= 0).AND.(i /= meshx(ireg))) then
!                     nodelist(2,node2)=inod
!                     nodelist(1,node1)=inod
!                     low=low+1
!                     znodbound(inod)=-2
!                 END IF
!             END IF
!
!             IF (j == meshy(ireg)) then
!                 IF((i /= 0).AND.(i /= meshx(ireg))) then
!                     nodelist(3,node3)=inod
!                     nodelist(4,node4)=inod
!                     up=up+1
!                     znodbound(inod)=-2
!                 END IF
!             END IF
!
!             IF (i == 0) then
!                 IF (j == 0) then  !corner node
!                     nodelist(1,node1)=inod
!                     low=low+1
!                     left=left+1
!                     znodbound(inod)=-3
!                 ELSE IF (j == meshy(ireg)) then !corner node
!                     nodelist(4,node4)=inod
!                     left=left+1
!                     up=up+1
!                     znodbound(inod)=-3
!                 ELSE   ! boundary node
!                     nodelist(4,node4)=inod
!                     nodelist(1,node1)=inod
!                     left=left+1
!                     znodbound(inod)=-1
!                 END IF
!             END IF
!
!             IF (i == meshx(ireg)) then
!                 IF (j == 0) then     !corner node
!                     nodelist(2,node2)=inod
!                     right=right+1
!                     low=low+1
!                     znodbound(inod)=-3
!                 ELSE IF (j == meshy(ireg)) then  !corner node
!                     nodelist(3,node3)=inod
!                     right=right+1
!                     up=up+1
!                     znodbound(inod)=-3
!                 ELSE   ! boundary node
!                     nodelist(3,node3)=inod
!                     nodelist(2,node2)=inod
!                     right=right+1
!                     znodbound(inod)=-1
!                 END IF
!             END IF
!
!             inod=inod+1
!         END DO
!     END DO
! END DO
!
! !need boundary arrays
! RETURN
! END SUBROUTINE geomcalc
!=================================================
END MODULE GEOM
