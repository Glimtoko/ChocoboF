
MODULE LAGSTEP
USE GlobalConstants
USE GEOM
! USE CONNECTARRAYS



IMPLICIT NONE

!        .. vectors ..
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::pre(:) !element pressures
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::rho(:)!element density
REAL (KIND=DP), PUBLIC, POINTER ::en(:)!element energy
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::cc(:) !element sound speed
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::qq(:) !element artificial viscosity
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::massel(:) !element mass
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::volel(:) !element volume
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::volelold(:) !element old vol
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::area(:) !element area
REAL (KIND=DP),DIMENSION(:),PUBLIC, POINTER ::uv !node x velocity
REAL (KIND=DP),DIMENSION(:),PUBLIC, POINTER ::vv!node y velocity

REAL (KIND=DP), PUBLIC, ALLOCATABLE ::uvold(:) !old x velocity
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::vvold(:) !old y velocity
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::uvbar(:) !av x velocity
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::vvbar(:) !av y velocity
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::xv05(:) !1/2dt node x
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::yv05(:) !1/2dt node y
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::pre05(:) !1/2dt pressures
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::rho05(:) !1/2dt density
REAL (KIND=DP), PUBLIC, POINTER ::en05(:)!1/2dt energy
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::volel05(:) !1/2dt volume
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::divint(:) !int divergence v
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::divvel(:) !int divergence v
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::dxtb(:),dxlr(:)
REAL (KIND=DP), PUBLIC,ALLOCATABLE ::dudxt(:),dudxb(:),dudxl(:),dudxr(:)
REAL (KIND=DP), PUBLIC,ALLOCATABLE ::phib(:), phit(:), phil(:), phir(:)
! used in momentum and hourgl
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::forcenodx(:) !force mass x
REAL (KIND=DP), PUBLIC, ALLOCATABLE ::forcenody(:) !force mass y

!fem
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::nint(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::elwtc(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::dndx(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::dndy(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::pdndx(:,:)
REAL (KIND=DP), DIMENSION(:,:),PUBLIC, ALLOCATABLE ::pdndy(:,:)
! cut off parameters
REAL(KIND=DP), PUBLIC, PARAMETER :: mindt=1.00000000000e-4
REAL(KIND=DP), PUBLIC, PARAMETER ::  zcut=1.00000000000e-8
REAL(KIND=DP), PUBLIC, PARAMETER :: dtminhg=0.000000008
REAL(KIND=DP), PUBLIC, PARAMETER ::  zerocut=1.0000000e-18
! stop divide by 0's
REAL(KIND=DP), PUBLIC, PARAMETER :: dencut=1.000000000e-6
! namelist tinp
! ideal gas gamma
REAL(KIND=DP), PUBLIC:: gamma
! artificial viscosity
REAL(KIND=DP), PUBLIC:: cq
REAL(KIND=DP), PUBLIC:: cl
!max time step
REAL(KIND=DP), PUBLIC:: maxallstep
! initial time step
REAL(KIND=DP), PUBLIC:: dtinit
INTEGER, PUBLIC :: dtoption ! timestep option
! time step groth factor
REAL(KIND=DP), PUBLIC :: t0  !initial time
REAL(KIND=DP), PUBLIC :: tf  !final time
REAL(KIND=DP), PUBLIC :: growth
! if 0 cartesian 1 axisymmetric
INTEGER, PUBLIC :: zaxis
INTEGER, PUBLIC :: zintdivvol   ! 0 not vol 1 indiv volume change
INTEGER, PUBLIC :: avtype ! type of artificial viscosity
INTEGER, PUBLIC :: zantihg  ! if 0 no hourglassing if 1 hourglassing
INTEGER, PUBLIC :: hgregtyp(1:maxreg)
REAL(KIND=DP), PUBLIC :: kappareg(1:maxreg)
!!!!!!!!!!!!!!end namelist tinp

REAL(KIND=DP), PUBLIC :: dt  !time step size
INTEGER, PUBLIC :: dtcontrol
REAL(KIND=DP), PUBLIC :: time  !time

REAL(KIND=DP), PUBLIC :: a1,a2,a3,b1,b2,b3
INTEGER, PUBLIC      ::prout, stepno
! total energy test
REAL(KIND=DP), PUBLIC :: totalenergy,totalke,totalie

CONTAINS
!===============================================
SUBROUTINE declare()
!allocate physical variables
IMPLICIT NONE
ALLOCATE (forcenodx(1:nnod),forcenody(1:nnod))
ALLOCATE (pre(1:nel),rho(1:nel),en(1:nel),cc(1:nel),qq(1:nel))
ALLOCATE (massel(1:nel),area(1:nel),volel(1:nel),volelold(1:nel))
ALLOCATE (uv(1:nnod),vv(1:nnod))
ALLOCATE (pre05(1:nel),rho05(1:nel),en05(1:nel),volel05(1:nel))
ALLOCATE (divint(1:nel),divvel(1:nel))
ALLOCATE (uvold(1:nnod),vvold(1:nnod),uvbar(1:nnod),vvbar(1:nnod))
ALLOCATE (xv05(1:nnod),yv05(1:nnod))
ALLOCATE (nint(1:4,1:nel),dndx(1:4,1:nel),dndy(1:4,1:nel))
ALLOCATE (elwtc(1:4,1:nel))
ALLOCATE (pdndx(1:4,1:nel),pdndy(1:4,1:nel))


! ALLOCATE (dxtb(1:nel),dxlr(1:nel))
! ALLOCATE (dudxb(1:nel),dudxt(1:nel),dudxl(1:nel),dudxr(1:nel))
! ALLOCATE (phib(1:nel),phit(1:nel),phil(1:nel),phir(1:nel))
RETURN
END SUBROUTINE declare

!========================================================
SUBROUTINE init()
!enter initial test problem configuration here
use sod_init
IMPLICIT NONE
REAL (KIND=DP)::   xbubble, ybubble
! SOD problem
pre=zero
rho=zero
uv=zero
vv=zero
! SOD circle quarter
pre=one/ten
rho=one/eight
uv=zero
vv=zero
Do iel=1,nel
    xbubble=xv(nodelist(1,iel))+xv(nodelist(2,iel))  &
            +xv(nodelist(3,iel))+xv(nodelist(4,iel))
    xbubble=xbubble/four
    ybubble=yv(nodelist(1,iel))+yv(nodelist(2,iel))  &
            +yv(nodelist(3,iel))+yv(nodelist(4,iel))
    ybubble=ybubble/four
    If ((xbubble)**2+(ybubble)**2 <= ((0.4_DP*(yf(1)-yo(1)))**2)+0.0000001_DP) then
        rho(iel)=one
        pre(iel)=one
    END IF
END DO

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
RETURN
END SUBROUTINE init



!=======================================================
SUBROUTINE volcalc(x,y,vol,area)
!calculate element area of quadralateral using jacobian
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  x(:),y(:)
REAL (KIND=DP), INTENT(OUT):: vol(:),area(:)

Do iel=1,nel
    a1=quarter*(-x(nodelist(1,iel))+x(nodelist(2,iel))   &
                +x(nodelist(3,iel))-x(nodelist(4,iel)))
    a2=quarter*(x(nodelist(1,iel))-x(nodelist(2,iel))   &
                +x(nodelist(3,iel))-x(nodelist(4,iel)))
    a3=quarter*(-x(nodelist(1,iel))-x(nodelist(2,iel))   &
                +x(nodelist(3,iel))+x(nodelist(4,iel)))
    b1=quarter*(-y(nodelist(1,iel))+y(nodelist(2,iel))   &
                +y(nodelist(3,iel))-y(nodelist(4,iel)))
    b2=quarter*(y(nodelist(1,iel))-y(nodelist(2,iel))   &
                +y(nodelist(3,iel))-y(nodelist(4,iel)))
    b3=quarter*(-y(nodelist(1,iel))-y(nodelist(2,iel))   &
                +y(nodelist(3,iel))+y(nodelist(4,iel)))
    vol(iel)=four*(a1*b3-a3*b1)
    area(iel)=vol(iel)
END DO

RETURN
END SUBROUTINE volcalc



!=======================================================
SUBROUTINE masscalc()
!calculate element mass=element density* element volume
IMPLICIT NONE
Do iel=1,nel
    massel(iel)=volel(iel)*rho(iel)
END DO
RETURN
END SUBROUTINE masscalc



!=======================================================
SUBROUTINE divv(u,v,pndx,pndy,diverg)
! not really divergence just a 2d mesure of gradient accross cell
! not axisymm
!calcs diverg for viscosity
! vsum rsum not used for artificial viscosity
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  u(:),v(:)
REAL (KIND=DP), INTENT(IN)::  pndx(:,:),pndy(:,:)
REAL (KIND=DP), INTENT(OUT):: diverg(:)
REAL (KIND=DP)             ::dterm
diverg=zero
DO iel=1,nel
    Do j=1,4
        dterm=u(nodelist(j,iel))*pndx(j,iel) + v(nodelist(j,iel))*pndy(j,iel)
        diverg(iel)=diverg(iel)+dterm
    END DO
END DO

RETURN
END SUBROUTINE divv



!=======================================================
! change 23rd June to add edge christensen
SUBROUTINE artificialviscosity(den,speed,diverg,viscosity)
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::diverg(:),den(:)
REAL (KIND=DP), INTENT(IN):: speed(:)
REAL (KIND=DP), INTENT(OUT):: viscosity(:)
REAL(KIND=DP)                    :: dudx

! Bulk q note our cl is Andy's/2
Do iel=1,nel
    IF(diverg(iel).LT.zero) THEN
        dudx=sqrt(area(iel))*diverg(iel)
        !dudx=u(nodelist(2,iel))-u(nodelist(1,iel))
        viscosity(iel)=(cq*den(iel)*(dudx)**2) +    &
        (cl*den(iel)*speed(iel)*ABS(dudx))
    ELSE
        viscosity(iel)=zero
    END IF
END DO
END SUBROUTINE artificialviscosity



!=======================================================
!=======================================================
SUBROUTINE soundspeed()
!calculate element sound speed
! ideal gas equation gamma is defined as parameter at init
IMPLICIT NONE
INTEGER :: jjjj
Do iel=1,nel
!     IF (rho(iel).LT.(0.0010000)) THEN
!         write(*,*)'rho small iel, rho, pre', iel, rho(iel), pre(iel)
!         Do jjjj=1,4
!             write(*,*)'u,v', uv(nodelist(jjjj,iel)),vv(nodelist(jjjj,iel))
!             write(*,*)'x,y', xv(nodelist(jjjj,iel)),yv(nodelist(jjjj,iel))
!         END DO
!     END IF
!
!     IF ((pre(iel).LT.(0.000000)).OR.(pre(iel).GT.(1000000.0))) THEN
!         write(*,*)'rho small iel, rho, pre', iel, rho(iel), pre(iel)
!         Do jjjj=1,4
!             write(*,*)'u,v', uv(nodelist(jjjj,iel)),vv(nodelist(jjjj,iel))
!             write(*,*)'x,y', xv(nodelist(jjjj,iel)),yv(nodelist(jjjj,iel))
!         END DO
!     END IF

    cc(iel)=sqrt(gamma*pre(iel)/rho(iel))
END DO
END SUBROUTINE soundspeed



!=======================================================
SUBROUTINE accn(step,x,y,u,v,xout,yout)
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  step,x(:),y(:),u(:),v(:)
REAL (KIND=DP), INTENT(OUT):: xout(:),yout(:)

Do inod=1,nnod
    xout(inod)=x(inod)+step*u(inod)
    yout(inod)=y(inod)+step*v(inod)
END DO
END SUBROUTINE accn



!=======================================================
SUBROUTINE femel(x,y,a1,a2,a3,b1,b2,b3,ni,ndx,ndy,pndx, pndy)
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  x(:),y(:)
REAL (KIND=DP), INTENT(OUT):: a1,a2,a3,b1,b2,b3
REAL (KIND=DP), INTENT(OUT):: ni(:,:),ndx(:,:),ndy(:,:)
REAL (KIND=DP), INTENT(OUT):: pndx(:,:),pndy(:,:)
REAL (KIND=DP)             ::jacob
INTEGER :: jjj

Do iel=1,nel
    a1=quarter*(-x(nodelist(1,iel))+x(nodelist(2,iel))   &
                +x(nodelist(3,iel))-x(nodelist(4,iel)))
    a2=quarter*(x(nodelist(1,iel))-x(nodelist(2,iel))   &
                +x(nodelist(3,iel))-x(nodelist(4,iel)))
    a3=quarter*(-x(nodelist(1,iel))-x(nodelist(2,iel))   &
                +x(nodelist(3,iel))+x(nodelist(4,iel)))
    b1=quarter*(-y(nodelist(1,iel))+y(nodelist(2,iel))   &
                +y(nodelist(3,iel))-y(nodelist(4,iel)))
    b2=quarter*(y(nodelist(1,iel))-y(nodelist(2,iel))   &
                +y(nodelist(3,iel))-y(nodelist(4,iel)))
    b3=quarter*(-y(nodelist(1,iel))-y(nodelist(2,iel))   &
                +y(nodelist(3,iel))+y(nodelist(4,iel)))
    ! for momentum
    ni(1,iel)=((three*b3-b2)*(three*a1-a2)-(three*a3-a2)*(three*b1-b2))/nine
    ni(2,iel)=((three*b3+b2)*(three*a1-a2)-(three*a3+a2)*(three*b1-b2))/nine
    ni(3,iel)=((three*b3+b2)*(three*a1+a2)-(three*a3+a2)*(three*b1+b2))/nine
    ni(4,iel)=((three*b3-b2)*(three*a1+a2)-(three*a3-a2)*(three*b1+b2))/nine

    ! integrated dndx's for intdiv and energy
    ndx(1,iel)=-b3+b1
    ndx(2,iel)= b3+b1
    ndx(3,iel)= b3-b1
    ndx(4,iel)=-b3-b1

    ndy(1,iel)= a3-a1
    ndy(2,iel)=-a3-a1
    ndy(3,iel)=-a3+a1
    ndy(4,iel)= a3+a1

    ! jacobian and dndx's for div and viscosity
    jacob=a1*b3-a3*b1
    Do jjj=1,4
        pndx(jjj,iel)=quarter*ndx(jjj,iel)/jacob
        pndy(jjj,iel)=quarter*ndy(jjj,iel)/jacob
    END DO

    ! calculate elwtc for SALE integral x over volume
    elwtc(1,iel)=ni(1,iel)
    elwtc(2,iel)=ni(2,iel)
    elwtc(3,iel)=ni(3,iel)
    elwtc(4,iel)=ni(4,iel)
    END DO
END SUBROUTINE femel



!=======================================================
SUBROUTINE density(mass,vol,den)
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  mass(:),vol(:)
REAL (KIND=DP), INTENT(OUT):: den(:)

Do iel=1,nel
    den(iel)=mass(iel)/vol(iel)
END DO
END SUBROUTINE density



!=======================================================
SUBROUTINE intdivv(step,vol,volold,u,v,ndx,ndy,intdiv)
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  u(:),v(:),step,vol(:),volold(:)
REAL (KIND=DP), INTENT(IN)::  ndx(:,:),ndy(:,:)
REAL (KIND=DP), INTENT(OUT):: intdiv(:)
REAL (KIND=DP)             ::eterm

IF (zintdivvol.EQ.0) THEN
    intdiv=zero
    DO iel=1,nel
        Do j=1,4
            eterm=u(nodelist(j,iel))*ndx(j,iel)         &
                    +v(nodelist(j,iel))*ndy(j,iel)
            intdiv(iel)=intdiv(iel)+eterm
        END DO
    END DO
ELSE
    stop 99
    Do iel=1,nel
        intdiv(iel)=(vol(iel)-volold(iel))/(step)
    END DO
END IF
END SUBROUTINE intdivv



!=======================================================
SUBROUTINE energy(step,press,visc,mass,ener,intdiv,enout)
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  step, press(:), visc(:)
REAL (KIND=DP), INTENT(IN)::  mass(:)
REAL (KIND=DP), POINTER :: ener(:)
REAL (KIND=DP), INTENT(IN)::  intdiv(:)
REAL (KIND=DP), POINTER :: enout(:)

DO iel=1,nel
    enout(iel)=ener(iel)-step*(press(iel)+visc(iel))*intdiv(iel)/mass(iel)
END DO
END SUBROUTINE energy



!=======================================================
SUBROUTINE hourglass(u,v,den,step,ndx,ndy,fx,fy)

! 26th june
!calculate anti-hourglass filters
!called in momentum subroutine

IMPLICIT NONE
REAL (KIND=DP), INTENT(INOUT)::fx(:),fy(:)
REAL (KIND=DP), INTENT(IN)::u(:),v(:),ndx(:,:),ndy(:,:)
REAL (KIND=DP), INTENT(IN)::den(:),step
REAL (KIND=DP)  :: ugam,vgam,temp,biibii,xdiff,ydiff
REAL (KIND=DP)  :: gam1,gam2,gam3,gam4,qx,qy, kap
INTEGER         :: hgtyp

! DO ireg=1,nreg   ! loop regions

  hgtyp=hgregtyp(1)

  IF (hgtyp.EQ.3) THEN
  ! DYNA hourglass type 2 with different l and no speed
    DO iel=1, nel
      ! calculate u  gamma
      !            ik     k
      ugam=u(nodelist(1,iel))-u(nodelist(2,iel))    &
          +u(nodelist(3,iel))-u(nodelist(4,iel))
      vgam=v(nodelist(1,iel))-v(nodelist(2,iel))    &
          +v(nodelist(3,iel))-v(nodelist(4,iel))

      ! add hg restoring force to nodal components
      kap=kappareg(1)
      temp = -1.0*kap*den(iel)*abs(area(iel))/max(step,dtminhg)

      fx(nodelist(1,iel)) = fx(nodelist(1,iel)) + temp*ugam
      fy(nodelist(1,iel)) = fy(nodelist(1,iel)) + temp*vgam
      fx(nodelist(2,iel)) = fx(nodelist(2,iel)) - temp*ugam
      fy(nodelist(2,iel)) = fy(nodelist(2,iel)) - temp*vgam
      fx(nodelist(3,iel)) = fx(nodelist(3,iel)) + temp*ugam
      fy(nodelist(3,iel)) = fy(nodelist(3,iel)) + temp*vgam
      fx(nodelist(4,iel)) = fx(nodelist(4,iel)) - temp*ugam
      fy(nodelist(4,iel)) = fy(nodelist(4,iel)) - temp*vgam

    END DO
  END IF

  IF (hgtyp.EQ.2) THEN
  ! as used in DYNA
  ! involves speed
  ! charcteristic length sqrt(area)
    DO iel=1, nel
      ! calculate u  gamma
      !            ik     k
      ugam=u(nodelist(1,iel))-u(nodelist(2,iel))    &
          +u(nodelist(3,iel))-u(nodelist(4,iel))
      vgam=v(nodelist(1,iel))-v(nodelist(2,iel))    &
          +v(nodelist(3,iel))-v(nodelist(4,iel))

      ! add hg restoring force to nodal components
      kap=kappareg(1)
      temp = -1.0*kap*den(iel)*sqrt(area(iel))*cc(iel)

      fx(nodelist(1,iel)) = fx(nodelist(1,iel)) + temp*ugam
      fy(nodelist(1,iel)) = fy(nodelist(1,iel)) + temp*vgam
      fx(nodelist(2,iel)) = fx(nodelist(2,iel)) - temp*ugam
      fy(nodelist(2,iel)) = fy(nodelist(2,iel)) - temp*vgam
      fx(nodelist(3,iel)) = fx(nodelist(3,iel)) + temp*ugam
      fy(nodelist(3,iel)) = fy(nodelist(3,iel)) + temp*vgam
      fx(nodelist(4,iel)) = fx(nodelist(4,iel)) - temp*ugam
      fy(nodelist(4,iel)) = fy(nodelist(4,iel)) - temp*vgam
    END DO
  END IF

  IF (hgtyp.EQ.1) Then
  ! type 1 is artificial damping and stiffness as described
  ! by Belystchko and Flanagan
    DO iel=1, nel

       a1=quarter*(-xv05(nodelist(1,iel))+xv05(nodelist(2,iel))   &
              +xv05(nodelist(3,iel))-xv05(nodelist(4,iel)))
       a2=quarter*(xv05(nodelist(1,iel))-xv05(nodelist(2,iel))   &
              +xv05(nodelist(3,iel))-xv05(nodelist(4,iel)))
       a3=quarter*(-xv05(nodelist(1,iel))-xv05(nodelist(2,iel))   &
              +xv05(nodelist(3,iel))+xv05(nodelist(4,iel)))
       b1=quarter*(-yv05(nodelist(1,iel))+yv05(nodelist(2,iel))   &
              +yv05(nodelist(3,iel))-yv05(nodelist(4,iel)))
       b2=quarter*(yv05(nodelist(1,iel))-yv05(nodelist(2,iel))   &
              +yv05(nodelist(3,iel))-yv05(nodelist(4,iel)))
       b3=quarter*(-yv05(nodelist(1,iel))-yv05(nodelist(2,iel))   &
              +yv05(nodelist(3,iel))+yv05(nodelist(4,iel)))
      biibii=four*(b3**2+b1**2+a3**2+a1**2)
      xdiff=xv05(nodelist(1,iel))-xv05(nodelist(2,iel))    &
            +xv05(nodelist(3,iel))-xv05(nodelist(4,iel))
      ydiff=yv05(nodelist(1,iel))-yv05(nodelist(2,iel))    &
            +yv05(nodelist(3,iel))-yv05(nodelist(4,iel))
      gam1= half-half*(ndx(1,iel)*xdiff+ndy(1,iel)*ydiff)/area(iel)
      gam2=-half-half*(ndx(2,iel)*xdiff+ndy(2,iel)*ydiff)/area(iel)
      gam3= half-half*(ndx(3,iel)*xdiff+ndy(3,iel)*ydiff)/area(iel)
      gam4=-half-half*(ndx(4,iel)*xdiff+ndy(4,iel)*ydiff)/area(iel)

      ! qx qy different forms according to
      ! whether artificial damping or stiffness are required

      kap=kappareg(1)
      qx= -kap*den(iel)*cc(iel)*sqrt(biibii)/ten  &!damping
          -half*kap*den(iel)*step*(cc(iel)**2)*biibii/area(iel) !stiffness

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

    END DO
  END IF
! END DO

end subroutine


!=======================================================
SUBROUTINE momentum(step,u,v,den, press,visc,ni, ndx, ndy, uout,vout)
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  step,press(:),den(:),visc(:),u(:),v(:)
REAL (KIND=DP), INTENT(IN)::  ni(:,:),ndx(:,:),ndy(:,:)
REAL (KIND=DP), INTENT(OUT):: uout(:),vout(:)
REAL (KIND=DP), ALLOCATABLE ::massnod(:) !node mass
REAL (KIND=DP), ALLOCATABLE ::massj(:) !node mass el
REAL (KIND=DP), ALLOCATABLE ::forcejx(:) !force mass el
REAL (KIND=DP), ALLOCATABLE ::forcejy(:) !force mass el
REAL (KIND=DP) :: temp

ALLOCATE (massnod(1:nnod),massj(1:nnod), forcejx(1:nnod))
ALLOCATE (forcejy(1:nnod))

massnod=zero
forcenodx=zero
forcenody=zero
Do inod=1,nnod
    Do iel=1,nel
        Do j=1,4
            ! this will calculate fine cells contrib to
            !disjoint nodes aswell
            IF (nodelist(j,iel).EQ.inod) THEN
                massnod(inod) = massnod(inod) + den(iel)*ni(j,iel)
                forcenodx(inod) = forcenodx(inod) + (press(iel)+visc(iel))*ndx(j,iel)
                forcenody(inod) = forcenody(inod) + (press(iel)+visc(iel))*ndy(j,iel)
            END IF
        END DO
    END DO
END DO

IF (zantihg.EQ.1) THEN
CALL hourglass(u,v,den,step,ndx,ndy,forcenodx,forcenody)
END IF

DO inod=1,nnod
    uout(inod)=u(inod)+step*forcenodx(inod)/massnod(inod)
    vout(inod)=v(inod)+step*forcenody(inod)/massnod(inod)

    IF ((znodbound(inod).eq.-1).OR.(znodbound(inod).eq.-3)) THEN
        uout(inod)=u(inod)
    END IF

    IF ((znodbound(inod).eq.-2).OR.(znodbound(inod).eq.-3)) THEN
        vout(inod)=v(inod)
    END IF
END DO

! print *, uout(12)

END SUBROUTINE momentum



!=======================================================
SUBROUTINE eos(ener,den,press)
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  ener(:),den(:)
REAL (KIND=DP), INTENT(OUT):: press(:)

Do iel=1,nel
    press(iel)=(gamma-one)*den(iel)*ener(iel)
END DO
END SUBROUTINE eos



!=======================================================
SUBROUTINE totalen(ener,den,u,v,toten,totke,totie)
IMPLICIT NONE
REAL (KIND=DP), INTENT(IN)::  ener(:),den(:)
REAL (KIND=DP), INTENT(IN)::  u(:),v(:)
REAL (KIND=DP), INTENT(OUT):: toten,totie,totke
REAL (KIND=DP)             :: elenergy
REAL (KIND=DP)             :: two_pi, tek

IF (zaxis.EQ.0) then
    two_pi=one
ELSE IF (zaxis.EQ.1) then
    two_pi=twopi
END IF

toten=zero
totke=zero
totie=zero
Do iel=1,nel
    tek=zero

    DO j=1,4
        tek= tek + half*den(iel)*elwtc(j,iel)*   &
            (u(nodelist(j,iel))**2+v(nodelist(j,iel))**2)*two_pi
    END DO
    elenergy=massel(iel)*ener(iel)*two_pi + tek
    totke=totke+tek
    totie=totie+ massel(iel)*ener(iel)*two_pi
    toten=toten+elenergy
END DO
END SUBROUTINE totalen
!=======================================================
!=======================================================
END MODULE LAGSTEP
!===============================================
