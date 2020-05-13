
! ===================================================================

MODULE GRAPHICS
! Simple ALE, reqional mesh relaxation.
! Mesh movement done using Winslow's equipotential solver
! Solution quantities advected using the ideas
! of Van Leer and Benson
USE GLOBALCONSTANTS
USE geom_data
use write_tio
USE LAGSTEP


IMPLICIT NONE
CHARACTER(LEN=10)::numbers
INTEGER::test_number,hundreds,tens,units
CHARACTER(LEN=32)::filename
CHARACTER(LEN=32)::filetwo
CHARACTER(LEN=32)::filethree
CHARACTER(LEN=32)::filefour
CHARACTER(LEN=32)::filefive
CHARACTER(LEN=32)::filesix
CHARACTER(LEN=32)::fileseven
CHARACTER(LEN=32)::fileeight
CHARACTER(LEN=32)::fileeigtwo
CHARACTER(LEN=32)::filenine
REAL(KIND=DP) ::xmidel,ymidel




CONTAINS
!===================================================================
SUBROUTINE output(stepcnt)
IMPLICIT NONE
integer :: stepcnt

numbers="0123456789"
filename="results/mypre_no"
filetwo="results/myvel_no"
filethree="results/myx_no"
filefour="results/nodnod_no"
fileseven="results/nodelist_no"
fileeigtwo="results/vol_no"

IF (time.GE.0.20) then
    test_number=prout
    hundreds=test_number/100
    test_number=test_number-100*hundreds
    tens=test_number/10
    test_number=test_number-10*tens
    units=test_number

!     write(*,*) "Output #",numbers(hundreds+1:hundreds+1)//numbers(tens+1:tens+1)// &
!         numbers(units+1:units+1), " @ ", time

    OPEN(UNIT=11,FILE=trim(fileseven)//numbers(hundreds+1:hundreds+1)//numbers(tens+1:tens+1)// &
        numbers(units+1:units+1)//".txt" )

    Do iel=1,nel
        WRITE(11,FMT='(I6,2X,I6,2X,I6,2X,I6,2X )') &
        nodelist(1,iel),nodelist(2,iel),nodelist(3,iel),nodelist(4,iel)
    END DO
    CLOSE(UNIT=11)


    OPEN(UNIT=12,FILE=trim(filename)//numbers(hundreds+1:hundreds+1)//numbers(tens+1:tens+1)// &
        numbers(units+1:units+1)//".txt" )
    OPEN(UNIT=82,FILE=trim(fileeigtwo)//numbers(hundreds+1:hundreds+1)//numbers(tens+1:tens+1)// &
        numbers(units+1:units+1)//".txt" )
    Do iel=1,nel
        WRITE(12,FMT='(I6,2X,F17.12,2X,F17.12,2x,F19.12 )') iel,pre(iel),rho(iel),en(iel)
        WRITE(82,FMT='(I6,2X,F17.12 )') iel,volel(iel)
    END DO
    CLOSE(UNIT=12)
    CLOSE(UNIT=82)


    OPEN(UNIT=13,FILE=trim(filetwo)//numbers(hundreds+1:hundreds+1)//numbers(tens+1:tens+1)// &
        numbers(units+1:units+1)//".txt" )
    Do inod=1,nnod
        WRITE(13,FMT='(F17.12,2x,F17.12,2x,F17.12 )') time,uv(inod),vv(inod)
    END DO
    CLOSE(UNIT=13)


    OPEN(UNIT=14,FILE=trim(filethree)//numbers(hundreds+1:hundreds+1)//numbers(tens+1:tens+1)// &
        numbers(units+1:units+1)//".txt" )
    Do inod=1,nnod
        WRITE(14,FMT='(F17.12,2x,F17.12 )') xv(inod),yv(inod)
    END DO
    CLOSE(UNIT=14)

    prout=prout+1
END IF


RETURN
END SUBROUTINE output
!+++++++++++++++++++++++
END MODULE GRAPHICS
