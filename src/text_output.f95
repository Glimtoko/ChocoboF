module text_output
use iso_fortran_env, only: int32, real64
implicit none

contains
subroutine output(problemname, stepcnt, nel, nnod, nodelist, time, stepno, xv, yv, rho, pre, en, uv, vv, volel, prout)
    implicit none
    character(len=*), intent(in) :: problemname
    integer(kind=int32), intent(in) :: stepcnt
    integer(kind=int32), intent(in) :: nel, nnod
    integer(kind=int32), dimension(:,:), intent(in) ::nodelist
    real(kind=real64), intent(in) :: time
    integer(kind=int32), intent(in) :: stepno
    real(kind=real64), dimension(:), intent(in) :: xv, yv, rho, pre, en, uv, vv, volel
    integer(kind=int32), intent(inout) :: prout

    integer(kind=int32) :: iel, inod
    integer(kind=int32) :: unitnum, unitnum2

    character(len=40) :: pressure_f
    character(len=40) :: velocity_f
    character(len=40) :: coord_f
    character(len=40) :: nodelist_f
    character(len=40) :: volume_f
    character(len=4) :: ext = ".txt"

    character(len=3) :: filenumber

    pressure_f = trim(problemname)//"/results/mypre_no"
    velocity_f = trim(problemname)//"/results/myvel_no"
    coord_f = trim(problemname)//"/results/myx_no"
    nodelist_f = trim(problemname)//"/results/nodelist_no"
    volume_f = trim(problemname)//"/results/vol_no"

    ! Generate the digits of the file number
    write(filenumber,"(I3.3)") prout

    ! Nodelist
    open(newunit=unitnum, file=trim(nodelist_f)//filenumber//ext)
    do iel=1,nel
        write(unitnum, fmt='(i6,2x,i6,2x,i6,2x,i6,2x )') &
            nodelist(1,iel),nodelist(2,iel),nodelist(3,iel),nodelist(4,iel)
    end do
    close(unit=unitnum)

    ! Cell-centred quantities
    open(newunit=unitnum, file=trim(pressure_f)//filenumber//ext)
    open(newunit=unitnum2, file=trim(volume_f)//filenumber//ext)
    do iel=1,nel
        write(unitnum, fmt='(i6,2x,f17.12,2x,f17.12,2x,f19.12 )') &
            iel,pre(iel),rho(iel),en(iel)
        write(unitnum2, fmt='(i6,2x,f17.12 )') iel,volel(iel)
    end do
    close(unit=unitnum)
    close(unit=unitnum2)

    ! Velocities
    open(newunit=unitnum, file=trim(velocity_f)//filenumber//ext)
    do inod=1,nnod
        write(unitnum, fmt='(f17.12,2x,f17.12,2x,f17.12 )') &
            time,uv(inod),vv(inod)
    end do
    close(unit=unitnum)

    ! Coordinates
    open(newunit=unitnum, file=trim(coord_f)//filenumber//ext)
    do inod=1,nnod
        write(unitnum, fmt='(f17.12,2x,f17.12 )') xv(inod),yv(inod)
    end do
    close(unit=unitnum)

    prout=prout+1
end subroutine output

end module text_output
