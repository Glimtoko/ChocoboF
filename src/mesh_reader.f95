module mesh_reader
use iso_fortran_env, only: int32, real64
implicit none

contains

subroutine get_mesh_size(problemname, nel, nnod, nreg)
character(len=*), intent(in) :: problemname
integer(kind=int32), intent(out) :: nel, nnod, nreg

character(len=:), allocatable :: meshfile
integer(kind=int32) :: status, unitnum

meshfile = trim(problemname) // "/mesh.chc"

! Open the mesh file
write(*,'("Reading headers from: ",a)') meshfile
open(newunit=unitnum, file=meshfile, status="OLD", iostat=status)
if (status /= 0) then
    write(*,*) "ERROR: File not found. Aborting"
    stop 98
end if

! Read the headers
read(unitnum, *, iostat=status) nnod
if (status /= 0) then
    write(*,*) "ERROR: Unable to read NNOD"
    stop 97
end if

read(unitnum, *, iostat=status) nel
if (status /= 0) then
    write(*,*) "ERROR: Unable to read NEL"
    stop 96
end if

read(unitnum, *, iostat=status) nreg
if (status /= 0) then
    write(*,*) "ERROR: Unable to read NREG"
    stop 95
end if

write(*,'("Number of nodes:",i5)') nnod
write(*,'("Number of cells:",i5)') nel
write(*,'("Number of regions:",i5)') nreg

close(unitnum)

end subroutine get_mesh_size


subroutine generate_mesh(problemname, nodelist, nodetype, region, material, regiontocell, xv, yv)
character(len=*), intent(in) :: problemname

integer(kind=int32), dimension(:,:), intent(out) :: nodelist
integer(kind=int32), dimension(:), intent(out) :: nodetype, region, material
integer(kind=int32), dimension(:,:), intent(out) :: regiontocell
real(kind=real64), dimension(:), intent(out) :: xv, yv

character(len=:), allocatable :: meshfile
integer(kind=int32) :: status, unitnum
integer(kind=int32) :: nel, nnod, nreg, i, d1, d2, d3

meshfile = trim(problemname) // "/mesh.chc"

! Open the mesh file
write(*,'("Reading mesh from: ",a)') meshfile
open(newunit=unitnum, file=meshfile, status="OLD", iostat=status)
if (status /= 0) then
    write(*,*) "ERROR: File not found. Aborting"
    stop 98
end if

! Read the headers
read(unitnum, *, iostat=status) nnod
read(unitnum, *, iostat=status) nel
read(unitnum, *, iostat=status) nreg

! Read coordinates
read(unitnum, *, iostat=status) xv
if (status /= 0) then
    write(*,*) "ERROR: Unable to read XV"
    stop 94
end if

read(unitnum, *, iostat=status) yv
if (status /= 0) then
    write(*,*) "ERROR: Unable to read YV"
    stop 93
end if

! Node type
read(unitnum, *, iostat=status) nodetype
if (status /= 0) then
    write(*,*) "ERROR: Unable to read NODETYPE"
    stop 92
end if


! Region and material arrays
read(unitnum, *, iostat=status) region
if (status /= 0) then
    write(*,*) "ERROR: Unable to read REGION"
    stop 91
end if

read(unitnum, *, iostat=status) material
if (status /= 0) then
    write(*,*) "ERROR: Unable to read MATERIAL"
    stop 90
end if

! Nodelist
read(unitnum, *, iostat=status) nodelist
if (status /= 0) then
    write(*,*) "ERROR: Unable to read NODELIST"
    stop 89
end if

! Region to cell iterator
do i = 1, nreg
    read(unitnum, *, iostat=status) d1, d2, d3
    if (status /= 0) then
        write(*,*) "ERROR: Unable to read RCI for region", i
        stop 89
    end if
    regiontocell(1,i) = d2
    regiontocell(2,i) = d2
end do

end subroutine


subroutine initial_conditions(problemname, nel, nmat, material, pre, rho, gamma)
implicit none
character(len=*), intent(in) :: problemname
integer(kind=int32), intent(in) :: nel, nmat
integer(kind=int32), dimension(:), intent(in) :: material

real(kind=real64), dimension(:), intent(out) :: pre, rho, gamma

character(len=:), allocatable :: eosfile
integer(kind=int32) :: status, unitnum, i, j
real(kind=real64) :: rhomat, premat, gammamat

eosfile = trim(problemname) // "/eos.dat"

! Open the mesh file
write(*,'("Reading EoS data from: ",a)') eosfile
open(newunit=unitnum, file=eosfile, status="OLD", iostat=status)
if (status /= 0) then
    write(*,*) "ERROR: File not found. Aborting"
    stop 88
end if

do i = 1, nmat
    read(unitnum, *, iostat=status) rhomat, premat, gammamat
    if (status /= 0) then
        write(*,*) "ERROR: Unable to read data for material", i
        stop 87
    end if
    write(*,'("Material ",i2)') i
    write(*,'("Initial density:",f8.4)') rhomat
    write(*,'("Initial pressure:",f8.4)') premat
    write(*,'("Gamma:",f8.4/)') gammamat

    do j = 1, nel
        if (material(j) == i) then
            rho(j) = rhomat
            pre(j) = premat
        end if
    end do
    gamma(i) = gammamat
end do



end subroutine

end module mesh_reader
