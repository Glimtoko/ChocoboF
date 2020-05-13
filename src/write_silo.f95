module write_silo
use lagstep
use mesh_data
use geom_data

implicit none

include "silo_f9x.inc"

integer, save :: silo_count = 0

contains
subroutine write_silo_file(state, name)
! Arguments
integer, intent(in) :: state
character(len=*), intent(in) :: name

! Filename
integer :: stubsize
character(len=:), allocatable :: filename

! SILO stuff
integer :: dbID, mesh_optlistID
integer :: status, status2

! Internal stuff
integer :: i, j, k
integer, dimension(4*nel) :: connectivity

silo_count = silo_count + 1

stubsize = len(name)
allocate(character(stubsize+9) :: filename)

write (filename, "(A6,I4.4,A5)") name, silo_count, ".silo"

! Create a silo file
status = dbcreate(filename, len(filename), DB_CLOBBER, DB_LOCAL, "", 0, DB_HDF5, dbID)

print *, dbID

! Create a zone list
k = 1
do i = 1, nel
    do j = 1, 4
        connectivity(k) = nodelist(j,i)
        k = k + 1
    end do
end do

status = dbputzl2(dbID, "ZL", 2, nel, 2, connectivity, 4*nel, 1, 0, 0,  &
                  [DB_ZONETYPE_QUAD], [4], [nel], 1, DB_F77NULL, status2)

! Set mesh options
status = dbmkoptlist(2, mesh_optlistID)
status = dbadddopt (mesh_optlistID, DBOPT_DTIME, time)
status = dbaddiopt (mesh_optlistID, DBOPT_CYCLE, stepno)

! ! Dump a mesh
status = dbputum(dbID, "Mesh", 4, 2, xv, yv, 0, "X", 1, "Y", 1, "", 0, &
                 DB_DOUBLE, nnod, nel, "ZL", 2, DB_F77NULLSTRING, 0, &
                 mesh_optlistID, status2)

! Dump variables
status = dbputuv1(dbID, "Density", 7, "Mesh", 4, rho, nel, 0, 0, DB_DOUBLE, &
                 DB_ZONECENT, DB_F77NULL, status2)

status = dbputuv1(dbID, "Pressure", 8, "Mesh", 4, pre, nel, 0, 0, DB_DOUBLE, &
                 DB_ZONECENT, DB_F77NULL, status2)

status = dbputuv1(dbID, "Energy", 6, "Mesh", 4, en, nel, 0, 0, DB_DOUBLE, &
                 DB_ZONECENT, DB_F77NULL, status2)

status = dbputuv1(dbID, "X_Velocity", 10, "Mesh", 4, uv, nnod, 0, 0, DB_DOUBLE, &
                 DB_NODECENT, DB_F77NULL, status2)

status = dbputuv1(dbID, "Y_Velocity", 10, "Mesh", 4, vv, nnod, 0, 0, DB_DOUBLE, &
                 DB_NODECENT, DB_F77NULL, status2)

! Close the file
status = dbclose(dbID)

end subroutine

end module
