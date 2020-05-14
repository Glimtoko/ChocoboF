module write_tio
use iso_fortran_env, only: int32, real64
use globalconstants
use core_input

use typhonio
implicit none

integer, save :: tio_count = 0

contains
subroutine write_tio_file(state, name, nel, nnod, nodelist, time, stepno, xv, yv, rho, pre, en, uv, vv)
! Arguments
integer(kind=int32), intent(in) :: state
character(len=*), intent(in) :: name
integer(kind=int32), intent(in) :: nel, nnod
integer(kind=int32), dimension(:,:), intent(in) ::nodelist
real(kind=real64), intent(in) :: time
integer(kind=int32), intent(in) :: stepno
real(kind=real64), dimension(:), intent(in) :: xv, yv, rho, pre, en, uv, vv

! Internal Stuff
integer :: i, j, k

! TIO Stuff
integer(kind=TIO_FILEK) :: fileID
integer(kind=TIO_OBJK) :: stateID
integer(kind=TIO_OBJK) :: meshID
integer(kind=TIO_ERRK) :: status

! Arguments to TIO_Create_f
integer :: stubsize
character(len=:), allocatable :: filename
character(len=:), allocatable :: code, version, date, title

! Arguments to TIO_Create_State_f
character(len=10) :: state_name  ! state_0001 - state_9999
character(len=:), allocatable :: units
integer(kind=TIO_STEPK) :: step_num
real(kind=TIO_TIMEK) :: step_time

! Arguments to TIO_Create_Mesh
character(len=:), allocatable :: mesh_name, group
integer(kind=TIO_SIZEK) :: order, n1, n2, n3, n4, nchunks
character(len=:), allocatable :: iunits, junits

! Arguments to TIO_Write_UnstrMesh_Chunk_f
integer(kind=TIO_SHAPEK), dimension(1) :: shapes = [TIO_SHAPE_QUAD4_F]
integer, dimension(1) :: ncells_per_shape
integer, dimension(nnod) :: nodeIDs
integer, dimension(nel) :: cellIDs
integer, dimension(4*nel) :: connectivity


! Set filename
if (tioonefile) then
    filename = name // "_all.h5"
else
    tio_count = tio_count + 1
    stubsize = len(name)
    allocate(character(stubsize+7) :: filename)

    write (filename, "(A6,I4.4,A3)") name, tio_count, ".h5"
end if


! Set code details
code = "gascode"
version = "0.1"
date = "N/A"
title = "Main Dump"

! Create or open state
if (state == 1 .or. .not. tioonefile) then
    ! First state, create a new file

    ! Create the file
    status = TIO_Create_f( &
        filename, &
        fileID, &
        TIO_ACC_REPLACE_F, &
        code, &
        version, &
        date, &
        title &
    )
else
    ! Open the file
    status = TIO_Open_f( &
        filename, &
        fileID, &
        TIO_ACC_READWRITE_F, &
        code, &
        version, &
        date, &
        title &
    )
    print *, code
    print *, title
end if


! Create a new state
write (state_name, "(A6,I4.4)") "state_", state
step_num = stepno
step_time = time
units = "s"

status = TIO_Create_State_f( &
    fileID, &
    state_name, &
    stateID, &
    step_num, &
    step_time, &
    units &
)

! Create a mesh
mesh_name = "Mesh"
group = "Dummy"
order = 1
n1 = nnod
n2 = nel
n3 = 1
n4 = 4*nel
nchunks = 1

status = TIO_Create_Mesh_f( &
    fileID = fileID, &
    stateID = stateID, &
    name = mesh_name, &
    meshID = meshID, &
    meshtype = TIO_MESH_UNSTRUCT_F, &
    coordtype = TIO_COORD_CARTESIAN_F, &
    isAMR = .false., &
    group = group, &
    order = order, &
    graph_datatype = TIO_INTEGER4_F, &
    coord_datatype = TIO_REAL8_F, &
    ndims = TIO_2D_F, &
    n1 = n1, &
    n2 = n2, &
    n3 = n3, &
    n4 = n4, &
    nchunks = nchunks, &
    iunits = iunits, &
    junits = junits &
)
print *, "MESH", status

status = TIO_Set_Unstr_Chunk_f(fileID, meshID, 1, 2, &
                      n1, n2, &
                      n3, n4, &
                      0,  0, &
                      0, 0, &
                      0, 0 )
print *, "CHUNK", status

! Write the mesh
ncells_per_shape(1) = nel

k = 1
do i = 1, nel
    cellIDs(i) = i
    do j = 1, 4
        connectivity(k) = nodelist(j,i)
        k = k + 1
    end do
end do

do i = 1, nnod
    nodeIDs(i) = i
end do

status = TIO_Write_UnstrMesh_Chunk_f( &
    fileID = fileID, &
    meshID = meshID, &
    idx = 1_TIO_SIZEK, &
    xfer = TIO_XFER_NULL_F, &
    graph_datatype = TIO_INTEGER4_F, &
    coord_datatype = TIO_REAL8_F, &
    nodeIDs = nodeIDs, &
    cellIDs = cellIDs, &
    shapes = shapes, &
    ncells_per_shape = ncells_per_shape, &
    connectivity = connectivity, &
    icoords = xv, &
    jcoords = yv &
)
print *, "STATUS", status


! Add some quants
call add_quant("Density", rho)
call add_quant("Pressure", pre)
call add_quant("Energy", en)
call add_quant("X Velocity", uv)
call add_quant("Y Velocity", vv)


! Close the mesh
status = TIO_Close_Mesh_f(fileID, meshID)

! Close the state
status = TIO_Close_State_f(fileID, stateID)

! Close the file
status = TIO_Close_f(fileID)

contains
subroutine add_quant(quantname, data)
character(len=*), intent(in) :: quantname
real(kind=real64), dimension(:), intent(in) :: data

! TIO Stuff
integer(kind=TIO_OBJK) :: quantID
integer(kind=TIO_IPK) :: centring

if (size(data) == nel) then
    print *, "Writing cell-centred quant: ", quantname
    centring = TIO_CENTRE_CELL_F
elseif (size(data) == nnod) then
    print *, "Writing node-centred quant: ", quantname
    centring = TIO_CENTRE_NODE_F
else
    print *, "Unknown quant centring"
    stop 98
end if

status = TIO_Create_Quant_f( &
    fileID, &
    meshID, &
    name = quantname, &
    quantID = quantID, &
    datatype = TIO_REAL8_F, &
    centring = centring, &
    nghosts = 0, &
    ismixed = .false., &
    units = "" &
)

status = TIO_Write_UnstrQuant_Chunk_f( &
    fileID, &
    quantID, &
    idx = 1, &
    xfer = TIO_XFER_NULL_F, &
    datatype = TIO_REAL8_F, &
    qdat = data&
)

status = TIO_Close_Quant_f(fileID, quantID)


end subroutine

end subroutine

end module
