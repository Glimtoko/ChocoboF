# ChocoboF
Chocobo CFD Code - Fortran Version

## Dependencies
The only required dependency for ChocoboF is a Fortran compiler. Currently, ChocoboF has only been tested under GFortran (version 9.2.1), but the intention exists to test it using Flang at some point.

ChocoboF also supports the following output libraries, which will be included in the executable if available on the build system:

* TyphonIO (https://github.com/UK-MAC/typhonio)
* SILO (https://wci.llnl.gov/simulation/computer-codes/silo)
* HDF5 (https://www.hdfgroup.org/solutions/hdf5/) (required if using either TyphonIO or SILO)

## Building ChocoboF
ChocoboF provides a makefile targetting the GFortran compiler. This makefile controls whether or not TyphonIO and/or SILO support is included in the executable. In order to use these optional dependencies, the `TIO_DIR`, `SILO_DIR` and `HDF5_DIR` variables should be set to the paths to the installation locations of the libraries. The standard way of using this makefile is:

    make; make install

By default, the installation path is `/prod/ChocoboF/bin`, with a simlink placed in `/prod/bin`, although these can be changed in the makefile.

In addition, the kweh folder contains a post-processor used to convert meshes into a format which can be read by ChocoboF. This too has a makefile which controls building and installation, which should be used in the same way as the main ChocoboF makefile.

## Running ChocoboF
**Note that, due to limitations of Fortran, the way ChocoboF is run is quite different to other flavours of Chocobo.**

To run ChocoboF, assuming its installation location is in your `$PATH`, simply type:

    chocobof

This will assume that the following input files are located in the current directory:

| File        | Purpose                                   |
|-------------|-------------------------------------------|
| param.dat   | Main input file - uses namelists          |
| mesh.chc    | Mesh file (in ChocoboF format, see below) |
| eos.dat     | EoS data file                             |

Unlike other versions of Chocobo, these file names cannot be changed. However, if placed in a directory, then this location can be passed to Chocobo via a single commandline argument. For example, say you had the following input files in a subfolder as follows:

```
.
+-- sod
    +-- param.dat
    +-- mesh.chc
    +-- eos.dat
```
This could be run in ChocoboF via:

    chocobof sod

## Control (input) file
ChocoboF's input takes the form of a *namelist*. The following details all the current input options:

    integer(kind=int32), public, parameter :: maxreg=10

    real(kind=real64), public:: gamma

    ! artificial viscosity
    real(kind=real64), public:: cq
    real(kind=real64), public:: cl

    !max time step
    real(kind=real64), public:: maxallstep

    ! initial time step
    real(kind=real64), public:: dtinit

    ! timestep option
    integer, public :: dtoption

    ! Problem time extent
    real(kind=real64), public :: t0  !initial time
    real(kind=real64), public :: tf  !final time

    ! time step growth factor
    real(kind=real64), public :: growth

    ! if 0 cartesian 1 axisymmetric
    integer(kind=int32), public :: zaxis

    ! 0 not vol 1 indiv volume change
    integer(kind=int32), public :: zintdivvol

    ! type of artificial viscosity
    integer(kind=int32), public :: avtype

    ! if 0 no hourglassing if 1 hourglassing
    integer(kind=int32), public :: zantihg

    ! Anti-hourglass filter type
    integer(kind=int32), public :: hgregtyp(1:maxreg)

    ! Kappa for A-H filter
    real(kind=real64), public :: kappareg(1:maxreg)

    ! SILO/TIO output controls
    real(kind=real64) :: dtsilo
    integer(kind=int32) :: h5type
    logical :: tioonefile

    ! Maximum number of timesteps - used for debugging
    integer(kind=int32) :: stepcnt = 0

## Mesh file
Meshes should be generated using a mesh generator which supports the following:

* Generation of **quad** meshes
* Output in LSDYNA Key format

The GMSH generator (https://gmsh.info/) has been tested and demonstrated to produce meshes of the form required. Note that by default it generates triangular meshes, and some care is needed to generate a good quad mesh. The following GMSH script will generate a mesh suitable for running the Sod shock tube problem:

    // Mesh physical size
    lc = 0.02;
    l = 1.0;
    h = 0.1;

    // Mesh resolution
    xc = 20;
    yc = 10;

    // Create mesh geometry from Points and Lines

    // Create geometry points
    Point(1) = {0, 0, 0, lc};
    Point(2) = {l/2, 0,  0, lc};
    Point(3) = {l/2, h, 0, lc};
    Point(4) = {0, h, 0, lc};
    Point(5) = {l, 0, 0, lc};
    Point(6) = {l, h, 0, lc};

    // Create geometry lines
    Line(1) = {1, 2};
    Line(2) = {2, 3};
    Line(3) = {3, 4};
    Line(4) = {4, 1};

    Line(5) = {2, 5};
    Line(6) = {5, 6};
    Line(7) = {6, 3};

    // Create surfaces. These produce regions in the mesh. Defined via loops of lines
    Curve Loop(1) = {4, 1, 2, 3};
    Curve Loop(2) = {-2, 5, 6, 7};
    Plane Surface(1) = {1};
    Plane Surface(2) = {2};

    // xc and yc refer to number of *nodes*, so add one to them
    xc++;
    yc++;

    // Create regions, each spanning half the X space
    Transfinite Surface {1} = {1, 2, 3, 4};
    Transfinite Line {1, 3} = xc/2.0;
    Transfinite Line {4, 2} = yc;

    Transfinite Surface {2} = {2, 5, 6, 3};
    Transfinite Line {5, 7} = xc/2.0;
    Transfinite Line {6} = yc;

    // Boundaries
    Physical Curve("YLOW") = {1, 5};
    Physical Curve("YHIGH") = {3, 7};
    Physical Curve("XLOW") = {4};
    Physical Curve("XHIGH") = {6};

    // Set meshing options
    Mesh.Algorithm = 8;
    Mesh.RecombinationAlgorithm = 3;
    Mesh.RecombineAll = 1;

    Mesh.SaveAll = 84;
    Mesh.SaveElementTagType = 1;
    Mesh.SaveTopology = 0;
    Mesh.SaveParametric = 0;
    Mesh.SaveGroupsOfElements = 1;
    Mesh.SaveGroupsOfNodes = 3;

    // Generate mesh
    Mesh 2;

    // Save mesh file
    Save "sod.key";

    Exit;

To correctly set the boundary conditions, the edges of the mesh **must** be included on `Physical Curve`s with the names *XLOW*, *XHIGH*, *YLOW* and *YHIGH*, as shown in the example.

**ChocoboF is, however, unable to directly read a .key file. Instead, the Kweh processor should be used.** The basic usage of Kweh is:

    kweh mesh.key

This will produce a file named mesh.chc. Alternatively, a different output file name may be provided (but remember ChocoboF will be unable to read it unless it is named mesh.chc):

    kweh mesh.key output.chc

## EoS file
The EoS file should consist of one line of data per material, with each line being a comma-separated list of the following parameters:

* Material Number, Initial Density, Initial Pressure, Gamma (for perfect gas EoS)

Note that, unlike ChocoboPy, comments are **not** supported. For example:

    1, 1, 1, 1.4
    2, 0.125, 0.1, 1.4

