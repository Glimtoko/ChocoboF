#!/prod/anaconda3/bin/python
import sys
import os

import data_structures as ds
import read_mesh as rm
import output
 
# Assume a standard set of boundary conditions
boundary_values = {
    "XLOW": -1,
    "XHIGH": -1,
    "YLOW": -2,
    "YHIGH": -2,
}

# Assume a 1:1 mapping of region to material number. We could change this later
material_list = {}
for i in range(1, 21):
    material_list[i] = i

try:
    mfile = sys.argv[1]
except IndexError:
    print("Usage: kweh <mesh filename> [output filename]")
    sys.exit(-1)
    
try:
    ofile = sys.argv[2]
except IndexError:
    ofile = "{}.chc".format(mfile)
    
print("KWEH - Mesh Convertor for ChocoboF\n")
print("Input key file: {}".format(mfile))
print("Output chc file: {}\n".format(ofile))

if not os.access(mfile, os.R_OK):
    print("ERROR: Input file not found")
    sys.exit(-1)

print("Input file exists!")

nnodes, ncells = rm.get_mesh_size(mfile)
print("\nMesh dimensions:")
print("  Number of nodes: {}".format(nnodes))
print("  Number of cells: {}".format(ncells))

mesh = ds.Mesh(nnodes, ncells)

print("\nReading mesh")
rm.read_mesh(mfile, boundary_values, material_list, mesh)

print("\nCreating output file")
with open(ofile, "w") as ofile:
    output.output(ofile, mesh)