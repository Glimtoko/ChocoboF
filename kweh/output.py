import numpy as np

def output(ofile, mesh):
    ofile.write("{}\n".format(mesh.nnodes))
    ofile.write("{}\n".format(mesh.ncells))
    ofile.write("{}\n".format(max(mesh.region)))
    
    # Node-centred
    np.savetxt(ofile, mesh.x, fmt="%.9e", newline="   ")
    ofile.write("\n")
    np.savetxt(ofile, mesh.y, fmt="%.9e", newline="   ")
    ofile.write("\n")
    np.savetxt(ofile, mesh.type, fmt="%3d", newline="   ")
    ofile.write("\n")
    
    # Cell-centred
    np.savetxt(ofile, mesh.region, fmt="%3d", newline="   ")
    ofile.write("\n")
    np.savetxt(ofile, mesh.material, fmt="%3d", newline="   ")
    ofile.write("\n")
    
    # Node list - convert to a Fortran-ordered 1D array
    nl = np.zeros(4*mesh.ncells)
    
    k = 0
    for j in range(mesh.ncells):
        for i in range(4):
            nl[k] = mesh.nodelist[j,i]
            
            k += 1
            
    np.savetxt(ofile, nl, fmt="%3d", newline="   ")
    ofile.write("\n")
    
    for i in range(max(mesh.region)):
        reg = i + 1
        fl = mesh.regioncelliterator[reg]
        ofile.write("{}   {}   {}\n".format(reg, fl[0], fl[1]))
        