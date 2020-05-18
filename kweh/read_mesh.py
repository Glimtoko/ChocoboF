import numpy as np

def get_mesh_size(filename):
    """
    Determines the number of nodes and cells in the mesh
    """
    nnodes = 0
    ncells = 0
    
    innodelist = False
    incelllist = False
    with open(filename, "r") as mfile:
        for line in mfile:
            line = line.strip()
            if "*" in line:
                innodelist = False
                incelllist = False
                
            if innodelist:
                nnodes += 1
                
            if incelllist:
                ncells += 1
                
            if "*NODE" in line:
                innodelist = True
                
            if "PhysicalSurface" in line:
                incelllist = True
                
    return nnodes, ncells
        

def read_mesh(filename, boundaries, material_list, mesh):
    """
    Reads the mesh
    """
    innodelist = False
    incelllist = False
    inbcs = False
    
    reg = None
    bc = None
    ignore_next = False
    index_cell = -1
    index_node = -1
    
    # Zero out BCs
    for i in range(mesh.nnodes):
        mesh.type[i] = 0.0 
    
    with open(filename, "r") as mfile:
        for line in mfile:
            line = line.strip()
            if "*" in line:
                innodelist = False
                incelllist = False
                
            # Nodal coordinates
            if innodelist:
                index_node += 1
                data = line.split(",")
                
                mesh.x[index_node] = float(data[1])
                mesh.y[index_node] = float(data[2])
                
            # Element-Node connectivity
            if incelllist:
                index_cell += 1
                
                data = line.split(",")
                
                n = [int(d) for d in data[2:6]]
                for i in range(4):
                    mesh.nodelist[index_cell,i] = n[i]
                mesh.region[index_cell] = reg
                
            # Boundary conditions
            if inbcs:
                if ignore_next:
                    ignore_next = False
                else:
                    if "#" in line:
                        bc = line.split("$# ")[1]
                        ignore_next = True
                    elif "*" not in line:
                        nodes = [int(n)-1 for n in line.split(",")]
                        for node in nodes:
                            mesh.type[node] += boundaries[bc]
                            
            # Detect start of nodelist
            if "*NODE" in line:
                print("  Node list found")
                innodelist = True
                
            # Detect start of region
            if "PhysicalSurface" in line:
                reg = int(line.split("PhysicalSurface")[1])/11
                print("  Found region {}".format(reg))
                incelllist = True
                
            # Detect node lists
            if "*SET_NODE_LIST" in line:
                inbcs = True
                        
    # Determine region/cell iterator forms
    nreg = max(mesh.region)
    print("  Mesh contains {} regions".format(nreg))
    
    canuserange = True  # Assume we can use simple range interators
    for reg in range(1,nreg+1):
        print("  Checking region {}".format(reg))
        all_elements = np.where(mesh.region == reg)
        first = all_elements[0][0]
        last = all_elements[0][-1]
              
        for i in range(first, last+1):
            if mesh.region[i] != reg:
                print("  Error: Non-consecutive region cell iteration not coded")
                canuserange = False
                
        if not canuserange:
            return False
        else:
            print("    Region is valid")
            mesh.regioncelliterator[reg] = (first+1, last+1)
            
    for i, reg in enumerate(mesh.region):
        mesh.material[i] = material_list[reg]
        
