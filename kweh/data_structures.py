import numpy as np
      
        
class Mesh():
    def __init__(self, nn, nc):
        self.ncells = nc
        self.nnodes = nn
 
        # Node centred       
        self.x = np.ndarray(nn, np.double)
        self.y = np.ndarray(nn, np.double)       
        self.type =  np.ndarray(nn, int)
               
        # Cell-centred
        self.region = np.ndarray(nc, int)
        self.material = np.ndarray(nc, int)
 
        # Special
        self.nodelist = np.ndarray((nc,4), int)       
        self.regioncelliterator = {}
