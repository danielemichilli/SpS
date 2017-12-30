from cython.parallel cimport prange

import numpy as np
#from matplotlib.mlab import griddata
from scipy import signal
from scipy.interpolate import griddata

def foo():
    cdef int i, j, n

    x = np.zeros((500, 200), float)
    
    
    def func(x, y):
      return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
    grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
    points = np.random.rand(1000, 2)
    values = func(points[:,0], points[:,1])
    
    
    n = x.shape[0]
    for i in prange(n, nogil=True):
        with gil:
            for j in range(100):
                #x[i,:] = np.cos(x[i,:])
                #x[i,:] = signal.detrend(x[i,:], type='constant')
                x[i,:] = griddata(points, values, (grid_x, grid_y), method='cubic')[30]  #Compare with matplotlib.mlab.griddata and matplotlib.tri.Triangulation
                

    return x
    
   
