import numpy as np



signal = #loaded from timeseries



def convolve(arr,max_length=1):
  arr = scipy.signal.detrend(arr, type='linear')
  chunklen = arr.size
  tmpchunk = arr.copy()
  tmpchunk.sort()
  stds = np.sqrt((tmpchunk[chunklen/40:-chunklen/40]**2.0).sum() / (0.95*chunklen))
  stds *= 1.148
  arr /= stds    

  if not isinstance(max_length,list):
    max_length = range(max_length)
  
  signal = np.zeros((chunklen,len(max_length)))
  
  for i,length in enumerate(max_length):
    signal[i] = np.convolve(arr, np.ones(length), mode='same') / np.sqrt(length)
  
  return np.sum(signal,axis=0)











#MAP MAKER

import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
from scipy import signal

fill = 10

ra = np.array([399, 399, 498, 497, 399, 302, 301, 399, 499, 598, 595, 594, 496,
       399, 303, 206, 203, 201, 300, 399, 500, 600, 698, 695, 692, 690,
       592, 495, 399, 304, 207, 110, 107, 103, 101, 199, 299, 399, 502,
       602, 701, 800, 796, 792, 788, 785, 687, 590, 494, 399, 305, 209,
       112,  14,  11,   7,   3,   0,  98, 198, 298])
 
dec = np.array([400, 500, 450, 350, 300, 350, 450, 600, 550, 500, 400, 300, 250,
       200, 250, 300, 400, 500, 550, 700, 650, 600, 550, 450, 350, 250,
       200, 150, 100, 150, 200, 250, 350, 450, 550, 600, 650, 800, 750,
       700, 650, 600, 500, 400, 300, 200, 150, 100,  50,   0,  50, 100,
       150, 200, 300, 400, 500, 600, 650, 700, 750])

val = np.random.randint(0,100,ra.size)
val[26-13] = 250
val[42-13] = 200
val[41-13] = 190

grid = np.linspace(0,800,8*fill)

beams_map = griddata(ra, dec, val, grid, grid, interp='nn')
beams_map[beams_map.mask==True] = 0

conv = signal.convolve2d(beams_map,np.ones((8*fill,8*fill))/fill,mode='same')

for i,n in enumerate(val):                                                                                               
    plt.annotate(str(n),xy=(ra[i],dec[i]),horizontalalignment='center',verticalalignment='center',color='k',size=18,weight='bold')
    
plt.contourf(grid, grid, conv, 10, cmap=plt.cm.hot_r)

plt.contourf(np.ones((8*fill,8*fill)))

plt.show()







#

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
from scipy import signal

lim1 = 10
lim2 = 

def time_fft(timeseries):
  downfacts = [2, 5, 10, 20, 50, 100, 200, 500, 1000]
  chunklen = timeseries.shape[1]
  for idx in range(timeseries.shape[0]):
    ts = timeseries[idx]
    ts = signal.detrend(ts, type='constant')
    tmpchunk = ts.copy()
    tmpchunk.sort()
    stds = np.sqrt((tmpchunk[chunklen/40:-chunklen/40]**2.0).sum() / (0.95*chunklen))
    stds *= 1.148
    ts /= stds
    
    down_series = np.zeros(ts.size)
    
    for downfact in downfacts:
      goodchunk = np.convolve(ts, np.ones(downfact) / np.sqrt(downfact), mode='same')
      down_series = np.max((goodchunk,down_series),axis=0)
      #store on a second array of ints the downfactor value for the highest SNR 
    
    ts = down_series
    
  return np.unravel_index(np.argpartition(timeseries, -lim1, axis=None)[-lim1:],timeseries.shape)[0]
  #return np.where(timeseries.max(axis=0)>lim)


def space_fft(timeseries,good_times):
  fill = 10

  ra = np.array([399, 399, 498, 497, 399, 302, 301, 399, 499, 598, 595, 594, 496,
        399, 303, 206, 203, 201, 300, 399, 500, 600, 698, 695, 692, 690,
        592, 495, 399, 304, 207, 110, 107, 103, 101, 199, 299, 399, 502,
        602, 701, 800, 796, 792, 788, 785, 687, 590, 494, 399, 305, 209,
        112,  14,  11,   7,   3,   0,  98, 198, 298])
  
  dec = np.array([400, 500, 450, 350, 300, 350, 450, 600, 550, 500, 400, 300, 250,
        200, 250, 300, 400, 500, 550, 700, 650, 600, 550, 450, 350, 250,
        200, 150, 100, 150, 200, 250, 350, 450, 550, 600, 650, 800, 750,
        700, 650, 600, 500, 400, 300, 200, 150, 100,  50,   0,  50, 100,
        150, 200, 300, 400, 500, 600, 650, 700, 750])

  grid = np.linspace(0,800,8*fill)
  
  good_spaces = []
  chunklen = timeseries.shape[1]
  for idx in good_times:
    ts = timeseries[:,idx]
    ts = signal.detrend(ts, type='constant')
    tmpchunk = ts.copy()
    tmpchunk.sort()
    stds = np.sqrt((tmpchunk[chunklen/40:-chunklen/40]**2.0).sum() / (0.95*chunklen))
    stds *= 1.148
    ts /= stds

    beams_map = griddata(ra, dec, ts, grid, grid, interp='nn')
    conv = signal.convolve2d(beams_map,np.ones((fill,fill)),mode='same') / fill
    
    if conv.max()>lim2:
      good_spaces.append(idx)
      
  return np.array(good_spaces)
  
    
    


obs = np.random.randint(0,255,(61,1000,500)).astype(np.float32)
#data = np.random.randint(0,255,(61,1000)).astype(np.float32)

for idx in data.shape[2]:
  data = obs[:,:,idx]

  good_times = time_fft(data.copy())

  good_spaces = space_fft(data.copy(),good_times)


#selezionare DM con segnale piu' grande tra tutti i good_spaces
  

#timeseries = np.ma.masked_less(timeseries,lim)
#timeseries.mask -= 1
#timeseries = np.ma.mask_cols(timeseries)
#timeseries.mask -= 1











