import os
import numpy as np
import glob
import multiprocessing as mp
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata


''' NEW

DM steps: 1

In standard pipeline:
- Applicare jumps removal
- Produrre 500 timeseries per ogni beam con -clip argument
- Chiamare beam_matrix
- Rimuovere timeseries
- Salvare matrici numpy
In LSPs:
- Chiamare datacube_creation
- Chiamare FRBs_finder




'''

#datacube(61,500,36000)

Number_of_candidates = 10


def beam_matrix():
  """
  Load the timesieries for a beam and transform the signal in SNR
  """

  #Create the matrix of the beam in DM-Time space and fill it
  beam = np.zeros((500,36000))
  for idx,DM in enumerate(np.arange(500)):
    try: file = glob.glob('*_DM{}.dat'.format(DM))[0]
    except IndexError: continue
    beam[idx] = np.fromfile(file,dtype=np.float32)[:36000]
  
  #Convert signal to SNR
  beam = signal.detrend(beam, type='constant', axis=1)
  stds = beam.copy()
  stds.sort(axis=1)
  stds = stds[:,900:-900]**2
  stds = np.sum(stds, axis=1)
  np.sqrt(stds / 34200., out=stds)
  beam /= stds[:,None]
  
  #Calculate the maximum value of the integrated signal for different downfactor windows in time
  downfacts = [2, 5, 10, 20, 50, 100, 200, 500, 1000]
  for downfact in downfacts:
    for idx,ts in enumerate(beam):
      conv = np.convolve(ts, np.ones(downfact) / np.sqrt(downfact), mode='same')
      ts = np.max((conv,ts),axis=0)
  
  np.save('../FRB_beam_SNR',beam)
  
  return
  
  
  
  
def datacube_creation():
  #Create the matrix of the datacube in beam-DM-Time space and fill it 
  datacube = np.zeros((61,500,36000))  #This will occupy ~9 GB
  for idx,beam in enumerate(np.arange(61)):
    try: beam[idx] = np.load('FRB_beam_SNR.npz')
    except IOError: continue
  
  #SNR along time to SNR along beam
  datacube = signal.detrend(datacube, type='constant', axis=0)
  stds = datacube.copy()
  stds.sort(axis=0)
  stds = stds[1:-1]**2
  stds = np.sum(stds, axis=0)
  np.sqrt(stds / 57.95, out=stds)
  datacube /= stds[:,None]
  
  datacube = datacube_parallel(datacube)
  
  
  
  
def datacube_parallel(datacube):
  CPUs = mp.cpu_count()
  idx_range = int(np.ceil(datacube.shape[-1]/float(CPUs)))  
  data_list = [datacube[...,i*idx_range:(i+1)*idx_range] for i in range(CPUs)]
  
  pool = mp.Pool()
  results = pool.map(space_fft, data_list)
  pool.close()
  pool.join()
  
  return np.dstack(results)
  
  

  
  
def space_fft(data):
  ra = np.array([ 8,  8, 10, 10,  8,  6,  6,  8, 10, 12, 12, 12, 10,  8,  6,  4,  4,
        4,  6,  8, 10, 12, 14, 14, 14, 14, 12, 10,  8,  6,  4,  2,  2,  2,
        2,  4,  6,  8, 10, 12, 14, 16, 16, 16, 16, 16, 14, 12, 10,  8,  6,
        4,  2,  0,  0,  0,  0,  0,  2,  4,  6])
  
  dec = np.array([ 8, 10,  9,  7,  6,  7,  9, 12, 11, 10,  8,  6,  5,  4,  5,  6,  8,
       10, 11, 14, 13, 12, 11,  9,  7,  5,  4,  3,  2,  3,  4,  5,  7,  9,
       11, 12, 13, 16, 15, 14, 13, 12, 10,  8,  6,  4,  3,  2,  1,  0,  1,
        2,  3,  4,  6,  8, 10, 12, 13, 14, 15])
  
  grid = np.arange(17)

  for i in range(data.shape[1]):
    for j in range(data.shape[2]):
      ts = data[:,i,j]

      beams_map = griddata(ra, dec, ts, grid, grid, interp='linear')
      beams_map[beams_map.mask==True] = 0
      conv = signal.convolve2d(beams_map,np.ones((fill,fill))/fill,mode='same')  #Check the statistics!


  
  
