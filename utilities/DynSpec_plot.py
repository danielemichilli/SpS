#!/usr/bin/env python

t0 = 20  #s
file_name = 'L337502_SAP0_BEAM0.fits'

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pyfits

#Open the raw file
fits = pyfits.open(file_name,memmap=True)

#Set the parameters for the spectrum
subint_duration = fits['SUBINT'].header['TBIN']*fits['SUBINT'].header['NSBLK']
subint_index = np.int(t0/subint_duration)
total_duration = 432  #s
frame_duration = np.int(36/subint_duration)*12

#Load the data
subint = fits['SUBINT'].data[subint_index:subint_index+frame_duration]['DATA']

#Downsample the data
N_channels = fits['SUBINT'].header['NCHAN']
N_spectra = frame_duration*fits['SUBINT'].header['NSBLK']
subint = subint.reshape(N_spectra/384,384,N_channels).mean(axis=1)
subint = subint[:,~np.all(subint == 0, axis=0)]

#Plot the spectrum
plt.figure(figsize=(20,10))
plt.xlabel('Time (s)')
plt.ylabel('Frequency (MHz)')
plt.imshow(subint.T,cmap=mpl.cm.hot_r,origin="lower",aspect='auto',interpolation='nearest',extent=[t0,t0+total_duration,119,151])
plt.colorbar()
x_lines = np.linspace(t0+36,t0+36*11,10)
y_lines = np.zeros(x_lines.shape[0])
plt.plot([x_lines,x_lines],[y_lines+118,y_lines+160],'b--')
plt.axis([t0,t0+total_duration,119,151])

for DM in [366,458,546]:
  freq = np.arange(160,117,-1,dtype=np.float)
  time = 4149 * DM * (np.power(freq,-2) - 160.**-2) + 5
  plt.plot(time,freq,'w-')
  plt.annotate(str(DM),(time[-1],122),color='white')

plt.show()
#plt.savefig('test.png',format='png',bbox_inches='tight')


