#!/usr/bin/env python

'''

Downsampler

Written by Daniele Michilli

'''


import os
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pyfits
import numpy as np
import single_pulse_search as sps

def dynspect_plot(idx,idL,filename,raw_file):
  
  #Constants declaration
  res = 36 #s  #bin resolution
  n_bins = 5  #number of bins to plot
  num_of_spectra = 1000  #number of spectra to actually plot

  #Raw data file
  fits = pyfits.open(raw_file,memmap=True)

  #Parameters for the spectrum
  subint_duration = fits['SUBINT'].header['TBIN']*fits['SUBINT'].header['NSBLK']
  total_duration = res*n_bins #s 
  frame_duration = np.int(total_duration/subint_duration)
 
  #Downsample the data
  N_channels = fits['SUBINT'].header['NCHAN']
  N_spectra = frame_duration*fits['SUBINT'].header['NSBLK']
  down_fact = N_spectra / num_of_spectra
  while N_spectra % down_fact != 0: down_fact -= 1 #find the closest integer divisor to average

  #Prepare the DM lines to plot
  freq = np.arange(151,117,-1,dtype=np.float)
  
  #Frequencies to remove
  freq_del = np.arange(0,N_channels,16)

  for ind0 in idx:
    #Set the start of the spectrum
    t0 = (ind0-1) * res
    subint_index = np.int(t0/subint_duration)

    #Load the data
    subint = fits['SUBINT'].data[subint_index:subint_index+frame_duration]['DATA']

    if subint_index+frame_duration > fits['SUBINT'].header['NAXIS2']:
      num_rows = subint_index+frame_duration - fits['SUBINT'].header['NAXIS2']
      end_rows = np.zeros((num_rows,subint.shape[1]))
      subint = np.append(subint,end_rows,axis=0)

    #Average and clean the spectrum
    subint = subint.reshape(N_spectra/down_fact,down_fact,N_channels).mean(axis=1)
    #subint = subint[:,~np.all(subint == 0, axis=0)]
    subint = np.delete(subint,np.arange(0,2592,16),axis=1)
    

    #Define the color range
    clean = subint[subint>0]
    min_element = clean.size/20
    max_element = clean.size*9/10
    vmin = np.partition(clean, min_element, axis=None)[min_element]   #invece di subint[subint>0] possibile subint[:-(num_rows/down_fact)]
    vmax = np.partition(clean, max_element, axis=None)[max_element]
    clean = 0
    
    #Plot the spectrum
    plt.figure(figsize=(20,10))
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (MHz)')
    plt.imshow(subint.T,cmap=mpl.cm.hot_r,origin="lower",aspect='auto',interpolation='nearest',extent=[t0,t0+total_duration,119,151],vmin=vmin,vmax=vmax)
    plt.colorbar()
    x_lines = np.linspace( t0+res, t0+res*(n_bins-1), n_bins-1 )
    y_lines = np.zeros(x_lines.shape[0])
    plt.plot([x_lines,x_lines],[y_lines+118,y_lines+151],'b--')
    plt.axis([t0,t0+total_duration,119,151])

    #plot the DM lines
    DM=300
    time = 4149 * DM * (np.power(freq,-2) - 151.**-2) + t0 + 1
    plt.plot(time,freq,'b-')
    plt.annotate(str(DM),(time[freq==120],120),color='b',horizontalalignment='left',fontweight='bold')

    DM=500
    time = 4149 * DM * (np.power(freq,-2) - 151.**-2) + t0 + 1
    plt.plot(time,freq,'b-')
    plt.annotate(str(DM),(time[freq==120],120),color='b',horizontalalignment='left',fontweight='bold')
 
    DM=700
    time = 4149 * DM * (np.power(freq,-2) - 151.**-2) + t0 + 1
    plt.plot(time,freq,'b-')
    plt.annotate(str(DM),(time[freq==120],120),color='b',horizontalalignment='left',fontweight='bold')

    plt.savefig('{}_{}_DynSpect_{}.pnx'.format(idL,filename,ind0),format='png',bbox_inches='tight',dpi=150)

  fits.close()
  
  



print 'Starting to downsample'

#Define the parser
parser = argparse.ArgumentParser()
parser.add_argument('idL', nargs=1)
parser.add_argument('filenm', nargs=1)
parser.add_argument('raw_file', nargs=1)
args = parser.parse_args()
filenm = args.filenm[0]
idL = args.idL[0]
raw_file = args.raw_file[0]

#Load the timeseries
timeseries = np.fromfile('{}_{}.dat'.format(idL,filenm),dtype=np.float32)



'''
#Set all the timeseries to start at the same time   ---   CHECK THIS!!
if filenm != 'DM0.00':
  #Check the starting time of the current DM timeseriess
  inf_file = open('{}_{}.inf'.format(idL,filenm))
  for line in inf_file:
    if line.startswith(" Epoch of observation"):
      mjd = line.split('=  ')[-1]
    #if line.startswith(" Width of each time series bin"):
      #res = line.split('=  ')[-1]
  inf_file.close()

  #Check the starting time of DM0 timeseriess
  inf_file = open('{}_DM0.00.inf'.format(idL))
  for line in inf_file:
    if line.startswith(" Epoch of observation"):
      mjd_DM0 = line.split('=  ')[-1]
  inf_file.close()
  
  #Calculate the time delay
  dt = (float(mjd_DM0) - float(mjd)) *86400 #s
  res = 0.000491519982460886*4
  didx = int(dt/res)
  
  #didx = (DM0[DM0!=DM0[-1]].size-DM[DM!=DM[-1]].size*4)
  
  #Apply the time delay
  timeseries = np.roll(timeseries,didx)
  
  #PROBLEMI: 1) Downsampling e' tenuto in considerazione? 2) timeseries inizia un fattore 1.6667 rispetto a quanto aspettato
'''



#Downsample to 36s time resloution
if timeseries.size == 7392000:  down_fact = 73920
elif timeseries.size == 3696000: down_fact = 36960
elif timeseries.size == 1848000: down_fact = 18480
else:
  print "Error: length of the timeseries unknown!"
  exit()

#Downsample to 36s
downsampled = np.mean(timeseries.reshape(100,down_fact),axis=1)  #The first value is the length after the downsample, the second the downsampling factor
downsampled.tofile('{}_{}_down_30s.ds'.format(idL,filenm))


#Plot bins that are some factor above the median and that have higher signal compared to DM0
#if filenm != 'DM0.00':
med = np.median(downsampled)
  #DM0 = np.fromfile('{}_DM0.00_down_30s.ds'.format(idL),dtype=np.float32)
  #idx = np.where( ( downsampled > 1.006 * med )&( downsampled > DM0 )&( downsampled > np.roll(DM0,-1) )&( downsampled > np.roll(DM0,-2) ))[0]  #roll(DM0,-1) e' ok!
  
  #idx = DM-np.maximum(DM0,np.roll(DM0,-1),np.roll(DM0,-2))  #da sviluppare!!
  

idx = np.where( downsampled > 1.002 * med )[0]
idx.tofile('{}_{}_down_30s(ind_cand).dx'.format(idL,filenm))

dynspect_plot(idx,idL,filenm,raw_file)
  
  
print 'Finished to downsample' 













def next2_to_n(x):
    #Return the first value of 2^n >= x.
    i = 1L
    while (i < x): i = i << 1
    return i

timeseries = np.fromfile('test_DM546.37.dat',dtype=np.float32)
downfacts = [2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100]

chunklen = 1789676  #number of points in the timeseries
timeseries = timeseries[:chunklen]
fftlen = int(next2_to_n(chunklen+max(downfacts)))
overlap = (fftlen - chunklen)/2





timeseries = scipy.signal.detrend(timeseries, type='linear')
tmpchunk = timeseries.copy()
tmpchunk.sort()

stds = np.sqrt((tmpchunk[detrendlen/40:-detrendlen/40]**2.0).sum() / (0.95*detrendlen))
stds *= 1.148
timeseries /= stds[:,np.newaxis]





extreme = np.zeros(overlap,np.float32)
timeseries = np.concatenate((extreme, timeseries, extreme))

fftd_chunk = sps.rfft(timeseries, -1)
fftd_kerns = sps.make_fftd_kerns(downfacts, fftlen)




for ii, downfact in enumerate(downfacts):
  goodchunk = sps.fft_convolve(fftd_chunk, fftd_kerns[ii], overlap, -overlap)
  
  hibins, hivals = prune_related1(hibins, hivals, downfact)










