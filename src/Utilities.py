import numpy as np
import pandas as pd
import os
import math
import logging
import sys
import matplotlib.pyplot as plt

try:
  import filterbank
  import sigproc
  import pyfits
  import presto
except ImportError: pass

from Parameters import *


def read_filterbank(filename,DM,bin_start,duration,offset,RFI_reduct=False,mask=False):
  if (not 'filterbank' in sys.modules) or (not 'sigproc' in sys.modules):
    logging.warning("Utilities - Additional modules missing")
    return
  
  #Read data from filterbank file
  head = filterbank.read_header(filename)[0]

  nbits = head['nbits']
  nchans = head['nchans']
  dtype = filterbank.get_dtype(nbits)

  filfile = open(filename, 'rb')
  filfile.seek(0)
  paramname = ""
  while (paramname != 'HEADER_END'): paramname, val = sigproc.read_hdr_val(filfile)
  header_size = filfile.tell()
  filfile.close()

  file_size = os.stat(filename)[6]
  data_size = file_size - header_size
  bytes_per_spectrum= nchans * nbits / 8
  nspec = data_size / bytes_per_spectrum
  
  bin_start -= offset
  if bin_start<0: bin_start = 0
  bin_end = bin_start + DM2delay(DM) + duration + 2*offset
  if bin_end >= nspec: bin_end = nspec
  
  spectrum = np.memmap(filename, dtype=dtype, mode='r', offset=header_size, shape=(nspec, nchans))
  spectrum = spectrum[bin_start:bin_end].copy()
  spectrum = np.fliplr(spectrum)
  
  if mask:
    #Zap the channels from the mask file
    zap_chans, zap_ints, chans_per_int, ptsperint = read_mask(mask,bin_start,bin_end)
    med_value = np.median(spectrum)
    
    spectrum[:,zap_chans] = med_value
    zap_ints = zap_ints[(zap_ints>=bin_start)&(zap_ints<=bin_end)]
    zap_ints -= bin_start
    spectrum[zap_ints,:] = med_value
    
    mask_start = int(np.floor(float(bin_start)/ptsperint)) * ptsperint - bin_start
    for i,chans in enumerate(chans_per_int):
      idx_start = mask_start + i * ptsperint
      idx_end = mask_start + (i + 1) * ptsperint
      spectrum[idx_start:idx_end,chans] = med_value

  if RFI_reduct:
      np.clip(spectrum - np.median(spectrum,axis=0) + 128, 0, 255, out=spectrum)    
  
  return spectrum
 


def read_fits(fits,DM,bin_start,duration,offset,RFI_reduct=False):
  bin_start = np.int(bin_start)
  
  duration = np.int(duration/RES+1)
  
  if isinstance(fits,str):
    try: fits = pyfits.open(fits,memmap=True)
    except NameError: 
      logging.warning("Utilities - Additional modules missing")
      return    
  header = fits['SUBINT'].header
  
  N_channels = header['NCHAN']
  N_spectra = header['NSBLK']
  spectra_per_file = N_spectra * header['NAXIS2']

  bin_start -= offset
  if bin_start<0: bin_start = 0
  subint_start = bin_start/N_spectra

  bin_end = bin_start + DM2delay(DM) + duration + 2*offset
  if bin_end >= spectra_per_file: 
    bin_end = spectra_per_file
    subint_end = bin_end/N_spectra
  else: subint_end = bin_end/N_spectra+1

  subint = fits['SUBINT'].data[subint_start:subint_end]['DATA']
  fits.close()
  subint = subint.reshape((subint_end-subint_start)*N_spectra,N_channels)
    
  #ATTENZIONE! Salvare .mask file nella pipeline!!
  try:
    #Zap the channels from the mask file
    filename = ''
    zap_chans, zap_ints, chans_per_int = read_mask(filename,subint_start,subint_end)
    med_value = np.median(subint)
    subint[:,zap_chans] = med_value
    zap_ints = zap_ints[(zap_ints>=bin_start)&(zap_ints<=bin_end)]
    subint[zap_ints,:] = med_value

    for i,chans in enumerate(chans_per_int):
      chunk = subint[i*N_spectra:(i+1)*N_spectra]
      chunk[:,chans] = med_value
      if RFI_reduct: np.clip(chunk - np.median(chunk,axis=0) + 128, 0, 255, out=chunk)
    
  except IOError:
    if RFI_reduct:
      for i in range(subint_end-subint_start):
        chunk = subint[i*N_spectra:(i+1)*N_spectra]
        np.clip(chunk - np.median(chunk,axis=0) + 128, 0, 255, out=chunk)
    
  bin_start = bin_start%N_spectra
  bin_end = (subint_end-subint_start-1)*N_spectra+bin_end%N_spectra+1
  subint = subint[bin_start%N_spectra:(subint_end-subint_start-1)*N_spectra+bin_end%N_spectra+1]
  
  return subint



def read_mask(filename,bin_start,bin_end):
  with open(filename,'r') as f:
    timesigma, freqsigma, mjd, dtint, lofreq, dfreq = np.fromfile(f,count=6,dtype=np.double)
    numchan, numint, ptsperint = np.fromfile(f,count=3,dtype=np.intc)
    
    num_zap_chans, = np.fromfile(f,count=1,dtype=np.intc)
    zap_chans = np.fromfile(f,count=num_zap_chans,dtype=np.intc)

    num_zap_ints, = np.fromfile(f,count=1,dtype=np.intc)
    zap_ints = np.fromfile(f,count=num_zap_ints,dtype=np.intc)
    
    num_chans_per_int = np.fromfile(f,count=numint,dtype=np.intc)
    num_chans_per_int[num_chans_per_int==numchan] = 0
    num_chans_cum = np.cumsum(num_chans_per_int)
    chans = np.fromfile(f,dtype=np.intc)

  subint_start = int(np.floor(float(bin_start)/ptsperint))
  subint_end = int(np.ceil(float(bin_end)/ptsperint))
  
  chans_per_int = []
  for i in range(subint_start,subint_end):
    if i==0:
      chans_per_int.append(chans[:num_chans_cum[i]])
    else:
      chans_per_int.append(chans[num_chans_cum[i-1]:num_chans_cum[i]])
    
  return zap_chans, zap_ints, chans_per_int, ptsperint
  


def clean(data):

  data_wso = data.astype(np.int)

  masked = np.ma.array(data_wso)
  masked[:,np.arange(0,2592,16)] = np.ma.masked


  #Set the mean of each spectra to zero
  #med_col = np.mean(masked,axis=1)
  med_col = np.min((\
   np.mean(masked[:,:2592/5],axis=1),\
   np.mean(masked[:,-2592/5:],axis=1)),\
   axis=0)

  data_wso = data_wso - med_col[:,np.newaxis]

  #Shift the mean to 128
  data_wso = data_wso + 128

  #Set the right proprieties to the data
  data_wso[:,np.arange(0,2592,16)] = 0
  dtype_min = np.iinfo(data.dtype).min
  dtype_max = np.iinfo(data.dtype).max
  np.clip(data_wso, dtype_min, dtype_max, out=data_wso)
  data_wso = np.around(data_wso)
  data = data_wso.astype(data.dtype)

  return data

 
 

def read_header(filename):
  name, ext = os.path.splitext(filename)
  if ext == '.fits':
    try: fits = pyfits.open(filename,memmap=True)
    except NameError: 
      logging.warning("Utilities - Additional modules missing")
      return None, None
    header = fits['SUBINT'].header + fits['PRIMARY'].header
    fits.close()
    return 'fits', header
  elif ext == '.fil':
    if not 'filterbank' in sys.modules:
      logging.warning("Utilities - Additional modules missing")
      return None, None
    return 'fil', filterbank.read_header(filename)[0]
  else: return None, None



def DM2delay(DM):
  delay = 4149. * (F_MIN**-2 - F_MAX**-2)
  return np.int(delay*DM/RES)


def delay4DM(DM):
  delay = 4149. * (F_MIN**-2 - F_MAX**-2)

  DM_steps = np.array((
          2.525,    5.055,    7.585,   10.115,   12.645,   15.175,
         17.705,   20.235,   22.765,   25.295,   27.825,   30.355,
         32.885,   35.415,   37.945,   40.475,   65.815,   91.115,
        116.415,  141.715,  242.965,  344.165,  445.365,  546.565))
  
  DM_n = (np.digitize([DM,], DM_steps) - 1)[0]
  
  if DM_n<9: DM_n = 0.253 * DM_n
  elif DM_n<15: DM_n = 0.253 * DM_n + 0.253
  elif DM_n<19: DM_n = 0.253 * 15 + 25.3 * delay * ( DM_n - 14 )
  else: DM_n = 0.253 * 15 + 25.3 * delay * 4 + 4 * 25.3 * delay * ( DM_n - 18 )
  
  return delay * DM / 2 + DM_n #s


def time2bin(time,DM):
  time -= delay4DM(DM)
  return time/RES


def span(n,k):
  return 54500.*(10./222./c(n,k)**(1./(k-1)))  #span: number of DMs
  
def c(n,k):
  return math.factorial(n)/math.factorial(k)/math.factorial((n-k))

#Probability that k elements of an ensamble n will group together in a space with n dimensions
#Formula approssimata, ricalcolare quella esatta
def p(n,k,dim = 5450.):
  return c(n,k)/dim**(k-1)



def d(n,k):                           
  return math.factorial(n)/math.factorial(n-k)
#def p(n,k,dim):
  #return 1-float(d(dim,(n-k+2)))*dim**(k-2)/dim**n
#float(Utilities.c(n,k)*dim**(n-k)-n+1)*dim/dim**n




#test of the probability formula
def test_p(n,k,tot,dim = 5450.):
  dim = float(dim)
  fav = 0
  for i in range(tot):
    test = pd.Series(np.random.randint(0,dim,n))
    test = test.groupby(by=test).size().max()
    fav += test>=k 
  return np.float(fav)/tot

def color_range(data):
  #Define the color range
  clean = data[data>0]
  min_element = clean.size/20
  max_element = clean.size*9/10
  vmin = np.partition(clean, min_element, axis=None)[min_element]   #invece di subint[subint>0] possibile subint[:-(num_rows/down_fact)]
  vmax = np.partition(clean, max_element, axis=None)[max_element]
  return vmin,vmax


def rrat_period(times, numperiods=20000):
  #Modified version of PRESTO
  ts = np.asarray(sorted(times))
  ps = (ts[1]-ts[0])/np.arange(1, numperiods+1)
  dts = np.diff(ts)
  xs = dts / ps[:,np.newaxis]
  metric = np.sum(np.fabs((xs - xs.round())), axis=1)
  pnum = metric.argmin()
  numrots = xs.round()[pnum].sum()
  p = (ts[-1] - ts[0]) / numrots
  rotations = dts / p
  diff_max = np.max(np.abs(np.round(rotations) - rotations))
  return p, diff_max
  


def closest_int_above(lenght, n=1, divisor=False):
  '''
  Return the closest int above the initial guess which is divisor of lenght.
  The initial guess is based either on a provided divisor or on the number n that is wanted as output
  '''
  lenght = int(lenght)
  n = int(n)
  if divisor: divisor = int(divisor)
  else: divisor = lenght / n
  while lenght % divisor != 0: divisor += 1 #find the closest integer divisor to average
  return divisor
  
  

def DynamicSpectrum(puls,filename,bary=True,mask=False,ax=False,res=RES,f_min=F_MIN,f_max=F_MAX):
  '''
  puls = {'DM':, 'Sigma':, 'Time':, 'Sample':, 'Downfact':, 'Downsample':, }
  '''
  if not os.path.isfile(filename): return
  
  sample = puls['Sample'] * puls['Downsample']
  
  if bary:
    file_type, header = Utilities.read_header(filename)
    if file_type == 'fil':
      MJD = header['tstart']
      RA = str(header['src_raj'])
      if header['src_raj'] >= 100000:
        RA = '{hh}:{mm}:{ss}'.format(hh=RA[:2], mm=RA[2:4], ss=RA[4:])
      else:
        RA = '{hh}:{mm}:{ss}'.format(hh=RA[:1], mm=RA[1:3], ss=RA[3:])
      DEC = str(header['src_dej'])
      if header['src_dej'] >= 100000:
        DEC = '{hh}:{mm}:{ss}'.format(hh=DEC[:2], mm=DEC[2:4], ss=DEC[4:])
      else:
        DEC = '{hh}:{mm}:{ss}'.format(hh=DEC[:1], mm=DEC[1:3], ss=DEC[3:])
    elif file_type == 'fits':
      MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400.
      RA = header['RA']
      DEC = header['DEC']
    else: raise KeyError('File format not supported')
    
    try: v = presto.get_baryv(RA,DEC,MJD,1800.,obs='LF')
    except NameError: 
      logging.warning("LSPplot - Additional modules missing")
      return
    sample += np.round(sample*v).astype(int)
    
  if isinstance(puls['Downfact'], float):
    duration = np.int(np.round(puls.Duration/res)) * puls['Downsample']
  else: duration = puls['Downsample'] * puls['Downfact']
  spectra_border = 20
  offset = duration*spectra_border
  
  #Load the spectrum
  if filename.endswith('.fits'):
    spectrum = read_fits(filename,puls['DM'],sample,duration,offset,RFI_reduct=True)
  elif filename.endswith('.fil'):
    spectrum = read_filterbank(filename,puls['DM'],sample,duration,offset,RFI_reduct=True,mask=mask)
  else:
    print "File type not recognised. This script only works with fits or fil files."
    return
  
  #De-dispersion
  freq = np.linspace(f_min,f_max,spectrum.shape[1])
  time = (4149 * puls['DM'] * (f_max**-2 - np.power(freq,-2)) / res).round().astype(np.int)
  for i in range(spectrum.shape[1]):
    spectrum[:,i] = np.roll(spectrum[:,i], time[i])
  spectrum = spectrum[:2*offset+duration]
  
  spectrum = np.mean(np.reshape(spectrum,[spectrum.shape[0]/duration,duration,spectrum.shape[1]]),axis=1)
  div = closest_int_above(spectrum.shape[1], n=64)
  spectrum = np.mean(np.reshape(spectrum,[spectrum.shape[0],spectrum.shape[1]/div,div]),axis=2)
  
  fig = plt.figure(figsize=(16,9))
  ax1 = plt.subplot2grid((4,1),(0,0),rowspan=3)
  ax2 = plt.subplot2grid((4,1),(3,0))

  #if not ax: ax = plt.subplot(111)
  extent = [(sample-offset)*res,(sample+duration+offset)*res,f_min,f_max]
  ax1.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent)
  ax1.scatter((sample+duration/2)*res,f_min+1,marker='^',s=1000,c='r')
  ax1.axis(extent)
  ax1.set_ylabel('Frequency (MHz)')
  
  flux =  spectrum.mean(axis=1)
  flux -= np.median(flux)
  flux /= flux.max()
  ax2.plot(flux,'k-')
  ax2.set_ylabel('Flux (rel.)')
  ax2.set_xlabel('Time (s)')
  
  return 






def dedisp(spectrum,DM,duration):
  freq = np.linspace(F_MIN,F_MAX,spectrum.shape[1])
  time = (4149 * DM * (F_MAX**-2 - np.power(freq,-2)) / RES).round().astype(np.int)
  for i in range(spectrum.shape[1]):
    spectrum[:,i] = np.roll(spectrum[:,i], time[i])
  duration = np.int(np.round(duration/RES))
  
  
