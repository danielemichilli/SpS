import argparse
import numpy as np
import single_pulse_search as sps
import scipy
import tarfile



'''OLD

DM steps: 0.45835940634304545
  DMs = dt / (4149 * (F_MIN**-2 - F_MAX**-2))
  dt = 0.05

In standard pipeline:
- Applicare jumps removal
- Produrre n timeseries in ogni beam #Trovare DM step#
- In LTA_fil.sh chiamare FRB.py per ogni DM step
  1.Rimuovere bad bins
    a.Basandosi su statistiche fits file
    b.Basandosi su statistiche timeseries
    c.Basandosi su DM0 timeseries
      + clip
      + rimuovere span in timeseries #Quando DM0>DMx rimuovere time span precedente#
      + rimuovere jumps estremi rispetto a bin precedente in DM0
  2.Downsample
  3.Rimuovere bad chunks
  4.Cercare candidati e appendere in DB (o txt)
  5.Salvare i DB di ogni beam in tar
in SP pipeline:
- Unire tutti i DB
- Applicare RFI removal
- Plottare dynamic spectrum e segnale in ogni beam dei candidati interessanti

'''




''' NEW

DM steps: 1

In standard pipeline:
- Applicare jumps removal
- Produrre 500 timeseries per ogni beam con -clip argument
- Caricare come tabella in memoria e applicare convolution per ogni DM
- Produrre array con segnale maggiore ad ogni tempo per ogni DM piu' DM e downfact corrispondenti



'''



def downsample(idL,filenm,raw_file,timeseries=None):
  print 'Starting to downsample'

  #Load the timeseries
  if timeseries is None:
    timeseries = np.fromfile('{}_{}.dat'.format(idL,filenm),dtype=np.float32)

  #Downsample to 0.5s time resloution
  if timeseries.size == 7392000:  down_fact = 100
  elif timeseries.size == 3696000: down_fact = 50
  elif timeseries.size == 1848000: down_fact = 25
  else:
    print "Error: length of the timeseries unknown!"
    exit()

  #Downsampled
  downsampled = np.mean(timeseries.reshape(73920,down_fact),axis=1)  #The first value is the length after the downsample, the second the downsampling factor
  
  #downsampled.tofile('{}_{}_down.ds'.format(idL,filenm))

  print 'Downsampled complete'
  
  return



def search_ts():

  #Plot bins that are some factor above the median and that have higher signal compared to DM0
  #if filenm != 'DM0.00':
  med = np.median(downsampled)
  #DM0 = np.fromfile('{}_DM0.00_down_30s.ds'.format(idL),dtype=np.float32)
  #idx = np.where( ( downsampled > 1.006 * med )&( downsampled > DM0 )&( downsampled > np.roll(DM0,-1) )&( downsampled > np.roll(DM0,-2) ))[0]  #roll(DM0,-1) e' ok!
  
  #idx = DM-np.maximum(DM0,np.roll(DM0,-1),np.roll(DM0,-2))  #da sviluppare!!
    

  #idx = np.where( downsampled > 1.002 * med )[0]
  #idx.tofile('{}_{}_down_30s(ind_cand).dx'.format(idL,filenm))


  threshold = 1.5


  timeseries = np.fromfile('test_DM546.37.dat',dtype=np.float32)


  downfacts = [2, 5, 10, 20, 50, 100, 200, 500, 1000]  #x0.05s

  chunklen = timeseries.size #73920  #number of points in the timeseries

  timeseries = scipy.signal.detrend(timeseries, type='linear')
  tmpchunk = timeseries.copy()
  tmpchunk.sort()
  stds = np.sqrt((tmpchunk[chunklen/40:-chunklen/40]**2.0).sum() / (0.95*chunklen))
  stds *= 1.148
  timeseries /= stds

  #Modifica a sps
  #timeseries = scipy.signal.detrend(timeseries, type='linear')
  for downfact in downfacts:
    goodchunk = np.convolve(timeseries, np.ones(downfact), mode='same') / np.sqrt(downfact)
    
    hibins = np.nonzero(goodchunk > threshold).tolist()
    hivals = goodchunk.tolist()
    
    hibins, hivals = sps.prune_related1(hibins, hivals, downfact)    
    

    #for bin, val in zip(hibins, hivals):
      #time = bin * dt
      #bisect.insort(dm_candlist, candidate(info.DM, val, time, bin, downfact))

















if __name__ == '__main__':
  #Define the parser
  parser = argparse.ArgumentParser()
  parser.add_argument('idL', nargs=1)
  parser.add_argument('filenm', nargs=1)
  parser.add_argument('raw_file', nargs=1)
  args = parser.parse_args()
  
  downsample(args.idL[0],args.filenm[0],args.raw_file[0])


  

import os
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pyfits
import numpy as np




  plot()
  
  
  
def plot():
  
  
  
  

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
  return


  


import pandas as pd
import os
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

os.system("module load python/2.7.9")

def ra_dec(folder,idL,sap):

  ra = np.zeros(61)
  dec = np.zeros(61)

  k = 0
  for beam in range(13,74):
    try:
      name = '{}_SAP{}_BEAM{}'.format(idL,sap,beam)
      path = 'SAP{}/{}/BEAM{}_sift/sp/'.format(sap,name,beam)
      events_path = '{}{}/{}{}_singlepulse.tgz'.format(folder,idL,path,name)
      tar_file = tarfile.open(events_path)
      inf_file = tar_file.extractfile(name+'.inf')
      
      for line in inf_file:
        if 'Right Ascension' in line:
          coo = line.split('=  ')[-1]
          h = coo.split(':')[0]
          m = coo.split(':')[1]
          s = coo.split(':')[2].split('.')[0]
          ra[k] = float(h)+float(m)/60+float(s)/3600
          
          
        if 'Declination' in line:
          coo = line.split('=  ')[-1]
          d = coo.split(':')[0]
          m = coo.split(':')[1]
          s = coo.split(':')[2].split('.')[0]
          dec[k] = float(d)+float(m)/60+float(s)/3600

      
    except IOError:
      ra[k] = np.nan
      dec[k] = np.nan
    
    k += 1

  return ra,dec




def plot_beams(idL,beams,name,ra,dec,t_min,t_max):
  plt.clf()
  plt.scatter(ra,dec,s=800,edgecolor='none',c=beams,cmap=mpl.cm.hot_r)
  bar = plt.colorbar()
  bar.set_label('Flux (Arb. Units)')
  plt.xlabel('RA (h)')
  plt.ylabel('DEC (deg)')
  [plt.annotate(str(i+13),(ra[i],dec[i]),horizontalalignment='center',verticalalignment='center') for i in range(0,61)]

  plt.annotate('DM = 546.37 pc/cm3\nTime = {} - {} s'.format(t_min,t_max), xy=(np.nanmin(ra)-0.2,np.nanmax(dec)+0.1))

  plt.xlim(np.nanmin(ra)-0.3,np.nanmax(ra)+0.3)
  plt.ylim(np.nanmin(dec)-0.2,np.nanmax(dec)+0.4)

  plt.savefig(name,format='png',bbox_inches='tight',dpi=200)


  folder = '/projects/lotaas/test/Daniele/out/'

  for sap in range(0,3):

    ra,dec = ra_dec(folder,idL,sap)

    if np.all(np.isnan(ra)) | np.all(np.isnan(dec)):

      ra = np.array([ 
          0.        ,  0.        ,  0.47083334,  0.46666667,  0.        ,
        -0.46666666, -0.47083333,  0.        ,  0.47916667,  0.95      ,
          0.9375    ,  0.92916667,  0.4625    ,  0.        , -0.45833333,
        -0.925     , -0.9375    , -0.94583333, -0.475     ,  0.        ,
          0.48333334,  0.95833334,  1.42916667,  1.4125    ,  1.4       ,
          1.3875    ,  0.92083334,  0.45833334,  0.        , -0.45416666,
        -0.91666666, -1.38333333, -1.39583333, -1.4125    , -1.425     ,
        -0.95416666, -0.47916666,  0.        ,  0.4875    ,  0.96666667,
          1.44166667,  1.9125    ,  1.89583334,  1.875     ,  1.85833334,
          1.84166667,  1.375     ,  0.9125    ,  0.45416667,  0.        ,
        -0.45      , -0.90833333, -1.37083333, -1.8375    , -1.85416666,
        -1.87083333, -1.89166666, -1.90833333, -1.4375    , -0.9625    ,
        -0.48333333])
  
      dec = np.array([ 
          0.        ,  0.28305556,  0.14166667, -0.14166666, -0.28305555,
        -0.14166666,  0.14166667,  0.56638889,  0.42472222,  0.28305556,
          0.        , -0.28305555, -0.42472222, -0.56638889, -0.42472222,
        -0.28305555,  0.        ,  0.28305556,  0.42472222,  0.84944445,
          0.70777778,  0.56638889,  0.42472222,  0.14166667, -0.14166666,
        -0.42472222, -0.56638889, -0.70777778, -0.84944444, -0.70777778,
        -0.56638889, -0.42472222, -0.14166666,  0.14166667,  0.42472222,
          0.56638889,  0.70777778,  1.1325    ,  0.99111111,  0.84944445,
          0.70777778,  0.56638889,  0.28305556,  0.        , -0.28305555,
        -0.56638889, -0.70777778, -0.84944444, -0.99111111, -1.13277778,
        -0.99111111, -0.84944444, -0.70777778, -0.56638889, -0.28305555,
          0.        ,  0.28305556,  0.56638889,  0.70777778,  0.84944445,
          0.99111111])

    signal = np.full([61,100],np.nan,dtype=np.float32)
    
    for beam in range(13,74):
      
      name = '{}_SAP{}_BEAM{}'.format(idL,sap,beam)
      path = '{}{}/SAP{}/BEAM{}_sift/sp/'.format(folder,idL,sap,beam)
      events_path = '{}{}_singlepulse_down.tgz'.format(path,name)
      num = beam-13
      
      try:
        #Open the file
        tar_file = tarfile.open(events_path)
        tar_file.extractall(path)
        tar_file.close()
        
        try: 
          signal[num] = np.fromfile(path+name+'_DM546.37_down_30s.ds',dtype=np.float32)
          [ os.remove(path+f) for f in os.listdir(path) if f.endswith(".ds") ]
        except IOError: pass
      
      except IOError: pass
    
    #median value for each bin over all the beams
    med = np.nanmedian(signal,axis=0)  
    
    #extract interesting signals
    for beam in range(13,74):
      name = '{}_SAP{}_BEAM{}'.format(idL,sap,beam)
      path = '{}{}/SAP{}/BEAM{}_sift/sp/'.format(folder,idL,sap,beam)
      num = beam-13

      try: 
        idx = np.fromfile(path+name+'_DM546.37_down_30s(ind_cand).dx',dtype=np.int64)
        os.remove(path+name+'_DM546.37_down_30s(ind_cand).dx')
      except IOError: continue

      for i in idx[ signal[num,idx] > med[num] ]:
        new_path = '{}{}/sp/ALERTS/'.format(folder,idL)
        if not os.path.exists(new_path): os.makedirs(new_path)
        plot = name+'_DM546.37_DynSpect_{}'.format(i)
        os.rename('{}{}.pnx'.format(path,plot),'{}{}.png'.format(new_path,plot))
        
        t_min = i * 36
        beams_name = new_path + name + '_DM546.37_{}s.png'.format(t_min)
        plot_beams(signal[:,i],beams_name,ra,dec,t_min,t_min+36)
        
      [ os.remove(path+f) for f in os.listdir(path) if f.endswith(".pnx") ]

  return

