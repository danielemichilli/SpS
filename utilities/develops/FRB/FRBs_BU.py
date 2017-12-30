#!/usr/bin/env python

'''

Fast Imager

Written by Daniele Michilli

'''


import pandas as pd
import os
import argparse
import tarfile
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




def plot_beams(beams,name,ra,dec,t_min,t_max):
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




parser = argparse.ArgumentParser()
parser.add_argument('filenm', nargs=1)
args = parser.parse_args()
idL = args.filenm[0]

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









