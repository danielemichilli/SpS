idL = 'L341466'
sap = 2
dm_min = 73.85
dm_max = 74.25
pulsar = 'J1937+15'

import pandas as pd
import os
import numpy as np
import tarfile
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


#ra = np.array([ 
        #499,  499,  622,  621,  499,  377,  376,  499,  624,  748,  744,
        #742,  620,  499,  379,  257,  254,  251,  375,  499,  625,  750,
        #873,  869,  865,  862,  740,  619,  499,  380,  259,  137,  134,
        #129,  126,  249,  374,  499,  627,  752,  876, 1000,  995,  990,
        #985,  981,  859,  738,  618,  499,  381,  261,  140,   18,   14,
          #9,    4,    0,  123,  247,  372])

 
#dec = np.array([ 
        #500,  625,  562,  437,  375,  437,  562,  750,  687,  625,  500,
        #375,  312,  250,  312,  375,  500,  625,  687,  875,  812,  750,
        #687,  562,  437,  312,  250,  187,  125,  187,  250,  312,  437,
        #562,  687,  750,  812, 1000,  937,  875,  812,  750,  625,  500,
        #375,  250,  187,  125,   62,    0,   62,  125,  187,  250,  375,
        #500,  625,  750,  812,  875,  937])

folder = '/projects/0/lotaas/data/out/new'

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
    print 'Incomplete observation'
    ra[k] = np.nan
    dec[k] = np.nan
  
  k += 1


  
#ra = ra[~np.isnan(ra)]
#dec = dec[~np.isnan(dec)]




pulses = pd.read_hdf('SinglePulses.hdf5',idL+'_pulses')
pulses = pulses[(pulses.SAP==sap)&(pulses.Pulse>0)&(pulses.Pulse<4)&(pulses.N_events>4)&(pulses.DM>dm_min)&(pulses.DM<dm_max)]

beams = pulses[pulses.BEAM!=12].groupby('BEAM').Sigma.sum()
 
ind = pd.Series(np.zeros(61))
ind.index += 13
beams = beams.reindex_like(ind)
#beams.dropna(inplace=True)

plt.clf()

plt.scatter(ra,dec,s=800,edgecolor='none',c=beams,cmap=mpl.cm.hot_r)
bar = plt.colorbar()
bar.set_label('Cumulative SNR')
plt.xlabel('RA (h)')
plt.ylabel('DEC (deg)')
[plt.annotate(str(i+13),(ra[i],dec[i]),horizontalalignment='center',verticalalignment='center') for i in range(0,61)]

plt.annotate('DM: {} - {}'.format(dm_min,dm_max), xy=(np.nanmin(ra)-0.02,np.nanmax(dec)+0.1))

plt.xlim(np.nanmin(ra)-0.025,np.nanmax(ra)+0.025)
plt.ylim(np.nanmin(dec)-0.25,np.nanmax(dec)+0.25)

plt.savefig('{}_beams.png'.format(pulsar),format='png',bbox_inches='tight',dpi=200)

