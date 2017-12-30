idL = 'L352194'

dm_min = 18 #18.95 77.1
dm_max = 20 #19.04 77.7
pulsar = 'J0301+20'

import pandas as pd
import os
import numpy as np
import tarfile
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

sap = 0

folder = '/projects/0/lotaas/data/out/confirmations/new/'

ra = np.zeros(128)
dec = np.zeros(128)

k = 0
for beam in range(0,128):
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
    print 'SAP',sap,'BEAM',beam,' missing'
    ra[k] = np.nan
    dec[k] = np.nan
  
  k += 1



pulses = pd.read_hdf('SinglePulses.hdf5','pulses')
pulses = pulses[(pulses.SAP==sap)&(pulses.Pulse>=0)&(pulses.Pulse<=2)&(pulses.DM>dm_min)&(pulses.DM<dm_max)]

beams = pulses.groupby('BEAM').Sigma.sum()

ind = pd.Series(np.zeros(128))
beams = beams.reindex_like(ind)

plt.clf()

plt.scatter(ra,dec,s=80,edgecolor='none',c=beams,cmap=mpl.cm.hot_r)
bar = plt.colorbar()
bar.set_label('Cumulative SNR')
plt.xlabel('RA (h)')
plt.ylabel('DEC (deg)')
[plt.annotate(str(i),(ra[i],dec[i]),horizontalalignment='center',verticalalignment='center') for i in range(0,128)]

plt.annotate('DM: {} - {}'.format(dm_min,dm_max), xy=(np.nanmin(ra)-0.02,np.nanmax(dec)+0.1))

plt.xlim(np.nanmin(ra)-0.025,np.nanmax(ra)+0.025)
plt.ylim(np.nanmin(dec)-0.25,np.nanmax(dec)+0.25)

#plt.savefig('{}_beams.png'.format(pulsar),format='png',bbox_inches='tight',dpi=200)

