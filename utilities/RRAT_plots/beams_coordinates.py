idL = 'L204720'
sap = 1

import os
import numpy as np
import tarfile

folder = '/projects/lotaas/data/out/'

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



