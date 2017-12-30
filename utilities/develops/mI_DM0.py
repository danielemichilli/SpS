#!/usr/bin/env python

'''


Written by Daniele Michilli

'''

folder = '/projects/lotaas/test/Daniele/raw/L204720_red/stokes/SAP2/BEAM38/L204720_SAP2_BEAM38.fits'

import numpy as np
import pyfits
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filenm', nargs=1)
args = parser.parse_args()
idL = args.filenm[0]

sap = 2
folder = '/projects/lotaas/test/Daniele/raw/'

for beam in range(13,74):
  name = '{}_SAP{}_BEAM{}'.format(idL,sap,beam)
  path = '{}{}_red/stokes/SAP{}/BEAM{}/'.format(folder,idL,sap,beam)
  file_name = '{}{}.fits'.format(path,name)


  fits = pyfits.open(file_name,memmap=True)
  subint = fits['SUBINT'].data['DATA']

  N_channels = fits['SUBINT'].header['NCHAN']
  N_spectra = fits['SUBINT'].header['NSBLK']*fits['SUBINT'].header['NAXIS2']

  subint = subint.reshape(N_spectra,N_channels)

  subint = subint[:,~np.all(subint == 0, axis=0)]
  N_channels = subint.shape[1]

  downfactor = 10
  N_spectra /= downfactor
  subint = subint.reshape(N_spectra,downfactor,N_channels).mean(axis=1)

  subint = np.vstack((subint,np.zeros((739200-N_spectra,N_channels))))
  N_spectra = 739200

  downfactor = 7392
  N_spectra /= downfactor
  subint = subint.reshape(N_spectra,downfactor,N_channels).mean(axis=1)

  lim = np.partition(subint,-20,axis=1)[:,-20]
  for idx,row in enumerate(subint):
    row[row>=lim[ind]] = np.median(row)

  downfactor = N_channels/3
  N_channels /= downfactor 
  subint = subint.reshape(N_spectra,N_channels,downfactor).mean(axis=2) 

  subint = subint.std(axis=1)

  #studiare max in osservazioni:
  subint.sort()
  print subint[-10:]

