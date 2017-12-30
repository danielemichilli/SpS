filename = 


import filterbank as fb
import os
import numpy as np
import sigproc

head = fb.read_header(filename)[0]
nbits = head['nbits']
nchans = head['nchans']
dtype = fb.get_dtype(nbits)
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

spectra = np.memmap(filename, dtype=dtype, mode='readwrite', offset=header_size, shape=(nspec, nchans))
spectra = np.fliplr(spectra)






filename = 


import pyfits

fits = pyfits.open(filename,memmap=True)
subint = fits['SUBINT'].data[]['DATA']
N_channels = fits['SUBINT'].header['NCHAN']
N_spectra = fits['SUBINT'].header['NSBLK']
subint = subint.reshape(N_spectra,N_channels)





