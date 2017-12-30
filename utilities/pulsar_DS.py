filename = 'L204720_SAP1_BEAM48.fil'
pulsar_period = 0.714519699726 #s
pulsar_DM = 26.7641


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import filterbank as fb
import os
import presto
import sigproc
import string

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

data = np.memmap(filename, dtype=dtype, mode='r', offset=header_size, shape=(nspec, nchans))
data = np.fliplr(data)


#chunks = 16
#data = data.reshape(chunks,data.size/chunks)


#Define the process parameters
tsamp = head['tsamp']
period = pulsar_period / tsamp
window = int(period)




ra = string.split(str(head['src_raj']),'.')[0]
ra = ra[-6:-4]+':'+ra[-4:-2]+':'+ra[-2:]
dec = string.split(str(head['src_dej']),'.')[0]
dec = dec[-6:-4]+':'+dec[-4:-2]+':'+dec[-2:]

v1 = presto.get_baryv(ra,dec,head['tstart'],0,obs='LF')
v2 = presto.get_baryv(ra,dec,head['tstart'],3600,obs='LF')
v = (v1 + v2) / 2.
bin_idx = np.abs(np.int(np.round(1./v)))


ind = np.arange(0,data.shape[0],bin_idx)
#mean = np.mean(data)
if v<0: data = np.insert(data,ind,0,axis=0)
else: data = np.delete(data,ind,axis=0)




#offset = period-window
#ind = []
#for i in range(1,data.shape[0]/window):
  #offset += offset
  #if offset > 1: 
    #ind.append(window*i)
    #offset -= 0.9999999999


ind = np.arange(data.shape[0]/window)*period
ind = ind.astype(np.int)%window
ind = np.where(np.abs(np.diff(ind))>0)[0]*window+window
ind += np.arange(ind.size)

data = np.delete(data,ind,axis=0)





data = data[:data.shape[0]/window*window]
data = np.mean(data.reshape(-1, window, nchans),axis=0)

ind = np.arange(0,2592,16)
data = np.delete(data,ind,axis=1)


#Plot 

def color_range(data):
  #Define the color range
  clean = data[data>0]
  min_element = clean.size/20
  max_element = clean.size*9/10
  vmin = np.partition(clean, min_element, axis=None)[min_element]   #invece di subint[subint>0] possibile subint[:-(num_rows/down_fact)]
  vmax = np.partition(clean, max_element, axis=None)[max_element]
  return vmin,vmax

vmin,vmax = color_range(data)



#plot the DM lines
freq = np.arange(151,119,-1,dtype=np.float)
time = 4149 * pulsar_DM * (np.power(freq,-2) - 151.**-2) + 1
plt.plot(time,freq,'b-')



plt.figure(figsize=(20,10))
plt.xlabel('Time (s)')
plt.ylabel('Frequency (MHz)')
plt.imshow(data.T,cmap=mpl.cm.hot_r,origin="lower",aspect='auto',interpolation='nearest',extent=[0,t0+total_duration,119,151],vmin=vmin,vmax=vmax)
plt.colorbar()






