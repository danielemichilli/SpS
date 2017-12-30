import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from RFIexcision import Pulse_Thresh

pulses = pd.read_hdf('/var/run/media/michilli/DATA/LOTAAS/pulsars/DB/DB.h5','pulses')
events = pd.read_hdf('/var/run/media/michilli/DATA/LOTAAS/pulsars/DB/DB.h5','events')
Pulse_Thresh(pulses,events)




plt.hist((pulses.Sigma[pulses.Pulse==0].tolist(),pulses.Sigma[pulses.Pulse==1].tolist(),pulses.Sigma[pulses.Pulse==2].tolist(),pulses.Sigma[pulses.Pulse>2].tolist()),bins=50,histtype='barstacked',log=True,label=['Excellent','Good','Poor','RFI'],color=['g','y','orange','red'])
plt.legend()
plt.xlim((6.5,12))
plt.xlabel('SNR')
plt.ylabel('counts')
plt.show()


plt.hist((pulses.Sigma[pulses.Pulse==0].tolist(),pulses.Sigma[pulses.Pulse==1].tolist(),pulses.Sigma[pulses.Pulse==2].tolist(),pulses.Sigma[pulses.Pulse>2].tolist()),bins=50,histtype='barstacked',label=['Excellent','Good','Poor','RFI'],color=['g','y','orange','red'])
plt.legend()
plt.xlim((6.5,12))
plt.xlabel('SNR')
plt.ylabel('counts')
plt.show()
