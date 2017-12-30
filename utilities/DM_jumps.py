idL = ''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

events = pd.read_hdf('SinglePulses.hdf5',idL)

y = np.array((0.0, 2.53, 5.06, 7.59, 10.12, 12.65, 15.18, 17.71, 20.24, 22.77, 25.3, 27.83, 30.36, 32.89, 35.42, 37.95, 40.52, 65.82, 91.12, 116.42, 141.77, 242.97, 344.17, 445.37))

sig = (events.Sigma/events.Sigma.max()*5)**4
plt.scatter(events.Time, events.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])
plt.plot([0,3600],[y,y],'r--')

plt.show()

exit()

