

import pandas as pd


pulses = pd.read_hdf('SinglePulses.hdf5',idL+'_pulses',where=['Pulse==0'])
events = pd.read_hdf('SinglePulses.hdf5',idL,where=['Pulse=pulses.index.tolist()'])

#long = pulses[(pulses.DM>(DM-0.2))&(pulses.DM<(DM+0.2))&(pulses.Sigma>6.5)]
#long.sort('N_events',inplace=True)
#print long[long.N_events<10].shape[0]/float(long.shape[0])*100


long = pulses[(pulses.DM>(DM-0.2))&(pulses.DM<(DM+0.2))&(pulses.N_events>=10)]
long.sort('Sigma',inplace=True)
#print long[long.N_events<10].shape[0]/float(long.shape[0])*100





import matplotlib.pyplot as plt

def sp_shape(pulses,events,store):

  fig = plt.figure()

  for i in range(0,pulses.shape[0]):
  
    puls = pulses.iloc[i]
    event = events.loc[events.Pulse==puls.name]

    sig = (event.Sigma/event.Sigma.max()*5)**4
  
    ax = plt.subplot2grid((2,5),(i/5,i%5))
    ax.scatter(event.Time, event.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])  
    ax.errorbar(puls.Time_c, puls.DM_c, xerr=puls.dTime, yerr=puls.dDM, fmt='none', ecolor='r')
    
    ax.set_title('Sigma = {0:.1f}, Rank = {1}'.format(event.Sigma.max(),i))
    
  
  # Set common labels
  fig.text(0.5, 0.05, 'Time (s)', ha='center', va='center', fontsize=8)
  fig.text(0.08, 0.5, 'DM (pc/cm3)', ha='left', va='center', rotation='vertical', fontsize=8)
  
  plt.savefig('{}'.format(store),format='png',bbox_inches='tight',dpi=200)
  
  return



def event_shape(dmlo,dmhi,tmlo,tmhi,events):
  
  selected = events[(events.DM<dmhi)&(events.DM>dmlo)&(events.Time<tmhi)&(events.Time>tmlo)]
  
  sig = (selected.Sigma/selected.Sigma.max()*5)**4


  plt.scatter(selected.Time, selected.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])    
  
  # Set common labels
  plt.text(0.5, 0.05, 'Time (s)', ha='center', va='center', fontsize=8)
  plt.text(0.08, 0.5, 'DM (pc/cm3)', ha='left', va='center', rotation='vertical', fontsize=8)
 
  return






def plot_pulse(idx):
  event = pd.read_hdf('SinglePulses.hdf5',idL,where=['Pulse==idx'])
  sig = (event.Sigma/event.Sigma.max()*10)**4
  plt.figure(figsize=(20,20))
  plt.xlim(event.Time.min()-0.5,event.Time.max()+0.5)
  plt.scatter(event.Time, event.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])
  print pulses.loc[idx]
  plt.show()
   
   
   
def plot_pulse(idx):
  event = events[events.Pulse==idx]
  sig = (event.Sigma/event.Sigma.max()*10)**4
  plt.figure(figsize=(20,20))
  plt.xlim(event.Time.min()-0.5,event.Time.max()+0.5)
  plt.scatter(event.Time, event.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])
  plt.show()
  

