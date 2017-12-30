idL = ''
astro = []
sap = 2
beam = 45



import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt



pulses = pd.read_hdf('SinglePulses.hdf5',idL+'_pulses')
pulses = pulses[(pulses.SAP==sap)&(pulses.BEAM==beam)]
rrat = pulses[pulses.index.isin(astro)]
rrat.sort('Sigma',ascending=False,inplace=True)
pulses = pulses[(pulses.Pulse==0) & (~ pulses.index.isin(rrat.index))]
meta_data = pd.read_hdf('SinglePulses.hdf5','meta_data')
meta_data = meta_data[(meta_data.SAP==str(sap))&(meta_data.BEAM==str(beam))]



def sp_plot(pulses,rrat,meta_data,store):

  col = pulses.Sigma
  cmap = plt.get_cmap('gist_heat_r')
  

  fig = plt.figure()
  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=4,rowspan=2)
  ax2 = plt.subplot2grid((3,4),(0,0))
  ax3 = plt.subplot2grid((3,4),(0,1))
  ax4 = plt.subplot2grid((3,4),(0,2))
  ax5 = plt.subplot2grid((3,4),(0,3))

  ax1.scatter(pulses.Time, pulses.DM, c=col, s=20., cmap=cmap, linewidths=[0.,], vmin=5, vmax=10) 
  
  
  
  mpl.rc('font',size=5)
  sig = (rrat.Sigma/rrat.Sigma.max()*30)**2
  ax1.scatter(rrat.Time, rrat.DM, s=sig, c=u'b', linewidths=[0.,])
  
  
  
  DM_m = rrat.DM.mean()
  ax1.plot([0,3600],[DM_m,DM_m],'r:')

  ax1.set_yscale('log')
  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('DM (pc/cm3)')
  ax1.axis([0,3600,5,550])
  
  mpl.rc('font', size=5)
  for i in range(0,rrat.shape[0]):
    ax1.annotate(i,xy=(rrat.Time.iloc[i],rrat.DM.iloc[i]*1.05),horizontalalignment='right',verticalalignment='center')
  
  all_pulses = pd.concat([pulses,rrat])
  ax2.hist(all_pulses.Sigma.tolist(),bins=100,histtype='step',color='k')
  ax2.set_xlabel('SNR')
  ax2.set_ylabel('Counts')
  ax2.set_yscale('log')

  hist = ax3.hist(all_pulses.DM.tolist(),bins=300,histtype='stepfilled',color=u'k')
  ax3.set_xscale('log')
  ax3.set_xlabel('DM (pc/cm3)')
  ax3.set_ylabel('Counts')
  ax3.set_xlim(5,550)
  ax3.plot([DM_m,DM_m],[0,hist[0].max()+10],'r:')
  
  ax4.scatter(pulses.DM,pulses.Sigma,c=col,s=3.,cmap=cmap,linewidths=[0.,],vmin=5,vmax=10)
  ax4.scatter(rrat.DM,rrat.Sigma,s=15.,linewidths=[0.,],c=u'g',marker='*')
  ax4.set_xscale('log')
  ax4.set_ylabel('SNR')
  ax4.set_xlabel('DM (pc/cm3)')
  limit = max(pulses.Sigma.max(),rrat.Sigma.max())
  ax4.axis([5,550,pulses.Sigma.min(),limit+3.])
  ax4.plot([DM_m,DM_m],[0,limit+3.],'r:')
  mpl.rc('font', size=3.5)
  for i in range(0,rrat.shape[0]):
    ax4.annotate(i,xy=(rrat.DM.iloc[i]/1.15,rrat.Sigma.iloc[i]),horizontalalignment='right',verticalalignment='center')

    
  mpl.rc('font', size=5)
  ax5.axis([0,10,0,7])
  ax5.annotate('File: '+meta_data.File.iloc[0], xy=(0,6))
  ax5.annotate('Telescope: '+meta_data.Telescope.iloc[0], xy=(0,5))
  ax5.annotate('Instrument: '+meta_data.Instrument.iloc[0], xy=(0,4))
  ax5.annotate('RA: '+meta_data.RA.iloc[0], xy=(0,3))
  ax5.annotate('DEC: '+meta_data.DEC.iloc[0], xy=(0,2))
  ax5.annotate('Epoch (MJD): '+meta_data.Epoch.iloc[0], xy=(0,1))
  ax5.axis('off')

  ax1.tick_params(which='both',direction='out')
  ax2.tick_params(which='both',direction='out')
  ax3.tick_params(which='both',direction='out')
  ax4.tick_params(which='both',direction='out')
  ax5.tick_params(which='both',direction='out')
  
  fig.tight_layout()
  mpl.rc('font',size=5)
  plt.savefig('{}'.format(store),format='png',bbox_inches='tight',dpi=200)
  
  return
