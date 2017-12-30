#############################
#
# LOTAAS Single Pulse plots
#
# Written by Daniele Michilli
#
#############################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tarfile
import os
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import logging
try: import presto
except ImportError: pass

import Utilities
from Parameters import *
import Paths
import RFIexcision

mpl.rc('font',size=5)


def sp_plot(pulses,rfi,meta_data,sap,beam,store):
  plt.clf()

  fill = u'b'
  square = u'g'
    
  fig = plt.figure()
  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=4,rowspan=2)
  ax2 = plt.subplot2grid((3,4),(0,0))
  ax3 = plt.subplot2grid((3,4),(0,1))
  ax4 = plt.subplot2grid((3,4),(0,2))
  ax5 = plt.subplot2grid((3,4),(0,3))

  scatter_beam(ax1,pulses,rfi=rfi)
  try: hist_DM(ax2,pulses)
  except ValueError: pass
  scatter_SNR(ax3,pulses)
  try: hist_SNR(ax4,pulses)
  except ValueError: pass
  meta_data_plot(ax5,meta_data)
  
  fig.tight_layout()
  plt.savefig(store,format='png',bbox_inches='tight',dpi=200)
  plt.close('all')
  return


def sp_shape(pulses,store,folder,idL):
  plt.clf()
  fig = plt.figure()
  
  events = pd.read_hdf('{}/SinglePulses.hdf5'.format(folder),'events',where=['Pulse==pulses.index.tolist()'])
  for i,(idx,puls) in enumerate(pulses.iterrows()):
    event = events[events.Pulse==idx]
    ax = plt.subplot2grid((2,5),(i/5,i%5))
    puls_DM_Time(ax,event,puls)
    ax.set_title('Sigma = {0:.1f}, Rank = {1}'.format(puls.Sigma,i))
    
  # Set common labels
  fig.text(0.5, 0.05, 'Time (s)', ha='center', va='center', fontsize=8)
  fig.text(0.08, 0.5, 'DM (pc/cm3)', ha='left', va='center', rotation='vertical', fontsize=8)
  fig.text(0.5, 0.95, '{} - SAP{}_BEAM{}'.format(idL,pulses.SAP.unique()[0],pulses.BEAM.unique()[0]), ha='center', va='center', fontsize=12)
  
  plt.savefig(store,format='png',bbox_inches='tight',dpi=200)
  plt.close('all')
  return


def obs_top_candidates(top_candidates,store,incoherent=False):
  plt.clf()
  
  if incoherent:
    col = top_candidates.SAP
    num = top_candidates.SAP.unique().size
  else:
    col = top_candidates.BEAM
    num = top_candidates.BEAM.unique().size
  
  if num > 1: cmap = discrete_cmap(num, 'spectral')
  else: cmap = plt.get_cmap('gist_heat_r')
  fig = plt.figure()
  plt.title("Best Candidates")
  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=4,rowspan=2)
  ax2 = plt.subplot2grid((3,4),(0,0))
  ax3 = plt.subplot2grid((3,4),(0,1))
  ax4 = plt.subplot2grid((3,4),(0,2))
  
  scatter_beam(ax1,top_candidates,cmap=cmap,col=col,legend=num)
  
  try: hist_DM(ax2,top_candidates)
  except ValueError: pass
  scatter_SNR(ax3,top_candidates,cmap=cmap,col=col,with_legend=True)
  try: hist_SNR(ax4,top_candidates)
  except ValueError: pass
  
  fig.tight_layout()
  plt.savefig(store,format='png',bbox_inches='tight',dpi=200)
  plt.close('all')
  return


def single_candidates(pulses,cands,meta_data,folder,idL):
  pulses_top = pulses[pulses.Pulse==0]
  rfi = pulses[pulses.Pulse>0]
  idx = pulses[pulses.Candidate.isin(cands.index)].groupby('Candidate',sort=False).head(1).index
  events = pd.read_hdf('{}/SinglePulses.hdf5'.format(folder),'events',where=['Pulse==idx.tolist()'])
  for idx,cand in cands.iterrows():
    puls = pulses[pulses.Candidate==idx]
    event = events[events.Pulse==puls.index[0]]
    sap = int(cand.SAP)
    beam = int(cand.BEAM)

    pulses_beam = pulses_top[(pulses_top.SAP==sap)&(pulses_top.BEAM==beam)]
    rfi_beam = rfi[(rfi.SAP==sap)&(rfi.BEAM==beam)]
    
    plt.clf()
    fig = plt.figure()

    ax1 = plt.subplot2grid((2,3),(0,0),colspan=2)
    ax2 = plt.subplot2grid((2,3),(0,2))
    ax3 = plt.subplot2grid((2,3),(1,0))
    ax4 = plt.subplot2grid((2,3),(1,1))
    ax5 = plt.subplot2grid((2,3),(1,2))

    scatter_beam(ax1,pulses_beam,rfi=rfi_beam)
    if cand.Rank>0: ax1.scatter(puls.Time, puls.DM, s=np.clip(np.log(puls.Sigma-5.5)*400+100,100,1200), linewidths=[0.,], c=u'b')
    ax1.scatter(puls.Time, puls.DM, s=300, linewidths=[0.,], marker='*', c='w')
    meta_data_puls(ax2,meta_data[(meta_data.SAP==sap)&(meta_data.BEAM==beam)],puls,cand)
    ax3.plot(event.DM, event.Sigma, 'k')
    ax3.set_xlabel('DM (pc/cm3)')    
    ax3.set_ylabel('SNR')
    puls_DM_Time(ax4,event,puls)
    DynamicSpectrum(ax5,puls.iloc[0].copy(),idL,sap,beam)
    
    plt.tight_layout()
    store = '{}/candidates'.format(folder)
    plt.savefig('{}/{}.png'.format(store,cand.id),format='png',bbox_inches='tight',dpi=200)
    plt.close('all')  #Maybe possible to reuse the figure without close it, should be faster
  return



def repeated_candidates(pulses,cands,meta_data,folder,idL):
  pulses_top = pulses[pulses.Pulse==0]
  rfi = pulses[pulses.Pulse>0]
  idx = pulses[pulses.Candidate.isin(cands.index)].groupby('Candidate',sort=False).head(2).index
  events = pd.read_hdf('{}/SinglePulses.hdf5'.format(folder),'events',where=['Pulse==idx.tolist()'])
  for idx,cand in cands.iterrows():
    puls1 = pulses[pulses.Candidate==idx].iloc[0]
    puls2 = pulses[pulses.Candidate==idx].iloc[1]
    event1 = events[events.Pulse==puls1.name]
    event2 = events[events.Pulse==puls2.name]
    sap = int(cand.SAP)
    beam = int(cand.BEAM)
    
    pulses_beam = pulses_top[(pulses_top.SAP==sap)&(pulses_top.BEAM==beam)]
    rfi_beam = rfi[(rfi.SAP==sap)&(rfi.BEAM==beam)]
    
    plt.clf()
    fig = plt.figure()

    gs = gridspec.GridSpec(2,4)
    ax1 = plt.subplot(gs[0,0:3])
    ax2 = plt.subplot(gs[0,3])
    gsA = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=gs[1,0], wspace=0.0)    
    axA1 = plt.subplot(gsA[0,0])
    axA2 = plt.subplot(gsA[0,1],sharey=axA1)
    axA2.label_outer()
    gsB = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=gs[1,1], hspace=0.0)    
    axB1 = plt.subplot(gsB[0,0])
    axB2 = plt.subplot(gsB[1,0],sharex=axB1)
    axB1.label_outer()
    gsC = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=gs[1,2], wspace=0.0, hspace=0.0)    
    axC1 = plt.subplot(gsC[0,0])
    axC2 = plt.subplot(gsC[0,1],sharey=axC1)
    axC2.label_outer()
    gsD = gridspec.GridSpecFromSubplotSpec(3,1,subplot_spec=gs[1,3], hspace=0.0)    
    axD1 = plt.subplot(gsD[0,0])
    axD2 = plt.subplot(gsD[1,0],sharex=axD1)
    axD3 = plt.subplot(gsD[2,0],sharex=axD1)
    axD1.label_outer() 
    axD2.label_outer() 
    
    puls_DM_Time(axA1,event1,puls1)
    puls_DM_Time(axA2,event2,puls2,sharey=True)
    axB1.plot(event1.DM, event1.Sigma, 'k')
    axB1.set_ylabel('SNR')
    axB2.plot(event2.DM, event2.Sigma, 'k')
    axB2.set_xlabel('DM (pc/cm2')
    axB2.set_ylabel('SNR')
    DynamicSpectrum(axC1,puls1.copy(),idL,sap,beam)
    DynamicSpectrum(axC2,puls2.copy(),idL,sap,beam,sharey=True)
    meta_data_repeat(ax2,meta_data[(meta_data.SAP==sap)&(meta_data.BEAM==beam)],cand,pulses[pulses.Candidate==idx])
    
    pulses_cand = pulses[(pulses.Pulse>0)&(pulses.Candidate==idx)]
    
    scatter_beam(ax1,pulses_beam,rfi=rfi_beam)
    scatter_beam(ax1,pulses_cand,col='b')
    ax1.scatter(pulses.Time[pulses.Candidate==idx], pulses.DM[pulses.Candidate==idx], s=300, linewidths=[0.,], marker='*', c='w')
    ax1.axhline(cand.DM,ls='--')
    if not pulses_beam.empty: scatter_SNR(axD2,pulses_beam)
    try: 
      hist_DM(axD1,pulses_beam)
      hist_SNR(axD3,pulses_beam)
    except ValueError: pass
    
    plt.tight_layout()
    store = '{}/candidates'.format(folder)
    plt.savefig('{}/{}.png'.format(store,cand.id),format='png',bbox_inches='tight',dpi=200)
    plt.close('all')  #Maybe possible to reuse the figure without close it, should be faster
  return






def scatter_beam(ax,pulses,cmap='gist_heat_r',col=None,rfi=False,legend=False):
  sig = np.clip(np.log(pulses.Sigma-5.5)*400+100,100,1200)
  if isinstance(col,type(None)): 
    col = sig
    vmin = col.min()/10.
  else: vmin = None
  if legend:
    main_plt = ax.scatter(pulses.Time, pulses.DM, c=col, s=sig, cmap=cmap, linewidths=[0.,],vmin=vmin)
    ticks = np.linspace(col.min(),col.max(),num=legend)
    bar = plt.colorbar(main_plt,ticks=ticks,ax=ax)
    bar.set_ticklabels(['{0:.0f}'.format(int(t)) for t in ticks])
    bar.ax.set_xlabel('sap',ha='left',labelpad=10)
    bar.update_ticks
    bar.ax.xaxis.set_ticks_position('top')
  else: ax.scatter(pulses.Time, pulses.DM, c=col, s=sig, cmap=cmap, linewidths=[0.,],vmin=vmin)
  if isinstance(rfi,pd.DataFrame): ax.scatter(rfi.Time, rfi.DM, s=5, c=u'k', marker='o', linewidths=[0.,])
  
  ax.axhline(40.48,c='k',ls='--',lw=.1)
  ax.axhline(141.68,c='k',ls='--',lw=.1)
  ax.set_yscale('log')
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('DM (pc/cm3)')
  ax.axis([0,3600,DM_MIN,550])
  
  top = pulses[pulses.Pulse==0].iloc[:10]
  for i in range(top.shape[0]):
    ax.annotate(str(i),xy=(top.Time.iloc[i],top.DM.iloc[i]),horizontalalignment='center',verticalalignment='center',color='dodgerblue',size=7,weight='bold')
  ax.tick_params(which='both',direction='out')
  
  return


def hist_DM(ax,pulses):
  ax.hist(pulses.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',range=(DM_MIN,550))
  ax.set_xscale('log')
  ax.set_xlabel('DM (pc/cm3)')
  ax.set_ylabel('Counts')
  ax.set_xlim(DM_MIN,550)
  ax.axvline(40.48,c='k',ls='--',lw=.1)
  ax.axvline(141.68,c='k',ls='--',lw=.1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return


def scatter_SNR(ax,pulses,cmap='gist_heat_r',col=None,with_legend=False):
  if isinstance(col,type(None)): 
    col = np.clip(np.log(pulses.Sigma-5.5)*400+100,100,1200)
    vmin = col.min()/10.
  else: vmin = None
  if with_legend: ax.scatter(pulses.DM,pulses.Sigma,c=col,s=6.,cmap=cmap,linewidths=[0.,],vmin=vmin)
  else: ax.scatter(pulses.DM,pulses.Sigma,c=col,s=3.,cmap=cmap,linewidths=[0.,],vmin=vmin)
  ax.set_xscale('log')
  ax.set_ylabel('SNR')
  ax.set_xlabel('DM (pc/cm3)')
  ax.axvline(40.48,c='k',ls='--',lw=.1)
  ax.axvline(141.68,c='k',ls='--',lw=.1)
  top = pulses.iloc[:10]
  for i in range(0,top.shape[0]):
    ax.annotate(i,xy=(top.DM.iloc[i]/1.15,top.Sigma.iloc[i]),horizontalalignment='right',verticalalignment='center',size=4,weight='medium')
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return


def hist_SNR(ax,pulses):
  ax.hist(pulses.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',weights=pulses.Sigma.tolist(),range=(DM_MIN,550))
  ax.set_xscale('log')
  ax.set_xlabel('DM (pc/cm3)')
  ax.set_ylabel('Cumulative SNR')
  ax.set_xlim(DM_MIN,550)
  ax.axvline(40.48,c='k',ls='--',lw=.1)
  ax.axvline(141.68,c='k',ls='--',lw=.1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return


def meta_data_plot(ax,meta_data):
  ax.axis([0,10,0,7])
  ax.annotate('File: '+meta_data.File.iloc[0], xy=(0,6))
  ax.annotate('Telescope: '+meta_data.Telescope.iloc[0], xy=(0,5))
  ax.annotate('Instrument: '+meta_data.Instrument.iloc[0], xy=(0,4))
  ax.annotate('RA: '+meta_data.RA.iloc[0], xy=(0,3))
  ax.annotate('DEC: '+meta_data.DEC.iloc[0], xy=(0,2))
  ax.annotate('Epoch (MJD): '+meta_data.Epoch.iloc[0], xy=(0,1))
  ax.axis('off')
  return

def meta_data_puls(ax,meta_data,puls,cand):
  ax.axis([0,10,0,8])
  ax.annotate('File: {}'.format(meta_data.File.iloc[0]), xy=(0,8))
  ax.annotate('RA, DEC: {0:.8}, {1:.8}'.format(meta_data.RA.iloc[0],meta_data.DEC.iloc[0]), xy=(0,7))
  ax.annotate('Epoch (MJD): {0:.11}'.format(meta_data.Epoch.iloc[0]), xy=(0,6))
  ax.annotate('DM (pc/cm2): {0:.2f}'.format(puls.DM.iloc[0]), xy=(0,5))
  ax.annotate('Time (s): {0:.2f}'.format(puls.Time.iloc[0]), xy=(0,4))
  ax.annotate('Sigma: {0:.1f}'.format(puls.Sigma.iloc[0]), xy=(0,3))
  ax.annotate('Duration (ms): {0:.0f}'.format(puls.Duration.iloc[0]*1000), xy=(0,2))
  ax.annotate('N. events: {}'.format(puls.N_events.iloc[0]), xy=(0,1))
  ax.annotate('Rank: {0:.0f}'.format(cand.Rank), xy=(0,0))
  ax.axis('off')
  return

def meta_data_repeat(ax,meta_data,cand,pulses):
  ax.axis([0,10,0,9])
  ax.annotate('File: {}'.format(meta_data.File.iloc[0]), xy=(0,9))
  ax.annotate('RA, DEC: {0:.8}, {1:.8}'.format(meta_data.RA.iloc[0],meta_data.DEC.iloc[0]), xy=(0,8))
  ax.annotate('Epoch (MJD): {0:.11}'.format(meta_data.Epoch.iloc[0]), xy=(0,7))
  ax.annotate('DM (pc/cm2): {0:.2f}'.format(cand.DM), xy=(0,6))
  ax.annotate('dDM (pc/cm2): {0:.2f}'.format(pulses.DM.max()-pulses.DM.min()), xy=(0,5))
  ax.annotate('Period (s): {0:.3f}'.format(cand.Period), xy=(0,4))
  ax.annotate('Period err. (s): {0:.3f}'.format(cand.Period_err), xy=(0,3))
  ax.annotate('Sigma (cum.): {0:.1f}'.format(cand.Sigma), xy=(0,2))
  ax.annotate('N. pulses: {0:.0f}'.format(cand.N_pulses), xy=(0,1))
  ax.annotate('Rank: {0:.0f}'.format(cand.Rank), xy=(0,0))
  ax.axis('off')
  return

def puls_DM_Time(ax,event,puls,sharey=False):
  #sig = (event.Sigma/event.Sigma.max()*5)**4
  #sig = np.clip(np.log(event.Sigma-5.5)*400+100,100,1200)
  #sig = event.Sigma/event.Sigma.max()*1000
  sig = np.clip(event.Sigma/event.Sigma.max()*1427-427,1,1000)
  ax.scatter(event.Time, event.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])  
  ax.errorbar(puls.Time, puls.DM, xerr=puls.dTime/2, yerr=puls.dDM/2, fmt='none', ecolor='r')
  ax.set_xlabel('Time (s)')  
  if not sharey:
    ax.set_ylabel('DM (pc/cm3)')
  return


def discrete_cmap(N, base_cmap):
  base = plt.cm.get_cmap(base_cmap)
  color_list = base(np.linspace(0, 1, N))
  cmap_name = base.name + str(N)
  return base.from_list(cmap_name, color_list, N)  
  
  

def DynamicSpectrum(ax1,puls,idL,sap,beam,sharey=False):
  if beam==12: stokes = 'incoherentstokes'
  else: stokes = 'stokes'
  filename = '{folder}/{idL}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=Paths.RAW_FOLDER,idL=idL,stokes=stokes,sap=sap,beam=beam)
  if not os.path.isfile(filename): return
  
  if puls.DM>141.71: sample = puls.Sample * 4
  elif puls.DM>40.47: sample = puls.Sample * 2
  else: sample = puls.Sample

  #controllare che questo vada dopo downsampling correction!
  filetype, header = Utilities.read_header(filename)
  MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400.
  try: v = presto.get_baryv(header['RA'],header['DEC'],MJD,1800.,obs='LF')
  except NameError: 
    logging.warning("LSPplot - Additional modules missing")
    return
  sample += np.round(sample*v).astype(int)
  
  duration = np.int(np.round(puls.Duration/RES))
  spectra_border = 20
  offset = duration*spectra_border
  
  #Load the spectrum
  spectrum = Utilities.read_fits(filename,puls.DM.copy(),sample.copy(),duration,offset,RFI_reduct=True)
  
  #De-dispersion
  freq = np.linspace(F_MIN,F_MAX,2592)
  time = (4149 * puls.DM * (F_MAX**-2 - np.power(freq,-2)) / RES).round().astype(np.int)
  for i in range(spectrum.shape[1]):
    spectrum[:,i] = np.roll(spectrum[:,i], time[i])
  spectrum = spectrum[:2*offset+duration]
  
  spectrum = np.mean(np.reshape(spectrum,[spectrum.shape[0]/duration,duration,spectrum.shape[1]]),axis=1)
  spectrum = np.mean(np.reshape(spectrum,[spectrum.shape[0],spectrum.shape[1]/32,32]),axis=2)
  
  extent = [(sample-offset)*RES,(sample+duration+offset)*RES,F_MIN,F_MAX]
  ax1.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent)
  ax1.scatter((sample+duration/2)*RES,F_MIN+1,marker='^',s=1000,c='r')
  ax1.axis(extent)
  ax1.set_xlabel('Time (s)')
  if not sharey:
    ax1.set_ylabel('Frequency (MHz)')
      
  return 
