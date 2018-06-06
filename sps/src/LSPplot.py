from glob import glob
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np

from waterfaller import waterfaller, psrfits, psr_utils


import matplotlib as mpl
mpl.rc('font',size=5)


def output(args, events, pulses, candidates, meta_data):
  #--------------------------------------------
  # Produce diagnostic plots for the candidates
  #--------------------------------------------

  pulses.sort('Sigma', ascending=False, inplace=True)
  candidates.sort('Sigma', ascending=False, inplace=True)
  
  plt.close('all')
  fig = plt.figure(figsize=(7,8))

  for i, (idx_c, cand) in enumerate(candidates.iterrows()):
    plot_name = os.path.splitext(args.plot_name)[0] + '_' + str(i) + os.path.splitext(args.plot_name)[1]
    with PdfPages(plot_name) as pdf:
      pulses_cand = pulses[pulses.Candidate == idx_c]
      cand_plot(pdf, cand, pulses_cand, pulses, events, meta_data)
      # Loop over 5 brightest pulses
      for j, (idx_p, puls) in enumerate(pulses_cand.head(5).iterrows()):
        puls_plot(pdf, puls, events, j, args)

  return


def cand_plot(pdf, cand, pulses, pulses_all, events, meta_data):
  #---------------
  # Candidate plot
  #---------------

  gs = gridspec.GridSpec(3, 6, wspace=.7, hspace=.3)

  # Meta data
  ax1 = plt.subplot(gs.new_subplotspec((0,0), 2, 1))
  meta_data_plot(ax1, pulses, cand, meta_data)

  # DM vs Time
  ax2 = plt.subplot(gs.new_subplotspec((0,1), 2, 5))
  scatter_beam(ax2, pulses, pulses_all, cand)

  # S/N vs DM
  if not pulses_beam.empty: 
    ax4 = plt.subplot(gs.new_subplotspec((2,2), 1, 2))
    scatter_SNR(ax4, events[events.Pulse.isin(pulses_all.index)], cand)

  # Histogram of DM
  ax3 = plt.subplot(gs.new_subplotspec((2,0), 1, 2))
  # Cumulative S/N vs DM
  ax5 = plt.subplot(gs.new_subplotspec((2,4), 1, 2))
  try: 
    hist_DM(ax3, pulses_all, cand)
    hist_SNR(ax5, pulses_all, cand)
  except ValueError: pass
  
  pdf.savefig(bbox_inches='tight', dpi=200)
  return


def puls_plot(pdf, puls, events, i, args):
  #-----------
  # Pulse plot
  #-----------

  ev = events[events.Pulse == puls.name]  # Events in the pulse

  if args.fits is not None: col = 4
  else: col = 3
  gs = gridspec.GridSpecFromSubplotSpec(2, col, wspace=0.5, hspace=0.2)

  # Meta data
  ax1 = plt.subplot(gs.new_subplotspec((0,0), 1, 1))
  puls_meta_data(ax1, puls, ev.Pulse.iloc[0], i)

  # DM vs Time (zoom on pulse)
  ax2 = plt.subplot(gs.new_subplotspec((0,1), 1, 1))
  puls_DM_Time(ax2, ev, events, puls)

  # Time series at different DMs and pulse profile
  if args.timeseries is not None:
    ax3 = plt.subplot(gs.new_subplotspec((1,0), 1, 2))
    ax4 = plt.subplot(gs.new_subplotspec((0,2), 1, 1))
    puls_dedispersed(ax3, ax4, puls, args)

  # S/N vs DM
  ax5 = plt.subplot(gs.new_subplotspec((1,2), 1, 1))
  puls_SNR_DM(ax5, ev)

  # Dynamic spectra
  if args.fits is not None:
    ax6 = plt.subplot(gs.new_subplotspec((1,3), 1, 1))
    ax7 = plt.subplot(gs.new_subplotspec((0,3), 1, 1), sharey=ax6)
    if args.timeseries is None: ax4 = plt.subplot(gs.new_subplotspec((0,2), 1, 1))
    else: ax4 = None
    puls_dynSpec(ax6, ax7, puls, ax_prof=ax4)

  pdf.savefig(bbox_inches='tight', dpi=200)
  return


def meta_data_plot(ax, pulses, cand, meta_data):
  #-------------------------
  # Plot candidate meta data
  #-------------------------

  ax.axis([0,10,0,10])
  ax.annotate('Cand ID: {}'.format(cand.name), xy=(0,10))
  ax.annotate('File: {}'.format(meta_data.File.iloc[0]), xy=(0,9))
  if 'RA' in meta_data.columns: ax.annotate('RA, DEC: {0:.8}, {1:.8}'.format(meta_data.RA.iloc[0],meta_data.DEC.iloc[0]), xy=(0,8))
  if 'Epoch' in meta_data.columns: ax.annotate('Epoch (MJD): {0:.11}'.format(meta_data.Epoch.iloc[0]), xy=(0,7))
  ax.annotate('DM (pc cm$^{-3}$): {0:.2f}'.format(cand.DM), xy=(0,6))
  ax.annotate('Sigma (cum.): {0:.1f}'.format(cand.Sigma), xy=(0,5))
  if cand.N_pulses == 1:
    ax.annotate('Time (s): {0:.2f}'.format(pulses.Time.iloc[0]), xy=(0,3))
    ax.annotate('Width (ms): {0:.0f}'.format(pulses.Duration.iloc[0]*1000), xy=(0,2))
    ax.annotate('N. events: {}'.format(pulses.N_events.iloc[0]), xy=(0,1))
  else:
    ax.annotate('dDM (pc cm$^{-3}$): {0:.2f}'.format(pulses.DM.max()-pulses.DM.min()), xy=(0,4))
    ax.annotate('Period (s): {0:.3f}'.format(cand.Period), xy=(0,3))
    ax.annotate('Period err. (s): {0:.3f}'.format(cand.Period_err), xy=(0,2))
    ax.annotate('N. pulses: {0:.0f}'.format(cand.N_pulses), xy=(0,1))
  if 'version' in meta_data.columns: ax.annotate('L-SpS {}'.format(meta_data.version.iloc[0]), xy=(0,0))
  ax.axis('off')
  return


def scatter_beam(ax, pulses, pulses_beam, cand):
  #--------------------------
  # Plot pulses in DM vs Time
  #--------------------------

  # Define circle sizes
  def circle_size(values):
    new_val = np.clip(values,6.5,20)
    m = 31.85
    q = -137.025
    return new_val * m + q

  # Plot all the pulses
  sig = circle_size(pulses_beam.Sigma)
  ax.scatter(pulses_beam.Time, pulses_beam.DM, c='k', s=sig, edgecolor='w', lw=.2, zorder=2)

  # Highlight pulses in the candidate
  sig = circle_size(pulses.Sigma)
  ax.scatter(pulses.Time, pulses.DM, s=sig, linewidths=[0.,], marker='*', c='w', zorder=3)

  # Candidate DM
  ax.axhline(cand.DM, c='r', ls='-', zorder=1)
  
  # Name 5 brightest pulses in the candidate
  top = pulses.iloc[:5]
  for i in range(top.shape[0]):
    ax.annotate(str(i),xy=(top.Time.iloc[i],top.DM.iloc[i]),horizontalalignment='center',verticalalignment='center',color='dodgerblue',size=5,fontweight='bold',zorder=4)
  ax.tick_params(which='both',direction='out')

  ax.set_yscale('log')
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('DM (pc/cc)')
  
  return


def scatter_SNR(ax, events_beam, cand):
  #-------------------------
  # Plot events in S/N vs DM
  #-------------------------

  # Plot all the events
  ax.scatter(events_beam.DM, events_beam.Sigma, c='k', s=5, linewidths=[0.,], zorder=2)

  # Candidate DM
  ax.axvline(cand.DM, c='r', ls='-', linewidth=.2, zorder=3)

  ax.set_xscale('log')
  ax.set_ylabel('S/N')
  ax.set_xlabel('DM (pc cm$^{-3}$)')

  ax.set_ylim((6.5, events_beam.Sigma.max()+1))
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')

  return


def hist_DM(ax, pulses, cand):
  #-------------------------
  # Plot histogram of pulses
  #-------------------------

  # Histogram of pulses
  ax.hist(pulses.DM.tolist(), bins='auto', histtype='stepfilled', color=u'k', zorder=2)

  # Candidate DM
  ax.axvline(cand.DM, c='r', ls='-', linewidth=.2, zorder=3)

  ax.set_xscale('log')
  ax.set_xlabel('DM (pc cm$^{-3}$)')
  ax.set_ylabel('Counts')
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return


def hist_SNR(ax,pulses,cand):
  #----------------------------------------
  # Plot histogram of pulses cumulative S/N
  #----------------------------------------

  # Histogram of pulses
  ax.hist(pulses.DM.tolist(), weights=pulses.Sigma.tolist(), bins='auto', histtype='stepfilled', color=u'k', zorder=2)

  # Candidate DM
  ax.axvline(cand.DM, c='r', ls='-', linewidth=.2, zorder=3)

  ax.set_xscale('log')
  ax.set_xlabel('DM (pc cm$^{-3}$)')
  ax.set_ylabel('Cumulative S/N')
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return


def puls_meta_data(ax, puls, idx, i):
  #---------------------
  # Plot pulse meta data
  #---------------------

  ax.axis([0,10,0,10])
  ax.annotate("Pulse n. {}".format(i), xy=(0,9)) 
  ax.annotate('Pulse code: {}'.format(idx), xy=(0,8))
  ax.annotate('DM (pc/cm2): {0:.2f}'.format(puls.DM), xy=(0,7))
  ax.annotate('dDM (pc/cm2): {0:.2f}'.format(puls.dDM), xy=(0,6))
  ax.annotate('Time (s): {0:.2f}'.format(puls.Time), xy=(0,5))
  ax.annotate('Sigma: {0:.1f}'.format(puls.Sigma), xy=(0,4))
  ax.annotate('Width (ms): {0:.0f}'.format(puls.Duration*1000), xy=(0,3))
  ax.annotate('N events: {0:.0f}'.format(puls.N_events), xy=(0,2))

  ax.axis('off')
  return


def puls_DM_Time(ax, event, all_events, puls):
  #--------------------------
  # Plot events in DM vs Time
  #--------------------------

  # Define circle sizes
  def circle_size(values):
    new_val = np.clip(values,6.5,20)
    m = 31.85
    q = -137.025
    return (new_val * m + q) / 2.

  # Plot eventsin the pulse
  sig = circle_size(event.Sigma)
  x_ev = (event.Time - puls.Time)
  ax.scatter(x_ev, event.DM, s=sig, marker='o', linewidths=[.5,], edgecolor='r', facecolors='none', zorder=0)
  
  # Select events in the window and plot them
  DM_min = puls.DM - puls.dDM * 3.
  DM_max = puls.DM + puls.dDM * 3.
  Time_min = puls.Time - puls.dt * 3.
  Time_max = puls.Time + puls.dt * 3.
  events = all_events[(all_events.DM > DM_min) & (all_events.DM < DM_max) & (all_events.Time > Time_min) & (all_events.Time < Time_max)]
  x_ev = (events.Time - puls.Time)
  ax.scatter(x_ev, events.DM, s=15., marker='o', c='k', linewidths=[0.,], zorder=1)

  # Plot pulse position
  ax.scatter(0, puls.DM, marker='x', s=30., color='dodgerblue', zorder=2)

  ax.set_xlabel('$\Delta$Time (s)')
  ax.set_ylabel('DM (pc cm$^{-3}$)')
  return


def puls_dedispersed(ax, prof_ax, puls, args):
  #----------------------------------------------------
  # Plot timeseries in a DM range and the pulse profile
  #----------------------------------------------------

  def load_ts(puls, filename=False):
    #---------------------------------------
    # Load timeseries from PRESTO .dat files
    #---------------------------------------

    # Select timeseries to load
    bin_peak = int(puls['Sample'])
    DM_peak = puls['DM']
    ts_list = sorted(glob(filename))
    DM_list = np.array([float(n[n.index('_DM') + 3 : a.index('_DM') + 8]) for n in ts_list])
    if np.diff(DM_list).unique().size > 1: raise ValueError('Spacing of timeseries is not uniform.')
    DM_range = np.where((DM_list > puls.DM - 2 * puls.dDM) & (DM_list < puls.DM + 2 * puls.dDM))[0]
    ts_list[DM_range[0] : DM_range[-1]]

    # Select number of pixels in the plot
    nProfBins = 3
    scrunch_fact = int(np.round(puls.Downfact / float(nProfBins)))
    if scrunch_fact < 1: scrunch_fact = 1
    nPlotBins = int(puls.Downsample * 10.)
    nBins = nPlotBins * nProfBins
    data = np.zeros((nDMs,nBins))
    bin_start = bin_peak - nBins/2 * scrunch_fact

    for j, ts in enumerate(ts_list):
      ts = np.memmap(ts, dtype=np.float32, mode='r', offset=bin_start*4, shape=(nBins*scrunch_fact,))
      ts = np.mean(np.reshape(ts, (nBins, scrunch_fact)), axis=1)
      data[j] = ts

    params = {'duration': nBins*scrunch_fact, 'DM_min': DM_range[0], 'DM_max': DM_range[-1], 'k': k}
  
    return data, params

  data, params = load_ts(puls, filename=args.timeseries)

  # Image plot
  extent=[-params['duration']/2.,params['duration']/2.,params['DM_min'], params['DM_max']
  ax.imshow(data,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent])
  ax.set_ylim((params['DM_min'], params['DM_max']))
  ax.set_ylabel('DM (pc cm$^{-3}$)')
  ax.set_xlim((-params['duration']/2.,params['duration']/2.))
  ax.set_xlabel('$\Delta$Time (s)')

  # Profile plot
  nPlotBins = 20
  nProfBins = 5
  nBins = nPlotBins * nProfBins
  ts = data[data.shape[0]/2+1]
  prof_ax.plot(ts, 'k')
  prof_ax.set_xlim((0,ts.size))
  prof_ax.set_xticks([])
  prof_ax.set_yticks([])

  return


def puls_SNR_DM(ax, event):
  #----------------------------------------------------------
  # Plot S/N and duration vs DM of the events forming a pulse
  #----------------------------------------------------------

  # Plot S/N vs DM
  ax.plot(event.DM, event.Sigma, 'k')
  ax.set_xlabel('DM (pc cm$^{-3}$)')
  ax.set_ylabel('S/N')
  
  # Plot duration vs DM
  ax2 = ax.twinx()
  ax2.plot(event.DM, event.Duration*1000, 'r')
  ax2.set_ylabel('Width (ms)', color='r')
  for tl in ax2.get_yticklabels():
    tl.set_color('r')

  ax.set_zorder(ax2.get_zorder()+1)
  ax.patch.set_visible(False)
  return


def puls_dynSpec(ax1, ax2, puls, args, prof_ax=None):
  #-----------------------------------------------------------
  # Plot dynamic spectra (dispersed and not) from a .fits file
  #-----------------------------------------------------------

  if maskfn is not None: mask = True  
  duration = 20
  df = int(puls.Downfact)
  psrfits_file = psrfits.PsrfitsFile(args.fits)

  # Plot de-dispersed dynamic spectrum
  start = puls.Time_org - puls.Duration * duration / 2.
  duration = puls.Duration * (duration + 1)
  ds, nbinsextra, nbins, start = waterfaller.waterfall(psrfits_file, start, duration, nsub=16, dm=puls.DM, width_bins=df, /
      maskfn=args.mask, mask=mask, scaleindep=False, bandpass_corr=True, zerodm=True)
  puls_t = - puls.Duration * duration / 2.
  waterfaller.plot_waterfall(ds, start, duration, ax_im=ax1, interactive=False, puls_t=puls_t)
  ax1.scatter(0., ax1.get_ylim[0], marker='^', s=50, c='r', lw=0.)

  # Plot pulse profile
  if prof_ax is not None:
    prof = ds.mean(axis=1)
    prof_ax.plot(prof, 'k')
    prof_ax.set_xlim((0, prof.size))
    prof_ax.set_xticks([])
    prof_ax.set_yticks([])

  # Plot dispersed dynamic spectrum
  DM_delay = psr_utils.delay_from_DM(puls.DM, F_MIN) - psr_utils.delay_from_DM(puls.DM, F_MAX)
  start = puls.Time_org - 0.1
  duration = puls.Duration + DM_delay + 0.2
  ds, nbinsextra, nbins, start = waterfaller.waterfall(psrfits_file, start, duration, dm=0., nsub=32*3, downsamp=df*3, /
      maskfn=maskfn, mask=mask, scaleindep=True, zerodm=True, bandpass_corr=True)
  waterfaller.plot_waterfall(ds, start, duration, ax_im=ax2, interactive=False, sweep_dms=[puls.DM], puls_t=-0.1)

  return



