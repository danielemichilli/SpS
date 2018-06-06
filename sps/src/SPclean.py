import os
import fnmatch
import multiprocessing as mp
import re
from glob import glob

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import Events
import Pulses
import RFIexcision
import LSPplot
import Candidates



def main(args, vers=None):
  
  # Load Events from HDF5 database of .singlepulse file
  try: 
    events = pd.read_hdf(args.filename, 'events')
    file_type = 'hdf5'
  except IOError: 
    events = Events.Loader(args)
    file_type = 'sp'
  if events.empty: 
    print "No events found. Exiting"
    return

  # Load meta data
  if args.meta_data is not None: meta_data = Events.meta_data_Loader(args.meta_data)
  if vers is not None: meta_data['version'] = vers
  meta_data['File'] = os.path.basename(args.filename)
  if args.no_store is not None: meta_data.to_hdf(args.store_name, 'meta_data')

  # If no_search is True, load from HDF5 file if exists
  if not args.no_search:
    # Load Pulses
    pulses = Pulses.Loader(events, args)
    pulses = pulses[pulses.Rank == 0]

    # Select pulses within the defined range
    if args.t_min is not None: pulses = pulses[pulses.Time >= args.t_min]
    if args.t_max is not None: pulses = pulses[pulses.Time <= args.t_max]
    if args.DM_min is not None: pulses = pulses[pulses.DM >= args.DM_min]
    if args.DM_max is not None: pulses = pulses[pulses.DM <= args.DM_min]
    if args.SNR_mim is not None: pulses = pulses[pulses.Sigma >= args.SNR_min]

    if pulses.empty: 
      print "No pulses found. Exiting"
      return

    # Load Candidates
    cands = Candidates.Loader(pulses, args)
    if cands.empty: 
      print "No candidates found. Exiting"
      return

  elif file_type == 'hdf5':
    # Load Pulses form database
    pulses = pd.read_hdf(args.filename, 'pulses')
    pulses.sort(['DM','Time'], inplace=True)

    # Load Candidates form database
    cands = pd.read_hdf(args.filename, 'candidates')
    cands.sort('Sigma', inplace=True, ascending=False)
    cands = cands[cands.main_cand == 0]

  else:
    print "Events stored into the HDF5 file. Exiting"
    return

  if cands.shape[0] > 100:
    print "{} candidates found, only the brightest 100 will be processed.".format(cands.shape[0])
    cands = cands.head(100)

  cands = cands[ ((cands.N_pulses == 1) & (cands.Sigma >= args.single_cand_SNR)) | ((cands.N_pulses > 1) & (cands.Sigma >= args.multiple_cand_SNR)) ]
  cands.sort('Sigma', inplace=True, ascending=False)

  #Produce the output
  if not args.no_plot:
    LSPplot.output(args, events, pulses, cands, meta_data)    

  return


