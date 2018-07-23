import os
import fnmatch
import multiprocessing as mp
import re
from glob import glob

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from tables.exceptions import HDF5ExtError

import Events
import Pulses
import LSPplot
import Candidates



def main(args, vers=None):
  
  # Load Events from HDF5 database of .singlepulse file
  try:
    events = pd.read_hdf(args.filename, 'events')
    with pd.HDFStore(args.filename) as store:
      db_keys = store.keys()
    file_type = 'hdf5'
  except (IOError, HDF5ExtError) as e:
    events = Events.Loader(args)
    file_type = 'sp'
 
  # Select events within the defined ranges
  if args.SNR_min is not None: events = events[events.Sigma >= args.SNR_min]
    
  if events.empty: 
    print "No events found. Exiting"
    return

  # Load meta data
  if args.meta_data is not None:
    meta_data = Events.meta_data_Loader(args.meta_data)
    if vers is not None: meta_data['version'] = vers
    meta_data['File'] = os.path.basename(args.filename)
    if args.no_store is not None: meta_data.to_hdf(args.store_name, 'meta_data')
  elif file_type == 'hdf5':
    if '/meta_data' in db_keys:
      meta_data = pd.read_hdf(args.filename, 'meta_data')
  else: meta_data = None

  # Load Pulses
  if not args.no_search: pulses = Pulses.Loader(events, args)
  elif file_type == 'hdf5':
   if '/pulses' in db_keys:
     pulses = pd.read_hdf(args.filename, 'pulses')
   else:
     print "Pulses not present in the database. Exiting"
     return
  else:
    print "Events have been loaded and stored into the HDF5 file. Exiting"
    return
    
  if not args.no_filter: pulses = pulses[pulses.Rank == 0]
 
  # Select pulses within the defined ranges
  if args.t_range is not None: pulses = pulses[(pulses.Time >= args.t_range[0]) &(pulses.Time <= args.t_range[1])]
  if args.DM_range is not None: pulses = pulses[(pulses.DM >= args.DM_range[0]) & (pulses.DM <= args.DM_range[1])]
  if args.SNR_min is not None: pulses = pulses[pulses.Sigma >= args.SNR_peak_min]
  if args.N_min is not None: pulses = pulses[pulses.N_events >= args.N_min]

  if pulses.empty: 
    print "No pulses found. Exiting"
    return

  # Load Candidates
  cands = Candidates.Loader(pulses, args)
  if not args.no_search or not args.no_store: cands.to_hdf(args.store_name, 'candidates')

  cands = cands[cands.main_cand == 0]
      
  if cands.empty: 
    print "No candidates found. Exiting"
    return

  if cands.shape[0] > 100:
    print "{} candidates found, only the brightest 100 will be processed.".format(cands.shape[0])
    cands = cands.head(100)

  cands = cands[ ((cands.N_pulses == 1) & (cands.Sigma >= args.single_cand_SNR)) | ((cands.N_pulses > 1) & (cands.Sigma >= args.multiple_cand_SNR)) ]
  cands.sort_values('Sigma', inplace=True, ascending=False)

  #Produce the output
  if not args.no_plot:
    LSPplot.output(args, events, pulses, cands, meta_data)    

  return


