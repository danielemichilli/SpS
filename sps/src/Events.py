import os
from glob import glob

import pandas as pd
import numpy as np

import C_Funct
from Parameters import *


def Loader(args):
  #--------------------------
  # Creates a table of events
  #--------------------------
  
  events = pd.DataFrame()
  sp_file = args.filename
  
  #Load the events and meta-data tables
  # events = pd.concat(pd.read_csv(f, delim_whitespace=True, dtype=np.float64) for f in sp_file if os.stat(f).st_size > 0)  # Parse a glob expression
  events = pd.read_csv(sp_file, delim_whitespace=True, dtype=np.float64, comment='#')
  events = events.ix[:,:5]
  events.columns = ['DM','Sigma','Time','Sample','Downfact']
  events.index.name = 'idx'

  events['Duration'] = events.Sampling * events.Downfact
  events.Duration = events.Duration.astype(np.float32)
  events.Sample = events.Sample.astype(np.int32)
  
  events['Pulse'] = 0
  events.Pulse = events.Pulse.astype(np.int64)

  #Correct for the time misalignment of events
  events['Time_org'] = events.Time.copy()
  events.sort(['DM','Time'],inplace=True)  #Needed by TimeAlign
  events.Time = TimeAlign(events.Time.copy(), events.DM, args)

  #Group the events
  events.sort(['DM','Time'],inplace=True) #Needed by Group
  C_Funct.Get_Group(events.DM.values, events.Sigma.values, events.Time.values, events.Duration.values, events.Pulse.values)

  #Store the events
  events.sort(['DM','Time'], inplace=True)
  if not args.no_store: events.to_hdf(args.store_name, 'events')
  events = events[events.Pulse > 0]

  return events


def TimeAlign(Time, DM, args):
  #-------------------------------------------------
  # Corrects for the time misalignment of the pulses
  #-------------------------------------------------
  
  # Quantifies the misalignment for a broad-band pulse
  # Only the extreme frequencies are taken into account
  k = 4148.808 #s-1
  delay = k * (args.F_MIN**-2 - args.F_MAX**-2)

  Time += np.float32( delay * DM / 2)
  
  return Time


def meta_data_Loader(inf_file):
  #------------------------------
  # Create a table with meta data
  #------------------------------

  inf = pd.read_csv(inf_file, sep="=", dtype=str, error_bad_lines=False, warn_bad_lines=False, header=None, skipinitialspace=True)
  inf = inf.iloc[[0,1,2,4,5,7],1]  # The selection of lines is now index-based, this should be substituted with label-based
  inf.iloc[0] = inf.iloc[0].replace("_rfifind","")
  inf = pd.DataFrame(inf).T
  inf.columns = ['File','Telescope','Instrument','RA','DEC','Epoch']
  inf = inf.astype(str)
  return inf


