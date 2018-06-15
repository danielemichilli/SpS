import pandas as pd
import numpy as np

import Pulses
import C_Funct


def Loader(pulses, args):
  #------------------------------
  # Creates a table of candidates
  #------------------------------
  
  pulses.sort_values(['Sigma','Rank'], ascending=[0,1], inplace=True)
  
  Repeated_candidates(pulses, args.DM_cand)
  
  # Identify bright pulses
  cands_unique = pulses[pulses.Candidate==-1].astype(np.int32)
  pulses.Candidate.loc[cands_unique.index.get_level_values('idx')] = 2 * np.arange(cands_unique.shape[0])  #Unique candidates have even ID
  
  if not pulses[pulses.Candidate>=0].empty:
    cands = candidates_generator(pulses[pulses.Candidate>=0].copy())
    cands['main_cand'] = 0
  
    #Unify the same repeated candidates in different beams
    cands.sort_values('Sigma',ascending=True,inplace=True)

    C_Funct.Compare_candidates(cands.DM.astype(np.float32).values,cands.Time.astype(np.float32).values,cands.index.values,cands.main_cand.values)
    
    cands.sort_values(['main_cand', 'Sigma'], ascending=[1,0], inplace=True)

  else: cands = pd.DataFrame()
  
  return cands


def Repeated_candidates(pulses, span):
  #------------------------------------------
  # Change candidate code for repeated pulses
  #------------------------------------------

  pulses.Candidate[:] = -1
  puls_beam = pulses.copy()
  
  # Find pulses within span DM
  def group_SNR(DM, pulses): return pulses.Sigma[(pulses.DM >= DM - span) & (pulses.DM <= DM + span)].sum()
  def group_count(DM, pulses): return pulses[(pulses.DM >= DM - span) & (pulses.DM <= DM + span)].shape[0]
  puls_beam['top_SNR'] = puls_beam.apply(lambda x: group_SNR(x.DM, puls_beam), axis=1)
  puls_beam['top_count'] = puls_beam.apply(lambda x: group_count(x.DM, puls_beam), axis=1)
  puls_beam = puls_beam[puls_beam.top_count > 1]
  
  # Assign the same candidate code to repeated pulses
  i = 1
  while puls_beam.shape[0] > 0:
    DM = puls_beam.DM[puls_beam.top_SNR.idxmax()]
    selected_pulses = puls_beam.Candidate[(puls_beam.DM >= DM - span) & (puls_beam.DM <= DM + span)]
    if selected_pulses.shape[0] > 1:
      pulses.Candidate.loc[selected_pulses.index] = 1 + 2 * i #Repeated candidates have odd ID
    puls_beam = puls_beam.drop(selected_pulses.index)
    i += 1
    
  return


def period(x):
  if x.size<=1: return 0
  else: return rrat_period(x)


def rrat_period(times, numperiods=20000):
  # Modified version of PRESTO utility
  ts = np.asarray(sorted(times))
  ps = (ts[1]-ts[0])/np.arange(1, numperiods+1)
  dts = np.diff(ts)
  xs = dts / ps[:,np.newaxis]
  metric = np.sum(np.fabs((xs - xs.round())), axis=1)
  pnum = metric.argmin()
  numrots = xs.round()[pnum].sum()
  p = (ts[-1] - ts[0]) / numrots
  return p


def candidates_generator(pulses):
  #-------------------------
  # Create a candidate table
  #-------------------------

  pulses['Period'] = pulses.Time
  pulses['Period_err'] = pulses.Time
  
  cands = pulses.groupby('Candidate',as_index=False,sort=False).agg({'Sigma':np.sum,'N_events':np.size,'DM':np.mean,'Time':np.min,'Period':period,'Rank':np.max})
  
  cands = cands.astype(np.float32)
  cands[['N_events','Rank']] = cands[['N_events','Rank']].astype(np.int16)
  
  cands.index = cands.Candidate.astype(int)
  cands.index.name = 'idx'
  cands.rename(columns={'N_events': 'N_pulses'}, inplace=True)
  cands = cands.drop('Candidate',axis=1)
  cands.Time[cands.N_pulses>1] = 0

  return cands
