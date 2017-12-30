import os
import multiprocessing as mp
import numpy as np

import Internet
import LSPplot
from Paths import *

def output(idL,pulses,meta_data,candidates):
  folder = '{}/sp/'.format(TEMP_FOLDER.format(idL))
  pulses.sort('Sigma',ascending=False,inplace=True)
  
  if not candidates.empty:
    candidates.sort(['Rank','Sigma'],ascending=False,inplace=True)
  
    #Repeated candidates
    if candidates[candidates.N_pulses>1].shape[0] > 0:
      LSPplot.repeated_candidates(pulses,candidates[candidates.N_pulses>1],meta_data,folder,idL)
  
    #Single candidates
    if candidates[candidates.N_pulses==1].shape[0] > 0:
      LSPplot.single_candidates(pulses,candidates[candidates.N_pulses==1],meta_data,folder,idL)

  beams_parallel(pulses,meta_data,folder,idL)

  pulses = pulses[pulses.Pulse==0]
  output_pointing(pulses,folder)
  
  return


def beams_parallel(pulses,meta_data,folder,idL):
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  dirs = [n for n in gb_puls.indices.iterkeys()]
  
  CPUs = mp.cpu_count()
  dirs_range = int(np.ceil(len(dirs)/float(CPUs)))

  pool = mp.Pool()
  dirs_CPU = [dirs[i*dirs_range:(i+1)*dirs_range] for i in range(CPUs) if len(dirs[i*dirs_range:(i+1)*dirs_range]) > 0]
  [pool.apply_async(output_beams, args=(pulses[pulses.SAP.isin(zip(*dir_CPU)[0])&pulses.BEAM.isin(zip(*dir_CPU)[1])],meta_data,folder,idL,dir_CPU)) for dir_CPU in dirs_CPU]
  pool.close()
  pool.join()
  return    


def output_beams(pulses,meta_data,folder,idL,dirs):
  for (sap,beam) in dirs:
    store = '{}/diagnostics'.format(folder)
    name = 'SAP{}_BEAM{}'.format(sap,beam)
    os.makedirs('{}/{}'.format(store,name))
    
    top = pulses[(pulses.SAP==sap) & (pulses.BEAM==beam) & (pulses.Pulse==0)]
    good = pulses[(pulses.SAP==sap) & (pulses.BEAM==beam) & (pulses.Pulse==1)]
    
    if beam == 12:
      if top.shape[0] > 0: 
        LSPplot.sp_shape(top.head(10),'{}/SAP{}_BEAM{}/top_candidates(0-9).png'.format(store,sap,beam),folder,idL)
        if top.shape[0] > 10:
          LSPplot.sp_shape(top.iloc[10:20],'{}/SAP{}_BEAM{}/top_candidates(10-19).png'.format(store,sap,beam),folder,idL)
          if top.shape[0] > 20: 
            LSPplot.sp_shape(top.iloc[20:30],'{}/SAP{}_BEAM{}/top_candidates(20-29).png'.format(store,sap,beam),folder,idL)
      LSPplot.sp_plot(top,good,meta_data,sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(store,sap,beam))
      
    else:
      if top.shape[0] > 0: 
        LSPplot.sp_shape(top.head(10),'{}/SAP{}_BEAM{}/top_candidates.png'.format(store,sap,beam),folder,idL)
        LSPplot.sp_plot(top,good,meta_data,sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(store,sap,beam))
    
  return
  
  

def output_pointing(pulses,folder):
  
  top_candidates = pulses[pulses.BEAM>12].groupby(['SAP','BEAM'],sort=False).head(10)
  top_candidates = top_candidates.append(pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30),ignore_index=False)
  top_candidates.sort(['SAP','BEAM','Sigma'],ascending=[True,True,False],inplace=True)
  top_candidates['code'] = top_candidates.index
  if not top_candidates.empty:
    a = top_candidates.groupby(['SAP','BEAM'],sort=False).apply(lambda x: range(len(x))).tolist()
    b = [val for sublist in a for val in sublist]
    top_candidates.index = b
  top_candidates.Duration *= 1000
  top_candidates['void'] = ''
  top_candidates.to_csv('{}/diagnostics/top_candidates.inf'.format(folder),sep='\t',float_format='%.2f',\
    columns=['SAP','BEAM','Sigma','DM','void','Time','void','Duration','void','code'],\
    header=['SAP','BEAM','Sigma','DM (pc/cm3)','Time (s)','Duration (ms)','code','','',''],index_label='rank',encoding='utf-8')
  
  puls = pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30)
  if not puls.empty: LSPplot.obs_top_candidates(puls,'{}/diagnostics/inc_top_candidates.png'.format(folder),incoherent=True)
  
  for sap in pulses.SAP.unique():
    puls = pulses[(pulses.SAP==sap)&(pulses.BEAM>12)].groupby('BEAM',sort=False).head(10)
    if not puls.empty: LSPplot.obs_top_candidates(puls,'{}/diagnostics/top_candidates(SAP{}).png'.format(folder,sap))

  return


 
