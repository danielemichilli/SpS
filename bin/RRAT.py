import argparse
import os
import pandas as pd
import numpy as np

folder_out = '/home/danielem/ALERTS/'

parser = argparse.ArgumentParser()
parser.add_argument('folder', nargs=1)
args = parser.parse_args()
folder = args.folder[0]
idL = os.path.basename(os.path.dirname(folder))

#if os.path.isfile('{}/ALERTS'.format(folder)):
  #print 'Observation {} already analyzed!\nExiting...'.format(idL)
  #exit()

pulses = pd.read_hdf('{}/SinglePulses.hdf5'.format(folder),idL+'_pulses',where=['(Pulse<=2)&(Pulse>=0)'])
pulses = pulses.loc[:,['SAP','BEAM','DM','Sigma','Pulse']]
pulses = pulses.loc[pulses.BEAM>12]
#pulses = pulses.loc[pulses.DM>5]
pulses.DM = pulses.DM.astype(np.float64)

file = 0

pulses.DM = pulses.DM.round(decimals=1)
noise_level = pulses.groupby(['SAP','BEAM','DM'],sort=False).Sigma.sum().astype(np.float64)
noise_level = noise_level.median(level=['SAP','DM'])

beams = pulses.groupby(['SAP','BEAM'],sort=False)
for ind,beam in beams:
  beam_output = 0
  val = beam.groupby('DM').Sigma.sum()
  lim = noise_level.loc[ind[0]]
  lim = lim[lim.index.isin(val.index)]
  ratio = val/lim
  ratio = ratio.loc[5.5:]
  ratio = ratio[ratio>7]
  if not ratio.empty:
    header = 0
    for i in ratio.index:
      counts = beam.groupby('DM',sort=False).Sigma.count().loc[i]
      val_max = beam.groupby('DM',sort=False).Sigma.max().loc[i]
      if (val_max > 7) | (counts > 10):
        if not file:
          file = open('{}{}.txt'.format(folder_out,idL),'w')
          file.write('Obs. {}\n'.format(idL))
        if beam_output == 0:
          file.write('\n\tSAP{}_BEAM{}\n'.format(ind[0],ind[1]))
          beam_output = 1
        if header == 0:
          file.write('\n\t\tStrong pulses compared to other beams:\n\n')
          file.write('\t\t\tDM\tRatio\tCounts\tSNR_max\n')
          header = 1
        file.write('\t\t\t{0:.2f}\t{1:.1f}\t{2}\t{3:.1f}\n'.format(i,ratio[i],counts,val_max))
  
  val = pd.DataFrame(val)
  val['med'] = val.index.round()
  val.med = val.groupby('med',sort=False).Sigma.median()
  val.med.fillna(method='ffill',inplace=True)
  val = val.loc[5.5:]
  #ratio = np.abs((val+5)/val.shift(10)).dropna()
  ratio = val.Sigma/val.med
  ratio = ratio[ratio>7]
  if not ratio.empty:
    header = 0
    for i in ratio.index:
      counts = beam.groupby('DM',sort=False).Sigma.count().loc[i]
      val_max = beam.groupby('DM',sort=False).Sigma.max().loc[i]
      if (val_max > 7) | (counts > 10):
        if not file:
          file = open('{}{}.txt'.format(folder_out,idL),'w')
          file.write('Obs. {}\n'.format(idL))
        if beam_output == 0:
          file.write('\n\tSAP{}_BEAM{}\n'.format(ind[0],ind[1]))
          beam_output = 1
        if header == 0:
          file.write('\n\t\tStrong pulses compared to other DMs:\n\n')
          file.write('\t\t\tDM\tRatio\tCounts\tSNR_max\n')
          header = 1
        file.write('\t\t\t{0:.2f}\t{1:.1f}\t{2}\t{3:.1f}\n'.format(i,ratio[i],counts,val_max))
  

if file:
  file.write('\n\n\n')
  file.close()

