#!/usr/bin/env python

import pandas as pd
import argparse

import Internet


def upload_website():
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description="The program upload existing candidates from a specific observation an the website.")
  parser.add_argument('idL', nargs=1, help='Observation ID')
  parser.add_argument('-folder', nargs=1, help='Path of the folder containig the observation')
  args = parser.parse_args()
  
  if args.folder: folder = args.folder[0] + '/'
  else: folder = Paths.OBS_FOLDER + '/'
  idL = args.idL[0]
  
  cands = pd.read_hdf('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'candidates')
  cands = cands[cands.main_cand==0].head(30)
  
  Internet.upload(cands,folder,idL)


if __name__ == '__main__':
  try: upload_website()
  except: print "ATTENTION! Website currently down. Try to upload the observation later"