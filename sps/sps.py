import time
import argparse
import subprocess
import os

from src import SPclean

import pandas as pd
pd.options.mode.chained_assignment = None

import warnings
#warnings.filterwarnings("error")
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


def parser():
  '''
  Command-line options
  '''
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description="Search for interesting single pulses in a dataset. See the README for detailed instructions.")
  parser.add_argument('filename', help="Name of the .singlepulse file or HDF5 database.")
  parser.add_argument('-meta_data', help="Name of the .inf file from PRESTO.", default=None)

  parser.add_argument('-no_search', help="Do not search for pulses and candidates.", action='store_true')
  parser.add_argument('-no_store', help="Do not store the output into an HDF5 database.", action='store_true')
  parser.add_argument('-store_name', help="Filename of the output HDF5 database.", default='SinglePulses.hdf5')
  parser.add_argument('-no_plot', help="Do not plot the diagnostics.", action='store_true')
  parser.add_argument('-plot_name', help="Name of the diagnostic plot file.", default='diagnostics.pdf')
  parser.add_argument('-no_filter', help="Do not filter pulses.", action='store_true')
  
  parser.add_argument('-F_range', help="Frequency range of the observation.", type=float, nargs=2)
  parser.add_argument('-t_range', help="Time range to search (seconds).", type=float, nargs=2)
  parser.add_argument('-DM_range', help="DM range to search (pc/cc).", type=float, nargs=2)
  parser.add_argument('-N_min', help="Minimum number of events forming the pulse.", default=5, type=int)
  parser.add_argument('-SNR_min', help="Minimum events S/N to search.", default=5, type=float)
  parser.add_argument('-SNR_peak_min', help="Minimum pulse S/N to search.", default=6.5, type=float)
  parser.add_argument('-events_dt', help="Time in seconds two consecutive events are considered part of the same pulse.", default=3e-2, type=float)
  parser.add_argument('-events_dDM', help="Number of DM steps two consecutive events are considered part of the same pulse.", default=10, type=int)
  
  parser.add_argument('-DM_cand', help="DM span to consider two pulses the same candidate.", default=0.1, type=float)
  parser.add_argument('-single_cand_SNR', help="Minimum S/N of single candidates to consider.", default=6.5, type=float)
  parser.add_argument('-multiple_cand_SNR', help="Minimum S/N of repeated candidates to consider.", default=6.5, type=float)
  
  #parser.add_argument('-timeseries', help="Glob expression of PRESTO .dat file for plotting.", default=None)
  #parser.add_argument('-fits', help="Name of a fits file for plotting.", default=None)
  #parser.add_argument('-mask', help="Name of a PRESTO .mask file for plotting.", default=None)
  args = parser.parse_args()
  return args


def main():
  print "The DataBase is being created"

  args = parser()
  if args.F_range is not None: args.F_range = sorted(args.F_range)
  if args.t_range is not None: args.t_range = sorted(args.t_range)
  if args.DM_range is not None: args.DM_range = sorted(args.DM_range)
  
  scriptFolder = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
  git_folder = os.path.join(scriptFolder, '.git')
  vers = subprocess.check_output(['git','--git-dir',git_folder,'describe','--tags','--abbrev=0','--always']).strip()
  print "L-SpS version used: {}".format(vers)
  
  time0 = time.time()

  SPclean.main(args, vers=vers)

  print "The DataBase has been created"
  print "Time spent: {:.2f} s".format(time.time() - time0)

  return


if __name__ == '__main__':
  main()
