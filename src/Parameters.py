#############################
#
# Parameters
#
# Written by Daniele Michilli
#
#############################


DM_MIN = 3  #pc/cm3
SNR_MIN = 6.5
#DURATION_MAX = 0.049152  #s


F_MIN = 119.43  #MHz
F_MAX = F_MIN + 31.64  #MHz
RES = 491.52e-6  #s


#Grouping
ERR_FLOAT = 0.001  #Error on float numbers
STEPS_GROUP = 20  #Tollerance on the number of DM steps without a candidate
DURAT_GROUP = 0.03  #Tollerance in seconds on the temporal shift of two consecutive events

SIGMA_TOLL = 4
SIGMA_TOLL_IB = 2

RFI_percent = 2
PULS_LENGHT = 15

DS_OFFSET = 100000 #bins


FILTERS = {
  'scattered': 0.09,
  'aligned': 0.045,
  'peak_central': 0.6,
  'duration_central': 1,
  'holes': 0.75,
  'variance': 0.06,
  'flat_duration': 0.9,
  'm': 23.684210526315788,
  'q': -60.86842105263156,
  'flat_SNR': 0.78,
  'flat_SNR_simmetric': 0.79,
  'flat_SNR_extremes': 0.68,
  'DM_extremes': 0.9,
  'sigma_min': 0.95,
  'cum_scatter': 0.0045,
  'std_scatter': 0.0068,
  'sigma_std': 0.46,
  'sigma_scatter': 0.02,
  'sigma_scatter_max': 0.47,
  'sigma_std_largest': 0.55,
  'flat_fit0': 0.72,
  'flat_fit1': 0.64,
  'sigma_std_largest_weak': 0.42,
  'flat_fit0_weak': 0.33,
  'flat_fit1_weak': 0.21,
  'bright_extremes_abs': 0.9,
  'bright_extremes_rel': 0.89,
  'pulse_simmetric': 0.155,
  'number_events': 1.,
  'monotonic': 0.0000001,
  'sigma_jumps': 0.3,
  'fit1_brightest': 0.09}

