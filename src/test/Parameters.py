#############################
#
# Parameters
#
# Written by Daniele Michilli
#
#############################


DM_MIN = 5.  #pc/cm3
SIGMA_MIN = 10.
#DURATION_MAX = 0.049152  #s


F_MIN = 119.43  #MHz
F_MAX = F_MIN + 31.64  #MHz

ERR_FLOAT = 0.001
STEPS_GROUP = 3
DURAT_GROUP = 1

SIGMA_TOLL = 2
SIGMA_TOLL_IB = 1.5

RFI_percent = 3

FILTERS = [
  [0.0344,0.0344,0.0344],  #duration
  [0.9,0.773,0.675],  #dDM/(N_events-1)/step
  [1,1,1],  #abs(puls.DM-puls.DM_c)/puls.dDM
  [1.03,1.05,1.08],  #puls.Sigma/Sigma_min
  [0.00633,0.00752,0.000341],  #puls.Sigma/Sigma_min**4
  [1.45,1.3,10.2],  #abs(Sigma_DM_max-Sigma_DM_min)
  [3.04,9.06,2.83],  #f (>)
  [-41,-36.3,-48.8],  #f (<)
  [0.0000361,0.0000824,0.00128]]  #gb.Time.apply(np.var)


FILTERS_BEST = [
  [0.667,0.636,0.537],  #dDM/(N_events-1)/step
  [0.667,0.636,0.714],  #abs(puls.DM-puls.DM_c)/puls.dDM
  [1.06,1.09,1.19],  #puls.Sigma/Sigma_min
  [0.00806,0.00836,0.00926],  #puls.Sigma/Sigma_min**4
  [0.54,0.41,3.39],  #abs(Sigma_DM_max-Sigma_DM_min)
  [-12.2,-14.2,-14.6],  #f (>)
  [-25.7,-26.1,-27.1],  #f (<)
  [0.00000309,0.00000728,0.0000634]]  #gb.Time.apply(np.var)
