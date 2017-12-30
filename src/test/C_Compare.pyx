cimport cython
from Parameters import *

@cython.boundscheck(False)
@cython.wraparound(False)


#-------------------------
# Compares different beams    #mettere duration invece di dT
#-------------------------
def Compare(long[::1] code1 not None,
            float[::1] DM_c_l not None,
            float[::1] dDM_l not None,
            float[::1] Time_c_l not None,
            float[::1] dTime_l not None,
            float[::1] Sigma_l not None,
            signed char[::1] Pulse_l not None,
            long[::1] code2 not None,
            float[::1] DM_c_r not None,
            float[::1] dDM_r not None,
            float[::1] Time_c_r not None,
            float[::1] dTime_r not None,
            float[::1] Sigma_r not None,
            signed char[::1] Pulse_r not None,
            int CB):

  cdef:
    unsigned int i, j, RFI
    unsigned int dim_l = len(DM_c_l)
    unsigned int dim_r = len(DM_c_r)
    unsigned int j_min = 0
    int TollSigma
    int rfi_limit = RFI_percent
    float DTime, Time, DM, DDM, sign
    float Duration_Max = FILTERS[0][2]
  
  # Uses different tollerances on sigma for coherent and incoherent beams
  if CB==int(1): TollSigma = SIGMA_TOLL
  else: TollSigma = SIGMA_TOLL_IB

  # Compare each pulse of the first group with each of the second
  # Assign different codes under certain conditions 
  for i in range(0,dim_l):
    
    RFI = 0
    if Pulse_l[i] >= rfi_limit: RFI = 1 
    
    j_flag = j_min
    
    for j in range(j_min, dim_r):
      
      if Pulse_r[j] >= rfi_limit: 
        
        if RFI > 0:
        
          continue
      
      Time = abs(Time_c_l[i]-Time_c_r[j])
      DTime = dTime_l[i]+dTime_r[j]
      
      sign = Time_c_l[i]-Time_c_r[j]
      
      if Time < 4.*DTime :
      
        if j_flag == j_min: j_min = j
        
        DM = abs(DM_c_l[i]-DM_c_r[j])
        DDM = dDM_l[i]+dDM_r[j]
        
        if DM < 2.*DDM :
          
          if CB==int(1):
        
            if abs(Sigma_l[i]-Sigma_r[j]) < TollSigma:
              
              Pulse_l[i] += 1
              Pulse_r[j] += 1
          
          else:
          
            if Sigma_l[i]/Sigma_r[j] < TollSigma:
            
              Pulse_l[i] += 1
              #Pulse_r[j] += 1  #non utilizzabile (se non con parametri molto ristretti): puls potrebbe essere dentro IB ma fuori CB e produrre Sigma simili
              
              
      elif sign < 0.: 
        
        break
      
  return