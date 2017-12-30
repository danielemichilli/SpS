#############################
#
# Cython Functions
#
# Functions written in cython
# are grouped here.
# 
# Written by Daniele Michilli
#
#############################

#Installation: sudo python setup.py build_ext --inplace

cimport cython
from Parameters import *

#---------------------------------
# Gives a pulse-code to each event
#---------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def Get_Group(float[::1] DM not None,
          float[::1] Sigma not None,
          float[::1] Time not None,
          float[::1] Duration not None,
          int[::1] Pulse not None):
  
  cdef:
    unsigned int i, j, k, j_min, j_max, empty, SNR_max
    unsigned int n_steps = STEPS_GROUP
    float durat = DURAT_GROUP
    unsigned int code = 0
    unsigned int dim = len(DM)
    float step, step_min, step_max, dDM, DM_min
    float DM_new = -1.
    float float_err = 0.001


  # Assign a code to each event.
  # Events must have been ordered in DM and then in Time
  for i in range(0,dim):
    
    #Remove close events at the same DM
    j = i+1
    if DM[i] == DM[j]:
    
      if abs(Time[i]-Time[j]) < durat:
        
        if j < dim : 
        
          if Sigma[i] < Sigma[j] : Pulse[i] = -1
          else : Pulse[j] = -1
  
    if Pulse[i]==-1: continue
    
    # Set a code to the events that aren't in a pulse
    if Pulse[i]==0: 
      code += 1
      Pulse[i] = code
      
    # Defines for a certain DM a range of events that can be grouped
    if DM[i] != DM_new:
      
      j_min = 0
      j_max = dim
      
      DM_new = DM[i]
      
      if DM_new < 40.49: step = 0.01
        
      elif DM_new < 141.69: step = 0.05
        
      else: step = 0.1
        
      step_min = step - float_err
      
      step_max = step * (n_steps + 1) + float_err
      
      
      #find the minimum and maximum event in the range
      for j in range(i+1,dim):
        
        dDM = DM[j] - DM_new

        if dDM > step_max:
          
          j_max = j
          
          break

        if dDM > step_min: 
          
          if j_min == 0: j_min = j
          
    empty = 0
    
    if j_min > 0:

      # Gives a code to the next event in the pulse
      for j in range(j_min,j_max):

        if abs(Time[i]-Time[j]) < durat:   #MAYBE add a condition on SNR (attention: dSNR depends on SNR!) 
          
          if Pulse[j] == -1: continue
          
          if Pulse[j] > 0: 
            
            Pulse[j] = -1
            continue
          
          if empty == 0:
            
            Pulse[j] = Pulse[i]
            SNR_max = j
            empty = 1
            DM_min = DM[j]
            
          else:
            
            if DM[j] > DM_min: break
                        
            if Sigma[j] > Sigma[SNR_max]:
              
              Pulse[SNR_max] = -1
              SNR_max = j
              Pulse[j] = Pulse[i]
              
            else:
              
              Pulse[j] = -1
         
  return
  
  
#-------------------------
# Compares different beams
#-------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def Compare(float[::1] DM_c_l not None,
            float[::1] dDM_l not None,
            float[::1] Time_c_l not None,
            float[::1] dTime_l not None,
            float[::1] Sigma_l not None,
            signed char[::1] Pulse_l not None,
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
          
          if abs(Sigma_l[i]-Sigma_r[j]) < TollSigma:
                        
            Pulse_l[i] += 1
            Pulse_r[j] += 1
          
      elif sign < 0.: 
        
        break
      
  return
  
  
  
  
  
  
  
  
#-------------------------
# Compares repeated pulses
#-------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def Compare_candidates(float[::1] DM not None,
            float[::1] Time not None,            
            long[::1] idx not None,            
            long[::1] cand not None):

  cdef:
    unsigned int i, j
    unsigned int dim = len(DM)
    float float_err = 0.0001
    
  # Compare each candidate
  for i in range(0,dim):
    
    for j in range(dim-1,i,-1):
      
      if abs(DM[j]-DM[i]) < 1.:
      
        if (abs(Time[j]-Time[i]) < 1.) | (Time[i]<float_err) | (Time[j]<float_err):
        
          cand[i] = idx[j]
          
          break

      
  return
  
  
  
  
@cython.boundscheck(False)
@cython.wraparound(False)
def time_span(float[::1] DM not None,
            float[::1] Time not None,          
            signed char[::1] cand not None):

  cdef:
    unsigned int i, j
    unsigned int dim = len(DM)

  # Compare each candidate
  for i in range(0,dim):
    
    for j in range(i+1,dim):
      
      if Time[j]-Time[i] < 10.:
      
        if abs(DM[j]-DM[i]) < 1.:
        
          cand[i] = 1
          cand[j] = 1
          
          break
        
      else: break
      
  return  


  

  
'''
  
#------------------
# RFI local filters
#------------------  
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef int pulses_apply(float[::1] Sigma_or,
                       float[::1] Time_or,
                       float[::1] DM_or
                       float[::1] Duration_or):
  
  cdef:
    int dim = Sigma_or.shape[0]
    float s1[dim], s2[dim], s[dim]
    float Sigma[?], Time[?], DM[?], Duration[?]
  
  s1 = Sigma_or - Sigma_or.shift(-1)
  s1.iloc[-1] = 0
  s2 = Sigma_or - Sigma_or.shift(1)
  s2.iloc[0] = 0
  s = pd.concat((s1[s1<s2],s2[s2<=s1]))

#modificare in type declaration
  Sigma = Sigma_or[s>-5]
  Time = Time_or[s>-5]
  DM = DM_or[s>-5]
  Duration = Duration_or[s>-5]

  if Sigma.shape[0] < 5: return 1
  return np.sum((
    np.mean(np.fabs( Sigma - Sigma.shift(-1) ) / Sigma) < FILTERS['sigma_scatter'],
    (np.mean(np.abs(Time - Time.shift(1))) > FILTERS['cum_scatter']) |
    (np.std(Time - Time.shift(1)) > FILTERS['std_scatter']),
    sigma_std_largest(Sigma) | fit0(DM,Sigma) | fit1(DM,Sigma),
    SNR_simmetric(DM,Sigma) > FILTERS['flat_SNR_simmetric'],
    bright_events_abs(DM, Sigma) > FILTERS['bright_extremes_abs'],
    bright_events_rel(DM,Sigma) > FILTERS['bright_extremes_rel'],
    pulse_simmetric(DM, Sigma) < FILTERS['pulse_simmetric'],
    flat_SNR_extremes(Sigma) < FILTERS['flat_SNR_extremes'],
    number_events(DM, Sigma, Duration) < FILTERS['number_events'],
    monotonic(Sigma) < FILTERS['monotonic'],
    sigma_jumps(Sigma) > FILTERS['sigma_jumps'],
    fit1_brightest(DM, Sigma) < FILTERS['fit1_brightest']))

cdef int sigma_std_largest(float Sigma):
  cdef float sigma_larg[Sigma.size*2/3]
  
  sigma_larg = Sigma.nlargest(Sigma.size*2/3)
  if sigma_larg.size < 20: return 0
  if sigma_larg.max() < 8: return np.std(sigma_larg) < FILTERS['sigma_std_largest_weak']
  else: return np.std(sigma_larg) < FILTERS['sigma_std_largest']

cdef int fit0(float x, float y):
  cdef double p[1]
  
  if x.size < 20: return 0
  p = np.polyfit(x, y, 0)
  if y.max() < 8: return np.sum((np.polyval(p, x) - y) ** 2) / x.size < FILTERS['flat_fit0_weak']
  else: return np.sum((np.polyval(p, x) - y) ** 2) / x.size < FILTERS['flat_fit0']

cdef int fit1(float x, float y):
  cdef double p[2]

  if x.size<20: return 0
  p = np.polyfit(x, y, 1)
  if y.max() < 8: return np.sum((np.polyval(p, x) - y) ** 2) / x.size < FILTERS['flat_fit1_weak']
  else: return np.sum((np.polyval(p, x) - y) ** 2) / x.size < FILTERS['flat_fit1']

cdef int pulse_simmetric(float DM, float Sigma):
  cdef float DM_c, xl[?], yl[?], xr[?], yr[?]
  cdef double ml, mr
  
  DM_c = DM.loc[Sigma.idxmax()]
  xl = DM[DM <= DM_c]
  yl = Sigma[DM <= DM_c]
  ml = np.polyfit(x, y, 1)[0]
  xr = DM[DM >= DM_c]
  yr = Sigma[DM >= DM_c]
  mr = np.polyfit(x, y, 1)[0]
  return np.min((-ml/mr, -mr/ml))

cdef int number_events(float DM, float Sigma, float Duration):
  cdef float sig[?], dm_c[?], sig_argmax, sig_max, lim_max, duration[?], dDM, y, diff_l, diff_r
  cdef int l, r
  cdef int dim = Sigma.shape[0] / 5
  
  if dim < 3: return 0
  sig = np.convolve(Sigma, np.ones(dim), mode='valid') / dim
  dm_c = DM.iloc[dim/2 : -int(dim-1.5)/2]
  sig_argmax = sig.argmax()
  sig_max = sig.max()
  try: lim_max = np.max((sig[:sig_argmax].min(), sig[sig_argmax:].min()))
  except ValueError: return 0
  lim_max += (sig_max - lim_max) / 5.
  l = np.where(sig[:sig_argmax] <= lim_max)[0][-1] + 1
  r = (np.where(sig[sig_argmax:] <= lim_max)[0] + sig_argmax)[0] - 1
  if (sig_argmax - l < 5) & (r - sig_argmax < 5): return 0
  duration = np.convolve(Duration, np.ones(dim), mode='valid') / dim
  duration = duration[sig_argmax]
  dDM = dm_c.iloc[sig_argmax] - dm_c.iloc[l]
  y = np.sqrt(np.pi)/2/(0.00000691*dDM*31.64/duration/0.13525**3)*special.erf(0.00000691*dDM*31.64/duration/0.13525**3)
  diff_l = lim_max / sig_max / y
  dDM = dm_c.iloc[r] - dm_c.iloc[sig_argmax]
  y = np.sqrt(np.pi)/2/(0.00000691*dDM*31.64/duration/0.13525**3)*special.erf(0.00000691*dDM*31.64/duration/0.13525**3)
  diff_r = lim_max / sig_max / y
  return np.nanmax((diff_l,diff_r))

cdef int monotonic(float Sigma):
  cdef int dim = Sigma.size
  cdef float sig[dim], l, r
  
  sig = np.convolve(Sigma, np.ones(dim / 5), mode='same') / dim * 5
  sig_max = sig.argmax()
  l = sig[:sig_max].size * 2/3
  r = sig[sig_max:].size * 2/3
  sig = sig[l:-r]
  if sig.size < 10: return 1
  sig_max = sig.argmax()
  sig = np.diff(sig)
  sig[sig_max:] *= -1
  return np.partition(sig, 1)[1]

cdef sigma_jumps(float Sigma):
  cdef float sig[?]
  
  sig = np.convolve(Sigma, np.ones(5), mode='same') / 5.
  sig_max = sig.argmax()
  sig = np.diff(sig)
  sig[sig_max:] *= -1
  return sig[sig < 0].size / float(sig.size)

cdef SNR_simmetric(float DM, float Sigma):
  cdef float DM_c
  
  DM_c = DM.loc[Sigma.idxmax()]
  l = Sigma[DM <= DM_c]
  r = Sigma[DM >= DM_c]
  return np.max((l.min(),r.min())) / Sigma.max()

cdef int fit1_brightest(float Sigma):
  cdef float sig[?], dm_cen[?], l[?], r[?], DM_c, lim_l, lim_r, y[?], x[?]
  cdef double p[2]
  
  sig = np.convolve(Sigma, np.ones(3), mode='valid') / 3.
  dm_cen = DM.iloc[3/2:-int(3-1.5)/2]
  sig = pd.Series(sig, index=dm_cen.index)
  DM_c = dm_cen.loc[sig.argmax()]
  l = sig[dm_cen <= DM_c]
  if l.size <= 4: return 10
  r = sig[dm_cen >= DM_c]
  if r.size < 4: return 10
  lim_l = l.min() + np.min((2., (l.max() - l.min()) / 4.))
  lim_r = r.min() + np.min((2., (r.max() - r.min()) / 4.))
  l = l[l > lim_l]
  r = r[r > lim_r]
  y = pd.concat((l, r))
  if y.size <= 5: return 10
  x = dm_cen.loc[y.index]
  p = np.polyfit(x, y, 1)
  return np.sum((np.polyval(p, x) - y) ** 2) / (x.size - 1)
  
#rimuove gli eventi piu' deboli a destra e sinistra.
cdef int bright_events_abs(float DM, float Sigma):
  cdef float DM_c, l[?], r[?], lim_l, lim_r
  
  DM_c = DM.loc[Sigma.idxmax()]
  l = Sigma[DM <= DM_c]
  if l.shape[0] <= 4: return 0
  r = Sigma[DM > DM_c]
  if r.shape[0] < 4: return 0
  lim_l = l.min() + np.min((2.,(l.max()-l.min())/4))
  lim_r = r.min() + np.min((2.,(r.max()-r.min())/4))
  l = l[l > lim_l]
  r = r[r > lim_r]
  Sigma = pd.concat((l,r))
  DM = DM.reindex_like(Sigma)
  if Sigma.shape[0] <= 5: return 0
  try: return np.max((Sigma[DM.argmin()], Sigma[DM.argmax()])) / Sigma.max()
  except ValueError: return 1

cdef int bright_events_rel(float DM, float Sigma):
  cdef float DM_c, Sigma_l[?], Sigma_r[?], DM_l[?], DM_r[?], lim_l, lim_r 
  
  if Sigma.max()<8: return 0
  DM_c = DM.loc[Sigma.idxmax()]
  Sigma_l = Sigma[DM <= DM_c]
  DM_l = DM[DM <= DM_c]
  if Sigma_l.shape[0] <= 4: return 0
  Sigma_r = Sigma[DM > DM_c]
  DM_r = DM[DM > DM_c]
  if Sigma_r.shape[0] < 4: return 0
  DM_r.sort(ascending=False)
  Sigma_r = Sigma_r.reindex_like(DM_r)
  lim_l = np.cumsum(Sigma_l - Sigma_l.iloc[0])
  lim_r = np.cumsum(Sigma_r - Sigma_r.iloc[0])
  Sigma_l = Sigma_l[lim_l >= Sigma.max() / 8.]
  Sigma_r = Sigma_r[lim_r >= Sigma.max() / 8.]
  Sigma = pd.concat((Sigma_l,Sigma_r))
  DM = DM.reindex_like(Sigma)
  if Sigma.shape[0] < 5: return 0
  else: return np.max((Sigma[DM.argmin()], Sigma[DM.argmax()])) / Sigma.max()

cdef int flat_SNR_extremes(float Sigma):                                            
  if Sigma.shape[0] < 30: return 0
  else: return np.max((Sigma.iloc[1],Sigma.iloc[-2])) / Sigma.max()




'''


