#############################
#
# Cython Functions
#
# Functions written in cython
# are grouped here.
#
#############################

cimport cython

#---------------------------------
# Gives a pulse-code to each event
#---------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def Get_Group(double[::1] DM not None,
          double[::1] Sigma not None,
          double[::1] Time not None,
          long[::1] Pulse not None,
          float durat,
          unsigned int n_steps):
  
  cdef:
    unsigned int i, j, k, j_min, j_max, empty, SNR_max
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
      
        if (abs(Time[j]-Time[i]) < 1.):
        
          cand[i] = idx[j]

          break
          
        elif Time[i] < float_err:
          
          cand[j] = idx[i]
          
          break
          
        elif Time[j] < float_err:
          
          cand[i] = idx[j]
          
          break
      
  return