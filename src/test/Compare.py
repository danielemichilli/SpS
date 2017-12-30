import pandas as pd
import numpy as np
import time

import C_Compare

from Parameters import *

idL = 'L204720'


puls = pd.read_hdf('sp/SinglePulses.hdf5',idL+'_pulses')
puls = puls[puls.Pulse<=RFI_percent]


print puls[(puls.SAP==1)&(puls.Pulse==0)&(puls.BEAM>12)&(puls.DM<26.80)&(puls.DM>26.74)].sort('Sigma').tail(20).Sigma
print puls[(puls.SAP==1)&(puls.Pulse==0)&(puls.BEAM==12)&(puls.DM<26.80)&(puls.DM>26.74)].sort('Sigma').tail(20).Sigma

exit()

puls.Pulse = 0
puls = puls.loc[puls.N_events >= 5]

puls.Pulse = puls.Pulse.astype(np.int8)

coh = puls[puls.BEAM>12]
incoh = puls[puls.BEAM==12]

a=coh[(coh.Pulse==0)&(coh.DM<26.80)&(coh.DM>26.74)].shape[0]
a0=coh[(coh.Pulse==0)&(coh.DM>35)|(coh.DM<15)].shape[0]


#sap0 = puls[(puls.SAP==0)&(puls.Pulse<=RFI_percent)].ix[:,['DM_c','dDM','Time_c','Duration','Sigma','Pulse']]
#sap0['Time_low'] = sap0.Time_c-sap0.Duration
#sap0.sort('Time_low',inplace=True)

#sap1 = puls[(puls.SAP==1)&(puls.Pulse<=RFI_percent)].ix[:,['DM_c','dDM','Time_c','Duration','Sigma','Pulse']]
#sap1['Time_low'] = sap1.Time_c-sap1.Duration
#sap1.sort('Time_low',inplace=True)

#sap2 = puls[(puls.SAP==2)&(puls.Pulse<=RFI_percent)].ix[:,['DM_c','dDM','Time_c','Duration','Sigma','Pulse']]
#sap2['Time_low'] = sap2.Time_c-sap2.Duration
#sap2.sort('Time_low',inplace=True)

#C_Compare.Compare(sap0.index.values,sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.Duration.values,sap0.Sigma.values,sap0.Pulse.values,\
                #sap1.index.values,sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.Duration.values,sap1.Sigma.values,sap1.Pulse.values,code)

#C_Compare.Compare(sap0.index.values,sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.Duration.values,sap0.Sigma.values,sap0.Pulse.values,\
                #sap2.index.values,sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.Duration.values,sap2.Sigma.values,sap2.Pulse.values,code)

#C_Compare.Compare(sap1.index.values,sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.Duration.values,sap1.Sigma.values,sap1.Pulse.values,\
                #sap2.index.values,sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.Duration.values,sap2.Sigma.values,sap2.Pulse.values,code)


#puls.Pulse.loc[(puls.SAP==0)&(puls.Pulse<=RFI_percent)]=sap0.Pulse

#puls.Pulse.loc[(puls.SAP==1)&(puls.Pulse<=RFI_percent)]=sap1.Pulse

#puls.Pulse.loc[(puls.SAP==2)&(puls.Pulse<=RFI_percent)]=sap2.Pulse


time0 = time.clock()
for sap in range(0,3):
  IB = incoh[(incoh.SAP==sap)&(incoh.Pulse<=RFI_percent)].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  IB['Time_low'] = IB.Time_c-IB.dTime
  IB.sort('Time_low',inplace=True)
  
  for beam in range(13,74):
    print "It's being compared sap{} beam{}".format(sap,beam)
    
    CB = coh[(coh.SAP==sap)&(coh.BEAM==beam)&(coh.Pulse<=RFI_percent)].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
    CB['Time_low'] = CB.Time_c-CB.dTime
    CB.sort('Time_low',inplace=True)

    C_Compare.Compare(CB.index.values,CB.DM_c.values,CB.dDM.values,CB.Time_c.values,CB.dTime.values,CB.Sigma.values,CB.Pulse.values,\
                      IB.index.values,IB.DM_c.values,IB.dDM.values,IB.Time_c.values,IB.dTime.values,IB.Sigma.values,IB.Pulse.values,np.int8(0))  
      
    coh.Pulse.loc[(coh.SAP==sap)&(coh.BEAM==beam)&(coh.Pulse<=RFI_percent)] = CB.Pulse
  incoh.Pulse.loc[(incoh.SAP==sap)&(incoh.Pulse<=RFI_percent)] = IB.Pulse

print 'Time for comparison: {} s'.format(time.clock()-time0)

puls = pd.concat([coh,incoh])


puls = puls[puls.Pulse==0]

b=coh[(coh.Pulse==0)&(coh.DM<26.80)&(coh.DM>26.74)].shape[0]
b0=coh[(coh.Pulse==0)&(coh.DM>35)|(coh.DM<15)].shape[0]

print 'RFI emilinato: ', 100.*(1.-b0/float(a0)), "%"
print 'Puls eliminati: ', 100.*(1.-b/float(a)), "%"

#print 'dimensione prima e dopo RFI:  ',a0,b0
#print 'dimensione prima e dopo Puls: ',a,b