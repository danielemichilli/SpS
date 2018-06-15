import numpy as Num

from psr_constants import *

def rotate(arr, bins):
    """
    rotate(arr, bins):
        Return an array rotated by 'bins' places to the left
    """
    bins = bins % len(arr)
    if bins==0:
        return arr
    else:
        return Num.concatenate((arr[bins:], arr[:bins]))

def delay_from_DM(DM, freq_emitted):
    """
    Return the delay in seconds caused by dispersion, given
    a Dispersion Measure (DM) in cm-3 pc, and the emitted
    frequency (freq_emitted) of the pulsar in MHz.
    """
    if (type(freq_emitted)==type(0.0)):
        if (freq_emitted > 0.0):
            return DM/(0.000241*freq_emitted*freq_emitted)
        else:
            return 0.0
    else:
        return Num.where(freq_emitted > 0.0,
                         DM/(0.000241*freq_emitted*freq_emitted), 0.0)




