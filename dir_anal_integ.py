##########################################################################
#                       ANALYTIC INTEGRALS direct ionization             #
##########################################################################
# Purpose:                                                               #
#          - A python module to calculate integrals of complex numbers.  #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer May 2018                                    #
##########################################################################

import numpy as np
#-------------------------------------------------------------------------
#   integration
# Integral 3
def integral_3(cdg, E_kin, TX, t):
    I = 1j * cdg * complex(0,1./E_kin) * (1 - np.exp(1j*(t-TX/2)*E_kin))
    return I

# Integral 5 and 8
def integral_5_8(cdg, E_kin, TX, TL, delta, t):
    dm = delta - TL/2
    I = 1j * cdg \
        * complex(0,1./E_kin) * (np.exp(-1j*(dm-TX/2)*E_kin) - 1)

# Integral 6
def integral_6(cdg, E_kin, TX, TL, delta, t):
    dm = delta - TL/2
    I = 1j * cdg \
        * complex(0,1./E_kin) * (np.exp(-1j*(t-dm)*E_kin) - 1)
    return I

# Integral 9
def integral_9(cdg, E_kin, TX, TL, delta):
    dm = delta - TL/2
    I = 1j * cdg \
        * complex(0,1./E_kin) * (np.exp(-1j*TL*E_kin) - 1)
    return I

# Integral 10
def integral_10(cdg, E_kin, TX, TL, delta, t):
    dm = delta - TL/2
    I = 1j * cdg \
        * complex(0,1./E_kin) * (np.exp(-1j*TL*E_kin) - np.exp(1j*(t-dm-TL/2)*E_kin))
    return I

