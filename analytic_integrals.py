##########################################################################
#                       ANALYTIC INTEGRALS                               #
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
def integral_3(Vr, rdg, E_kin, TX, res, res_kin, t):
    I = Vr * rdg / res_kin * (1./res * (np.exp((t-TX/2) * res) - 1)
                              + complex(0,1./E_kin) * (np.exp(-(t-TX/2) * res) - 1))
    return I

def integral_4(Vr, rdg, E_kin, TX, res, res_kin, t):
    I = Vr * rdg / (res * reskin) \
        * (1 - np.exp((t-TX/2) * reskin) - np.exp(-(t-TX/2) * res)
           + np.exp(1j*(t-TX/2) * E_kin))
