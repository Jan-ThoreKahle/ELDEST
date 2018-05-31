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
                              - complex(0,1./E_kin) * (np.exp(-(t-TX/2) * res) - 1))
    return I
