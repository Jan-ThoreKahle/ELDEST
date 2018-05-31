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
    I = - Vr * rdg / res_kin * (1./res * (np.exp((t-TX/2) * res) - 1)
                              + complex(0,1./E_kin)
                                * (np.exp(-1j*(t-TX/2)*E_kin) - 1))
    return I

# Integral 4
def integral_4(Vr, rdg, E_kin, TX, res, res_kin, t):
    I = - Vr * rdg / (res * res_kin) \
          * (1 - np.exp((t-TX/2) * res_kin) - np.exp(-(t-TX/2) * res)
             + np.exp(1j*(t-TX/2) * E_kin))

# Integral 6 and 12 (they are the same)
def integral_6_12(Vr, rdg, E_kin, TX, TL, delta, res, res_kin):
    dm = delta - TL/2
    I = - Vr * rdg / res_kin * (1./res * (np.exp((dm-TX/2) * res) - 1)
                            + complex(0,1./E_kin) * (np.exp(-(dm-TX/2) * res) - 1))
    return I

# Integral 7 and 13
def integral_7_13(Vr, rdg, E_kin, TX, TL, delta, res, res_kin):
    dm = delta - TL/2
    I = - Vr * rdg / (res * res_kin) \
          * (1 - np.exp((dm-TX/2) * res_kin) - np.exp(-(dm-TX/2) * res)
             + np.exp(1j*(dm-TX/2) * E_kin))
    return I

# Integral 8
def integral_8(Vr, rdg, E_kin, TX, TL, delta, res, res_kin, t):
    dm = delta - TL/2
    I = - Vr * rdg / res_kin \
        * (1./res * (np.exp((t-TX/2) * res) - np.exp((dm-TX/2) * res))
           + complex(0,1./E_kin) * (np.exp(-(t-TX/2) * res)
                                    - np.exp(-(dm-TX/2) * E_kin)))
    return I

# Integral 9
def integral_9(Vr, rdg, E_kin, TX, TL, delta, res, res_kin, t):
    dm = delta - TL/2
    I = - Vr * rdg / (res * res_kin) \
          * (np.exp((t-dm) * res)
             - np.exp((t-TX/2) * res) * np.exp(1j*(dm-TX/2) * E_kin)
             - 1
             + np.exp((dm-TX/2) * res_kin))
    return I

# Integral 10 but without the extra phase facot coming from the second pulse
def integral_10(Vr, rdg, E_kin, TX, TL, delta, res, res_kin, t):
    dm = delta - TL/2
    I = - Vr * rdg / (res * res_kin) \
          * (np.exp(1j*(t-dm) * E_kin)
             - np.exp((t-dm) * res)
             - np.exp(-(t-dm) * res_kin)
             + 1)
    return I

# Integral 15
def integral_15(Vr, rdg, E_kin, TX, TL, delta, res, res_kin):
    dm = delta - TL/2
    dp = delta + TL/2
    I = - Vr * rdg / (res * res_kin) \
          * (np.exp(TL * res)
             - np.exp((dp-TX/2) * res) * np.exp(1j*(dm-TX/2) * E_kin)
             - 1
             + np.exp((dm-TX/2) * res_kin))
    return I
