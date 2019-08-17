##########################################################################
#                                   Potentials                           #
##########################################################################
# Purpose:                                                               #
#          - A python module for the conversion of numbers between       #
#            different unit systems.                                     #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer July 2019                                   #
##########################################################################

import scipy.constants as constants
import numpy as np
import sciconv as sc
#-------------------------------------------------------------------------
#   unit conversion

a_0 = constants.hbar / (constants.m_e * constants.c * constants.fine_structure)
E_h = constants.fine_structure**2 * constants.m_e * constants.c**2


def expr6(a,b,c,d,r_au):
    r = sc.bohr_to_angstrom(r_au)
    V = a* np.exp(-b*r) + c / r**6 + d
    V_au = sc.ev_to_hartree(V)
    return V_au

def hyperbel(a,b,r_au):
    r = sc.bohr_to_angstrom(r_au)
    V = a / r + b
    V_au = sc.ev_to_hartree(V)
    return V_au

def gammar6(a,r_au):
    Gamma_au = a / r_au**6
    Vr_au = np.sqrt(Gamma_au / 2 / np.pi)
    return Vr_au
    

##-------------------------------------------------------------------------
##      Distances, Areas, Volume
#def bohr_to_angstrom(d_au):
#    d_AA = d_au *a_0 * 1E10
#    return d_AA
