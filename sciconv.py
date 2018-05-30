##########################################################################
#                                   SCICONV                              #
##########################################################################
# Purpose:                                                               #
#          - A python module for the conversion of numbers between       #
#            different unit systems.                                     #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer May 2018                                    #
##########################################################################

import scipy.constants as constants
#-------------------------------------------------------------------------
#   unit conversion
#      Distances, Areas, Volume

a_0 = constants.hbar / (constants.m_e * constants.c * constants.fine_structure)

def bohr_to_angstrom(d_au):
    d_AA = d_au *a_0 * 1E10
    return d_AA

def angstrom_to_bohr(d_AA):
    d_au = d_AA / a_0 * 1E-10
    return d_au

def bohr_to_meter(d_au):
    d_m = d_au * a_0
    return d_m

def meter_to_bohr(d_m):
    d_au = d_m / a_0
    return d_au

def megabarn_to_sqmeter(sigma_mb):
    sigma_m2 = sigma_mb * 1E-22
    return sigma_m2

def sqmeter_to_megabarn(sigma_m2):
    sigma_mb = sigma_m2 * 1E22

#      Energies
#      Time
