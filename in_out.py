##########################################################################
#                     INPUT AND OUTPUT ROUTINES                          #
##########################################################################
# Purpose:                                                               #
#          - A python module handling input and output.                  #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer May 2018                                    #
##########################################################################

import sciconv
import numpy as np
from sys import exit

#-------------------------------------------------------------------------
#   input
def check_input(Er, E_kin, E_fin, Gamma,
                Omega_min, Omega_max, TX, n_X, A0X,
                omega, TL, A0L, delta_t,
                tmax, timestep, Omega_step):
    print 'Input Check'

    if (TL/2 > (delta_t + TX/2)):
        exit('Warning: TL/2 > delta_t + TX/2' + '\n'
             + 'Stopping Script')

    if (Omega_min > Omega_max):
        exit('Warning: Omega_min > Omega_max' + '\n'
             + 'Stopping Script')

    if (TX < n_X * 2 * np.pi / Omega_min):
        exit('Warning: TX is too short to cover for the minimum energy requested'
             + '\n'
             + 'Stopping Script')

    print 'Input fullfills requirements'

    return 0
    

#-------------------------------------------------------------------------
#   output
def prep_output(I, Omega_au, t_au):
    square = np.absolute(I)**2
    Omega_eV = sciconv.hartree_to_ev(Omega_au)
    t_s = sciconv.atu_to_second(t_au)
    string = str(Omega_eV) + '   ' + format(t_s, '.18f') + '   ' + str(square)
    return string

def doout(t_au, outlines):
    # output filename will give the time in ps
    t_s = sciconv.atu_to_second(t_au)
    t_ps = t_s * 1E12
    if t_ps < 0.0:
        t_ps = np.absolute(t_ps)
        filename = 'm' + format(t_ps, '.8f') + '.dat'
    else:
        filename = format(t_ps, '.8f') + '.dat'
    outfile = open(filename, mode='w')
    res_lines = '\n'.join(outlines)
    outfile.write(res_lines)
    outfile.close

def doout_1f(filename, outlines):
    res_lines = '\n'.join(outlines)
    res_lines = res_lines + '\n' + '' + '\n'
    filename.write(res_lines)
    #outfile.close
