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
def read_input(inputfile, outfile):
    f = open(inputfile, 'r')
    
    lines = f.readlines()
    
    for line in lines:
        words = line.split()
        if (words[0] == 'rdg_au'):
            rdg_au = float(words[2])
            print 'rdg_au = ', rdg_au
            outfile.write('rdg_au = ' + str(rdg_au) + '\n')
    
    # energy parameters of the system
        elif (words[0] == 'Er_eV'):
            Er_eV = float(words[2])
            print 'Er_eV = ', Er_eV
            outfile.write('Er_eV = ' + str(Er_eV) + '\n')
        elif (words[0] == 'E_fin_eV'):
            E_fin_eV = float(words[2])
            print 'E_fin_eV = ', E_fin_eV
            outfile.write('E_fin_eV = ' + str(E_fin_eV) + '\n')
        elif (words[0] == 'tau_s'):
            tau_s = float(words[2])
            print 'tau_s = ', tau_s
            outfile.write('tau_s = ' + str(tau_s) + '\n')
    
    # exciting laser parameters
        elif (words[0] == 'Omega_eV'):
            Omega_eV = float(words[2])
            print 'Omega_eV = ', Omega_eV
            outfile.write('Omega_eV = ' + str(Omega_eV) + '\n')
        elif (words[0] == 'n_X'):
            n_X = float(words[2])
            print 'n_X = ', n_X
            outfile.write('n_X = ' + str(n_X) + '\n')
        elif (words[0] == 'I_X'):
            I_X = float(words[2])
            print 'I_X = ', I_X
            outfile.write('I_X = ' + str(I_X) + '\n')
        elif (words[0] == 'X_shape'):
            if (words[2] == 'sinsq'):
                X_sinsq = True
                X_gauss = False
                print 'X_shape = Sin**2'
                outfile.write('X_shape = Sin**2 \n')
            elif (words[2] == 'gauss'):
                X_sinsq = False
                X_gauss = True
                print 'X_shape = Gauss'
                outfile.write('X_shape = Gauss \n')
            else:
                print 'no XUV pulse shape selected'
                outfile.write('no XUV pulse shape selected \n')
        elif (words[0] == 'Xshape'):
            if (words[2] == 'infinite'):
                Xshape = 'infinite'
                outfile.write('Infinite XUV pulse selected \n')
                print 'Infinite XUV pulse selected'
            elif (words[2] == 'convoluted'):
                Xshape = 'convoluted'
                outfile.write('Convoluted XUV pulse selected \n')
                print 'Convoluted XUV pulse selected'
    
    # dressing laser parameters
        elif (words[0] == 'omega_eV'):
            omega_eV = float(words[2])
            print 'omega_eV = ', omega_eV
            outfile.write('omega_eV = ' + str(omega_eV) + '\n')
        elif (words[0] == 'n_L'):
            n_L = float(words[2])
            print 'n_L = ', n_L
            outfile.write('n_L = ' + str(n_L) + '\n')
        elif (words[0] == 'I_L'):
            I_L = float(words[2])
            print 'I_L = ', I_L
            outfile.write('I_L = ' + str(I_L) + '\n')
        elif (words[0] == 'delta_t_s'):
            delta_t_s = float(words[2])
            print 'delta_t_s = ', delta_t_s
            outfile.write('delta_t_s = ' + str(delta_t_s) + '\n')
        elif (words[0] == 'shift_step_s'):
            shift_step_s = float(words[2])
            print 'shift_step_s = ', shift_step_s
            outfile.write('shift_step_s = ' + str(shift_step_s) + '\n')
        elif (words[0] == 'phi'):
            phi = float(words[2])
            print 'phi = ', phi
            outfile.write('phi = ' + str(phi) + '\n')
        elif (words[0] == 'q'):
            q = float(words[2])
            print 'q = ', q
            outfile.write('q = ' + str(q) + '\n')
    
    # parameters of the simulation
        elif (words[0] == 'tmax_s'):
            tmax_s = float(words[2])
            print 'tmax_s = ', tmax_s
            outfile.write('tmax_s = ' + str(tmax_s) + '\n')
        elif (words[0] == 'timestep_s'):
            timestep_s = float(words[2])
            print 'timestep_s = ', timestep_s
            outfile.write('timestep_s = ' + str(timestep_s) + '\n')
        elif (words[0] == 'E_step_eV'):
            E_step_eV = float(words[2])
            print 'E_step_eV = ', E_step_eV
            outfile.write('E_step_eV = ' + str(E_step_eV) + '\n')
    
        elif (words[0] == 'E_min_eV'):
            E_min_eV = float(words[2])
            print 'E_min_eV = ', E_min_eV
            outfile.write('E_min_eV = ' + str(E_min_eV) + '\n')
        elif (words[0] == 'E_max_eV'):
            E_max_eV = float(words[2])
            print 'E_max_eV = ', E_max_eV
            outfile.write('E_max_eV = ' + str(E_max_eV) + '\n')

        elif (words[0] == 'integ'):
            if (words[2] == 'romberg'):
                integ = 'romberg'
                print 'Integration Scheme of the inner integral = Romberg'
                outfile.write('Integration Scheme of the inner integral = Romberg \n')
            elif (words[2] == 'quadrature'):
                integ = 'quadrature'
                print 'Integration Scheme of the inner integral = Gaussian Quadrature'
                outfile.write('Integration Scheme of the inner integral = Gaussian Quadrature \n')
            elif (words[2] == 'analytic'):
                integ = 'analytic'
                print 'Integration Scheme of the inner integral = analytic'
                outfile.write('Integration Scheme of the inner integral = analytic \n')
            else:
                print 'no integration scheme selected'
                outfile.write('no integration scheme selected \n')
            
        elif (words[0] == 'integ_outer'):
            if (words[2] == 'romberg'):
                integ_outer = 'romberg'
                print 'Integration Scheme of the outer integral = Romberg'
                outfile.write('Integration Scheme of the outer integral = Romberg \n')
            elif (words[2] == 'quadrature'):
                integ_outer = 'quadrature'
                print 'Integration Scheme of the outer integral = Gaussian Quadrature'
                outfile.write('Integration Scheme of the outer integral = Gaussian Quadrature \n')
            else:
                print 'no integration scheme selected'
                outfile.write('no integration scheme selected \n')
    
    f.close()
    return (rdg_au,
            Er_eV, E_fin_eV, tau_s,
            Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
            omega_eV, n_L, I_L, delta_t_s, shift_step_s, phi, q,
            tmax_s, timestep_s, E_step_eV,
            E_min_eV, E_max_eV,
            integ, integ_outer)

#-------------------------------------------------------------------------
def check_input(Er, E_fin, Gamma,
                Omega, TX, n_X, A0X,
                omega, TL, A0L, delta_t,
                tmax, timestep, E_step):
    print 'Input Check'

    if (E_fin > Omega):
        exit('Warning: E_fin > Omega' + '\n'
             + 'Stopping Script')

    print 'Input fullfills requirements'

    return 0
    

#-------------------------------------------------------------------------
#   output
def prep_output(I, Omega_au, t_au):
#    square = np.absolute(I)**2
    #print I
    Omega_eV = sciconv.hartree_to_ev(Omega_au)
    t_s = sciconv.atu_to_second(t_au)
    string = str(Omega_eV) + '   ' + format(t_s, '.18f') + '   ' + format(I, '.15e')
    return string

def prep_output_int(I, Omega_au, integer):
#    square = np.absolute(I)**2
    #print I
    Omega_eV = sciconv.hartree_to_ev(Omega_au)
    string = str(Omega_eV) + '   ' + str(integer) + '   ' + format(I, '.15e')
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
