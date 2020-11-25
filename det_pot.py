#!/usr/bin/python

import numpy as np
import sciconv as sc
import wellenfkt as wf
import sys

mu = wf.red_mass_au(20.1797,20.1797)

another = raw_input('Do you want another value? (y/n)    ')
while another == 'y':
    De_eV = raw_input('How deep should the potential be in eV?    ')
    De_eV = float(De_eV)
    De_au = sc.ev_to_hartree(De_eV)
    n_n = raw_input('How many bound states do you want?    ')
    n_n = int(n_n)
    alpha_min = 2* np.sqrt(2*mu*De_au) / (1 + 2*n_n)
    print 'alpha_min = ', alpha_min
    alpha = raw_input('Choose alpha:   ')
    alpha = float(alpha)
    lambda_param = np.sqrt(2*mu*De_au) / alpha
    n_max = int(lambda_param - 0.5)
    print "n_max = ", n_max
    for n in range (0,n_max+1):
        ev = wf.eigenvalue(n,De_au,alpha,mu)
        print "Eigenvalue = ", ev, "n = ", n
else:
    sys.exit()

#def red_mass_au(mass1,mass2):
#    red_mass_gmol = mass1 * mass2 / (mass1 + mass2)
#    red_mass_au = sc.gmol_to_me(red_mass_gmol)
#    return red_mass_au
#
#def eigenvalue(n, De, alpha, red_mass):
#    En =   (n + 0.5) * alpha * np.sqrt(2*De / red_mass) \
#         - (n + 0.5)**2 * alpha**2 / (2* red_mass)
#    return En
#
#def sqrt_fact(real):
#    if (np.absolute(1 - real) < 1.0E-7):
#        return np.sqrt(real)
#    elif real < 1.0:
#        return np.sqrt(scipy.misc.factorial(real) )
#    else:
#        return np.sqrt(real) * sqrt_fact(real-1)


