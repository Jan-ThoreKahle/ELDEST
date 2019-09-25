#!/usr/bin/python

##########################################################################
#                             HARMOSZI                                   #
##########################################################################
# Purpose:                                                               #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer August 2019                                 #
##########################################################################

import scipy.misc
import scipy.integrate as integrate
#from scipy.signal import argrelextrema
#from scipy.special import erf
from scipy.special import gamma
from scipy.special import eval_genlaguerre
import numpy as np
import sciconv as sc
#import complex_integration as ci
#import in_out
import sys
import warnings
import potentials


#infile = sys.argv[1]
#print infile

#-------------------------------------------------------------------------
# open outputfile
#outfile = open("eldest.out", mode='w')

#outfile.write("The results were obtained with time_measure_split.py \n")
#-------------------------------------------------------------------------

n = 0

#-------------------------------------------------------------------------
# Parameters for potentials
#-------------------------------------------------------------------------
gs_de      =  0.003013
gs_a       =  1.105308
gs_Req     =  5.918877
gs_const   = -0.005044

u_de       =  0.096727
u_a        =  0.447152
u_Req      =  4.523888
u_const    = -256.233965

g_de       = 0.079073
g_a        = 0.143584
g_Req      = 10.003604
g_const    = -256.182536

mass1 = 20.1797
mass2 = 20.1797

red_mass_gmol = mass1 * mass2 / (mass1 + mass2)
red_mass = sc.gmol_to_me(red_mass_gmol)
print "redmass = ", red_mass

def eigenvalue(n, omega):
    En =   (n + 0.5) * omega 
    return En

def sqrt_fact(real):
    if (np.absolute(1 - real) < 1.0E-7):
        return np.sqrt(real)
    elif real < 1.0:
        return np.sqrt(scipy.misc.factorial(real) )
    else:
        return np.sqrt(real) * sqrt_fact(real-1)

print "Grundzustand:"
print "Potentialtiefe", gs_de
print eigenvalue(n, gs_de, gs_a, red_mass)


De = gs_de
alpha = gs_a
Req = gs_Req
n_gs = 0

lambda_param_gs = np.sqrt(2*red_mass*gs_de) / gs_a
z_gs = lambda R: 2* lambda_param_gs * np.exp(-gs_a * (R - gs_Req))
Nn_gs = np.sqrt(scipy.misc.factorial(n_gs) * (2*lambda_param_gs - 2*n_gs - 1)
                / gamma(2*lambda_param_gs - n_gs))

psi_n_gs = lambda R: Nn_gs *( 1
           * z_gs(R)**(lambda_param_gs - n_gs - 0.5) 
           * np.exp(-z_gs(R)/2) 
           * eval_genlaguerre(n_gs, 2*lambda_param_gs - 2*n_gs -1, z_gs(R))
            )

De = u_de
alpha = u_a
Req = u_Req

#De = g_de
#alpha = g_a
#Req = g_Req

lambda_param = np.sqrt(2*red_mass*De) / alpha
z = lambda R: 2* lambda_param * np.exp(-alpha * (R - Req))
#Nn = np.sqrt(scipy.misc.factorial(n) * (2*lambda_param - 2*n - 1)
#                / gamma(2*lambda_param - n))
Nn = np.sqrt(scipy.misc.factorial(n) * (2*lambda_param - 2*n - 1)) \
                / sqrt_fact(2*lambda_param - n - 1)
print "Nn = ", Nn

#psi_n = lambda R: Nn *( 1
#           * z(R)**(lambda_param - n - 0.5) 
#           * np.exp(-z(R)/2) 
#           * eval_genlaguerre(n, 2*lambda_param - 2*n -1, z(R))
#            )

# try eq 50 for n = 0
s = 2*lambda_param - 2*n - 1
print "s = ", s
print "sqrt = ", 1.0 / sqrt_fact(2*lambda_param-2)
psi_n = lambda R:  ( 1.0 
                        / sqrt_fact(2*lambda_param-2)
                        * z(R)**(s/4)
                        * z(R)**(s/4)
                        * np.exp(-z(R) / 2)
                        )

print "z(Req) = ", z(Req)
print "psi_n = ", psi_n(Req)
#print "psi_n = ", psi_n

func = lambda R: psi_n(R) * psi_n(R)

R_min = sc.angstrom_to_bohr(0.0)
R_max = sc.angstrom_to_bohr(30.0)

print "R_min = ", R_min
print "R_max = ", R_max

FC = integrate.quad(func, R_min, R_max)
print "FC = ", FC[0]
