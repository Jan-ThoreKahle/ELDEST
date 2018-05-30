#!/usr/bin/python

##########################################################################
#                                    ELDEST                              #
#        Investigating Electronic Decay Processes with Streaking         #
##########################################################################
# Purpose:                                                               #
#          - A program to simulate the streaking process of electronic   #
#            decay processes.                                            #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer May 2018                                    #
##########################################################################

import scipy
import scipy.integrate as integrate
import numpy as np
import sciconv
import complex_integration as ci

#-------------------------------------------------------------------------
# Input parameters

rdg           = 0.5           # transition dipole moment into the resonant state
cdg           = 0.5           # transition dipole moment into any continuum state

# parameters of the investigated system
# the ground state energy is being defined as Eg = 0
Er_eV         = 14.0          # resonance energy in eV
E_kin_eV      = 2.0           # kinetic energy of secondary electron
E_fin_eV      = 12.0          # final state energy in eV

Gamma_eV      = 0.5           # electronic decay width of the resonant state

# laser parameters
Omega_min_eV  = 10.0          # scanning XUV pulse from Omega_min-eV to
Omega_max_eV  = 18.0          #
TX_s          = 100E-18       # duration of the XUV pulse in seconds
A0X           = 1.0           # amplitude of the XUV pulse

omega_eV      = 1.0           # IR pulse
TL_s          = 1.0E-12       # duration of the IR streaking pulse
A0L           = 1.0           # amplitude of the IR pulse
delta_t_s     = 1.0E-15       # time difference between the maxima of the two pulses

# parameters of the simulation
tmax_s        = 1.0E-10       # simulate until time tmax in seconds
timestep_s    = 10E-18        # evaluate expression every timestep_s seconds 
Omega_step_eV = 0.1           # energy difference between different evaluated Omegas
#-------------------------------------------------------------------------

print 'Hello World'
Omega= 13.5

#-------------------------------------------------------------------------
# Definitions of reusable functions
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_au          = sciconv.ev_to_hartree(Er_eV)
Er_kin_au      = sciconv.ev_to_hartree(E_kin_eV)
Er_fin_au      = sciconv.ev_to_hartree(E_fin_eV)

Gamma_au       = sciconv.ev_to_hartree(Gamma_eV)

# laser parameters
Omega_min_au  = sciconv.ev_to_hartree(Omega_min_eV)
Omega_max_au  = sciconv.ev_to_hartree(Omega_max_eV)
TX_au         = sciconv.second_to_atu(TX_s)

omega_au      = sciconv.ev_to_hartree(omega_eV)
TL_au         = sciconv.second_to_atu(TL_s)
delta_t_au    = sciconv.second_to_atu(delta_t_s)

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
Omega_step_au = sciconv.ev_to_hartree(Omega_step_eV)

#-------------------------------------------------------------------------
# physical defintions of functions
# XUV pulse
f  = lambda t1: 1./4 (np.exp(2j*np.pi*t1/TX_au) + 2 + np.exp(-2j*np.pi*t1/TX_au) )
fp = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi*t1/TX_au)
                                     + np.exp(-2j*np.pi*t1/TX_au) )
FX = lambda t1: - A0X * np.cos(Omega_au * t1) * fp(t1) \
                + A0X * Omega_au * np.sin(Omega_au * t1) * f(t1)

#-------------------------------------------------------------------------
## very important: The first Variable in the definition of the function marks the inner
## integral, while the second marks the outer integral.
## If any limit is replaced by the integration variable of the outer integral,
## this is always specified as x, never as the actual name of the variable.
f = lambda t2, t1: np.exp(t1 * complex(Gamma_eV/2,Er_eV)) * np.exp(t2 * complex(Gamma_eV/2, Er_eV + E_kin_eV))
##f = lambda t2, t1: np.exp(t1 * Gamma_eV/2) * np.exp(t2 * Gamma_eV/2)
#


fun1 = lambda t1: np.exp(t1 * complex(Gamma_eV/2,Er_eV))
fun2 = lambda t2: np.exp(t2 * complex(Gamma_eV/2, Er_eV + E_kin_eV))

I = ci.complex_double_quadrature(fun1,fun2, -TX_s/2, TX_s/2, lambda x: x, lambda x: TX_s/2)

print I

print sciconv.bohr_to_angstrom(2)
