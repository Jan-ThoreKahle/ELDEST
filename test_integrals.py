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
import analytic_integrals as ai
import in_out

#-------------------------------------------------------------------------
# Input parameters

rdg_au        = 0.5           # transition dipole moment into the resonant state
cdg           = 0.5           # transition dipole moment into any continuum state

# parameters of the investigated system
# the ground state energy is being defined as Eg = 0
Er_eV         = 14.0          # resonance energy in eV
E_kin_eV      = 2.0           # kinetic energy of secondary electron
E_fin_eV      = 12.0          # final state energy in eV

Gamma_eV      = 0.5           # electronic decay width of the resonant state

# laser parameters
Omega_min_eV  = 42.0          # scanning XUV pulse from Omega_min-eV to
Omega_max_eV  = 46.0          #
TX_s          = 100E-18       # duration of the XUV pulse in seconds
A0X           = 1.0           # amplitude of the XUV pulse

omega_eV      = 1.0           # IR pulse
TL_s          = 1.0E-14       # duration of the IR streaking pulse
A0L           = 1.0           # amplitude of the IR pulse
delta_t_s     = 1.0E-14       # time difference between the maxima of the two pulses

# parameters of the simulation
tmax_s        = 5.0E-14       # simulate until time tmax in seconds
timestep_s    = 2E-15        # evaluate expression every timestep_s seconds 
Omega_step_eV = 0.2           # energy difference between different evaluated Omegas
#-------------------------------------------------------------------------

in_out.check_input(Er_eV, E_kin_eV, E_fin_eV, Gamma_eV,
                   Omega_min_eV, Omega_max_eV, TX_s, A0X,
                   omega_eV, TL_s, A0L, delta_t_s,
                   tmax_s, timestep_s, Omega_step_eV)

#-------------------------------------------------------------------------
# Definitions of reusable functions
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_au          = sciconv.ev_to_hartree(Er_eV)
E_kin_au       = sciconv.ev_to_hartree(E_kin_eV)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)

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

p_au          = np.sqrt(2*E_kin_au)
VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))

res_kin = complex(Gamma_au/2,Er_au + E_kin_au)
res     = complex(Gamma_au/2,Er_au)

#-------------------------------------------------------------------------
# physical defintions of functions
# XUV pulse
f  = lambda t1: 1./4 * ( np.exp(2j * np.pi * t1 / TX_au)
                      + 2
                      + np.exp(-2j * np.pi * t1 /TX_au) )

fp = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi*t1/TX_au)
                                     + np.exp(-2j*np.pi*t1/TX_au) )

FX = lambda t1: - A0X * np.cos(Omega_au * t1) * fp(t1) + A0X * Omega_au * np.sin(Omega_au * t1) * f(t1)

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) / TL_au)**2
integ_IR = lambda t3: (p_au + A_IR(t3))**2

#-------------------------------------------------------------------------
# technical defintions of functions
fun_inf_TX2_1 = lambda t1: np.exp(t1 * complex(Gamma_au/2,Er_au)) * FX(t1)
fun_inf_TX2_2 = lambda t2: np.exp(t2 * complex(Gamma_au/2, Er_au + E_kin_au))

fun_TX2_delta_1 = lambda t1: np.exp(t1 * complex(Gamma_au/2,Er_au))
fun_TX2_delta_2 = lambda t2: np.exp(t2 * complex(Gamma_au/2, Er_au + E_kin_au))

#-------------------------------------------------------------------------
# very important: The first Variable in the definition of the function marks the inner
# integral, while the second marks the outer integral.
# If any limit is replaced by the integration variable of the outer integral,
# this is always specified as x, never as the actual name of the variable.
#

#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2


#-------------------------------------------------------------------------
# constant integrals, they are independent of both Omega and t
# Integral 6 - 12
integral_6_12 = ai.integral_6_12(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin)
J = integral_6_12
print 'integral 6 12 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 TX_au/2, delta_t_au - TL_au/2,
                                 lambda x: x, lambda x: TX_au/2)
I_TX2_delta_t1 = I[0] * np.exp(1j * E_kin_au * TX_au/2) * rdg_au * VEr_au
print 'integral 6 quadrature = ', I_TX2_delta_t1

# Integral 7 - 13
integral_7_13 = ai.integral_7_13(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin)
J = integral_7_13
print 'integral 7 12 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 TX_au/2, delta_t_au - TL_au/2,
                                 lambda x: TX_au/2,
                                 lambda x: delta_t_au -TL_au/2)
I_TX2_delta_TX2 = (I[0] * np.exp(1j * E_kin_au * (delta_t_au - TL_au/2))
                   * rdg_au * VEr_au)
print 'integral 7 quadrature = ', I_TX2_delta_TX2

# Integral 14
integral_14 = ai.integral_14(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin)
J = integral_14
print 'integral 14 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au - TL_au/2, delta_t_au + TL_au/2,
                                 lambda x: x,
                                 lambda x: TX_au/2)
I_deltam_deltap_t1 = (I[0] * np.exp(1j * E_kin_au * TX_au/2)
                * rdg_au * VEr_au)
print 'integral 14 quadrature = ', I_deltam_deltap_t1

# Integral 15
integral_15 = ai.integral_15(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin)
J = integral_15
print 'integral 15 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au - TL_au/2, delta_t_au + TL_au/2,
                                 lambda x: TX_au/2,
                                 lambda x: delta_t_au - TL_au/2)
I_deltam_deltap_TX2 = (I[0] * np.exp(1j * E_kin_au * (delta_t_au - TL_au/2))
                * rdg_au * VEr_au)
print 'integral 15 quadrature = ', I_deltam_deltap_TX2

# Integral 16
integral_16 = ai.integral_16(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin)
J = integral_16
print 'integral 16 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au - TL_au/2, delta_t_au + TL_au/2,
                                 lambda x: delta_t_au - TL_au/2,
                                 lambda x: delta_t_au + TL_au/2)
I_deltam_deltap_deltam = (I[0] * np.exp(1j * E_kin_au * (delta_t_au - TL_au/2))
#                          * np.exp(1j/2 * I_IR[0])
                          * rdg_au * VEr_au)
print 'integral 16 quadrature = ', I_deltam_deltap_deltam



#-------------------------------------------------------------------------
# time-dependent integrals
#-------------------------------------------------------------------------
# between the pulses
t_au = TX_au
# integrals 3 and 4 are independent of Omega, they are therefore
# evaluated before integral 2 and especially outside the loop
#Integral 3
integral_3 = ai.integral_3(VEr_au, rdg_au, E_kin_au, TX_au, res, res_kin, t_au)
J = integral_3
print 'integral 3 analytic = ', J

I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 TX_au/2, t_au,
                                 lambda x: x, lambda x: TX_au/2)
I_TX2_t_t1 = I[0] * np.exp(1j * E_kin_au * TX_au/2) * rdg_au * VEr_au
print 'integral 3 quadrature = ', I_TX2_t_t1

#Integral 4
integral_4 = ai.integral_3(VEr_au, rdg_au, E_kin_au, TX_au, res, res_kin, t_au)
J = integral_4
print 'integral 4 analytic = ', J


I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 TX_au/2, t_au,
                                 lambda x: TX_au/2, lambda x: t_au)
I_TX2_t_TX2 = I[0] * np.exp(1j * E_kin_au * t_au) * rdg_au * VEr_au
print 'integral 4 quadrature = ', I_TX2_t_TX2


#-------------------------------------------------------------------------
# during the IR pulse
t_au = delta_t_au -TL_au/2 + TX_au
#-------------------------------------------------------------------------
# integrals, that are independent of Omega
# Integral 8
integral_8 = ai.integral_8(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                           TX=TX_au, TL=TL_au, delta=delta_t_au,
                           res=res, res_kin=res_kin, t=t_au)
J = integral_8
print 'integral 8 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au - TL_au/2, t_au,
                                 lambda x: x,
                                 lambda x: TX_au/2)
I_delta_t_t1 = (I[0] * np.exp(1j * E_kin_au * TX_au/2)
                * rdg_au * VEr_au)
print 'integral 8 quadrature = ', I_delta_t_t1

# Integral 9
integral_9 = ai.integral_9(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                           TX=TX_au, TL=TL_au, delta=delta_t_au,
                           res=res, res_kin=res_kin, t=t_au)
J = integral_9
print 'integral 9 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au - TL_au/2, t_au,
                                 lambda x: TX_au/2,
                                 lambda x: delta_t_au - TL_au/2)
I_delta_t_TX2 = (I[0] * np.exp(1j * E_kin_au * (delta_t_au - TL_au/2))
                * rdg_au * VEr_au)
print 'integral 9 quadrature = ', I_delta_t_TX2

# Integral 10
integral_10 = ai.integral_10(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin, t=t_au)
J = integral_10
print 'integral 10 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au - TL_au/2, t_au,
                                 lambda x: delta_t_au - TL_au/2,
                                 lambda x: t_au)
#I_IR = integrate.quad(integ_IR, delta_t_au - TL_au/2, t_au)
I_delta_t_delta = (I[0] * np.exp(1j * E_kin_au * (delta_t_au - TL_au/2))
#                   * np.exp(1j/2 * I_IR[0])
                   * rdg_au * VEr_au)
print 'integral 10 quadrature = ', I_delta_t_delta
#I_IR = integrate.quad(integ_IR, delta_t_au - TL_au/2, t_au)



#-------------------------------------------------------------------------
# after the second pulse
t_au = delta_t_au + TL_au/2 + TX_au
#-------------------------------------------------------------------------
print 'after the second pulse'

# Omega independent integrals
#Integral 17
integral_17 = ai.integral_17(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin, t=t_au)
J = integral_17
print 'integral 17 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au + TL_au/2, t_au,
                                 lambda x: x,
                                 lambda x: TX_au/2)
I_deltap_t_t1 = (I[0] * np.exp(1j * E_kin_au * TX_au/2)
                 * rdg_au * VEr_au)
print 'integral 17 quadrature = ', I_deltap_t_t1

#Integral 18
integral_18 = ai.integral_18(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,

                             res=res, res_kin=res_kin, t=t_au)
J = integral_18
print 'integral 18 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au + TL_au/2, t_au,
                                 lambda x: TX_au/2,
                                 lambda x: delta_t_au - TL_au/2)
I_deltap_t_TX2 = (I[0] * np.exp(1j * E_kin_au * (delta_t_au - TL_au/2))
                  * rdg_au * VEr_au)
print 'integral 18 quadrature = ', I_deltap_t_TX2

#Integral 19
integral_19 = ai.integral_19(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin, t=t_au)
#integral_19 = integral_19 * I_IR[0]
J = integral_19
print 'integral 19 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au + TL_au/2, t_au,
                                 lambda x: delta_t_au - TL_au/2,
                                 lambda x: delta_t_au + TL_au/2)
I_deltap_t_deltam = (I[0] * np.exp(1j * E_kin_au * (delta_t_au - TL_au/2))
#                          * np.exp(1j/2 * I_IR[0])
                          * rdg_au * VEr_au)
print 'integral 19 quadrature = ', I_deltap_t_deltam


#Integral 20
integral_20 = ai.integral_20(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin, t=t_au)
#integral_20 = integral_20 * I_IR[0]
J = integral_20
print 'integral 20 analytic = ', J
I = ci.complex_double_quadrature(fun_TX2_delta_1,fun_TX2_delta_2,
                                 delta_t_au + TL_au/2, t_au,
                                 lambda x: delta_t_au + TL_au/2,
                                 lambda x: t_au)
I_deltap_t_deltam = (I[0] * np.exp(1j * E_kin_au * (delta_t_au - TL_au/2))
#                          * np.exp(1j/2 * I_IR[0])
                          * np.exp(1j * E_kin_au * (t_au
                                                    - (delta_t_au + TL_au/2)))
                          * rdg_au * VEr_au)
print 'integral 20 quadrature = ', I_deltap_t_deltam

