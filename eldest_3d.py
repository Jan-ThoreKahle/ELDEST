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
import res_anal_integ as aires
import dir_anal_integ as aidir
import in_out

#-------------------------------------------------------------------------
# Input parameters

rdg_au        = 0.5           # transition dipole moment into the resonant state
cdg_au        = 0.5           # transition dipole moment into any continuum state

# parameters of the investigated system
# the ground state energy is being defined as Eg = 0
Er_eV         = 44.0          # resonance energy in eV
E_kin_eV      = 2.0           # kinetic energy of secondary electron
E_fin_eV      = 12.0          # final state energy in eV

#Gamma_eV      = 0.5           # electronic decay width of the resonant state
tau_s         = 2.0E-15       # lifetime

# laser parameters
Omega_min_eV  = 20.0          # scanning XUV pulse from Omega_min-eV to
Omega_max_eV  = 78.0          #
TX_s          = 700E-18       # duration of the XUV pulse in seconds
n_X           = 3
I_X           = 5.0E11        # intensity of the XUV pulse in W/cm^2
#A0X           = 1.0           # amplitude of the XUV pulse

omega_eV      = 1.0           # IR pulse
TL_s          = 1.0E-15       # duration of the IR streaking pulse
print TL_s
n_L           = 4
I_L           = 1.0E09        # intensity of the IR pulse in W/cm^2
#A0L           = 1.0           # amplitude of the IR pulse
delta_t_s     = 1.0E-14       # time difference between the maxima of the two pulses
print delta_t_s
phi           = 0

# parameters of the simulation
tmax_s        = 1.0E-15       # simulate until time tmax in seconds
timestep_s    = 50E-18        # evaluate expression every timestep_s seconds 
Omega_step_eV = 1.0           # energy difference between different evaluated Omegas
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# Definitions of reusable functions
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_au          = sciconv.ev_to_hartree(Er_eV)
E_kin_au       = sciconv.ev_to_hartree(E_kin_eV)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)

#Gamma_au       = sciconv.ev_to_hartree(Gamma_eV)
tau_au         = sciconv.second_to_atu(tau_s)
Gamma_au       = 1. / tau_au

# laser parameters
Omega_min_au  = sciconv.ev_to_hartree(Omega_min_eV)
Omega_max_au  = sciconv.ev_to_hartree(Omega_max_eV)
TX_au         = sciconv.second_to_atu(TX_s)
#TX_au         = n_X * 2 * np.pi / Omega_min_au
I_X_au        = sciconv.Wcm2_to_aiu(I_X)
E0X           = np.sqrt(I_X_au)
A0X           = E0X / Omega_min_au # this could be wrong and might have
                                   # to be evaluated for each Omega

omega_au      = sciconv.ev_to_hartree(omega_eV)
TL_au         = sciconv.second_to_atu(TL_s)
TL_au         = n_L * 2 * np.pi / omega_au
print TL_au/2
I_L_au        = sciconv.Wcm2_to_aiu(I_L)
E0L           = np.sqrt(I_L_au)
A0L           = E0L / omega_au
delta_t_au    = sciconv.second_to_atu(delta_t_s)
print delta_t_au

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
Omega_step_au = sciconv.ev_to_hartree(Omega_step_eV)

p_au          = np.sqrt(2*E_kin_au)
VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))

res_kin = complex(Gamma_au/2,Er_au + E_kin_au)
res     = complex(Gamma_au/2,Er_au)

#-------------------------------------------------------------------------
in_out.check_input(Er_au, E_kin_au, E_fin_au, Gamma_au,
                   Omega_min_au, Omega_max_au, TX_au, n_X, A0X,
                   omega_au, TL_au, A0L, delta_t_au,
                   tmax_au, timestep_au, Omega_step_au)
#-------------------------------------------------------------------------
# open outputfile
outfile = open("eldest.out", mode='w')
pure_out = open('full.dat', mode='w')
#-------------------------------------------------------------------------
# physical defintions of functions
# XUV pulse
f_t  = lambda tau: 1./4 * ( np.exp(2j * np.pi * (t_au - tau) / TX_au)
                      + 2
                      + np.exp(-2j * np.pi * (t_au - tau) /TX_au) )

fp_t = lambda tau: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (t_au - tau)/TX_au)
                                     + np.exp(-2j*np.pi* (t_au - tau) /TX_au) )

FX_t = lambda tau: - A0X * np.cos(Omega_au * (t_au - tau)) * fp_t(tau) + A0X * Omega_au * np.sin(Omega_au * (t_au - tau)) * f_t(tau)
#Variante mit TX
f_TX = lambda tau: 1./4 * ( np.exp(2j * np.pi * (TX_au/2 - tau) / TX_au)
                      + 2
                      + np.exp(-2j * np.pi * (TX_au/2 - tau) /TX_au) )

fp_TX = lambda tau: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (TX_au/2 - tau)/TX_au)
                                     + np.exp(-2j*np.pi* (TX_au/2 - tau) /TX_au) )

FX_TX = lambda tau: - A0X * np.cos(Omega_au * (TX_au/2 - tau)) * fp_TX(tau) + A0X * Omega_au * np.sin(Omega_au * (TX_au/2 - tau)) * f_TX(tau)

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) * omega_au / TL_au
                               + phi)**2
integ_IR = lambda t3: (p_au + A_IR(t3))**2

#-------------------------------------------------------------------------
# technical defintions of functions

# probiere Umschreiben der Integrationsvariable
fun_t_1 = lambda tau: np.exp(-tau * res) * FX_t(tau)
fun_t_2 = lambda tau: np.exp(complex(0,E_kin_au) * tau) * FX_t(tau)

fun_TX2_1 = lambda tau: np.exp(-tau * res) * FX_TX(tau)
fun_TX2_2 = lambda tau: np.exp(complex(0,E_kin_au) * tau) * FX_TX(tau)


#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2

outfile.write(' '.join(('TX/2                 = ',
                        str(sciconv.atu_to_second(TX_au/2)), 's', '\n')))
outfile.write(' '.join(('TL/2                 = ',
                        str(sciconv.atu_to_second(TL_au/2)), 's', '\n')))
outfile.write(' '.join(('delta_t_au - TL_au/2 = ',
                        str(sciconv.atu_to_second(delta_t_au - TL_au/2)), 's', '\n')))
outfile.write(' '.join(('delta_t_au + TL_au/2 = ',
                        str(sciconv.atu_to_second(delta_t_au + TL_au/2)), 's', '\n')))
outfile.write(' '.join(('tmax                 = ',
                        str(sciconv.atu_to_second(tmax_au)), 's', '\n')))

#-------------------------------------------------------------------------
# constant integrals, they are independent of both Omega and t
integral_6_12 = aires.integral_6_12(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin)
res_integral_6_12 = integral_6_12 * (-VEr_au) * rdg_au
indir_integral_6_12 = integral_6_12 * 1j * np.pi / VEr_au**2 * cdg_au

integral_7_13 = aires.integral_7_13(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin)
res_integral_7_13 = integral_7_13 * (-VEr_au) * rdg_au
indir_integral_7_13 = integral_7_13 * 1j * np.pi / VEr_au**2 * cdg_au

integral_14 = aires.integral_14(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin)
res_integral_14 = integral_14 * (-VEr_au) * rdg_au
indir_integral_14 = integral_14 * 1j * np.pi / VEr_au**2 * cdg_au

integral_15 = aires.integral_15(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin)
res_integral_15 = integral_15 * (-VEr_au) * rdg_au
indir_integral_15 = integral_15 * 1j * np.pi / VEr_au**2 * cdg_au

integral_16 = aires.integral_16(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin)
res_integral_16 = integral_16 * (-VEr_au) * rdg_au
indir_integral_16 = integral_16 * 1j * np.pi / VEr_au**2 * cdg_au

# direct ionization
dir_integral_5_8 = aidir.integral_5_8(cdg=cdg_au, E_kin=E_kin_au,
                                      TX=TX_au, TL=TL_au, delta=delta_t_au, t=t_au)
dir_integral_9 = aidir.integral_9(cdg=cdg_au, E_kin=E_kin_au,
                                      TX=TX_au, TL=TL_au, delta=delta_t_au)

# sums of constant terms
res_const_after = (res_integral_6_12 + res_integral_7_13 + res_integral_14
                   + res_integral_15)
indir_const_after = (indir_integral_6_12 + indir_integral_7_13 + indir_integral_14
                   + indir_integral_15)



#-------------------------------------------------------------------------
while ((t_au <= TX_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the first pulse \n')

    outlines = []
    Omega_au = Omega_min_au
    
    while (Omega_au < Omega_max_au):

# integral 1
# other integration variable
        I1 = ci.complex_quadrature(fun_t_1, (t_au + TX_au/2), 0)
        I2 = ci.complex_quadrature(fun_t_2, (t_au + TX_au/2), 0)

        res_J = - rdg_au * VEr_au / res_kin * (I1[0] - I2[0])
        indir_J = 1j * np.pi / VEr_au**2 * cdg_au / res_kin * (I1[0] - I2[0])
        dir_J = 1j * cdg_au * I2[0]

        J = res_J + indir_J + dir_J

        string = in_out.prep_output(J, Omega_au, t_au)
        outlines.append(string)
        
        Omega_au = Omega_au + Omega_step_au
    
    
    in_out.doout_1f(pure_out, outlines)

    t_au = t_au + timestep_au




#-------------------------------------------------------------------------
while (t_au >= TX_au/2 and t_au <= (delta_t_au - TL_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('between the pulses \n')

    Omega_au = Omega_min_au
    outlines = []

    # integrals 3 and 4 are independent of omega, they are therefore
    # evaluated before integral 2 and especially outside the loop
    #integral 3
    integral_3 = aires.integral_3(VEr_au, rdg_au, E_kin_au, TX_au, res, res_kin, t_au)
    res_integral_3 = integral_3 * (-VEr_au) * rdg_au
    indir_integral_3 = integral_3 * 1j * np.pi / VEr_au**2 * cdg_au
    dir_integral_3 = aidir.integral_3(cdg=cdg_au, E_kin=E_kin_au, TX=TX_au, t=t_au)

    K = res_integral_3 + indir_integral_3 + dir_integral_3

    #integral 4
    integral_4 = aires.integral_3(VEr_au, rdg_au, E_kin_au, TX_au, res, res_kin, t_au)
    res_integral_4 = integral_4 * (-VEr_au) * rdg_au
    indir_integral_4 = integral_4 * 1j * np.pi / VEr_au**2 * cdg_au
    K = K + res_integral_4 + indir_integral_4
    
    
    while (Omega_au < Omega_max_au):
        # integral 2

# other integration variable
        I1 = ci.complex_quadrature(fun_TX2_1, (TX_au/2 + TX_au/2), 0)
        I2 = ci.complex_quadrature(fun_TX2_2, (TX_au/2 + TX_au/2), 0)

        res_J = - rdg_au * VEr_au / res_kin * (I1[0] - I2[0])
        indir_J = 1j * np.pi / VEr_au**2 * cdg_au / res_kin * (I1[0] - I2[0])
        dir_J = 1j * cdg_au * I2[0]

        J = res_J + indir_J + dir_J

        L = K + J

        string = in_out.prep_output(L, Omega_au, t_au)
        outlines.append(string)
        
        Omega_au = Omega_au + Omega_step_au
    
    
    in_out.doout_1f(pure_out,outlines)

    t_au = t_au + timestep_au



#-------------------------------------------------------------------------
# during the ir pulse
while (t_au >= (delta_t_au - TL_au/2)
       and t_au <= (delta_t_au + TL_au/2)
       and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the second pulse \n')
    # integrals, that are independent of omega
    integral_8 = aires.integral_8(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                               TX=TX_au, TL=TL_au, delta=delta_t_au,
                               res=res, res_kin=res_kin, t=t_au)
    res_integral_8 = integral_8 * (-VEr_au) * rdg_au
    indir_integral_8 = integral_8 * 1j * np.pi / VEr_au**2 * cdg_au

    integral_9 = aires.integral_9(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                               TX=TX_au, TL=TL_au, delta=delta_t_au,
                               res=res, res_kin=res_kin, t=t_au)
    res_integral_9 = integral_9 * (-VEr_au) * rdg_au
    indir_integral_9 = integral_9 * 1j * np.pi / VEr_au**2 * cdg_au

    integral_10 = aires.integral_10(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    res_integral_10 = integral_10 * (-VEr_au) * rdg_au
    indir_integral_10 = integral_10 * 1j * np.pi / VEr_au**2 * cdg_au

    dir_integral_6 = aidir.integral_6(cdg=cdg_au, E_kin=E_kin_au, TX=TX_au, TL=TL_au,
                                      delta=delta_t_au, t=t_au)

    I_IR = integrate.quad(integ_IR, delta_t_au - TL_au/2, t_au)
    res_I10  = res_integral_10 * I_IR[0]
    indir_I10  = indir_integral_10 * I_IR[0]
    dir_I6 = dir_integral_6 * I_IR[0]

    K = (res_integral_8 + res_integral_9 + res_I10 + res_integral_7_13
         + indir_integral_8 + indir_integral_9 + indir_I10 + indir_integral_7_13
         + dir_integral_5_8 + dir_I6)

    Omega_au = Omega_min_au
    outlines = []
    
    while (Omega_au < Omega_max_au):
        # integral 5 = integral 2
# other integration variable
        I1 = ci.complex_quadrature(fun_TX2_1, (TX_au/2 + TX_au/2), 0)
        I2 = ci.complex_quadrature(fun_TX2_2, (TX_au/2 + TX_au/2), 0)

        res_J = - rdg_au * VEr_au / res_kin * (I1[0] - I2[0])
        indir_J = 1j * np.pi / VEr_au**2 * cdg_au / res_kin * (I1[0] - I2[0])
        dir_J = 1j * cdg_au * I2[0]

        J = res_J + indir_J + dir_J

        L = J + K 
        

        string = in_out.prep_output(L, Omega_au, t_au)
        outlines.append(string)
        
        Omega_au = Omega_au + Omega_step_au
    
    
    in_out.doout_1f(pure_out,outlines)

    t_au = t_au + timestep_au



#-------------------------------------------------------------------------
# after the second pulse
while (t_au >= (delta_t_au + TL_au/2)
       and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('after the second pulse')

    # omega independent integrals
    #integral 16
    I_IR = integrate.quad(integ_IR, delta_t_au - TL_au/2, delta_t_au + TL_au/2)
    res_integral_16_p = integral_16 * I_IR[0]
    indir_integral_16_p = integral_16 * I_IR[0]
    #integral 17
    integral_17 = aires.integral_17(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    res_integral_17 = integral_17 * (-VEr_au) * rdg_au
    indir_integral_17 = integral_17 * 1j * np.pi / VEr_au**2 * cdg_au
    #integral 18
    integral_18 = aires.integral_18(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    res_integral_18 = integral_18 * (-VEr_au) * rdg_au
    indir_integral_18 = integral_18 * 1j * np.pi / VEr_au**2 * cdg_au
    #integral 19
    integral_19 = aires.integral_19(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    res_integral_19 = integral_19 * (-VEr_au) * rdg_au
    indir_integral_19 = integral_19 * 1j * np.pi / VEr_au**2 * cdg_au
    res_integral_19_p = res_integral_19 * I_IR[0]
    indir_integral_19_p = indir_integral_19 * I_IR[0]
    #integral 20
    integral_20 = aires.integral_20(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    res_integral_20 = integral_20 * (-VEr_au) * rdg_au
    indir_integral_20 = integral_20 * 1j * np.pi / VEr_au**2 * cdg_au
    res_integral_20_p = res_integral_20 * I_IR[0]
    indir_integral_20_p = indir_integral_20 * I_IR[0]

    dir_integral_9_p = dir_integral_9 * I_IR[0]
    dir_integral_10 = aidir.integral_10(cdg=cdg_au, E_kin=E_kin_au, TX=TX_au, TL=TL_au,
                                        delta=delta_t_au, t=t_au)
    dir_integral_10_p = dir_integral_10 * I_IR[0]

    K = (res_integral_16_p + res_integral_17 + res_integral_18 + res_integral_19_p
         + res_integral_20_p + res_const_after
         + indir_integral_16_p + indir_integral_17 + indir_integral_18
         + indir_integral_19_p + indir_integral_20_p + indir_const_after
         + dir_integral_5_8 + dir_integral_9_p)

    Omega_au = Omega_min_au
    outlines = []
    
    while (Omega_au < Omega_max_au):
        # integral 11 = integral 5 = integral 2
# other integration variable
        I1 = ci.complex_quadrature(fun_TX2_1, (TX_au/2 + TX_au/2), 0)
        I2 = ci.complex_quadrature(fun_TX2_2, (TX_au/2 + TX_au/2), 0)

        res_J = - rdg_au * VEr_au / res_kin * (I1[0] - I2[0])
        indir_J = 1j * np.pi / VEr_au**2 * cdg_au / res_kin * (I1[0] - I2[0])
        dir_J = 1j * cdg_au * I2[0]

        J = res_J + indir_J + dir_J

        L = J + K


        string = in_out.prep_output(L, Omega_au, t_au)
        outlines.append(string)
        
        Omega_au = Omega_au + Omega_step_au
    
    
    in_out.doout(pure_out,outlines)

    t_au = t_au + timestep_au

outfile.close
pure_out.close
