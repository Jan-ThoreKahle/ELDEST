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

rdg_au        = 0.30          # transition dipole moment into the resonant state
cdg_au        = 0.5           # transition dipole moment into any continuum state

# parameters of the investigated system
# the ground state energy is being defined as Eg = 0
Er_eV         = 150.0         # resonance energy in eV
E_fin_eV      =  70.0         # final state energy in eV

tau_s         =  400.0E-18       # lifetime

# laser parameters
Omega_eV      = 150.0          #
TX_s          = 250.0E-18       # duration of the XUV pulse in seconds
n_X           = 10
I_X           = 5.0E16        # intensity of the XUV pulse in W/cm^2

omega_eV      = 1.6           # IR pulse
n_L           = 10
I_L           = 5.0E12        # intensity of the IR pulse in W/cm^2
delta_t_s     = 15000.0E-18       # time difference between the maxima of the two pulses
phi           = 0
q             = 5

# parameters of the simulation
tmax_s        = 8.0E-14       # simulate until time tmax in seconds
timestep_s    = 500.0E-18     # evaluate expression every timestep_s seconds 
E_step_eV     = 3.0           # energy difference between different evaluated Omegas

E_min_eV      =  50.0
E_max_eV      = 110.0
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_au          = sciconv.ev_to_hartree(Er_eV)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)

tau_au         = sciconv.second_to_atu(tau_s)
Gamma_au       = 1. / tau_au

# laser parameters
Omega_au      = sciconv.ev_to_hartree(Omega_eV)
TX_au         = sciconv.second_to_atu(TX_s)
TX_au         = n_X * 2 * np.pi / Omega_au
print 'end of the first pulse = ', sciconv.atu_to_second(TX_au)
I_X_au        = sciconv.Wcm2_to_aiu(I_X)
print 'I_X = ', I_X
print 'I_X_au = ', I_X_au
E0X           = np.sqrt(I_X_au)
A0X           = E0X / Omega_au
print 'A0X = ', A0X

omega_au      = sciconv.ev_to_hartree(omega_eV)
TL_au         = n_L * 2 * np.pi / omega_au
print 'start of IR pulse = ', delta_t_s - sciconv.atu_to_second(TL_au/2)
print 'end of IR pulse = ', delta_t_s + sciconv.atu_to_second(TL_au/2)
I_L_au        = sciconv.Wcm2_to_aiu(I_L)
print 'I_L = ', I_L
print 'I_L_au = ', I_L_au
E0L           = np.sqrt(I_L_au)
print 'E0L', E0L
A0L           = E0L / omega_au
print 'A0L = ', A0L
delta_t_au    = sciconv.second_to_atu(delta_t_s)

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
E_step_au = sciconv.ev_to_hartree(E_step_eV)

E_min_au = sciconv.ev_to_hartree(E_min_eV)
E_max_au = sciconv.ev_to_hartree(E_max_eV)

VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))
print 'VEr_au = ', VEr_au

#test q=1
cdg_au = rdg_au / ( q * np.pi * VEr_au)
print 'cdg_au = ', cdg_au


#-------------------------------------------------------------------------
in_out.check_input(Er_au, E_fin_au, Gamma_au,
                   Omega_au, TX_au, n_X, A0X,
                   omega_au, TL_au, A0L, delta_t_au,
                   tmax_au, timestep_au, E_step_au)
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

# functions for the norm
f_t1  = lambda t1: 1./4 * ( np.exp(2j * np.pi * t1 / TX_au)
                      + 2
                      + np.exp(-2j * np.pi * t1 /TX_au) )

fp_t1 = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* t1/TX_au)
                                     + np.exp(-2j*np.pi* t1 /TX_au) )

FX_t1 = lambda t1: - A0X * np.cos(Omega_au * t1) * fp_t(t1) + A0X * Omega_au * np.sin(Omega_au * (t1)) * f_t(t1)

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) / TL_au)**2 \
                      * np.cos(omega_au * t3 + phi)
integ_IR = lambda t3: (p_au + A_IR(t3))**2
integ_IR_off = lambda t3, t1: (p_au + A_IR(t3) - A_IR(t1))**2

#-------------------------------------------------------------------------
# technical defintions of functions

# probiere Umschreiben der Integrationsvariable
fun_t_1 = lambda tau: np.exp(-tau * res) * FX_t(tau)
fun_t_2 = lambda tau: np.exp(complex(0,E_kin_au) * tau) * FX_t(tau)

fun_TX2_1 = lambda tau: np.exp(-tau * res) * FX_TX(tau)
fun_TX2_2 = lambda tau: np.exp(complex(0,E_kin_au) * tau) * FX_TX(tau)

#direct ionization
#fun_t_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
#                                   * np.exp(1j * E_kin_au * (t1-t_au))
fun_t_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                   * np.exp(1j * (p_au - A_IR(t1))**2
                                               * (t1-t_au) / 2)
#fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
#                                   * np.exp(1j * E_kin_au * (t1-TX_au/2))
fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                     * np.exp(1j * (p_au - A_IR(t1))**2
                                     * (t1-TX_au/2) / 2)

dress_I = lambda t1: integrate.quad(integ_IR,t1,t_au)[0]
dress_I1 = lambda t1: integrate.quad(integ_IR_off,delta_t_au - TL_au/2, t_au,
                                     args=(t1,))[0]
dress = lambda t1: np.exp(-1j/2 * dress_I(t1))
dress = lambda t1: np.exp(-1j/2 * dress_I1(t1))

#dress_I_after = lambda t1: integrate.quad(integ_IR,t1,(delta_t_au + TL_au/2))[0]
dress_I_after = lambda t1: integrate.quad(integ_IR_off,(delta_t_au-TL_au/2),
                                         (delta_t_au + TL_au/2), args=(t1,))[0]
dress_after = lambda t1: np.exp(-1j/2 * dress_I_after(t1))
#fun_dress_after = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
#                              * np.exp(1j * E_kin_au * ((delta_t_au + TL_au/2)-t_au)) \
#                              * dress_after(t1)
fun_dress_after = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                  * np.exp(1j * (p_au - A_IR(t1))**2
                                  * (t1 + TL_au/2 - t_au) / 2) \
                                  * dress_after(t1)

#fun_IR_dir = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
#                                  * dress(t1)
fun_IR_dir = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                  * np.exp(1j * (p_au - A_IR(t1))**2
                                  * (t1 - delta_t_au + TL_au/2) / 2) \
                                  * dress(t1)
                                  #* np.exp(-1j/2 *
                                  #  integrate.quad((p_au+A_IR-A_IR(t1))**2,
                                  #  delta_t_au-TL_au/2, t_au))

#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2


#-------------------------------------------------------------------------
# constants / prefactors
#res_kin = complex(Gamma_au/2,Er_au + E_kin_au)
res     = complex(Gamma_au/2,Er_au)
print 'res = ', res

prefac_res = - VEr_au * rdg_au
prefac_indir = 1j * np.pi * VEr_au**2 * cdg_au
#prefac_indir = 0
prefac_dir = 1j * cdg_au

print 'prefac_res', prefac_res
print 'prefac_indir', prefac_indir


#-------------------------------------------------------------------------
while ((t_au <= TX_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the first pulse \n')
    print 'during the first pulse'

    outlines = []
    E_kin_au = E_min_au
    
    print 't_s = ', sciconv.atu_to_second(t_au)
    while (E_kin_au <= E_max_au):

        p_au = np.sqrt(2*E_kin_au)

# integral 1
        I = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), t_au)

        dir_J = prefac_dir * I[0]

        J = (0
             + dir_J
             )

        square = np.absolute(J)**2

        string = in_out.prep_output(square, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au
    
    
    in_out.doout_1f(pure_out, outlines)

    t_au = t_au + timestep_au




#-------------------------------------------------------------------------
while (t_au >= TX_au/2 and t_au <= (delta_t_au - TL_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('between the pulses \n')
    print 'between the pulses'

    outlines = []
    E_kin_au = E_min_au
    
    print 't_s = ', sciconv.atu_to_second(t_au)
    #while (E_kin_au <= Omega_au):
    while (E_kin_au <= E_max_au):

        p_au = np.sqrt(2*E_kin_au)

# integral 1
        I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)

        dir_J = prefac_dir * (I1[0]
                              )

        J = (0
             + dir_J
             )

        square = np.absolute(J)**2

        string = in_out.prep_output(square, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au

    
    
    in_out.doout_1f(pure_out,outlines)

    t_au = t_au + timestep_au


#-------------------------------------------------------------------------
# during the ir pulse
while (t_au >= (delta_t_au - TL_au/2)
       and t_au <= (delta_t_au + TL_au/2)
       and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the IR pulse \n')
    print 'during the IR pulse'

    outlines = []
    E_kin_au = E_min_au
    
    print 't_s = ', sciconv.atu_to_second(t_au)
    #while (E_kin_au <= Omega_au):
    while (E_kin_au <= E_max_au):

        p_au = np.sqrt(2*E_kin_au)


        I1 = ci.complex_quadrature(fun_IR_dir, (-TX_au/2), TX_au/2)

        dir_J = prefac_dir * (I1[0]
                              )

        J = (0
             + dir_J
             )

        square = np.absolute(J)**2

        string = in_out.prep_output(square, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au

    
    
    in_out.doout_1f(pure_out,outlines)

    t_au = t_au + timestep_au



#-------------------------------------------------------------------------
# after the second pulse
while (t_au >= (delta_t_au + TL_au/2)
       and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('after the IR pulse \n')
    print 'after the IR pulse'

    outlines = []
    E_kin_au = E_min_au
    
    print 't_s = ', sciconv.atu_to_second(t_au)
    #while (E_kin_au <= Omega_au):
    while (E_kin_au <= E_max_au):

        p_au = np.sqrt(2*E_kin_au)

        I1 = ci.complex_quadrature(fun_dress_after, (-TX_au/2), TX_au/2)

        dir_J = prefac_dir * (I1[0]
                              )

        J = (0
             + dir_J
             )

        square = np.absolute(J)**2

        string = in_out.prep_output(square, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au

    
    
    in_out.doout_1f(pure_out,outlines)

    t_au = t_au + timestep_au



outfile.close
pure_out.close
