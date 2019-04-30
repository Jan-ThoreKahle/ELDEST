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
from scipy.signal import argrelextrema
from scipy.special import erf
import numpy as np
import sciconv
import complex_integration as ci
import pulses
import in_out
import sys

infile = sys.argv[1]
print infile

#-------------------------------------------------------------------------
# open outputfile
outfile = open("eldest.out", mode='w')
pure_out = open('full.dat', mode='w')

outfile.write("The results were obtained with dblexc.py \n")
#-------------------------------------------------------------------------
# set some defaults
Xshape = 'convoluted'


#-------------------------------------------------------------------------
# read inputfile
(rdg_au, cdg_au,
 Er_a_eV, Er_b_eV, tau_a_s, tau_b_s, E_fin_eV, tau_s, E_fin_eV_2, tau_s_2,
 interact_eV,
 Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
 omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q, sigma_L,
 tmax_s, timestep_s, E_step_eV,
 E_min_eV, E_max_eV,
 integ, integ_outer) = in_out.read_input(infile, outfile)


#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_a_au        = sciconv.ev_to_hartree(Er_a_eV)
Er_b_au        = sciconv.ev_to_hartree(Er_a_eV)
Er_au          = Er_a_au + Er_b_au
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)

tau_au         = sciconv.second_to_atu(tau_s)
Gamma_au       = 1. / tau_au
Gamma_eV       = sciconv.hartree_to_ev(Gamma_au)
print "Gamma_eV = ", Gamma_eV
outfile.write('Gamma_eV = ' + str(Gamma_eV) + '\n')

# laser parameters
Omega_au      = sciconv.ev_to_hartree(Omega_eV)
if (X_sinsq):
    TX_au     = n_X * 2 * np.pi / Omega_au
elif(X_gauss):
    sigma     = np.pi * n_X / (Omega_au * np.sqrt(np.log(2)))
    FWHM      = 2 * np.sqrt( 2 * np.log(2)) * sigma
    TX_au     = 5 * sigma
    print 'sigma = ', sciconv.atu_to_second(sigma)
    print 'FWHM = ', sciconv.atu_to_second(FWHM)
    outfile.write('sigma = ' + str(sciconv.atu_to_second(sigma)) + '\n')
    outfile.write('FWHM = ' + str(sciconv.atu_to_second(FWHM)) + '\n')
print 'end of the first pulse = ', sciconv.atu_to_second(TX_au)
outfile.write('end of the first pulse = ' + str(sciconv.atu_to_second(TX_au)) + '\n')
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


##-------------------------------------------------------------------------
#in_out.check_input(Er_au, E_fin_au, Gamma_au,
#                   Omega_au, TX_au, n_X, A0X,
#                   omega_au, TL_au, A0L, delta_t_au,
#                   tmax_au, timestep_au, E_step_au)
#-------------------------------------------------------------------------
# physical defintions of functions
# functions for the shape of the XUV pulse
# use RPA and assume very slowly varying pulse shape
if (X_sinsq):
    print 'use sinsq function'
    f_t1  = lambda t1: 1./4 * ( np.exp(2j * np.pi * (t1 + TX_au/2) / TX_au)
                          + 2
                          + np.exp(-2j * np.pi * (t1 + TX_au/2) /TX_au) )

    fp_t1 = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (t1 + TX_au/2) / TX_au)
                                         + np.exp(-2j*np.pi* (t1 + TX_au/2) / TX_au) )
elif (X_gauss):
    print 'use gauss function'
    f_t1  = lambda t1: ( 1./ np.sqrt(2*np.pi * sigma**2)
                       * np.exp(-t1**2 / (2*sigma**2)))
    fp_t1 = lambda t1: ( -t1 / np.sqrt(2*np.pi) / sigma**3
                       * np.exp(-t1**2 / (2*sigma**2)))
else:
    print 'no pulse shape selected'

if (Xshape == 'convoluted'):
    FX_t1 = lambda t1: (
                        0
                        - (A0X * 2j
                           * Omega_au
                           * np.exp(-1j*Omega_au * (t1))
                           * f_t1(t1)
                          )
                       )
elif (Xshape == 'infinite'):
    FX_t1 = lambda t1: + A0X * Omega_au * np.cos(Omega_au * t1)
    #FX_t1 = lambda t1: - A0X * np.sin(Omega_au * t1)
                       

#-------------------------------------------------------------------------
# technical defintions of functions

#direct ionization
fun_t_dir = lambda t1: - np.exp(-1j*t_au * (E_kin_au + E_fin_au)) \
                       * np.exp(-sigma**2/2 * (Er_a_au - Omega_au)**2) \
                       * FX_t1(t1) * np.exp(-1j * t1 * (Er_a_au - E_kin_au - E_fin_au)) \
                       * erf(1./sigma/np.sqrt(2) * (t1 - 1j * sigma**2 * (Er_a_au - Omega_au)))
fun_TX2_dir = lambda t1: - np.exp(-1j*t_au * (E_kin_au + E_fin_au)) \
                       * np.exp(-sigma**2/2 * (Er_a_au - Omega_au)**2) \
                       * FX_t1(t1) * np.exp(-1j * t1 * (Er_a_au - E_kin_au - E_fin_au)) \
                       * erf(1./sigma/np.sqrt(2) * (t1 - 1j * sigma**2 * (Er_a_au - Omega_au)))



res_inner_fun = lambda t2: np.exp(-t2 * (np.pi * VEr_au**2 + 1j*(Er_au))) \
                           * IR_during(t2)

if (integ == 'romberg'):
    res_inner = lambda t1: integrate.romberg(res_inner_fun, t1, t_au)
elif (integ == 'quadrature'):
    res_inner = lambda t1: integrate.quad(res_inner_fun, t1, t_au)[0]
elif (integ == 'analytic'):
# analytic inner integral
    res_inner = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au) - np.pi * VEr_au**2)
                            * (np.exp(t_au * (1j*(E_kin_au + E_fin_au - Er_au) - np.pi * VEr_au**2))
                              - np.exp(t1 * (1j*(E_kin_au + E_fin_au - Er_au) - np.pi * VEr_au**2)))
                            * np.exp(-1j*t_au * (E_kin_au + E_fin_au))
                           )

res_outer_fun = lambda t1: FX_t1(t1) * np.exp(t1 * (np.pi* VEr_au**2 + 1j*Er_au)) \
                           * res_inner(t1)

res_int1 = lambda t1: FX_t1(t1) * np.exp(t1 * (np.pi * VEr_au**2 + 1j* (Er_au-Er_a_au)))
res_int2 = lambda t1: FX_t1(t1) * np.exp(t1 * (np.pi * VEr_au**2 + 1j* (Er_au-Er_a_au))) \
                      * erf(1./np.sqrt(2)/sigma * (t1 - 1j*sigma**2 * (Er_a_au - Omega_au)))
res_int3 = lambda t1: FX_t1(t1) * np.exp(-t1 * 1j* (Er_a_au - E_kin_au - E_fin_au))
res_int4 = lambda t1: FX_t1(t1) * np.exp(-t1 * 1j* (Er_a_au - E_kin_au - E_fin_au)) \
                      * erf(1./np.sqrt(2)/sigma * (t1 - 1j*sigma**2 * (Er_a_au - Omega_au)))

#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2

# construct list of energy points
Ekins = []
E_kin_au = E_min_au
while (E_kin_au <= E_max_au):
    Ekins.append(sciconv.hartree_to_ev(E_kin_au))
    E_kin_au = E_kin_au + E_step_au


#-------------------------------------------------------------------------
# constants / prefactors
# jdg is assumed to be 1
prefac_res = -VEr_au * rdg_au * A0X * Omega_au/ 2
prefac_indir = 1j * np.pi * VEr_au**2 * cdg_au * A0X * Omega_au/ 2
prefac_dir = 1j * cdg_au * A0X * Omega_au/ 2


#-------------------------------------------------------------------------
while ((t_au <= TX_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the first pulse \n')
    print 'during the first pulse'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')

    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)

        dir_const = np.exp(-1j* t_au * (E_kin_au + E_fin_au)) \
                    * np.exp(-sigma**2/2 * (Er_a_au - Omega_au)**2) \
                    * erf(1./np.sqrt(2)/sigma * (t_au - 1j*sigma**2 * (Er_a_au - Omega_au)) )
        res_p1 = np.exp(-t_au * (1j*(Er_au - E_fin_au - E_kin_au) + np.pi * VEr_au**2)) \
                 * erf(1./np.sqrt(2)/sigma * (t_au - 1j*sigma**2 * (Er_a_au - Omega_au)))
        res_p2 = - np.exp(-t_au * (1j*(Er_au - E_fin_au - E_kin_au) + np.pi * VEr_au**2))
        res_p3 = - erf(1./np.sqrt(2)/sigma * (t_au - 1j*sigma**2 * (Er_a_au - Omega_au)))
        res_p4 = 1
# integral 1
        if (integ_outer == "quadrature"):
            I = ci.complex_quadrature(fun_t_dir, (-TX_au/2), t_au)
            res_1 = ci.complex_quadrature(res_int1, (-TX_au/2), t_au)
            res_2 = ci.complex_quadrature(res_int2, (-TX_au/2), t_au)
            res_3 = ci.complex_quadrature(res_int3, (-TX_au/2), t_au)
            res_4 = ci.complex_quadrature(res_int4, (-TX_au/2), t_au)

            dir_J = prefac_dir * (I[0] + dir_const)
            res_J = prefac_res * (
                      res_1[0] * res_p1
                    + res_2[0] * res_p2
                    + res_3[0] * res_p3
                    + res_4[0] * res_p4
                    )
            indir_J = prefac_indir * (
                      res_1[0] * res_p1
                    + res_2[0] * res_p2
                    + res_3[0] * res_p3
                    + res_4[0] * res_p4
                    )

        elif (integ_outer == "romberg"):
            I = ci.complex_romberg(fun_t_dir, (-TX_au/2), t_au)
            #res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), t_au)
            res_1 = ci.complex_romberg(res_int1, (-TX_au/2), t_au)
            res_2 = ci.complex_romberg(res_int2, (-TX_au/2), t_au)
            res_3 = ci.complex_romberg(res_int3, (-TX_au/2), t_au)
            res_4 = ci.complex_romberg(res_int4, (-TX_au/2), t_au)

            dir_J = prefac_dir * (I + dir_const)
            res_J = prefac_res * (
                      res_1 * res_p1
                    + res_2 * res_p2
                    + res_3 * res_p3
                    + res_4 * res_p4
                    )
            indir_J = prefac_indir * (
                      res_1 * res_p1
                    + res_2 * res_p2
                    + res_3 * res_p3
                    + res_4 * res_p4
                    )

        J = (0
             + dir_J
             + res_J
             + indir_J
             )

        square = np.absolute(J)**2
        dir_term = np.absolute(dir_J)**2
        squares = np.append(squares, square)

        string = in_out.prep_output_comp(square, dir_term, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au
    
    
    in_out.doout_1f(pure_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print Ekins[max_pos[i]], squares[max_pos[i]]
            outfile.write(str(Ekins[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')
    

    t_au = t_au + timestep_au




#-------------------------------------------------------------------------
while (t_au >= TX_au/2 and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('between the pulses \n')
    print 'between the pulses'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')
    while (E_kin_au <= E_max_au):

        dir_const = np.exp(-1j* TX_au/2 * (E_kin_au + E_fin_au)) \
                    * np.exp(-sigma**2/2 * (Er_a_au - Omega_au)**2) \
                    * erf(1./np.sqrt(2)/sigma * (TX_au/2 - 1j*sigma**2 * (Er_a_au - Omega_au)) )
        res_p1 = np.exp(-t_au * (1j*(Er_au - E_fin_au - E_kin_au) + np.pi * VEr_au**2)) \
                 * erf(1./np.sqrt(2)/sigma * (TX_au/2 - 1j*sigma**2 * (Er_a_au - Omega_au)))
        res_p2 = - np.exp(-t_au * (1j*(Er_au - E_fin_au - E_kin_au) + np.pi * VEr_au**2))
        res_p3 = - erf(1./np.sqrt(2)/sigma * (TX_au/2 - 1j*sigma**2 * (Er_a_au - Omega_au)))
        res_p4 = 1

# integral 1
        if (integ_outer == "quadrature"):
            I1 = ci.complex_quadrature(fun_TX2_dir, (-TX_au/2), TX_au/2)
            #res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)
            res_1 = ci.complex_quadrature(res_int1, (-TX_au/2), TX_au/2)
            res_2 = ci.complex_quadrature(res_int2, (-TX_au/2), TX_au/2)
            res_3 = ci.complex_quadrature(res_int3, (-TX_au/2), TX_au/2)
            res_4 = ci.complex_quadrature(res_int4, (-TX_au/2), TX_au/2)
 
            dir_J = prefac_dir * (I1[0] + dir_const)
            res_J = prefac_res * (
                      res_1[0] * res_p1
                    + res_2[0] * res_p2
                    + res_3[0] * res_p3
                    + res_4[0] * res_p4
                    )
            indir_J = prefac_indir * (
                      res_1[0] * res_p1
                    + res_2[0] * res_p2
                    + res_3[0] * res_p3
                    + res_4[0] * res_p4
                    )
        
        elif (integ_outer == "romberg"):
            I1 = ci.complex_romberg(fun_TX2_dir, (-TX_au/2), TX_au/2)
            res_1 = ci.complex_romberg(res_int1, (-TX_au/2), TX_au/2)
    
            dir_J = prefac_dir * (I1 + dir_const)
            res_J = prefac_res * (
                      res_1 * res_p1
                    + res_2 * res_p2
                    + res_3 * res_p3
                    + res_4 * res_p4
                    )
            indir_J = prefac_indir * (
                      res_1 * res_p1
                    + res_2 * res_p2
                    + res_3 * res_p3
                    + res_4 * res_p4
                    )

        J = (0
             #+ dir_J
             + res_J
             + indir_J
             )

        square = np.absolute(J)**2
        dir_term = np.absolute(dir_J)**2
        squares = np.append(squares, square)

        string = in_out.prep_output_comp(square, dir_term, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au

    
    
    in_out.doout_1f(pure_out,outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print Ekins[max_pos[i]], squares[max_pos[i]]
            outfile.write(str(Ekins[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')

    t_au = t_au + timestep_au



outfile.close
pure_out.close

print "end of script"
sys.exit(0)
