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
import numpy as np
import sciconv
import complex_integration as ci
import pulses
import in_out
import sys
import warnings


# don't print warnings unless python -W ... is used
if not sys.warnoptions:
    warnings.simplefilter("ignore")

infile = sys.argv[1]
print infile

#-------------------------------------------------------------------------
# open outputfile
outfile = open("eldest.out", mode='w')
pure_out = open('full.dat', mode='w')

outfile.write("The results were obtained with photoelectron.py \n")
#-------------------------------------------------------------------------
# set some defaults
Xshape = 'convoluted'


#-------------------------------------------------------------------------
# read inputfile
(rdg_au, cdg_au,
 Er_a_eV, Er_b_eV, tau_a_s, tau_b_s, E_fin_eV, tau_s, E_fin_eV_2, tau_s_2,
 interact_eV,
 Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
 omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q,
 tmax_s, timestep_s, E_step_eV,
 E_min_eV, E_max_eV,
 integ, integ_outer) = in_out.read_input(infile, outfile)


#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_a_au          = sciconv.ev_to_hartree(Er_a_eV)
Er_b_au          = sciconv.ev_to_hartree(Er_b_eV)
tau_a_au         = sciconv.second_to_atu(tau_a_s)
tau_b_au         = sciconv.second_to_atu(tau_b_s)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)
E_fin_au_1     = sciconv.ev_to_hartree(E_fin_eV)

interact_au    = sciconv.ev_to_hartree(interact_eV)

tau_au         = sciconv.second_to_atu(tau_s)

Gamma_a_au       = 1. / tau_a_au
Gamma_b_au       = 1. / tau_b_au

Gamma_au       = 1. / tau_au
Gamma_eV       = sciconv.hartree_to_ev(Gamma_au)
outfile.write('Gamma_eV = ' + str(Gamma_eV) + '\n')

# second final state
E_fin_au_2       = sciconv.ev_to_hartree(E_fin_eV_2)
tau_au_2         = sciconv.second_to_atu(tau_s_2)
Gamma_au_2       = 1. / tau_au_2

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
print 'end of the first pulse = ', sciconv.atu_to_second(TX_au/2)
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

VEr_a_au      = np.sqrt(Gamma_a_au/ (2*np.pi))
VEr_b_au      = np.sqrt(Gamma_b_au/ (2*np.pi))

#test q=1
cdg_au_V = rdg_au / ( q * np.pi * VEr_a_au)
#cdg_au_W = rdg_au / ( q * np.pi * WEr_au)


#-------------------------------------------------------------------------
in_out.check_input(Er_au, E_fin_au, Gamma_au,
                   Omega_au, TX_au, n_X, A0X,
                   omega_au, TL_au, A0L, delta_t_au,
                   tmax_au, timestep_au, E_step_au)
#-------------------------------------------------------------------------
# physical defintions of functions
# functions for the shape of the XUV pulse
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
                        - (A0X
                           * np.cos(Omega_au * t1)
                           * fp_t1(t1)
                          )
                        + (A0X
                           * Omega_au
                           * np.sin(Omega_au * (t1))
                           * f_t1(t1)
                          )
                       )
elif (Xshape == 'infinite'):
    FX_t1 = lambda t1: + A0X * Omega_au * np.cos(Omega_au * t1)
    #FX_t1 = lambda t1: - A0X * np.sin(Omega_au * t1)
                       

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) / TL_au)**2 \
                      * np.cos(omega_au * t3 + phi)
integ_IR = lambda t3: (p_au + A_IR(t3))**2

IR_during = lambda t2:  np.exp(-1j * (E_kin_au + E_fin_au) * (t_au - t2))# \

IR_after = lambda t2:  np.exp(-1j * E_kin_au * (t_au - t2)) #\

#-------------------------------------------------------------------------
# technical defintions of functions

#direct ionization
fun_t_dir_1 = lambda t1: FX_t1(t1)   * np.exp(1j * E_fin_au * (t1-t_au)) \
                                     * np.exp(1j * E_kin_au * (t1-t_au))
fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
                                     * np.exp(1j * E_kin_au * (t1-t_au))

dress_I = lambda t1: integrate.quad(integ_IR,t1,t_au)[0]
dress = lambda t1: np.exp(-1j/2 * dress_I(t1))

dress_I_after = lambda t1: integrate.quad(integ_IR,t1,(delta_t_au + TL_au/2))[0]
dress_after = lambda t1: np.exp(-1j/2 * dress_I_after(t1))
fun_dress_after = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                              * np.exp(1j * E_kin_au * ((delta_t_au + TL_au/2)-t_au)) \
                              * dress_after(t1)

fun_IR_dir = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                  * dress(t1)



res_inner_fun = lambda t2: np.exp(-t2 * (np.pi * (VEr_au**2 + WEr_au**2) + 1j*(Er_au))) \
                           * IR_during(t2)

if (integ == 'romberg'):
    res_inner = lambda t1: integrate.romberg(res_inner_fun, t1, t_au)
elif (integ == 'quadrature'):
    res_inner = lambda t1: integrate.quad(res_inner_fun, t1, t_au)[0]
elif (integ == 'analytic'):
# analytic inner integral
    res_inner = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au)
                                    - np.pi * (VEr_au**2 + WEr_au**2))
                            * (np.exp(t_au * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  - np.pi * (VEr_au**2 + WEr_au**2)))
                              - np.exp(t1 * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  - np.pi * (VEr_au**2 + WEr_au**2))))
                            * np.exp(-1j*t_au * (E_kin_au + E_fin_au))
                           )
# check formula for res_inner!

res_outer_fun = lambda t1: FX_t1(t1) \
                           * np.exp(t1 * (np.pi* (VEr_au**2 + WEr_au**2) + 1j*Er_au)) \
                           * res_inner(t1)

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
# diagonalized properties
E_plus  = (Er_a_au + Er_b_au) / 2 + np.sqrt( (Er_a_au - Er_b_au)**2/4 + interact_au )
E_minus = (Er_a_au + Er_b_au) / 2 - np.sqrt( (Er_a_au - Er_b_au)**2/4 + interact_au )
# transformation matrix
A_a_plus  = - ((Er_a_au - E_plus + interact_au) /
              (np.sqrt((Er_a_au - E_plus + interact_au)**2
                     + (Er_b_au - E_plus + interact_au)**2) ))
A_b_plus  =  (1 /
             (np.sqrt(1 + (Er_b_au - E_plus + interact_au)**2
                        / (Er_a_au - E_plus + interact_au)**2) ))
A_a_minus = - ((Er_b_au - E_minus + interact_au) /
             (np.sqrt((Er_a_au - E_minus + interact_au)**2
                    + (Er_b_au - E_minus + interact_au)**2) ))
A_b_minus =  (1 /
            (np.sqrt(1 + (Er_b_au - E_minus + interact_au)**2
                       / (Er_a_au - E_minus + interact_au)**2) ))
#diagonalized lifetimes
V_plus  = VEr_a_au * A_a_plus + VEr_b_au * A_b_plus
V_minus = VEr_a_au * A_a_minus + VEr_b_au * A_b_minus

r1 = np.sqrt(((E_plus - E_minus)**2 - (np.pi*V_plus**2 + np.pi*V_minus**2)**2)**2
             + 4*np.pi**2 * (E_plus - E_minus)**2 * (V_plus**2 + V_minus**2)**2)
phi1 = np.arctan(2*np.pi * (E_plus - E_minus) * (V_plus**2 - V_minus**2)
                / ((E_plus - E_minus)**2 - np.pi**2 * (V_plus**2 + V_minus**2)))

# auxiliary energies
E2 = (E_plus + E_minus)/2 + 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
     - np.sqrt(r1)/2 * np.cos(phi1/2) - 1j*np.pi/2 np.sin(phi1/2)
E2 = (E_plus + E_minus)/2 - 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
     + np.sqrt(r1)/2 * np.cos(phi1/2) - 1j*np.pi/2 np.sin(phi1/2)
E4 = (E_plus + E_minus)/2 - 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
     - np.sqrt(r1)/2 * np.cos(phi1/2) + 1j*np.pi/2 np.sin(phi1/2)

#aV = VEr_au / np.sqrt(VEr_au**2 + WEr_au**2)
#aW = WEr_au / np.sqrt(VEr_au**2 + WEr_au**2)
#
#prefac_res1 = aV * VEr_au * rdg_au
#prefac_res2 = aW * WEr_au * rdg_au
#prefac_indir1 = -1j * np.pi * VEr_au * (VEr_au + WEr_au) * cdg_au_V
#prefac_indir2 = -1j * np.pi * WEr_au * (VEr_au + WEr_au) * cdg_au_W
##prefac_indir = 0
#prefac_dir1 = 1j * aV * cdg_au_V
#prefac_dir2 = 1j * aW * cdg_au_W


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

# integral 1
        if (integ_outer == "quadrature"):
            E_fin_au = E_fin_au_1

            I1 = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), t_au)
            res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), t_au)

            dir_J1 = prefac_dir1 * I1[0]
            res_J1 = prefac_res1 * res_I[0]
            indir_J1 = prefac_indir1 * res_I[0]

            E_fin_au = E_fin_au_2

            I1 = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), t_au)
            res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), t_au)

            dir_J2   = prefac_dir2 * I1[0]
            res_J2   = prefac_res2 * res_I[0]
            indir_J2 = prefac_indir2 * res_I[0]

        elif (integ_outer == "romberg"):
            E_fin_au = E_fin_au_1

            I1 = ci.complex_romberg(fun_t_dir_1, (-TX_au/2), t_au)
            res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), t_au)
    
            dir_J1 = prefac_dir1 * I1
            res_J1 = prefac_res1 * res_I
            indir_J1 = prefac_indir1 * res_I

            E_fin_au = E_fin_au_2

            I1 = ci.complex_romberg(fun_t_dir_1, (-TX_au/2), t_au)
            res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), t_au)
    
            dir_J2   = prefac_dir2 * I1
            res_J2   = prefac_res2 * res_I
            indir_J2 = prefac_indir2 * res_I

        J = (0
             + dir_J1 + dir_J2
             + res_J1 + res_J2
             + indir_J1 + indir_J2
             )

        square = np.absolute(J)**2
        squares = np.append(squares, square)

        string = in_out.prep_output(square, E_kin_au, t_au)
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
    outfile.write('after the XUV pulse \n')
    print 'after the XUV pulse'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)

# integral 1
        if (integ_outer == "quadrature"):
            E_fin_au = E_fin_au_1

            I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
            res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)

            dir_J1 = prefac_dir1 * I1[0]
            res_J1 = prefac_res1 * res_I[0]
            indir_J1 = prefac_indir1 * res_I[0]

            E_fin_au = E_fin_au_2

            I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
            res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)

            dir_J2   = prefac_dir2 * I1[0]
            res_J2   = prefac_res2 * res_I[0]
            indir_J2 = prefac_indir2 * res_I[0]
        
        elif (integ_outer == "romberg"):
            E_fin_au = E_fin_au_1

            I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
            res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
    
            dir_J1 = prefac_dir1 * I1
            res_J1 = prefac_res1 * res_I
            indir_J1 = prefac_indir1 * res_I

            E_fin_au = E_fin_au_2

            I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
            res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
    
            dir_J2   = prefac_dir2 * I1
            res_J2   = prefac_res2 * res_I
            indir_J2 = prefac_indir2 * res_I

        J = (0
             + dir_J1 + dir_J2
             + res_J1 + res_J2
             + indir_J1 + indir_J2
             )

        square = np.absolute(J)**2
        squares = np.append(squares, square)

        string = in_out.prep_output(square, E_kin_au, t_au)
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

