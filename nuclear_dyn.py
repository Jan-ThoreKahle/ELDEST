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
# written by: Elke Fasshauer November 2020                               #
##########################################################################

import scipy
import scipy.integrate as integrate
from scipy.signal import argrelextrema
from scipy.special import erf
import numpy as np
import sciconv
import complex_integration as ci
#import pulses
import in_out
import sys
import warnings
import wellenfkt as wf


# don't print warnings unless python -W ... is used
if not sys.warnoptions:
    warnings.simplefilter("ignore")

infile = sys.argv[1]
print infile

#-------------------------------------------------------------------------
# open outputfile
outfile = open("eldest.out", mode='w')
pure_out = open('full.dat', mode='w')
movie_out = open('movie.dat', mode='w')
popfile = open("pop.dat", mode='w')

outfile.write("The results were obtained with nuclear_dyn.py \n")
#-------------------------------------------------------------------------
# set some defaults
Xshape = 'convoluted'

#-------------------------------------------------------------------------

(rdg_au, cdg_au,
 Er_a_eV, Er_b_eV, tau_a_s, tau_b_s, E_fin_eV, tau_s, E_fin_eV_2, tau_s_2,
 interact_eV,
 Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
 omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q, FWHM_L,
 tmax_s, timestep_s, E_step_eV,
 E_min_eV, E_max_eV,
 integ, integ_outer,
 mass1, mass2, grad_delta, R_eq_AA,
 gs_de, gs_a, gs_Req, gs_const,
 res_de, res_a, res_Req, res_const,
 fin_a, fin_b, fin_c, fin_d, fin_pot_type
 ) = in_out.read_input(infile, outfile)

#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_a_au        = sciconv.ev_to_hartree(Er_a_eV)
Er_b_au        = sciconv.ev_to_hartree(Er_b_eV)
Er_au          = Er_a_au
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)
E_fin_au_1     = sciconv.ev_to_hartree(E_fin_eV)

tau_au_1       = sciconv.second_to_atu(tau_s)
tau_au         = tau_au_1
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
FWHM_L_au     = sciconv.second_to_atu(FWHM_L)
sigma_L_au    = FWHM_L_au / np.sqrt(8 * np.log(2))
a             = 5./2 * sigma_L_au
print "FWHM_L = ", sciconv.atu_to_second(FWHM_L_au)
print "sigma_L = ", sciconv.atu_to_second(sigma_L_au)
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

VEr_au_1      = VEr_au

#test q=1
cdg_au_V = rdg_au / ( q * np.pi * VEr_au)

#-------------------------------------------------------------------------
# Potential details
# vibrational energies of Morse potentials
print
print '----------------------------------------------------------------'
outfile.write('\n' + "--------------------------------------------------------" + '\n')
red_mass = wf.red_mass_au(mass1,mass2)
print "red_mass = ", red_mass
#ground state
print "Ground state"
lambda_param_gs = np.sqrt(2*red_mass*gs_de) / gs_a
n_gs_max = int(lambda_param_gs - 0.5)
E_kappas = []
print "n_gs_max = ", n_gs_max
for n in range (0,n_gs_max+1):
    ev = wf.eigenvalue(n,gs_de,gs_a,red_mass)
    E_kappas.append(ev)
#resonant state
print "Resonant state"
lambda_param_res = np.sqrt(2*red_mass*res_de) / res_a
n_res_max = int(lambda_param_res - 0.5)
E_lambdas = []
print "n_res_max = ", n_res_max
for n in range (0,n_res_max+1):
    ev = wf.eigenvalue(n,res_de,res_a,red_mass)
    E_lambdas.append(ev)
#final state
print "Final state"
if (fin_pot_type == 'morse'):
    fin_de    = fin_a
    fin_a     = fin_b
    fin_Req   = fin_c
    fin_const = fin_d
    lambda_param_fin = np.sqrt(2*red_mass*fin_de) / fin_a
    n_fin_max = int(lambda_param_fin - 0.5)
    E_mus = []
    print "n_fin_max = ", n_fin_max
    for n in range (0,n_fin_max+1):
        ev = wf.eigenvalue(n,fin_de,fin_a,red_mass)
        E_mus.append(ev)

#-------------------------------------------------------------------------
# Franck-Condon factors
#-------------------------------------------------------------------------
# ground state - resonant state <lambda|kappa>
gs_res =  []
gs_fin =  []
res_fin = []
R_min = sciconv.angstrom_to_bohr(1.5)
R_max = sciconv.angstrom_to_bohr(30.0)
print
print '----------------------------------------------------------------'
print "Franck Condon overlaps between ground and resonant state"
outfile.write('\n' + "--------------------------------------------------------" + '\n')
outfile.write("Franck Condon overlaps between ground and resonant state" + '\n')
outfile.write('n_gs  ' +'n_res  ' + '<res|gs>' + '\n')
for i in range (0,n_gs_max+1):
    tmp = []
    for j in range (0,n_res_max+1):
        FC = wf.FC(j,res_a,res_Req,res_de,red_mass,
                   i,gs_a,gs_Req,gs_de,R_min,R_max)
        #tmp.append(wf.FC(j,res_a,res_Req,res_de,red_mass,
        #                 i,gs_a,gs_Req,gs_de,R_min,R_max))
        tmp.append(FC)
        outfile.write('{:4d}  {:4d}  {:14.10E}\n'.format(i,j,FC))
        print ('{:4d}  {:4d}  {:14.10E}'.format(i,j,FC))
    gs_res.append(tmp)
#print "gs_res"
#print gs_res
    
# ground state - final state <mu|kappa>
print
print '----------------------------------------------------------------'
print "Franck Condon overlaps between ground and final state"
outfile.write('\n' + "--------------------------------------------------------" + '\n')
outfile.write("Franck Condon overlaps between ground and final state" + '\n')
outfile.write('n_gs  ' +'n_fin  ' + '<fin|gs>' + '\n')
if fin_pot_type == 'morse':
    for i in range (0,n_gs_max+1):
        tmp = []
        for j in range (0,n_fin_max+1):
            FC = wf.FC(j,fin_a,fin_Req,fin_de,red_mass,
                       i,gs_a,gs_Req,gs_de,R_min,R_max)
            tmp.append(FC)
            outfile.write('{:4d}  {:4d}  {:14.10E}\n'.format(i,j,FC))
            print ('{:4d}  {:4d}  {:14.10E}'.format(i,j,FC))
        gs_fin.append(tmp)
#    print "gs_fin"
#    print gs_fin

# resonant state - final state <mu|lambdaa>
print
print '----------------------------------------------------------------'
print "Franck Condon overlaps between final and resonant state"
outfile.write('\n' + "--------------------------------------------------------" + '\n')
outfile.write("Franck Condon overlaps between final and resonant state" + '\n')
outfile.write('n_res  ' +'n_fin  ' + '<fin|res>' + '\n')
if fin_pot_type == 'morse':
    for i in range (0,n_res_max+1):
        tmp = []
        for j in range (0,n_fin_max+1):
            FC = wf.FC(j,fin_a,fin_Req,fin_de,red_mass,
                       i,res_a,res_Req,res_de,R_min,R_max)
            tmp.append(FC)
            outfile.write('{:4d}  {:4d}  {:14.10E}\n'.format(i,j,FC))
            print ('{:4d}  {:4d}  {:14.10E}'.format(i,j,FC))
        res_fin.append(tmp)
#    print "res_fin"
#    print res_fin

# sum over mup of product <lambda|mup><mup|kappa>
indir_FCsums = []
for i in range (0,n_res_max+1):
    indir_FCsum = 0
    for j in range (0,n_fin_max+1):
        tmp = np.conj(res_fin[i][j]) * gs_fin[0][j]
        indir_FCsum = indir_FCsum + tmp
    indir_FCsums.append(indir_FCsum)
#print indir_FCsums
print
print '----------------------------------------------------------------'
outfile.write('\n' + "--------------------------------------------------------" + '\n')

#-------------------------------------------------------------------------
# determine total decay width matrix element
print "determine total decay width"
if fin_pot_type == 'morse':
    W_lambda = []
    for i in range (0,n_res_max+1):
        tmp = 0
        for j in range (0,n_fin_max+1):
            tmp = tmp + VEr_au**2 * res_fin[i][j]
        W_lambda.append(tmp)



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
                       

#-------------------------------------------------------------------------
# technical defintions of functions
#direct ionization
fun_t_dir_1 = lambda t1: FX_t1(t1)   * np.exp(1j * E_fin_au * (t1-t_au)) \
                                     * np.exp(1j * E_kin_au * (t1-t_au))
fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
                                     * np.exp(1j * E_kin_au * (t1-t_au))

#res_inner_fun = lambda t2: np.exp(-t2 * (np.pi * W_au + 1j*(Er_au))) \
#                           * IR_during(t2)

if (integ == 'romberg'):
    res_inner = lambda t1: integrate.romberg(res_inner_fun, t1, t_au)
elif (integ == 'quadrature'):
    res_inner = lambda t1: integrate.quad(res_inner_fun, t1, t_au)[0]
elif (integ == 'analytic'):
# analytic inner integral
    res_inner = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au - E_lambda)
                                    - np.pi * W_au)
                            * (np.exp(t_au * (1j*(E_kin_au + E_fin_au
                                                  - Er_au - E_lambda)
                                                  - np.pi * W_au))
                              - np.exp(t1 * (1j*(E_kin_au + E_fin_au
                                                 - Er_au - E_lambda)
                                                  - np.pi * W_au)))
                            * np.exp(-1j*t_au * (E_kin_au + E_fin_au))
                           )

res_outer_fun = lambda t1: FX_t1(t1) \
                           * np.exp(t1 * (np.pi* (W_au) + 1j*(Er_au+E_lambda))) \
                           * res_inner(t1)


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2


Ekins = []
E_kin_au = E_min_au
while (E_kin_au <= E_max_au):
    Ekins.append(sciconv.hartree_to_ev(E_kin_au))
    E_kin_au = E_kin_au + E_step_au


#-------------------------------------------------------------------------
# constants / prefactors
prefac_res1 = VEr_au * rdg_au
prefac_indir1 = -1j * np.pi * VEr_au**2 * cdg_au_V
prefac_dir1 = 1j * cdg_au_V


#-------------------------------------------------------------------------
while ((t_au <= TX_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the first pulse \n')
    print 'during the first pulse'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')
    t_s = sciconv.atu_to_second(t_au)
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)

        sum_square = 0
        for nmu in range (0,n_fin_max+1):
        
            J = 0
            for nlambda in range (0,n_res_max+1):
                E_fin_au = E_fin_au_1 + E_mus[nmu]
                Er_au = Er_a_au
                E_lambda = E_lambdas[nlambda]
                W_au = W_lambda[nlambda]
                if (integ_outer == "quadrature"):
                    I1 = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), t_au)
                    res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), t_au)
    
                    dir_J1 = prefac_dir1 * I1[0] * gs_fin[0][nmu]
                    res_J1 = (prefac_res1 * res_I[0]
                              * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                    indir_J1 = (prefac_indir1 * res_I[0]
                                * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
    
                elif (integ_outer == "romberg"):
                    I1 = ci.complex_romberg(fun_t_dir_1, (-TX_au/2), t_au)
                    res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), t_au)
                
                    dir_J1 = prefac_dir1 * I1 * gs_fin[0][nmu]
                    res_J1 = (prefac_res1 * res_I
                              * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                    indir_J1 = (prefac_indir1 * res_I
                                * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
    
                J = (J
                     + dir_J1
                     + res_J1
                     + indir_J1
                     )
    
            square = np.absolute(J)**2
            sum_square = sum_square + square
        squares = np.append(squares, sum_square)

        string = in_out.prep_output(sum_square, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au
    
    
    in_out.doout_1f(pure_out, outlines)
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print Ekins[max_pos[i]], squares[max_pos[i]]
            outfile.write(str(Ekins[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')
    

    t_au = t_au + timestep_au




#-------------------------------------------------------------------------
while (t_au >= TX_au/2 and (t_au <= (delta_t_au - a)) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('between the pulses \n')
    print 'between the pulses'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)

        sum_square = 0
        for nmu in range (0,n_fin_max+1):

            J = 0
            for nlambda in range (0,n_res_max+1):
                E_fin_au = E_fin_au_1 + E_mus[nmu]
                Er_au = Er_a_au
                E_lambda = E_lambdas[nlambda]
                W_au = W_lambda[nlambda]
      
                if (integ_outer == "quadrature"):
                    I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                    res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)
    
                    dir_J1 = prefac_dir1 * I1[0] * gs_fin[0][nmu]
                    res_J1 = (prefac_res1 * res_I[0]
                              * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                    indir_J1 = (prefac_indir1 * res_I[0]
                                * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
                
                elif (integ_outer == "romberg"):
                    I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                    res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)

                    dir_J1 = prefac_dir1 * I1 * gs_fin[0][nmu]
                    res_J1 = (prefac_res1 * res_I
                              * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                    indir_J1 = (prefac_indir1 * res_I
                                * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
    
                J = (J
                     + dir_J1
                     + res_J1
                     + indir_J1
                     )
    
            square = np.absolute(J)**2
            sum_square = sum_square + square
        squares = np.append(squares, sum_square)

        string = in_out.prep_output(sum_square, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au

    
    
    in_out.doout_1f(pure_out,outlines)
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print Ekins[max_pos[i]], squares[max_pos[i]]
            outfile.write(str(Ekins[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')

    t_au = t_au + timestep_au







outfile.close
pure_out.close
movie_out.close
