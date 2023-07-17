#!/usr/bin/python

##########################################################################
#                                    ELDEST                              #
#        Investigating Electronic Decay Processes with Streaking         #
##########################################################################
# Purpose:                                                               #
#          - A program to simulate the time-resolved RICD spectroscopy   #
#            including quantum nuclear dynamics.                         #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer November 2020                               #
# extended by: Alexander Riegel July 2023                                #
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
print(infile)

#-------------------------------------------------------------------------
# open outputfile
outfile = open("eldest.out", mode='w')
pure_out = open('full.dat', mode='w')
movie_out = open('movie.dat', mode='w')
#popfile = open("pop.dat", mode='w')

outfile.write("The results were obtained with nuclear_dyn.py \n")
#-------------------------------------------------------------------------
# set some defaults
Xshape = 'convoluted'

#-------------------------------------------------------------------------
# read inputfile
# (see next section for explanations of most symbols)
# ( * X_sinsq, X_gauss are simply Booleans, created by in_out from X_shape)
# ( * phi is the phase for the IR pulse potential cosine-oscillation, a remnant from PRA 2020)
# ( * integ, integ_outer are integration schemes: [analytic,] quadrature, romberg)
# (currently NOT in use: cdg_au, tau_a_s, tau_b_s interact_eV, Lshape, shift_step_s, phi, grad_delta, R_eq_AA, gs_const, res_const)
# ( * Er_b_eV and E_fin_eV_2 will be converted to au, but these will not be used afterwards)
# ( * tau_s_2 will be converted to au at this to Gamma, but this will not be used afterwards)
# ( * omega_eV will be converted to au, from which TL and A0L are calculated, but other than being used for needless printing and for check_input, they will not be used afterwards)
# ( * n_L and I_L only lead to related qnts like TL, E0L and A0L, for which above holds)
# ( * FWHM_L will be converted to au and this printed, but not be used afterwards)
# ( * fin_d will be used to bind fin_const for Morse final potential, but both will not be used afterwards)

# (q is explicit input, not calced as q = rdg / (cdg pi VEr) = sqrt(2 tau / pi) rdg / cdg )

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
Er_a_au        = sciconv.ev_to_hartree(Er_a_eV)     # resonance E for RICD + AI
#Er_b_au        = sciconv.ev_to_hartree(Er_b_eV)     # resonance E for ICD
Er_au          = Er_a_au        # ? One could delete Er_a_au altogether
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)    # (same as for Er)
E_fin_au_1     = sciconv.ev_to_hartree(E_fin_eV)    # final E for sRICD

tau_au_1       = sciconv.second_to_atu(tau_s)       # lifetime for sRICD res. st.
tau_au         = tau_au_1                           # (same as for Er)
Gamma_au       = 1. / tau_au
Gamma_eV       = sciconv.hartree_to_ev(Gamma_au)
outfile.write('Gamma_eV = ' + str(Gamma_eV) + '\n')

# second final state
#E_fin_au_2       = sciconv.ev_to_hartree(E_fin_eV_2)
#tau_au_2         = sciconv.second_to_atu(tau_s_2)
#Gamma_au_2       = 1. / tau_au_2

# laser parameters
Omega_au      = sciconv.ev_to_hartree(Omega_eV)
if (X_sinsq):
    TX_au     = n_X * 2 * np.pi / Omega_au
elif(X_gauss):
    sigma     = np.pi * n_X / (Omega_au * np.sqrt(np.log(2)))
    FWHM      = 2 * np.sqrt( 2 * np.log(2)) * sigma
    TX_au     = 5 * sigma
    print('sigma = ', sciconv.atu_to_second(sigma))
    print('FWHM = ', sciconv.atu_to_second(FWHM))
    outfile.write('sigma = ' + str(sciconv.atu_to_second(sigma)) + '\n')
    outfile.write('FWHM = ' + str(sciconv.atu_to_second(FWHM)) + '\n')
print('end of the first pulse = ', sciconv.atu_to_second(TX_au/2))
outfile.write('end of the first pulse = ' + str(sciconv.atu_to_second(TX_au/2)) + '\n')
I_X_au        = sciconv.Wcm2_to_aiu(I_X)
print('I_X = ', I_X)
print('I_X_au = ', I_X_au)
E0X           = np.sqrt(I_X_au)
A0X           = E0X / Omega_au
print('A0X = ', A0X)

omega_au      = sciconv.ev_to_hartree(omega_eV)
#FWHM_L_au     = sciconv.second_to_atu(FWHM_L)
#sigma_L_au    = FWHM_L_au / np.sqrt(8 * np.log(2))      # assume Gaussian envelope for second pulse
#a             = 5./2 * sigma_L_au       # half duration of IR pulse (delta_t - a, delta_t + a); in PRA 2020: small-delta t
#print("FWHM_L = ", FWHM_L)
#print("sigma_L = ", sciconv.atu_to_second(sigma_L_au))
TL_au         = n_L * 2 * np.pi / omega_au
#print('start of IR pulse = ', delta_t_s - sciconv.atu_to_second(TL_au/2))
#print('end of IR pulse = ', delta_t_s + sciconv.atu_to_second(TL_au/2))
I_L_au        = sciconv.Wcm2_to_aiu(I_L)
#print('I_L = ', I_L)
#print('I_L_au = ', I_L_au)
E0L           = np.sqrt(I_L_au)
##print('E0L = ', E0L)
A0L           = E0L / omega_au
#print('A0L = ', A0L)
delta_t_au    = sciconv.second_to_atu(delta_t_s)        # t diff between the maxima of the two pulses

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
E_step_au = sciconv.ev_to_hartree(E_step_eV)

E_min_au = sciconv.ev_to_hartree(E_min_eV)
E_max_au = sciconv.ev_to_hartree(E_max_eV)

VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))
print('VEr_au = ', VEr_au)

#VEr_au_1      = VEr_au      # (same as for Er)

#test q=1
cdg_au_V = rdg_au / ( q * np.pi * VEr_au)

#-------------------------------------------------------------------------
# Potential details
# vibrational energies of Morse potentials
print()
print('-----------------------------------------------------------------')
outfile.write('\n' + '-----------------------------------------------------------------' + '\n')
red_mass = wf.red_mass_au(mass1,mass2)
print("red_mass = ", red_mass)

#ground state
print()
print("Ground state")
print('-----------------------------------------------------------------')
print("Energies of vibrational states of the ground state")
outfile.write('\n' + '-----------------------------------------------------------------' + '\n')
outfile.write("Energies of vibrational states of the ground state" + '\n')
lambda_param_gs = np.sqrt(2*red_mass*gs_de) / gs_a
n_gs_max = int(lambda_param_gs - 0.5)
print("n_gs_max = ", n_gs_max)
E_kappas = []   # collects vibr energies of GS
print('n_gs  ' + 'E [au]            ' + 'E [eV]')
outfile.write('n_gs  ' + 'E [au]            ' + 'E [eV]' + '\n')
for n in range (0,n_gs_max+1):
    ev = wf.eigenvalue(n,gs_de,gs_a,red_mass)   # ev stands for eigenvalue, not for electronvolt (it is, in fact, in au!)
    E_kappas.append(ev)
    outfile.write('{:4d}  {:14.10E}  {:14.10E}\n'.format(n,ev,sciconv.hartree_to_ev(ev)))
    print('{:4d}  {:14.10E}  {:14.10E}'.format(n,ev,sciconv.hartree_to_ev(ev)))

#resonant state
print()
print("Resonant state")
print('-----------------------------------------------------------------')
print("Energies of vibrational states of the resonant state")
outfile.write('\n' + '-----------------------------------------------------------------' + '\n')
outfile.write("Energies of vibrational states of the resonant state" + '\n')
lambda_param_res = np.sqrt(2*red_mass*res_de) / res_a
n_res_max = int(lambda_param_res - 0.5)
print("n_res_max = ", n_res_max)
E_lambdas = []
outfile.write('n_res  ' + 'E [au]            ' + 'E [eV]' + '\n')
print('n_res  ' + 'E [au]            ' + 'E [eV]')
for n in range (0,n_res_max+1):
    ev = wf.eigenvalue(n,res_de,res_a,red_mass)
    E_lambdas.append(ev)
    outfile.write('{:5d}  {:14.10E}  {:14.10E}\n'.format(n,ev,sciconv.hartree_to_ev(ev)))
    print('{:5d}  {:14.10E}  {:14.10E}'.format(n,ev,sciconv.hartree_to_ev(ev)))

#final state
print()
print("Final state")
print('-----------------------------------------------------------------')
print("Energies of vibrational states of the final state")
outfile.write('\n' + '-----------------------------------------------------------------' + '\n')
outfile.write("Energies of vibrational states of the final state" + '\n')
if (fin_pot_type == 'morse'):
    fin_de    = fin_a
    fin_a     = fin_b
    fin_Req   = fin_c
    fin_const = fin_d
    lambda_param_fin = np.sqrt(2*red_mass*fin_de) / fin_a
    n_fin_max = int(lambda_param_fin - 0.5)     # Maximum quantum number = n_fin_max -> number of states = n_fin_max + 1
    print("n_fin_max = ", n_fin_max)
    E_mus = []
    print('n_fin  ' + 'E [au]            ' + 'E [eV]')
    outfile.write('n_fin  ' + 'E [au]            ' + 'E [eV]' + '\n')
    for n in range (0,n_fin_max+1):
        ev = wf.eigenvalue(n,fin_de,fin_a,red_mass)
        E_mus.append(ev)
        outfile.write('{:5d}  {:14.10E}  {:14.10E}\n'.format(n,ev,sciconv.hartree_to_ev(ev)))
        print('{:5d}  {:14.10E}  {:14.10E}'.format(n,ev,sciconv.hartree_to_ev(ev)))
if (fin_pot_type == 'hyperbel'):
    print('Final state is repulsive')
    outfile.write('Final state is repulsive' + '\n')
    fin_hyp_a = fin_a
    fin_hyp_b = fin_b
    E_fin_au = fin_hyp_b        # Since for an all-repulsive state there is no minimum (E_fin), E_fin is set to the final potential at infinite distance, i.e. fin_hyp_b
    E_fin_au_1 = fin_hyp_b
    E_hyp_step = fin_c
    n_fin_max_list = []              # Max quantum number considered for each lambda (all vibr fin states above the resp res state are discarded)
    for E_l in E_lambdas:
        n_fin_max_list.append(int((Er_a_au + E_l - E_fin_au) / E_hyp_step))
    outfile.write('Continuous vibrational states of the final state are discretized,\nstep width: {:5E} au   = {:5E} eV\n'.format(
        E_hyp_step, sciconv.hartree_to_ev(E_hyp_step)))
    outfile.write('Thus, up to {} final vibrational states will be considered\n'.format(n_fin_max_list[-1]))
    print('Continuous vibrational states of the final state are discretized,\nstep width: {:5E} au   = {:5E} eV'.format(
        E_hyp_step, sciconv.hartree_to_ev(E_hyp_step)))
    print('Thus, up to {} final vibrational states will be considered'.format(n_fin_max_list[-1]))
    E_mus = []
    E_mu = fin_hyp_b
    for n in range (0,n_fin_max_list[-1]+1):
        E_mus.append(E_mu)
        E_mu = E_mu + E_hyp_step

#-------------------------------------------------------------------------
# Franck-Condon factors
#-------------------------------------------------------------------------
gs_res =  []    # collects lists of FC: [<l1|k1>, <l2|k1>, ...], [<l1|k2, <l2|k2>, ...], ...
gs_fin =  []
res_fin = []
R_min = sciconv.angstrom_to_bohr(1.5)
R_max = sciconv.angstrom_to_bohr(30.0)

# ground state - resonant state <lambda|kappa>
print()
print('-----------------------------------------------------------------')
print("Franck-Condon overlaps between ground and resonant state")
print('n_gs  ' + 'n_res  ' + '<res|gs>')
outfile.write('\n' + '-----------------------------------------------------------------' + '\n')
outfile.write("Franck-Condon overlaps between ground and resonant state" + '\n')
outfile.write('n_gs  ' + 'n_res  ' + '<res|gs>' + '\n')
for i in range (0,n_gs_max+1):
    tmp = []
    for j in range (0,n_res_max+1):
        FC = wf.FC(j,res_a,res_Req,res_de,red_mass,
                   i,gs_a,gs_Req,gs_de,R_min,R_max)
        tmp.append(FC)
        outfile.write('{:4d}  {:5d}  {:14.10E}\n'.format(i,j,FC))
        print(('{:4d}  {:5d}  {:14.10E}'.format(i,j,FC)))
    gs_res.append(tmp)
    
# ground state - final state <mu|kappa>
print()
print('-----------------------------------------------------------------')
print("Franck-Condon overlaps between ground and final state")
outfile.write('\n' + '-----------------------------------------------------------------' + '\n')
outfile.write("Franck-Condon overlaps between ground and final state" + '\n')
if (fin_pot_type == 'hyperbel'):
    n_fin_max = n_fin_max_list[-1]      # Calculate FC for all final states ever to be considered, so get highest n_fin_max
print('n_gs  ' +'n_fin  ' + '<fin|gs>')
outfile.write('n_gs  ' +'n_fin  ' + '<fin|gs>' + '\n')
for i in range (0,n_gs_max+1):
    tmp = []
    for j in range (0,n_fin_max+1):
        if (fin_pot_type == 'morse'):
            FC = wf.FC(j,fin_a,fin_Req,fin_de,red_mass,
                i,gs_a,gs_Req,gs_de,R_min,R_max)
            tmp.append(FC)
            outfile.write('{:4d}  {:5d}  {:14.10E}\n'.format(i,j,FC))
            print(('{:4d}  {:5d}  {:14.10E}'.format(i,j,FC)))
        elif (fin_pot_type == 'hyperbel'):
            if (j == 0):
                R_start = np.inf        # Technically that is the vibr ground state for the elec final state
            else:
                R_start = fin_hyp_a / (E_mus[j] - fin_hyp_b)
            FC = wf.FCmor_freehyp(i,gs_a,gs_Req,gs_de,red_mass,
                fin_hyp_a,fin_hyp_b,R_start,R_min,R_max,limit=500)
            tmp.append(FC)
            if (j == 0 or j == n_fin_max-1 or j == n_fin_max):      # Don't print all the FC, just the first two and last two (per GS vibr state)
                outfile.write('{:4d}  {:5d}  {: 14.10E}\n'.format(i,j,FC))
                print(('{:4d}  {:5d}  {: 14.10E}'.format(i,j,FC)))
            elif (j == 1):
                outfile.write('{:4d}  {:5d}  {: 14.10E}\n'.format(i,j,FC))
                outfile.write('   ...\n')
                print(('{:4d}  {:5d}  {: 14.10E}'.format(i,j,FC)))
                print('   ...')
    gs_fin.append(tmp)

# resonant state - final state <mu|lambda>
print()
print('-----------------------------------------------------------------')
print("Franck-Condon overlaps between final and resonant state")
outfile.write('\n' + '-----------------------------------------------------------------' + '\n')
outfile.write("Franck-Condon overlaps between final and resonant state" + '\n')
print('n_res  ' +'n_fin  ' + '<fin|res>')
outfile.write('n_res  ' +'n_fin  ' + '<fin|res>' + '\n')
for i in range (0,n_res_max+1):
    tmp = []
    if (fin_pot_type == 'hyperbel'):
        n_fin_max = n_fin_max_list[i]
    for j in range (0,n_fin_max+1):
        if (fin_pot_type == 'morse'):
            FC = wf.FC(j,fin_a,fin_Req,fin_de,red_mass,
                       i,res_a,res_Req,res_de,R_min,R_max)
            tmp.append(FC)
            outfile.write('{:5d}  {:5d}  {:14.10E}\n'.format(i,j,FC))
            print(('{:5d}  {:5d}  {:14.10E}'.format(i,j,FC)))
        elif (fin_pot_type == 'hyperbel'):
            if (j == 0):
                R_start = np.inf
            else:
                R_start = fin_hyp_a / (E_mus[j] - fin_hyp_b)
            FC = wf.FCmor_freehyp(i,res_a,res_Req,res_de,red_mass,
                fin_hyp_a,fin_hyp_b,R_start,R_min,R_max,limit=500)
            tmp.append(FC)
            if (j == 0 or j == n_fin_max-1 or j == n_fin_max):
                outfile.write('{:5d}  {:5d}  {: 14.10E}\n'.format(i,j,FC))
                print(('{:5d}  {:5d}  {: 14.10E}'.format(i,j,FC)))
            elif (j == 1):
                outfile.write('{:5d}  {:5d}  {: 14.10E}\n'.format(i,j,FC))
                outfile.write('    ...\n')
                print(('{:5d}  {:5d}  {: 14.10E}'.format(i,j,FC)))
                print('    ...')
    res_fin.append(tmp)

# sum over mup of product <lambda|mup><mup|kappa>       where mup means mu prime
indir_FCsums = []
for i in range (0,n_res_max+1):
    indir_FCsum = 0
    if (fin_pot_type == 'hyperbel'):
        n_fin_max = n_fin_max_list[i]
    for j in range (0,n_fin_max+1):
        tmp = np.conj(res_fin[i][j]) * gs_fin[0][j]     # = <mu_j|lambda_i>* <mu_j|kappa_0> = <lambda_i|mu_j><mu_j|kappa_0> = <li|mj><mj|k0>
        indir_FCsum = indir_FCsum + tmp                 # = sum_j <li|mj><mj|k0>
    indir_FCsums.append(indir_FCsum)                    # = [sum_j <l1|mj><mj|k0>, sum_j <l2|mj><mj|k0>, ...]
print()
print('-----------------------------------------------------------------')
outfile.write('\n' + '-----------------------------------------------------------------' + '\n')

#-------------------------------------------------------------------------
# determine total decay width matrix element
print('Effective decay widths in eV and lifetimes in s:')
outfile.write('Effective decay widths in eV and lifetimes in s:' + '\n')
W_lambda = []   # [W_(l=0), W_(l=1), ...]
for i in range (0,n_res_max+1):
    tmp = 0
    if (fin_pot_type == 'hyperbel'):
        n_fin_max = n_fin_max_list[i]       # To each lambda their own n_fin_max (v.s.)
    for j in range (0,n_fin_max+1):
        tmp = tmp + VEr_au**2 * np.abs(res_fin[i][j])**2      # W_li = sum_j ( VEr**2 |<mj|li>|**2 )
    W_lambda.append(tmp)
    ttmp = 1./ (2 * np.pi * tmp)        # lifetime tau_l = 1 / (2 pi W_l)
    print(sciconv.hartree_to_ev(tmp), sciconv.atu_to_second(ttmp))
    outfile.write( str(sciconv.hartree_to_ev(tmp)) + ' '
                 + str(sciconv.atu_to_second(ttmp)) + '\n')
print()


#-------------------------------------------------------------------------
in_out.check_input(Er_au, E_fin_au, Gamma_au,
                   Omega_au, TX_au, n_X, A0X,
                   omega_au, TL_au, A0L, delta_t_au,
                   tmax_au, timestep_au, E_step_au)
#-------------------------------------------------------------------------
# physical definitions of functions
# functions for the shape of the XUV pulse
if (X_sinsq):
    print('use sinsq function')
    f_t1  = lambda t1: 1./4 * ( np.exp(2j * np.pi * (t1 + TX_au/2) / TX_au) # There should be a minus sign in front [from (1/2i)**2]
                          + 2                                               # & this 2 be negative [sin**2 = (exp - exp*)**2 = exp**2 - 2 exp exp* + (exp*)**2]
                          + np.exp(-2j * np.pi * (t1 + TX_au/2) /TX_au) )
    # fp_t1 = f'(t1)
    fp_t1 = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (t1 + TX_au/2) / TX_au)  # Accordingly, these signs must be flipped
                                         + np.exp(-2j*np.pi* (t1 + TX_au/2) / TX_au) )
elif (X_gauss):
    print('use gauss function')
    f_t1  = lambda t1: ( 1./ np.sqrt(2*np.pi * sigma**2)
                       * np.exp(-t1**2 / (2*sigma**2)))
    # fp_t1 = f'(t1)
    fp_t1 = lambda t1: ( -t1 / np.sqrt(2*np.pi) / sigma**3
                       * np.exp(-t1**2 / (2*sigma**2)))
else:
    print('no pulse shape selected')

print()

if (Xshape == 'convoluted'):    # Calculate field strength EX = -(AX fX)'
    FX_t1 = lambda t1: (
                        0
                        - (A0X
                           * np.cos(Omega_au * t1)
                           * fp_t1(t1)
                          )
                        + (A0X
                           * Omega_au
                           * np.sin(Omega_au * t1)
                           * f_t1(t1)
                          )
                       )
elif (Xshape == 'infinite'):
    FX_t1 = lambda t1: + A0X * Omega_au * np.cos(Omega_au * t1)
    #FX_t1 = lambda t1: - A0X * np.sin(Omega_au * t1)
                       

#-------------------------------------------------------------------------
# technical definitions of functions (remember: FX is the field strength EX)
#direct ionization
fun_t_dir_1 = lambda t1: FX_t1(t1)   * np.exp(1j * E_fin_au * (t1-t_au)) \
                                     * np.exp(1j * E_kin_au * (t1-t_au))        # Note: before any of these fncts are called, E_fin is redefined to also include E_mu
fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
                                     * np.exp(1j * E_kin_au * (t1-t_au))        # Same as fun_t_dir_1 - why keep ?

#res_inner_fun = lambda t2: np.exp(-t2 * (np.pi * W_au + 1j*(Er_au))) \
#                           * IR_during(t2)

if (integ == 'romberg'):                                                        # numerical inner int not possible (res_inner_fun deactivated) ?
    res_inner = lambda t1: integrate.romberg(res_inner_fun, t1, t_au)
elif (integ == 'quadrature'):
    res_inner = lambda t1: integrate.quad(res_inner_fun, t1, t_au)[0]
elif (integ == 'analytic'):
    # analytic inner integral
    res_inner = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au - E_lambda)     # See the above note on E_fin also including E_mu
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
                           * np.exp(t1 * (np.pi* W_au + 1j*(Er_au + E_lambda))) \
                           * res_inner(t1)


#-------------------------------------------------------------------------
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
prefac_res1 = VEr_au * rdg_au
prefac_indir1 = -1j * np.pi * VEr_au**2 * cdg_au_V
prefac_dir1 = 1j * cdg_au_V

if (fin_pot_type == 'hyperbel'):
    threshold = fin_d

########################################
# now follow the integrals themselves, for the temporal phases:
# 'during the first pulse' (-TX/2, TX/2)
# 'between the pulses' (TX/2, tmax)

#-------------------------------------------------------------------------
while ((t_au <= TX_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the first pulse \n')
    print('during the first pulse')

    outlines = []       # will contain lines containing triples of E_kin, time and signal intensity
    squares = np.array([])  # signal intensity ( = |amplitude|**2 = |J|**2 )
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    print('t_s = ', t_s)
    outfile.write('t_s = ' + str(t_s) + '\n')
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)
        print(sciconv.hartree_to_ev(E_kin_au)) #?

        if (fin_pot_type == 'morse'):
            sum_square = 0      # Total spectrum |J @ E_kin|**2 = sum_mu |J_mu @ E_kin|**2      (sum of contributions of all final states with E_kin)
            for nmu in range (0,n_fin_max+1):           # loop over all mu, calculate J_mu = J_dir,mu + J_nondir,mu
                E_fin_au = E_fin_au_1 + E_mus[nmu]      # E_fin_au_1: inputted electronic E_fin_au, E_mus: Morse vibrational eigenvalues of fin state
    #            Er_au = Er_a_au
                
                # Direct term
                if (integ_outer == "quadrature"):
                    I1 = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), t_au)
                    dir_J1 = prefac_dir1 * I1[0] * gs_fin[0][nmu]        # [0] of quad integ result = integral (rest is est error & info); FC = <mu_n|kappa_0>
        
                elif (integ_outer == "romberg"):
                    I1 = ci.complex_romberg(fun_t_dir_1, (-TX_au/2), t_au)
                    dir_J1 = prefac_dir1 * I1 * gs_fin[0][nmu]           # romberg returns only the integral, so no [0] necessary

                # J_nondir,mu = sum_lambda J_nondir,mu,lambda = sum_lambda (J_res,mu,lambda + J_indir,mu,lambda)
                J = 0
                for nlambda in range (0,n_res_max+1):
                    E_lambda = E_lambdas[nlambda]
                    W_au = W_lambda[nlambda]
                    if (integ_outer == "quadrature"):
                        res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), t_au)
        
                        res_J1 = (prefac_res1 * res_I[0]
                                  * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                        indir_J1 = (prefac_indir1 * res_I[0]
                                    * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
        
                    elif (integ_outer == "romberg"):
                        res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), t_au)
                    
                        res_J1 = (prefac_res1 * res_I
                                  * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                        indir_J1 = (prefac_indir1 * res_I
                                    * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
        
                    J = (J
                         + res_J1
                         + indir_J1
                         )
        
                square = np.absolute(J + dir_J1)**2     # |J_mu|**2
                sum_square = sum_square + square        # |J|**2 = sum_mu |J_mu|**2
        

        if (fin_pot_type == 'hyperbel'):
            sum_square = 0      # replace integral dE_mu |J_(E_mu) @ E_kin|**2  by  E_hyp_step * sum_(n_fin) |J_dir,mu + J_nondir,mu|**2
            dir_J1 = threshold + 5  # Initialize J_dir arbitrarily above threshold (to pass initial subsequent "while" check)
            nmu = 0
            while (nmu <= n_fin_max_list[-1]) or (np.abs(dir_J1) >= threshold): # Cap sum if n_fin > max n_fin_max & J_dir fell sub-threshold
                if (nmu > n_fin_max_list[-1]):           # If precalced list of gs-fin-FCF is exhausted but still J_dir >= threshold:
                    E_mus.append(E_mus[-1] + E_hyp_step) # a) calc new E_mu (above even highest |res>|lambda>), 1 step above prev max E_mu
                    for i in range (0,n_gs_max+1):       # b) calc new gs-fin-FCF (one for each gs vibr state kappa)
                        R_start = fin_hyp_a / (E_mus[nmu] - fin_hyp_b)
                        FC = wf.FCmor_freehyp(i,gs_a,gs_Req,gs_de,red_mass,
                                               fin_hyp_a,fin_hyp_b,R_start,R_min,R_max,limit=500)
                        gs_fin[i].append(FC)             # Append each k-list of <k|m> by one value (with the current mu)
                        #print(nmu, E_mus[-1], gs_fin[0][-1])   #?
                E_fin_au = E_fin_au_1 + E_mus[nmu]       # Again: E_fin_au_1 is the potential at infinity
    #            Er_au = Er_a_au
                
                # Direct term
                if (integ_outer == "quadrature"):
                    I1 = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), t_au)
                    dir_J1 = prefac_dir1 * I1[0] * gs_fin[0][nmu]
                    #print(nmu, np.abs(dir_J1)) #?
        
                elif (integ_outer == "romberg"):
                    I1 = ci.complex_romberg(fun_t_dir_1, (-TX_au/2), t_au)
                    dir_J1 = prefac_dir1 * I1 * gs_fin[0][nmu]
                 
                # J_nondir,mu = sum_lambda J_nondir,mu,lambda = sum_lambda (J_res,mu,lambda + J_indir,mu,lambda)
                J = 0
                for nlambda in range (0,n_res_max+1):
                    if nmu > n_fin_max_list[nlambda]:   # J_nondir,mu,lambda = 0 if |fin>|mu> lies higher than |res>|lambda>
                        continue                        # -> Also no new res-FCF must be calced, since if list is exhausted, all J_nondir = 0
                    E_lambda = E_lambdas[nlambda]
                    W_au = W_lambda[nlambda]
                    if (integ_outer == "quadrature"):
                        res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), t_au)
        
                        res_J1 = (prefac_res1 * res_I[0]
                                  * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                        indir_J1 = (prefac_indir1 * res_I[0]
                                    * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
        
                    elif (integ_outer == "romberg"):
                        res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), t_au)
                    
                        res_J1 = (prefac_res1 * res_I
                                  * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                        indir_J1 = (prefac_indir1 * res_I
                                    * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
        
                    J = (J
                         + res_J1
                         + indir_J1
                         )
        
                square = np.absolute(J + dir_J1)**2     # |J_mu|**2
                sum_square = sum_square + square        # |J|**2 = sum_mu |J_mu|**2
                nmu = nmu + 1

        squares = np.append(squares, sum_square)

        string = in_out.prep_output(sum_square, E_kin_au, t_au)     # returns str: E_kin_eV, t_s, sum_square = intensity
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au     # @ t = const.
    
    
    in_out.doout_1f(pure_out, outlines)     # writes each (E_kin, t = const, |J|**2) triple in a sep line into output file
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]      # finds position of relative (i. e. local) maxima of |J|**2 in an array
    if (len(max_pos > 0)):                               # if there are such:
        for i in range (0, len(max_pos)):
            print(Ekins[max_pos[i]], squares[max_pos[i]])      # print all loc max & resp E_kin
            outfile.write(str(Ekins[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')
    

    t_au = t_au + timestep_au



#-------------------------------------------------------------------------
while (t_au >= TX_au/2\
#        and (t_au <= (delta_t_au - a))\
        and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('between the pulses \n')
    print('between the pulses')

    # all equal to during-1st-pulse section, except for integrating over entire XUV pulse now

    outlines = []       # will contain lines containing triples of E_kin, time and signal intensity
    squares = np.array([])  # signal intensity ( = |amplitude|**2 = |J|**2 )
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    print('t_s = ', t_s)
    outfile.write('t_s = ' + str(t_s) + '\n')
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)
        print(sciconv.hartree_to_ev(E_kin_au)) #?
        

        if (fin_pot_type == 'morse'):
            sum_square = 0      # Total spectrum |J @ E_kin|**2 = sum_mu |J_mu @ E_kin|**2      (sum of contributions of all final states with E_kin)
            for nmu in range (0,n_fin_max+1):           # loop over all mu, calculate J_mu = J_dir,mu + J_nondir,mu
                E_fin_au = E_fin_au_1 + E_mus[nmu]      # E_fin_au_1: inputted electronic E_fin_au, E_mus: Morse vibrational eigenvalues of fin state
    #            Er_au = Er_a_au
                
                # Direct term
                if (integ_outer == "quadrature"):
                    I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)      # same function as fun_t_dir_1 before,
                    dir_J1 = prefac_dir1 * I1[0] * gs_fin[0][nmu]
    
                elif (integ_outer == "romberg"):
                    I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                    dir_J1 = prefac_dir1 * I1 * gs_fin[0][nmu]

                # J_nondir,mu = sum_lambda J_nondir,mu,lambda = sum_lambda (J_res,mu,lambda + J_indir,mu,lambda)
                J = 0
                for nlambda in range (0,n_res_max+1):
                    E_lambda = E_lambdas[nlambda]
                    W_au = W_lambda[nlambda]
                    if (integ_outer == "quadrature"):
                        res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)
        
                        res_J1 = (prefac_res1 * res_I[0]
                                  * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                        indir_J1 = (prefac_indir1 * res_I[0]
                                    * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
        
                    elif (integ_outer == "romberg"):
                        res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
                    
                        res_J1 = (prefac_res1 * res_I
                                  * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                        indir_J1 = (prefac_indir1 * res_I
                                    * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
        
                    J = (J
                         + res_J1
                         + indir_J1
                         )
        
                square = np.absolute(J + dir_J1)**2     # |J_mu|**2
                sum_square = sum_square + square        # |J|**2 = sum_mu |J_mu|**2
        

        if (fin_pot_type == 'hyperbel'):
            sum_square = 0      # replace integral dE_mu |J_(E_mu) @ E_kin|**2  by  E_hyp_step * sum_(n_fin) |J_dir,mu + J_nondir,mu|**2
            dir_J1 = threshold + 5  # Initialize J_dir arbitrarily above threshold (to pass initial subsequent "while" check)
            nmu = 0
            while (nmu <= n_fin_max_list[-1]) or (np.abs(dir_J1) >= threshold): # Cap sum if n_fin > max n_fin_max & J_dir fell sub-threshold
                if (nmu > n_fin_max_list[-1]):           # If precalced list of gs-fin-FCF is exhausted but still J_dir >= threshold:
                    E_mus.append(E_mus[-1] + E_hyp_step) # a) calc new E_mu (above even highest |res>|lambda>), 1 step above prev max E_mu
                    for i in range (0,n_gs_max+1):       # b) calc new gs-fin-FCF (one for each gs vibr state kappa)
                        R_start = fin_hyp_a / (E_mus[nmu] - fin_hyp_b)
                        FC = wf.FCmor_freehyp(i,gs_a,gs_Req,gs_de,red_mass,
                                               fin_hyp_a,fin_hyp_b,R_start,R_min,R_max,limit=500)
                        gs_fin[i].append(FC)             # Append each k-list of <k|m> by one value (with the current mu)
                        #print(nmu, E_mus[-1], gs_fin[0][-1])   #?
                E_fin_au = E_fin_au_1 + E_mus[nmu]       # Again: E_fin_au_1 is the potential at infinity
    #            Er_au = Er_a_au
                
                # Direct term
                if (integ_outer == "quadrature"):
                    I1 = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), TX_au/2)
                    dir_J1 = prefac_dir1 * I1[0] * gs_fin[0][nmu]
                    #if (nmu > n_fin_max_list[-1]): print(nmu, dir_J1)  #?
        
                elif (integ_outer == "romberg"):
                    I1 = ci.complex_romberg(fun_t_dir_1, (-TX_au/2), TX_au/2)
                    dir_J1 = prefac_dir1 * I1 * gs_fin[0][nmu]
                 
                # J_nondir,mu = sum_lambda J_nondir,mu,lambda = sum_lambda (J_res,mu,lambda + J_indir,mu,lambda)
                J = 0
                for nlambda in range (0,n_res_max+1):
                    if nmu > n_fin_max_list[nlambda]:   # J_nondir,mu,lambda = 0 if |fin>|mu> lies higher than |res>|lambda>
                        continue                        # -> Also no new res-FCF must be calced, since if list is exhausted, all J_nondir = 0
                    E_lambda = E_lambdas[nlambda]
                    W_au = W_lambda[nlambda]
                    if (integ_outer == "quadrature"):
                        res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)
        
                        res_J1 = (prefac_res1 * res_I[0]
                                  * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                        indir_J1 = (prefac_indir1 * res_I[0]
                                    * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
        
                    elif (integ_outer == "romberg"):
                        res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
                    
                        res_J1 = (prefac_res1 * res_I
                                  * gs_res[0][nlambda] * res_fin[nlambda][nmu])
                        indir_J1 = (prefac_indir1 * res_I
                                    * indir_FCsums[nlambda] * res_fin[nlambda][nmu])
        
                    J = (J
                         + res_J1
                         + indir_J1
                         )
        
                square = np.absolute(J + dir_J1)**2     # |J_mu|**2
                sum_square = sum_square + square        # |J|**2 = sum_mu |J_mu|**2
                nmu = nmu + 1

        squares = np.append(squares, sum_square)

        string = in_out.prep_output(sum_square, E_kin_au, t_au)     # returns str: E_kin_eV, t_s, sum_square = intensity
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au     # @ t = const.
    
    
    in_out.doout_1f(pure_out, outlines)     # writes each (E_kin, t = const, |J|**2) triple in a sep line into output file
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]      # finds position of relative (i. e. local) maxima of |J|**2 in an array
    if (len(max_pos > 0)):                               # if there are such:
        for i in range (0, len(max_pos)):
            print(Ekins[max_pos[i]], squares[max_pos[i]])      # print all loc max & resp E_kin
            outfile.write(str(Ekins[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')
    

    t_au = t_au + timestep_au








outfile.close
pure_out.close
movie_out.close
