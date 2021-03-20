#!/usr/bin/python

import numpy as np
import sciconv as sc
import wellenfkt as wf
import sys
import warnings
import math

# don't print warnings unless python -W ... is used
if not sys.warnoptions:
    warnings.simplefilter("ignore")

#--------------------------------------------------------------------------
# input parameters
mu = wf.red_mass_au(20.1797,20.1797)

gs_de     = 0.0001102
gs_a      = 1.5
#gs_a      = 3.8
gs_Req    = 6.0
gs_const  = 0.0

nmax_res = 1
nmax_fin = 1
n_vis_resstates = 2
n_vis_finstates = 4

De_min_eV = 0.01
De_max_eV = 6.0
#De_max_eV = 0.2
alpha_max = 30.0

De_step_eV = 0.1
alpha_step = 0.5
#alpha_in_step = 0.05

FCmin = 0.4
FCmin_fin = 0.2
FCmin_res = 0.2
FCmax = 1.0E-3

minEdiff_eV = 0.19

res_Req = 6.0
fin_Req = 6.0
R_min = sc.angstrom_to_bohr(1.5)
R_max = sc.angstrom_to_bohr(30.0)
#--------------------------------------------------------------------------
#convert units
De_min = sc.ev_to_hartree(De_min_eV)
De_max = sc.ev_to_hartree(De_max_eV)
minEdiff = sc.ev_to_hartree(minEdiff_eV)
De_step = sc.ev_to_hartree(De_step_eV)
#--------------------------------------------------------------------------
#definition of functions

def alpha_min(mu,De,nmax):
    return 2* np.sqrt(2*mu*De) / (1 + 2*(nmax+1))

def lambda_param(mu,De,alpha):
    return np.sqrt(2*mu*De) / alpha

def is_correct_nmax(l_param,nmax):
    calc_n = int(l_param - 0.5)
    if (calc_n == nmax):
        return True
    else:
        return False
#--------------------------------------------------------------------------

#--------------------------------------------------------------------
# 2 lambda, viele (10) mu
#--------------------------------------------------------------------
#start with resonant state
De_res = De_min
while (De_res < De_max):

    res_alpha = alpha_min(mu,De_res,nmax_res)
    
    while (is_correct_nmax(lambda_param(mu,De_res,res_alpha+alpha_step),nmax_res)
           and (res_alpha <= alpha_max)):
        res_alpha = res_alpha + alpha_step

        # consider only cases with minimum energy gap between vibrational states
        res_evs = []
        for i in range(0,nmax_res+1):
            tmp = wf.eigenvalue(i,De_res,res_alpha,mu)
            res_evs.append(tmp)
        #if (nmax_res == 1):
        #    if ((res_evs[1]-res_evs[0]) < minEdiff):
        #        continue

            
        # consider only cases, where the overlap of both res vib states with the
        # ground state is sufficient
        FCres = []
        for i in range(0,nmax_res+1):
            FC_tmp = wf.FC(i,res_alpha,res_Req,De_res,mu,
                      0,gs_a,gs_Req,gs_de,R_min,R_max)
            FCres.append(FC_tmp)
        if any(math.isnan(x) for x in FCres):
            break
        
        n_realistic_res_states = sum(abs(x) >= FCmin_res for x in FCres)
        if (n_realistic_res_states < n_vis_resstates):
        #if not all(abs(x) >= FCmin_fin for x in FCfins):
            continue


        De_fin = De_min
        while (De_fin < De_max):

            # consider only cases with minimum energy gap between vibrational states
            fin_evs = []

            fin_alpha = alpha_min(mu,De_fin,nmax_fin)
            while (is_correct_nmax(lambda_param(mu,De_fin,fin_alpha+alpha_step),nmax_fin) and (fin_alpha <= alpha_max)):
                fin_alpha = fin_alpha + alpha_step

                for i in range(0,nmax_fin+1):
                    tmp = wf.eigenvalue(i,De_fin,fin_alpha,mu)
                    fin_evs.append(tmp)
                if (nmax_fin == 1):
                    if ((fin_evs[1]-fin_evs[0]) < minEdiff):
                        continue

                # consider only cases, where the overlap of both res vib states with the
                # ground state is sufficient
                FCfins = []
                for i in range(0,nmax_res+1):
                    for j in range(0,nmax_fin+1):
                        FC_tmp = wf.FC(i,res_alpha,res_Req,De_res,mu,
                                  j,fin_alpha,fin_Req,De_fin,R_min,R_max)
                        FCfins.append(FC_tmp)
                if any(math.isnan(x) for x in FCfins):
                    break
                #print FCfins
                n_realistic_fin_states = sum(abs(x) >= FCmin_fin for x in FCfins)
                if (n_realistic_fin_states < n_vis_finstates):
                #if not all(abs(x) >= FCmin_fin for x in FCfins):
                    continue
                print "n_realistic_res_states = ", n_realistic_res_states
                print "n_realistic_fin_states = ", n_realistic_fin_states

                print "YES!!!!!!!!!!!!!!!!!!!!!!!!"
                print 'resonant', sc.hartree_to_ev(De_res), De_res, res_alpha
                print 'final   ', sc.hartree_to_ev(De_fin), De_fin, fin_alpha
                print 'FCres', FCres
                print 'FCfins', FCfins
                if (nmax_res+1 >= 2):
                    print 'Ediff_res', sc.hartree_to_ev(res_evs[1]-res_evs[0])
                if (nmax_fin+1 >= 2):
                    print 'Ediff_fin', sc.hartree_to_ev(fin_evs[1]-fin_evs[0])
                print '---------------------------------------------------'
                
            

            De_fin = De_fin + De_step


    De_res = De_res + De_step


