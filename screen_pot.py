#!/usr/bin/python

import numpy as np
import sciconv as sc
import wellenfkt as wf
import sys
import warnings

# don't print warnings unless python -W ... is used
if not sys.warnoptions:
    warnings.simplefilter("ignore")

mu = wf.red_mass_au(20.1797,20.1797)

gs_de     = 0.0001102
gs_a      = 1.5
gs_Req    = 6.0
gs_const  = 0.0

nmax_res = 1
nmax_fin = 0

De_min_eV = 0.1
De_max_eV = 8.0

De_step_eV = 0.01
alpha_step = 0.01

FCmin = 0.1
FCmax = 1.0E-5

minEdiff_eV = 0.25

Req = 6.0
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

#start with resonant state
De_res = De_min
while (De_res < De_max):

    res_alpha = alpha_min(mu,De_res,nmax_res)# + alpha_step
    
    while (is_correct_nmax(lambda_param(mu,De_res,res_alpha+alpha_step),nmax_res)):
        res_alpha = res_alpha + alpha_step

        # consider only cases with minimum energy gap between vibrational states
        ev1 = wf.eigenvalue(1,De_res,res_alpha,mu)
        ev0 = wf.eigenvalue(0,De_res,res_alpha,mu)
        #print "<<<>>>"
        #print ev1-ev0
        #print minEdiff
        #print "<<<>>>"
        if ((ev1-ev0) < minEdiff):
            continue

        #print "positive case"
        #print minEdiff
        #print ev1-ev0
        #print "---------------------------------------------------"
            
        # consider only cases, where the overlap of both res vib states with the
        # ground state is sufficient
        FCres1 = wf.FC(1,res_alpha,Req,De_res,mu,
                     0,gs_a,gs_Req,gs_de,R_min,R_max)
        FCres0 = wf.FC(0,res_alpha,Req,De_res,mu,
                     0,gs_a,gs_Req,gs_de,R_min,R_max)
        if not (FCres1 >= FCmin and FCres0 >= FCmin):
            continue

        #print FCres0, FCres1
        #print "---------------------------------------------------"

        
        #for j in range(0,nmax_res+1):
        #    FC_res_gs = wf.FC(j,res_alpha,Req,De_res,mu,
        #                 0,gs_a,Req,gs_de,R_min,R_max)

        De_fin = De_min
        while (De_fin < De_max):

            fin_alpha = alpha_min(mu,De_fin,nmax_fin)
            while (is_correct_nmax(lambda_param(mu,De_fin,fin_alpha+alpha_step),nmax_fin)):
                fin_alpha = fin_alpha + alpha_step
                # consider only cases, where the overlap of both res vib states with the
                # ground state is sufficient
                FCfin1 = wf.FC(1,res_alpha,Req,De_res,mu,
                             0,fin_alpha,Req,De_fin,R_min,R_max)
                FCfin0 = wf.FC(0,res_alpha,Req,De_res,mu,
                             0,fin_alpha,Req,De_fin,R_min,R_max)
                if not (FCfin1 >= FCmin and FCfin0 >= FCmin):
                    continue

                print 'resonant', De_res, res_alpha
                print 'final   ', De_fin, fin_alpha
                print FCres0, FCres1
                print FCfin0, FCfin1
                print sc.hartree_to_ev(ev1-ev0)
                print '---------------------------------------------------'
                
            

            De_fin = De_fin + De_step


    De_res = De_res + De_step






