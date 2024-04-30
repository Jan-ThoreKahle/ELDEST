#!/usr/bin/python

##########################################################################
#                             WELLENFKT                                  #
##########################################################################
# Purpose: Determine Franck-Condon factors and eigenenergies of Morse    #
#          potentials based on analytical solutions                      #
#          of Morse-Potentials (JCP 88, 4535 (1988).)                    #
##########################################################################
# written by: Elke Fasshauer August 2019                                 #
##########################################################################

import scipy.integrate as integrate
#from scipy.special import gamma
from scipy.special import factorial
#from scipy.special import eval_genlaguerre
#from scipy.special import genlaguerre
import numpy as np
import sciconv as sc
#import in_out
#import sys
#import warnings
import potentials
import complex_integration as ci
from mpmath import coulombf, coulombg

#-------------------------------------------------------------------------

#n = 0
#n1 = n
#n2 = n

#-------------------------------------------------------------------------
# Parameters for potentials
#-------------------------------------------------------------------------
#Bedenke, dass diese Parameter Grottenschlecht sind
#gs_de      =  0.0001107257
#gs_a       =  1.105308
#gs_Req     =  5.918877
#gs_const   = -0.00018536
#
#u_de       =  0.096727
#u_a        =  0.447152
#u_Req      =  4.523888
#u_const    = -256.233965
#
#g_de       = 0.079073
#g_a        = 0.143584
#g_Req      = 10.003604
#g_const    = -256.182536
#
#mass1 = 20.1797
#mass2 = 20.1797
#-------------------------------------------------------------------------

#red_mass_gmol = mass1 * mass2 / (mass1 + mass2)
#red_mass = sc.gmol_to_me(red_mass_gmol)
#print "redmass = ", red_mass

def red_mass_au(mass1,mass2):
    red_mass_gmol = mass1 * mass2 / (mass1 + mass2)
    red_mass_au = sc.gmol_to_me(red_mass_gmol)
    return red_mass_au


def eigenvalue(n, De, alpha, red_mass):
    En =   (n + 0.5) * alpha * np.sqrt(2*De / red_mass) \
         - (n + 0.5)**2 * alpha**2 / (2* red_mass)
    return En

def sqrt_fact(real):
    if (np.absolute(1 - real) < 1.0E-7):
        return np.sqrt(real)
    elif real < 1.0:
        return np.sqrt(factorial(real) )
    else:
        return np.sqrt(real) * sqrt_fact(real-1)

#print "Grundzustand:"
#print "Potentialtiefe", gs_de
#print eigenvalue(n, gs_de, gs_a, red_mass)
#
#print "--------------------------------"



#####
# Provide everything in atomic units from here on!
# I. e. all energies (De, V_hyp_b) in Hartree, all lengths (R, Req, Rmin, Rmax) in Bohr,
# all inverse lengths (alpha) in inverse Bohr, all masses (red_mass) in electron masses,
# all composite quantities (V_hyp_a) in the respective atomic units.
#####

def const_s_psi(R,n,s,alpha,Req,lambda_param):
    z = 2* lambda_param * np.exp(-alpha * (R - Req))
    if (n == 0):
        psi_0 = ( 1.0 
                     * np.sqrt(alpha) # not mentioned part of norm. fact.
                     # needed due to different type of Laguerre polyomials
                     #/ sqrt_fact(2*lambda_param-2) alternative normalization
                     # factor based on Eq. (41)
                     * np.sqrt(s) * sqrt_fact(n) / sqrt_fact(s+n)
                     * z**(s/4) #improves numerical stability to split
                     * np.exp(-z / 2)
                     * z**(s/4)
                     )
        return psi_0
    elif (n == 1):
        psi_1 = ( 1.0 
                     * (2*n + s -1 -z)
                     * np.sqrt(1./(n*(s + n)))
                     )
        return psi_1 * const_s_psi(R,n-1,s,alpha,Req,lambda_param)
    else:
        prefac  =  np.sqrt(1./(n*(s + n)))
        prefac1 =  (2 * n + s -1 - z)
        prefac2 = np.sqrt((n-1) * (n + s - 1))
        return prefac * (  prefac1 * const_s_psi(R,n-1,s,alpha,Req,lambda_param)
                         - prefac2 * const_s_psi(R,n-2,s,alpha,Req,lambda_param)  )


def psi_n(R,n,alpha,Req,red_mass,De):
    lambda_param = np.sqrt(2*red_mass*De) / alpha
    s = 2*lambda_param - 2*n - 1
    psi = const_s_psi(R,n,s,alpha,Req,lambda_param)
    return psi
    

def FC(n1,alpha1,Req1,De1,red_mass,n2,alpha2,Req2,De2,R_min,R_max):
    func = lambda R: (np.conj(psi_n(R,n1,alpha1,Req1,red_mass,De1))
                            * psi_n(R,n2,alpha2,Req2,red_mass,De2) )
    tmp = integrate.quad(func, R_min, R_max)
    FC = tmp[0]
    return FC


def psi_freehyp(R,a,b,red_mass,R_start,phase=0):    # model: free particle with energy corresponding to a point (at R_start) on a hyperbola, psi = 0 for section left of R_start
    a_eV = sc.hartree_to_ev(a)
    b_eV = sc.hartree_to_ev(b)
    E_au = potentials.hyperbel(a_eV,b_eV,R_start) - b 
    if (R <= R_start):
        psi = 0
    else:
        K_au = np.sqrt(2 * red_mass * E_au)              # momentum of the nuclear system
        norm = np.sqrt(red_mass / (2 * np.pi * K_au))    # normalization factor for energy-normalized plane wave
        psi = norm * np.exp(1.j * (K_au * (R - R_start) + phase))
    return psi


def FCmor_freehyp(n1,alpha1,Req1,De1,red_mass,V2a,V2b,R_start,R_min,R_max,**kwargs):
    lim = kwargs.get("limit", 50)
    phi = kwargs.get("phase", 0)
    func = lambda R: (np.conj(psi_n(R,n1,alpha1,Req1,red_mass,De1))
                            * psi_freehyp(R,V2a,V2b,red_mass,R_start,phi) )
    tmp = ci.complex_quadrature(func, R_min, R_max, limit=lim)
    FC = tmp[0]
    return FC


def FCmor_freehyp_R6(n1,alpha1,Req1,De1,red_mass,V2a,V2b,R_start,R_min,R_max,**kwargs):
    lim = kwargs.get("limit", 50)
    phi = kwargs.get("phase", 0)
    func = lambda R: (np.conj(psi_n(R,n1,alpha1,Req1,red_mass,De1))
                            * psi_freehyp(R,V2a,V2b,red_mass,R_start,phi) * R**(-3) )
    tmp = ci.complex_quadrature(func, R_min, R_max, limit=lim)
    FC = tmp[0]
    return FC


def FCfreehyp_freehyp(V1a,V1b,R_start1,red_mass,V2a,V2b,R_start2,R_min,R_max,**kwargs):
    lim = kwargs.get("limit", 50)
    phi1 = kwargs.get("phase1", 0)
    phi2 = kwargs.get("phase2", 0)
    func = lambda R: (np.conj(psi_freehyp(R,V1a,V1b,red_mass,R_start1,phi1))
                            * psi_freehyp(R,V2a,V2b,red_mass,R_start2,phi2) )
    tmp = ci.complex_quadrature(func, R_min, R_max, limit=lim)
    FC = tmp[0]
    return FC


def psi_free(R,a,b,red_mass,R_start,phase=0):    # model: free particle with energy corresponding to a point (at R_start) on a hyperbola
    a_eV = sc.hartree_to_ev(a)
    b_eV = sc.hartree_to_ev(b)
    E_au = potentials.hyperbel(a_eV,b_eV,R_start) - b 
    K_au = np.sqrt(2 * red_mass * E_au)              # momentum of the nuclear system
    norm = np.sqrt(red_mass / (2 * np.pi * K_au))    # normalization factor for energy-normalized plane wave
    psi = norm * np.exp(1.j * (K_au * (R - R_start) + phase))
    return psi

def norm_free(R,a,b,red_mass,R_start,**kwargs):
    lim = kwargs.get("limit", 50)
    phi = kwargs.get("phase", 0)
    func = lambda R: (np.conj(psi_free(R,a,b,red_mass,R_start,phi))
                            * psi_free(R,a,b,red_mass,R_start,phi) )
    tmp = ci.complex_quadrature(func, R_min, R_max, limit=lim)
    FC = tmp[0]
    return FC


def psi_hyp(R,a,b,red_mass,R_start):        # model: particle in a hyperbolic potential a/R + b, energy chosen to be a/R_start
    a_eV = sc.hartree_to_ev(sc.bohr_to_angstrom(a))
    b_eV = sc.hartree_to_ev(b)
    E_au = potentials.hyperbel(a_eV,b_eV,R_start) - b
    K_au = np.sqrt(2 * red_mass * E_au)             # momentum of the nuclear system
    norm = np.sqrt(2 * red_mass / (np.pi * K_au))   # normalization factor for energy-normalized plane wave (F_kL = sqrt(2/pi) * F_L is k-normalized, F_EL = sqrt(m/k) * F_kL)
    eta = a * red_mass / K_au
    z = K_au * R
    func = coulombf(l = 0, eta = eta, z = z, maxterms=10**6)      # so that psi->sin[K*Theta(x)] for x->inf and regular (for x->0)
    #func = coulombg(l = 0, eta = eta, z = z) + 1.j * coulombf(l = 0, eta = eta, z = z)      # lin comb chosen so that psi->exp[iKx] for x->inf (up to a constant phase shift)
    psi = norm * func
    return float(psi)


def FCmor_hyp(n1,alpha1,Req1,De1,red_mass,V2a,V2b,R_start,R_min,R_max,**kwargs):
    lim = kwargs.get("limit", 50)
    func = lambda R: (np.conj(psi_n(R,n1,alpha1,Req1,red_mass,De1))
                            * psi_hyp(R,V2a,V2b,red_mass,R_start) )
    tmp = ci.complex_quadrature(func, R_min, R_max, limit=lim)
    FC = tmp[0]
    return FC


def FCmor_hyp_R6(n1,alpha1,Req1,De1,red_mass,V2a,V2b,R_start,R_min,R_max,**kwargs):
    lim = kwargs.get("limit", 50)
    func = lambda R: (np.conj(psi_n(R,n1,alpha1,Req1,red_mass,De1))
                            * psi_hyp(R,V2a,V2b,red_mass,R_start) * R**(-3) )
    tmp = ci.complex_quadrature(func, R_min, R_max, limit=lim)
    FC = tmp[0]
    return FC


def FCmor_mor_R6(n1,alpha1,Req1,De1,red_mass,n2,alpha2,Req2,De2,R_min,R_max):
    func = lambda R: (np.conj(psi_n(R,n1,alpha1,Req1,red_mass,De1))
                            * psi_n(R,n2,alpha2,Req2,red_mass,De2) * R**(-3) )
    tmp = integrate.quad(func, R_min, R_max)
    FC = tmp[0]
    return FC

    
#R_min = sc.angstrom_to_bohr(1.5)
#R_max = sc.angstrom_to_bohr(30.0)
#
#print "R_min = ", R_min
#print "R_max = ", R_max
#
#FCO = FC(n1,u_a,u_Req,u_de,red_mass,n2,u_a,u_Req,u_de,R_min,R_max)
#
#print "FC = ", FCO



