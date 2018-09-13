#!/usr/bin/python

##########################################################################
#                    PULSE DEFINITION ROUTINES                           #
#        Investigating Electronic Decay Processes with Streaking         #
##########################################################################
# Purpose:                                                               #
#          - Handle definitions of pulses in one place                   #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer May 2018                                    #
##########################################################################

import numpy as np

##-------------------------------------------------------------------------
##Variante mit TX
#f_TX = lambda tau: 1./4 * ( np.exp(2j * np.pi * (TX_au/2 - tau) / TX_au)
#                      + 2
#                      + np.exp(-2j * np.pi * (TX_au/2 - tau) /TX_au) )
#
#fp_TX = lambda tau: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (TX_au/2 - tau)/TX_au)
#                                     + np.exp(-2j*np.pi* (TX_au/2 - tau) /TX_au) )
#
#FX_TX = lambda tau: - A0X * np.cos(Omega_au * (TX_au/2 - tau)) * fp_TX(tau) + A0X * Omega_au * np.sin(Omega_au * (TX_au/2 - tau)) * f_TX(tau)
#
## functions for the norm
#if (X_sinsq):
#    print 'use sinsq function'
#    f_t1  = lambda t1: 1./4 * ( np.exp(2j * np.pi * t1 / TX_au)
#                          + 2
#                          + np.exp(-2j * np.pi * t1 /TX_au) )
#    
#    fp_t1 = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* t1/TX_au)
#                                         + np.exp(-2j*np.pi* t1 /TX_au) )
#elif (X_gauss):
#    print 'use gauss function'
#    f_t1  = lambda t1: ( 1./ np.sqrt(2*np.pi * sigma**2)
#                       * np.exp(-t1**2 / (2*sigma**2)))
#    fp_t1 = lambda t1: ( -t1 / np.sqrt(2*np.pi) / sigma**3
#                       * np.exp(-t1**2 / (2*sigma**2)))
#else:
#    print 'no pulse shape selected'
#
#FX_t1 = lambda t1: (- A0X * np.cos(Omega_au * t1) * fp_t1(t1)
#                    + A0X * Omega_au * np.sin(Omega_au * (t1)) * f_t1(t1)
#                   )

## IR pulse
##-------------------------------------------------------------------------
#A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) / TL_au)**2 \
#                      * np.cos(omega_au * t3 + phi)
#integ_IR = lambda t3: - 1j*p_au * A_IR(t3)
##IR_after = lambda t1: ci.complex_romberg(integ_IR, t1, delta_t_au + TL_au/2)
#
##integ_IR = lambda t3: -1j*p_au**2/2 - 1j*p_au * A_IR(t3)
##IR_after = lambda t1: ci.complex_romberg(integ_IR, t1, delta_t_au + TL_au/2) \
##                      -1j * p_au**2/2 * (t_au - delta_t_au - TL_au/2) 
#
#IR_during = lambda t1:  np.exp(-1j * p_au**2/2 * (t_au - t1)) \
#                        * np.exp( -1j * p_au * A0L / 4
#                        * (np.sin(2*np.pi/TL_au * (t_au - delta_t_au) - omega_au * t_au
#                                  - phi)
#                            / (2*np.pi/TL_au - omega_au)
#                           + np.sin(-2*np.pi/TL_au * (t1 - delta_t_au) + omega_au * t1
#                                  + phi) 
#                            / (2*np.pi/TL_au - omega_au)
#                           + np.sin(2*np.pi/TL_au * (t_au - delta_t_au) + omega_au * t_au
#                                  + phi) 
#                            / (2*np.pi/TL_au + omega_au)
#                           + np.sin(-2*np.pi/TL_au * (t1 - delta_t_au) - omega_au * t1
#                                  - phi) 
#                            / (2*np.pi/TL_au + omega_au)
#                           + 2./omega_au * np.sin(omega_au * t_au + phi)
#                           - 2./omega_au * np.sin(omega_au * t1 + phi)
#                          )
#                       )
#
##IR_after = lambda t1:  np.exp(-1j * p_au**2/2 * (t_au - t1)) \
##                       * np.exp( -1j * p_au * A0L / 4
##                       * (np.sin(np.pi - omega_au * (delta_t_au + TL_au/2)
##                                 - phi)
##                           / (2*np.pi/TL_au - omega_au)
##                          + np.sin(-2*np.pi/TL_au * (t1 - delta_t_au) + omega_au * t1
##                                 + phi) 
##                           / (2*np.pi/TL_au - omega_au)
##                          + np.sin(np.pi + omega_au * (delta_t_au + TL_au/2)
##                                 + phi) 
##                           / (2*np.pi/TL_au + omega_au)
##                          + np.sin(-2*np.pi/TL_au * (t1 - delta_t_au) - omega_au * t1
##                                 - phi) 
##                           / (2*np.pi/TL_au + omega_au)
##                          + 2./omega_au * np.sin(omega_au * (delta_t_au + TL_au/2) + phi)
##                          - 2./omega_au * np.sin(omega_au * t1 + phi)
##                         )
##                      )
#
##IR_after = lambda t1:  -1j * p_au**2/2 * (t_au - t1) \
##                        -1j * p_au * A0L / 4 \
##                       * (np.sin(np.pi - omega_au * (delta_t_au + TL_au/2)
##                                 - phi)
##                           / (2*np.pi/TL_au - omega_au)
##                          + np.sin(-2*np.pi/TL_au * (t1 - delta_t_au) + omega_au * t1
##                                 + phi) 
##                           / (2*np.pi/TL_au - omega_au)
##                          + np.sin(np.pi + omega_au * (delta_t_au + TL_au/2)
##                                 + phi) 
##                           / (2*np.pi/TL_au + omega_au)
##                          + np.sin(-2*np.pi/TL_au * (t1 - delta_t_au) - omega_au * t1
##                                 - phi) 
##                           / (2*np.pi/TL_au + omega_au)
##                          + 2./omega_au * np.sin(omega_au * (delta_t_au + TL_au/2) + phi)
##                          - 2./omega_au * np.sin(omega_au * t1 + phi)
##                         )

def IR_after(t1, p_au,
             A0L, omega_au , delta_t_au, TL_au, phi):
    IR_after = -1j * p_au * A0L / 4 \
               * (np.sin(np.pi - omega_au * (delta_t_au + TL_au/2)
                         - phi)
                   / (2*np.pi/TL_au - omega_au)
                  + np.sin(-2*np.pi/TL_au * (t1 - delta_t_au) + omega_au * t1
                         + phi) 
                   / (2*np.pi/TL_au - omega_au)
                  + np.sin(np.pi + omega_au * (delta_t_au + TL_au/2)
                         + phi) 
                   / (2*np.pi/TL_au + omega_au)
                  + np.sin(-2*np.pi/TL_au * (t1 - delta_t_au) - omega_au * t1
                         - phi) 
                   / (2*np.pi/TL_au + omega_au)
                  + 2./omega_au * np.sin(omega_au * (delta_t_au + TL_au/2) + phi)
                  - 2./omega_au * np.sin(omega_au * t1 + phi)
                 )
    return IR_after

###-------------------------------------------------------------------------
### technical defintions of functions
##
###direct ionization
##fun_t_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
##                                   * np.exp(1j * p_au**2/2 * (t1-t_au))
##fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
##                                   * np.exp(1j * p_au**2/2 * (t1-TX_au/2))
##
##dress_I = lambda t1: integrate.quad(integ_IR,t1,t_au)[0]
##dress = lambda t1: np.exp(-1j/2 * dress_I(t1))
##
##dress_I_after = lambda t1: integrate.quad(integ_IR,t1,(delta_t_au + TL_au/2))[0]
##dress_after = lambda t1: np.exp(-1j/2 * dress_I_after(t1))
##fun_dress_after = lambda t1: (FX_t1(t1)
##                              * np.exp(1j * E_fin_au * (t1-t_au)) \
##                              * IR_after(t1)
##                             )
##
##fun_IR_dir = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
##                                  * IR_during(t1)
##
##
###-------------------------------------------------------------------------
## resonant state functions
#inner_prefac = lambda x,y:  np.exp(-1j * y * (p_au**2/2 + E_fin_au)) \
#                        * np.exp(-1j * p_au * A0L / (4*(2*np.pi/TL_au - omega_au))
#                                 *np.sin(2*np.pi/TL_au * (x - delta_t_au) 
#                                         - omega_au * x - phi) ) \
#                        * np.exp(-1j * p_au * A0L / (4*(2*np.pi/TL_au + omega_au))
#                                 *np.sin(2*np.pi/TL_au * (x - delta_t_au) 
#                                         + omega_au * x + phi) ) \
#                        * np.exp(-1j * p_au * A0L / omega_au
#                                 *np.sin(omega_au * x + phi) )
#
#inner_int_part = lambda x,y: 1./(complex(-np.pi * VEr_au**2, p_au**2/2 + E_fin_au - Er_au)
#                              +1j*p_au*A0L/4
#                                 * np.cos(-2*np.pi/TL_au * (x-delta_t_au) + omega_au * x + phi)
#                              +1j*p_au*A0L/4 
#                                 * np.cos(-2*np.pi/TL_au * (x-delta_t_au) - omega_au * x - phi)
#                              +1j*A0L*p_au 
#                                 * np.cos(omega_au * x + phi)
#                              ) \
#                           *(np.exp(y*(complex(-np.pi * VEr_au**2, p_au**2/2 + E_fin_au - Er_au)))
#                           *np.exp(-1j*A0L*p_au /(4*(2*np.pi/TL_au - omega_au))
#                                  * np.sin(-2*np.pi/TL_au * (x - delta_t_au) 
#                                        + omega_au * x + phi) )
#                           *np.exp(-1j*A0L*p_au /(4*(2*np.pi/TL_au + omega_au))
#                                  * np.sin(-2*np.pi/TL_au * (x - delta_t_au) 
#                                        - omega_au * x - phi) )
#                           *np.exp(1j*A0L*p_au / omega_au
#                                  * np.sin(omega_au * x + phi) )
#                           )
#
#res_inner_fun = lambda t2: np.exp(-t2 * (np.pi * VEr_au**2 + 1j*(Er_au))) \
#                           * IR_during(t2)
#
#if (integ == 'romberg'):
#    res_inner = lambda t1: ci.complex_romberg(res_inner_fun, t1, t_au)
#elif (integ == 'quadrature'):
#    res_inner = lambda t1: ci.complex.quadrature(res_inner_fun, t1, t_au)[0]
#elif (integ == 'analytic'):
#    res_inner = lambda t1: inner_prefac(t_au,t_au) * \
#                           (inner_int_part(t_au,t_au) - inner_int_part(t1,t1))
#
#
#res_outer_fun = lambda t1: FX_t1(t1) * np.exp(t1 * (np.pi* VEr_au**2 + 1j*Er_au)) \
#                           * res_inner(t1)
#
## after the pulse
#res_inner_after = lambda t2: np.exp(-t2 * (np.pi * VEr_au**2 + 1j*(Er_au))) \
#                             * IR_after(t2)
#
#if (integ == 'romberg'):
#    res_inner_a = lambda t1: ci.complex_romberg(res_inner_after, t1, t_au)
#elif (integ == 'quadrature'):
#    res_inner_a = lambda t1: ci.complex_quadrature(res_inner_after, t1, t_au)[0]
#elif (integ == 'analytic'):
#    res_inner_a = lambda t1: inner_prefac(delta_t_au + TL_au/2,t_au) * \
#                           (inner_int_part(delta_t_au + TL_au/2,t_au) - inner_int_part(t1,t1))
#
#res_outer_after = lambda t1: FX_t1(t1) * np.exp(t1 * (np.pi* VEr_au**2 + 1j*Er_au)) \
#                           * res_inner_a(t1)
