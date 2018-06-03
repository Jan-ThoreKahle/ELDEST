#!/usr/bin/python

import scipy
import scipy.integrate as integrate
import numpy as np
import sciconv
import complex_integration as ci
import analytic_integrals as ai
import in_out


f = lambda x: np.exp(2*x)
g = lambda y: np.exp(1j*y)

I = ci.complex_double_quadrature(f,g,0,np.pi,lambda x: 0, lambda x: x)
print I
print "Result: ((107.29833110495294+481.8424899722883j), (1.1640029546930332e-11,))"

I = ci.complex_double_quadrature(f,g,0,np.pi, 0, np.pi)
print I

print '-------------------'
f = lambda x: np.exp(2*x)
g = lambda y: np.exp(complex(2,1) * y)
I = ci.complex_double_quadrature(f,g,0,np.pi, 0, np.pi)
print I
#((-57350.06262733066+28675.031313665328j), (1.3068783974514146e-09,))
print '-------------------'

print '-------------------'
f = lambda x: np.exp(complex(2,1) * x)
g = lambda y: np.exp(-complex(1,1) * y)
I = ci.complex_double_quadrature(f,g,0,np.pi, 0, np.pi)
print I
#((-57350.06262733066+28675.031313665328j), (1.3068783974514146e-09,))
print '-------------------'

#-------------------------------------------------------
alpha = 2
beta = 1
gamma = 2
delta = 2

f = lambda x: np.exp(complex(alpha,beta)*x)
g = lambda y: np.exp(complex(gamma,delta)*y)

I = ci.complex_double_quadrature(f,g,0,np.pi,lambda x: x, lambda x: np.pi)
print I

