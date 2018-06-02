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

I = ci.complex_double_quadrature(f,g,0,np.pi, 0, np.pi)
print I

print "Result: ((107.29833110495294+481.8424899722883j), (1.1640029546930332e-11,))"
