f = lambda x: np.exp(2*x)
g = lambda y: np.exp(1j*y)

I = complex_double_quadrature(f,g,0,np.pi,lambda x: 0, lambda x: x)

Result: 
((107.29833110495294+481.8424899722883j), (1.1640029546930332e-11,))
