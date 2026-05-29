from sympy import nsolve, var
import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol

def SepDistance(values):
	sep_output = np.zeros([values.shape[0]])
	h_output = np.zeros([values.shape[0]])
	Pp_output = np.zeros([values.shape[0]])
	NPRd = 37.5
	Po1 = 100000.0*37.5
	xt = 0.1
	
	for i in range(len(sep_output)):
		error = 1
		tol = 10**-3
		M1guess = 2.3
		ls_old = 0
		gamma1 = values[i][0]
		gammaj = values[i][1]
		NPR = values[i][2]
		SPR = values[i][3]
		Cd = values[i][4]
		alpha = values[i][5]
		xm = values[i][6]
		rth = values[i][7]
		Dj = values[i][8]
		print i, Cd, rth, Dj
		if NPR >= NPRd:
			Pav = Po1/NPRd
		if NPR < NPRd:
			Pav = Po1/NPR

		Poj = SPR*Po1
		Pj = Pav
		while np.abs(error) > tol:
			xs = xm - ls_old
			rs = rth + xs*np.tan(alpha*np.pi/180)
			area_ratio = (rs/rth)**2
			M = var('M')

			# Solving A/A* for M
			eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio
			M1 = float(nsolve(eq, M, M1guess))

			# Getting dynamic pressure, assuming constant Po1 before the shock
			P1 = Po1/((1 + 0.5*(gamma1-1)*M1**2)**(gamma1/(gamma1-1)))
			
			# calculating Cp; these dummy variables separate the computation for visualization
			num = ((gamma1+1)**2)*(M1**2)
			den = 4*gamma1*(M1**2) - 2*(gamma1-1)
			dummy1 = num/den
			dummy2 = (2*gamma1*M1**2 - (gamma1-1))/(gamma1+1)
			Cp = (dummy2*(dummy1**(gamma1/(gamma1-1))) - 1)*(2./(gamma1*(M1**2)))

			q1 = 0.5*gamma1*P1*(M1**2)
			term1 = (1./(gammaj**2 - 1))*(1 - (Pj/Poj)**((gammaj-1)/gammaj))
			
			h = ((Cd*(Dj**2)*gammaj*Poj*((2./(gammaj+1))**(1./(gammaj-1)))*term1**0.5)/(P1 - Pav + 0.5*Cp*q1))**0.5

			rise = (1 + 0.5*M1)**-1
			# rise = 0.582*(1 + 0.5*(gamma1 - 1)*M1**2)**(-0.1197*gamma1/(gamma1-1))
			
			Pp = P1*1./rise
			
			# oblique shock relations
			beta = np.arcsin((((1./rise)*(gamma1+1) + (gamma1-1))/(2*gamma1*M1**2))**0.5)
			epsilon = np.arctan((M1**2*np.sin(2*beta) - 2/np.tan(beta))/(2 + M1**2*(gamma1 + np.cos(2*beta))))
			ls_new = h*(1./np.sin(epsilon) - 1)

			error = ls_new - ls_old
			ls_old = ls_new

		sep_output[i] = xs
		h_output[i] = h
		Pp_output[i] = Pp

	return sep_output, h_output, Pp_output

