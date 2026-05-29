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

def Forces(values):
	delta_output = np.zeros([values.shape[0]])
	Nxs = 50
	NPRd = 37.5
	Po1 = 100000.0*37.5

	for i in range(len(delta_output)):
		print i
		Fxb = 0
		Fxd = 0
		Fyb = 0 
		Fyd = 0
		M1guess = 2.3

		gamma1 = values[i][0]
		gammaj = values[i][1]
		NPR = values[i][2]
		SPR = values[i][3]
		Cd = values[i][4]
		alpha = values[i][5]
		xm = values[i][6]
		rth = values[i][7]
		Dj = values[i][8]
		xs = values[i][9]
		psi_max = values[i][10]
		Pp = values[i][11]
		xd = values[i][12]
		if NPR >= NPRd:
			Pav = Po1/NPRd
		if NPR < NPRd:
			Pav = Po1/NPR

		Poj = SPR*Po1
		Aj = Cd*(np.pi/4.)*Dj**2
		Pc = Po1/((1 + 0.5*(gamma1-1))**(gamma1/(gamma1-1)))
		Pj = Poj/((1 + 0.5*(gammaj-1))**(gammaj/(gammaj-1)))
		Fc = np.pi*rth**2*(1 + gamma1)*Pc
		Fxj = Aj*(1 + gammaj)*Pj*np.sin(alpha*np.pi/180)
		Fyj = Aj*(1 + gammaj)*Pj*np.cos(alpha*np.pi/180)

		x_edge = np.linspace(xs, xm, Nxs)
		radius_edge = np.zeros(Nxs)
		x_down = np.linspace(xm, xd, Nxs)
		radius_down = np.zeros(Nxs)
		
		# populating the radius at each x location
		for j in range(len(x_edge)):
			radius_edge[j] = rth + x_edge[j]*np.tan(alpha*np.pi/180)  

		# integrating the forces in the x-direction	
		for j in range(1, len(x_edge)):
		 	area_ratio = (radius_edge[j]/rth)**2
			M = var('M')
			# Solving A/A* for M
			eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio
			M1_edge = float(nsolve(eq, M, M1guess))
			# Getting dynamic pressure, assuming constant Po1 before the shock
			P1_edge = Po1/((1 + 0.5*(gamma1-1)*M1_edge**2)**(gamma1/(gamma1-1)))

			Pb = (x_edge[j] - xs)*(-Pp + P1_edge)/(xm - xs) + (Pp - P1_edge)
			angle = (x_edge[j] - xs)*(psi_max)/(xm - xs) 

			Fxb += 2*Pb*angle*radius_edge[j]*(x_edge[j]-x_edge[j-1])*np.sin(alpha*np.pi/180)
			Fyb += 2*Pb*angle*radius_edge[j]*(x_edge[j]-x_edge[j-1])*np.cos(alpha*np.pi/180)

		# OVER-EXPANDED CASE
		if NPR < NPRd: 
			# populating the radius at each x location
			for k in range(len(x_down)):
				radius_down[k] = rth + x_down[k]*np.tan(alpha*np.pi/180)
			# integrating the forces in the x-direction	
			for k in range(1, len(x_down)):
			 	area_ratio = (radius_down[k]/rth)**2
				M = var('M')
				# Solving A/A* for M
				eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio
				M1_down = float(nsolve(eq, M, M1guess))
				# Getting dynamic pressure, assuming constant Po1 before the shock
				P1_down = Po1/((1 + 0.5*(gamma1-1)*M1_down**2)**(gamma1/(gamma1-1)))

				Fxd += 2*(Pav-P1_down)*psi_max*radius_down[k]*(x_down[k] - x_down[k-1])*np.sin(alpha*np.pi/180)
				Fyd += 2*(Pav-P1_down)*psi_max*radius_down[k]*(x_down[k] - x_down[k-1])*np.cos(alpha*np.pi/180)


		Fx_total = Fc + Fxj + Fxb + Fxd
		Fy_total = Fyj + Fyb + Fyd
		delta_output[i] = np.arctan(Fy_total/Fx_total)

	return delta_output


