from sympy import nsolve, var
import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol


def Forces(values):
	delta_output = np.zeros([values.shape[0]])
	Nxs = 50
	NPRd = 37.5

	for i in range(len(delta_output)):
		print i
		Fxb = 0
		Fxd = 0
		Fyb = 0 
		Fyd = 0
		M1guess = 2.3

		rth = values[i][0]
		re = values[i][1]
		Dj = values[i][2]
		Po1 = values[i][3]
		Poj = values[i][4]
		gamma1 = values[i][5]
		gammaj = values[i][6]
		xm = values[i][7]
		xt = values[i][8]
		Cd = values[i][9]
		xs = values[i][10]
		xd = values[i][11]
		Pp = values[i][12]
		psi_max = np.pi*values[i][13]/180.

		NPR = 15.0
		alpha = np.arctan((re - rth)/xt)
		
		if NPR >= NPRd:
			Pav = Po1/NPRd
		if NPR < NPRd:
			Pav = Po1/NPR

		Aj = Cd*(np.pi/4.)*Dj**2
		Pc = Po1/((1 + 0.5*(gamma1-1))**(gamma1/(gamma1-1)))
		Pj = Poj/((1 + 0.5*(gammaj-1))**(gammaj/(gammaj-1)))
		Fc = np.pi*rth**2*(1 + gamma1)*Pc
		Fxj = Aj*(1 + gammaj)*Pj*np.sin(alpha)
		Fyj = Aj*(1 + gammaj)*Pj*np.cos(alpha)

		x_edge = np.linspace(xs, xm, Nxs)
		radius_edge = np.zeros(Nxs)
		
		# populating the radius at each x location
		for j in range(len(x_edge)):
			radius_edge[j] = rth + x_edge[j]*np.tan(alpha)  

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

			Fxb += 2*Pb*angle*radius_edge[j]*(x_edge[j]-x_edge[j-1])*np.sin(alpha)
			Fyb += 2*Pb*angle*radius_edge[j]*(x_edge[j]-x_edge[j-1])*np.cos(alpha)

		# OVER-EXPANDED CASE
		if NPR < NPRd: 
			x_down = np.linspace(xm, xd, Nxs)
			radius_down = np.zeros(Nxs)
			# populating the radius at each x location
			for k in range(len(x_down)):
				radius_down[k] = rth + x_down[k]*np.tan(alpha)
			# integrating the forces in the x-direction	
			for k in range(1, len(x_down)):
			 	area_ratio = (radius_down[k]/rth)**2
				M = var('M')
				# Solving A/A* for M
				eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio
				M1_down = float(nsolve(eq, M, M1guess))
				# Getting dynamic pressure, assuming constant Po1 before the shock
				P1_down = Po1/((1 + 0.5*(gamma1-1)*M1_down**2)**(gamma1/(gamma1-1)))

				Fxd += 2*(Pav-P1_down)*psi_max*radius_down[k]*(x_down[k] - x_down[k-1])*np.sin(alpha)
				Fyd += 2*(Pav-P1_down)*psi_max*radius_down[k]*(x_down[k] - x_down[k-1])*np.cos(alpha)


		Fx_total = Fc + Fxj + Fxb + Fxd
		Fy_total = Fyj + Fyb + Fyd
		delta_output[i] = np.arctan(Fy_total/Fx_total)

	return delta_output


