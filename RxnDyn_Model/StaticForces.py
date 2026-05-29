from sympy import nsolve, var
import numpy as np

Nxs = 50
NPRd = 37.5
NPR = 15.0

Fxb = 0
Fxd = 0
Fyb = 0 
Fyd = 0
M1guess = 2.3

rth = 9.72/1000
re = 20.13/1000
Dj = 6./1000
Po1 = 3.75*10**6
Poj = Po1
gamma1 = 1.4
gammaj = 1.4
xt = 100./1000
Cd = 0.85

if NPR == 15:
	Pp = 0.146*Po1
	psi_max = np.pi*74.3/180.
	xs = 0.581*xt
	xm = 95.4/1000
	xd = 90./1000

if NPR == 37.5:
	Pp = 0.143*Po1
	psi_max = np.pi*73.9/180.
	xd = 0.
	xm = 95.4/1000
	xs = 0.592*xt

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

print Fxd, Fyd

Fx_total = Fc + Fxj + Fxb + Fxd
Fy_total = Fyj + Fyb + Fyd
delta = np.arctan(Fy_total/Fx_total)

print delta*180/np.pi