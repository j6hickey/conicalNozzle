from sympy import nsolve, var
import numpy as np
import matplotlib.pyplot as plt

gamma1 = 1.4
gammaj = 1.4
Dj = 6.15*(10**-3)
Cd = 0.95
NPR = 37.5
SPR = 1
xt = 100*(10**-3)
xm = 0.9*xt 
rth = 9.72*(10**-3)
exp_ratio = 4.24
Me = 3.0
Po1 = 300*(10**3)
Poj = SPR*Po1
Pav = Po1/NPR
alpha = 5.42

error = 1
tol = 10**-6
ls_old = 0

M1guess = 2.3
Zukoski = 0
Summerfield = 0
Schilling = 1
Schmucker = 0
Ost_Kur = 0
i = 0

if Summerfield == 1:
	P1 = 0.4*Pav

if Schilling == 1:
	P1 = 0.541*Pav*(Po1/Pav)**-0.136

M1 = ((((Po1/P1)**((gamma1-1)/gamma1)) - 1)*2./(gamma1-1))**0.5

print M1

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
	den = 4*gamma1*M1**2 - 2*(gamma1-1)
	dummy1 = num/den
	dummy2 = (2*gamma1*M1**2 - (gamma1-1))/(gamma1+1)
	Cp = (dummy2*(dummy1**(gamma1/(gamma1-1))) - 1)*(2./(gamma1*M1**2))

	# term1 is a part of h under the root
	q1 = 0.5*gamma1*P1*M1**2
	term1 = (1./(gammaj**2 - 1))*(1 - (Pj/Poj)**((gammaj-1)/gammaj))
	
	h = ((Cd*Dj**2*gammaj*Poj*((2./(gammaj+1))**(1./(gammaj-1)))*term1**0.5)/((P1 - Pav) + 0.5*Cp*q1))**0.5
	
	if Schilling == 1:
		rise = 0.582*(1 + 0.5*(gamma1 - 1)*M1**2)**(-0.1197*gamma1/(gamma1-1))
	if Zukoski == 1:
		rise = (1 + 0.5*M1)**-1
	Pp = P1*1./rise
	
	Psep = Pinfty*0.541*(Po1/Pinfty)**-0.136
	# oblique shock relations
	beta = np.arcsin((((1./rise)*(gamma1+1) + (gamma1-1))/(2*gamma1*M1**2))**0.5)
	epsilon = np.arctan((M1**2*np.sin(2*beta) - 2/np.tan(beta))/(2 + M1**2*(gamma1 + np.cos(2*beta))))
	ls_new = h*(1./np.sin(epsilon) - 1)

	error = ls_new - ls_old
	ls_old = ls_new
	# xp is not very important at the moment; it will be interpolated later
	xp = xs + h
	i += 1
	print i
	print "M1  | Psep  |  Pp   | xs/xt | xp/xt | Beta | Eps. |   h   |"
	print '{:.1f}'.format(round(M1,1)),"|", '{:.3f}'.format(round(Psep/Po1,3)),"|",'{:.3f}'.format(round(Pp/Po1,3)),"|",'{:.3f}'.format(round(xs,3)),"|", \
	'{:.3f}'.format(round(xp/xt,3)),"|", '{:.1f}'.format(round(beta*180/np.pi,1)),"|",'{:.1f}'.format(round(epsilon*180/np.pi,1)),"|",'{:.3f}'.format(round(h,3)),"|"
	print "------------------------------------------------------"
