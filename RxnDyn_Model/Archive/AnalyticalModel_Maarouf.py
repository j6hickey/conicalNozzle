from sympy import nsolve, var
import numpy as np
import matplotlib.pyplot as plt

#############################
### INITIAL CONDITIONS ######
#############################
gamma1 = 1.4
gammaj = 1.4
Dj = 6.15*(10**-3)
Cd = 0.1
Pinfty = 8*(10**3)
NPR = 37.5
SPR = 1
xt = 100*(10**-3)
xm = 0.9*xt - Dj/2
rth = 9.72*(10**-3)
exp_ratio = 0.236
Me = 3.0
Po1 = NPR*Pinfty
Poj = SPR*Po1
alpha = 5.42

#############################
### OPERATING CONDITIONS ####
#############################
over_exp = 0
adapted = 1
under_exp = 0

if over_exp == 1:
	Pj = Pinfty

if adapted == 1 or under_exp == 1:
	Pj = Po1/((1 + 0.5*(gamma1-1)*Me**2)**(gamma1/(gamma1-1)))

Pav = Pj
error = 1

tol = 10**-6

M1guess = 2.3
ls_old = 0
Zukoski = 0
Schillig = 1
i = 0
j = 0

#############################
###### SEP DISTANCE #########
#############################
xs = xm
xsprime = xm	
while np.abs(error) > tol:
	errorp = 1
	# rs = rth + rth*(1./(exp_ratio**0.5) - 1)*length_ratio
	rs = rth + xs*np.tan(alpha*np.pi/180)
	area_ratio = (rs/rth)**2
	M = var('M')
	eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio
	M1 = float(nsolve(eq, M, M1guess))
	P1 = Po1/((1 + 0.5*(gamma1-1)*M1**2)**(gamma1/(gamma1-1)))
	
	num = ((gamma1+1)**2)*(M1**2)
	den = 4*gamma1*M1**2 - 2*(gamma1-1)
	dummy1 = num/den
	dummy2 = (2*gamma1*M1**2 - (gamma1-1))/(gamma1+1)
	Cp = (dummy2*dummy1**(gamma1/(gamma1-1)) - 1)*(2./(gamma1*M1**2))
	term1 = (1./(gammaj**2 - 1))*(1 - (Pj/Poj)**((gammaj-1)/gammaj))
	q1 = 0.5*gamma1*P1*M1**2
	h = ((2*Cd*Dj**2*gammaj*Poj*((2./(gammaj+1))**(1./(gammaj-1)))*term1**0.5)/(2*(P1 - Pav) + Cp*q1))**0.5

	while np.abs(errorp) > tol:	
		rsprime = rth + xsprime*np.tan(alpha*np.pi/180)
		area_ratiop = (rsprime/rth)**2
		eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratiop
		M1prime = float(nsolve(eq, M, M1guess))

		Psep = 0.541*P1*(Po1/Pinfty)**0.136
	
		if Schillig == 1:
			rise = 0.582*(1 + 0.5*(gamma1 - 1)*M1prime**2)**(-0.1197*gamma1/(gamma1-1))
		if Zukoski == 1:
			rise = (1 + 0.5*M1prime)**-1
		Pp = P1*1./rise
	
		beta = np.arcsin((((1./rise)*(gamma1+1) + (gamma1-1))/(2*gamma1*M1prime**2))**0.5)
		epsilon = np.arctan((M1prime**2*np.sin(2*beta) - 2/np.tan(beta))/(2 + M1prime**2*(gamma1 + np.cos(2*beta))))
		# epsilon = np.arctan(1./((0.5*(gamma1+1)*M1**2/(M1**2*np.sin(beta)*np.sin(beta)-1)-1)*np.tan(beta)))
		ls_new = h*(1./np.sin(epsilon) - 1)
		xsprime = xm - ls_new
		errorp = ls_new - ls_old
		ls_old = ls_new
		j += 1 
		print j
		print "M1p  | Psep  |  Pp   | xsp/xt | Beta | Eps. |   h   |"
		print '{:.1f}'.format(round(M1prime,1)),"|", '{:.3f}'.format(round(P1/Po1,3)),"|",'{:.3f}'.format(round(Pp/Po1,3)),"|",'{:.3f}'.format(round(xsprime/xt,3)),"|", \
		'{:.1f}'.format(round(beta*180/np.pi,1)),"|",'{:.1f}'.format(round(epsilon*180/np.pi,1)),"|",'{:.3f}'.format(round(h,3)),"|"
		print "------------------------------------------------------"

	error = xs - xsprime
	xs = xsprime
	xp = xs + h
	i += 1
	print i
	print "M1  | Psep  |  Pp   | xs/xt | xp/xt | Beta | Eps. |   h   |"
	print '{:.1f}'.format(round(M1,1)),"|", '{:.3f}'.format(round(P1/Po1,3)),"|",'{:.3f}'.format(round(Pp/Po1,3)),"|",'{:.3f}'.format(round(xs/xt,3)),"|", \
	'{:.3f}'.format(round(xp/xt,3)),"|", '{:.1f}'.format(round(beta*180/np.pi,1)),"|",'{:.1f}'.format(round(epsilon*180/np.pi,1)),"|",'{:.3f}'.format(round(h,3)),"|"
	print "------------------------------------------------------"

#############################
#######  BOW SHOCK ##########
#############################
Rc = h*1.143*np.exp(0.54/(M1 - 1)**1.2)
Delta = h*0.143*np.exp(3.24/M1**2)
mu = np.arcsin(1/M1)

Rs = ls_new + Rc - Delta

cotan = 1/np.tan(mu)

y = np.linspace(0, Rs, 50)
x = h + Delta - Rs*cotan**2*((1 + y**2*(np.tan(mu))**2/Rs**2)**0.5 - 1)

y1 = np.linspace(0, h, 50)
x1 = (h**2 - y1**2)**0.5 

# y2 = [h]*50
# x2 = np.linspace(h + Delta - Rs*cotan**2*((1 + (np.tan(mu))**2)**0.5 - 1), 0, 50)

# y3 = np.linspace(0, Rc, 50)
# x3 = h + Delta - Rc*cotan**2*((1 + y3**2*(np.tan(mu))**2/Rc**2)**0.5 - 1)

# plt.figure()
# plt.xlim(0, 0.002)
# plt.plot(x, y)
# plt.plot(x1, y1)
# plt.plot(x2, y2)
# plt.plot(x3, y3)
# plt.show()

#############################
#######  FORCES #############
#############################

# Dynalpy Flux (due to flow onto throat)
Pc = Po1/((1 + 0.5*(gamma1-1))**(gamma1/(gamma1-1)))
Pj = Poj/((1 + 0.5*(gammaj-1))**(gammaj/(gammaj-1)))
Fc = np.pi*rth**2*(1 + gamma1)*Pc
Fjx = Cd*np.pi*(Dj/2)**2*(1 + gammaj)*Pj*np.sin(alpha*np.pi/180)
Fjy = Cd*np.pi*(Dj/2)**2*(1 + gammaj)*Pj*np.cos(alpha*np.pi/180)

print Fjy

