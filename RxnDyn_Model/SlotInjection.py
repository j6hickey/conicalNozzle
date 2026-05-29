from sympy import *
import numpy as np

#############################
###### GLOBAL NUMERICS ######
#############################
error = 1
tol = 10**-3
M1guess = 2.3
ls_old = 0
nozzle_profile = "Masuya"
gamma1 = 1.4
gammaj = 1.4

#############################
### INITIAL CONDITIONS ######
#############################
if nozzle_profile == 'Wing':
	NPRd = 8.26
	Pinfty = 100000.0
	Po1 = NPRd*Pinfty
	NPR = 3.0
	SPR = 1.0

if nozzle_profile == 'Masuya':
	# This was calculated based on the
	# exp. ratio of this nozzle.
	# From Ko 2002, exp = 4.64
	NPRd = 42.1
	Pinfty = 100000.0
	Po1 = NPRd*Pinfty
	NPR = 20.0
	SPR = 1.0

if nozzle_profile == 'Sellam':
	NPRd = 37.5
	Pinfty = 100000.0
	Po1 = NPRd*Pinfty
	NPR = 37.5
	SPR = 1.0


#############################
######## GEOMETRY ###########
#############################
if nozzle_profile == 'Wing':
	# The length of the entire nozzle is 2 inches
	# from the throat.
	xt = 2.0*0.0254
	# We compute alpha using the expansion ratio:
	# 1.738 (given).
	alpha = 9.85
	# This was obtained from the area of the throat
	# in Wing 1997.
	rth = 0.0277

if nozzle_profile == 'Masuya':
	# These were all obtained from 
	# Ko 2002.
	xt = 0.09
	alpha = 9.6
	rth = 0.013

if nozzle_profile == 'Sellam':
	# These were all obtained from 
	# Sellam 2012 and Deng 2014.
	xt = 0.1
	alpha = 5.42
	rth = 9.72*10**-3


#############################
####### INJECTION ###########
#############################
if nozzle_profile == 'Wing':
	Cd = 1.0
	# The injection slot starts at 1.15 inches.
	xm = 0.575*xt
	Poj = SPR*Po1
	# Slot width is 1.375 - 1.15 inches.
	b = 0.225*0.0254
	# Psi is given as +/- 30 deg.
	psi = 60.0*np.pi/180
	Rj = rth + xm*np.tan(alpha*np.pi/180)
	Aj = Cd*b*psi*Rj

if nozzle_profile == 'Masuya':
	Cd = 0.9
	# xm 0.05 is strictly from Maarouf's thesis.
	xm = 0.05
	Poj = SPR*Po1
	# This yields a mass-flow ratio of 5%, 
	# as quoted in Maarouf Thesis.	
	b = 1.4*10**-3
	# Psi is given as +/- 30 deg.
	psi = 60.0*np.pi/180
	Rj = rth + xm*np.tan(alpha*np.pi/180)
	Aj = Cd*b*psi*Rj

if nozzle_profile == 'Sellam':
	Cd = 1.0
	xm = 0.92*xt
	Poj = SPR*Po1
	# This yields a mass-flow ratio of 5%, 
	# as quoted in Maarouf Thesis.	
	b = 2.5*10**-3
	psi = 30.0*np.pi/180
	Rj = rth + xm*np.tan(alpha*np.pi/180)
	Aj = Cd*b*psi*Rj


#############################
### OPERATING CONDITIONS ####
#############################
if NPR >= NPRd:
	Pav = Po1/NPRd
	print "Ideal or under-expanded regime."
if NPR < NPRd:
	Pav = Po1/NPR
	print "Over-expanded regime." 

# Pj is the exit pressure
# Pav is the ambient pressure
Pj = Pav


#############################
######## SEP. CRITERIA ######
#############################
Zukoski = 1
Schilling = 0
i = 0


#############################
###### SEP DISTANCE #########
#############################
print "------------------------------------------------------"
while np.abs(error) > tol:
	xs = xm - ls_old
	if xs < 0:
		xs = tol*10./i
		print ""
		print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		print "WARNING -- Entering throat zone."
		print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		print ""
	if nozzle_profile == 'Wing' or nozzle_profile == 'Masuya'\
	or nozzle_profile == 'Sellam':
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

	delta = ((2*q1*Cp/3. + P1 - Pav)**2)*(Rj**2)*(psi**2) + 8*Aj*gammaj*Poj*(np.pi*q1*Cp/4. + 0.5*(P1 - Pav)*np.pi)*(term1**0.5)*((2./(gammaj+1))**(1./(gammaj-1)))
	
	h = (-(2.*q1*Cp/3 + P1 - Pav)*Rj*psi + delta**0.5)/(2.*(np.pi*q1*Cp/4. + 0.5*(P1 - Pav)*np.pi))
	
	if Schilling == 1:
		rise = 0.582*(1 + 0.5*(gamma1 - 1)*M1**2)**(-0.1197*gamma1/(gamma1-1))
	if Zukoski == 1:
		rise = (1 + 0.5*M1)**-1
	
	Pp = P1*1./rise
	
	# oblique shock relations
	beta = np.arcsin((((1./rise)*(gamma1+1) + (gamma1-1))/(2*gamma1*M1**2))**0.5)
	epsilon = np.arctan((M1**2*np.sin(2*beta) - 2/np.tan(beta))/(2 + M1**2*(gamma1 + np.cos(2*beta))))
	ls_new = h*(1./np.sin(epsilon) - 1)

	error = ls_new - ls_old
	ls_old = ls_new
	i += 1
	
	print "M1  | Psep  |  Pp   | xs/xt  | Beta | Eps. |   h   |"
	print '{:.1f}'.format(round(M1,1)),"|", '{:.3f}'.format(round(P1/Po1,3)),"|",'{:.3f}'.format(round(Pp/Po1,3)),"|",'{:.3f}'.format(round(xs/xt,3)),"|", \
	'{:.1f}'.format(round(beta*180/np.pi,1)),"|",'{:.1f}'.format(round(epsilon*180/np.pi,1)),"|",'{:.3f}'.format(round(h,3)),"|"
	print "------------------------------------------------------"


#############################
########## DYNALPY ########## 
#############################
Pc = Po1/((1 + 0.5*(gamma1-1))**(gamma1/(gamma1-1)))
Pj = Poj/((1 + 0.5*(gammaj-1))**(gammaj/(gammaj-1)))
Fc = np.pi*rth**2*(1 + gamma1)*Pc
Fxj = Aj*(1 + gammaj)*Pj*np.sin(alpha*np.pi/180)
Fyj = Aj*(1 + gammaj)*Pj*np.cos(alpha*np.pi/180)


#############################
#######  FORCES #############
#############################
Fxs = 0
Fys = 0
Fxb = 0
Fxd = 0
Fxr = 0
Fyb = 0 
Fyd = 0
Fyr = 0

Nxs = 50
x_up = np.linspace(xs, xm, Nxs)
radius_up = np.zeros(Nxs)


#############################
####### UPSTREAM ############ 
#############################
# populating the radius at each x location
for i in range(len(x_up)):
	radius_up[i] = rth + x_up[i]*np.tan(alpha*np.pi/180)  

# integrating the forces in the x-direction	
for i in range(1, len(x_up)):
	area_ratio = (radius_up[i]/rth)**2
	M = var('M')
	# Solving A/A* for M
	eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio
	M1_up = float(nsolve(eq, M, M1guess))
	# Getting dynamic pressure, assuming constant Po1 before the shock
	P1_up = Po1/((1 + 0.5*(gamma1-1)*M1_up**2)**(gamma1/(gamma1-1)))

	Fxs += (Pp - P1_up)*psi*radius_up[i]*(x_up[i]-x_up[i-1])*np.sin(alpha*np.pi/180)
	Fys += (Pp - P1_up)*psi*radius_up[i]*(x_up[i]-x_up[i-1])*np.cos(alpha*np.pi/180)


#############################
########## EDGES ############ 
#############################
x_edge = x_up
radius_edge = radius_up

# bow shock
Rc = h*1.143*np.exp(0.54/(M1 - 1)**1.2)
Delta = h*0.143*np.exp(3.24/M1**2)
mu = np.arcsin(1/M1)
Rs = ls_new + Rc - Delta
cotan = 1/np.tan(mu)
y = var('y')
Billig_eq = (h + Delta - Rs*cotan**2*((1 + (y)**2*(np.tan(mu))**2/Rs**2)**0.5 - 1))*-1 + xs - xm
y1 = float(nsolve(Billig_eq, y, 0.1))
psi_max = np.arctan(y1/radius_edge[-1]) + psi/2
# psi_max = y1/radius_edge[-1]

print "Psi max is:", np.round(psi_max*180/np.pi)

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
	angle = (x_edge[j] - xs)*(psi_max - psi/2)/(xm - xs) 

	Fxb += 2*Pb*angle*radius_edge[j]*(x_edge[j]-x_edge[j-1])*np.sin(alpha*np.pi/180)
	Fyb += 2*Pb*angle*radius_edge[j]*(x_edge[j]-x_edge[j-1])*np.cos(alpha*np.pi/180)


#############################
####### DOWNSTREAM ########## 
#############################

# OVER-EXPANDED CASE
if NPR < NPRd: 
	x_down = np.linspace(xm, xt, Nxs)
	radius_down = np.zeros(Nxs)
	FSS = False
	Zukoski = 1
	Summerfield = 0
	Schmucker = 0

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

		# Check whether pressure post-shock drops below Pinfty (using Zukoski's Criteria)
		# If True, then shock wave will form to equalize the pressure, break from loop
		# If False, then no shock will form, carry on integration as normal
		if Zukoski == 1:
			ratio = P1_down*(1 + 0.5*M1_down)
		if Summerfield == 1:
			ratio = P1_down*1./0.4
		if Schmucker == 1:
			ratio = P1_down*(1.88*M1_down - 1)**0.64
		# print P1_down, ratio, Pav 
		if ratio < Pav:
			print "Secondary shock appears at", np.round(x_down[k], 3)
			FSS = True
			break

		Fxd += 2*(Pav-P1_down)*psi_max*radius_down[k]*(x_down[k] - x_down[k-1])*np.sin(alpha*np.pi/180)
		Fyd += 2*(Pav-P1_down)*psi_max*radius_down[k]*(x_down[k] - x_down[k-1])*np.cos(alpha*np.pi/180)

	# Still not completely sure about this one
	# The calculation of Fxd, Fyd = 0 assumes that
	# the shock forms immediately downstream of injection
	# and thus the pressures equalize right away.
	if FSS == False:
		print "No secondary shock forms inside the nozzle."
		Fxd = 0
		Fyd = 0


# UNDER-EXPANDED CASE
if NPR >= NPRd:
	# Calculate the Mach post first oblique shock (using Green's criteria)
	M2 = M1*0.76
	# Calculate psi_r (angle of reattachment) using experimental results 
	# conducted by Mangin at the ONERA lab
	psi_r = 32.6 - 29.2/M2

	xr = xm + h/np.tan(psi_r*np.pi/180) + Rc - Delta
	print "Reattachment takes place at", np.round(xr, 3)
	x_down = np.linspace(xm, xr, Nxs)
	radius_down = np.zeros(Nxs)

	# The value of M3 is calculated by solving the Prandtl-Meyer Expansion Fan eq.
	# Assuming that epsilon is the turning angle
	# The value of M4 is calculated similarly, by assuming that psi_r is the turning angle
	fan_angle1 = ((gamma1+1)/(gamma1-1))**0.5*np.arctan(((gamma1-1)/(gamma1+1)*(M2**2 - 1))**0.5) - np.arctan((M2**2 - 1)**0.5)
	fan_angle2 = epsilon + fan_angle1
	fan_angle3 = fan_angle2 + psi_r*np.pi/180
	fan_eq = ((gamma1+1)/(gamma1-1))**0.5*sp.atan(((gamma1-1)/(gamma1+1)*(M**2 - 1))**0.5) - sp.atan((M**2 - 1)**0.5) - fan_angle3
	M4 = float(nsolve(fan_eq, M, M2))

	first_term = (((gamma1+1)*M1**2*(np.sin(beta))**2)/((gamma1-1)*M1**2*(np.sin(beta))**2 + 2.))**(gamma1/(gamma1-1))
	second_term = ((gamma1+1)/(2*gamma1*M1**2*(np.sin(beta))**2 - (gamma1-1)))**(1./(gamma1-1))
	oblique_shock = first_term*second_term
	Po2 = Po1*oblique_shock
	Pr_init = Po2/((1 + 0.5*(gamma1-1)*M4**2)**(gamma1/(gamma1-1)))

	# populating the radius at each x location
	for k in range(len(x_down)):
		radius_down[k] = rth + x_down[k]*np.tan(alpha*np.pi/180)

	# integrating the forces in the x-direction	
	for k in range(1, len(x_down)-1):
	 	area_ratio = (radius_down[k]/rth)**2
		M = var('M')
		# Solving A/A* for M
		eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio
		M1_down = float(nsolve(eq, M, M1guess))
		# Getting dynamic pressure, assuming constant Po1 before the shock
		P1_down = Po1/((1 + 0.5*(gamma1-1)*M1_down**2)**(gamma1/(gamma1-1)))
		
		Pr = (x_down[k] - xm)*(-Pr_init + P1_down)/(xr - xm) + (Pr_init - P1_down)
		Fxr += (Pr_init - P1_down)*psi*radius_down[k]*(x_down[k+1] - x_down[k-1])*np.sin(alpha*np.pi/180)
		Fyr += (Pr_init - P1_down)*psi*radius_down[k]*(x_down[k+1] - x_down[k-1])*np.cos(alpha*np.pi/180)


#############################
###### POST-PROCESSING ###### 
#############################
Fx_total = Fc + Fxj + Fxs + Fxb + Fxd + Fxr
Fy_total = Fyj + Fys + Fyb + Fyd + Fyr
def_angle = np.arctan(Fy_total/Fx_total) 
fm = SPR*(Aj/np.pi/rth**2)

#############################
########## OUTPUT ########### 
#############################
print "Sum of forces in x, Fx (N):", np.round(Fx_total)
print "Sum of forces in y, Fy (N):", np.round(Fy_total)
print "Net deflection angle (deg.):", np.round(180*def_angle/np.pi, 1)
print "Mass-flow ratio (%):", np.round(fm*100, 2)
print ""
print "SUMMARY OF FORCES IN X"
print "------------------------------------------------------"
print "Fc     | Fxj | Fxs  | Fxb | Fxd |  Fxr |"
print '{:.1f}'.format(round(Fc,1)),"|", '{:.1f}'.format(round(Fxj,1)),"|",'{:.1f}'.format(round(Fxs,3)),"|", \
'{:.1f}'.format(round(Fxb,1)),"|", '{:.1f}'.format(round(Fxd,1)),"|",'{:.1f}'.format(round(Fxr,1)),"|"
print "------------------------------------------------------"

print ""
print "SUMMARY OF FORCES IN Y"
print "------------------------------------------------------"
print "Fyj  |  Fys  | Fyb  | Fyd |  Fyr  |"
print '{:.1f}'.format(round(Fyj,1)),"|", '{:.1f}'.format(round(Fys,1)),"|",'{:.1f}'.format(round(Fyb,3)),"|", \
'{:.1f}'.format(round(Fyd,1)),"|", '{:.1f}'.format(round(Fyr,1)),"|"
print "------------------------------------------------------"