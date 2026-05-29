from sympy import *
import numpy as np

#############################
###### GLOBAL NUMERICS ######
#############################
error = 1
tol = 10**-3
M1guess = 2.3
ls_old = 0
nozzle_profile = "NASA2D"
gamma1 = 1.4
gammaj = 1.4

#############################
### INITIAL CONDITIONS ######
#############################
if nozzle_profile == 'NASA2D':
	NPRd = 8.78
	Pinfty = 100000.0
	Po1 = NPRd*Pinfty
	NPR = 4.6
	SPR = 0.7

if nozzle_profile == 'Li':
	# NPRd is calculated from Mach 2 exit,
	# from Li 2016.
	NPRd = 7.82
	Pinfty = 100000.0
	Po1 = NPRd*Pinfty
	NPR = 5.0
	SPR = 0.5


#############################
######## GEOMETRY ###########
#############################
if nozzle_profile == 'NASA2D':
	# The length of the entire nozzle is 4.55 inches
	# When converted to m and subtracted from conv
	# section; end up with 0.0577.
	xt = 0.0577
	alpha = 11.01
	# This was obtained from the area of the throat
	# in Waithe 2003.
	rth = 0.02977

if nozzle_profile == 'Li':
	# xt is calculated from fig. 4.
	xt = 90.0*10**-3
	# alpha is computed using xt and 
	# area ratio (height of exit/throat).
	alpha = 3.814
	# height of throat given.
	rth = 20.0*10**-3


#############################
####### INJECTION ###########
#############################

# CONFIGURATION 1
if nozzle_profile == 'NASA2D':
	Cd = 0.95
	# The injection slot starts at 4.1 inches.
	# After subtracting the converging section, 
	# and dividing by xt; end up with 0.8.
	xm = 0.8*xt
	Poj = SPR*Po1
	# Slot width is 0.08 inches.
	b = 0.08*0.0254
	w = 3.49*0.0254
	Aj = Cd*b*w

if nozzle_profile == 'Li':
	Cd = 1.0
	# The injection slot is 5 mm from exit;
	# this yields 0.94 of total length.
	xm = 0.94*xt
	Poj = SPR*Po1
	# Slot width is 1 mm.
	b = 1.0*10**-3
	# Nozzle width given as 40 mm.
	w = 40.0*10**-3
	Aj = Cd*b*w


#############################
### OPERATING CONDITIONS ####
#############################
if NPR >= NPRd:
	Pav = Po1/NPRd
	print "Ideal or under-expanded regime."
if NPR < NPRd:
	Pav = Po1/NPR
	print "Over-expanded regime." 

Pj = Pav


#############################
######## SEP. CRITERIA ######
#############################
Zukoski = 1
Schilling = 0
i = 0


#############################
###### SEP. DISTANCE ########
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
	if nozzle_profile == 'NASA2D' or nozzle_profile == 'Li':
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
	CD = 2.*Cp/3.

	q1 = 0.5*gamma1*P1*M1**2
	term1 = (1./(gammaj**2 - 1))*(1 - (Pj/Poj)**((gammaj-1)/gammaj))

	h = 2.*(Cd/CD)*b*Poj*gammaj*term1**0.5*(2./(gammaj+1))**(1./(gammaj-1))/((P1 - Pav) + 1./3.*P1*gamma1*M1**2*Cp) 

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
if nozzle_profile == 'NASA2D':
	print '{:.1f}'.format(round(M1,1)),"|", '{:.3f}'.format(round(P1/Po1,3)),"|",'{:.3f}'.format(round(Pp/Po1,3)),"|",'{:.3f}'.format(round(xs/xt + 1,3)),"|", \
	'{:.1f}'.format(round(beta*180/np.pi,1)),"|",'{:.1f}'.format(round(epsilon*180/np.pi,1)),"|",'{:.3f}'.format(round(h,3)),"|"
if nozzle_profile == 'Li':
	print '{:.1f}'.format(round(M1,1)),"|", '{:.3f}'.format(round(P1/Po1,3)),"|",'{:.3f}'.format(round(Pp/Po1,3)),"|",'{:.3f}'.format(round(xs/xt,3)),"|", \
	'{:.1f}'.format(round(beta*180/np.pi,1)),"|",'{:.1f}'.format(round(epsilon*180/np.pi,1)),"|",'{:.3f}'.format(round(h,3)),"|"
print "------------------------------------------------------"


#############################
########## DYNALPY ##########
#############################
Pc = Po1/((1 + 0.5*(gamma1-1))**(gamma1/(gamma1-1)))
Pj = Poj/((1 + 0.5*(gammaj-1))**(gammaj/(gammaj-1)))
Fc = np.pi*rth**2*(1 + gamma1)*Pc
Fxj = Cd*b*w*(1 + gammaj)*Pj*np.sin(alpha*np.pi/180)
Fyj = Cd*b*w*(1 + gammaj)*Pj*np.cos(alpha*np.pi/180)


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
	Fxs += (Pp - P1_up)*w*(x_up[i]-x_up[i-1])*np.sin(alpha*np.pi/180)
	Fys += (Pp - P1_up)*w*(x_up[i]-x_up[i-1])*np.cos(alpha*np.pi/180)


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
			print "Secondary shock appears at", np.round(x_down[k]/xt + 1, 3)
			FSS = True
			break

		Fxd += (Pav-P1_down)*w*(x_down[k] - x_down[k-1])*np.sin(alpha*np.pi/180)
		Fyd += (Pav-P1_down)*w*(x_down[k] - x_down[k-1])*np.cos(alpha*np.pi/180)

	if FSS == False:
		print "No secondary shock forms inside the nozzle."
		Fxd = 0
		Fyd = 0

# UNDER-EXPANDED CASE
if NPR >= NPRd:
	print "Ignoring reattachment."


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
print "Fc    | Fxj | Fxs | Fxb | Fxd | Fxr |"
print '{:.1f}'.format(round(Fc,1)),"|", '{:.1f}'.format(round(Fxj,1)),"|",'{:.1f}'.format(round(Fxs,3)),"|", \
'{:.1f}'.format(round(Fxb,1)),"|", '{:.1f}'.format(round(Fxd,1)),"|",'{:.1f}'.format(round(Fxr,1)),"|"
print "------------------------------------------------------"

print ""
print "SUMMARY OF FORCES IN Y"
print "------------------------------------------------------"
print "Fyj  | Fys | Fyb | Fyd | Fyr |"
print '{:.1f}'.format(round(Fyj,1)),"|", '{:.1f}'.format(round(Fys,1)),"|",'{:.1f}'.format(round(Fyb,3)),"|", \
'{:.1f}'.format(round(Fyd,1)),"|", '{:.1f}'.format(round(Fyr,1)),"|"
print "------------------------------------------------------"
