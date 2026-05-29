from sympy import nsolve, var
import numpy as np
import matplotlib.pyplot as plt

#############################
### INITIAL CONDITIONS ######
#############################
Masuya = 1
Comparison = 0

if Masuya == 1:
	gamma1 = 1.4
	gammaj = 1.4
	psi = 60.*np.pi/180
	NPRd = 42.1
	NPR = 42
	SPR = 1.0
	xt = 0.09
	xm = 0.05
	rth = 0.013
	exp_ratio = 4.6
	Me = 3.1
	Pinfty = 100*10**3
	Po1 = NPRd*Pinfty
	Poj = SPR*Po1
	alpha = 9.6
	Pav = Po1/NPR
	Pj = Pav
	# yj and Aj are not explicity given 
	# I calculated them based on a digitized version 
	# of Fig. 5-39 in Maarouf's Thesis. 
	# Everything else is explicitly mentioned.
	yj = 0.0086
	Aj = 12.6*(10**-6)

if Comparison == 1:
	gamma1 = 1.4
	gammaj = 1.4
	psi = 30.*np.pi/180
	NPRd = 42.1
	NPR = 20
	SPR = 1.0
	xt = 0.09
	xm = 0.047
	rth = 0.013
	exp_ratio = 4.6
	Me = 3.1
	Pinfty = 100*10**3
	Po1 = NPRd*Pinfty
	Poj = SPR*Po1
	alpha = 9.6
	Pav = Po1/NPR
	Pj = Pav
	# yj and Aj are not explicity given 
	# I calculated them based on a digitized version 
	# of Fig. 5-51 in Maarouf's Thesis. 
	# Everything else is explicitly mentioned.
	yj = 0.0086*2
	Aj = 12.6*(10**-6)

error = 1
tol = 10**-6

M1guess = 2.3
ls_old = 0
Zukoski = 1
Schilling = 0
i = 0

#############################
###### SEP DISTANCE #########
#############################
while np.abs(error) > tol:
	xs = xm - ls_old
	length_ratio = xs/xt
	# rs = rth + rth*(exp_ratio**0.5 - 1)*length_ratio
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

	delta = ((2*q1*Cp/3. + P1 - Pav)**2)*(yj**2)*(psi**2) + 8*Aj*gammaj*Poj*(np.pi*q1*Cp/4. + 0.5*(P1 - Pav)*np.pi)*(term1**0.5)*((2./(gammaj+1))**(1./(gammaj-1)))
	
	h = (-(2*q1*Cp/3 + P1 - Pav)*yj*psi + delta**0.5)/(2.*(np.pi*q1*Cp/4. + 0.5*(P1 - Pav)*np.pi))
	
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
	
	# xp is not very important at the moment; it will be interpolated later
	xp = xs + h

	print "M1  | Psep  |  Pp   | xs/xt | xp/xt | Beta | Eps. |   h   |"
	print '{:.1f}'.format(round(M1,1)),"|", '{:.3f}'.format(round(P1/Po1,3)),"|",'{:.3f}'.format(round(Pp/Po1,3)),"|",'{:.3f}'.format(round(xs,3)),"|", \
	'{:.3f}'.format(round(xp/xt,3)),"|", '{:.1f}'.format(round(beta*180/np.pi,1)),"|",'{:.1f}'.format(round(epsilon*180/np.pi,1)),"|",'{:.3f}'.format(round(h,3)),"|"
	print "------------------------------------------------------"

#############################
#######  FORCES #############
#############################
Fxs = 0
Fys = 0
Fxb = 0
Fxd = 0
Fyb = 0 
Fyd = 0

Nxs = 20
x_up = np.linspace(xs, xm, Nxs)
radius_up = np.zeros(Nxs)
x_edge = np.zeros(Nxs)
radius_edge = np.zeros(Nxs)
x_down = np.linspace(xm, xt, Nxs)
radius_down = np.zeros(Nxs)

#############################
# Forces upstream of sep. 
#############################

# populating the radius at each x location
for i in range(len(x_up)):
	radius_up[i] = rth + x_up[i]*np.tan(alpha*np.pi/180)  

# integrating the forces in the x-direction	
for i in range(1, len(x_up)-1):
	area_ratio = (radius_up[i]/rth)**2
	M = var('M')
	# Solving A/A* for M
	eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio
	M1_up = float(nsolve(eq, M, M1guess))
	# Getting dynamic pressure, assuming constant Po1 before the shock
	P1_up = Po1/((1 + 0.5*(gamma1-1)*M1_up**2)**(gamma1/(gamma1-1)))
	 
	Fxs += (Pp - P1_up)*psi*radius_up[i]*0.5*(radius_up[i+1]-radius_up[i-1])
	Fys += 2*(Pp - P1_up)*radius_up[i]*0.5*(x_up[i+1] - x_up[i-1])*np.sin(psi)

#############################
# Forces at edges of sep. 
#############################

# bow shock
Rc = h*1.143*np.exp(0.54/(M1 - 1)**1.2)
Delta = h*0.143*np.exp(3.24/M1**2)
mu = np.arcsin(1/M1)

Rs = ls_new + Rc - Delta
print Rs
cotan = 1/np.tan(mu)
y_bow = 0
y = var('y')
x_bow = h + Delta - Rs*cotan**2*((1 + y_bow**2*(np.tan(mu))**2/Rs**2)**0.5 - 1)
eq = h + Delta - Rs*cotan**2*((1 + y**2*(np.tan(mu))**2/Rs**2)**0.5 - 1) - 0
y1 = float(nsolve(eq, y, 0))
print y1/radius_up[-1], psi

# populating the radius at each x location
for j in range(len(x_up)):
	x_edge[j] = x_up[j] + x_bow
	radius_edge[j] = rth + x_edge[j]*np.tan(alpha*np.pi/180)

# integrating the forces in the x-direction	
for j in range(1, len(x_edge)-1):
 	area_ratio = (radius_edge[j]/rth)**2
	M = var('M')
	# Solving A/A* for M
	eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio
	M1_edge = float(nsolve(eq, M, M1guess))
	# Getting dynamic pressure, assuming constant Po1 before the shock
	P1_edge = Po1/((1 + 0.5*(gamma1-1)*M1_edge**2)**(gamma1/(gamma1-1)))

	Fxb += 2*(Pp - P1_edge)*radius_edge[j]*0.5*(radius_edge[j+1] - radius_edge[j-1])
	Fyb += 2*(Pp - P1_edge)*radius_edge[j]*0.5*(x_edge[j+1] - x_edge[j-1])*np.sin(psi)

#############################
# Forces downstream of sep. 
#############################

# populating the radius at each x location
for k in range(len(x_up)):
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

	Fxd += 2*(Pav - P1_down)*radius_down[k]*0.5*(radius_down[k+1] - radius_down[k-1])
	Fyd += 2*(Pav - P1_down)*radius_down[k]*0.5*(x_down[k+1] - x_down[k-1])*np.sin(psi)

# Dynalpy Flux (due to flow onto throat)
Pc = Po1/((1 + 0.5*(gamma1-1))**(gamma1/(gamma1-1)))
Pj = Poj/((1 + 0.5*(gammaj-1))**(gammaj/(gammaj-1)))
Fc = np.pi*rth**2*(1 + gamma1)*Pc
Fxj = Aj*(1 + gammaj)*Pj*np.sin(alpha*np.pi/180)
Fyj = -Aj*(1 + gammaj)*Pj*np.cos(alpha*np.pi/180)


Fx_total = Fc + Fxs + Fxj + Fxb + Fxd 
Fy_total = Fys + Fyj + Fyb + Fyd

print Fx_total, Fy_total, 180*np.arctan(Fy_total/Fx_total)/np.pi

y_axis1 = [Fc, Fxs, Fxj, Fxb, Fxd]
y_axis2 = [Fyj, Fys, Fyb, Fyd]

print y_axis1, y_axis2
x_axis1 = range(len(y_axis1))
x_axis2 = range(len(y_axis2))

plt.figure()
plt.bar(x_axis1, y_axis1)
plt.title('Forces in x')

plt.figure()
plt.bar(x_axis2, y_axis2)
plt.title('Forces in y')
# plt.show()