from sympy import nsolve, var
import numpy as np
import matplotlib.pyplot as plt

#############################
### INITIAL CONDITIONS ######
#############################
gamma1 = 1.4
Pinfty = 8*(10**3)
NPR = 37.5
xt = 100*(10**-3)
xm = 0.9*xt
rth = 9.72*(10**-3)
exp_ratio = 0.236
Po1 = NPR*Pinfty
plot = 1

N = 51
x_global = np.linspace(0, xt, N)
P_global = np.zeros(N)
M_global = np.zeros(N)
M_global[0] = 1.0
P_global[0] = Po1*(2./(gamma1+1))**(gamma1/(gamma1-1))

for i in range(1, len(x_global)):
	print i
	length_ratio = x_global[i]/xt
	r = rth + rth*(1./(exp_ratio**0.5) - 1)*length_ratio
	area_ratio = (r/rth)**2
	M = var('M')
	eq = (1/M)*(((2 + (gamma1-1)*M**2)/(gamma1+1))**(0.5*(gamma1+1)/(gamma1-1))) - area_ratio

	Mguess = M_global[i-1]
	M_global[i] = float(nsolve(eq, M, Mguess))
	P_global[i] = Po1/((1 + 0.5*(gamma1-1)*M_global[i]**2)**(gamma1/(gamma1-1)))

# function to calculate the cumulative trapezoidal rule
def CumTrapz(y, x, *args):
	integral = np.zeros(len(y))
	for i in range(0, len(y)-1):
		integral[i+1] = integral[i] + np.trapz([y[i], y[i+1]], [x[i], x[i+1]])

	return integral 

import AnalyticalModel
xs = AnalyticalModel.xs
Dj = AnalyticalModel.Dj
sep_index = np.where(x_global < xs)[-1][-1]
end_index = np.where(x_global > xm + Dj/2)[0][0]

Fx1 = 2*np.pi*np.tan(5.42*np.pi/180)*CumTrapz((rth + x_global[0:sep_index]*np.tan(5.42*np.pi/180))*P_global[0:sep_index], x_global[0:sep_index])
Fx4 = 2*np.pi*np.tan(5.42*np.pi/180)*CumTrapz((rth + x_global[end_index:-1]*np.tan(5.42*np.pi/180))*P_global[end_index:-1], x_global[end_index:-1])
print Fx1[-1]
print Fx4[-1]

print M_global[end_index]

if plot == 1:	
	plt.figure(1, figsize=(5, 4))
	plt.plot(x_global/xt, P_global/Po1)
	# plt.xlim(0.5, x_global[-/xt)
	plt.ylabel('$\\frac{P}{Po}$', fontname='Times New Roman', fontsize=20)
	plt.xlabel('$\\chi$', fontname='Times New Roman', fontsize=20)
	plt.xticks(fontname='Times New Roman', fontsize=20)
	plt.yticks(fontname='Times New Roman', fontsize=20)
	plt.tick_params('both', length=5, width=1.5, which='major')
	plt.tight_layout()
	plt.show()