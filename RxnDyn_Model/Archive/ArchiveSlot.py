if Masuya == 1:
	gamma1 = 1.4
	gammaj = 1.4
	R1 = 287.058
	Rj = 287.058
	psi = 60.*np.pi/180
	NPR = 20.
	SPR = 1.0
	alpha = 9.6
	xt = 0.09
	xm = 0.05
	rth = 0.013
	exp_ratio = 4.6
	Me = 3.09
	Pinfty = 100.*10**3
	Po1 = NPR*Pinfty
	Poj = SPR*Po1
	Tinfty = 298.
	To1 = Tinfty*(Po1/Pinfty)**((gamma1-1)/gamma1)
	Toj = Tinfty*(Poj/Pinfty)**((gammaj-1)/gammaj)
	Pe = Po1/((1 + 0.5*(gamma1-1)*Me**2)**(gamma1/(gamma1-1)))
	NPR_critical = Po1/Pe

	print ""
	if NPR >= NPR_critical:
		Pav = Pe
		print "Ideal or under-expanded regime."
	if NPR < NPR_critical:
		Pav = Pinfty
		print "Over-expanded regime."

	# yj and Aj are not explicity given 
	# I calculated them based on a digitized version 
	# of Fig. 5-39 in Maarouf's Thesis. 
	# Everything else is explicitly mentioned.
	yj = 0.0086
	Aj = 12.6*(10**-6)
	b = Aj/yj/psi


# The Wing test case is not currently
# working because the injection is close
# to the throat that it yields a larger
# separation distance than the length of
# the nozzle.
if Wing == 1:
	gamma1 = 1.4
	gammaj = 1.4
	R1 = 287.058
	Rj = 287.058
	psi = 60.*np.pi/180
	NPR = 3.
	SPR = 1.
	xt = 5./100
	xm = 3./100 
	rth = 2.77/100
	exp_ratio = 1.74
	alpha = 10.0
	Me = 2.0
	Mj_exit = ((((SPR*NPR)**((gammaj-1)/gammaj)) - 1)*2./(gammaj-1))**0.5
	Po1 = 300.*(10**3)
	Poj = SPR*Po1
	Pinfty = Po1/NPR
	To1 = 300.
	Toj = 300.
	Pe = Po1/((1 + 0.5*(gamma1-1)*Me**2)**(gamma1/(gamma1-1)))
	NPR_critical = Po1/Pe

	print ""
	if NPR >= NPR_critical:
		Pav = Pe
		print "Ideal or under-expanded regime."
	if NPR < NPR_critical:
		Pav = Pinfty
		print "Over-expanded regime."

	yj = 2.*rth
	Aj = 7.*25.0*(10**-6)
	b = Aj/yj/psi

if Comparison == 1:
	gamma1 = 1.4
	gammaj = 1.4
	R1 = 287.058
	Rj = 287.058
	psi = 30.*np.pi/180
	NPR = 20.
	SPR = 1.0
	xt = 0.09
	xm = 0.047
	rth = 0.013
	exp_ratio = 4.6
	Me = 3.09
	Mj_exit = ((((SPR*NPR)**((gammaj-1)/gammaj)) - 1)*2./(gammaj-1))**0.5
	Po1 = 2.*10**6
	Poj = SPR*Po1
	alpha = 9.6
	Pinfty = Po1/NPR
	Tinfty = 298.
	To1 = Tinfty*(Po1/Pinfty)**((gamma1-1)/gamma1)
	Toj = Tinfty*(Poj/Pinfty)**((gammaj-1)/gammaj)
	Pe = Po1/((1 + 0.5*(gamma1-1)*Me**2)**(gamma1/(gamma1-1)))
	NPR_critical = Po1/Pe

	print ""
	if NPR >= NPR_critical:
		Pav = Pe
		print "Ideal or under-expanded regime."
	if NPR < NPR_critical:
		Pav = Pinfty
		print "Over-expanded regime."

	# yj and Aj are not explicity given 
	# I calculated them based on a digitized version 
	# of Fig. 5-39 in Maarouf's Thesis. 
	# Everything else is explicitly mentioned.
	yj = 0.0086*2
	# Aj = 20.0*(10**-6)
	Aj = 12.6*(10**-6)
	b = Aj/yj/psi

if Conical == 1:
	gamma1 = 1.4
	gammaj = 1.4
	R1 = 287.058
	Rj = 287.058
	psi = 30.*np.pi/180
	NPR = 37.5
	SPR = 1.
	xt = 100.*(10**-3)
	xm = 0.9*xt 
	rth = 9.72*(10**-3)
	exp_ratio = 4.234
	Me = 3.0
	Mj_exit = ((((SPR*NPR)**((gammaj-1)/gammaj)) - 1)*2./(gammaj-1))**0.5
	Po1 = 300.*(10**3)
	Poj = SPR*Po1
	alpha = 5.42
	Pinfty = Po1/NPR
	To1 = 243.
	Toj = 260.
	Pe = Po1/((1 + 0.5*(gamma1-1)*Me**2)**(gamma1/(gamma1-1)))
	NPR_critical = Po1/Pe

	print ""
	if NPR >= NPR_critical:
		Pav = Pe
		print "Ideal or under-expanded regime."
	if NPR < NPR_critical:
		Pav = Pinfty
		print "Over-expanded regime."

	yj = 2.0*rth
	Aj = 25.1*(10**-6)
	b = Aj/yj/psi


Pj = Pav
error = 1
tol = 10**-6

M1guess = 2.3
ls_old = 0
Zukoski = 1
Schilling = 0
i = 0

#############################
######### PLOTTING ########## 
#############################
y_axis1 = [Fc, Fxj, Fxs, Fxb, Fxd, Fxr]
y_axis2 = [Fyj, Fys, Fyb, Fyd, Fyr]
x_axis1 = range(len(y_axis1))
x_axis2 = range(len(y_axis2))

NPR_Maarouf = [10, 12, 15, 20, 23, 30, 42.5]
Maarouf = [11.0, 9.5, 8.01, 5.5, 5.42, 5.31, 5.22]
NPR_CFD = [11, 15, 20, 30, 42.5]
CFD = [12.6, 8.9, 6.01, 5.44, 5.22]
NPR_Younes = [10, 20, 23, 30, 42.5]
Younes1 = [9.9, 11.0, 11.6, 6.7, 6.7, 6.6, 5.5]
Younes2 = [9.9, 9.9, 8.5, 7.0, 6.5, 5.7, 5.5]

SPR = [0.3, 0.5, 0.7, 1.0, 1.2]
# NPR = 42.5
SPR_Angle = [1.0, 2.1, 3.4, 5.5, 7.0]

# NPR = 11
SPR_Angle2 = [3.3, 5.5, 7.6, 10.6, 12.6]

# NPR = 20
SPR_Angle3 = [1.7, 3.1, 4.5, 6.7, 8.2]

# plt.figure()
# plt.plot(NPR_CFD, CFD, linestyle='-', linewidth=1.5, marker='x', label='Maarouf CFD', color='k')
# plt.plot(NPR_Maarouf, Maarouf, linestyle='--', linewidth=1.25, marker='o', label='Maarouf Model', color='g')
# plt.plot(NPR_Maarouf, Younes1, linestyle=':', linewidth=1.5, marker='^', label='Younes Model', color='b')
# plt.plot(NPR_Maarouf, Younes2, linestyle='-.', linewidth=1.75, marker='v', label='Younes Model 2', color='r')
# plt.ylabel('$\\delta$', fontname='Times New Roman', fontsize=16)
# plt.xlabel('NPR', fontname='Times New Roman', fontsize=16)
# plt.xlim(10, 42.5)
# plt.xticks()
# plt.ylim(3, 15)
# plt.xticks([10, 20, 30, 40], ["10", "20", "30", "40"], fontname='Times New Roman', fontsize=16)
# plt.yticks([3, 6, 9, 12, 15], ["3", "6", "9", "12", "15"], fontname='Times New Roman', fontsize=16)
# plt.tick_params('both', which='major', direction='in', right=True, top=True, pad=10)
# plt.tight_layout()
# plt.legend(loc='upper right',prop={'size':14,"family":"Times New Roman"})
# plt.show()

# plt.figure()
# plt.plot(SPR, SPR_Angle, linestyle='-', linewidth=1.5, marker='s', label='Younes Model', color='k')
# plt.ylabel('$\\delta$', fontname='Times New Roman', fontsize=16)
# plt.xlabel('SPR', fontname='Times New Roman', fontsize=16)
# plt.xlim(0.3, 1.2)
# plt.ylim(1, 7)
# plt.xticks(fontname='Times New Roman', fontsize=16)
# plt.yticks([1, 3, 5, 7], ["1", "3", "5", "7"], fontname='Times New Roman', fontsize=16)
# plt.tick_params('both', which='major', direction='in', right=True, top=True, pad=10)
# plt.tight_layout()
# plt.legend(loc='lower right',prop={'size':14,"family":"Times New Roman"})
# plt.show()

# plt.figure()
# plt.bar(x_axis1, y_axis1)
# plt.title('Forces in x')

# plt.figure()
# plt.bar(x_axis2, y_axis2)
# plt.title('Forces in y')
# plt.show()