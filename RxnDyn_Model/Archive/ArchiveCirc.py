if Conical == 1:
	gamma1 = 1.4
	gammaj = 1.4
	R1 = 287.058
	Rj = 287.058
	Dj = 6.*(10**-3)
	Aj = (np.pi/4.)*Dj**2
	NPR = 37.5
	SPR = 1.
	xt = 100*(10**-3)
	xm = 0.9*xt 
	rth = 9.72*(10**-3)
	exp_ratio = 4.234
	Me = 3.0
	Mj_exit = ((((SPR*NPR)**((gammaj-1)/gammaj)) - 1)*2./(gammaj-1))**0.5
	Po1 = 300.*(10**3)
	Poj = SPR*Po1
	alpha = 5.42
	# Cd is the only parameter that is
	# not explicitly mentioned. 
	Cd = 0.12
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
	
if Masuya == 1:
	gamma1 = 1.4
	gammaj = 1.4
	R1 = 287.058
	Rj = 287.058
	Dj = 6.*(10**-3)
	Cd = 0.12
	Aj = 0.25*np.pi*Dj**2
	NPR = 20.
	SPR = 1.
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
	# Tinfty = 298.
	# To1 = (1 + 0.5*(gamma1-1)*Me**2)*Tinfty
	# Toj = (1 + 0.5*(gammaj-1)*Mj_exit**2)*Tinfty
	To1 = 616.
	Toj = 616.
	Pe = Po1/((1 + 0.5*(gamma1-1)*Me**2)**(gamma1/(gamma1-1)))
	NPR_critical = Po1/Pe

	print ""
	if NPR >= NPR_critical:
		Pav = Pe
		print "Ideal or under-expanded regime."
	if NPR < NPR_critical:
		Pav = Pinfty
		print "Over-expanded regime."

if TIC == 1:
	gamma1 = 1.4
	gammaj = 1.4
	R1 = 287.058
	Rj = 287.058
	Dj = 5.8*(10**-3)
	Aj = (np.pi/4.)*Dj**2
	NPR = 37.5
	SPR = 1.0
	xt = 68.*(10**-3)
	xm = 0.88*xt 
	rth = 10.*(10**-3)
	exp_ratio = 4.87
	Me = 3.03
	Mj_exit = ((((SPR*NPR)**((gammaj-1)/gammaj)) - 1)*2./(gammaj-1))**0.5
	Po1 = 300.*(10**3)
	Poj = SPR*Po1
	alpha = 10.06
	# Cd is the only parameter that is
	# not explicitly mentioned. 
	Cd = 0.12
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


#############################
######### PLOTTING ########## 
#############################
y_axis1 = [Fc, Fxj, Fxs, Fxb, Fxd, Fxr]
y_axis2 = [Fyj, Fys, Fyb, Fyd, Fyr]
x_axis1 = range(len(y_axis1))
x_axis2 = range(len(y_axis2))

NPR = [10, 15, 20, 30, 40]
Maarouf = [6.96, 5.5, 4.2, 3.94, 3.85]
CFD = [9.24, 6.6, 4.6, 4.0, 3.9]
Younes = [7.4, 5.5, 4.5, 3.6, 3.1]

# plt.figure()
# plt.plot(NPR, CFD, linestyle='-', linewidth=1.5, marker='x', label='Maarouf CFD', color='k')
# plt.plot(NPR, Maarouf, linestyle='--', linewidth=1.25, marker='o', label='Maarouf Model', color='g')
# plt.plot(NPR, Younes, linestyle=':', linewidth=1.5, marker='^', label='Younes Model', color='b')
# plt.ylabel('$\\delta$', fontname='Times New Roman', fontsize=16)
# plt.xlabel('NPR', fontname='Times New Roman', fontsize=16)
# plt.xlim(10, 40)
# plt.ylim(3, 12)
# plt.xticks([10, 20, 30, 40], ["10", "20", "30", "40"], fontname='Times New Roman', fontsize=16)
# plt.yticks([3, 6, 9, 12], ["3", "6", "9", "12"], fontname='Times New Roman', fontsize=16)
# plt.tick_params('both', which='major', direction='in', right=True, top=True, pad=10)
# plt.tight_layout()
# plt.legend(loc='upper right',prop={'size':14,"family":"Times New Roman"})
# plt.show()

# plt.figure()
# plt.bar(x_axis1, y_axis1)
# plt.title('Forces in x')

# plt.figure()
# plt.bar(x_axis2, y_axis2)
# plt.title('Forces in y')
# plt.show()
