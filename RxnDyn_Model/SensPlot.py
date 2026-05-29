import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 14.0
mpl.rcParams['font.family'] = 'Times New Roman'
# mpl.rcParams['hatch.linewidth'] = 1.75
# mpl.rcParams['hatch.color'] = 'k'

error = 0
bar = 0
contour = 0


labels = ('$r_{th}$', '$r_e$', '$D_j$', 
    '$P_{o1}$', '$P_{oj}$', 
    '$\\gamma_1$', '$\\gamma_j$', 
    '$x_m$', '$x_t$',
    '$C_d$', 
    '$x_s$', '$x_d$',
    '$P_p$', '$\\psi$')
y_pos = np.arange(len(labels))

if error == 1:
	import glob
	filelist = glob.glob("SensitivityResults/NPRd/Convergence/*.txt")
	filelist.sort(key=lambda f: int(filter(str.isdigit, f)))
	total_error = np.zeros(len(filelist))
	Ns = [1, 2, 4, 8, 16, 32, 64, 128, 256]

	for i in range(1, len(filelist)):
		base = np.loadtxt(filelist[i-1])
		new = np.loadtxt(filelist[i])
		for j in range(len(base)):
			total_error[i] += np.abs(new[j] - base[j])

	total_error[0] = 1.*len(base)

	plt.figure(figsize=(5, 4))
	plt.plot(Ns, total_error/len(base), linewidth=1.0, linestyle='-', marker=' ', mfc='none', color='k')

	filelist = glob.glob("SensitivityResults/NPR15/Convergence/*.txt")
	filelist.sort(key=lambda f: int(filter(str.isdigit, f)))
	total_error = np.zeros(len(filelist))
	Ns = [1, 2, 4, 8, 16, 32, 64, 128, 256]

	for i in range(1, len(filelist)):
		base = np.loadtxt(filelist[i-1])
		new = np.loadtxt(filelist[i])
		for j in range(len(base)):
			total_error[i] += np.abs(new[j] - base[j])

	total_error[0] = 1.*len(base)

	plt.plot(Ns, total_error/len(base), linewidth=1.0,linestyle='--', marker=' ', mfc='none', color='k', dashes=(15,10))
	# plt.xlim(0, 13)
	# plt.ylim(-6, 0)
	plt.yscale('log')
	plt.ylim(10**-4, 1.5)
	plt.ylabel('$\\epsilon$')
	plt.xlabel('$N_s$')
	plt.tight_layout()
	# plt.yticks([1, 10^-1, 10^-2, 10^-3], ['1', '$10^{-1}$', '$10^{-2}$', '$10^{-3}$'])
	# plt.tick_params('both', length=5, width=1.5, which='major')
	plt.savefig('SensitivityResults/' + 'Error.pdf', bbox_inches='tight')

if bar == 1:
	delta_S1_NPRd = np.loadtxt('SensitivityResults/NPRd/delta_S1.txt')
	delta_ST_NPRd = np.loadtxt('SensitivityResults/NPRd/delta_ST.txt')
	# print y_pos

	width = 0.5
	plt.figure(figsize=(5, 4))
	# plt.barh(y_pos, delta_S1, height=0.5, alpha=0.3, hatch='...', linewidth=0.5, edgecolor='k')
	plt.bar(y_pos, delta_S1_NPRd, width=width, color='grey', linewidth=1.0, edgecolor='k')
	plt.bar(y_pos+width, delta_ST_NPRd, width=width, color=(0.1, 0.1, 0.1, 0.1), linewidth=1.25, edgecolor='k')
	plt.ylim(0, 0.6)
	plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6], ['0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6'])
	plt.ylabel('$\\delta$ Sensitivity Indices')
	plt.xlim(-0.5, 14)
	plt.xticks(y_pos+0.5*width, labels)
	# plt.gca().invert_yaxis() 
	# plt.tick_params(axis='x', length=0, width=0, which='major')
	plt.tight_layout()
	plt.savefig('SensitivityResults/' + 'delta_NPRd.pdf', bbox_inches='tight')

	delta_S1_NPR15 = np.loadtxt('SensitivityResults/NPR15/delta_S1.txt')
	delta_ST_NPR15 = np.loadtxt('SensitivityResults/NPR15/delta_ST.txt')
	plt.figure(figsize=(5, 4))
	# plt.barh(y_pos, delta_S1, height=0.5, alpha=0.3, hatch='...', linewidth=0.5, edgecolor='k')
	plt.bar(y_pos, delta_S1_NPR15, width=width, color='grey', linewidth=1.0, edgecolor='k')
	plt.bar(y_pos+width, delta_ST_NPR15, width=width, color=(0.1, 0.1, 0.1, 0.1), linewidth=1.25, edgecolor='k')
	plt.ylim(0, 0.8)
	plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], ['0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8'])
	plt.ylabel('$\\delta$ Sensitivity Indices')
	plt.xlim(-0.5, 14)
	plt.xticks(y_pos+0.5*width, labels)
	# plt.gca().invert_yaxis() 
	# plt.tick_params(axis='x', length=0, width=0, which='major')
	plt.tight_layout()
	plt.savefig('SensitivityResults/' + 'delta_NPR15.pdf', bbox_inches='tight')

if contour == 1:
	delta_S2_NPRd = np.loadtxt('SensitivityResults/NPRd/delta_S2.txt')
	plt.figure(figsize=(5,4))
	plt.contourf(delta_S2_NPRd, cmap='Greys')
	# plt.xticks([0.2, 1., 1.9, 3., 4., 5., 6, 7, 8], labels)
	# plt.xlim(0, 13)
	# plt.ylim(-0.5, 14)
	plt.xticks(y_pos, labels)
	plt.yticks(y_pos, labels)
	# plt.xlabel('$\\delta$ Heatmap')
	plt.gca().invert_yaxis()
	# plt.tick_params(axis='both', length=0, width=0, which='major')
	plt.colorbar()
	# plt.clim(-0.016, 0.016)
	plt.tight_layout()
	plt.savefig('SensitivityResults/' + 'heat_NPRd.pdf', bbox_inches='tight')

	delta_S2_NPR15 = np.loadtxt('SensitivityResults/NPR15/delta_S2.txt')
	plt.figure(figsize=(5,4))
	plt.contourf(delta_S2_NPR15, cmap='Greys')
	# plt.xticks([0.2, 1., 1.9, 3., 4., 5., 6, 7, 8], labels)
	# plt.xlim(0, 13)
	# plt.ylim(-0.5, 14)
	plt.xticks(y_pos, labels)
	plt.yticks(y_pos, labels)
	# plt.xlabel('$\\delta$ Heatmap')
	plt.gca().invert_yaxis()
	# plt.tick_params(axis='both', length=0, width=0, which='major')
	plt.colorbar()
	# plt.clim(-0.06, 0.045)
	plt.tight_layout()
	plt.savefig('SensitivityResults/' + 'heat_NPR15.pdf', bbox_inches='tight')


xs = [-6., 0., 6.]
delta_xs_NPRd = [3.63, 0., -3.88]
delta_xm_NPRd = [-8.77, 0., 9.15]

delta_xd_NPR15 = [-15.6, 0., 17.2]
delta_xm_NPR15 = [6.71, 0., -8.17]

# print (xs-xs[1])/xs[1]
plt.figure(figsize=(5,4))
plt.plot(xs, delta_xs_NPRd, linestyle='-', linewidth=1.0, color='k', marker='o', markevery=2, mfc='none')
plt.plot(xs, delta_xm_NPRd, linestyle='-', linewidth=1.0, color='k', marker='x', markevery=2)

plt.plot(xs, delta_xd_NPR15, linestyle='--', linewidth=1.0, dashes=(15,10), markevery=2, color='k', marker='s', mfc='none')
plt.plot(xs, delta_xm_NPR15, linestyle='--', linewidth=1.0, dashes=(15,10), markevery=2, color='k', marker='x')
plt.xlim(-6.5, 6.5)
plt.xticks([-6, 0, 6], ['$-$6', '0', '$+$6'])
plt.xlabel('$\\Delta x$ (%)')
plt.ylim(-20, 20)
plt.ylabel('$\\Delta\\delta$ (%)')
plt.tight_layout()
plt.savefig('SensitivityResults/' + 'delta_opt.pdf', bbox_inches='tight')
plt.show()


plt.show()