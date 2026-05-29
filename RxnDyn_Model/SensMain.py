from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np

SubroutineOne = 0
SubroutineTwo = 1

if SubroutineOne == 1:
	import SensSep
	# Define the model inputs
	problem = {
	    'num_vars': 9,
	    'names': ['gamma1', 'gammaj', 'NPR', 'SPR', 'Cd', 'alpha', 'xm', 'rth', 'Dj'],
	    'bounds': [[1.3, 1.4], [1.4, 1.66], [10.0, 40.0], [0.9, 1.1], [0.765, 0.935], 
	    [4.9, 6.0], [0.81*0.1, 0.99*0.1], [8.7*10**-3, 11.0*10**-3], 
	    [5.4*10**-3, 6.6*10**-3]]
	}

	# Generate samples
	param_values = saltelli.sample(problem, 500)
	print len(param_values)

	dummy = SensSep.SepDistance(param_values)
	xs = dummy[0]
	h =  dummy[1]
	Pp = dummy[2]

	# Perform analysis
	Si1 = sobol.analyze(problem, xs, print_to_console=True)
	Si2 = sobol.analyze(problem, h)
	Si3 = sobol.analyze(problem, Pp)

	np.savetxt('SensitivityResults/Xs_S1.txt', Si1['S1'])
	np.savetxt('SensitivityResults/h_S1.txt', Si2['S1'])
	np.savetxt('SensitivityResults/Pp_S1.txt', Si3['S1'])

	np.savetxt('SensitivityResults/Xs_ST.txt', Si1['ST'])
	np.savetxt('SensitivityResults/h_ST.txt', Si2['ST'])
	np.savetxt('SensitivityResults/Pp_ST.txt', Si3['ST'])

	np.savetxt('SensitivityResults/Xs_S2.txt', Si1['S2'])
	np.savetxt('SensitivityResults/h_S2.txt', Si2['S2'])
	np.savetxt('SensitivityResults/Pp_S2.txt', Si3['S2'])


if SubroutineTwo == 1:
	import SensForces
	# Define the model inputs
	problem = {
	    'num_vars': 14,
	    'names': ['rth', 're', 'Dj', 
	    'Po1', 'Poj', 
	    'gamma1', 'gammaj', 
	    'xm', 'xt',
	    'Cd', 
	    'xs', 'xd',
	    'Pp', 'psi'],
	    'bounds': [[9.67/1000, 9.77/1000], [20.08/1000, 20.18/1000], [5.95/1000, 6.05/1000], 
	    [3.749*10**6, 3.751*10**6], [3.749*10**6, 3.751*10**6], 
	    [1.32, 1.401], [1.32, 1.401], 
	    [0.085, 0.095], [0.095, 0.105], 
	    [0.8415, 0.8585], 
	    [0.07426, 0.08374], [0.0846, 0.0954], 
	    [3.525*10**5, 3.975*10**5], [69.5, 78.3]]
	    }

	# Generate samples
	Ns = 256
	param_values = saltelli.sample(problem, Ns)
	print len(param_values)

	delta = SensForces.Forces(param_values)

	# Perform analysis
	Si1 = sobol.analyze(problem, delta, print_to_console=True)

	np.savetxt('SensitivityResults/NPR15/Convergence/delta_ST_' + \
	str(Ns) + '.txt', Si1['ST'])
	np.savetxt('SensitivityResults/NPR15/delta_S1.txt', Si1['S1'])
	np.savetxt('SensitivityResults/NPR15/delta_S2.txt', Si1['S2'])
	np.savetxt('SensitivityResults/NPR15/delta_ST.txt', Si1['ST'])
	
