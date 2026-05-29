#############################
# KHALED YOUNES
# University of Waterloo
# kyounes@uwaterloo.ca
#############################

import numpy as np
import matplotlib.pyplot as plt
import time
from datetime import timedelta
import pdb

start_time = time.time()

#############################
###### GLOBAL NUMERICS ######
#############################
error = 1
tol = 10**-6
M1guess = 2.3
ls_old = 0

# Accepted options are Bell and conical.
nozzle_profile = "Bell"

# Currently only circular injection is supported.
inj_type = "circular"

# Load user-defined injection geometry and properties.
geometry = np.loadtxt('RD_Properties.txt', delimiter=',', comments='#')

#############################
### INITIAL CONDITIONS ######
#############################
gamma1 = geometry[0]
gammaj = geometry[1]
Po1 = geometry[2]*10**6
SPR = geometry[3]
NPR = geometry[4]
xt = geometry[5]/1000.
xm = geometry[6]*xt
rth = geometry[7]/1000.
Me = geometry[8]
Cd = geometry[9]

# Run appropriate solver.
if inj_type == "2D":
	import RD_TwoDInjection

if inj_type == "slot":
	import RD_SlotInjection

if inj_type == "circular":
	Dj = geometry[10]/1000.
	import RD_CircInjection

elapsed_time_secs = time.time() - start_time

print "Execution of solver took: %s secs (Wall clock time)" % timedelta(seconds=round(elapsed_time_secs))