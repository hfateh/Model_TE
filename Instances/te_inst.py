# distribution modules
import matplotlib.pyplot as plt
import os
import sys
import numpy as np

# local user modules
cmd_folder = os.path.dirname(os.path.abspath('../Modules/hx.py'))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import te_pair
reload(te_pair)

te_pair = te_pair.TE_Pair()

te_pair.Ntype.material = 'MgSi'
te_pair.Ptype.material = 'HMS'

te_pair.T_c_conv = 300.
te_pair.T_h_conv = 800.

te_pair.U_cold = 800000.
te_pair.U_hot = 20000.

# te_pair.solve_te_pair()



R_load = np.linspace(0.001, 0.03, 30)
R_internal = np.zeros(R_load.size)
P = np.zeros(R_load.size)
for i in range(R_load.size):
    te_pair.R_load = R_load[i]
    te_pair.solve_te_pair()
    #print "Actual J is ", te_pair.J
    P[i] = te_pair.P
    print "Power is ", P[i]
    R_internal[i] = te_pair.R_internal

# Plot configuration
FONTSIZE = 14
plt.rcParams['axes.labelsize'] = FONTSIZE
plt.rcParams['axes.titlesize'] = FONTSIZE
plt.rcParams['legend.fontsize'] = FONTSIZE
plt.rcParams['xtick.labelsize'] = FONTSIZE
plt.rcParams['ytick.labelsize'] = FONTSIZE
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['lines.markersize'] = 10

plt.close()

plt.figure()

plt.plot(R_load, P, 's')
# plt.plot(R_load, R_internal, 'o')
plt.grid()
plt.xlabel('R_load (ohms)')
plt.ylabel('Power (W)')
#plt.ylim(0,0)
#plt.xlim(0,0)

plt.show()









# print "\nLoad resistance is ", te_pair.R_load
# print "\nNtype T distribution is \n", te_pair.Ntype.T_x
# print "\nPtype T distribution is \n", te_pair.Ptype.T_x

# print "\nte_pair q_h is ", te_pair.q_h
# print "\nte_pair q_c is ", te_pair.q_c
# print "\nPower output for te_pair is ", te_pair.P
