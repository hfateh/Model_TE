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

te_pair.T_c_conv = 323.
te_pair.T_h_conv = 443.

te_pair.U_cold = 800000.
te_pair.U_hot = 20000.

te_pair.solve_te_pair()



R_load_total = np.linspace(0.1, 1.0, 30)
R_internal = np.zeros(R_load_total.size)
P = np.zeros(R_load_total.size)
for i in range(R_load_total.size):
    te_pair.R_load_total = R_load_total[i]
    # it has to be called for this to work
    # te_pair.set_R_load()
    te_pair.solve_te_pair()
    P[i] = te_pair.P
    R_internal[i] = te_pair.R_internal
    print "\nR_load_total is ", R_load_total[i]
    print "R_internal is ", R_internal[i]
    print "Power is ", P[i]

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

plt.plot(R_load_total, P, 's')
# plt.plot(R_load, R_internal, 'o')
plt.grid()
plt.xlabel('R_load_total (ohms)')
plt.ylabel('Power (W)')
#plt.ylim(0,0)
#plt.xlim(0,0)

plt.show()









# R_load = np.linspace(0.001, 1.0, 30)
# R_internal = np.zeros(R_load.size)
# P = np.zeros(R_load.size)
# for i in range(R_load.size):
#     te_pair.R_load = R_load[i]
#     te_pair.solve_te_pair()
#     #print "Actual J is ", te_pair.J
#     P[i] = te_pair.P
#     print "Power is ", P[i]
#     R_internal[i] = te_pair.R_internal

# # Plot configuration
# FONTSIZE = 14
# plt.rcParams['axes.labelsize'] = FONTSIZE
# plt.rcParams['axes.titlesize'] = FONTSIZE
# plt.rcParams['legend.fontsize'] = FONTSIZE
# plt.rcParams['xtick.labelsize'] = FONTSIZE
# plt.rcParams['ytick.labelsize'] = FONTSIZE
# plt.rcParams['lines.linewidth'] = 1.5
# plt.rcParams['lines.markersize'] = 10

# plt.close()

# plt.figure()

# plt.plot(R_load, P, 's')
# # plt.plot(R_load, R_internal, 'o')
# plt.grid()
# plt.xlabel('R_load (ohms)')
# plt.ylabel('Power (W)')
# #plt.ylim(0,0)
# #plt.xlim(0,0)

# plt.show()

