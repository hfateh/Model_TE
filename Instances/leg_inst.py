import matplotlib.pyplot as plt
import numpy as np
import os
import sys

cmd_folder = os.path.dirname('../Modules/')
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import leg
reload(leg)

leg = leg.Leg()

# leg.length = 1.e-4
# leg.material = 'HMS'
# leg.nodes = 10 
# leg.T_h_conv = 750.
# leg.T_c_conv = 360.
# leg.U_hot = 54.e2
# leg.U_cold = 253.e3

# R_load = 0.0033
# leg.R_load = R_load
# leg.set_constants()
# leg.set_I()
# leg.solve_leg_for_real()

leg.solve_leg_for_real()





# R_load = np.linspace(0.001, 0.007, 12)
# R_internal = np.zeros(R_load.size)
# P = np.zeros(R_load.size)
# for i in range(R_load.size):
#     leg.R_load = R_load[i]
#     leg.set_I()
#     leg.set_constants()
#     leg.solve_leg_for_real()
#     P[i] = leg.P[0]
#     R_internal[i] = leg.R_internal

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
