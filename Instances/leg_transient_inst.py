# distribution modules
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# User Defined Modules
cmd_folder = os.path.dirname('../Modules/')
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
import leg
reload(leg)

leg = leg.Leg()
leg.length = 3.56e-4
# leg.I = 13.
leg.material = 'HMS'
leg.nodes = 20
leg.t_array = np.linspace(0, 1, 50)

leg.T_h_conv = 400.
leg.U_hot = 54e3
leg.T_c_conv = 300.
leg.U_cold = 253e3

#leg.t_array = np.logspace(np.log10(0.001), np.log10(1), 10)
leg.set_constants()

leg.solve_leg()

leg.T_h_conv += 300.
leg.solve_leg_transient_once()
#leg.T_h_conv -= 300.
#leg.T_c_conv += 15. 
#leg.solve_leg_transient_once()

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

plt.plot(leg.t_array, leg.I_transient)
# plt.plot(leg.t_array, leg.Power_transient)
# for i in range(leg.t_array.size):
#     j = i + leg.t_array.size
#     plt.plot(leg.x * 1e3, leg.Vsxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.Rxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.T_xt[j, :])

# plt.grid()
# plt.xlabel('Position (mm)')
# plt.ylabel('Temperature (K)')
# plt.ylim(leg.Rxt.min() - 25., leg.Rxt.max() + 25.)
# plt.xlim(-0.05, 0.40)
# plt.subplots_adjust(left=0.15)

# #plt.savefig('../Plots/leg_instance/transient.pdf')

plt.show()