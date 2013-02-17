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

# Numeric code variables
te_pair.nodes = 50
te_pair.t_array = np.linspace(0., 5., 100)

# Materials 
te_pair.Ntype.material = 'MgSi'
te_pair.Ptype.material = 'HMS'

# Number of pairs
te_pair.pairs = 127

# Total load resistance
te_pair.R_load_total = 1.0

# height of te_pair legs
te_pair.length = 5.e-3

# leg_area_ratio
te_pair.leg_area_ratio = .7

# fill fraction
te_pair.fill_fraction = 0.3

# Ptype area
te_pair.Ptype.area = (5.e-3) ** 2 

# BC's
te_pair.T_c_conv = 300.
te_pair.T_h_conv = 500.

# Should come from heat exchanger
te_pair.U_cold = 500.
te_pair.U_hot = 500.

# Steady state solution used as BS for transient
te_pair.solve_te_pair()

# te_pair.optimize()

# # Change in BC's
# te_pair.T_h_conv += 20.
# te_pair.T_c_conv += 0.

# # Transient te_pair solution
# te_pair.solve_te_pair_transient_once()


# #==========================================
# #==========================================
# # Plots
# T_xn = te_pair.Ptype.qxt

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

# # plt.plot(leg.t_array, leg.I_transient)
# # plt.plot(leg.t_array, leg.Power_transient)
# # plt.plot(leg.t_array, leg.R_internal_transient)
# for i in range(te_pair.t_array.size):
#     #j = i + te_pair.t_array.size
#     plt.plot(te_pair.Ntype.x * 1e3, T_xn[i, :])
#     #plt.plot(te_pair.Ptype.x * 1e3, te_pair.Ptype.qxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.Rxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.T_xt[j, :])

# plt.grid()
# plt.xlabel('Position (mm)')
# plt.ylabel('Temperature (K)')
# #plt.ylim(T_xn.min() - 10., T_xn.max() + 10.)
# #plt.ylim(te_pair.Ptype.Txt.min() - 10., te_pair.Ptype.Txt.max() + 10.)
# #plt.ylim(te_pair.Ptype.Txt.min() - 10., te_pair.Ptype.Txt.max() + 10.)
# #plt.xlim(-0.05, 1.05)
# plt.subplots_adjust(left=0.15)

# #plt.savefig('../Plots/leg_instance/transient.pdf')

# plt.show()
# #=============================================
# #=============================================

