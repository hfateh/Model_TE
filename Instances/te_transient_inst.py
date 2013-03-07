# To get some transient results for meeting
# ========================================
# working code set 1
# ========================================
# ========================================
# perfectly working code is given below
# ========================================
# ========================================

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
te_pair.nodes = 10
te_pair.t_array = np.linspace(0., 2., 10)

# Materials 
te_pair.Ntype.material = 'MgSi'
te_pair.Ptype.material = 'HMS'

# Number of pairs
te_pair.pairs = 1.

# Total load resistance
# te_pair.R_load_total = 1.0/127.
# te_pair.R_load_total = 0.1171587
#te_pair.R_load_total = 0.0158
te_pair.R_load_total = 0.004

# height of te_pair legs
te_pair.length = 1.e-3

# leg_area_ratio
te_pair.leg_area_ratio = .7

# fill fraction
te_pair.fill_fraction = 0.15

# total te_pair area
#te_pair.area = (10.e-3) ** 2

# Ptype area
te_pair.Ptype.area = (3.e-3) ** 2 

# BC's
te_pair.T_c_conv = 300.
te_pair.T_h_conv = 680.

# Should come from heat exchanger
te_pair.U_cold = 8000.
te_pair.U_hot = 2000.

te_pair.solve_te_pair()

# fill_fraction = np.linspace(0.1, 1.0, 20)
# power = np.zeros(fill_fraction.size)
# efficiency = np.zeros(fill_fraction.size)

# for i in range(fill_fraction.size):
#     te_pair.fill_fraction = fill_fraction[i]
#     te_pair.set_constants()
#     te_pair.solve_te_pair()
#     power[i] = te_pair.P
#     efficiency[i] = te_pair.eta_thermal
#     print "\nfill fraction is ", te_pair.fill_fraction 
#     print "Power is ", te_pair.P
#     print "Power flux is ", te_pair.P_flux
#     print "Ntype T_h is ", te_pair.Ntype.T_h
#     print "Ptype T_h is ", te_pair.Ptype.T_h

# Steady state solution used as BS for transient


# # te_pair.optimize()

# # Change in BC's
# te_pair.T_h_conv += 100.
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

# # plt.plot(fill_fraction * 1.e2, p_flux_output)
# # plt.plot(leg.t_array, leg.Power_transient)
# # plt.plot(leg.t_array, leg.R_internal_transient)
# for i in range(te_pair.t_array.size):
#     j = i + te_pair.t_array.size
#     plt.plot(te_pair.Ntype.x * 1e3, T_xn[i, :])
#     #plt.plot(te_pair.Ptype.x * 1e3, te_pair.Ptype.qxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.Rxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.T_xt[j, :])

# plt.grid()
# plt.xlabel('Leg height (mm)')
# plt.ylabel('Heat flux (W/m2)')
# #plt.ylim(T_xn.min() - 10., T_xn.max() + 10.)
# #plt.ylim(te_pair.Ptype.Txt.min() - 10., te_pair.Ptype.Txt.max() + 10.)
# #plt.ylim(te_pair.Ptype.Txt.min() - 10., te_pair.Ptype.Txt.max() + 10.)
# #plt.xlim(-0.05, 1.05)
# plt.subplots_adjust(left=0.15)

# #plt.savefig('../Plots/leg_instance/transient.pdf')

# plt.show()
# #=============================================
# #=============================================

































































# trying to stack transient modules
# # distribution modules
# import matplotlib.pyplot as plt
# import os
# import sys
# import numpy as np

# # local user modules
# cmd_folder = os.path.dirname(os.path.abspath('../Modules/hx.py'))
# if cmd_folder not in sys.path:
#     sys.path.insert(0, cmd_folder)

# import te_pair
# reload(te_pair)

# te_pair = te_pair.TE_Pair()
# nodes = 20

# # Numeric code variables
# te_pair.nodes = nodes


# # Materials 
# te_pair.Ntype.material = 'MgSi'
# te_pair.Ptype.material = 'HMS'

# # Number of pairs
# te_pair.pairs = 127

# # Total load resistance
# te_pair.R_load_total = 1.0

# # height of te_pair legs
# te_pair.length = 5.e-3

# # leg_area_ratio
# te_pair.leg_area_ratio = .7

# # fill fraction
# te_pair.fill_fraction = 0.3

# # Ptype area
# te_pair.Ptype.area = (3.e-3) ** 2 

# # BC's
# te_pair.T_c_conv = 300.
# # te_pair.T_h_conv = 680.

# # Should come from heat exchanger
# te_pair.U_cold = 8000.
# te_pair.U_hot = 2000.

# # ===============================
# # for transient stacking solution
# # ===============================

# input_array = np.array([[0., 700.],
#                         [1., 600.],
#                         [2., 940.]])

# runs_needed = (input_array.size / 2) - 1
# T_h_conv_array = input_array[:,1]
# t_array = input_array[:,0]
# t_size = 10

# # + 1 since SS solution is also included
# row_size = 1 + (t_size * runs_needed)
# column_size = nodes

# # I know the size now, lets initiate arrays
# NTxt_array = np.zeros([row_size, column_size])
# PTxt_array = np.zeros([row_size, column_size])
# Nqxt_array = np.zeros([row_size, column_size])
# Pqxt_array = np.zeros([row_size, column_size])

# print "run this ..."
# # set SS BC
# te_pair.t_array = np.linspace(0., 5., 100)
# te_pair.T_h_conv = T_h_conv_array[0]

# # solve the pair once (SS)
# te_pair.solve_te_pair()

# # # get steady state profiles
# NTxt_array[0,:] = te_pair.Ntype.T_x
# PTxt_array[0,:] = te_pair.Ptype.T_x
# Nqxt_array[0,:] = te_pair.Ntype.q_x
# Pqxt_array[0,:] = te_pair.Ptype.q_x

# # Now get transient profiles
# for i in range(runs_needed):
#     te_pair.T_h_conv = T_h_conv_array[i+1]
#     print "T_h_conv is set to ", te_pair.T_h_conv
#     te_pair.t_array = np.linspace(t_array[i], t_array[i+1], 10)
#     te_pair.set_t_array()
#     te_pair.set_constants()
#     te_pair.solve_te_pair_transient_once()


#     start_row = 1 + i * te_pair.t_array.size
#     end_row = 1 + (i + 1) * te_pair.t_array.size

#     print "start_row ", start_row
#     print "end_row ", end_row

#     print "\nte_pair T is "
#     print te_pair.Ntype.Txt
#     print te_pair.Ntype.Txt.size
#     print te_pair.Ntype.Txt.shape


#     print "\nMy array is "
#     print NTxt_array[start_row:end_row,:]
#     print NTxt_array[start_row:end_row,:].size
#     print NTxt_array[start_row:end_row,:].shape


#     NTxt_array[start_row:end_row,:] = te_pair.Ntype.Txt
#     PTxt_array[start_row:end_row,:] = te_pair.Ptype.Txt


# # Steady state solution used as BS for transient
# # te_pair.solve_te_pair()

# # te_pair.optimize()

# # # Change in BC's
# # te_pair.T_h_conv += 20.
# # te_pair.T_c_conv += 0.

# # # Transient te_pair solution
# # te_pair.solve_te_pair_transient_once()


# #==========================================
# #==========================================
# # Plots
# #T_xn = te_pair.Ptype.Txt
# T_xn = NTxt_array

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
# for i in range(row_size):
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




































































# # This plots fill fraction vs power flux
# # ========================================
# # working code set 1
# # ========================================
# # ========================================
# # perfectly working code is given below
# # ========================================
# # ========================================



# # distribution modules
# import matplotlib.pyplot as plt
# import os
# import sys
# import numpy as np

# # local user modules
# cmd_folder = os.path.dirname(os.path.abspath('../Modules/hx.py'))
# if cmd_folder not in sys.path:
#     sys.path.insert(0, cmd_folder)

# import te_pair
# reload(te_pair)

# te_pair = te_pair.TE_Pair()

# # Numeric code variables
# te_pair.nodes = 10
# te_pair.t_array = np.linspace(0., 5., 100)

# # Materials 
# te_pair.Ntype.material = 'MgSi'
# te_pair.Ptype.material = 'HMS'

# # Number of pairs
# te_pair.pairs = 1.

# # Total load resistance
# # te_pair.R_load_total = 1.0/127.
# # te_pair.R_load_total = 0.1171587
# te_pair.R_load_total = 0.0158

# # height of te_pair legs
# te_pair.length = 4.e-3

# # leg_area_ratio
# te_pair.leg_area_ratio = .7

# # fill fraction
# # te_pair.fill_fraction = 0.22

# # Ptype area
# te_pair.Ptype.area = (3.e-3) ** 2 

# # BC's
# te_pair.T_c_conv = 300.
# te_pair.T_h_conv = 800.

# # Should come from heat exchanger
# te_pair.U_cold = 8000.
# te_pair.U_hot = 2000.

# fill_fraction = np.linspace(0.2, 1.0, 10)
# p_flux_output = np.zeros(fill_fraction.size)

# for i in range(fill_fraction.size):
#     te_pair.fill_fraction = fill_fraction[i]
#     te_pair.solve_te_pair()
#     p_flux_output[i] = te_pair.P_flux

# # Steady state solution used as BS for transient
# #te_pair.solve_te_pair()




# # te_pair.optimize()

# # # Change in BC's
# # te_pair.T_h_conv += 20.
# # te_pair.T_c_conv += 0.

# # # Transient te_pair solution
# # te_pair.solve_te_pair_transient_once()


# #==========================================
# #==========================================
# # Plots
# # T_xn = te_pair.Ptype.Txt

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

# plt.plot(fill_fraction * 1.e2, p_flux_output)
# # plt.plot(leg.t_array, leg.Power_transient)
# # plt.plot(leg.t_array, leg.R_internal_transient)
# # for i in range(te_pair.t_array.size):
#     #j = i + te_pair.t_array.size
#     #plt.plot(te_pair.Ntype.x * 1e3, T_xn[i, :])
#     #plt.plot(te_pair.Ptype.x * 1e3, te_pair.Ptype.qxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.Rxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.T_xt[j, :])

# plt.grid()
# plt.xlabel('Fill fraction (%)')
# plt.ylabel('Power_flux (kW/m2)')
# #plt.ylim(T_xn.min() - 10., T_xn.max() + 10.)
# #plt.ylim(te_pair.Ptype.Txt.min() - 10., te_pair.Ptype.Txt.max() + 10.)
# #plt.ylim(te_pair.Ptype.Txt.min() - 10., te_pair.Ptype.Txt.max() + 10.)
# #plt.xlim(-0.05, 1.05)
# plt.subplots_adjust(left=0.15)

# #plt.savefig('../Plots/leg_instance/transient.pdf')

# plt.show()
# #=============================================
# #=============================================





















































































# # This plots current and voltage
# # ========================================
# # working code set 1
# # ========================================
# # ========================================
# # perfectly working code is given below
# # ========================================
# # ========================================



# # distribution modules
# import matplotlib.pyplot as plt
# import os
# import sys
# import numpy as np

# # local user modules
# cmd_folder = os.path.dirname(os.path.abspath('../Modules/hx.py'))
# if cmd_folder not in sys.path:
#     sys.path.insert(0, cmd_folder)

# import te_pair
# reload(te_pair)

# te_pair = te_pair.TE_Pair()

# # Numeric code variables
# te_pair.nodes = 10
# te_pair.t_array = np.linspace(0., 5., 100)

# # Materials 
# te_pair.Ntype.material = 'MgSi'
# te_pair.Ptype.material = 'HMS'

# # Number of pairs
# te_pair.pairs = 1.

# # Total load resistance
# # te_pair.R_load_total = 1.0/127.
# # te_pair.R_load_total = 0.1171587
# te_pair.R_load_total = 0.0158

# # height of te_pair legs
# te_pair.length = 5.e-3

# # leg_area_ratio
# te_pair.leg_area_ratio = .7

# # fill fraction
# # te_pair.fill_fraction = 0.22

# # Ptype area
# te_pair.Ptype.area = (5.e-3) ** 2 

# # BC's
# te_pair.T_c_conv = 300.
# te_pair.T_h_conv = 680.

# # Should come from heat exchanger
# te_pair.U_cold = 8000.
# te_pair.U_hot = 2000.

# # fill_fraction = np.linspace(0.2, 1.0, 10)
# # fill_fraction = np.linspace(0.2, 1.0, 10)

# R_load_total = np.linspace(0.5, 2.0, 10) * te_pair.R_load_total
# p_flux_output = np.zeros(R_load_total.size)
# current = np.zeros(R_load_total.size)
# voltage = np.zeros(R_load_total.size)

# for i in range(R_load_total.size):
#     te_pair.R_load_total = R_load_total[i]
#     te_pair.solve_te_pair()
#     p_flux_output[i] = te_pair.P_flux
#     current[i] = te_pair.I
#     voltage[i] = te_pair.V

# # Steady state solution used as BS for transient
# #te_pair.solve_te_pair()




# # te_pair.optimize()

# # # Change in BC's
# # te_pair.T_h_conv += 20.
# # te_pair.T_c_conv += 0.

# # # Transient te_pair solution
# # te_pair.solve_te_pair_transient_once()


# #==========================================
# #==========================================
# # Plots
# # T_xn = te_pair.Ptype.Txt

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

# plt.plot(R_load_total, current)
# plt.xlabel('R_load (ohm)')
# plt.ylabel('Current (A)')

# plt.plot(R_load_total, voltage)
# # plt.plot(leg.t_array, leg.Power_transient)
# # plt.plot(leg.t_array, leg.R_internal_transient)
# # for i in range(te_pair.t_array.size):
#     #j = i + te_pair.t_array.size
#     #plt.plot(te_pair.Ntype.x * 1e3, T_xn[i, :])
#     #plt.plot(te_pair.Ptype.x * 1e3, te_pair.Ptype.qxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.Rxt[i, :])
#     #plt.plot(leg.x * 1e3, leg.T_xt[j, :])

# plt.grid()
# plt.xlabel('Fill fraction (%)')
# plt.ylabel('Power_flux (kW/m2)')
# #plt.ylim(T_xn.min() - 10., T_xn.max() + 10.)
# #plt.ylim(te_pair.Ptype.Txt.min() - 10., te_pair.Ptype.Txt.max() + 10.)
# #plt.ylim(te_pair.Ptype.Txt.min() - 10., te_pair.Ptype.Txt.max() + 10.)
# #plt.xlim(-0.05, 1.05)
# plt.subplots_adjust(left=0.15)

# #plt.savefig('../Plots/leg_instance/transient.pdf')

# plt.show()
# #=============================================
# #=============================================

