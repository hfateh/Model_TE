# Haiyan Fateh
# Created on February 22, 2013

# Distribution Modules
import numpy as np
import os
import sys
import time

# User Defined Modules
cmd_folder = os.path.dirname(os.path.abspath('../Modules/hx.py'))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import te_pair
reload(te_pair)

# instantiate a te_design object
te_design = te_pair.TE_Pair()

te_design.nodes = 10

# declare materials to be used for property calculations
te_design.Ntype.material = 'MgSi'
te_design.Ptype.material = 'HMS'

te_design.pairs = 1.
te_design.Ptype.area = (3.e-3) ** 2 

fill_fraction = 0.6
length = 4.00e-3
area_ratio = 0.7
R_load_total = 0.009

te_design.fill_fraction = fill_fraction
te_design.length = length
te_design.leg_area_ratio = area_ratio
te_design.R_load_total = R_load_total

te_design.set_constants()

te_design.T_c_conv = 300.  # cold side convection temperature (K)
te_design.T_h_conv = 680.  # hot side convection temperature (K)

te_design.U_cold = 8000.
# cold side overall heat transfer coeffcient (kW / (m ** 2 * K))
te_design.U_hot = 2000.
# hot side overall heat transfer coeffcient (kW / (m ** 2 * K))

te_design.optimize()

# ================================
# postprocessing after optimization
# ================================

R_internal = te_design.R_internal
R_load_total = te_design.R_load_total
leg_area_ratio = te_design.leg_area_ratio
length = te_design.length
R_ratio = te_design.R_ratio

SIZE = 50

R_load_total_array = (
    np.linspace(0.5, 2, SIZE) * te_design.R_load_total
    )
R_ratio_array = (
    np.linspace(0.5, 2, SIZE) * R_ratio
    )
leg_area_ratio_array = (
    np.linspace(0.5, 2, SIZE) * leg_area_ratio
    )
length_array = (
    np.linspace(0.5, 2, SIZE) * length
    )

# Variables paired
power_R_area = np.zeros(
    [R_load_total_array.size, leg_area_ratio_array.size]
    )
power_R_length = np.zeros(
    [R_load_total_array.size, length_array.size]
    )
power_area_length = np.zeros(
    [leg_area_ratio_array.size, length_array.size]
    )

#print "R_area ", power_R_area.size, power_R_area.shape
#print "R_length ", power_R_length.size, power_R_length.shape
#print "area_length ", power_area_length.size, power_area_length.shape

t0 = time.clock()
for index in np.ndindex(R_load_total_array.size, leg_area_ratio_array.size):
    i = index[0]
    j = index[1]
    if j == 0:
        print "i =", i

    te_design.R_load_total = R_load_total_array[i]
    te_design.leg_area_ratio = leg_area_ratio_array[j]
    te_design.set_constants()

    te_design.solve_te_pair()

    power_R_area[i, j] = te_design.P_flux

te_design.R_load_total = R_load_total
te_design.leg_area_ratio = leg_area_ratio
te_design.length = length

t1 = time.clock() - t0
print "t1 =", t1
t0 = time.clock()
for index in np.ndindex(R_load_total_array.size, length_array.size):
    j = index[0]
    k = index[1]
    if k == 0:
        print "j =", j

    te_design.R_load_total = R_load_total_array[j]
    te_design.length = length_array[k]
    te_design.set_constants()

    te_design.solve_te_pair()

    power_R_length[j, k] = te_design.P_flux

te_design.R_load_total = R_load_total
te_design.leg_area_ratio = leg_area_ratio
te_design.length = length

t2 = time.clock() - t0
print "t2 =", t2
for index in np.ndindex(leg_area_ratio_array.size, length_array.size):
    k = index[0]
    i = index[1]
    if i == 0:
        print "k =", k

    te_design.leg_area_ratio = leg_area_ratio_array[k]
    te_design.length = length_array[i]
    te_design.set_constants()

    te_design.solve_te_pair()

    power_area_length[k, i] = te_design.P_flux

te_design.R_load_total = R_load_total
te_design.leg_area_ratio = leg_area_ratio
te_design.length = length

data_dir = '../Output/te_design/'
np.save(data_dir + 'power_R_area', power_R_area)
np.save(data_dir + 'power_R_length', power_R_length)
np.save(data_dir + 'power_area_length', power_area_length)
np.save(data_dir + 'R_load_total_array', R_load_total_array)
np.save(data_dir + 'leg_area_ratio_array', leg_area_ratio_array)
np.save(data_dir + 'length_array', length_array)
np.save(data_dir + 'R_ratio_array', R_ratio_array)

print "\nProgram finished."
print "\nPlotting..."

execfile('plot_te_design.py')

















































# This code works fine for R_load_total instead of R_ratio




# # Haiyan Fateh
# # Created on February 22, 2013

# # Distribution Modules
# import numpy as np
# import os
# import sys
# import time

# # User Defined Modules
# cmd_folder = os.path.dirname(os.path.abspath('../Modules/hx.py'))
# if cmd_folder not in sys.path:
#     sys.path.insert(0, cmd_folder)

# import te_pair
# reload(te_pair)

# # instantiate a te_design object
# te_design = te_pair.TE_Pair()

# te_design.nodes = 10

# # declare materials to be used for property calculations
# te_design.Ntype.material = 'MgSi'
# te_design.Ptype.material = 'HMS'

# te_design.pairs = 1.
# te_design.Ptype.area = (3.e-3) ** 2 

# fill_fraction = 0.6
# length = 4.00e-3
# area_ratio = 0.7
# R_load_total = 0.009

# te_design.fill_fraction = fill_fraction
# te_design.length = length
# te_design.leg_area_ratio = area_ratio
# te_design.R_load_total = R_load_total

# te_design.set_constants()

# te_design.T_c_conv = 300.  # cold side convection temperature (K)
# te_design.T_h_conv = 680.  # hot side convection temperature (K)

# te_design.U_cold = 8000.
# # cold side overall heat transfer coeffcient (kW / (m ** 2 * K))
# te_design.U_hot = 2000.
# # hot side overall heat transfer coeffcient (kW / (m ** 2 * K))

# te_design.optimize()

# # ================================
# # postprocessing after optimization
# # ================================

# R_internal = te_design.R_internal
# R_load_total = te_design.R_load_total
# leg_area_ratio = te_design.leg_area_ratio
# length = te_design.length

# R_ratio = (te_design.R_load_total/te_design.R_internal)

# SIZE = 50
# R_load_total_array = (
#     np.linspace(0.5, 2, SIZE) * te_design.R_load_total
#     )
# leg_area_ratio_array = (
#     np.linspace(0.5, 2, SIZE) * te_design.leg_area_ratio
#     )
# length_array = (
#     np.linspace(0.5, 2, SIZE) * te_design.length
#     )

# # Variables paired
# power_R_area = np.zeros(
#     [R_load_total_array.size, leg_area_ratio_array.size]
#     )
# power_R_length = np.zeros(
#     [R_load_total_array.size, length_array.size]
#     )
# power_area_length = np.zeros(
#     [leg_area_ratio_array.size, length_array.size]
#     )

# print "R_area ", power_R_area.size, power_R_area.shape
# print "R_length ", power_R_length.size, power_R_length.shape
# print "area_length ", power_area_length.size, power_area_length.shape


# t0 = time.clock()
# for index in np.ndindex(R_load_total_array.size, leg_area_ratio_array.size):
#     i = index[0]
#     j = index[1]
#     if j == 0:
#         print "i =", i

#     te_design.R_load_total = R_load_total_array[i]
#     te_design.leg_area_ratio = leg_area_ratio_array[j]
#     te_design.set_constants()

#     te_design.solve_te_pair()

#     power_R_area[i, j] = te_design.P_flux


# te_design.R_load_total = R_load_total
# te_design.leg_area_ratio = leg_area_ratio
# te_design.length = length

# t1 = time.clock() - t0
# print "t1 =", t1
# t0 = time.clock()
# for index in np.ndindex(R_load_total_array.size, length_array.size):
#     j = index[0]
#     k = index[1]
#     if k == 0:
#         print "j =", j

#     te_design.R_load_total = R_load_total_array[j]
#     te_design.length = length_array[k]
#     te_design.set_constants()

#     te_design.solve_te_pair()

#     power_R_length[j, k] = te_design.P_flux

# te_design.R_load_total = R_load_total
# te_design.leg_area_ratio = leg_area_ratio
# te_design.length = length

# t2 = time.clock() - t0
# print "t2 =", t2
# for index in np.ndindex(leg_area_ratio_array.size, length_array.size):
#     k = index[0]
#     i = index[1]
#     if i == 0:
#         print "k =", k

#     te_design.leg_area_ratio = leg_area_ratio_array[k]
#     te_design.length = length_array[i]
#     te_design.set_constants()

#     te_design.solve_te_pair()

#     power_area_length[k, i] = te_design.P_flux

# te_design.R_load_total = R_load_total
# te_design.leg_area_ratio = leg_area_ratio
# te_design.length = length

# data_dir = '../Output/te_design/'
# np.save(data_dir + 'power_R_area', power_R_area)
# np.save(data_dir + 'power_R_length', power_R_length)
# np.save(data_dir + 'power_area_length', power_area_length)
# np.save(data_dir + 'R_load_total_array', R_load_total_array)
# np.save(data_dir + 'leg_area_ratio_array', leg_area_ratio_array)
# np.save(data_dir + 'length_array', length_array)

# print "\nProgram finished."
# print "\nPlotting..."

# execfile('plot_te_design.py')
