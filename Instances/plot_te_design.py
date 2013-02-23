# Haiyan Fateh
# Created on February 23, 2013

# Distribution Modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import os
import sys

# User Defined Modules
cmd_folder = os.path.dirname(os.path.abspath('../Modules/hx.py'))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
data_dir = '../Output/te_design/'
power_R_area = np.load(data_dir + 'power_R_area.npy')
power_R_length = np.load(data_dir + 'power_R_length.npy')
power_area_length = np.load(data_dir + 'power_area_length.npy')
R_load_total_array = np.load(data_dir + 'R_load_total_array.npy')
leg_area_ratio_array = np.load(data_dir + 'leg_area_ratio_array.npy')
length_array = np.load(data_dir + 'length_array.npy')
R_ratio_array = np.load(data_dir + 'R_ratio_array.npy')


# Plot configuration
FONTSIZE = 14
plt.rcParams['axes.labelsize'] = FONTSIZE
plt.rcParams['axes.titlesize'] = FONTSIZE
plt.rcParams['legend.fontsize'] = FONTSIZE
plt.rcParams['xtick.labelsize'] = FONTSIZE
plt.rcParams['ytick.labelsize'] = FONTSIZE
plt.rcParams['lines.linewidth'] = 1.5

plt.close()
save_dir = "../Plot/plot_te_design/"

# fig = plt.figure()
# ax = fig.gca(projection='3d')

X1, Y1 = np.meshgrid(R_ratio_array, leg_area_ratio_array)
Z1 = power_R_area.T
# cset = ax.contourf(X1, Y1, Z1, zdir='z', offset=length_array[0])

X2, Y2 = np.meshgrid(R_ratio_array, length_array * 1.e3)
Z2 = power_R_length.T
# cset = ax.contourf(X2, Y2, Z2, zdir='x', offset=current_array[0])

X3, Y3 = np.meshgrid(leg_area_ratio_array, length_array * 1.e3)
Z3 = power_area_length.T
# cset = ax.contourf(X3, Y3, Z3, zdir='x', offset=fill_array[-1])

# ax.set_xlabel('Current (A)')
# ax.set_xlim(current_array[0], current_array[-1])
# ax.set_ylabel('Fill Fraction')
# ax.set_ylim(fill_array[0], fill_array[-1])
# ax.set_zlabel('Leg Height (mm)')
# ax.set_zlim(length_array[0], length_array[-1])

RIGHT = 0.75
BOTTOM = 0.17
FRACTION = 0.15
MAX_ARRAY = np.array([Z1.max(), Z2.max(), Z3.max()])
MAX = MAX_ARRAY.max()
LEVELS = (
    (MAX + 0.1 - np.logspace(np.log10(0.1), np.log10(MAX), 12))[::-1]
    )
TICKS = LEVELS
FORMAT = '%0.2f'

fig1 = plt.figure('R_load/R_int and leg area ratio')
cset1 = plt.contourf(X1, Y1, Z1, levels=LEVELS)
CB1 = plt.colorbar(
    cset1, orientation='vertical', format=FORMAT, fraction=FRACTION,
    ticks=TICKS
    )
CB1.set_label(r'Power Flux (kW/m$^2$)')
#plt.xticks(rotation=40)
plt.xlabel('R_load/R_int')
plt.ylabel('Leg area ratio (N/P)')
plt.subplots_adjust(bottom=BOTTOM, right=RIGHT)
plt.grid()
fig1.savefig(save_dir + "power_R_area.pdf")


fig2 = plt.figure('R_load/R_int and length')
cset2 = plt.contourf(X2, Y2, Z2, levels=LEVELS)
CB2 = plt.colorbar(
    cset2, orientation='vertical', format=FORMAT, fraction=FRACTION,
    ticks=TICKS
    )
CB2.set_label(r'Power Flux (kW/m$^2$)')
#plt.xticks(rotation=40)
plt.xlabel('R_load/R_int')
plt.ylabel('Length (mm)')
plt.subplots_adjust(bottom=BOTTOM, right=RIGHT)
plt.grid()
fig2.savefig(save_dir + "power_R_length.pdf")


fig3 = plt.figure('leg area ratio and length')
cset3 = plt.contourf(X3, Y3, Z3.T, levels=LEVELS)
CB3 = plt.colorbar(
    cset3, orientation='vertical', format=FORMAT, fraction=FRACTION,
    ticks=TICKS 
    )
CB3.set_label(r'Power Flux (kW/m$^2$)')
#plt.xticks(rotation=40)
plt.xlabel('Leg area ratio (N/P)')
plt.ylabel('Length (mm)')
plt.subplots_adjust(bottom=BOTTOM, right=RIGHT)
plt.grid()
fig3.savefig(save_dir + "power_area_length.pdf")


# fig4 = plt.figure('wackiness')
# cset4 = plt.contourf(X2 / Y2, Y3, Z3.T, levels=LEVELS)
# CB4 = plt.colorbar(
#     cset4, orientation='vertical', format=FORMAT, fraction=FRACTION,
#     ticks=TICKS 
#     )
# CB3.set_label(r'Power Flux (kW/m$^2$)')
# plt.xticks(rotation=40)
# plt.xlabel('Length / Fill Fraction (mm)')
# plt.ylabel('Current (A)')
# plt.subplots_adjust(bottom=BOTTOM, right=RIGHT)
# plt.grid()
# fig3.savefig(save_dir + "power_wacky_I.pdf")

plt.show()
