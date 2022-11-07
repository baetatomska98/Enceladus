#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 08:00:27 2019

@author: akaragiannis
"""

import numpy as np
import matplotlib.pyplot as plt

filename = "baseline/baseline.csv"
#filename1 = "L15/L15.csv"
filename2 = "L150/L150.csv"
filename_isentr = "baseline/baseline_isentr.csv"

# Set normalisation parameters
# P_0 = 25000.0 # [Pa], stagnation pressure at inlet
# x_length = 0.15 # [m], nozzle length

# Read experimental and simulation data from respective files
data = np.transpose(np.genfromtxt(filename, skip_header = 1, delimiter=','))
#data1 = np.transpose(np.genfromtxt(filename1, skip_header = 1, delimiter=','))
data2 = np.transpose(np.genfromtxt(filename2, skip_header = 1, delimiter=','))
data_isentr = np.transpose(np.genfromtxt(filename_isentr, skip_header = 1, delimiter=','))

# Assign experimental pressure data and x coordinate to separate arrays
Mach_array = data[0, :]
P_array = data[4, :]
x_array = data[1, :]/1.5
#Mach1_array = data1[0, :]
#P1_array = data1[4, :]
#x1_array = data1[1, :]/15
Mach2_array = data2[0, :]
P2_array = data2[4, :]
x2_array = data2[1, :]/150

# Same for simulation data
# Adjust simulation data to match format of the Moore et al (1973) paper plots
# -0.25 to correct for different coordinate starting points
isentr_Mach_array = data_isentr[0, :]
isentr_P_array = data_isentr[4, :]
isentr_x_array = data_isentr[1, :]/1.5


###############
### PLOTTING
###############

fig, ax1 = plt.subplots()

# Color of data in the plot
colors = ['k'] # Black

# plot - multiphase pressure
plot_pressure, = ax1.plot(x_array, P_array, color = 'k',
                               label = '$L=1.5 m$ ($P/P_0$)')
#plot_pressure1, = ax1.plot(x1_array, P1_array, color = 'k', linestyle=":",
                               #label = '$L=15 m$ ($P/P_{0}$)')
plot_pressure2, = ax1.plot(x2_array, P2_array, color = 'k', linestyle="-.",
                               label = '$L=150 m$ ($P/P_0$)')
# plot - isentropic pressure
plot_pressure_isentr, = ax1.plot(isentr_x_array, isentr_P_array, color = 'k', linestyle="--", label = 'Isentropic model ($P/P_0$)')

# Plot frame
ax1.set_xlim((0, 1))
ax1.set_ylim((0.05, 1.0))
x0,x1 = ax1.get_xlim()
y0,y1 = ax1.get_ylim()

# Text box on figure's upper left corner
# props = dict(boxstyle='round', facecolor='white', alpha=1.0)
# textstr = '\n'.join((
    # r'Baseline channel',
    # r'$P_{0} = 611.2 $ Pa',
    # r'$T_{0} = 273.16 $ K',
    # r'Supersonic outlet'))
# ax1.text(0.03, 0.97, textstr, transform=ax1.transAxes, fontsize=10,
        # verticalalignment='top', bbox=props)

# Axes labels
ax1.set_xlabel('$X/Length$ [-]', fontsize = 14)
ax1.set_ylabel('$P/P_{0}$ [-]', fontsize = 14)
ax1.tick_params(axis = 'x', labelsize = 12, which = 'both',
                width = 1.0, length = 5.0)
ax1.tick_params(axis = 'y', labelsize = 12, which = 'both',
                width = 1.0, length = 5.0)


# Grid
ax1.grid(b=True, which='major', color='k', linestyle='--', alpha = 0.2)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
plot_Mach, = ax2.plot(x_array, Mach_array, 
                        color = 'g', label = '$L=1.5 m$ ($Mach$)')
#plot_Mach1, = ax2.plot(x1_array, Mach1_array, 
                        #color = 'g', linestyle=":", label = '$L=15 m$ ($Mach$)')
plot_Mach2, = ax2.plot(x2_array, Mach2_array, 
                        color = 'g', linestyle="-.", label = '$L=150 m$ ($Mach$)')
plot_Mach_isentr, = ax2.plot(isentr_x_array, isentr_Mach_array, 
                        color = 'g', linestyle="--", label = 'Isentropic model ($Mach$)')
ax1.set_adjustable("datalim")

# Switch to ax2 for the radius plot on the right vertical axis
ax2.set_xlim((0, 1))
ax2.set_ylim((0.2, 2))
x2,x3 = ax2.get_xlim()
y2,y3 = ax2.get_ylim()

ax2.set_ylabel('Mach number [-]', fontsize = 14)
ax2.tick_params(axis = 'y', labelsize = 12, which = 'major', 
                width = 1.0, length = 5.0)

# Legend on figure's upper right corner
plots = [plot_pressure, plot_Mach, plot_pressure2, plot_Mach2, plot_pressure_isentr, plot_Mach_isentr]
ax1.legend(plots, [plots_.get_label() for plots_ in plots],
           loc= 'lower left', edgecolor='k', fontsize = 10)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

# Output figure as scaled .pdf
plt.savefig('L150/L150.jpg', bbox_inches = 'tight', pad_inches = 0)
plt.show()
