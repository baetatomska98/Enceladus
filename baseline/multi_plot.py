# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 15:42:25 2022

@author: tomas
"""

import numpy as np
import matplotlib.pyplot as plt

filename = "L150/L150_multi.csv"
#filename1 = "L15/L15_multi.csv"
filename_baseline = "baseline/baseline_multi.csv"

# Read experimental and simulation data from respective files
data = np.transpose(np.genfromtxt(filename, skip_header = 1, delimiter=','))
#data1 = np.transpose(np.genfromtxt(filename1, skip_header = 1, delimiter=','))
data_baseline = np.transpose(np.genfromtxt(filename_baseline, skip_header = 1, delimiter=','))

# Assign experimental pressure data and x coordinate to separate arrays
J = data[0, :]
N = data[2, :]
x = data[3, :]/150
S = data[6, :]
Y = data[7, :]
drdt = data[8, :]
r = data[9, :]
#J1 = data1[0, :]
#N1 = data1[2, :]
#x1 = data1[3, :]/15
#S1 = data1[6, :]
#Y1 = data1[7, :]
#drdt1 = data1[8, :]
#r1 = data1[9, :]
J_baseline = data_baseline[0, :]
N_baseline = data_baseline[2, :]
x_baseline = data_baseline[3, :]/1.5
S_baseline = data_baseline[6, :]
Y_baseline = data_baseline[7, :]
drdt_baseline = data_baseline[8, :]
r_baseline = data_baseline[9, :]

plt.figure()
plt.yscale('log')
plt.ylim((1e13, 1e23))
plt.plot(x,J,label=r'$L=150 m$')
#plt.plot(x1,J1,label=r'$L=15 m$')
plt.plot(x_baseline,J_baseline,label='Baseline')
plt.xlabel(r'$x/Length [-]$')
plt.ylabel(r'$J$ [# $m^{-3} s^{-1}]$')
plt.legend()

plt.figure()
plt.yscale('log')
plt.ylim((1e12, 1e22))
plt.plot(x,N,label=r'$L=150 m$')
#plt.plot(x1,N1,label=r'$L=15 m$')
plt.plot(x_baseline,N_baseline,label='Baseline')
plt.xlabel(r'$x/Length [-]$')
plt.ylabel(r'$N$ [#$/kg]$')
plt.legend()

plt.figure()
plt.plot(x,S,label=r'$L=150 m$')
#plt.plot(x1,S1,label=r'$L=15 m$')
plt.plot(x_baseline,S_baseline,label='Baseline')
plt.xlabel(r'$x/Length [-]$')
plt.ylabel(r'$S [-]$')
plt.legend()

plt.figure()
plt.plot(x,Y,label=r'$L=150 m$')
#plt.plot(x1,Y1,label=r'$L=15 m$')
plt.plot(x_baseline,Y_baseline,label='Baseline')
plt.xlabel(r'$x/Length [-]$')
plt.ylabel(r'$Y [-]$')
plt.legend()

plt.figure()
plt.plot(x,drdt,label=r'$L=150 m$')
#plt.plot(x1,drdt1,label=r'$L=15 m$')
plt.plot(x_baseline,drdt_baseline,label='Baseline')
plt.xlabel(r'$x/Length [-]$')
plt.ylabel(r'$\frac{dr_{droplet}}{dt} [m/s]$')
plt.legend()

plt.figure()
plt.yscale('log')
plt.ylim((1e-10, 1e-6))
plt.plot(x,r,label=r'$L=150 m$')
#plt.plot(x1,r1,label=r'$L=15 m$')
plt.plot(x_baseline,r_baseline,label='Baseline')
plt.xlabel(r'$x/Length [-]$')
plt.ylabel(r'$r_{droplet} [m]$')
plt.legend()

plt.show()
