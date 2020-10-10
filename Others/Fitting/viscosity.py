#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 12:50:23 2020
This script is to predict the energy data from the fitted coefficients in the SI for 
the GLJ paper
@author: simon
"""
import numpy as np

import multid_fit_liquid_jup as mf
import matplotlib.pyplot as plt
import copy

coefficient_file = 'fit_coefficients_liquid_rc2.0.txt'
values_file = 'viscosity_liquid_rc2.0.txt'


# Reading the coefficients and exponents
m_values, n_values, coefficients = mf.read_coeff_file(coefficient_file)

# Reading the energy
header, data = mf.read_file(values_file)

poly_eta = mf.polynomial(n_values, m_values,[1],[1],[1,0],[1,0])

x_e = data[:,0] # Density
y_e = data[:,1] # Temperature
z_e = data[:,2] # Viscosity

## Getting the volume of the system
#
#n_part = 1000
##vol = n_part/y_e[0]

variables = copy.deepcopy(data[:,0:2])

# Getting beta
variables[:,1] = 1/variables[:,1]

er_results_eta = mf.test_prediction(coefficients, np.transpose(variables), data[:,2], poly_eta)

z_predicted = er_results_eta[:,1]

mf.set_plot_appearance()

#plt.close('all')
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter( y_e, x_e, z_e, zdir='z',marker='.',label="Simulation",color='r')
ax1.scatter(y_e, x_e, z_predicted, zdir='z',marker='.',label="Fitting",color='b')
ax1.set_xlabel(r'$\rho$', labelpad = 5)
ax1.set_ylabel(r'$T$', labelpad = 5)
ax1.set_zlabel(r'$\eta$', labelpad = 5)




# Testing the polynomial by hand



