#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 19:48:06 2020
This script is to predict the Thermal conductivity
data from the fitted coefficients in the SI 
@author: simon
"""
import numpy as np
import matplotlib.pyplot as plt
import copy
import fitting_functions as ff
import os 
import sys
Utilities_path = os.path.join(os.path.dirname(__file__), '../../')
sys.path.append(Utilities_path) #This falls into Utilities path
import Lammps.core_functions as cf

# =============================================================================
# Input parameters
# =============================================================================
coefficient_file = 'fit_coefficients_thermal_conductivity_liquid_rc2.0.txt'
values_file = 'thermal_conductivity_liquid_rc2.0.txt'

# =============================================================================
# Reading files
# =============================================================================
# Reading the coefficients and exponents
m_values, n_values, coefficients = ff.read_coeff_file(coefficient_file)

# Reading the energy
header, data = ff.read_file(values_file)

# Definining the polynomial
poly_d = ff.Polynomial(n_values, m_values,[1],[1],[1,0],[1,0])

# Check the header varible but in general the variables are like this
x_e = data[:,0] # Density
y_e = data[:,1] # Temperature
z_e = data[:,2] # Diffusion

variables = copy.deepcopy(data[:,0:2])

# Getting beta
variables[:,1] = 1/variables[:,1]

# Evaluating the polynomial at the given ponts
z_predicted = poly_d.evaluate(np.transpose(variables), coefficients)

# =============================================================================
# Plotting the points and predicted
# =============================================================================
cf.set_plot_appearance()
plt.close('all')
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter( y_e, x_e, z_e, zdir='z',marker='.',label="Simulation",color='r')
ax1.scatter(y_e, x_e, z_predicted, zdir='z',marker='.',label="Fitting",color='b')
ax1.set_xlabel(r'$T$', labelpad = 5)
ax1.set_ylabel(r'$\rho$', labelpad = 5)
ax1.set_zlabel(r'$\kappa$', labelpad = 5)
fig1.show()