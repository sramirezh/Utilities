#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 12:50:23 2020
This script is to predict the energy data from the fitted coefficients in the SI for 
the GLJ paper
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
plot_dir = "plots"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)    

coefficient_file = 'fit_coefficients_liquid_rc2.0.txt'
values_file = 'viscosity_liquid_rc2.0.txt'

# =============================================================================
# Reading files
# =============================================================================
# Reading the coefficients and exponents
m_values, n_values, coefficients = ff.read_coeff_file(coefficient_file)

# Reading the energy
header, data = ff.read_file(values_file)

# Definining the polynomial
poly_eta = ff.Polynomial(n_values, m_values,[1],[1],[1,0],[1,0])

# Check the header varible but in general the variables are like this
x_e = data[:,0] # Density
y_e = data[:,1] # Temperature
z_e = data[:,2] # Viscosity

variables = copy.deepcopy(data[:,0:2])

# Getting beta
variables[:,1] = 1/variables[:,1]

# Evaluating the polynomial at the given ponts
z_predicted = poly_eta.evaluate(np.transpose(variables), coefficients)

# =============================================================================
# Plotting the points and predicted
# =============================================================================
cf.set_plot_appearance()
plt.close('all')
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter( y_e, x_e, z_e, zdir='z',marker='.',label="Simulation",color='r')
ax1.scatter(y_e, x_e, z_predicted, zdir='z',marker='.',label="Fitting",color='b')
ax1.set_xlabel(r'$T$', labelpad = 10)
ax1.set_ylabel(r'$\rho$', labelpad = 5)
ax1.set_zlabel(r'$\eta$', labelpad = 5)
ax1.legend(loc=(0.15,0.65),frameon=0)
# To reduce the empty space around the figure
fig1.subplots_adjust(left=0, right=1, bottom=0, top=1) 
fig1.show()
fig1.savefig("%s/eta.pdf"%plot_dir)
logger.info("plotted %s/eta.pdf"%plot_dir)

