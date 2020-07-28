#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 11:53:49 2020
Script to create the density distributions using NEMD and EMD
Run inside Benchmark/x_0.2

The EMD analysis is base on plot_v_theo_sim
@author: simon
"""

import os
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.DO.EMD.density_analysis as da



        
# =============================================================================
# Main
# =============================================================================
plot_dir = "plots/2.rho_dist"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   

#folder for the small bin simulations
dir_f1= "4.1.Applying_force_0.125_bin_0.1/1"
dir_theo = "3.Measuring"

#logger.info("The data for the main data is from %s" %dir_normal_bin)
#logger.info("The data for the insert data is from %s" %dir_small_bin)

# Loading the NEMD data
fluid_f1 = da.PropertyDistribution("properties_short.dat", directory = dir_f1) 


# Loading the EMD data

solution = da.DensityDistribution("properties_short.dat", "rBulk", directory = dir_theo) 

## Loading the data for the bin 0.25 \sigma
#fluid_s = da.DensityDistribution("properties_short.dat", "rBulk", directory = dir_small_bin) 

    
# ========================================================================
# Plotting 
# =============================================================================
cf.set_plot_appearance()
plt.close("all")
# =============================================================================
# Plotting Density distribution
# =============================================================================
fig1, ax1 = plt.subplots()

ax1.plot(solution.positions, solution.data_frame["density/mass"], label = "Equilibrium")
ax1.plot(fluid_f1.positions, fluid_f1.data_frame["density/mass"], ls = '--', label = r'$\nabla \mu_s = -0.125$')


ax1.set_xlim(0, 5)
#ymin, ymax = ax1.get_ylim()
ax1.set_ylim(0, None)
ax1.set_xlabel(r'$z[\sigma] $')
ax1.set_ylabel(r'$c(z)$')
ax1.legend(loc = 'upper right')

fig1.tight_layout()
fig1.savefig('%s/densities.pdf'%plot_dir)


