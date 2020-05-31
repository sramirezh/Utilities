#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 11:53:49 2020
Script to create the plot of the velocity profile with the insert close to the wall
Run this inside Benchmark/x_0.2
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
logger = cf.log(__file__, os.getcwd())    


#folder for the small bin simulations
dir_small_bin = "5.Applying_force_p_0.001_bin_0.1/1"
dir_normal_bin = "5.Applying_force_p_0.001"


#dir_small_bin = "4.Applying_force_0.125_bin_0.1/1"
#dir_normal_bin = "4.Applying_force_0.125/1"

logger.info("The data for the main data is from %s" %dir_normal_bin)
logger.info("The data for the insert data is from %s" %dir_small_bin)

# Loading the data for the bin 0.1 \sigma
fluid_n = da.DensityDistribution("properties_short.dat", "rBulk", directory = dir_normal_bin) 
solute_n = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = dir_normal_bin) 
solvent_n = da.DensityDistribution("Fproperties_short.dat", "rBulk", directory = dir_normal_bin) 

# Loading the data for the bin 0.25 \sigma
fluid_s = da.DensityDistribution("properties_short.dat", "rBulk", directory = dir_small_bin) 
solute_s = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = dir_small_bin) 
solvent_s = da.DensityDistribution("Fproperties_short.dat", "rBulk", directory = dir_small_bin) 

# =============================================================================
# Plotting the velocities
# =============================================================================
cf.set_plot_appearance()

plt.close("all")
fig1, ax1 = plt.subplots()


solute_n.plot_property_dist("vx", ax = ax1)
solvent_n.plot_property_dist("vx", ax = ax1)
fluid_n.plot_property_dist("vx", ax = ax1)


#ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')


ax1.set_xlim(0, fluid_n.positions[-1])
ax1.set_ylim(0, None)
ax1.set_xlabel(r'$z[\sigma] $')
ax1.set_ylabel(r'$v_x(z)$')
ax1.legend(["Solute", "Solvent", "Fluid"], loc = 'upper left')

## Adding the insert

left, bottom, width, height = [0.55, 0.30, 0.4, 0.30]
ax2 = fig1.add_axes([left, bottom, width, height])

solute_s.plot_property_dist("vx", ax = ax2)
solvent_s.plot_property_dist("vx", ax = ax2)
fluid_s.plot_property_dist("vx", ax = ax2)


xmin = 0
xmax = 2
ax2 = cf.plot_zoom(ax2, [xmin, xmax])
ax2.tick_params(axis='both', which='major', labelsize = 12)
#ax2.set_ylabel(r'$V(r)$',fontsize =17, labelpad=-5)
#ax2.set_xlabel(r'$r$' ,fontsize =17, labelpad=-5)


fig1.tight_layout()
fig1.savefig('vprofile_insert.pdf')


