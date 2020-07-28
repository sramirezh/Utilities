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
import Lammps.lammps_utilities as lu


class SimulationInsert(lu.SimulationType):
    """
    Class to analyse the velocity profiles and create an insert"""

    def __init__(self, name, small_bin, normal_bin, multiplier, legend):
        """
        Args:
            name: Identifier of the simulation
            small_bin Direction for the simulation with small binning
            normal_bin Direction for the simulation with normal binning
            multiplier: velocity multiplier, as the gradient in some simulations
                        has a minus sign
            legend: localisation
        Other atributes:
            sampling_interval
        """
        self.name = name
        self.small_bin = small_bin
        self.normal_bin = normal_bin
        self.multiplier = multiplier
        self.legend = legend

# =============================================================================
# Main
# =============================================================================
cwd = os.getcwd()
plot_dir = "plots/0.v_profile_insert"


if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)    


# =============================================================================
# Simulation type definitions
# =============================================================================


# pressure driven simulations with p_dist
pdist = SimulationInsert("pdist", "6.1.Applying_force_p_dist_0.00063_bin_0.1/1",
                   "6.Applying_force_p_dist_0.00063/1", -1,'upper left' ) 

# Gravitation-like force
pgrav = SimulationInsert("pgrav", "5.1.Applying_force_p_0.00063_bin_0.1/1",
                   "5.Applying_force_p_0.00063/1", 1, 'upper left') 
# Difusio-osmosis
do = SimulationInsert("do", "4.1.Applying_force_0.125_bin_0.1/1",
                   "4.Applying_force_0.125/1", 1,'upper right') 



# =============================================================================
# Main
# =============================================================================
# Define the type of simulation
sim = pdist.copy
sim.print_params(logger)

logger.info("The data for the main data is from %s" %sim.normal_bin)
logger.info("The data for the insert data is from %s" %sim.small_bin)

# Loading the data for the bin 0.1 \sigma
fluid_n = da.DensityDistribution("properties_short.dat", "rBulk", directory = sim.normal_bin) 
solute_n = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = sim.normal_bin) 
solvent_n = da.DensityDistribution("Fproperties_short.dat", "rBulk", directory = sim.normal_bin) 

# Loading the data for the bin 0.25 \sigma
fluid_s = da.DensityDistribution("properties_short.dat", "rBulk", directory = sim.small_bin) 
solute_s = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = sim.small_bin) 
solvent_s = da.DensityDistribution("Fproperties_short.dat", "rBulk", directory = sim.small_bin) 


# Applying the multiplier
solute_n.data_frame['vx'] = sim.multiplier * solute_n.data_frame['vx']
solvent_n.data_frame['vx'] = sim.multiplier * solvent_n.data_frame['vx']
fluid_n.data_frame['vx'] = sim.multiplier * fluid_n.data_frame['vx']
solute_s.data_frame['vx'] = sim.multiplier * solute_s.data_frame['vx']
solvent_s.data_frame['vx'] = sim.multiplier * solvent_s.data_frame['vx']
fluid_s.data_frame['vx'] = sim.multiplier * fluid_s.data_frame['vx']


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
ax1.legend(["Solute", "Solvent", "Fluid"], loc = sim.legend)

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
fig1.savefig('%s/vprofile_insert_%s.pdf'%(plot_dir,sim.name))


