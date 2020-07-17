#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Scripts that computes the force distribution and generates the forces to 
be applied in LAMMMPS


The quantities for solute, solvent and fluid where measured in chunks with
"bin/1d z lower 0.1"

TODO this code assumes the surface at z = 0

@author: sr802
"""

import numpy as np
import sys
import os
from shlex import split
import matplotlib.pyplot as plt
import linecache
import re
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.General.thermo_analyser as ta
import Lammps.lammps_utilities as lu
import density_analysis as da


def force_dist_mu(solvent, solute, fluid):
    """
    Average force  per chunks due to a chemical potential gradient 
    Args:
        
    solute: Instance of the DensityDistribution with all the solute properties
    solvent: Instance of the DensityDistribution with all the solvent properties
    solution: Instance of the DensityDistribution with all the solution properties
    
    Returns
    The Force distribution as a file Force.dat

    """
    
    cs_b = solute.rho_bulk
    cf_b = solvent.rho_bulk
    fs = 1
    ff = -fs * cs_b / cf_b
    n, m = fluid.data_frame.shape
    force = np.zeros((n, 2))
    force[:, 0] = fluid.data_frame['Coord1']
    
    rho_fluid = fluid.data_frame['density/mass'].values
    
    # The force has rho_fluid in the denominator, so we need to avoid the zeros
    i_n_zero = np.where( rho_fluid != 0 )
    
    f_dist_s = fs * solute.data_frame['density/mass'].values
    f_dist_f = ff * solvent.data_frame['density/mass'].values
    
    
    force[i_n_zero, 1] = (f_dist_s[i_n_zero] + f_dist_f[i_n_zero]) / rho_fluid[i_n_zero]

    np.savetxt("Force.dat", force)
    print("\n*********************Getting the force profile*********************\n")
    print("The force on the Solutes is %f, on the Solvents %f" %(fs,ff))
    print("Created the force distribution File Force.dat ")

    return force


def force_dist_p(fluid):
    """
    Average force  per chunks due to a pressure gradient, that is assumed to be
    1
    Args:
        
    solute: Instance of the DensityDistribution with all the solute properties
    solvent: Instance of the DensityDistribution with all the solvent properties
    solution: Instance of the DensityDistribution with all the solution properties
    
    Returns
    The Force distribution as a file Force.dat

    """

    grad_p = 1
    rho_fluid = fluid.data_frame['density/mass'].values
    
     # The force has rho_fluid in the denominator, so we need to avoid the zeros
    i_n_zero = np.where( rho_fluid != 0 )
    n, m = fluid.data_frame.shape
    force = np.zeros((n, 2))
    
    force[:, 0] = fluid.data_frame['Coord1']
    
    force[i_n_zero, 1] = -grad_p / rho_fluid[i_n_zero]
        

    np.savetxt("Force_p.dat", force)
    print("\n*********************Getting the force profile for grad_p*****\n")
    print("Created the force distribution File Force_p.dat ")

    return force

def plot_force_dist(force, sim, fluid, ax):
    """
    Plots the force distribution
    showing the bulk and the limit where the chunks for the force will finish
    """
    z_min = sim.limits[0][sim.index]
    z_max = sim.limits[1][sim.index]
    x_tick_distance = 5
    # Force distribution

    ax1.plot(force[:,0],force[:,1])
    ax1.axvspan(fluid.limits_b[0], fluid.limits_b[1], alpha=0.5, color='green')
    ax1.set_xlim(z_min, z_max)
    ax1.set_xticks(np.arange(z_min, z_max, x_tick_distance))
    ax1.axhline(y=0, xmin=0, xmax=1,ls='-.',c='black')




# =============================================================================
# Main
# =============================================================================

cwd = os.getcwd()
plot_dir = "plots/0.Force_distributions"


if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   


sim = lu.Simulation("log.lammps")
sim.plot_dir = plot_dir


solute = da.DensityDistribution("Sproperties_short.dat", [10,20])
solvent = da.DensityDistribution("Fproperties_short.dat", [10,20])
fluid = da.DensityDistribution("properties_short.dat", [10,20])

#ns = sim.thermo_ave.at['v_cBSolu','Average']
#nf = sim.thermo_ave.at['v_cBSolv','Average']

# =============================================================================
#  Getting the concentration in the bulk, using the results from log.lammps
# =============================================================================


if sim.is_2d == True:
    logger.info("\nThis is a 2d Simulation\n")
    #Added the 0.5 because of the 2D correction
    v_b = (sim.volume * (fluid.h_b / (sim.limits[1][1] - sim.limits[0][1]))) / 0.5 
    index  = 1
    
else:
    logger.info("\nThis is a 3d Simulation\n")
    v_b = (sim.volume * (fluid.h_b/ (sim.limits[1][2] - sim.limits[0][2]))) 
    index = 2
    
# Tells which dimension is perpendicular to the wall
sim.index = index
    
# =============================================================================
# Creating the force distribution for the chemical potential gradient
# =============================================================================

force_mu = force_dist_mu(solvent, solute, fluid)
# =============================================================================
# Finding the cutoff for the force
# =============================================================================
# To find the plateau after the maximum, as before there might be a plateau
imax = np.argmax(force_mu[:,1])

cut_off = imax + min(cf.plateau_finder(force_mu[imax::,1]))[-1]

zmax_force = fluid.data_frame['Coord1'][cut_off]

# The original bins we indicated with the lower limit of the bin
z_pos = force_mu[:cut_off+1, 0]
# This is the force applied to all the particles inside the bin
f_mu = np.transpose(force_mu[:cut_off, 1])

logger.info("The Force Cut-off is %f, this is where the region of applied forces finishes"%np.max(z_pos))
logger.info("Creating the Files to iterate in Lammps")
logger.info("In the bulk, the solute concentration is %s, the solvent is %s"% (solute.rho_bulk, solvent.rho_bulk))
np.savetxt("Zpos_iterate.dat", z_pos)
np.savetxt("Force_iterate.dat", f_mu)



# =============================================================================
# Creating the force distribution for the Pressure gradient
# =============================================================================
force_p =  force_dist_p(fluid)

# Finding the nearest point after the begining of the bulk
imax = min(np.where(force_p[:,0]>=10)[0])

force_p_red = force_p[:imax, :]

z_pos_p = np.append(force_p_red[:,0], sim.limits[1,2])
f_p = force_p_red[:,1]

logger.info("The Force Cut-off is %f, this is where the region of applied forces finishes"%np.max(z_pos))
logger.info("Creating the Files to iterate in Lammps")
logger.info("In the bulk, the solute concentration is %s, the solvent is %s"% (solute.rho_bulk, solvent.rho_bulk))
np.savetxt("Zpos_iterate_p.dat", z_pos_p)
np.savetxt("Force_iterate_p.dat", f_p)

# Modifying the force distribution array to have one last chunk after the 

# =============================================================================
# Plots
#Plots without shifting
# =============================================================================
cf.set_plot_appearance()

# Force due to the chemical potential
fig1, ax1=plt.subplots()
plot_force_dist(force_mu, sim, fluid, ax1)
ax1.set_ylabel(r"$F^{\mu}_{\text{ave}}$")
ax1.set_xlabel(r'$z[\sigma]$')
ax1.set_xlim(0, 30)
ax1.axvline(x = z_pos[-1], ymin=0, ymax=1,ls=':',c='black')
fig1.tight_layout()
fig1.savefig("%s/Force_mu_dist.pdf"%plot_dir)


# Force due to the pressure gradient
fig1 ,ax1 = plt.subplots()

# Because I applied the wrong sign of the gradient
force_p[:,1] = -force_p[:,1]

plot_force_dist(force_p, sim, fluid, ax1)
ax1.set_ylabel(r"$F^{P}$")
ax1.set_xlabel(r'$z[\sigma]$')
ax1.set_yscale('symlog')
#ax1.set_ylim(-1, 10)
ax1.set_xlim(0, 30)
ax1.axhline(y=1/fluid.rho_bulk, xmin=0, xmax=1,ls='-.',c='black')
ax1.axvline(x = z_pos_p[-2], ymin=0, ymax=1,ls=':',c='black')
fig1.tight_layout()
fig1.savefig("%s/Force_p_dist.pdf"%plot_dir)



