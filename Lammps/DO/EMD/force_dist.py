#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Scripts that computes the force distribution and generates the forces to 
be applied in LAMMMPS

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
import Lammps.General.Thermo_Analyser as ta
import Lammps.lammps_utilities as lu
import density_analysis as da


def force_dist(ns, nf, solvent, solute, fluid):
    """
    Force calculation
    :param Ns: Solute Number in the Bulk
    :param Nf: Solvent Number in the Bulk
    :param SProperties: Solute Properties
    :param FProperties: Solvent Properties
    :param AProperties: Fluid Properties
    :return:
    The Force distribution 

    """

    fs = 1
    ff = -fs * ns / nf
    n, m = fluid.data_frame.shape
    force = np.zeros((n, 2))
    force[:, 0] = fluid.data_frame['Coord1']
    
    rho_fluid = fluid.data_frame['density/mass'].values
    i_n_zero = np.where( rho_fluid != 0 )
    
    f_dist_s = fs * solute.data_frame['density/mass'].values
    f_dist_f = ff * solvent.data_frame['density/mass'].values
    
    
    force[i_n_zero, 1] = (f_dist_s[i_n_zero] + f_dist_f[i_n_zero]) / rho_fluid[i_n_zero]

    np.savetxt("Force.dat", force)
    print("\n*********************Getting the force profile*********************\n")
    print("The force on the Solutes is %f, on the Solvents %f" %(fs,ff))
    print("Created the force distribution File Force.dat ")

    return force

def plot_force_dist(force, z_min, z_max):
    cf.set_plot_appearance()
    
    x_tick_distance = 5
    # Force distribution
    fig1,ax1=plt.subplots()
    ax1.plot(force[:,0],force[:,1])
    ax1.set_ylabel(r"$F$")
    ax1.set_xlabel(r'$d[\sigma]$')
    ax1.set_xlim(z_min, z_max)
    ax1.set_xticks(np.arange(z_min, z_max, x_tick_distance))
    ax1.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
    ax1.axvline(x=z_pos[-1], ymin=0, ymax=1,ls='-.',c='black')
    ax1.axvline(x=10, ymin=0, ymax=1,ls='-.',c='b')
    ax1.axvline(x=20, ymin=0, ymax=1,ls='-.',c='b')
    plt.tight_layout()
    fig1.savefig("Force_dist.pdf")

# =============================================================================
# Main
# =============================================================================

sim = lu.Simulation("log.lammps")

solute = da.DensityDistribution("Sproperties_short.dat")
solvent = da.DensityDistribution("Fproperties_short.dat")
fluid = da.DensityDistribution("properties_short.dat")

ns = sim.thermo_ave.at['v_cBSolu','Average']
nf = sim.thermo_ave.at['v_cBSolv','Average']

# =============================================================================
#  Getting the concentration in the bulk, using the results from log. lammps
# =============================================================================

h_b, limits_b = lu.read_region_height("rBulk") # Getting bulk limits

if sim.is_2d:
    print("\nThis is a 2d Simulation\n")
    
    #Added the 0.5 because of the 2D correction
    v_b = (sim.volume * (h_b / (sim.limits[1][1] - sim.limits[0][1]))) / 0.5 
    index  = 1
    
else:
    print("\nThis is a 3d Simulation\n")
    v_b = (sim.volume * (h_b / (sim.limits[1][2] - sim.limits[0][2]))) 
    index = 2

# Concentration of species in the bulk
cs_b = ns / v_b
cf_b = nf / v_b


force = force_dist(ns, nf, solvent, solute, fluid)


# =============================================================================
# Finding the cutoff for the force
# =============================================================================

# To find the plateau after the maximum, as before there might be a plateau
imax =np.argmax(force[:,1])

cut_off = imax + min(cf.plateau_finder(force[imax::,1]))[-1]

zmax_force = fluid.data_frame.values[cut_off,1]


# =============================================================================
# Creating the force distribution
# =============================================================================
z_pos = force[:cut_off+1, 0]
f_mu = np.transpose(force[:cut_off, 1])
print("The Force Cut-off is %f, this is where the region of applied forces finishes"%np.max(z_pos))
print("Creating the Files to iterate in Lammps")
print ("In the bulk, the solute concentration is %s, the solvent is %s"% (cs_b, cf_b))
np.savetxt("Zpos_iterate.dat", z_pos)
np.savetxt("Force_iterate.dat", f_mu)



# =============================================================================
# Plots
# =============================================================================

#Plots without shifting
z_min = 0
z_max = sim.limits[1][index]
x_tick_distance = 5


plot_force_dist(force, z_min, z_max)


