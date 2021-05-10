#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 23 11:53:49 2020
Script to focus on the profiles, specially, to analyse the simulatios of DO and
pressure driven
with small binning close to the surface
@author: simon
"""


import numpy as np
import pandas as pd
import argparse
import os
import sys
from io import StringIO
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.General.chunk_utilities as cu
import Lammps.core_functions as cf
import Lammps.DO.EMD.density_analysis as da


# =============================================================================
# Main
cwd = os.getcwd()
plot_dir = "plots/0.studying_profiles"


if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)      

fluid = da.DensityDistribution("properties_short.dat", "rBulk") 
solute = da.DensityDistribution("Sproperties_short.dat", "rBulk") 
solvent = da.DensityDistribution("Fproperties_short.dat", "rBulk") 

# =============================================================================
# Plotting the velocities
# =============================================================================
cf.set_plot_appearance()

plt.close("all")
fig1, ax1 = plt.subplots()


solute.plot_property_dist("vx", ax = ax1)
solvent.plot_property_dist("vx", ax = ax1)
fluid.plot_property_dist("vx", ax = ax1)


#ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')


ax1.set_xlim(0, fluid.positions[-1])
ax1.set_ylim(0, None)
ax1.set_xlabel(r'$z[\sigma] $')
ax1.set_ylabel(r'$v_x(z)$')

# Getting the index closest to the limit
xmin = 0
xmax = 30
ax1 = cf.plot_zoom(ax1, [xmin, xmax])

fig1.tight_layout()
ax1.legend(["Solute", "Solvent", "Fluid" ])
fig1.savefig('%s/vprofile.pdf'%plot_dir)

ymin = fluid.data_frame['vx'].min()


# Getting the index closest to the limit
xmin = 0
xmax = 2
ax1 = cf.plot_zoom(ax1, [xmin, xmax])
fig1.tight_layout()
fig1.savefig('%s/vprofile_zoom.pdf'%plot_dir)



xmin = 0
xmax = 5
ax1 = cf.plot_zoom(ax1, [xmin, xmax])
fig1.tight_layout()
fig1.savefig('%s/vprofile_zoom2.pdf'%plot_dir)

## Adding the insert
#
#left, bottom, width, height = [0.55, 0.25, 0.4, 0.30]
#ax2 = fig1.add_axes([left, bottom, width, height])
#ax2.set_ylabel(r'$V(r)$',fontsize =17, labelpad=-5)
#ax2.set_xlabel(r'$r$' ,fontsize =17, labelpad=-5)
#
#
#fig1.tight_layout()
#fig1.savefig('vprofile_insert.pdf')

# =============================================================================
# Plotting the densities
# =============================================================================

fig2, ax2 = plt.subplots()

solute.plot_property_dist("density/mass", ax = ax2)
solvent.plot_property_dist("density/mass", ax = ax2)
fluid.plot_property_dist("density/mass", ax = ax2)


#ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
#ax2.axvspan(fluid.limits_b[0], fluid.limits_b[1], alpha=0.5, color='green')


ax2.set_xlim(0, fluid.positions[-1])
ax2.set_ylim(0, None)
ax2.set_xlabel(r'$z[\sigma] $')
ax2.set_ylabel(r'$c(z)$')

fig2.tight_layout()
ax2.legend(["Solute", "Solvent", "Fluid" ], loc = 'upper right')
# Getting the index closest to the limit
xmin = 0
xmax = 30
ax2 = cf.plot_zoom(ax2, [xmin, xmax])
fig2.tight_layout()
fig2.savefig('%s/rhoprofile.pdf'%plot_dir)

ymin = fluid.data_frame['vx'].min()


# Getting the index closest to the limit
xmin = 0
xmax = 2
ax2 = cf.plot_zoom(ax2, [xmin, xmax])
fig2.tight_layout()
fig2.savefig('%s/rhoprofile_zoom.pdf'%plot_dir)


xmin = 0
xmax = 5

ax2 = cf.plot_zoom(ax2, [xmin, xmax])
fig2.tight_layout()
fig2.savefig('%s/rhoprofile_zoom2.pdf'%plot_dir)



