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





def plot_zoom(ax,fig, xlim, ylim, name):
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    fig.tight_layout()
    fig.savefig("%s.pdf"%name)
    
def get_ymax(xmin, xmax, fluid, solute, solvent, prop_name):
    """
    Returns the closest x to xmax (from below) and the bounding values in 
    y for that range
    
    Assumes that the property grows
    """
    
    # Assuming that all objects have the same x coordinates
    ind_x_max = fluid.data_frame['Coord1'].iloc[(fluid.data_frame['Coord1']- xmax ).abs().argsort()[:2]].idxmin()
    xmax = fluid.data_frame['Coord1'].iloc[ind_x_max]
    
    ind_x_min= fluid.data_frame['Coord1'].iloc[(fluid.data_frame['Coord1']- xmin ).abs().argsort()[:2]].idxmin()
    xmin = fluid.data_frame['Coord1'].iloc[ind_x_min]
    
    # getting the maximum values in that range
    ymax = max(fluid.data_frame[prop_name].iloc[ind_x_min:ind_x_max].max(), 
               solute.data_frame[prop_name].iloc[ind_x_min:ind_x_max].max(), 
               solvent.data_frame[prop_name].iloc[ind_x_min:ind_x_max].max())
    
    ymin = min(fluid.data_frame[prop_name].iloc[ind_x_min:ind_x_max].min(), 
           solute.data_frame[prop_name].iloc[ind_x_min:ind_x_max].min(), 
           solvent.data_frame[prop_name].iloc[ind_x_min:ind_x_max].min())

    return xmin, xmax, ymin, ymax
 

logger = cf.log(__file__, os.getcwd())    

fluid = da.DensityDistribution("properties_short.dat", "rBulk") 
solute = da.DensityDistribution("Sproperties_short.dat", "rBulk") 
solvent = da.DensityDistribution("Fproperties_short.dat", "rBulk") 

# =============================================================================
# Plotting the velocities
# =============================================================================
cf.set_plot_appearance()

plt.close("all")
fig, ax = plt.subplots()

fluid.plot_property_dist("vx", ax = ax)
solvent.plot_property_dist("vx", ax = ax)
solute.plot_property_dist("vx", ax = ax)
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')


ax.set_xlim(0, fluid.positions[-1])
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$v_x(z)$')

fig.tight_layout()
ax.legend(["Fluid", "Solvent", "Solute"])
fig.savefig('vprofile.pdf')

ymin = fluid.data_frame['vx'].min()


# Getting the index closest to the limit
xmin = 0
xmax = 2
xmin, xmax, ymin, ymax = get_ymax(xmin, xmax, fluid, solute, solvent, "vx" )
plot_zoom(ax, fig,  [0,xmax], [ymin, ymax], "vprofile_zoom")


xmin = 0
xmax = 5
xmin, xmax, ymin,ymax = get_ymax(xmin, xmax, fluid, solute, solvent, "vx" )
plot_zoom(ax, fig,  [0, xmax], [ymin, ymax], "vprofile_zoom2")


# =============================================================================
# Plotting the densities
# =============================================================================

fig, ax = plt.subplots()

fluid.plot_property_dist("density/mass", ax = ax)
solvent.plot_property_dist("density/mass", ax = ax)
solute.plot_property_dist("density/mass", ax = ax)
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
ax.axvspan(fluid.limits_b[0], fluid.limits_b[1], alpha=0.5, color='green')


ax.set_xlim(0, fluid.positions[-1])
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$c(z)$')

fig.tight_layout()
ax.legend(["Fluid", "Solvent", "Solute"], loc = 'upper right')
fig.savefig('rhoprofile.pdf')

ymin = fluid.data_frame['vx'].min()


# Getting the index closest to the limit
xmin = 0
xmax = 2
xmin, xmax, ymin,ymax = get_ymax(xmin, xmax, fluid, solute, solvent, "density/mass" )
plot_zoom(ax, fig,  [0,xmax], [ymin, ymax], "rhoprofile_zoom")

xmin = 0
xmax = 5
xmin, xmax, ymin,ymax = get_ymax(xmin,xmax, fluid, solute, solvent, "density/mass" )
plot_zoom(ax, fig,  [0, xmax], [ymin, ymax], "rhoprofile_zoom2")


