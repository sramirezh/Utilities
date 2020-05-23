#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 23 11:53:49 2020
Script to focus at the interface, specially, to analyse the simulatios of DO
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
    
def get_ymax(xmax, fluid, solute, solvent, prop_name):
    """
    Returns the closest x to xmax (from below) and the largest value of the property 
    at that point.
    
    Assumes that the property grows
    """
    
    ind_x = fluid.data_frame['Coord1'].iloc[(fluid.data_frame['Coord1']- xmax ).abs().argsort()[:2]].idxmin()
    xmax = fluid.data_frame['Coord1'].iloc[ind_x]
    ymax = max(fluid.data_frame[prop_name].iloc[ind_x], solute.data_frame[prop_name].iloc[ind_x], solvent.data_frame[prop_name].iloc[ind_x] )
    
    return xmax, ymax
 
fluid = da.PropertyDistribution("properties_short.dat") 
solute = da.PropertyDistribution("Sproperties_short.dat") 
solvent = da.PropertyDistribution("Fproperties_short.dat") 

# =============================================================================
# Plotting
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
xmax = 2
xmax, ymax = get_ymax(xmax, fluid, solute, solvent, "vx" )
plot_zoom(ax, fig,  [0,xmax], [ymin, ymax], "vprofile_zoom")


xmax = 5
xmax, ymax = get_ymax(xmax, fluid, solute, solvent, "vx" )
plot_zoom(ax, fig,  [0, xmax], [ymin, ymax], "vprofile_zoom2")
