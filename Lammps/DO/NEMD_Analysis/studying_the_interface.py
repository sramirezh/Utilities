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


fluid = da.PropertyDistribution("properties_short.dat") 
solute = da.PropertyDistribution("Sproperties_short.dat") 
solvent = da.PropertyDistribution("Fproperties_short.dat") 

# =============================================================================
# Plotting
# =============================================================================
cf.set_plot_appearance()

plt.close("all")
fig, ax = plt.subplots()

ax = fluid.plot_property_dist("vx", ax = ax)

#fluid.data_frame.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'line', legend = False)
#solvent.data_frame.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'line', legend = False)
#solute.data_frame.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'line', legend = False)
#ax.axvline(x = lz_min_half, ls=':',c='black')
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$v_x(z)$')
fig.tight_layout()
ax.legend(["Fluid", "Solvent", "Solute"])
fig.savefig('vprofile.pdf')


ymin, ymax = plt.ylim()

# Detail at the interface
ax.set_xlim(0,2)
ax.set_ylim(ymin,0.005)
fig.tight_layout()
fig.savefig('vprofile_zoom.pdf')


# Detail at the interface
ax.set_xlim(0,5)
ax.set_ylim(ymin,0.03)
fig.tight_layout()
fig.savefig('vprofile_zoom2.pdf')