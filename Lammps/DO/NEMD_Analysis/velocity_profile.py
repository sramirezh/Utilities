#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 13:40:22 2020
Script to plot the velocity profiles
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
import density_analysis as da


fluid = da.PropertyDistribution("properties_short.dat") 
solute = da.PropertyDistribution("Sproperties_short.dat") 
solvent = da.PropertyDistribution("Fproperties_short.dat") 

# =============================================================================
# Plotting
# =============================================================================
cf.set_plot_appearance()

plt.close("all")
fig, ax = plt.subplots()
fluid.data_frame.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'line', legend = False)
solvent.data_frame.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'line', legend = False)
solute.data_frame.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'line', legend = False)
#ax.axvline(x = lz_min_half, ls=':',c='black')
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$v_x(z)$')
fig.tight_layout()
ax.legend(["Fluid", "Solvent", "Solute"])
fig.savefig('vprofile.pdf')


# Detail at the interface
ax.set_xlim(0,2)
ax.set_ylim(0,0.005)
fig.savefig('vprofile_zoom.pdf')