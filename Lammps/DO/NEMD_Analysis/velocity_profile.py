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
import Lammps.General.chunk_utilities as cu
from io import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import matplotlib.pyplot as plt



chunks = cf.read_data_file("properties_short.dat") 

# =============================================================================
# Plotting
# =============================================================================
cf.set_plot_appearance()

plt.close("all")
fig, ax = plt.subplots()
chunks.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'line', legend = False,)
#ax.axvline(x = lz_min_half, ls=':',c='black')
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$v_x(z)$')
fig.tight_layout()
fig.savefig('vprofile.pdf')

