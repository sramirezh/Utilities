#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 17:27:17 2020

Script to animate the velocity profile
https://towardsdatascience.com/animations-with-matplotlib-d96375c5442c

@author: simon
"""

import numpy as np
import pandas as pd
import argparse
import os
import sys
from io import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.General.chunk_utilities as cu
import Lammps.core_functions as cf
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation



# =============================================================================
# Reading data
# =============================================================================
results = cu.chunk_reader("properties.dat") 

cf.set_plot_appearance()

fig, ax = plt.subplots()


def init():
    line.set_data([], [])
    return line,


fig = plt.figure()
ax = plt.axes(xlim=(0, 4), ylim=(-2, 2))
line, = ax.plot([], [], lw=3)


def animate(i, results):
    x = results[i].data['Coord1'].values
    y = results[i].data['vx'].values
#    x = res[i].data['vx'].values
#    y = res[i].data['vx'].values
    line.set_data(x, y)
    return line,



indexes = np.arange(0,len(results),1000)

anim = FuncAnimation(fig, animate, init_func=init,
                               frames = len(indexes), interval=20, blit=True, fargs =[results])


anim.save('sine_wave.gif', writer='imagemagick')


#
#chunks.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'line', legend = False,)
##ax.axvline(x = lz_min_half, ls=':',c='black')
#ax.set_xlabel(r'$z[\sigma] $')
#ax.set_ylabel(r'$v_x(z)$')
#fig.tight_layout()
#fig.savefig('vprofile.pdf')

