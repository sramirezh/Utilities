# -*- coding: utf-8 -*-
"""
Spyder Editor


Script to analyse lammps chunk files 
At this point it is based on P
"""

import numpy as np
import pandas as pd
import argparse
import os
import sys


sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import matplotlib.pyplot as plt


cf.set_plot_appearance
chunks = cf.read_data_file("velocity.dat")

plt.close("all")

# Plotting
fig, ax = plt.subplots()
chunks.plot(x = 'Coord1',y = 'vx', ax = ax)

plt.show()

#class Lammps_chunk(object):
#    def