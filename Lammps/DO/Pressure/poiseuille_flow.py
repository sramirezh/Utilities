# -*- coding: utf-8 -*-
"""
Spyder Editor
Script to analyse lammps chunk files from pressure driven
Poiseuille_flow/Todd1995_NoseHoover

"""

import numpy as np
import pandas as pd
import argparse
import os
import sys


sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import matplotlib.pyplot as plt



# =============================================================================
#  Problem parameters and simple calculations
# =============================================================================
lz_wall = 1.5625
lz_fluid = 25
n_atoms = 563
lx = 7.4217
fx = 0.01


rho_fluid = n_atoms/(lx**2*lz_fluid)
lz_min = -(lz_wall+lz_fluid/2)
lz_min_half = -(lz_wall+lz_fluid)/2



# =============================================================================
# Reading data
# =============================================================================

chunks = cf.read_data_file("velocity_all.dat")

data = chunks.values


# =============================================================================
# Plotting
# =============================================================================

plt.close("all")
fig, ax = plt.subplots()
chunks.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'scatter', legend = False,)
ax.axvline(x = lz_min_half, ls=':',c='black')
plt.show()




# =============================================================================
# Data manipulation
# =============================================================================
#Now I need to move the data before the vertical line to the end, working with numpy

indexes = np.where(data[:,1]<= lz_min_half)[0]

data[indexes,1] = (data[indexes,1]-lz_min)+lz_fluid/2
data =  np.roll(data, -len(indexes), axis = 0)

chunks_new = pd.DataFrame(data[:,1:], columns = chunks.columns[1:])


fig, ax = plt.subplots()
chunks_new.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'scatter', legend = False,)
plt.show()



# =============================================================================
#  Fitting
# =============================================================================
from scipy import optimize


def model(p,x):
    return p[0] + p[1] * x**2

def residual(p,x,y):
    return y - model(p,x)
    

pinit = [1.0, -1.0]
out = optimize.leastsq(residual, pinit,
                       args=(chunks_new['Coord1'],chunks_new['vx']), full_output=1)
  
pfinal = out[0]

viscosity = -((rho_fluid*fx)/(2*pfinal[1]))

vx_predicted = model(pfinal,chunks_new['Coord1']) 


print ("The viscosity is %s for a fluid with average density %s"%(viscosity, rho_fluid))


fig, ax = plt.subplots()
chunks_new.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'scatter', legend = False)
ax.plot(chunks_new['Coord1'],vx_predicted, c = 'k')
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$v_x(z)$')
ax.legend(["NEMD","NS"])


plt.show()

