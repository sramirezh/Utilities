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
import chunk_reader as cr
from io import StringIO

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
results = cr.chunk_reader("velocity_all_steps.dat") 


# Creating all the frames

f = open( results.filename)

series = []

byte_pos = results.offsets

for i, start in enumerate(byte_pos[:-1]):

    f.seek(start)
    step, n_chunks, _ = f.readline().split() # Timestep Number-of-chunks Total-count
    start = f.tell()
    stuff = f.read(byte_pos[i+1] - start)
    data = pd.read_csv(StringIO(stuff),sep=" ",header=None).dropna(axis=1,how='all')
    data.columns = results.header
    series.append(cr.timestep(step, n_chunks, data))
    
    data['vx'] get this one and append






#chunks = cf.read_data_file()

data = chunks.values



results = chunk_reader("velocity_all.dat") 
        
    




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

