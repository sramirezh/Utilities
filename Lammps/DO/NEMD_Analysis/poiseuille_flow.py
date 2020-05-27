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

from io import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.General.chunk_utilities as cu
import Lammps.core_functions as cf
import matplotlib.pyplot as plt
import Others.Statistics.FastAverager as stat
from scipy import optimize


def err_func(p, x, y, err):
    return (y - model(p, x)) / err

def model(p,x):
    return p[0] + p[1] * x**2

def residual(p,x,y):
    return y - model(p,x)

class SimulationType(object):
    """Class to define the type of simulation analysed, for example N2 or
    Octane
    """
    def __init__(self, name, lz_wall, h_fluid, lx, fx, n_atoms):
        """
        Args:
            name: to identify the simulation type, eg. N2
            lz_wall: wall heigth 
            h_fluid: fluid region heigth
            lx: box side in x and y
            fx: force applied on the particles
            n_atoms: number of fluid atoms
        """
        self.name = name
        self.lz_wall= lz_wall  # Time in femptoseconds
        self.h_fluid = h_fluid   # sampling_interval myDump
        self.lx = lx  
        self.fx = fx
        self.n_atoms = n_atoms
        
        
    def print_params(self):
        print("Using the parameters from %s"%self.name)
 


# =============================================================================
# main
# =============================================================================
logger = cf.log(__file__, os.getcwd())    


# =============================================================================
#  Problem parameters and simple calculations
# =============================================================================

Todd1997 = SimulationType("Todd1997", 37.0193942870738-34.5, 69, 4.6840, 0.005, 1278 )


sim = Todd1997

sim.print_params()


lz_fluid = sim.h_fluid - 1 # available for the fluid excluding the volume See Todd1995, Todd1997 and Todd2017)                        
rho_fluid = sim.n_atoms/(sim.lx**2*lz_fluid)
lz_min = -(sim.lz_wall+lz_fluid/2)
lz_min_half = -(sim.lz_wall+lz_fluid)/2



# =============================================================================
# Reading data
# =============================================================================
results = cu.chunk_reader("velocity_all_steps.dat") 
# Creating all the frames

f = open( results.filename)

series = []
velocities = []
byte_pos = results.offsets

for i, start in enumerate(byte_pos[:-1]):

    f.seek(start)
    step, n_chunks, _ = f.readline().split() # Timestep Number-of-chunks Total-count
    start = f.tell()
    stuff = f.read(byte_pos[i+1] - start)
    data = pd.read_csv(StringIO(stuff),sep=" ",header=None).dropna(axis=1,how='all')
    data.columns = results.header
    series.append(cu.timestep(step, n_chunks, data))
    velocities.append(data['vx'].values)
    # This could be done inside the time step maybe

from scipy.stats import sem 
velocities = np.array(velocities)
ave_velocities = np.average(velocities, axis = 0 )
error_velocities = sem(velocities, axis = 0)



stat_data = []
error_correlation =[]
for i in range(int(n_chunks)):
    results = stat.fast_averager(velocities[:,i])[0]
    error_correlation.append(results[2])
    
    
# Creating a data frame with the average velocity and error

data['vx'] = ave_velocities
data['vx_error'] = error_correlation

chunks = data

data = chunks.values

# =============================================================================
# Plotting
# =============================================================================

cf.set_plot_appearance()


plt.close("all")
fig, ax = plt.subplots()
#chunks.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'scatter', legend = False,)
ax.plot(data[:,1], ave_velocities)
ax.fill_between(data[:,1], ave_velocities-error_correlation, ave_velocities+error_correlation , alpha=0.4)
ax.axvline(x = lz_min_half, ls=':',c='black')
ax.axvline(x = -0.5*lz_fluid, ls='-.',c='black')
ax.axvline(x = 0.5*lz_fluid, ls='-.',c='black')
ax.axvline(x = 0, ls=':',c='black')
fig.savefig("raw_data.pdf")

# =============================================================================
# Data manipulation OLD
# =============================================================================
# Only taking into account the data inside the fluid region

indexes = np.where((data[:,1]<= 0.5*lz_fluid)&(data[:,1]>= -0.5*lz_fluid))[0]
data = data[indexes, :]

#data = data[:-5,:]
# New df with only the data from fluid region
chunks_new = pd.DataFrame(data[:,1:], columns = chunks.columns[1:])

# Ploting the interes data
fig, ax = plt.subplots()
chunks_new.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'scatter', legend = False,)
ax.axvline(x = 0, ls=':',c='black')
#plt.show()


# =============================================================================
#  Fitting
# =============================================================================

pinit = [1.0, -1.0]
out = optimize.leastsq(err_func, pinit,
                       args=(chunks_new['Coord1'], chunks_new['vx'],chunks_new['vx_error'] ), full_output=1)
  
pfinal = out[0]

viscosity = -((rho_fluid*sim.fx)/(2*pfinal[1]))

vx_predicted = model(pfinal, chunks_new['Coord1']) 


logger.info("The viscosity is %s for a fluid with average density %s"%(viscosity, rho_fluid))


ave_velocities = chunks_new['vx']
error_correlation = chunks_new['vx_error']

fig, ax = plt.subplots()
chunks_new.plot(x = 'Coord1',y = 'vx', ax = ax, kind = 'line', legend = False)
ax.fill_between(chunks_new['Coord1'], ave_velocities - error_correlation, ave_velocities + error_correlation , alpha = 0.4)
ax.plot(chunks_new['Coord1'],vx_predicted, c = 'k')
#ax.axvline(x = 0, ls='--', c='black')
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$v_x(z)$')
ax.set_ylim(0,None)
ax.legend(["NEMD","NS"])
fig.tight_layout()
#plt.show()
fig.savefig("poiseuille.pdf")



