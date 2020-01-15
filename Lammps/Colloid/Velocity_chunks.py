#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 11:31:23 2020
Script to analyse the chunks
@author: sr802
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path


from Lammps.PDP.Plots.LJ_plotter import LJ
import Lammps.core_functions as cf


from ovito.modifiers import *
import ovito.io as ov


header = cf.read_data_file("Sproperties_t.dat").columns.values


solutes_t = cf.read_data_file("Sproperties_t.dat").values
solvents_t = cf.read_data_file("Fproperties_t.dat").values

solutes_b = cf.read_data_file("Sproperties_b.dat").values
solvents_b = cf.read_data_file("Fproperties_b.dat").values

sol = cf.read_data_file("properties_b.dat").values

solutes = 0.5*(solutes_t+solutes_b)
solvents = 0.5*(solvents_t+solvents_b)

vx_total = sol[:,5]
# determining the space between the reservoirs

node = ov.import_file("equil.dat", multiple_frames = False)
data = node.compute()
list(data.particle_properties.keys())


# =============================================================================
# Some calculations
# =============================================================================

#Getting the properties of the box from the equilibrated configuration
#Assuming that the box does not change size
box = node.compute(0).cell # getting the properties of the box
L = np.diag(box.matrix) 
center = L/2.0

lattice_constant = (4/0.8)**(1/3)

indexes = np.where((solutes[:,1]> 6*lattice_constant) & (solutes[:,1]<18*lattice_constant))[0]

mu_high_c = -1.7478856237659908 
mu_low_c = -3.1341799848858813


grad_mu = (-1.7478856237659908 -(-3.1341799848858813))/(12*lattice_constant)

print ("The average velocity of the solution in the bulk is %s" %np.average(vx_total[indexes]))
print ("The chemical potential gradient is %s" %grad_mu )
fig,ax = plt.subplots()

r_colloid  = 3.23

ax.plot(solutes[indexes,1],solutes[indexes,5], label='Solutes')
ax.plot(solvents[indexes,1],solvents[indexes,5], label='Solvents')
ax.plot(solutes[indexes,1], vx_total[indexes], label='Total')  
#ax.axhline(y=0.6, xmin = 0, xmax=1, ls=':',c='black')
#ax.axhline(y=0.15, xmin = 0, xmax=1, ls=':',c='black')
#ax.axvline(x=L[0]/2, ymin = 0, ymax=1, ls=':',c='black')
#ax.set_xlim(0, L[0])   
#ax.axvspan(L[0]/2-r_colloid,L[0]/2+r_colloid, alpha=0.5, color='red')
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$v_x(x)$")
ax.legend()
plt.tight_layout()
plt.savefig("velocities.pdf")

