#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 14:23:36 2019
Analyse the particle density distribution using ovito!
    
This scripts plots the density distribution of the solvents and solutes using a lammps trajectory

Useful links

https://ovito.org/manual/python/modules/ovito_data.html
https://www.ovito.org/manual/python/modules/ovito_modifiers.html


An important modifier is 
class ovito.modifiers.ComputePropertyModifier

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
import numpy as np

nbins = 100

node = ov.import_file("all.atom", multiple_frames = True)



data = node.compute(0)
list(data.particle_properties.keys())


#Getting the properties of the box

#Assuming that the box does not change size
box = node.compute(0).cell # getting the properties of the box
L = np.diag(box.matrix) 
center = L/2.0


def density_distribution(pos, center, nbins = 40, rmin = 0, rmax = 8):
    """
    Get the density distribution of particles measured from the center
    Args:
        pos: numpy arry with the particles position
        center: vector with the center position [X,Y,Z]
        
    returns:
        positions
        density
    """
    r = np.linalg.norm(pos-center,axis = 1 )
    hist, bin_edges = np.histogram(r, bins = nbins,range = [rmin,rmax])
    radii = bin_edges[:-1]
    radii_right = bin_edges[1:]
    factor = 4./3. * np.pi
    rho_dist = hist / (factor * (radii_right**3 - radii**3))
    
    return radii, rho_dist


m_density_solu = []
m_density_solv = []


discard = 1000


n_frames = node.source.num_frames

# To use either total number or percentage to discard
if discard< 1:
    discard = discard*n_frames

for frame in range(n_frames):
    
    if frame >= discard:
    
        print ("Analysing frame %s of %s)"%(frame,n_frames))
        
        
        data = node.compute(frame)
        pos = data.particle_properties.position.array
        types = data.particle_properties.particle_type.array
        print ("The number of solutes is %s"%len(np.where(types == 2)[0]))
        ind_solute = np.where(types == 2)[0]
        ind_solvent = np.where(types == 1)[0]
        
        pos_solute = pos[ind_solute]
        pos_solvent = pos[ind_solvent]
        
        radii, density_solu = density_distribution(pos_solute, center)
        radii, density_solv = density_distribution(pos_solvent, center)
        
        m_density_solu.append(density_solu)
        m_density_solv.append(density_solv)
    

rho_solu = np.average(m_density_solu, axis = 0)
rho_solv = np.average(m_density_solv, axis = 0)
rho_total = rho_solu+rho_solv

# Getting the averages concentrations in the bulk assumed to be outside r=6

ind_bulk =np.where(radii>6)

cs_b = np.average(rho_solu[ind_bulk])
ct_b = np.average(rho_total[ind_bulk])

#Plotting the results

cf.set_plot_appearance()

fig,ax = plt.subplots()




ax.plot(radii,rho_solu, label='Solutes')
ax.plot(radii,rho_solv, label='Solvents')
ax.plot(radii,rho_solu+rho_solv, label='Total')
ax.set_xlim(0, 8)   
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$c(r)$")
ax.axvline(x=3.23, ymin = 0, ymax=1, ls=':',c='black')
ax.legend()
plt.tight_layout()
plt.savefig("conc.pdf")


#To analyse with colmobility need to add some columns to not change the col_mobility_anderson.py

n = len(radii)
filler = np.zeros(n)
data_save = np.column_stack((filler,radii,filler,rho_solu))
np.savetxt("prof_u_atom.dat", data_save )




