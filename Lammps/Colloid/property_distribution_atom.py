#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 17:49:52 2020
BAsed on density_distribuion_from_atom.py, computes the property distributions from the dump files from lammps using ovito
@author: sr802
"""

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


def property_chunk(prop,region_limits,n_bins = 50, v_weight = True ):
    
    """
    computes the property distribution for now only on the x direction
    Args:
        pos: numpy arry with the particles position
        region_limits array with the definitions of the region [xmin,xmax,ymin,ymax...]
        v_weight True if we want to weight the property by the volume of the chunck
    
    """
    indexes = np.where((prop[:,0]> region_limits[0] )&(prop[:,0] < region_limits[1])
              &(prop[:,1]> region_limits[2] )&(prop[:,1] < region_limits[3])
              &(prop[:,2]> region_limits[4] )&(prop[:,2] < region_limits[5]))

    hist, bin_edges = np.histogram(prop[indexes,0], bins = nbins,range = [region_limits[0],region_limits[1]])
    left_edge = bin_edges[:-1] 
    right_edge = bin_edges[1:]
    
    vol = 1
    if v_weight == True:
        vol = ((right_edge-left_edge)*(region_limits[3]-region_limits[2])*(region_limits[5]-region_limits[4]))
    
    prop_dist = hist/ vol
    
    
    return left_edge, prop_dist
    


m_density_solu = []
m_density_solv = []


m_density_c_solu = []
m_density_c_solv = []

discard = 0.3
n_frames = node.source.num_frames

region = [-L[0]/2, L[0]/2, -L[1]/2, L[1]/2, -L[2]/2, L[2]/2] #for swol type
region = [0, L[0], 0, L[1], 25, L[2]]


ave_vel_m = []
ave_vel_solu_m = []
for frame in range(n_frames):
    print ("Analysing frame %s of %s)"%(frame,n_frames))
    
    if frame>discard*n_frames:
    
        data = node.compute(frame)
        pos = data.particle_properties.position.array
        velocity = data.particle_properties.velocity.array[:,0]
        types = data.particle_properties.particle_type.array
    
        ind_solute = np.where(types == 2)[0]
        ind_solvent = np.where(types == 1)[0]
        
        # Velocities
        ave_vel = (np.sum(velocity[ind_solute])+np.sum(velocity[ind_solvent]))/(len(velocity[ind_solute])+len(velocity[ind_solvent]))
        ave_vel_m.append(ave_vel)
        ave_vel_solu_m.append(np.average(velocity[ind_solute]))
        
        pos_solute = pos[ind_solute]
        pos_solvent = pos[ind_solvent]
        
        x,density_c_solu = property_chunk(pos_solute, region)
        x,density_c_solv = property_chunk(pos_solvent, region)
        
        radii, density_solu = density_distribution(pos_solute, center)
        radii, density_solv = density_distribution(pos_solvent, center)
        
        m_density_solu.append(density_solu)
        m_density_solv.append(density_solv)
        
        m_density_c_solu.append(density_c_solu)
        m_density_c_solv.append(density_c_solv)
    

# Velocities

print ("The velocity for the solution in all the system is %s"%np.average(ave_vel_m))
print ("The velocity for the solutes in all the system is %s"%np.average(ave_vel_solu_m))


rho_solu = np.average(m_density_solu, axis = 0)
rho_solv = np.average(m_density_solv, axis = 0)
rho_total = rho_solu+rho_solv

# From density chunks
rho_solu_c = np.average(m_density_c_solu, axis = 0)
rho_solv_c = np.average(m_density_c_solv, axis = 0)
rho_total_c = rho_solu_c+rho_solv_c

# Getting the averages concentrations in the bulk assumed to be outside r=6

ind_bulk =np.where(radii>6)

cs_b = np.average(rho_solu[ind_bulk])
ct_b = np.average(rho_total[ind_bulk])

#Plotting the results for the radial

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


#Plotting the results for the chunk
lattice_constant = (4/0.8)**(1/3)

fig,ax = plt.subplots()

r_colloid  = 3.23

ax.plot(x,rho_solu_c, label='Solutes')
ax.plot(x,rho_solv_c, label='Solvents')
ax.plot(x, rho_total_c, label='Total')  
ax.axhline(y=0.6, xmin = 0, xmax=1, ls=':',c='black')
ax.axhline(y=0.15, xmin = 0, xmax=1, ls=':',c='black')
ax.axvline(x=L[0]/2, ymin = 0, ymax=1, ls=':',c='black')
ax.set_xlim(0, L[0])   
ax.axvspan(L[0]/2-r_colloid,L[0]/2+r_colloid, alpha=0.5, color='green')
ax.axvspan(6*lattice_constant,9*lattice_constant, alpha=0.5, color='blue')
ax.axvspan(21*lattice_constant,24*lattice_constant, alpha=0.5, color='red')

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$c(x)$")
ax.legend()
plt.tight_layout()
plt.savefig("conc_chunks.pdf")

