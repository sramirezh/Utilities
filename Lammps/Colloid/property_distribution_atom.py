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
from joblib import Parallel, delayed
import multiprocessing

from Lammps.PDP.Plots.LJ_plotter import LJ
import Lammps.core_functions as cf


from ovito.modifiers import *
import ovito.io as ov





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



def reduce_to_region(pos,region_limits):
    """
    Gets the indexes for particles inside a given region
    """
    ind_region = np.where((pos[:,0]> region_limits[0] )&(pos[:,0] < region_limits[1])
          &(pos[:,1]> region_limits[2] )&(pos[:,1] < region_limits[3])
          &(pos[:,2]> region_limits[4] )&(pos[:,2] < region_limits[5]))
    
    return ind_region
    
    
    
def property_chunk(pos, prop ,region_limits, n_bins = 50, v_weight = True ):
    global prop_dist, count,bin_edges,hist, edges, ind_region
    """
    computes the property distribution for now only on the x direction
    Args:
        pos: numpy arry with the particles position
        property : property to be measured, BE careful that if you want the number density, the property is a ones array
        region_limits array with the definitions of the region [xmin,xmax,ymin,ymax...]
        v_weight True if we want to weight the property by the volume of the chunck
    
    """

    ind_region = reduce_to_region(pos,region_limits)
    
    # Reducing the arrays to the values withing the region
    prop = prop[ind_region]
    pos = pos[ind_region]
    
    bin_edges = np.linspace(region_limits[0],region_limits[1],n_bins+1)
    left_edge = bin_edges[:-1] 
    right_edge = bin_edges[1:]
    
    hist, edges = np.histogram(pos[:,0],bins =n_bins, range=[region_limits[0],region_limits[1]])
    prop_dist = np.zeros(n_bins)
    count = np.zeros(n_bins)

    for i,left in enumerate(left_edge):
        indexes = np.where((pos[:,0]>left) & (pos[:,0] <= right_edge[i]))[0]
        count[i] = len(indexes)
        prop_dist[i] = np.sum(prop[indexes])
        
    
    if v_weight == True: #Normalising with the volume
        norm = ((right_edge-left_edge)*(region_limits[3]-region_limits[2])*(region_limits[5]-region_limits[4]))
    else:
        count[count == 0 ] = 1 # To avoid dividing by zero
        norm = count
    
    prop_dist = prop_dist/ norm
    
    
    return left_edge, prop_dist
    


# =============================================================================
# MAIN
# =============================================================================

num_cores = multiprocessing.cpu_count()


#External parameters 
lattice_constant = (4/0.8)**(1/3)
r_colloid  = 3.23

node = ov.import_file("all.atom", multiple_frames = True)
#node = ov.import_file("equil.dat", multiple_frames = False)

#Getting the properties of the box

#Assuming that the box does not change size

data = node.compute(0)
list(data.particle_properties.keys())

box = node.compute(0).cell # getting the properties of the box
L = np.diag(box.matrix) 
center = L/2.0




region = [-L[0]/2, L[0]/2, -L[1]/2, L[1]/2, -L[2]/2, L[2]/2] #for swol type
region = [0, L[0], 0, L[1], 25, L[2]] # In the top bulk



# density Arrays
m_density_solu = []
m_density_solv = []

m_density_c_solu = []
m_density_c_solv = []


# velocity arrays
m_velocity_c_solu = []
m_velocity_c_solv = []
m_velocity_c_sol = []

#Global analysis
vel_m = []
vel_solu_m = []
vel_solv_m = []

# For the region
vel_r = []
vel_solu_r = []
vel_solv_r = []


#Parallel(n_jobs=num_cores,verbose=10)(delayed(compute_one_time)(pos_init,fil) for fil in input_files)










#
#def analyse_one_frame(frame,region):
#    
#    """
#    Analyses the radial distribution of densities and chunk distribution of velocities and densities in a chunk
#    
#    Args:
#        Frame is just the 
#    """
    

discard = 0


n_frames = node.source.num_frames

# To use either total number or percentage to discard
if discard< 1:
    discard = discard*n_frames 
    
for frame in range(n_frames):
    print ("Analysing frame %s of %s)"%(frame,n_frames))
    
    if frame >= discard:
    
        data = node.compute(frame)
        pos = data.particle_properties.position.array
        velocity = data.particle_properties.velocity.array[:,0]
        types = data.particle_properties.particle_type.array
    
        ind_solute = np.where(types == 2)[0]
        ind_solvent = np.where(types == 1)[0]
        

        pos_solute = pos[ind_solute]
        pos_solvent = pos[ind_solvent]
        
        # Chunk density Analysis 
        x,density_c_solu = property_chunk(pos_solute, np.ones(len(ind_solute)), region)
        x,density_c_solv = property_chunk(pos_solvent, np.ones(len(ind_solvent)), region)
        
        
        # Radial density Analysis
        radii, density_solu = density_distribution(pos_solute, center)
        radii, density_solv = density_distribution(pos_solvent, center)
        
        m_density_solu.append(density_solu)
        m_density_solv.append(density_solv)
        
        m_density_c_solu.append(density_c_solu)
        m_density_c_solv.append(density_c_solv)
        
        
        # Chunk velocity analysis
        
        velocity_solute = velocity[ind_solute]
        velocity_solvent = velocity[ind_solvent]

        
        x,velocity_c_solu = property_chunk(pos_solute, velocity_solute, region, v_weight = False)
        x,velocity_c_solv = property_chunk(pos_solvent, velocity_solvent, region, v_weight = False)
        x,velocity_c_sol = property_chunk(pos, velocity, region, v_weight = False)
        
        
        m_velocity_c_solu.append(velocity_c_solu)
        m_velocity_c_solv.append(velocity_c_solv)
        m_velocity_c_sol.append(velocity_c_sol)  
        
        # Velocities for the entire system
        vel = (np.sum(velocity[ind_solute])+np.sum(velocity[ind_solvent]))/(len(velocity[ind_solute])+len(velocity[ind_solvent]))
        vel_m.append(vel)
        vel_solu_m.append(np.average(velocity[ind_solute]))
        
        #Velocity for the top region
        ind_solute_r = reduce_to_region(pos_solute,region)
        ind_solvent_r =reduce_to_region(pos_solute,region)
        ind_sol_r = reduce_to_region(pos, region)
        
        vel_solu_r.append(np.average(velocity[ind_solute_r]))
        vel_solv_r.append(np.average(velocity[ind_solvent_r]))
    
        ind_solution_r = np.append(ind_solute_r[0],ind_solvent_r[0])
        vel_r.append(np.average(velocity[ind_sol_r]))
        

  
#    

# Velocities

print ("The velocity for the solution in all the system is %s"%np.average(vel_m))
print ("The velocity for the solutes in all the system is %s"%np.average(vel_solu_m))

print ("The velocity for the solution in the region is %s"%np.average(vel_r))
print ("The velocity for the solutes in the region is %s"%np.average(vel_solu_r))
print ("The velocity for the solvents in the region is %s"%np.average(vel_solv_r))


rho_solu = np.average(m_density_solu, axis = 0)
rho_solv = np.average(m_density_solv, axis = 0)
rho_total = rho_solu+rho_solv

# From density chunks
rho_solu_c = np.average(m_density_c_solu, axis = 0)
rho_solv_c = np.average(m_density_c_solv, axis = 0)
rho_total_c = rho_solu_c+rho_solv_c


# From velocity chunks

vel_solu_c = np.average(m_velocity_c_solu, axis = 0)
vel_solv_c = np.average(m_velocity_c_solv, axis = 0)
vel_sol_c = np.average(m_velocity_c_sol, axis = 0)


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
ax.axvline(x = r_colloid, ymin = 0, ymax=1, ls=':',c='black')
ax.legend()
plt.tight_layout()
plt.savefig("conc.pdf")


#Plotting the densities for the chunk



fig,ax = plt.subplots()



ax.plot(x,rho_solu_c, label='Solutes', ls = '--')
ax.plot(x,rho_solv_c, label='Solvents', ls = '--')
#ax.plot(x, rho_total_c, label='Total')  
ax.axhline(y=0.6, xmin = 0, xmax=1, ls='--',c='black')
ax.axhline(y=0.15, xmin = 0, xmax=1, ls='--',c='black')
#ax.axvline(x=L[0]/2, ymin = 0, ymax=1, ls=':',c='black')
ax.set_xlim(0, L[0])   
ax.axvspan(L[0]/2-r_colloid,L[0]/2+r_colloid, alpha=0.5, color='green')
#ax.axvspan(6*lattice_constant,9*lattice_constant, alpha=0.5, color='blue')
#ax.axvspan(21*lattice_constant,24*lattice_constant, alpha=0.5, color='red')

ymin,ymax=ax.get_ylim()
ax.set_ylim(0,ymax*1.3)

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$c^B(x)$")
ax.legend()
plt.tight_layout()
plt.savefig("conc_chunks.pdf")


epsilon =os.getcwd().split('/')[-1]
np.savetxt("distribution_%s.dat"%epsilon, np.column_stack((x,rho_solu_c,rho_solv_c)), header = "x rho_solu rho_solv" )

#Plotting the velocities for the chunk


fig,ax = plt.subplots()



ax.plot(x,vel_solu_c, label='Solutes')
ax.plot(x,vel_solv_c, label='Solvents')
ax.plot(x, vel_sol_c, label='Total')  
ax.axvline(x=L[0]/2, ymin = 0, ymax=1, ls=':',c='black')
ax.set_xlim(0, L[0])   
ax.axvspan(L[0]/2-r_colloid,L[0]/2+r_colloid, alpha=0.5, color='green')
#ax.axvspan(6*lattice_constant,9*lattice_constant, alpha=0.5, color='blue')
#ax.axvspan(21*lattice_constant,24*lattice_constant, alpha=0.5, color='red')

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$v_x(x)$")
ax.legend()
plt.tight_layout()
plt.savefig("velocity_chunks.pdf")


#Plotting the flux for the chunk


fig,ax = plt.subplots()



ax.plot(x,vel_solu_c*rho_solu_c, label='Solutes')
ax.plot(x,vel_solv_c*rho_solv_c, label='Solvents')
ax.plot(x, vel_sol_c*rho_total_c, label='Total')  
ax.axvline(x=L[0]/2, ymin = 0, ymax=1, ls=':',c='black')
ax.set_xlim(0, L[0])   
ax.axvspan(L[0]/2-r_colloid,L[0]/2+r_colloid, alpha=0.5, color='green')
ax.axvspan(6*lattice_constant,9*lattice_constant, alpha=0.5, color='blue')
ax.axvspan(21*lattice_constant,24*lattice_constant, alpha=0.5, color='red')

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$J_x(x)$")
ax.legend()
plt.tight_layout()
plt.savefig("flux_chunks.pdf")


