#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:04:24 2020

@author: simon
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.transformations as tr
import MDAnalysis.analysis.rdf as rdf
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
from scipy.spatial.distance import pdist,squareform



# =============================================================================
# Input parameters
# =============================================================================

r_max = 1.1+8.0

print ("Using r_max=%s"%r_max)


u = mda.Universe("system.data", "dcd_nvt.dcd")  # The universe for all the atoms
v = mda.Universe("system.data","velocities.dat", format = "LAMMPSDUMP" ) # Reading all the velocities

scale = v.dimensions[0]


n_molecules = u.atoms.n_residues
molecules = u.atoms.residues # I can select them here and then iterate later through the trajectory



time_steps = u.trajectory.n_frames

centroids_traj = np.empty((time_steps, n_molecules, 3 ))

distances = []

# =============================================================================
#  Unwtrapping the coordinates
# =============================================================================
ag = u.atoms
transform = tr.unwrap(ag)
u.trajectory.add_transformations(transform)


f = open("colliding.dat","w")

list_ni = [] # has the time step and the list of non-interacting molecules
ratio = []
for i,ts in enumerate(u.trajectory):
    
    v.trajectory[i]
    
    #Adding velocities to the universe 
    ts.has_velocities = True      
    # TODO mda reads lammpsdumps and scales them
#    print (v.atoms.positions/scale)
#    print ("***********************************")
    
    u.atoms.velocities = v.atoms.positions/ scale        
    
#    print (u.atoms.velocities)
        
    list_ni_t = ["Frame %s"%ts.time]
    # Computing the centroid
    centroid_pos = molecules.atoms.centroid(compound='residues')
    centroids_traj[i, :,:]  = centroid_pos
    
    # Measuring the positions for all the molecules
    pos_m = squareform(pdist(centroid_pos))
    indexes = np.where(pos_m>r_max)
    
    # seting all the non interacting particles to zero
    pos_m[indexes] = 0
    
    suma = np.sum(pos_m, axis = 1)
    
    ni_molecules = np.where(suma == 0)[0] # Non-interacting molecules
    
    
    # Writing the properties
    f.write("%s\n"%len(ni_molecules))
    for j in ni_molecules:
        res = u.atoms.residues[j]
        for atom in res.atoms:
            f.write("%s %s %s %s %s %s %s %s \n"%(atom.id,j,atom.position[0],atom.position[1],atom.position[2],atom.velocity[0],atom.velocity[1],atom.velocity[2]))
        
        
    
    
    ratio.append(len(ni_molecules)/n_molecules) # Ratio of non interacting molecules
    
    list_ni_t.extend(ni_molecules)
    
    list_ni.extend(np.array(list_ni_t))
    
f.close()
array_dist = np.transpose(np.array(distances))

print ("The average percentage of non-interacting molecules is: %s" %np.average(ratio))

arr = np.ravel(np.array(list_ni))

print ("created the list of nonpinteracting molecules for each step in ni_list.dat" )
np.savetxt("ni_list.dat", arr , fmt = "%s")


# =============================================================================
# Creating an universe for the molecules
# =============================================================================

uc = mda.Universe.empty(n_molecules, trajectory = True)
uc.dimensions = u.dimensions

uc.load_new(centroids_traj)


positions = []
with mda.coordinates.XYZ.XYZWriter("trajectory.xyz") as W:
    for ts in uc.trajectory:
        ts.dimensions = u.dimensions
        uc.atoms.pack_into_box(box = u.dimensions)
        W.write(uc.atoms)
        positions.append(uc.atoms.positions)
        
pos_tensor = np.array(positions)



# =============================================================================
# Computing the RDF
# =============================================================================
#
g_r=rdf.InterRDF(uc.atoms,uc.atoms, exclusion_block=(1, 1))
g_r.run()

counts = g_r.count
plt.plot(g_r.bins, g_r.rdf)
plt.savefig("gr.pdf")

    

# Integrating to get the average number of particles at a given radius


volume = uc.trajectory[0].volume # Assuming a constant volume
rho_ave = uc.atoms.n_atoms/volume

integrand = 4*np.pi*g_r.bins**2*g_r.rdf



N_interacting = cf.integrate(g_r.bins,integrand,0,r_max)*rho_ave

print("The average number of molecules interacting with each molecule is %s" %N_interacting)






