#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:04:24 2020

@author: simon
"""

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt

import MDAnalysis.analysis.rdf as rdf


u = mda.Universe("system.data", "dcd_nvt.dcd", format="LAMMPS")  # The universe for all the atoms





n_molecules = u.atoms.n_residues
molecules = u.atoms.residues # I can select them here and then iterate later through the trajectory



time_steps = u.trajectory.n_frames

centroids_traj = np.empty((time_steps, n_molecules, 3 ))

distances = []

# =============================================================================
#  Unwtrapping the coordinates
# =============================================================================
ag = u.atoms
transform = mda.transformations.unwrap(ag)
u.trajectory.add_transformations(transform)

for i,ts in enumerate(u.trajectory):
    distances_t = []
    for j,mol in enumerate(molecules):
        centroids_traj[i, j,:] = mol.atoms.centroid(compound='residues')
        pos = mol.atoms.positions
        distances_t.append(np.linalg.norm(pos[1]-pos[0]))
    
    distances.append(distances_t)


array_dist = np.transpose(np.array(distances))

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
plt.show()

    

# Integrating to get the average number of particles at a given radius

r_max = 1.1+8.0

