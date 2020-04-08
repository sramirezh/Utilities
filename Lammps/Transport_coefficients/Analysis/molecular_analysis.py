#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:04:24 2020

@author: simon
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations

import MDAnalysis.analysis.rdf as rdf
import matplotlib.pyplot as plt

u = mda.Universe("system.data", "dcd_nvt.dcd", format="LAMMPS")  # The universe for all the atoms





n_molecules = u.atoms.n_residues
molecules = u.atoms.residues # I can select them here and then iterate later through the trajectory



time_steps = u.trajectory.n_frames

centroids_traj = np.empty((time_steps, n_molecules, 3 ))

distances = []

# Unwtrapping the coordinates
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

uc.load_new(centroids_traj)


with mda.coordinates.XYZ.XYZWriter("trajectory.xyz") as W:
    for ts in uc.trajectory:
        W.write(uc.atoms)


#I can write an XYZ and then read it again or just analyse on the flight

#    # analyze frame
#    if take_this_frame == True:
#    with mda.Writer('frame.data') as W:
#            W.write(u.atoms)
#         break

#u= mda.Universe("1.cxyz",format='xyz')
#
#monomers=u.select_atoms("type 3")
#solutes= u.select_atoms("type 2")
#
#g_r=rdf.InterRDF(monomers,solutes)
#g_r.run()
#
#plt.plot(g_r.bins, g_r.rdf)
#plt.show()

    
#g_r = rdf.InterRDF()