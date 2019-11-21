#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 10:43:05 2019

@author: sr802
"""


import numpy as np
import hoomd
import hoomd.hpmc
import sys
import os
from scipy.spatial.distance import pdist,squareform
import glob
sys.path.append(os.path.join(os.path.dirname(__file__), '../')) #This falls into Utilities path

import Lammps.core_functions as cf
# =============================================================================
# Input parameters
# =============================================================================

sphere_radius = 0.5
N = 100
vol_part = 4/3*np.pi*sphere_radius**3
ff = 0.5


# =============================================================================
# Simple computations
# =============================================================================

box_side = (N*vol_part/ff)**(1/3) 
n = int(np.round(N**(1/3))) #Number of particles per dimension


# getting the FCC lattice, following the jupyther nb
rho = vol_part/ff

vol_cell = 4*rho
a = vol_cell**(1/3)



n = int(np.ceil((N/4)**(1/3))) #Number of unitary boxes per dimension

def generate_configuration(steps, seed):
    """
    saves the configuration as trajectory.gsd
    
    Args:
        steps number of steps to run
        seed for the MC
    """

    hoomd.context.initialize("--mode=cpu --nthreads=12")


    #Creating the HCP
    #http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/hcp/hcp.php
    
    
    #a = 1
    #b = a
    #c = 1.633*a
    
    uc = hoomd.lattice.fcc(a = a )
    
    #uc = hoomd.lattice.unitcell(N = 4, 
    #                            a1 = [a/2,-3**(0.5)/2,0],
    #                            a2 = [a/2,3**(0.5)/2,0],
    #                            a3 = [0,0,c],
    #                            position = [[0.0,0.0,0.0],[0.5,0.5,0.0],[0.5,5.0/6.0,0.5],[0.0,1.0/3.0,0.5]],
    #                            dimensions = 3);
                                
                                
    system = hoomd.init.create_lattice(unitcell=uc, n = n)
    
    #system.box = hoomd.data.boxdim(L = box_side)
    
    mc = hoomd.hpmc.integrate.sphere(d=0.2, seed=seed)
    
    mc.shape_param.set('A', diameter= 2*sphere_radius )
    print ("The system is made of %s" %mc.get_type_shapes())
    
    
    d = hoomd.dump.gsd("trajectory.gsd", period=10, group=hoomd.group.all(), overwrite=True)
    
    
    hoomd.run(steps)
    
    snap = system.take_snapshot(all=True)
    
    
   #cf.save_instance(snap, "snap.pkl")
    
    
    
    return snap


# =============================================================================
# Creating or loading the configuration
# =============================================================================

#if len(glob.glob("*.pkl"))>0:
#    snap = cf.load_instance("snap.pkl")
#else:
    
snap = generate_configuration(10000,1231432)
positions = snap.particles.position

real_radius = 150 *10**-3 #Radius in microns
def convert_into_real(length,radius,dim):
    """
    Converts into microns
    Args:
        quantity: coulde be an array or single number, in length, area, volume
        real_radius in microns
        dim 1 if length, 3 if volume
    """
    scaling = radius*2
    length = length * (scaling)**dim
    
    return length
    
positions = convert_into_real(positions, real_radius,1)

np.savetxt("real_positions.txt",positions)

pos_m = squareform(pdist(positions))
np.fill_diagonal(pos_m,100)

# =============================================================================
# Basic checks in real units
# =============================================================================

print ("\nThe minimum distance between spheres is %s" %np.min(pos_m))
print ("The box side is %f"%convert_into_real(snap.box.Lx,real_radius,1))

new_ff = snap.particles.N*vol_part/snap.box.get_volume()

print ("The new ff is %lf" %new_ff)


# =============================================================================
# With MDAAnalysis
# =============================================================================

import MDAnalysis
import MDAnalysis as mda
from MDAnalysis.coordinates.GSD import GSDReader
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt


simulation = GSDReader('trajectory.gsd')


#for ts in simulation:
#    print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, simulation.trajectory.time))
    
    
#I can acces the positions with ts.positions
#to see how are the atoms
    
u = mda.Universe("trajectory.gsd")
group = u.select_atoms("type A")   #Checked  u.atoms[0] to see how are the atoms

com1 = group.center_of_mass()

discard_perc = 0.5

discard_frames = int(discard_perc*u.trajectory.n_frames)

rdf = InterRDF(group,group,exclusion_block = (1,1),range = [0,2.5],verbose = True, start = discard_frames)

rdf.run()

plt.plot(rdf.bins, rdf.rdf)
plt.savefig("gr.pdf")
plt.show()