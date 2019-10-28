#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 10:43:05 2019

@author: sr802
"""


import numpy as np
import hoomd
import hoomd.hpmc

from scipy.spatial.distance import pdist,squareform

# =============================================================================
# Input parameters
# =============================================================================

sphere_radius = 0.5
N = 100
vol_part = 4/3*np.pi*sphere_radius**3
ff = 0.6


# =============================================================================
# Simple computations
# =============================================================================

box_side = (N*vol_part/ff)**(1/3) 
n = int(np.round(N**(1/3))) #Number of particles per dimension


# getting the FCC lattice, following the jupyther nb
rho = vol_part/ff

vol_cell = 4*rho
a = vol_cell**(1/3)




# =============================================================================
# Simple computations
# =============================================================================

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
                            
                            
system = hoomd.init.create_lattice(unitcell=uc, n = 3)

#system.box = hoomd.data.boxdim(L = box_side)

mc = hoomd.hpmc.integrate.sphere(d=0.2, seed=1)

mc.shape_param.set('A', diameter=1)
print ("The system is made of %s" %mc.get_type_shapes())


d = hoomd.dump.gsd("trajectory.gsd", period=10, group=hoomd.group.all(), overwrite=True)


hoomd.run(10000)

snap = system.take_snapshot(all=True)

positions = snap.particles.position

pos_m = squareform(pdist(positions))
np.fill_diagonal(pos_m,100)
print ("\nThe minimum distance between spheres is %s" %np.min(pos_m))

