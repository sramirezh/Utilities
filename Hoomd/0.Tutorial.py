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
N = 108
vol_part = 4/3*np.pi*sphere_radius**3
ff = 0.6


# =============================================================================
# Simple computations
# =============================================================================

box_side = (N*vol_part/ff)**(1/3) 
n = int(np.round(N**(1/3))) #Number of particles per dimension

hoomd.context.initialize("--mode=cpu --nthreads=12")


system = hoomd.init.create_lattice(unitcell=hoomd.lattice.fcc(a=2), n=3)

system.box = hoomd.data.boxdim(L = box_side)

mc = hoomd.hpmc.integrate.sphere(d=0.2, seed=1)

mc.shape_param.set('A', diameter=0.5)
print ("The system is made of %s" %mc.get_type_shapes())


d = hoomd.dump.gsd("trajectory.gsd", period=10, group=hoomd.group.all(), overwrite=True)


hoomd.run(1000)

mc.shape_param.set('A', diameter=1.0)

hoomd.run(1000)


snap = system.take_snapshot(all=True)

positions = snap.particles.position

pos_m = squareform(pdist(positions))


