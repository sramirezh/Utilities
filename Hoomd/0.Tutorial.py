#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 10:43:05 2019

@author: sr802
"""



import hoomd
import hoomd.hpmc

hoomd.context.initialize("--mode=cpu")

system = hoomd.init.create_lattice(unitcell=hoomd.lattice.sc(a=1.05), n=10)

mc = hoomd.hpmc.integrate.sphere(d=0.2, seed=1)

mc.shape_param.set('A', diameter=1.0)



d = hoomd.dump.gsd("trajectory.gsd", period=10, group=hoomd.group.all(), overwrite=True)


hoomd.run(100)



import ex_render
ex_render.display_movie(ex_render.render_disk_frame, 'trajectory.gsd');