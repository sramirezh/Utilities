#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 18:24:28 2018
Test of LJ melting
@author: simon
"""
import pandas as pd
import numpy as np
from lammps import IPyLammps
import matplotlib.pyplot as plt
# 3d Lennard-Jones melt

L = IPyLammps() #Creates the object
L.units("lj")
L.atom_style("atomic")
#L.atom_modify("map array")

L.region("box block", 0, 6, 0, 6, 0, 6)
L.create_box(2, "box")
L.create_atoms(1, "single 0 0 0")
L.create_atoms(1, "single 0 1.12 0")
L.create_atoms(2, "single 3 0 0")
L.create_atoms(2, "single 3 1.12 0")
L.mass('*', 1.0)

#L.velocity("all create", 1.44, 87287, "loop geom")

L.pair_style("hybrid","frenkel",2.5, "lj/cut", 2.5)
L.pair_coeff(1, 1,"frenkel", 1.0, 1.0, 2.5, 4)
L.pair_coeff(2, 2,"lj/cut", 1.0, 1.0, 2.5)
L.pair_coeff(1, 2,"lj/cut", 1.0, 1.0, 2.5)
L.command("pair_modify pair lj/cut shift yes")
L.thermo(1)
L.variable("potential_energy equal pe")
L.neighbor(0.3, "bin")
L.neigh_modify("delay 0 every 20 check no")
L.dump("id", "all", "atom", "50", "dump.melt")
L.fix("1 all nve")
L.command("compute peratom all pe/atom")
L.command("dump properties all custom 1 atom_properties.dat id type fx fy fz x y z c_peratom")
L.dump("3 all xyz 5 trajectory.xyz")
L.run(100)

L.write_script("input.lmp")

dir(L.runs[0].thermo) #Lists all the thermo parameters and also you can access them easily
