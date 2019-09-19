# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 16:15:07 2019

@author: sr802
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

l_box = 20
L.lattice("fcc", 0.8)
L.command("region box block  0 %f 0 %f 0 %f units box"%(l_box,l_box, l_box))
L.create_box(1, "box")
L.create_atoms(1, "box")
L.mass(1, 1.0)
L.image()

# This has to be defined to be able to run
L.pair_style("lj/cut", 2.5)
L.pair_coeff(1, 1, 1.0, 1.0, 2.5)

n_atoms = L.system.natoms
v = L.eval ('vol')

print(("The density of the system is:%f"%(n_atoms/v)))


L.command("region r_nocolloid sphere 10 10 10 3 side out")

L.group("gdel region r_nocolloid")

L.command("delete_atoms region r_nocolloid")

L.command("dump all custom 1 trajectory.atom id type x y z ix iy iz")

L.run(0)