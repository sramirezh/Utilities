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

L.lattice("fcc", 0.8442)
L.region("box block", 0, 4, 0, 4, 0, 4)
L.create_box(1, "box")
L.create_atoms(1, "box")
L.mass(1, 1.0)

L.velocity("all create", 1.44, 87287, "loop geom")

L.pair_style("lj/cut", 2.5)
L.pair_coeff(1, 1, 1.0, 1.0, 2.5)

L.neighbor(0.3, "bin")
L.neigh_modify("delay 0 every 20 check no")

L.fix("1 all nve")

L.variable("fx atom fx")
L.write_script("input.lmp")
L.run(100)

L.dump("id", "all", "atom", "50", "dump.melt")
L.command("compute gr all rdf 100 1 1 ")
L.fix(" 2 all ave/time 1 100 100 c_gr[*] file gr.dat mode vector")
L.variable("g_r vector f_2")
L.variable("f_x atom fx")
fx=L.variables['fx'].value
L.dump(" da2 all custom 50 config*.dat id type x y z ix iy iz vx vy vz")
L.run(1001)


Data=pd.read_csv("gr.dat",sep=" ",skiprows=4,dtype=np.float64,header=None).values

plt.plot(Data[:,1],Data[:,2],label="lammps")
plt.legend()
plt.show()

###Comments

"""

in theory you can define the atom positions as shown in the MC example.
 L.atoms[i].position = (x+dx, y+dy)
"""