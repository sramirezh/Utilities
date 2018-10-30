#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 09:29:22 2018
This scripts computes the mean-squared distance both from lammps and from the trajectories
@author: sr802
"""

import pandas as pd
import numpy as np
from lammps import IPyLammps 
import matplotlib.pyplot as plt


L=IPyLammps()
L.units("lj")
L.atom_style("atomic")

L.variable("L equal 20")
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
L.command("compute msd all msd")
L.variable("meansqrdisp equal m")
L.command("thermo 10")
L.command("thermo_style custom step temp press c_msd[4]")

L.run(1000)


