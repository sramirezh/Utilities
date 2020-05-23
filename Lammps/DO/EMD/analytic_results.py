#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 20 13:42:54 2020
Script to get the theoretical estimation of the diffusio-osmotic velocity

@author: simon
"""

import numpy as np
import sys
import os
from shlex import split
import matplotlib.pyplot as plt
import linecache
import re
import density_analysis as da
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

# =============================================================================
# Main
# =============================================================================


solute = da.DensityDistribution("Sproperties_short.dat",'rBulk')
solvent = da.DensityDistribution("Fproperties_short.dat",'rBulk')
fluid = da.DensityDistribution("properties_short.dat",'rBulk')

# Parameters
Kb = 1
T = 1
eta =  1.57
grad_mu_s = - 0.125


z_wall = 0

solute.compute_all_properties(solute.lower_limit)
solvent.compute_all_properties(solvent.lower_limit)

# For the solutes
grad_c_s = grad_mu_s * solute.rho_bulk
v_0_s = - (Kb* T)* solute.l*solute.k/eta
vx_s = v_0_s * grad_c_s


# For the solvents
grad_mu_f = - (solute.rho_bulk / solvent.rho_bulk)* grad_mu_s
grad_c_f = grad_mu_f * solvent.rho_bulk
v_0_f = - (Kb* T)* solvent.l*solvent.k/eta
vx_f = v_0_f * grad_c_f



vx_total = vx_s + vx_f


print (vx_s)
print (vx_f)


print ("total is %s" % (vx_total))


# =============================================================================
# Computing the velocity using both species
# =============================================================================
print("The property distributions from the solutes are:" )
print(solute.property_dist)
#solute.plot_property_dist()


solute.get_k(0)
solute.plot_property_dist("integrand_k", 0,8)





#  Integral



cf.integrate()


