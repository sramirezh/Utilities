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
import Lammps.lammps_utilities as lu


class SimulationEMD(lu.SimulationType):
    def __init__(self, name, T, eta, grad_mu_s):
        self.name = name
        self.T = T
        self.eta = eta
        self.grad_mu_s = grad_mu_s
    
        
            

# =============================================================================
# Main
# =============================================================================
logger = cf.log(__file__, os.getcwd())    

# Write all the new types of simulations here
Yoshida = SimulationEMD("Yoshida0125", 1, 1.57, -0.125)

# Reading all the final averaged chunks
solute = da.DensityDistribution("Sproperties_short.dat",'rBulk')
solvent = da.DensityDistribution("Fproperties_short.dat",'rBulk')
fluid = da.DensityDistribution("properties_short.dat",'rBulk')

# Loading the data from the simualtion
sim_log = lu.Simulation("log.lammps")
z_wall = sim_log.limits[0, 2]
Kb = 1

#  define the type of simulation
sim = Yoshida
sim.print_params(logger)

# =============================================================================
# Computing the analytic properties
# =============================================================================
solute.compute_all_properties(solute.lower_limit)
solvent.compute_all_properties(solvent.lower_limit)




# =============================================================================
# Getting the velocities as in Anderson1984
# =============================================================================
# For the solutes
grad_c_s = sim.grad_mu_s * solute.rho_bulk
v_0_s = - (Kb* sim.T)* solute.l*solute.k/sim.eta
vx_s = v_0_s * grad_c_s


# For the solvents
grad_mu_f = - (solute.rho_bulk / solvent.rho_bulk)* sim.grad_mu_s
grad_c_f = grad_mu_f * solvent.rho_bulk
v_0_f = - (Kb* sim.T)* solvent.l*solvent.k/sim.eta
vx_f = v_0_f * grad_c_f

vx_total = vx_s + vx_f

logger.info("The contribution to the velocity due to the solutes is %s"%vx_s)
logger.info("The contribution to the velocity due to the solvents is %s"%vx_f)
logger.info("The predicted DO velocity is  %s" % (vx_total))

## =============================================================================
## Computing the velocity using both species
## =============================================================================
#print("The property distributions from the solutes are:" )
#print(solute.property_dist)
##solute.plot_property_dist()
#
#
#solute.get_k(0)
#solute.plot_property_dist("integrand_k", 0,8)
#
#
#
#
#
##  Integral
#
#
#
#cf.integrate()


