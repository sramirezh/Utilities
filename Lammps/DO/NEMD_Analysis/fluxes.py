#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:42:01 2019

This computes the solute and total fluxes for diffusio-osmotic simulations

TODO created copying flux gatherer on porous

@author: simon
"""

# TODO This line is added such that spyder works!
import os
import sys
Utilities_path=os.path.join(os.path.dirname(__file__), '../../../')
sys.path.append(Utilities_path) #This falls into Utilities path
from Lammps.Pore.qsub import simulation_results as sr
from uncertainties import ufloat,unumpy
import Lammps.lammps_utilities as lu

cwd = os.getcwd() #current working directory


## =============================================================================
## Chemical potential simulations
## =============================================================================


#  TODO Below is an specific case using pore/qsub simulation results
root_pattern="."
directory_pattern='.'
parameter_id='mu'


dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}


# Running the required analysis 
cna = sr.check_n_analyse(root_pattern,directory_pattern)
cna.check_finished("vdata.dat")
cna.check_stat("statistics.dat")
cna.stat_analysis("vdata.dat")
cna.check_thermo("thermo.dat")
cna.thermo_analysis("log.lammps")


# Creating the simulation instance
sim = sr.simulation("grad_mu",0.5)
properties = sr.read_properties(root_pattern,["statistics.dat","thermo.dat"])
sim.add_properties(properties)



# Getting the dimensions of the simulation box and the regions
box_volume, box_limits = lu.read_box_limits('log.lammps')

bulk_heigth  = lu.read_region_height('rBulk')

#This is a region defined as the bulk+interface
sys_heigth  = lu.read_region_height('rSystem') 

vol_bulk = (box_limits[1][0]-box_limits[0][0])*(box_limits[1][1]-box_limits[0][1])*bulk_heigth[0]
vol_sys = (box_limits[1][0]-box_limits[0][0])*(box_limits[1][1]-box_limits[0][1])*sys_heigth[0]


# =============================================================================
# Quantities in the bulk
# =============================================================================
n_s_bulk = sim.get_property('cBSolu',True)[1][0][0]
n_f_bulk = sim.get_property('cBSolv',True)[1][0][0]
n_fluid_bulk = n_f_bulk + n_s_bulk

rho_bulk = n_fluid_bulk/vol_bulk
cs_bulk = n_s_bulk/vol_bulk

vx_s_bulk = ufloat(sim.get_property('vxB_Solu')[1][0][0],sim.get_property('vxB_Solu')[1][0][1])
J_s_bulk = cs_bulk*vx_s_bulk
Q_bulk = ufloat(sim.get_property('vxB_Sol',exact=True)[1][0][0],sim.get_property('vxB_Sol',exact=True)[1][0][1])

exc_s_bulk = J_s_bulk- cs_bulk*Q_bulk

sim.add_property('Q_bulk',Q_bulk)
sim.add_property('J_s_bulk',J_s_bulk)
sim.add_property('J_s_exc_bulk',exc_s_bulk)



# =============================================================================
# Quantities in the entire Volume
# =============================================================================
n_s = sim.get_property('cSolu',True)[1][0][0]
n_f = sim.get_property('cSolv',True)[1][0][0]
n_fluid = n_f + n_s

rho = n_fluid/vol_sys
cs = n_s/vol_sys

vx_s = ufloat(sim.get_property('vx_Solu')[1][0][0],sim.get_property('vx_Solu')[1][0][1])
J_s = cs*vx_s
Q = ufloat(sim.get_property('vx_Sol',exact=True)[1][0][0],sim.get_property('vx_Sol',exact=True)[1][0][1])

exc_s = J_s- cs*Q

sim.add_property('Q',Q)
sim.add_property('J_s',J_s)
sim.add_property('J_s_exc',exc_s)


# =============================================================================
# Quantities in the entire Volume weighted by the bulk
# =============================================================================

exc_s_hyb = J_s- cs_bulk*Q
sim.add_property('J_s_exc_hyb',exc_s_hyb)


# =============================================================================
# Quantities difference of velocities
# =============================================================================

exc_s_vel = cs_bulk*(vx_s-Q)
exc_s_HY = J_s-(n_s_bulk/(n_s_bulk+n_f_bulk)*rho)*Q

sim.add_property('J_s_exc_vel',exc_s_vel)
sim.add_property('J_s_exc_HY',exc_s_HY)



# Printing results to write on excel
sim.get_property('J_s',text = True)
sim.get_property('Q',text = True)





