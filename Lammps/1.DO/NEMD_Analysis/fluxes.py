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
import Lammps.core_functions as cf
import Others.Statistics.FastAverager as stat
import Lammps.General.Log_Analysis.Thermo_Analyser as thermo
import matplotlib.pyplot as plt
import numpy as np
import copy
import re
from scipy import optimize
from Lammps.Pore.qsub import simulation_results as sr
import pickle as pickle
from uncertainties import ufloat,unumpy
import glob
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


cna = sr.check_n_analyse(root_pattern,directory_pattern)
cna.check_finished("vdata.dat")
cna.check_stat("statistics.dat")
cna.stat_analysis("vdata.dat")
cna.check_thermo("thermo.dat")
cna.thermo_analysis("log.lammps")

sim = sr.simulation("grad_mu",0.5)

properties=sr.read_properties(root_pattern,["statistics.dat","thermo.dat"])
sim.add_properties(properties)



# Getting the dimensions of the simulation box and the regions
box_volume, box_limits = lu.read_box_limits('log.lammps')

##The following two parameters are obtained from the knowledge of the bulk propertiesbox_volume=8000 
rho_bulk= 0.752375
cs_bulk=0.375332


n_solutes=sim.get_property('cSolu')[1][0][0]


n_solvents=sim.get_property('cSolv')[1][0][0]
n_total=n_solutes+n_solvents
vx_solu=ufloat(sim.get_property('vx_Solu')[1][0][0],sim.get_property('vx_Solu')[1][0][1])
J_s=n_solutes/box_volume*vx_solu
Q=ufloat(sim.get_property('vx_Sol',exact=True)[1][0][0],sim.get_property('vx_Sol',exact=True)[1][0][1])
exc_sol_flux=J_s-cs_bulk*Q
sim.add_property('Q',Q)
sim.add_property('J_s',J_s)
sim.add_property('J_s_exc',exc_sol_flux)