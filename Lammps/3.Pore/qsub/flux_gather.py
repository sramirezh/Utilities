#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 17:35:55 2019

Gathers the flow from diffusio-osmotic simulations, both from Pressure and chemical potential simulations
Uses the classes in simulation_results.py

@author: sr802
"""

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
import simulation_results as sr

try:
    from uncertainties import ufloat
except ImportError as err2:
    print err2

cwd = os.getcwd() #current working directory


def specific_plot_all(sim_bundle,fit=True):
    """
    Very specific function
    Plots all the properties but only fits to a line the velocities
    Args:
        simulation_bundle
    """
    
    #Copy of plot_all_properties
    for i,prop in enumerate(sim_bundle.simulations[-1].property_names):
        
        if prop!="time" and i>0:
            print "\ncreating the plot of %s"%prop
            
            if "vx" in prop: 
                
                sim_bundle.plot_property(prop,fit=fit)
            
            else:
                sim_bundle.plot_property(prop)


        
# =============================================================================
# main 
# =============================================================================
plt.close('all')



## =============================================================================
## Chemical potential simulations
## =============================================================================
#
root_pattern="mu_force*"
directory_pattern='[0-9]*'
parameter_id='mu'

dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}


bundles_mu=sr.initialise_sim_bundles(root_pattern,parameter_id,directory_pattern,dictionary)
final_mu=sr.simulation_bundle(bundles_mu,parameter_id,3,cwd,dictionary=dictionary)

specific_plot_all(final_mu)



# =============================================================================
# Pressure simulations
# =============================================================================

root_pattern="p_force*"
directory_pattern='[0-9]*'
parameter_id='p'

dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}

# This function could be included inside the class simulation_bundle
bundles_p=sr.initialise_sim_bundles(root_pattern,parameter_id,directory_pattern,dictionary)
final_p=sr.simulation_bundle(bundles_p,parameter_id,3,cwd,dictionary=dictionary)

specific_plot_all(final_p,fit=True)