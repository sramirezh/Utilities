#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 17:35:55 2019

Gathers the flow from diffusio-osmotic simulations, both from Pressure and chemical potential simulations

@author: sr802
"""

import os
import sys
import glob
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
from simulation_results import *

try:
    from uncertainties import ufloat
except ImportError as err2:
    print err2

cwd = os.getcwd() #current working directory

def initialise_sim_bundles(root_pattern,directory_pattern,dictionary):
    
    #Needed parameters
    roots=glob.glob(root_pattern)
    mu=cf.extract_digits(roots, sort=False)
    
    
    bundles=[]

    for i,root in enumerate(roots):
        directories=glob.glob('%s/%s'%(root,directory_pattern))
        # =============================================================================
        # Checking if the simulation and the required files are on each folder
        # =============================================================================
        
        print "\nChecking if the simulations finished with vdata\n"
        dir_fin=filter_directories(directories,"vdata.dat")
        
        print "\nChecking if the simulations finished with statistics\n"
        dir_stat=filter_directories(dir_fin,"statistics.dat")
        
        print "\nChecking if the simulations finished with thermo\n"
        dir_thermo=filter_directories(dir_fin,"thermo.dat")
    
    
    
        # =============================================================================
        # Running the necessary analysis
        # =============================================================================
        
        #directories to run statistics
        dir_run_vdata=[x for x in dir_fin if x  not in dir_stat]
        run_stat_file(dir_run_vdata,"vdata.dat",0.3,"statistics.dat")
        
        
        
        #Directories to run thermo analysis
        dir_run_thermo=[x for x in dir_fin if x  not in dir_thermo]
        run_thermo_directories(dir_run_thermo,"log.lammps",0.3)
        
        
        
        # Creating the statistic summary (MAYBE Get rid of this)
        array=gather_statistics(dir_fin,'Time',root)
    
    
        #Building the simulations
        times=construct_simulations(directories)
    
        #Creating the directory for the plots
        
        directory="%s/plots/"%root
        if not os.path.exists(directory):
            os.mkdir(directory)
    
        
        #Creating the bundle
        bundles.append(simulation_bundle(times,"mu",mu[i],root,dictionary=dictionary))
    
        #Plot for all the properties
    #    for prop in bundles[-1].simulations[-1].property_names:
    #        print "\ncreating the plot of %s"%prop
    #        if prop!="time":
    #            bundles[-1].plot_property(prop)
        bundles[-1].plot_all_properties()
        
    return bundles


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
                sim_bundle.plot_property(prop,fit)
            
            else:
                sim_bundle.plot_property(prop)


        
# =============================================================================
# main 
# =============================================================================
plt.close('all')



# =============================================================================
# Chemical potential simulations
# =============================================================================

#root_pattern="mu_force*"
#directory_pattern='[0-9]*'
#
#
#dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}
#
#
#bundles_mu=initialise_sim_bundles(root_pattern,directory_pattern,dictionary)
#final_mu=simulation_bundle(bundles_mu,'mu',3,cwd,dictionary=dictionary)
#
#specific_plot_all(final_mu)



# =============================================================================
# Pressure simulations
# =============================================================================

root_pattern="p_force*"
directory_pattern='[0-9]*'


dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}


bundles_p=initialise_sim_bundles(root_pattern,directory_pattern,dictionary)
final_p=simulation_bundle(bundles_p,'p',3,cwd,dictionary=dictionary)

specific_plot_all(final_p,fit=False)