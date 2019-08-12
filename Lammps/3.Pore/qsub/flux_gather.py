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
import cPickle as pickle
from uncertainties import ufloat,unumpy

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

def load_structure(file_name):
    """
    Loads the data structure to be used later
    """
    file1 = open(file_name, 'rb')
    instance = pickle.load(file1)
    file1.close()

    return instance

### =============================================================================
### Chemical potential simulations
### =============================================================================
##
#root_pattern="mu_force*"
#directory_pattern='[0-9]*'
#parameter_id='mu'
#
#dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}
#
#
#bundles_mu=sr.initialise_sim_bundles(root_pattern,parameter_id,directory_pattern,dictionary)
#final_mu=sr.simulation_bundle(bundles_mu,parameter_id,3,cwd,dictionary=dictionary)
#
#specific_plot_all(final_mu)
#
#
#
## =============================================================================
## Pressure simulations
## =============================================================================
#

#root_pattern="p_force*"
#directory_pattern='[0-9]*'
#parameter_id='p'
#
#dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}
#
## This function could be included inside the class simulation_bundle
#bundles_p=sr.initialise_sim_bundles(root_pattern,parameter_id,directory_pattern,dictionary)
#final_p=sr.simulation_bundle(bundles_p,parameter_id,3,cwd,dictionary=dictionary)
#
#
#specific_plot_all(final_p)

#writing the structure
#afile = open(r'p.pkl', 'wb')
#pickle.dump(final_p, afile)
#afile.close()




# =============================================================================
# Calculations for the pressure driven to get the excess solute flux
# =============================================================================
    
fitfunc1 = lambda p, x: p * x  #Fitting to a line that goes through the origin
errfunc1 = lambda p, x, y, err: (y - fitfunc1(p, x)) / err #To include the error in the least squares

final_p=load_structure("p.pkl")
#The following two parameters are obtained from the knowledge of the bulk properties
box_volume=8000 
rho_bulk=0.752375
cs_bulk=0.375332

f_p=[]
exc_solute=[]
for bund in final_p.simulations:
    
    #Getting the applied forces
    f_p.extend(bund.get_property('p',exact=True)[1])
    print f_p[-1]
    #Getting the solute excess
    exc_sol_array=[]
    for sim in bund.simulations:

        n_solutes=sim.get_property('cSolu')[1][0][0]
        vx_solu=ufloat(sim.get_property('vx_Solu')[1][0][0],sim.get_property('vx_Solu')[1][0][1])
        J_s=n_solutes/box_volume*vx_solu
        Q=ufloat(sim.get_property('vx_Sol',exact=True)[1][0][0],sim.get_property('vx_Sol',exact=True)[1][0][1])
        exc_sol_flux=J_s-cs_bulk*Q
        
        exc_sol_array.append(exc_sol_flux)
    
    exc_solute.append(sum(exc_sol_array)/len(exc_sol_array))
    
grad_p=rho_bulk*np.array(f_p)
y=[i.n for i in exc_solute]
y_error=[i.s for i in exc_solute]









cf.set_plot_appearance()

fig,ax=plt.subplots()

plt.errorbar(grad_p,y,yerr=y_error,xerr=None,fmt='o')



pinit=[1.0]
out = optimize.leastsq(errfunc1, pinit, args=(grad_p, y, y_error), full_output=1)
pfinal = out[0] #fitting coefficients
grad_p=np.insert(grad_p,0,0)
ax.plot(np.unique(grad_p),fitfunc1(pfinal,np.unique(grad_p)),linestyle='--')



ax.set_xlabel(r'$\nabla P$')
ax.set_ylabel(r'$J_s-c_s^BQ$')
plt.tight_layout()