#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 17:35:55 2019

Gathers the flow from simulations where the  
chemical potential on each one of the species

@author: sr802
"""

import os
import sys
Utilities_path=os.path.join(os.path.dirname(__file__), '../../../')
sys.path.append(Utilities_path) #This falls into Utilities path
import Lammps.core_functions as cf
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import Lammps.Pore.qsub.simulation_results as sr
from uncertainties import ufloat,unumpy
import glob
import argparse
import re

cwd = os.getcwd() #current working directory


# def specific_plot_all(sim_bundle,fit=True):
#     """
#     Very specific function
#     Plots all the properties but only fits to a line the velocities
#     Args:
#         simulation_bundle
#     """
    
#     #Copy of plot_all_properties
#     for i,prop in enumerate(sim_bundle.simulations[-1].property_names):
        
#         if prop!="time" and i>0:
#             print("\ncreating the plot of %s"%prop)
            
#             if "vx" in prop: 
                
#                 sim_bundle.plot_property(prop,fit=fit)
            
#             else:
#                 sim_bundle.plot_property(prop)


# fitfunc1 = lambda p, x: p * x  #Fitting to a line that goes through the origin
# errfunc1 = lambda p, x, y, err: (y - fitfunc1(p, x)) / err #To include the error in the least squares





def build_bundle(root_pattern, directory_pattern, box_volume, rho_bulk, cs_bulk, sim_type):
    """
    SPECIFIC FUNCTION FOR THIS PROBLEM
    
    
    Constructs a bundle of simulations based on a directory tree.
    suppose a directory tree as follows:
         |-mus_force_*
         | |-202000
         | |-203000
         | |-...
         |-muf_force_*
         | |-202000
         | |-203000
         | |-...
         
    And creates the specific parameters for this analysis
         
    Args: 
        root_pattern: in the example above "4.Applying_force_*"
        directory_pattern: the runs inside each root, in the example "[0-9]*"
        rho_bulk: Density of the fluid in the bulk
        cs_bulk: Density on solutes in the bulk
        sim_type: identifier of the simulations in the bundle, could be "s", "f", etc
        
    
    Returns:
        final_mu: an instance of sr.simulation_bundle containing all the 
                    information from the simulation bundle
        
    """

    dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}

    #If the object was not saved
    if not glob.glob("mu%s.pkl"%sim_type):
        bundles_mu = sr.initialise_sim_bundles(root_pattern,'mu',directory_pattern,dictionary)
        final_mu = sr.simulation_bundle(bundles_mu,'mu_bundle',3,cwd,dictionary = dictionary, ave = False)
        final_mu.save("mu%s"%sim_type)

    # =============================================================================
    # Calculations of the transport coefficients from DP problem
    # =============================================================================

    final_mu = cf.load_instance("mu%s.pkl"%sim_type) 
    
    final_mu.average = False
    count = 0
    for bund in final_mu.simulations:
        
        #This is just to compute the excess of solute flux
        
        for sim in bund.simulations:
            n_solutes = sim.get_property('cSolu')[1][0][0]
            n_solvents = sim.get_property('cSolv')[1][0][0]
            c_total = (n_solutes + n_solvents) / box_volume
            c_solvents = n_solvents / box_volume
            c_solutes = n_solutes / box_volume
            
            vx_solu = ufloat(sim.get_property('vx_Solu')[1][0][0],sim.get_property('vx_Solu')[1][0][1])
            vx_solv = ufloat(sim.get_property('vx_Solv')[1][0][0],sim.get_property('vx_Solv')[1][0][1])
            
            J_s = c_solutes * vx_solu
            J_f = c_solvents * vx_solv
            Q = ufloat(sim.get_property('vx_Sol',exact=True)[1][0][0],sim.get_property('vx_Sol',exact=True)[1][0][1])
            pref = c_total * cs_bulk/rho_bulk
            exc_sol_flux = J_s - pref * Q
            
            
            #Adding specific properties to the individual simulations
            
            sim.add_property('Q',Q)
            sim.add_property('Js',J_s)
            sim.add_property('Jf',J_f)
            sim.add_property('Js_exc',exc_sol_flux)

            count+=1
        
        #Adding specific properties to the secondary bundle
        
        bund.add_property('grad_mu',float(bund.get_property('mu',exact=True)[1][0]))
        bund.update_properties() # To update the bundle


    final_mu.update_properties()
    
    
    return final_mu

# =============================================================================
#     Main
# =============================================================================

plot_dir = "plots/3.fluxes"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

logger = cf.log(__file__, os.getcwd(),plot_dir)   

ms_path = 'mus_force_*'
mf_path = 'muf_force_*'
ms_dir = '[0-9]*'
mf_dir = '[0-9]*'


logger.info("Analysing %s and %s"%(ms_path, mf_path))

# =============================================================================
# Assumptions and external parameters
# =============================================================================
box_volume = 20**3
rho_bulk =  0.752375
cs_bulk = 0.375332

logger.info("I am using the following parameters:\n box volume = %f\n rho_bulk = %f\n cs_bulk = %f\n"%(box_volume, rho_bulk, cs_bulk ))

# =============================================================================
#     Loading all the results for different gradients
# =============================================================================

final_mus = build_bundle(ms_path, ms_dir, box_volume,rho_bulk,cs_bulk,"s")
final_muf = build_bundle(mf_path, mf_dir, box_volume,rho_bulk,cs_bulk,"f")

final_mus.logger = logger
final_muf.logger = logger

# =============================================================================
#  Plotting
# =============================================================================
cf.set_plot_appearance()
plt.close("all")


# plot_properties(final_mus, 'grad_mu','Js', x_label = r'$-\nabla \mu_s$', y_label = r'$J_s$', plot_name ="ss" )
#     plot_properties(final_mus, 'grad_mu','Jf', x_label = r'$-\nabla \mu_s$',  y_label = r'$J_f$', plot_name ="fs" )
#     plot_properties(final_muf, 'grad_mu','Js',x_label = r'$-\nabla \mu_f$', y_label = r'$J_s$' ,plot_name ="sf" )
#     plot_properties(final_muf, 'grad_mu','Jf', x_label = r'$-\nabla \mu_f$', y_label = r'$J_f$', plot_name ="ff"  )

# Js vs grad mu_s
fig1, ax1 = plt.subplots()
final_mus.plot_property( ax1,'Js','grad_mu', fit = True) 
ax1.set_ylabel(r'$J_s$')
ax1.set_xlabel(r"$-\nabla \mu_{s}$")
#ax.legend(loc = 'upper right')
#ax.set_ylim(0, None)
#ax.set_xlim(0, 8)
fig1.tight_layout()
fig1.savefig('%s/Js_vs_grad_mu_s.pdf'%(plot_dir), Transparent = True)


# Jf vs grad mu_s
fig1, ax1 = plt.subplots()
final_mus.plot_property( ax1,'Jf','grad_mu', fit = True) 
ax1.set_ylabel(r'$J_s$')
ax1.set_xlabel(r"$-\nabla \mu_{s}$")
#ax.legend(loc = 'upper right')
#ax.set_ylim(0, None)
#ax.set_xlim(0, 8)
fig1.tight_layout()
fig1.savefig('%s/Jf_vs_grad_mu_s.pdf'%(plot_dir), Transparent = True)


# vs vs grad mu_s
fig1, ax1 = plt.subplots()
final_muf.plot_property( ax1,'Js','grad_mu', fit = True) 
ax1.set_ylabel(r'$J_s$')
ax1.set_xlabel(r"$-\nabla \mu_{f}$")
#ax.legend(loc = 'upper right')
#ax.set_ylim(0, None)
#ax.set_xlim(0, 8)
fig1.tight_layout()
fig1.savefig('%s/Js_vs_grad_mu_f.pdf'%(plot_dir), Transparent = True)

# vs vs grad mu_s
fig1, ax1 = plt.subplots()
final_muf.plot_property( ax1,'Jf','grad_mu', fit = True) 
ax1.set_ylabel(r'$J_f$')
ax1.set_xlabel(r"$-\nabla \mu_{f}$")
#ax.legend(loc = 'upper right')
#ax.set_ylim(0, None)
#ax.set_xlim(0, 8)
fig1.tight_layout()
fig1.savefig('%s/Jf_vs_grad_mu_f.pdf'%(plot_dir), Transparent = True)


# # TODO change group pattern to something more flexible as just file_names
# if __name__ == "__main__":
#     """
#     THIS IS VERY SPECIFIC
#     The arguments of this depend on the application
#     """
#     parser = argparse.ArgumentParser(description='Launch simulations from restart',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#     parser.add_argument('-ms_pat', metavar='ms_pat',help='Generic name of the directories with mu_s gradient',default='mus_force_*')
#     parser.add_argument('-mf_pat', metavar='mf_pat',help='Generic name of the directories with mu_f gradient',default='muf_force_*')
#     parser.add_argument('-ms_dir', metavar='ms_dir',help='Patter of the files inside, in this case restart are like 202000',default='[0-9]*')
#     parser.add_argument('-mf_dir', metavar='mf_dir',help='Patter of the files inside, in this case restart are like 202000',default='[0-9]*')
#     args = parser.parse_args()

    
#     main(args.ms_pat,args.mf_pat,args.ms_dir,args.mf_dir)
    
    
    

# # =============================================================================
# # Peclet number computation
# # =============================================================================

# D = 0.066 # Diffusion coefficient
# L = 2
# pe_p = []
# for bund in final_p.simulations:
#     pe_p.append(bund.get_property('Q')[1][0][0]/D * L)
#     pe_p.sort()
    
    
    
    
# # =============================================================================
# # Peclet number computation
# # =============================================================================

# pe_mu = []
# for bund in final_mu.simulations:
#     pe_mu.append(bund.get_property('Q')[1][0][0]/D * L)
#     pe_mu.sort()