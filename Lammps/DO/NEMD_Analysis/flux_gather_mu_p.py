#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 17:35:55 2019

Gathers the flow from diffusio-osmotic simulations, both from chemical 
potentials of solutes and solvents
Uses the classes in simulation_results.py

Some details about the simulations:
    1. There are two volumes where the fluxes are measured in both grad mu and
    grad p simulations. rBulk (z between 10-25) and rSystem (z between 0-25)
    

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
import re
import Lammps.lammps_utilities as lu
import Lammps.DO.EMD.density_analysis as da

cwd = os.getcwd() #current working directory


def plot_properties(instance, x_name, y_name, ax, plot_fit = True):
    """
    TODO this could be part of the class simulation_bundle, actually I could 
    make the get_property to return either the average or the list
    This is a partly based on simulation_bundle.plot_property 
    
    Args:
        instance: instance of simulation bundle
        x_name: property in the x axis
        y_name: property in the y axis
        x_label: label of the x axis in latex format, eg: r'\nabla p'
        y_label: label of the y axis in latex format
        plot_name: name of the pdf file without extension
        ax: an axis could be passed
    
    Returns:
        ax:
    """
    
    y = [i.n for i in instance.get_property(y_name, exact = True)[1][0]]
    y_error=[i.s for i in instance.get_property(y_name, exact = True)[1][0]]

    x = [i.n for i in instance.get_property(x_name, exact = True)[1][0]]
    

    ax.errorbar(x, y, yerr=y_error, fmt='o')
    pinit = [1.0]
    out = optimize.leastsq(errfunc1, pinit, args=(x, y, y_error), full_output=1)
    pfinal = out[0] #fitting coefficients

    
    error = np.sqrt(out[1])
    logger.info("The transport coefficient \Gamma_{%s%s} is %.6f +/- %.6f"%(y_name, x_name, pfinal[0],error[0][0]))
    
    if plot_fit == True:
        x.sort()
        y_fit = fitfunc1(pfinal, x)
        ax.plot(x, y_fit, ls = '--', c = ax.lines[-1].get_color())
    
    plt.tight_layout()
    
    return ax

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
            print("\ncreating the plot of %s"%prop)
            
            if "vx" in prop: 
                
                sim_bundle.plot_property(prop,fit=fit)
            
            else:
                sim_bundle.plot_property(prop)


fitfunc1 = lambda p, x: p * x  #Fitting to a line that goes through the origin
errfunc1 = lambda p, x, y, err: (y - fitfunc1(p, x)) / err #To include the error in the least squares




# TODO This function could be included inside the class simulation_bundle
def build_bundle(root_pattern, directory_pattern, rho_bulk, cs_bulk, sim_type):
    """
    SPECIFIC FUNCTION FOR THIS PROBLEM
    
    
    Constructs a bundle of simulations based on a directory tree.
    suppose a directory tree as follows:
         |-4.Applying_force_0.025
         | |-1
         | |-2
         | |-3
         | |-4
         |-4.Applying_force_0.030
         | |-1
         | |-2
         | |-3
         | |-4
         
    And creates the specific parameters for this analysis
         
    Args: 
        root_pattern: in the example above "4.Applying_force_*"
        directory_pattern: the runs inside each root, in the example "[0-9]*"
        rho_bulk: Density of the fluid in the bulk
        cs_bulk: Density on solutes in the bulk
        sim_type: identifier of the simulations in the bundle, could be "p", "mu_s", etc
        
    
    Returns:
        final_mu: an instance of sr.simulation_bundle containing all the 
                    information from the simulation bundle
        
    """

    dictionary = {'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}
# =============================================================================
#     Create or read the bundle
# =============================================================================
    if not glob.glob("%s"%sim_type):
        logger.info("\nAnalysing the mu grad simulations\n")
        bundles_mu = sr.initialise_sim_bundles(root_pattern,sim_type,directory_pattern,dictionary, plot = False)
        final_mu = sr.simulation_bundle(bundles_mu, "bundle" , 3,cwd,dictionary = dictionary, ave = False)
        final_mu.save("%s"%sim_type)
    else:
        final_mu = cf.load_instance("%s"%sim_type) 
        
# =============================================================================
#     Obtaining some geometry parameters
# =============================================================================
    #         TODO improve this 
    root = glob.glob(root_pattern)[0]
    folder = glob.glob("%s/%s"%(root,directory_pattern))[0]
    
    # Getting the dimensions of the simulation box and the regions
    box_volume, box_limits = lu.read_box_limits('%s/log.lammps'%folder)
    
    bulk_heigth  = lu.read_region_height('rBulk','%s/in.geom'%folder)
    
    #This is a region defined as the bulk+interface
    sys_heigth  = lu.read_region_height('rSystem','%s/in.geom'%folder) 
    
    vol_bulk = (box_limits[1][0]-box_limits[0][0])*(box_limits[1][1]-box_limits[0][1])*bulk_heigth[0]
    vol_sys = (box_limits[1][0]-box_limits[0][0])*(box_limits[1][1]-box_limits[0][1])*sys_heigth[0]

    # =============================================================================
    # Calculations of the transport coefficients from DP problem
    # =============================================================================
    final_mu.average = False
    count = 0
    for bund in final_mu.simulations:
        
        #This is just to compute the excess of solute flux
        
        for sim in bund.simulations:
            
# =============================================================================
#             # Quantities in the bulk 
# =============================================================================
            n_s_bulk = sim.get_property('cBSolu',True)[1][0][0]
            n_f_bulk = sim.get_property('cBSolv',True)[1][0][0]
            n_fluid_bulk = n_f_bulk + n_s_bulk
            rho_bulk = n_fluid_bulk/vol_bulk
            cs_bulk = n_s_bulk/vol_bulk
            
            vx_s_bulk = ufloat(sim.get_property('vxB_Solu')[1][0][0],sim.get_property('vxB_Solu')[1][0][1])
            J_s_bulk = cs_bulk*vx_s_bulk
            Q_bulk = ufloat(sim.get_property('vxB_Sol',exact=True)[1][0][0],sim.get_property('vxB_Sol',exact=True)[1][0][1])
            
            exc_s_bulk = J_s_bulk- cs_bulk * Q_bulk
            
            sim.add_property('Q_bulk',Q_bulk)
            sim.add_property('J_s_bulk',J_s_bulk)
            sim.add_property('J_s_exc_bulk',exc_s_bulk)
            
# =============================================================================
#             # Quantities in the entire system
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
            
            
        # To add all the properties that we just added in the sub bundle
            bund.update_properties()


            count+=1
        

    final_mu.update_properties()
    
    
    return final_mu

# =============================================================================
#     Main
# =============================================================================

plot_dir = "plots/3.fluxes"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   

ms_path = '4.Applying_force_[0-9]*'
mp_path = '5.Applying_force_p_[0-9]*'
ms_dir = '[0-9]*'
mp_dir = '[0-9]*'
# =============================================================================
# Assumptions and external parameters
# =============================================================================

# Obtaining some basic properties from the bulk, assuming that the NEDM
# run does not alter the densities distribution

dir_measure = "3.Measuring"
fluid = da.DensityDistribution("properties_short.dat", "rBulk", directory = dir_measure) 
solute = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = dir_measure) 

rho_bulk = fluid.rho_bulk
cs_bulk = solute.rho_bulk

logger.info("Using the following parameters:\n rho_bulk = %f\n cs_bulk = %f\n"%(rho_bulk, cs_bulk ))

# =============================================================================
#     Loading all the results for different gradients
# =============================================================================
final_mus = build_bundle(ms_path, ms_dir, rho_bulk, cs_bulk,"grad_mu_s")
final_p = build_bundle(mp_path, mp_dir, rho_bulk, cs_bulk,"grad_p")


# =============================================================================
# Getting the grad s primed
# =============================================================================
for bund in final_mus.simulations:
        
        bund.add_property('grad_mu_primed',rho_bulk/(rho_bulk-cs_bulk)*float(bund.get_property('grad_mu_s',exact=True)[1][0]))
        bund.update_properties() # To update the bundle

final_mus.update_properties()


# =============================================================================
#  Plotting
# =============================================================================
cf.set_plot_appearance()
plt.close("all")

# Q vs grad mu primed
fig1, ax1 = plt.subplots()
plot_properties(final_mus, 'grad_mu_primed','Q',ax1) 
ax1.set_ylabel(r'$Q$')
ax1.set_xlabel(r"$-\nabla \mu'_{s}$")
#ax.legend(loc = 'upper right')
#ax.set_ylim(0, None)
#ax.set_xlim(0, 8)
fig1.tight_layout()
fig1.savefig('%s/q_vs_grad_mu.pdf'%(plot_dir), Transparent = True)
                
# Js vs grad P
fig2, ax2 = plt.subplots()
plot_properties(final_p, 'grad_p','J_s_exc_HY', ax2)
ax2.set_ylabel(r'$J_s-c^*_sQ$')
ax2.set_xlabel(r'$-\nabla P$')
#ax.legend(loc = 'upper right')
ax2.set_ylim(0, None)
ax2.set_xlim(0, None)
fig2.tight_layout()
fig2.savefig('%s/Js_exc_HY_vs_grad_P.pdf'%(plot_dir), Transparent = True)

# Checking the different exess fluxes

fig3, ax3 = plt.subplots()
plot_properties(final_p, 'grad_p', 'J_s_exc_bulk', ax3)
plot_properties(final_p, 'grad_p', 'J_s_exc', ax3)
plot_properties(final_p, 'grad_p', 'J_s_exc_hyb', ax3)
plot_properties(final_p, 'grad_p', 'J_s_exc_vel', ax3)
plot_properties(final_p, 'grad_p', 'J_s_exc_HY', ax3)
ax3.set_ylabel(r'$J_s$')
ax3.set_xlabel(r'$-\nabla P$')
ax3.legend([r"$J'^B_s$", r"$J_s'$", r"$J_s'^{hyb}$", 
            r"$J_s'^{vel}$", r"$J_s'^{HY}$"], loc = 'upper left', ncol = 2)
#ax.set_ylim(0, None)
#ax.set_xlim(0, 8)
fig3.tight_layout()
fig3.savefig('%s/Js_exc_comparison.pdf'%(plot_dir), Transparent = True)

 
    



