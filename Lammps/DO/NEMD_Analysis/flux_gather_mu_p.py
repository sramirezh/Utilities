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
import argparse
import re
import Lammps.lammps_utilities as lu
import Lammps.DO.EMD.density_analysis as da

cwd = os.getcwd() #current working directory


def plot_properties(instance, x_name, y_name, x_label = None, y_label = None, 
                    plot_name = None, plot_fit = True):
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
    
    Returns:
        fig:
        ax:
    """
    
    cf.set_plot_appearance()
    
    y = [i.n for i in instance.get_property(y_name, exact = True)[1][0]]
    y_error=[i.s for i in instance.get_property(y_name, exact = True)[1][0]]

    x = [i.n for i in instance.get_property(x_name, exact = True)[1][0]]
    
    fig, ax = plt.subplots()

    plt.errorbar(x, y,yerr=y_error,xerr=None,fmt='o')
    pinit=[1.0]
    out = optimize.leastsq(errfunc1, pinit, args=(x, y, y_error), full_output=1)
    pfinal = out[0] #fitting coefficients

    
    error = np.sqrt(out[1]) 
    print("The transport coefficient \Gamma_{%s%s} is %.6f +/- %.6f"%(y_name,x_name, pfinal[0],error[0][0]))
    
    if plot_fit == True:
        x.sort()
        y_fit = fitfunc1(pfinal, x)
        ax.plot(x, y_fit, ls = '--', c = ax.lines[-1].get_color())
    
    #Checking if there are axis names
    if x_label == None:
        x_label = re.sub('_',' ', x_name)
        ax.set_xlabel(x_label) #The zeroth-property is the param_id
    else:
        ax.set_xlabel(x_label)
        
    if y_label == None:

        if y_name in list(instance.dictionary.keys()):
            y_label = instance.dictionary[y_name]
        else:
            y_label = re.sub('_',' ',y_name)
        ax.set_ylabel(y_label)
    else:
        ax.set_ylabel(y_label)
     
    plt.tight_layout()
    
    fig.savefig("%s.pdf"%plot_name, transparent=True)
    cf.save_instance(ax,"%s"%plot_name)
    
    return fig, ax

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
def build_bundle_mu(root_pattern, directory_pattern, rho_bulk, cs_bulk, sim_type):
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
        sim_type: identifier of the simulations, could be "p", "mu_s", etc
    
    Returns:
        final_mu: an instance of sr.simulation_bundle containing all the 
                    information from the simulation bundle
        
    """

    dictionary = {'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}
# =============================================================================
#     Create or read the bundle
# =============================================================================
    if not glob.glob("%s"%sim_type):
        print("\nAnalysing the mu grad simulations\n")
        bundles_mu = sr.initialise_sim_bundles(root_pattern,'mu',directory_pattern,dictionary, plot = False)
        final_mu = sr.simulation_bundle(bundles_mu,sim_type,3,cwd,dictionary = dictionary, ave = False)
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
            n_s = sim.get_property('cSolu')[1][0][0]
            n_f = sim.get_property('cSolv')[1][0][0]
            n_s_bulk = sim.get_property('cBSolu',True)[1][0][0]
            n_f_bulk = sim.get_property('cBSolv',True)[1][0][0]
            c_total =  (n_s+n_f)/vol_sys
#            c_solvents = n_solvents/sim.volume
            c_s = n_s / vol_sys
#            
            vx_solu = ufloat(sim.get_property('vx_Solu')[1][0][0],sim.get_property('vx_Solu')[1][0][1])
#            vx_solv = ufloat(sim.get_property('vx_Solv')[1][0][0],sim.get_property('vx_Solv')[1][0][1])
            
            J_s = c_s * vx_solu
            Q = ufloat(sim.get_property('vx_Sol',exact=True)[1][0][0],sim.get_property('vx_Sol',exact=True)[1][0][1])
            pref = (n_s_bulk/(n_s_bulk+n_f_bulk)*c_total)
            exc_sol_flux = J_s - pref * Q
            
            
            #Adding specific properties to the individual simulations
            
            sim.add_property('Q', Q)
            sim.add_property('Js',J_s)
            sim.add_property('Js_exc',exc_sol_flux)

            count+=1
        
        #Adding specific properties to the secondary bundle
        
        bund.add_property('grad_mu',rho_bulk/(rho_bulk-cs_bulk)*float(bund.get_property('mu',exact=True)[1][0]))
        bund.update_properties() # To update the bundle


    final_mu.update_properties()
    
    
    return final_mu



def build_bundle_p(root_pattern, directory_pattern, rho_bulk, cs_bulk, sim_type):
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
        sim_type: identifier of the simulations, could be "p", "mu_s", etc
    
    Returns:
        final_mu: an instance of sr.simulation_bundle containing all the 
                    information from the simulation bundle
        
    """

    dictionary = {'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}
    
# =============================================================================
#     Create or read the bundle
# =============================================================================
    if not glob.glob("%s"%sim_type):
        print("\nAnalysing the mu grad simulations\n")
        bundles_p = sr.initialise_sim_bundles(root_pattern,'p',directory_pattern,dictionary, plot = False)
        final_p = sr.simulation_bundle(bundles_p,sim_type,3,cwd,dictionary = dictionary, ave = False)
        final_p.save("%s"%sim_type)
    else:
        final_p = cf.load_instance("%s"%sim_type) 
        
        
        
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
    final_p.average = False
    count = 0
    for bund in final_p.simulations:
        
        #This is just to compute the excess of solute flux
        
        for sim in bund.simulations:
            n_s = sim.get_property('cSolu')[1][0][0]
            n_f = sim.get_property('cSolv')[1][0][0]
            n_s_bulk = sim.get_property('cBSolu',True)[1][0][0]
            n_f_bulk = sim.get_property('cBSolv',True)[1][0][0]
            c_total =  (n_s+n_f)/vol_sys
#            c_solvents = n_solvents/sim.volume
            c_s = n_s / vol_sys
#            
            vx_solu = ufloat(sim.get_property('vx_Solu')[1][0][0],sim.get_property('vx_Solu')[1][0][1])
#            vx_solv = ufloat(sim.get_property('vx_Solv')[1][0][0],sim.get_property('vx_Solv')[1][0][1])
            
            J_s = c_s * vx_solu
            Q = ufloat(sim.get_property('vx_Sol',exact=True)[1][0][0],sim.get_property('vx_Sol',exact=True)[1][0][1])
            pref = (n_s_bulk/(n_s_bulk+n_f_bulk)*c_total)
            exc_sol_flux = J_s - pref * Q
            
            
            #Adding specific properties to the individual simulations
            
            sim.add_property('Q', Q)
            sim.add_property('Js',J_s)
            sim.add_property('Js_exc',exc_sol_flux)

            count+=1
        
        # It has to be updated as new properties were added to the simulations
        bund.update_properties()
        
#        bund.add_property('grad_mu',rho_bulk/(rho_bulk-cs_bulk)*float(bund.get_property('mu',exact=True)[1][0]))
#        .update_properties() # To update the bundle

        # It has to be updated as new properties were added to the subbundle
        final_p.update_properties()
    
    
    return final_p

# =============================================================================
#     Main
# =============================================================================

#def main(ms_pat, mp_pat, ms_dir, mp_dir):

logger = cf.log(__file__, os.getcwd())

ms_path = '4.Applying_force_[0-9]*'
mp_path = '5.Applying_force_p_[0-9]*'
ms_dir = '[0-9]*'
mp_dir = '[0-9]*'
# =============================================================================
# Assumptions and external parameters
# =============================================================================

#Obtaining some basic properties from the bulk, assuming that the NEDM
#Run does not alter the densities distribution

dir_measure = "3.Measuring"
fluid = da.DensityDistribution("properties_short.dat", "rBulk", directory = dir_measure) 
solute = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = dir_measure) 

rho_bulk = fluid.rho_bulk
cs_bulk = solute.rho_bulk

logger.info("Using the following parameters:\n rho_bulk = %f\n cs_bulk = %f\n"%(rho_bulk, cs_bulk ))

# =============================================================================
#     Loading all the results for different gradients
# =============================================================================
final_mus = build_bundle_mu(ms_path, ms_dir, rho_bulk, cs_bulk,"grad_mu_s")
final_p = build_bundle_p(mp_path, mp_dir, rho_bulk, cs_bulk,"grad_p")


plot_properties(final_mus, 'grad_mu','Q', x_label = r"$-\nabla \mu'_{s}$", 
                y_label = r'$Q$', plot_name ="q_vs_grad_mu" )
plot_properties(final_p, 'p','Js_exc', x_label = r'$-\nabla P$', 
                y_label = r'$J_s-c^*_sQ$', plot_name ="Js_exc_vs_grad_P"  )
    


## TODO change group pattern to something more flexible as just file_names
#if __name__ == "__main__":
#    """
#    THIS IS VERY SPECIFIC
#    The arguments of this depend on the application
#    """
#    parser = argparse.ArgumentParser(description='Launch simulations from restart',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#    parser.add_argument('-ms_pat', metavar='ms_pat',help='Generic name of the directories with mu_s gradient',default='4.Applying_force_[0-9]*')
#    parser.add_argument('-mp_pat', metavar='mp_pat',help='Generic name of the directories with p gradient',default='5.Applying_force_p_[0-9]*')
#    parser.add_argument('-ms_dir', metavar='ms_dir',help='Patter of the files inside, in this case restart are like 202000',default='[0-9]*')
#    parser.add_argument('-mp_dir', metavar='mp_dir',help='Patter of the files inside, in this case restart are like 202000',default='[0-9]*')
#    args = parser.parse_args()
#    
#    logger = cf.log(__file__, os.getcwd()) 
#    print = logger.info
#    
#    main(args.ms_pat,args.mp_pat,args.ms_dir,args.mp_dir)
#    