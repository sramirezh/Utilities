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
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import simulation_results as sr
from uncertainties import ufloat,unumpy
import glob
import argparse
import re

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
            print("\ncreating the plot of %s"%prop)
            
            if "vx" in prop: 
                
                sim_bundle.plot_property(prop,fit=fit)
            
            else:
                sim_bundle.plot_property(prop)


fitfunc1 = lambda p, x: p * x  #Fitting to a line that goes through the origin
errfunc1 = lambda p, x, y, err: (y - fitfunc1(p, x)) / err #To include the error in the least squares





def mu_simulations(root_pattern, directory_pattern, box_volume,rho_bulk,cs_bulk,T_qs, T_ss):

    dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}

    # TODO This function could be included inside the class simulation_bundle
    #If the object was not saved
    if not glob.glob("mu.pkl"):
        print("\nAnalysing the mu grad simulations\n")
        bundles_mu = sr.initialise_sim_bundles(root_pattern,'mu',directory_pattern,dictionary)
        final_mu = sr.simulation_bundle(bundles_mu,'mu_bundle',3,cwd,dictionary = dictionary, ave = False)
        final_mu.save("mu")

    # =============================================================================
    # Calculations of the transport coefficients from DP problem
    # =============================================================================

    final_mu=cf.load_instance("mu.pkl") 
    
    final_mu.average = False
    count = 0
    for bund in final_mu.simulations:
        
        #This is just to compute the excess of solute flux
        
        for sim in bund.simulations:
            n_solutes=sim.get_property('cSolu')[1][0][0]
            n_solvents=sim.get_property('cSolv')[1][0][0]
            c_total = (n_solutes+n_solvents)/box_volume
            c_solutes = n_solutes/box_volume
            
            vx_solu=ufloat(sim.get_property('vx_Solu')[1][0][0],sim.get_property('vx_Solu')[1][0][1])
            J_s=c_solutes*vx_solu
            Q=ufloat(sim.get_property('vx_Sol',exact=True)[1][0][0],sim.get_property('vx_Sol',exact=True)[1][0][1])
            pref = c_total*cs_bulk/rho_bulk
            exc_sol_flux = J_s - pref * Q
            
            
            #Adding specific properties to the individual simulations
            
            sim.add_property('Q',Q)
            sim.add_property('Js',J_s)
            sim.add_property('Js_exc',exc_sol_flux)

            count+=1
        
        #Adding specific properties to the secondary bundle
        
        bund.add_property('grad_mu',rho_bulk/(rho_bulk-cs_bulk)*float(bund.get_property('mu',exact=True)[1][0]))
        bund.add_upd_properties() # To update the bundle


    final_mu.add_upd_properties()
    
    
    return final_mu


def p_simulations(root_pattern, directory_pattern, box_volume,rho_bulk,cs_bulk,T_sq,T_qq):



    dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}

    if not glob.glob("p.pkl"):

        bundles_p = sr.initialise_sim_bundles(root_pattern,'p',directory_pattern,dictionary)
        final_p = sr.simulation_bundle(bundles_p,'p_bundle',3,cwd,dictionary = dictionary, ave = False)
        final_p.save("p")


    # =============================================================================
    # Calculations for the pressure driven flow
    # =============================================================================

    final_p=cf.load_instance("p.pkl") 

    for bund in final_p.simulations:
        

        for sim in bund.simulations:

            n_solutes=sim.get_property('cSolu')[1][0][0]
            n_solvents=sim.get_property('cSolv')[1][0][0]
            c_total = (n_solutes+n_solvents)/box_volume
            c_solutes = n_solutes/box_volume
            
            vx_solu=ufloat(sim.get_property('vx_Solu')[1][0][0],sim.get_property('vx_Solu')[1][0][1])
            J_s=c_solutes*vx_solu
            Q=ufloat(sim.get_property('vx_Sol',exact=True)[1][0][0],sim.get_property('vx_Sol',exact=True)[1][0][1])
            pref = c_total*cs_bulk/rho_bulk
            exc_sol_flux = J_s - pref * Q
            
            sim.add_property('Q',Q)
            sim.add_property('Js',J_s)
            sim.add_property('Js_exc',exc_sol_flux)
                
        #Adding specific properties to the secondary bundle
        bund.add_property("grad_p", bund.get_property('p',exact=True)[1][0])
        
        bund.add_upd_properties() # To update the bundle
        
            
    final_p.add_upd_properties()
    
    
    return final_p


def plot_properties(instance, x_name, y_name, x_label = None, y_label = None, plot_name=None):
    """
    TODO this could be part of the class simulation_bundle, actually I could make the get_property to return either the average or the list
    This is a partly based on simulation_bundle.plot_property that 
    """
    y=[i.n for i in instance.get_property(y_name, exact = True)[1][0]]
    y_error=[i.s for i in instance.get_property(y_name, exact = True)[1][0]]

    x = [i.n for i in instance.get_property(x_name, exact = True)[1][0]]


    fig,ax=plt.subplots()

    plt.errorbar(x,y,yerr=y_error,xerr=None,fmt='o')



    pinit=[1.0]
    out = optimize.leastsq(errfunc1, pinit, args=(x, y, y_error), full_output=1)
    pfinal = out[0] #fitting coefficients
    error = np.sqrt(out[1]) 
    print("The transport coefficient \Gamma_{%s%s} is %.6f +/- %.6f"%(y_name,x_name, pfinal[0],error[0][0]))
        
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




def main(m_pat, p_pat, m_dir, p_dir):


    # =============================================================================
    # Assumptions and external parameters
    # =============================================================================


    ##The following two parameters are obtained from the knowledge of the bulk properties
    # TODO this can be obtained with lammps_utilities.py, add this property to class

    box_volume = 20**3
    rho_bulk =  0.752375
    cs_bulk = 0.375332

    print("I am using the following parameters:\n box volume = %f\n rho_bulk = %f\n cs_bulk = %f\n"%(box_volume, rho_bulk, cs_bulk ))


    # Transport coefficients from GK

    # # E_3.0
    # T_sq = 0.0601964
    # T_qq = 0.3412958
    # T_qs = 0.0594056
    # T_ss = 0.0167278

    # E_1.5
    T_qq = 0.8860249969715818
    T_sq = 0.05275808785481095
    T_qs = 0.05275808785481095
    T_ss = 0.0154820979540676

    final_p = p_simulations(p_pat, p_dir, box_volume,rho_bulk,cs_bulk,T_sq,T_qq)
    final_mu = mu_simulations(m_pat, m_dir, box_volume,rho_bulk,cs_bulk,T_qs,T_ss)
    
    plot_properties(final_p, 'grad_p','Q', x_label = r'$-\nabla P$', y_label = r'$Q$', plot_name ="qq" )
    plot_properties(final_p, 'grad_p','Js_exc', x_label = r'$-\nabla P$',  y_label = r'$J_s-c_s^BQ$', plot_name ="sq" )
    plot_properties(final_mu, 'grad_mu','Q',x_label = r'$-\nabla \mu^\prime_s$', y_label = r'$Q$' ,plot_name ="qs" )
    plot_properties(final_mu, 'grad_mu','Js_exc', x_label = r'$-\nabla \mu^\prime_s$', y_label = r'$J_s-c_s^BQ$', plot_name ="ss"  )
    


# TODO change group pattern to something more flexible as just file_names
if __name__ == "__main__":
    """
    THIS IS VERY SPECIFIC
    The arguments of this depend on the application
    """
    parser = argparse.ArgumentParser(description='Launch simulations from restart',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m_pat', metavar='m_pat',help='Generic name of the directories with mu gradient',default='mu_force_*')
    parser.add_argument('-p_pat', metavar='p_pat',help='Generic name of the directories with p gradient',default='p_force_*')
    parser.add_argument('-m_dir', metavar='m_dir',help='Patter of the files inside, in this case restart are like 202000',default='[0-9]*')
    parser.add_argument('-p_dir', metavar='p_dir',help='Patter of the files inside, in this case restart are like 202000',default='[0-9]*')
    args = parser.parse_args()

    
    main(args.m_pat,args.p_pat,args.m_dir,args.p_dir)
    
    
    

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