#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 12:26:36 2020
This script returns the viscosity as predicted by Meyer 2018
and retrieves all the results from the folder, it has to be run inside 
Meyer_2018
@author: simon
"""

import os
import sys
import numpy as np
import glob
import matplotlib.pyplot as plt
from  scipy.stats import sem
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

def eta_meyer(rho, temp):
    """
    returns the eta as given by
    Meyer, N., Wax, J. F. and Xu, H. (2018) 
    ‘Viscosity of Lennard-Jones mixtures: A systematic study and 
    empirical law’, Journal of Chemical Physics, 148(23).
    
    Args:
        rho  reduced density
        temp reduced temperature
    
    Returns:
        eta: viscosity of the LJ system
    """
    
    a = [0.007]
    b = [-2.2621, 3.1567]
    c = [-2.6643, 5.3625]

    eta = (a[0]*temp**2 + (c[0]+c[1]*rho))*np.exp((b[0]+b[1]*rho)/temp)
    
    return eta



def retrieve_results():
    """
    Goes to each temperature T_X and then inside each subfolder
    getting all the viscosisies.
    
    
    Returns:
        results: object containing an entry for each Temperature like
        [[T_X], [viscosities list]]
    """
    cwd = os.getcwd() #current working directory
    folders = glob.glob("T*")
    
    results = []
    for f in folders:
        logger.info("Inside folder %s"%f)
        viscosity = []
        temperature = float(f.split('_')[-1])
        results_temperature = [temperature]
        os.chdir(f)
        subfolders = glob.glob('[1-9]*')
        for s in subfolders:
            os.chdir(s)
            line, error = cf.bash_command("""tail -2 log.lammps| head -1""")
            viscosity.append(cf.extract_digits(line.decode('utf-8'))[0])
            os.chdir('../')
        results_temperature.append(viscosity)
        
        np.savetxt("viscosities.dat", np.array(viscosity))
        results.append(results_temperature)
        os.chdir(cwd)
        
    cf.save_instance(results, "results")
    
    return results


# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    logger = cf.log(__file__, os.getcwd())
    
    
    if not os.path.exists("results.pkl"):
        logger.info("\nRunning the analysis")
        results = retrieve_results()
        
    else: 
        logger.info("\nReading the results from results.pkl")
        results = cf.load_instance("results.pkl")
        
    
    
    # =============================================================================
    # Analysing the results
    # =============================================================================
    
    ave_results = np.zeros((len(results),3))  
    
    for i,res in enumerate(results):
        ave_results[i,0] = res[0] #  Getting the temperature
        ave_results[i,1] = np.average(res[1]) # Average
        ave_results[i,2] = sem(res[1])  # Statistic mean error
        
    
    
    # =============================================================================
    # Getting results from the empirical results in Meyer
    # =============================================================================
    rho = 0.8 
    
    logger.info("Obtaining results using the law in Meyer et al, with rho = %s" %rho)
    temperature = np.linspace(np.min(ave_results[:,0]), np.max(ave_results[:,0]))
    viscosity = eta_meyer(rho,temperature)
    # =============================================================================
    # Plot results
    # =============================================================================
         
    plt.close('all')
    cf.set_plot_appearance()
    
    #For c11
    fig,ax = plt.subplots()
    
    
    
    ax.plot(temperature, viscosity, label = 'Meyer et al')
    ax.errorbar(ave_results[:,0], ave_results[:,1], yerr = ave_results[:,2], label="GK", fmt='o')
    ax.set_xlabel(r'$T^*$')
    ax.set_ylim([0,3])
    ax.set_ylabel(r'$\eta^*$')
    plt.legend( loc = 'upper right')
    plt.tight_layout()
    plt.savefig("validating_eta.pdf")
    