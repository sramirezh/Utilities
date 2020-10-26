#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Plots the concentration evolution found in the log file
@author: sr802
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import warnings
import argparse
warnings.filterwarnings("ignore")

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.General.Log_Analysis.Thermo_Analyser as thermo

cwd = os.getcwd() #current working directory

def run_thermo_directories(directories,log_name,dmin):
    """ 
    Taken from flux gather
    THIS IS A VERY GENERAL FUNCTION (Kind of replaces compute statistics from PDP)
    
    Runs the thermo analysis inside the specific directories
    
    Args:
        directories: list of directories to run the analysis
        dmin: the minimum of steps to be discarded as defined in fast_averager
    """

    for folder in directories:
        os.chdir(folder)
        thermo.thermo_analyser(log_name,dmin)
        os.chdir(cwd)

def plot_concentrations(log_file):
    #Input parameters for ideal case
    
    T=1.0
    mu=-3.2105263157894735
    target=np.exp(mu/T) #Only for the ideal case
    
    
    run_thermo_directories(".",log_file,0) #a special case of the function only in the current directory
    thermo_data=cf.read_data_file("Parameters.dat")
    dic_legends={'Density':'Total','v_rhoSolv':"Solvent",'v_rhoSolu':"Solute"}
    
    names=thermo_data.columns.values
    data=thermo_data.values
    
    cf.set_plot_appearance()
    
    ind_steps=cf.parameter_finder(names,"Step")
    
    #For single component
    indexes=cf.parameter_finder(names,"Dens")
    
    #If there are more components
    indexes.extend(cf.parameter_finder(names,"rhoS"))
    
    plt.close("all")
    fig1,ax1=plt.subplots()
    
    for i in indexes:
        ax1.plot(data[:,ind_steps],data[:,i],label=dic_legends[names[i]])
    
    #ax1.axhline(y=target, xmin=0, xmax=1,ls=':',c='red')
    ax1.set_xlabel(r'$steps $')
    ax1.set_ylabel(r'$c[\sigma^{-3}] $')
    #ax1.set_xlim(0,1)
    #ax1.set_ylim(-5,10)
    ax1.legend(loc='upper_left')
    fig1.tight_layout()
    fig1.savefig("Concentration.pdf")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script extracs the denisties of the different components and generate a plots",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-log_name',help='name of the log file',default="log.lammps")
    args = parser.parse_args()
    plot_concentrations(args.log_name)
