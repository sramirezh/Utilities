#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Plots the concentrations found in the log file
@author: sr802
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import warnings
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


#Input parameters

T=2.0
mu=-2.0
target=np.exp(mu/T)


run_thermo_directories(".","log.mu-2.0",0)
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

ax1.axhline(y=target, xmin=0, xmax=1,ls=':',c='red')
ax1.set_xlabel(r'$steps $')
ax1.set_ylabel(r'$c[\sigma^{-3}] $')
#ax1.set_xlim(0,1)
#ax1.set_ylim(-5,10)
ax1.legend(loc='upper_left')
fig1.tight_layout()
fig1.savefig("Concentration.pdf")
plt.show()