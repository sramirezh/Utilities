#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:56:37 2020
To plot the comparison between FDMD and GCMC
@author: sr802
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
from joblib import Parallel, delayed
import multiprocessing

from Lammps.PDP.Plots.LJ_plotter import LJ
import Lammps.core_functions as cf


dcv = np.loadtxt("DCVMC.dat",skiprows = 1 )
fdmd_free = np.loadtxt("FDMD_free.dat", skiprows = 1)
fdmd_fixed = np.loadtxt("FDMD_fixed.dat", skiprows = 1)
theory = np.loadtxt("Theory.dat", skiprows = 1)
theory_r = np.loadtxt("Theory_rh.dat", skiprows = 1)

cf.set_plot_appearance()
fig, ax = plt.subplots()

ax.errorbar(dcv[:,0], -dcv[:,1], yerr= dcv[:,2], label="Explicit", ls = '-', fmt='o')
ax.errorbar(fdmd_free[:,0],-fdmd_free[:,1],yerr= fdmd_free[:,2],label="Implicit", fmt='o', ls = '-')
#ax.errorbar(fdmd_fixed[:,0],fdmd_fixed[:,1],yerr= fdmd_fixed[:,2],label="EFMD fixed", fmt='o')
#ax.scatter(theory[:,0], theory[:,1], label="Theory v0")
ax.plot(theory[:,0], theory[:,2], label="Theory", marker = "s")
#ax.plot(theory_r[:,0], theory_r[:,2], label="Theory $r=R_h$")
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xlabel(r'$\varepsilon_{cs} $')
ax.set_ylabel(r'$v_x$')
ymin,ymax=plt.ylim()
ax.set_ylim(ymin, ymax*1.3)
ax.legend(loc='upper right',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
       frameon=True, fancybox=False, edgecolor='k',ncol = 1, fontsize = 10)
fig.tight_layout()
plt.savefig("Result_comparison.pdf")



# The following is to plot the evolution by doing the 

dcv = np.loadtxt("DCVMC2.dat",skiprows = 1 )

for end in range(len(dcv[:,0])):

    cf.set_plot_appearance()
    fig, ax = plt.subplots()
    
    
    
    ax.errorbar(dcv[:end+1,0], -dcv[:end+1,1], yerr= dcv[:end+1,2], fmt='o', c = 'k')
    
    
    ax.set_xlabel(r'$\varepsilon_{cs} $', fontsize = 40)
    ax.set_ylabel(r'$v_x$', fontsize = 40)
    
    ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
    ax.set_ylim(-0.01, 0.025)
    ax.set_xlim(0, 5.5)
    ax.set_xticks(np.arange(0, 6, 1.0))
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    fig.tight_layout()
    plt.savefig("DCVM_result_%s.pdf"%end, transparent = True)
    
    
    # Peclet number
    fig, ax = plt.subplots()
    
    
    
    ax.scatter(dcv[:end+1,0], dcv[:end+1,3])
    
    
    ax.set_xlabel(r'$\varepsilon_{cs} $')
    ax.set_ylabel(r'$Pe$')
    
    ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
    ax.set_ylim(-10, 100)
    ax.set_xlim(0, 5.5)
    ax.set_xticks(np.arange(0, 6, 1.0))
#    plt.xticks(fontsize=30)
#    plt.yticks(fontsize=30)
    fig.tight_layout()
    plt.savefig("Peclet_%s.pdf"%end, transparent = True)