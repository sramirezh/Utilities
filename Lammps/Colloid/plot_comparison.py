#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:56:37 2020
To plot the comparison between FDMD and GCMC, it needs to be run where all the
compiled data files are:
    
    DCVMC double control volume
    DCVMC2 same as before but with the Peclet number
    FDMD_free Field driven with the colloid free to move
    FDMD_fixed with the colloid fixed
    Theory for the theoretical results using the read radius of the colloid
    Theory_rh for the hydrodynamic radius
@author: sr802
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf



# =============================================================================
# Main
# =============================================================================
plot_dir = "plots/0.plot_comparison"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)    

# Loading all the compiled data
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
ax.set_ylabel(r'$v_c^x$')
ymin,ymax=plt.ylim()
ax.set_ylim(ymin, ymax*1.3)
ax.legend(loc='upper right',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
       frameon=True, fancybox=False, edgecolor='k',ncol = 1, fontsize = 10)
fig.tight_layout()
plt.savefig("%s/Result_comparison.pdf"%plot_dir)



# The following is to plot the evolution

dcv = np.loadtxt("DCVMC2.dat",skiprows = 1 )

for end in range(len(dcv[:,0])):

    cf.set_plot_appearance()
    fig, ax = plt.subplots()
    
    
    
    ax.errorbar(dcv[:end+1,0], -dcv[:end+1,1], yerr= dcv[:end+1,2], fmt='o', c = 'k')
    
    
    ax.set_xlabel(r'$\varepsilon_{cs} $')
    ax.set_ylabel(r'$v_c^x$')
    
    ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
    ax.set_ylim(-0.01, 0.025)
    ax.set_xlim(0, 5.5)
    ax.set_xticks(np.arange(0, 6, 1.0))

    fig.tight_layout()
    plt.savefig("%s/DCVM_result_%s.pdf"%(plot_dir,end), transparent = True)
    
    
    # Peclet number
    fig1, ax1 = plt.subplots()
    
    
    
    ax1.scatter(dcv[:end+1,0], dcv[:end+1,3])
    
    
    ax1.set_xlabel(r'$\varepsilon_{cs} $')
    ax1.set_ylabel(r'$Pe$')
    
    ax1.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
    ax1.set_ylim(-10, 100)
    ax1.set_xlim(0, 5.5)
    ax1.set_xticks(np.arange(0, 6, 1.0))
    fig1.tight_layout()
    plt.savefig("%s/Peclet_%s.pdf"%(plot_dir,end), transparent = True)
    
# The final plot The velocities including the negative value
    
ax.set_ylim(-0.02, 0.025)
fig.tight_layout()
fig.savefig("%s/DCVM_result_all.pdf"%(plot_dir), transparent = True)
    
    
    
    
    