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
logger.info("All the files load here were created by hand, while gathering all the data in  the excel Colloid results")
# Loading all the compiled data
dcv = np.loadtxt("DCVMC.dat",skiprows = 1 )
fdmd_free = np.loadtxt("FDMD_free.dat", skiprows = 1)
fdmd_fixed = np.loadtxt("FDMD_fixed.dat", skiprows = 1)
theory = np.loadtxt("Theory.dat", skiprows = 1)
theory_r = np.loadtxt("Theory_rh.dat", skiprows = 1)
theory_r_ns = np.loadtxt("Theory_rh_non_slip.dat", skiprows = 1)
rh = np.loadtxt("rh.dat", skiprows = 1)
diffusion = np.loadtxt("D.dat", skiprows = 1)
cf.set_plot_appearance()



# Ploting the comparison between BD-NEMD, FD-NEMD and Theory
fig, ax = plt.subplots()

ax.errorbar(dcv[:,0], -dcv[:,1], yerr= dcv[:,2], label="BD-NEMD", ls = '--', fmt='o')
ax.errorbar(fdmd_free[:,0],-fdmd_free[:,1],yerr= fdmd_free[:,2],label="FD-NEMD", fmt='o', ls = '--')
#ax.errorbar(fdmd_fixed[:,0],fdmd_fixed[:,1],yerr= fdmd_fixed[:,2],label="EFMD fixed", fmt='o')
#ax.scatter(theory[:,0], theory[:,1], label="Theory v0")
ax.plot(theory[:,0], theory[:,2], label="Theory", marker = "s", ls = '--')
#ax.plot(theory_r[:,0], theory_r[:,2], label="Theory $r=R_h$")
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xlabel(r'$\varepsilon_{cs} $')
ax.set_ylabel(r'$v_c^x$')
ymin,ymax=plt.ylim()
ax.set_ylim(ymin, ymax*1.3)
ax.legend(loc='upper right',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
       frameon=True, fancybox=False, edgecolor='k',ncol = 1, fontsize = 10)
fig.tight_layout()
plt.savefig("%s/Result_comparison_1.pdf"%plot_dir)
logger.info("plotted %s/Result_comparison_1.pdf"%plot_dir)

# Ploting the comparison between FD-NEMD fixed and mobile a


fig, ax = plt.subplots()


ax.errorbar(fdmd_free[:,0],-fdmd_free[:,1],yerr= fdmd_free[:,2],label="Mobile", fmt='o', ls = '--')
ax.errorbar(fdmd_fixed[:,0],fdmd_fixed[:,1],yerr= fdmd_fixed[:,2],label="Fixed", fmt='o',ls = '--')
#ax.scatter(theory[:,0], theory[:,1], label="Theory v0")
#ax.plot(theory[:,0], theory[:,2], label="Theory", marker = "s")
#ax.plot(theory_r[:,0], theory_r[:,2], label="Theory $r=R_h$")
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xlabel(r'$\varepsilon_{cs} $')
ax.set_ylabel(r'$v_c^x$')
ymin,ymax=plt.ylim()
ax.set_ylim(ymin, ymax*1.3)
ax.legend(loc='upper right',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
       frameon=True, fancybox=False, edgecolor='k',ncol = 1, fontsize = 10)
fig.tight_layout()
plt.savefig("%s/Result_comparison_2.pdf"%plot_dir)
logger.info("plotted %s/Result_comparison_2.pdf"%plot_dir)

# Ploting the comparison between Theoretical approximations


fig, ax = plt.subplots()


ax.errorbar(fdmd_free[:,0],-fdmd_free[:,1],yerr= fdmd_free[:,2],label="FD-NEMD", fmt='o', ls = '--')
#ax.plot(theory[1:,0], theory[1:,1], label=r"Theory $v_0$", marker = "s", ls = '--')
ax.plot(theory[:,0], theory[:,2], label=r"$\text{Theory} \quad a=\sigma_{cs}$", marker = "s", ls = '--')
ax.plot(theory_r[:,0], theory_r[:,2], label=r"$\text{Theory} \quad a=R_H^{\text{slip}}$", marker = "s", ls = '--')
ax.plot(theory_r_ns[:,0], theory_r_ns[:,2], label=r"$\text{Theory} \quad a=R_H^{\text{non-slip}}$", marker = "s", ls = '--')
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xlabel(r'$\varepsilon_{cs} $')
ax.set_ylabel(r'$v_c^x$')
ymin,ymax=plt.ylim()
ax.set_ylim(ymin, ymax*1.3)
ax.legend(loc='upper right',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
       frameon=True, fancybox=False, edgecolor='k',ncol = 1, fontsize = 10)
fig.tight_layout()
plt.savefig("%s/Result_comparison_3.pdf"%plot_dir)
logger.info("plotted %s/Result_comparison_3.pdf"%plot_dir)



# Ploting the radius


fig, ax = plt.subplots()

#ax.plot(theory[1:,0], theory[1:,1], label=r"Theory $v_0$", marker = "s", ls = '--')
ax.plot(rh[:,0], rh[:,2], label=r"$R_h^{\text{slip}}$", marker = "s", ls = '--')
ax.plot(rh[:,0], rh[:,1], label=r"$R_h^{\text{non-slip}}$", marker = "s", ls = '--')
ax.axhline(y=3.23, xmin=0, xmax=1,ls='--',c='black', label = r"$r =\sigma_{cs}$" )

ax.set_xlabel(r'$\varepsilon_{cs} $')
ax.set_ylabel(r'$r$')
ymin,ymax=plt.ylim()
ax.set_ylim(ymin, ymax*1.3)
ax.legend(loc='upper right',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
       frameon=True, fancybox=False, edgecolor='k',ncol = 1, fontsize = 10)
fig.tight_layout()
plt.savefig("%s/R_H.pdf"%plot_dir)
logger.info("plotted %s/R_H.pdf"%plot_dir)


# Ploting the Diffusion coefficient

fig, ax = plt.subplots()

#ax.plot(theory[1:,0], theory[1:,1], label=r"Theory $v_0$", marker = "s", ls = '--')
ax.plot(diffusion[:,0], diffusion[:,1], marker = "s", ls = '--')

ax.set_xlabel(r'$\varepsilon_{cs} $')
ax.set_ylabel(r'$D$')
ymin,ymax=plt.ylim()
ax.set_ylim(ymin, ymax*1.3)
#ax.legend(loc='upper right',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
#       frameon=True, fancybox=False, edgecolor='k',ncol = 1, fontsize = 10)
fig.tight_layout()
plt.savefig("%s/D.pdf"%plot_dir)
logger.info("plotted %s/D.pdf"%plot_dir)


# =============================================================================
# Plot of the v or Pe vs interaction 
# =============================================================================
# The following is to plot the evolution

dcv = np.loadtxt("DCVMC2.dat",skiprows = 1 )

D_s = 0.13036074653166377 # Diffusion coefficient for the solutes in the bulk
L = 20.51976 # Distance between the source and the Sink

pref = L/D_s

for end in range(len(dcv[:,0])):

    cf.set_plot_appearance()
    fig, ax = plt.subplots()
    
    
    
    ax.errorbar(dcv[:end+1,0], -dcv[:end+1,1], yerr= dcv[:end+1,2], fmt='o', c = 'k')
    
    
    ax.set_xlabel(r'$\varepsilon_{cs} $')
    ax.set_ylabel(r'$v^x$')
    
    ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
    ax.set_ylim(-0.01, 0.025)
    ax.set_xlim(0, 5.5)
    ax.set_xticks(np.arange(0, 6, 1.0))

    fig.tight_layout()
    plt.savefig("%s/DCVM_result_%s.pdf"%(plot_dir,end), transparent = True)
    
    
    # Peclet number Pe = L*v/D_c D_c for colloid
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
    
    
    # Peclet number Pe = L*v/D_s
    
    fig3, ax3 = plt.subplots()
    
    
    
    ax3.errorbar(dcv[:end+1,0], -pref*dcv[:end+1,1], yerr= pref*dcv[:end+1,2], fmt='o', c = 'k')
    
    
    ax3.set_xlabel(r'$\varepsilon_{cs} $')
    ax3.set_ylabel(r'$Pe$')
    
    ax3.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
#    ax.set_ylim(-0.01, 0.025)
    ax3.set_xlim(0, 5.5)
    ax3.set_xticks(np.arange(0, 6, 1.0))

    fig3.tight_layout()
    plt.savefig("%s/Peclet2_%s.pdf"%(plot_dir,end), transparent = True)
    
    
# The final plot The velocities including the negative value
    
ax.set_ylim(-0.02, 0.025)
fig.tight_layout()
fig.savefig("%s/DCVM_result_all.pdf"%(plot_dir), transparent = True)
    
    
    
    
    