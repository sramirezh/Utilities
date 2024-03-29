#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 13:50:16 2020

Plots all the density distributions generated by density distribution from atom
that are in 16. DCV-SGCMC

The distributions are inside each folder in 2.DCV

Got the distributions ready to scp as follows:
rsync -avz --include="distribution_E_*" --include="*/" --exclude="*" --progress
--partial 2.DCV ~/distributions
@author: sr802
"""

import warnings
import sys
import os
import glob
import numpy as np

warnings.filterwarnings("ignore")


sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf

import matplotlib.pyplot as plt

# =============================================================================
# Main
# =============================================================================
plot_dir = "plots/1.DCV_concentration_distribution"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)  



# Input data
Lx = 51.29928 
#External parameters 
lattice_constant = (4/0.8)**(1/3)
r_colloid  = 3.23


files = glob.glob("E*/distribution*")

cf.set_plot_appearance()

fig,ax = plt.subplots()

markers = ['v','^','<','>','s','o','D','p']
ax.axhline(y=0.6, xmin = 0, xmax=1, ls='--',c='black')
ax.axhline(y=0.15, xmin = 0, xmax=1, ls='--',c='black')

# sorting the files by the n
a = cf.extract_digits(files, sort = False)
a = np.c_[a, files]
sort_a = a[a[:,0].argsort()]
files = sort_a[:,-1]

# Only include the following epsilons
epsilon_included = [0.5, 1, 1.5]

for i, file in enumerate(files):

    data = cf.read_data_file(file).values
    x = data[:,0]
    solute_dist  = data[:,1]
    epsilon = file.strip('.dat').split('_')[-1]
    print(epsilon)
    if float(epsilon) in epsilon_included:
        print (epsilon)
        ax.plot(x,solute_dist, label=r'$\varepsilon_{cs} = %s$'%epsilon, ls = '--', marker = markers[i], ms = 5)


ax.axvspan(Lx/2-r_colloid,Lx/2+r_colloid, alpha=0.3, color='green')
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$c_s^B(x)$")
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlim(0, x[-1]) 
ax.set_ylim(0,0.87)
# ax.axvspan(6*lattice_constant,9*lattice_constant, alpha=0.5, color='blue')
# ax.axvspan(21*lattice_constant,24*lattice_constant, alpha=0.5, color='red')

ax.legend(loc='upper left', labelspacing=0.5, borderpad=0.4, scatteryoffsets=[0.6], 
           frameon=True, fancybox=False, edgecolor='k',ncol=1, fontsize = 10)

# For more information about the legendloc
# https://stackoverflow.com/questions/39803385/what-does-a-4-element-tuple-argument-for-bbox-to-anchor-mean-in-matplotlib
#ax.legend(loc= 1,labelspacing=0.2,borderpad=0.4,scatteryoffsets=[0.2],
#           frameon=True, fancybox=False, edgecolor='k',ncol=3, fontsize = 14, 
#           columnspacing = 0.5, bbox_to_anchor=(0.025, 0.42, 0.95, 0.6), mode="expand")
plt.tight_layout()
plt.savefig("%s/conc_chunks_%s.pdf"%(plot_dir,i+1), transparent = True)
logger.info("Created the plot %s/conc_chunks_%s.pdf"%(plot_dir,i+1))


#ax.plot(x,rho_solu_c, label='Solutes', ls = '--')
#ax.plot(x,rho_solv_c, label='Solvents', ls = '--')
##ax.plot(x, rho_total_c, label='Total')  
#ax.axhline(y=0.6, xmin = 0, xmax=1, ls='--',c='black')
#ax.axhline(y=0.15, xmin = 0, xmax=1, ls='--',c='black')
##ax.axvline(x=L[0]/2, ymin = 0, ymax=1, ls=':',c='black')
#ax.set_xlim(0, L[0])   
#
##ax.axvspan(6*lattice_constant,9*lattice_constant, alpha=0.5, color='blue')
##ax.axvspan(21*lattice_constant,24*lattice_constant, alpha=0.5, color='red')
#
#ymin,ymax=ax.get_ylim()
#ax.set_ylim(0,ymax*1.3)
#

#ax.legend()
#plt.tight_layout()
#plt.savefig("conc_chunks.pdf")
#
#
#epsilon =os.getcwd().split('/')[-1]