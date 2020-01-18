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


from ovito.modifiers import *
import ovito.io as ov





dcv = np.loadtxt("DCVMC.dat",skiprows = 1 )
fdmd = np.loadtxt("FDMD.dat", skiprows = 1)

cf.set_plot_appearance()
fig, ax = plt.subplots()

ax.errorbar(dcv[:,0],np.abs(dcv[:,1]),yerr= dcv[:,2],label="DCV", fmt='o')
ax.errorbar(fdmd[:,0],np.abs(fdmd[:,1]),yerr= fdmd[:,2],label="EFMD", fmt='o')

ax.set_xlabel(r'$\varepsilon $')
ax.set_ylabel(r'$v_x$')
ymin,ymax=plt.ylim()
ax.set_ylim(0,ymax*1.2)
ax.legend(loc='upper right',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
       frameon=True, fancybox=False, edgecolor='k',ncol=1)
fig.tight_layout()
plt.savefig("Result_comparison.pdf")

