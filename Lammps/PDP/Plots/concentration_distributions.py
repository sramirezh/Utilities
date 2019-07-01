#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 15:47:20 2019
Plots all the concentration files inside the folders
@author: sr802
"""

import sys
import os
import numpy as np
import warnings
import glob
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf



Rh=[2.758531,2.755166] #LJ,GLJ
cf.set_plot_appearance()


directories=glob.glob('*/')

name={"u":"Solute","v":"Solvent","t":"Solution"}
colors={"u":"red","v":"blue","t":"black"}
ltype={"u":"-","v":"-","t":"--"}

plt.close('all')

for counter,directory in enumerate(directories):
    files=glob.glob('%s/prof_*.dat'%directory)
    fig,ax=plt.subplots()
    files=sorted(files)[::-1]
    
    for f in files:
        key=f.split('_')[1][0]
        data=cf.read_data_file(f).values
        plt.plot(data[:,1],data[:,3],color=colors[key],label=name[key],linestyle=ltype[key])
        
    plt.legend(loc='upper right')
    ax.set_ylabel(r'$c[\sigma^{-3}]$')
    ax.set_xlabel(r'$r[\sigma]$')
    ax.set_xlim(0,9)
    ax.set_ylim(0,ax.get_ylim()[1])
    ax.axvline(x=Rh[counter], ymin=0, ymax=1,ls='-.',c='black')
    fig.tight_layout()
    fig.savefig('%s.pdf'%directory.split('/')[0])

plt.show()