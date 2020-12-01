#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 28/02/2019
Plots the Mobility vs the interactions for several input files generated perhaps from different potentials.
@author: sr802
"""

import sys
import os
import argparse
import numpy as np
from general_plotter import general_plotter,pre_processing
import pandas as pd
import matplotlib.pyplot as plt
import Lammps.PDP.trajectory_analysis.first_n_analysis as fna

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf


parser = argparse.ArgumentParser(description='This script plots several files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file_name', metavar='InputFile',help='Input filename',nargs='+',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-columns', metavar='columns',help='Properties to plot',nargs=2,default=[0,1],type=int)
parser.add_argument('-path_name',metavar='path_name',help='depth in tree to define name',default=3,type=int )
parser.add_argument('-plot_name',metavar='plot_name',help='Name of the pdf file generated', default='plotter.pdf',type=str)
args = parser.parse_args()
files=args.file_name
columns=args.columns
path_name=args.path_name


# Creating the logger
logger = cf.log(__file__, os.getcwd())   

logger.info("Using the following arguments for the paser")
logger.info(args)


"""
This is the general structure of anything in a file
"""

data,names=pre_processing(files,path_name)


"""
converting into pandas df
"""

for i in range(len(data)):
    data[i]=pd.DataFrame(data=data[i][1:,1:],index=data[i][1:,0], columns=data[i][0,1:])



"""
Data pre-processing
"""

index_mobility=cf.parameter_finder(data[0].columns,"mobility")

data_to_plot=[]
for dat in data:
    epsilon=cf.extract_digits(list(dat.index))[:,0]
    data_to_plot.append(np.column_stack((epsilon,np.array(dat.values[:,index_mobility],dtype=float))))




"""
Plotting
"""


names = ["GLJ","LJ"]
plt.close('all')
cf.set_plot_appearance()
ax,fig=general_plotter(data_to_plot,yerror=2)



ax.set_xlabel(r'$\varepsilon_{ms} $')
ax.set_ylabel(r'$M_{ps} [\tau/m]$')
plt.legend(names, loc='upper left')
ax.axhline(y=0.07155442595335286, xmin=0, xmax=1,ls='--',c='black')
ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
ax.set_xlim(0,ax.get_xlim()[1])
plt.xticks(np.arange(0,14,2))


# =============================================================================
# Creating the insert
# =============================================================================
insert= True
if insert == True:
    rmin = 0.985
    rmax = 2.5
    potentials = fna.get_potentials(rmin,rmax,epsilon_lj=1)
    left, bottom, width, height = [0.55, 0.25, 0.4, 0.30]
    ax2 = fig.add_axes([left, bottom, width, height])
    
    ax2.set_ylabel(r'$V(r)$',fontsize =17, labelpad=-5)
    ax2.set_xlabel(r'$r$' ,fontsize =17, labelpad=-5)
    ax2.plot(potentials[:,0],potentials[:,1],label="GLJ")
    ax2.plot(potentials[:,0],potentials[:,2],label="LJ")
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
    
    

    
    
    
fig.tight_layout()
fig.savefig(args.plot_name, transparent = True)
#plt.show()

