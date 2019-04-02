#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 28/02/2019
Plots the Mobility vs the interactions for several input files generated perhaps from different potentials.
@author: sr802
"""
from __future__ import division
import sys
import os
import argparse
import numpy as np
from general_plotter import general_plotter,pre_processing
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
except ImportError as err:
    print err

import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='This script plots several files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file_name', metavar='InputFile',help='Input filename',nargs='+',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-columns', metavar='columns',help='Properties to plot',nargs=2,default=[0,1],type=int)
parser.add_argument('-path_name',metavar='path_name',help='depth in tree to define name',default=3,type=int )
parser.add_argument('-plot_name',metavar='plot_name',help='Name of the pdf file generated', default='plotter.pdf',type=str)
args = parser.parse_args()
files=args.file_name
columns=args.columns
path_name=args.path_name

"""
This is the general structure of anything in a file
"""

data,names=pre_processing(files,path_name)


"""
converting into pandas df
"""

for i in xrange(len(data)):
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
names=["SRLJ","LJ"]
plt.close('all')
cf.set_plot_appearance()
ax,fig=general_plotter(data_to_plot,yerror=2)



ax.set_xlabel(r'$\epsilon_{ms} $')
ax.set_ylabel(r'$\Gamma_{ps} [\tau/m]$')
plt.legend(names,loc='upper_left')
ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
ax.set_xlim(0,ax.get_xlim()[1])
plt.xticks(np.arange(0,14,2))
fig.tight_layout()
fig.savefig(args.plot_name)
plt.show()
