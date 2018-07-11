#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 11:58:15 2018
Rg Vs N
@author: sr802
"""
from __future__ import division
import sys
import os
import argparse
import numpy as np
from general_plotter import general_plotter,pre_processing


sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
except ImportError as err:
    print err

parser = argparse.ArgumentParser(description='This script plots several files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file_name', metavar='InputFile',help='Input filename',nargs='+',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-columns', metavar='columns',help='Properties to plot',nargs=2,default=[0,1],type=int)
parser.add_argument('-path_name',metavar='path_name',help='depth in tree to define name',default=3,type=int )
parser.add_argument('-plot_name',metavar='plot_name',help='Name of the pdf file generated', default='plotter.pdf',type=str)
args = parser.parse_args()
files=args.file_name
columns=args.columns
path_name=args.path_name





data,names=pre_processing(files,path_name)
"""
###############################################################################
Specific processing of the data
###############################################################################
"""
data=np.log(data)

x=data[0][:,0]
y=3/5*x

deltay=data[0][0,1]-y[0]
y=y+deltay
"""
###############################################################################
Specific processing of the data
###############################################################################
"""

ax,fig=general_plotter(data,columns)

ax.plot(x,y,'--')

"""Legend"""

handles, labels = ax.get_legend_handles_labels() #This is to check the references
new_names=['$\epsilon=0.5, \sigma=0.8$','$\epsilon=1.0, \sigma=1,5$','$\epsilon=1.5, \sigma=1.5$','$R_g \sim N^{3/5}$']
plt.legend(new_names,fontsize=17,s)

"""Axis"""
ax.tick_params(labelsize=18)
ax.set_ylabel(r'$\log(R_g)$',fontsize=20)
ax.set_xlabel(r'$\log(N)$',fontsize=20)
ymin,ymax=plt.ylim()
ax.set_ylim(ymin,ymax*1.25)  #To add 20% more in the y direction to fit the legend
#ax.set_xlim(0,1)

"""Lines"""


"""General"""
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["text.usetex"] =True
ax.grid(False)
plt.tight_layout()
fig.savefig(args.plot_name)

plt.close()
