#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 12:36:57 2018
This plots the results from different methods, moving

Just run it like run ~/dev/Utilities/Lammps/PDP/Plots/size_analysis.py *.dat
@author: sr802
"""

import sys
import os
import argparse
import numpy as np
from general_plotter import pre_processing, general_plotter
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

logger = cf.log(__file__, os.getcwd())  


parser = argparse.ArgumentParser(description='This script plots the method files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file_name', metavar='InputFile',help='Input filename',nargs='+',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-columns', metavar='columns',help='Properties to plot',nargs=2,default=[0,1],type=int)
parser.add_argument('-path_name',metavar='path_name',help='depth in tree to define name',default=3,type=int )
parser.add_argument('-plot_name',metavar='plot_name',help='Name of the file generated, including the extension', default='method_comparison.pdf',type=str)
args = parser.parse_args()
files=args.file_name
columns=args.columns
path_name=args.path_name
file_name=args.plot_name


"""
This is the general structure of anything in a file
"""

data,names = pre_processing(files, path_name)


"""
Preprocessing data
"""
index = cf.parameter_finder(names,"free")
data[1][:,1]*=-1 #As the free polymer moves in the opposite direction of the flow
# Very bad way of getting the names but I will not spend time
for i,file in enumerate(files):
    names[i] = file.split('_')[-1].split('.')[-2].capitalize()

"""
Plot and modifications
"""

ax,fig = general_plotter(data, yerror=2)

#"""
#Changing properties of the lines
#"""
#for i,line in enumerate(ax.lines):
#    line.set_color(colors[i]) #The line colors
#    ax.collections[i].set_color(colors[i]) #The error bars




"""Legend"""
plt.legend(names,loc='upper left',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
   frameon=True, fancybox=False, edgecolor='k')

"""Axis"""

ax.set_xlabel(r'$F_{s}^{\mu}=-\nabla_x \mu_s [\varepsilon/\sigma]$')
ax.tick_params(direction='in',top=True, right=True)

ax.set_ylabel(r'$|v_{\text{dp}}^x|[\sigma/\tau]$')

ymin,ymax=plt.ylim()
deltay=ymax-ymin
ax.set_ylim(0,0.0125)

xmin,xmax=plt.xlim()
deltax=xmax-xmin
ax.set_xlim(0,0.12)



plt.xticks(np.arange(0.02,0.12,0.02))
ax.spines["top"].set_visible(True)
ax.spines["right"].set_visible(True)


"""Lines"""
if ymin*ymax<0:
    ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')

"""General"""

plt.tight_layout()
plt.savefig(file_name)
plt.close()
