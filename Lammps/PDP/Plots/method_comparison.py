#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 12:36:57 2018
This plots the results from different methods, moving
@author: sr802
"""

import sys
import os
import argparse
import numpy as np
from .general_plotter import pre_processing, general_plotter

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf


try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt


except ImportError as err:
    print(err)

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

data,names=pre_processing(files,path_name)


"""
Preprocessing data
"""
index=cf.parameter_finder(names,"free")
data[1][:,1]*=-1 #As the free polymer moves in the opposite direction of the flow

for i,name in enumerate(names):
    names[i]=name.split('_')[-2].capitalize()

"""
Plot and modifications
"""

ax,fig=general_plotter(data,yerror=2)

#General plot parameters
axis_font=24
tick_font=20
legend_font=18
xoffset=0.1
yoffset=0.1
error_cap=4
colors=['r','b','k']



#"""
#Changing properties of the lines
#"""
#for i,line in enumerate(ax.lines):
#    line.set_color(colors[i]) #The line colors
#    ax.collections[i].set_color(colors[i]) #The error bars




"""Legend"""
plt.legend(names,fontsize=legend_font,loc='upper left',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
   frameon=True, fancybox=False, edgecolor='k')

"""Axis"""

ax.set_xlabel(r'$F_{s}^{\mu}=-\nabla_x \mu_s [\epsilon/\sigma]$',fontsize=axis_font)
ax.tick_params(labelsize=tick_font,direction='in',top=True, right=True)

ax.set_ylabel(r'$|V_p^x|[\sigma/\tau]$',fontsize=axis_font)

ymin,ymax=plt.ylim()
deltay=ymax-ymin
ax.set_ylim(0,ymax+deltay*0.45)

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

plt.grid(False)
try:
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["text.usetex"] = True
except:
    pass
plt.tight_layout()
plt.savefig(file_name)
plt.close()
