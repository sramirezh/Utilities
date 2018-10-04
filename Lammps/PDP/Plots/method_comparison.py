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
from general_plotter import pre_processing

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf


try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
except ImportError as err:
    print err

parser = argparse.ArgumentParser(description='This script plots the method files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file_name', metavar='InputFile',help='Input filename',nargs='+',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-columns', metavar='columns',help='Properties to plot',nargs=2,default=[0,1],type=int)
parser.add_argument('-path_name',metavar='path_name',help='depth in tree to define name',default=3,type=int )
parser.add_argument('-plot_name',metavar='plot_name',help='Name of the file generated, including the extension', default='plotter.png',type=str)
args = parser.parse_args()
files=args.file_name
columns=args.columns
path_name=args.path_name


"""
This is the general structure of anything in a file
"""

data,names=pre_processing(files,path_name)



colors=['r','b','k','g']
dic_legenend={'force':r'$C_s^B [\sigma^{-3}]$','vx_poly':r'$V_p^x[\sigma/\tau]$','rg_ave':r'$R_g [\sigma]$','rRg2':r'$R_{g}^2 [\sigma^2]$'}


#def general_plotter(data,columns=[0,1],marker="."):
#    """
#    Plot a property from several files
#    Args:
#        data that wants to be plotted
#        columns an array [x.y] containing the positions you want to plot
#        path_name automatically gives a name for the legends of the plot taken from the path, so this is the index of the
#        name counting from the back that follows the name of the file.
#
#    Returns:
#        fig,ax to be handled and later customized.
#
#    """
#    x_index=columns[0]
#    y_index=columns[1]
#
#    fig,ax=plt.subplots()
#
#    for dat in data:
#
#
#        ax.plot(dat[:,x_index],dat[:,y_index],marker=marker)
#
#    ax.legend()
#
#
#    return ax,fig



