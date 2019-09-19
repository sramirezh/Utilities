#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  7 13:47:44 2018
General plotter for several files, it helps to analyse results fast
@author: simon
"""
import sys
import os
import argparse
import numpy as np


sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
except ImportError as err:
    print(err)


"""
###############################################################################
Function definition
###############################################################################
"""

def get_name(file_name,path_name):
    path=os.path.abspath(file_name).split("/")
    name_no_extension=os.path.splitext(file_name)[0]
    name='%s_%s'%(name_no_extension,path[-path_name])
    return name


def general_plotter(data,columns=[0,1],xerror=None,yerror=None,marker=None):
    """
    Plot a property from several files
    Args:
        data that wants to be plotted
        columns an array [x.y] containing the positions you want to plot
        xerror,yerror give the position in the data of the errors
        path_name automatically gives a name for the legends of the plot taken from the path, so this is the index of the
        name counting from the back that follows the name of the file.

    Returns:
        fig,ax to be handled and later customized.

    """

    colors=['r','b','k','g']
    error_cap=4

    x_index=columns[0]
    y_index=columns[1]

    fig,ax=plt.subplots()

    for i,dat in enumerate(data):
        if yerror == None: yerr_vals=yerror
        else: yerr_vals=dat[:,yerror]

        if xerror == None: xerr_vals=xerror
        else: xerr_vals=dat[:,xerror]

        #Defining the first colors from array and the rest by random numbers
        if i<len(colors):color=colors[i]
        else: color=np.random.rand(3)

        ax.errorbar(dat[:,x_index],dat[:,y_index],xerr=xerr_vals,yerr=yerr_vals,fmt='o',color=color,capsize=error_cap)




    return ax,fig

def pre_processing(files,path_name):
    """
    Added this to have an intermediate step before plotting
    Returns:
        data The data contained in the files
        names a reference name, that contains the file name without extension and a reference folder 2 levels above (see -path_name)

    """
    data=[]
    names=[]
    for fil in files:
        data.append(cf.read_data_file(fil).values)
        names.append(get_name(fil,path_name))

    return data,names

"""
###############################################################################
Main
###############################################################################
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script plots several files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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

    """
    Process the data here before plotting
    """

    ax,fig=general_plotter(data,columns)
    plt.legend(names)
    ymin,ymax=plt.ylim()
    ax.set_ylim(ymin,ymax*1.2)  #To add 20% more in the y direction to fit the legend
    plt.tight_layout()
    ax.grid()
    fig.savefig(args.plot_name)

    plt.close()
