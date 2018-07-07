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



sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
except ImportError as err:
    print err

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script plots several files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file_name', metavar='InputFile',help='Input filename',nargs='+',type=lambda x: cf.is_valid_file(parser, x))
    parser.add_argument('-columns', metavar='columns',help='Properties to plot',nargs=2,default=[0,1],type=int)
    parser.add_argument('-path_name',metavar='path_name',help='depth in tree to define name',default=3,type=int )
    parser.add_argument('-plot_name',metavar='plot_name',help='Name of the pdf file generated', default='plotter.pdf',type=str)
    args = parser.parse_args()
    files=args.file_name
    columns=args.columns
    path_name=args.path_name


def get_name(file_name,path_name ):
    path=os.path.abspath(file_name).split("/")
    name_no_extension=path[-1].split(".")[-2]
    name='%s_%s'%(name_no_extension,path[-path_name])
    return name


def general_plotter(files,columns,path_name=3):
    """
    Plot a property from several files
    Args:
        files the files that contain the data
        columns an array [x.y] containing the positions you want to plot
        path_name automatically gives a name for the legends of the plot taken from the path, so this is the index of the 
        name counting from the back that follows the name of the file. 
    
    Returns:
        fig,ax to be handled and later customized.
    
    """
    
    x_index=columns[0]
    y_index=columns[1]
    
    fig,ax=plt.subplots()
    
    for fil in files:
        data=cf.read_data_file(fil).values
        name=get_name(fil)
        
        ax.plot(data[:,x_index],data[:,y_index],label=name)
    
    ax.legend()
    
    return ax,fig

ax,fig=general_plotter(files,columns)  
ymin,ymax=plt.ylim()
ax.set_ylim(ymin,ymax*1.2)  #To add 20% more in the y direction to fit the legend
plt.tight_layout()
ax.grid()
fig.savefig(args.plot_name)

plt.close()


