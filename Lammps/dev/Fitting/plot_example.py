#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 17:45:06 2019

@author: simon
"""

import numpy as np
import argparse
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 


def read_data(fname,rho_ref,beta_ref):
    """
    Reads the data from the input file
    Args:
        fname: name of the file containing the data, run build example to see the structure 
        of the input_example.dat
        rho_ref: reference density
        beta_ref
        
    Returns:
        variables
    """
    data=np.loadtxt(fname)
    
    
    rho=data[:,0]-rho_ref
    temperature=data[:,1]
    prope=data[:,2]
    sigma_prope=data[:,3]
    beta=1/temperature-beta_ref
    
    
    
    return rho,beta, prope, sigma_prope

input_file="input_example.dat"
x,y, z, zerr=read_data(input_file,0,0)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, zdir='z',marker='.',label="Simulation",color='r')
ax.plot_wireframe(np.reshape(x,(30,30)),np.reshape(y,(30,30)),np.reshape(z,(30,30)),color='b',label="Fitting")
fig.legend()
fig.savefig("3Dplot.pdf")