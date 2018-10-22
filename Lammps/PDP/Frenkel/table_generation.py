#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 12:56:45 2018
This files creates the table for the frenkel potential to be used in Lammps
@author: simon
"""

from __future__ import division
import argparse
import pandas as pd
import numpy as np
import warnings
import sys 
import os
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")


sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path


import Lammps.PDP.Plots.LJ_plotter as ljplot

rmin=0.95
rc=1.6
n_points=500
epsilon=2.5
sigma=1.0



def generate_points(epsilon, sigma,n,rc,n_points):
    """
    generates the points for the potential
    """
    
    data=np.zeros((n_points,3))
    data[:,0]=-(rc-rmin)*np.cos(np.linspace(0,np.pi/2,num=n_points))+rc
    data[:,1]=ljplot.frenkel(data[:,0],epsilon,sigma,rc,4)
    data[:,2]=ljplot.force_frenkel(data[:,0],epsilon,sigma,rc,4)
    
    return data



def write_potential(f,epsilon, sigma,n,rc,n_points):
    """
    writing the file as described in https://lammps.sandia.gov/doc/pair_table.html
    """
    
    
    data=generate_points(epsilon,sigma,n,rc,n_points)
    f.write("\n")
    f.write("FRENKEL_N_%s_E_%s\n"%(n,epsilon))
    f.write("N %d\n"%(n_points)) #Did not use R or any other option here, cus the points are crerated with a sinusoidal distribution
    f.write("\n")
    for i,point in enumerate(data):
        f.write("%d %f %f %f\n"%(i+1,data[i,0],data[i,1],data[i,2]) )
        
        
    return f,data
      


def run_lammps():
    """
    Short lammps test
    """
    from lammps import IPyLammps
    
    L = IPyLammps() #Creates the object
    L.units("lj")
    L.atom_style("atomic")
    #L.atom_modify("map array")
    
    L.region("box block", 0, 6, 0, 6, 0, 6)
    L.create_box(2, "box")
    L.create_atoms(1, "single 0 0 0")
    L.create_atoms(1, "single 0 1.12 0")
    L.create_atoms(2, "single 3 0 0")
    L.create_atoms(2, "single 3 1.12 0")
    L.mass('*', 1.0)
    L.neighbor(0.3, "bin")
    L.command("pair_style table linear 1000")
    L.pair_coeff(" 1 1 frenkel.dat FRENKEL_N_4_E_2.5")
    L.command("shell rm frenkel_lammps.dat")
    L.command("pair_write 1 1 600 r %f %f  frenkel_lammps.dat Frenkel_lammps" %(rmin,rc) )
    L.write_script("input.lmp")
    L.run(0)

def plot_test(data_generated):
    """
    reading the lammps data
    """
    data=np.genfromtxt("frenkel_lammps.dat",skip_header=5)
    
    """
    plotting everything
    """
    plt.close('all')
    plt.plot(data_generated(), phi,'*',label="Generated potential")
    
    plt.plot(data[:,1],data[:,2],label="Lammps Potential")
    
    plt.legend()
    
    
def main():
    
    f=open("frenkel.dat",'w')
    f,data=write_potential(f,2.5,1.0,4,1.6,500)
    f.close()
    
main()