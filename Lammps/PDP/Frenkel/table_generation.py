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





def generate_points(epsilon, sigma,n,rc,n_points,rmin):
    """
    generates the points for the potential
    """
    data=np.zeros((n_points,3))
    data[:,0]=-(rc-rmin)*np.cos(np.linspace(0,np.pi/2,num=n_points))+rc  #Harmonic dist
    #data[:,0]=np.linspace(rmin,rc,n_points) #Linear dist
    data[:,1]=ljplot.frenkel(data[:,0],epsilon,sigma,rc,4)
    data[:,2]=ljplot.force_frenkel(data[:,0],epsilon,sigma,rc,4)
    
    return data



def write_potential(f,epsilon, sigma,n,rc,n_points,rmin):
    """
    writing the file as described in https://lammps.sandia.gov/doc/pair_table.html
    """
    
    
    data=generate_points(epsilon,sigma,n,rc,n_points,rmin)
    f.write("\n")
    f.write("FRENKEL_N_%s_E_%s\n"%(n,epsilon))
    f.write("N %d\n"%(n_points)) #Did not use R or any other option here, cus the points are crerated with a sinusoidal distribution
    f.write("\n")
    for i,point in enumerate(data):
        f.write("%d %f %f %f\n"%(i+1,data[i,0],data[i,1],data[i,2]) )
        
        
    return f,data
      


def run_lammps(rmin,rc,n_points):
    """
    Short lammps test cheking the table
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
    L.command("pair_style table linear %s" %n_points)
    L.pair_coeff(" 1 1 frenkel.dat FRENKEL_N_4_E_2.5")
    L.pair_coeff(" 1 2 frenkel.dat FRENKEL_N_4_E_1.0")
    L.pair_coeff("2 2 frenkel.dat FRENKEL_N_4_E_1.0")
    L.command("shell rm frenkel_lammps11.dat")
    L.command("shell rm frenkel_lammps22.dat")
    
    L.command("pair_write 1 1 600 r %f %f  frenkel_lammps11.dat Frenkel_lammps" %(rmin,rc) )
    L.command("pair_write 2 2 600 r %f %f  frenkel_lammps22.dat Frenkel_lammps" %(rmin,rc) )
    L.write_script("input.lmp")
    L.run(0)
    
    
def run_lammps_hybrid(rmin,rc, n_points):
    """
    Test with 3 types of particles and two potentials: LJ and frenkel
    """
    from lammps import IPyLammps
    
    L = IPyLammps() #Creates the object
    L.units("lj")
    L.atom_style("atomic")
    #L.atom_modify("map array")
    
    L.region("box block", 0, 6, 0, 6, 0, 6)
    L.create_box(3, "box")
    L.create_atoms(1, "single 0 0 0")
    L.create_atoms(1, "single 0 1.12 0")
    L.create_atoms(2, "single 3 0 0")
    L.create_atoms(3, "single 3 1.12 0")
    L.mass('*', 1.0)
    L.neighbor(0.3, "bin")
    L.command("pair_style lj/cut 2.5")
    L.command("pair_coeff * * 1.0 1.0 2.5")
    L.command("pair_style hybrid table linear %s lj/cut 2.5" %n_points)
    L.pair_coeff(" 1 1 table frenkel.dat FRENKEL_N_4_E_2.5")
    L.pair_coeff(" 1 2 table frenkel.dat FRENKEL_N_4_E_1.0")
    L.pair_coeff("2 2 table frenkel.dat FRENKEL_N_4_E_1.0")
    L.pair_coeff("* 3 lj/cut 1.0 1.0 ")
    L.command("shell rm frenkel_lammps11.dat")
    L.command("shell rm frenkel_lammps22.dat")
    
    L.command("pair_write 1 1 600 r %f %f  frenkel_lammps11.dat Frenkel_lammps" %(rmin,rc) )
    L.command("pair_write 2 2 600 r %f %f  frenkel_lammps22.dat Frenkel_lammps" %(rmin,rc) )
    L.write_script("table.lmp")
    L.run(0)

def plot_test(data11,data22):
    """
    reading the lammps data
    """
    lammps11=np.genfromtxt("frenkel_lammps11.dat",skip_header=5)
    lammps22=np.genfromtxt("frenkel_lammps22.dat",skip_header=5)
    
    """
    plotting everything
    """
    plt.close('all')
    fig1,(ax1,ax12)=plt.subplots(2,1)
    fig2,(ax2,ax22)=plt.subplots(2,1)
    
    ax1.plot(data11[:,0] , data11[:,1],'*',label="Generated potential 11")
    ax12.plot(data11[:,0] , data11[:,1],'*',label="Generated potential 11")

    ax2.plot(data22[:,0] , data22[:,1],'*',label="Generated potential 22")
    ax22.plot(data22[:,0] , data22[:,1],'*',label="Generated potential 11")

    
    ax1.plot(lammps11[:,1],lammps11[:,2],label="Lammps Potential")
    ax2.plot(lammps22[:,1],lammps22[:,2],label="Lammps Potential")
    
    ax12.plot(lammps11[:,1],lammps11[:,2],label="Lammps Potential")
    ax22.plot(lammps22[:,1],lammps22[:,2],label="Lammps Potential")
    
    ax22.set_xlim(1,1.6)
    ax22.set_ylim(-1.0,0)
    ax12.set_xlim(1,1.6)
    ax12.set_ylim(-2.5,0)
    
    ax1.legend()
    ax2.legend()
    fig1.savefig("Potential1.png")
    fig2.savefig("Potential2.png")
    
    
    
def main():
    
    rmin=0.1
    rc=1.6
    n_points=2000
    n=4
    
    f=open("frenkel.dat",'w')
    f,data1=write_potential(f,2.5,1.0,n,rc,n_points,rmin)
    f,data2=write_potential(f,1.0,1.0,n,rc,n_points,rmin)
    f.close()
    
    run_lammps(rmin,rc,n_points)
    
    plot_test(data1,data2)
    
    run_lammps_hybrid(rmin,rc,n_points)
    print ("\n****************************************************************")
    print ("Always check that the generated potential in Lammps is accurate!")
    print ("****************************************************************")
    
    
main()
