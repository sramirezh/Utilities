#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:16:38 2020
Copied launch_restart for automatation 
@author: sr802
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell
Utilities_path=os.path.join(os.path.dirname(__file__), '../../')
sys.path.append(Utilities_path) #This falls into Utilities path
import Lammps.core_functions as cf
import shutil
from Lammps.Pore.qsub.simulation_utilities import simulation
import argparse

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script


#class atom():

# TODO abstract this in a class
def create_init(box_limits,T, x , v, n = 2, ntypes = 2):
    """
    Creates the input file that LAMMPS read with read_data following Lammps guidelines in https://lammps.sandia.gov/doc/read_data.html
    
    Notice that The first molecule is centered and has zero velocity
    
    Args:
        box_limits vector containing the box dimensions [xlo, xhi, ylo, yhi, zlo, zhi]
        T temperature to get velocities from the Maxwell Boltzmann distribution
        x relative position [Vector]
        v relative velocitie [Vector]
        n number of atoms
    """
    
    f = open("initial.dat",'w')
    f.write("# This file was created using my own algorithm in the Transport coefficients project\n\n")
    f.write("%s atoms\n"%n)
    f.write("%s atom types\n"%ntypes)
    f.write('\n')
    
    #Box limits
    f.write("%s %s xlo xhi\n"%(box_limits[0], box_limits[1]))
    f.write("%s %s ylo yhi\n"%(box_limits[2], box_limits[3]))
    f.write("%s %s zlo zhi\n"%(box_limits[4], box_limits[5]))
    f.write('\n')
    
    # Masses
    f.write('Masses\n')
    f.write('\n')
    f.write('1 1\n')
    f.write('2 1\n')
    f.write('\n')
    
    #Atom (Check the atom section in lammps)
    # e.g atomic type atom-ID atom-type x y z
    f.write('Atoms\n')
    f.write('\n')
    f.write('1 1 0 0 0\n')
    f.write('2 2 %s %s %s\n'%(x[0], x[1], x[2]))
    f.write('\n')
    
    #Velocities
    f.write('Velocities\n')
    f.write('\n')
    f.write('1 0 0 0\n')
    f.write('2 %s %s %s\n'%(v[0], v[1], v[2]))
    f.write('\n')
    
    f.close()
    
    


##Main
#    
#T = 1
#m = [1,1] # Vector with all the masses of the particles
#
#m_red =  m[0]*m[1]/(m[0]+m[1]) #Reduced mass
#scale = np.sqrt(T/m_red)
#
#b = -6 #Impact parameter
#vx = maxwell.rvs(size = 1, scale = scale)[0]
#
#v = [vx, 0 , 0]
#x = [b, 0, 0]
#
#
#L = 30 # Half lenght of the simulation box
#create_init([-L, L, -L, L, -L, L,], 1, x, v)
#    





# =============================================================================
# Main
# =============================================================================
def main(name, root, template, n_vel , Temperature , b, run):       
    #Getting the path to all the restart files
    #files=glob.glob('%s/*'%conf_folder)
    
    home=root+'/'+name
    
    # =============================================================================
    #     Additional calculations
    # =============================================================================
    L = 30 # Half lenght of the simulation box
    m = [1,1]
    m_red =  m[0]*m[1]/(m[0]+m[1]) #Reduced mass
    x = [b, 0, 0]
    scale = np.sqrt(Temperature/m_red)
    velocities =  maxwell.rvs(size = n_vel, scale = scale)
    
    for i,vel in enumerate(velocities):
          
    # =============================================================================
    #     Creating the simulation instance 
    # =============================================================================
        
        sim = simulation(home,template,i)
        sim.create_folder()
        sim.create_qsub('short',1,16,24,'input.lmp')
    # =============================================================================
    #     #Mofications to the files here (THIS IS SPECIFIC)
    # =============================================================================
        os.chdir(sim.folder)
        v = [vel, 0 , 0]
        create_init([-L, L, -L, L, -L, L,], 1, x, v)
        os.chdir(cwd)
        
        
    # =============================================================================
    #     Running the simulation
    # =============================================================================
        if run == True:
            sim.run_simulation()
        os.chdir(cwd)
            
if __name__ == "__main__":
    """
    THIS IS VERY SPECIFIC
    The arguments of this depend on the application
    """
    parser = argparse.ArgumentParser(description='Launch simulations from restart',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-name', metavar='name',help='Name of the folder to keep all the simulations',default='b')
    parser.add_argument('-template', metavar='path_template',help='Directory to take as template',default=cwd+'/Template')
    parser.add_argument('-root', metavar='root directory',help='Directory to create the folder for the simulations',default=cwd)
    parser.add_argument('-n_vel',metavar='n conf',help='number of velocities',default=5,type = int)
    parser.add_argument('-Temperature',metavar='Temperature',help='Temperature of the system',default=1.0,type=float)
    parser.add_argument('-b',metavar='b',help='inpact parameter', default = 6.0, type=float)
    parser.add_argument('-run',metavar='run',help='Define if run simulations or not. If not, just creates the folder structure',default = False,type=cf.str2bool)
    
    args = parser.parse_args()

    
    main(args.name,args.root,args.template,args.n_vel,args.Temperature,args.b,args.run)


    