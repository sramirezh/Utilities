#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 16:16:17 2019
Widom insertion method using pylammps,
be careful that 'pe' gives the potential energy per particle!
@author: sr802
"""


from __future__ import division
import numpy as np



from lammps import IPyLammps

import pandas as pd
import argparse
import os
import sys
import time
import random
import re

sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

def print_lvar(L,var):
    """
    Prints a lammps variable
    Args:
        L pylamms instance
        var variable name as defined in the lammps instance
    Return:
        x a float containing the value of the variable.
    """
    
    if sys.version_info < (3, 0):
        # In Python 2 'print' is a restricted keyword, which is why you have to use the lmp_print function instead.
        x = float(L.lmp_print('"${%s}"'%var))
    else:
        # In Python 3 the print function can be redefined.
        # x = float(L.print('"${a}"')")
    
        # To avoid a syntax error in Python 2 executions of this notebook, this line is packed into an eval statement
        x = float(eval("L.print('\"${$s}\"')"%var))
    
    return x

def read_interactions(interaction_file):
    """
    reads the relevant lines from the in.interaction file 
    
    Returns:
        interaction file lines
    """
    
    f=open(interaction_file,"r")
    #avoiding commented or blank lines
    interactions=[]
    for line in f:
        if (line.startswith('#')) or (len(line.strip())==0):pass        
        else:
            line=line.split("#")[0].strip('\n')  #To remove and end of lines
            interactions.append(line.replace("\t",''))
    return interactions


def include_interactions(L,interacions):
   """
   includes the interactions in the Lammps object
       Args:
           L the lammps instance 
           interactions are the lines from the interaction file
           
       Returns: 
           L the lammps objects with the loaded interactions
   """
   for line in interactions:
       L.command("%s"%line)
      
        
   return 

def getting_properties(file_name,interactions):
    """
    Returns the intial energy of the system and the box size
    """
    box=[]
    L=IPyLammps()
    L.command("read_data %s"%file_name) #Not necessary to create a new type of particles
    include_interactions(L,interactions)
    box.append(L.system.xlo)
    box.append(L.system.xhi)
    box.append(L.system.ylo)
    box.append(L.system.yhi)
    box.append(L.system.zlo)
    box.append(L.system.zhi)
    L.run(0)
    ene_ini=L.eval('pe')
    #This has to be per species
    natoms=L.atoms.natoms
    L.__del__()
    vol=(box[1]-box[0])*(box[3]-box[2])*(box[5]-box[4])
    return ene_ini,box,natoms,vol

def random_position(box):
    """
    Returns a random position inside the box
    """
    x=random.random()*(box[1]-box[0])+box[0]
    y=random.random()*(box[3]-box[2])+box[0]
    z=random.random()*(box[5]-box[4])+box[0]
    
    return x,y,z

def mu_id(temp,rho):
    return temp*np.log(rho)

def main():
    
    # =============================================================================
    # Main
    # =============================================================================
        
    
    parser = argparse.ArgumentParser(description="This script computes the chemical potential for a given configuration",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file', metavar='input_file',help='configuration file generated in Lammps via write_data',type=lambda x: cf.is_valid_file(parser, x))
    parser.add_argument('interaction_file', metavar='interaction_file',help='interaction file lammps style',type=lambda x: cf.is_valid_file(parser, x))
    parser.add_argument('-temp',help='Temperature of the simulation',default=2.0,type=float)
    parser.add_argument('-nin',help='Number of particles insertions',default=100,type=int)
    
    args = parser.parse_args()
    file_name=args.input_file
    interaction_file=args.interaction_file
    
    
    interactions=read_interactions(interaction_file)
    ene_i,box,natoms,vol=getting_properties(file_name,interactions)
    temperature=args.temp
    n_trials=args.nin
    
    #Basic computations
    beta=1/temperature
    rho=natoms/vol    #Has to be per species
    mu_id=mu_id(temperature,rho )
    
    
    #Particle insertion
    
    Boltzmann=np.zeros(n_trials)
    for i in xrange(n_trials):
        L=IPyLammps()
        
        L.atom_style("atomic") #Necessary to create atom maps
        L.command("atom_modify map yes") 
        L.command("read_data %s"%file_name) #Not necessary to create a new type of particles
        include_interactions(L,interactions)  
        L.create_atoms(1, "single %f %f %f"%(random_position(box)))
        L.run(0)
        ene_f=L.eval('pe')  #Final energy
        Boltzmann[i]=np.exp(-beta*((natoms+1)*ene_f-natoms*ene_i))
        L.__del__()
        
        
    
    
    ave_bol=np.average(Boltzmann)
    mu_ex=-temperature*np.log(ave_bol)
    mu=mu_ex+mu_id
    print mu,mu_ex,mu_id



#def test():

file_name="final_conf_2.0_binary.dat"

interaction_file="in.interaction_binary"
interactions=read_interactions(interaction_file)
#ene_i,box,natoms,vol=getting_properties(file_name,interactions)



L=IPyLammps()
L.atom_style("atomic") #Necessary to create atom maps
L.command("atom_modify map yes") 
L.command("read_data %s"%file_name) #Not necessary to create a new type of particles
include_interactions(L, interactions)
L.run(0)

#
#ene_i=L.eval('pe')
#L.create_atoms(1, "single %f %f %f"%(random_position(box)))
#
#print "single %f %f %f"%(random_position(box))
#L.run(0)
#ene_f=L.eval('pe')
#delta_energy=ene_f-ene_i
#boltzmann=np.exp(-beta*delta_energy)
#mu_ex=-temperature*np.log(boltzmann)
#
#
#print  ene_f,ene_i,delta_energy,mu_ex
#
#print np.exp(-0.5*(ene_f-ene_i))
#
#L.image()
