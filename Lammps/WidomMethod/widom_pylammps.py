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

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script


parser = argparse.ArgumentParser(description="This script gets the results created by dp_poly and the averages of vdata.dat, computes relevant quantities and generates plots, It has to be run inside every N_X",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-folder_names',help='Folders where the analysis is going to run i.e E_*',default=glob.glob('E_*'),nargs='+')
    parser.add_argument('-s','--source',choices=['read','run','gather'], default="read",help='Decides if the if the file Statistics_summary.dat needs to be read, run, gather  ')
    parser.add_argument('--vdatamin', help='Number of samples to be discarded in vdata.dat', default=1000, type=int)
    parser.add_argument('--dpolymin', help='Number of samples to be discarded in DPpoly', default=100, type=int)


def getting_properties(file_name):
    """
    Returns the intial energy of the system and the box size
    """
    box=[]
    L=IPyLammps()
    L.command("read_data %s"%file_name) #Not necessary to create a new type of particles
    L.pair_style("lj/cut", 3.0)
    L.pair_coeff(1, 1, 1.0, 1.0)
    L.command("pair_modify     tail no")
#    L.command("pair_modify     shift yes")
    L.command("mass * 1.0")
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
    
    return ene_ini,box,natoms

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

# =============================================================================
# Main
# =============================================================================
t=time.time()

file_name="final_conf_-4.0.dat"
ene_i,box,natoms=getting_properties(file_name)
temperature=2.0
beta=1/temperature
vol=(box[1]-box[0])*(box[3]-box[2])*(box[5]-box[4])
rho=natoms/vol    #Has to be per species

mu_id=mu_id(temperature,rho )

n_trials=1
Boltzmann=np.zeros(n_trials)



for i in xrange(n_trials):
    L=IPyLammps()
    
    L.atom_style("atomic") #Necessary to create atom maps
    L.command("atom_modify map yes") 
    L.command("read_data %s"%file_name) #Not necessary to create a new type of particles
    L.pair_style("lj/cut",3.0)
    L.pair_coeff(1, 1, 1.0, 1.0)
    L.command("pair_modify     tail no")
#    L.command("pair_modify     shift yes")
    L.command("mass * 1.0")
    #L.command("include in.interaction")   
    L.create_atoms(1, "single %f %f %f"%(random_position(box)))
    L.run(0)
    ene_f=L.eval('pe')  #Final energy
    Boltzmann[i]=np.exp(-beta*((natoms+1)*ene_f-natoms*ene_i))
    L.__del__()
    
    
print time.time()-t


ave_bol=np.average(Boltzmann)
mu_ex=-temperature*np.log(ave_bol)
mu=mu_ex+mu_id
print mu,mu_ex,mu_id


L=IPyLammps()
L.atom_style("atomic") #Necessary to create atom maps
L.command("atom_modify map yes") 
L.command("read_data %s"%file_name) #Not necessary to create a new type of particles
L.command("include in.interaction")
L.run(0)

ene_i=L.eval('pe')
L.create_atoms(2, "single %f %f %f"%(random_position(box)))

print "single %f %f %f"%(random_position(box))
L.run(0)
ene_f=L.eval('pe')
delta_energy=ene_f-ene_i
boltzmann=np.exp(-beta*delta_energy)
mu_ex=-temperature*np.log(boltzmann)


print  ene_f,ene_i,delta_energy,mu_ex

print np.exp(-0.5*(ene_f-ene_i))

L.image()


#Try to paralellize
#[x for x in dir(L.atoms[0]) if not x.startswith('__')] this tells me the properties of the atoms