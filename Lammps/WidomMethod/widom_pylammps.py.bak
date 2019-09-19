#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 16:16:17 2019
 insertion method using pylammps,
be careful that 'pe' gives the potential energy per particle!
For now it just reads a restart file

PENDING:
    
    read_restart could be a variable that could also be called as read_data, depending on the input file
@author: sr802
"""


from __future__ import division
from lammps import IPyLammps
import numpy as np
import argparse
import os
import sys
import random
from joblib import Parallel, delayed
import multiprocessing
import Others.Statistics.FastAverager as stat
from uncertainties import unumpy,ufloat

sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script



# =============================================================================
# Function definition
# =============================================================================

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




def include_interactions(L,interactions):
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


def particle_insertion(initial,atom_type):
    """
    One particle insertion
    Args:
        initial the instance of the system to analyse
        atom_type the atom to insert
    """
    L=IPyLammps()
    L.atom_style("atomic") #Necessary to create atom maps
    L.command("atom_modify map yes") 
    L.command("read_restart %s"%initial.file_name) #Not necessary to create a new type of particles
    include_interactions(L,initial.interactions)  
    L.create_atoms(atom_type, "single %f %f %f"%(random_position(initial.box)))
    L.run(0)
    ene_f=L.eval('pe')  #Final energy
    Boltzmann=np.exp(-initial.beta*((np.sum(initial.natoms)+1)*ene_f-np.sum(initial.natoms)*initial.ene))
    L.__del__()
    return Boltzmann



def _method(initial,n_trials,atom_type):
    num_cores = multiprocessing.cpu_count()
    Boltzmann=Parallel(n_jobs=num_cores,verbose=10)(delayed(particle_insertion)(initial,atom_type) for i in xrange(n_trials))
    return Boltzmann


    

def visualise(initial, name="snapshot.png"):
    """
    visualise the system
    Args: 
        initial: the instance of the system to visualise
    """
    L=IPyLammps()
        
    L.atom_style("atomic") #Necessary to create atom maps
    L.command("atom_modify map yes") 
    L.command("read_restart %s"%initial.file_name) #Not necessary to create a new type of particles
    include_interactions(L,initial.interactions)  
    L.create_atoms(atom_type, "single %f %f %f"%(random_position(initial.box)))
    L.run(0)
    L.image(filename=name)
    print "created a file called %s" %name
    L.__del__()
    

# =============================================================================
# Class definition
# =============================================================================
class system(object):
    
    def __init__(self,file_name,interaction_file,temp):
        self.file_name=file_name
        self.interaction_file=interaction_file
        self.temperature=temp
        self.beta=1/self.temperature
        
    def read_interactions(self):
        """
        reads the relevant lines from the in.interaction file 
        
        Returns:
            interaction file lines
        """
        f=open(self.interaction_file,"r")
        #avoiding commented or blank lines
        interactions=[]
        for line in f:
            if (line.startswith('#')) or (len(line.strip())==0):pass        
            else:
                line=line.split("#")[0].strip('\n')  #To remove and end of lines
                interactions.append(line.replace("\t",''))
        self.interactions=interactions
        
    def get_natoms(self):
        """
        gets the number of atom per species
        """
        L=IPyLammps()
        L.command("read_restart %s"%self.file_name) #Not necessary to create a new type of particles
        include_interactions(L,self.interactions)
        L.run(0)
        self.n_types=L.system.ntypes
        numbers=[]
        for j in xrange(self.n_types):
            i=j+1
            L.command("group g%s type %s"%(i,i))
            L.variable('a%s equal count(g%s)'%(i,i))
            x=print_lvar(L,"a%s"%i)
            numbers.append(int(x))
        self.natoms=numbers
        L.__del__()

    def get_properties(self):
        """
        Returns the intial energy of the system and the box size
        """
        box=[]
        L=IPyLammps()
        L.command("read_restart %s"%self.file_name) #Not necessary to create a new type of particles
        include_interactions(L,self.interactions)
        box.append(L.system.xlo)
        box.append(L.system.xhi)
        box.append(L.system.ylo)
        box.append(L.system.yhi)
        box.append(L.system.zlo)
        box.append(L.system.zhi)
        L.run(0)
        self.ene=L.eval('pe')
        #This has to be per species
        L.__del__()
        self.vol=(box[1]-box[0])*(box[3]-box[2])*(box[5]-box[4])
        self.box=box
        print "\nFound a simulation box with volume=%f, a total of %d particles and %d species"%(self.vol,np.sum(self.natoms),self.n_types)
        
        rho_vector=[]
        mu_id_vector=[]
        for i in xrange(self.n_types):
            rho=(self.natoms[i]/self.vol)
            rho_vector.append(rho)
            mu_id_vector.append(mu_id(self.temperature,rho))
            
        self.rho=rho_vector
        self.mu_id=mu_id_vector
        
        
#def main():
    
    # =============================================================================
    # Main
    # =============================================================================
        
    
parser = argparse.ArgumentParser(description="This script computes the chemical potential for a given configuration",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input_file', metavar='input_file',help='configuration file generated in Lammps via write_data',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('interaction_file', metavar='interaction_file',help='interaction file lammps style',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-type',help='Atom type to analyse as in Lammps',default=1,type=int)
parser.add_argument('-temp',help='Temperature of the simulation',default=2.0,type=float)
parser.add_argument('-nin',help='Number of particles insertions',default=100,type=int)

#Getting the initial parameters

args = parser.parse_args()
file_name=args.input_file
interaction_file=args.interaction_file
temperature=args.temp
n_trials=args.nin
atom_type=args.type



initial_system=system(file_name,interaction_file,temperature)
initial_system.read_interactions()
initial_system.get_natoms()
initial_system.get_properties()



print "\nAnalysing the chemical potential for particles of species %d"%atom_type   
Boltzmann=_method(initial_system,n_trials,atom_type)    



ave_bol=np.array(stat.fast_averager(Boltzmann)[0])

#Standard deviation as in Daan's book
#std_mu_ex=np.std(Boltzmann,ddof=1)/np.sum(Boltzmann)/initial_system.beta**2

#computing the error propagation with ufloats [NOTICE THAT IT IS NOT THE STANDARD DEVIATION]

mu_ex=-temperature*log(ufloat(ave_bol[1],ave_bol[3])) #the average and the error given by the blocking analysis
mu=mu_ex+initial_system.mu_id[atom_type-1]


#Creating the log file
array=[initial_system.rho[atom_type-1],mu,mu_ex.n,initial_system.mu_id[atom_type-1],mu_ex.s]
array=np.array(array).reshape((1,len(array)))
header="Number of particle insertions %s for species %s\n"%(n_trials,atom_type) + "rho mu mu_ex mu_id std_mu_ex"
name="_%s_%s.log"%(n_trials,atom_type)
np.savetxt(name,array,fmt="%3.12f",header=header)

print "Created the file %s with all the information"%name

#main()

#def test():
#
#file_name="restart_-2.0.dat"
#
#interaction_file="in.interaction"
##interactions=read_interactions(interaction_file)
##ene_i,box,natoms,vol=getting_properties(file_name,interactions)
#
#L=IPyLammps()
#L.atom_style("atomic") #Necessary to create atom maps
#L.command("atom_modify map yes") 
#L.command("read_restart %s"%file_name) #Not necessary to create a new type of particles
#L.command("group gSolv type 1")
#L.command("group gSolu type 2")
#L.variable('a equal count(gSolu)')
#if sys.version_info < (3, 0):
#    # In Python 2 'print' is a restricted keyword, which is why you have to use the lmp_print function instead.
#    x = float(L.lmp_print('"${a}"'))
#else:
#    # In Python 3 the print function can be redefined.
#    # x = float(L.print('"${a}"')")
#    
#    # To avoid a syntax error in Python 2 executions of this notebook, this line is packed into an eval statement
#    x = float(eval("L.print('\"${a}\"')"))
#include_interactions(L, interactions)
#L.run(0)

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
