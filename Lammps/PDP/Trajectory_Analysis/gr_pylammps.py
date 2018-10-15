#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 15:37:04 2018
This scripts analyses the pair correlation function from the *.gz files
@author: simon
"""


from __future__ import division
import numpy as np
import pandas as pd
import argparse
import os
import sys
import glob
import poly_analysis as pa
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

from scipy.spatial.distance import pdist,squareform
import time

from joblib import Parallel, delayed
import multiprocessing
import matplotlib.pyplot as plt
from lammps import IPyLammps



cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script



def call_pylammps(fil,p_types):
    """
    Calls pylammps to read the configuration and compute the g(r)
    fil is the name of the dump file
    p_types is an array with the particle types.
    """
    
    file_time=int(filter(str.isdigit,fil))
    
    L=IPyLammps()
    L.region("box block", 0, 20, 0, 1.9469999999999999e+01, 0, 30)
    L.create_box(len(p_types), "box")
    
    L.command("read_dump %s %s x y z ix iy iz vx vy vz box yes add yes"%(fil,file_time))
    L.command("mass * 1.0")
    L.pair_style("lj/cut", 2.5)
    L.command("pair_coeff * * 1.0 1.0 2.5")
    
    
    for i in p_types:
        for j in p_types:
            if j>=i:
                compute_name="gr_%s%s"%(i,j)
                fix_name="f_%s" %compute_name
                L.command("compute %s all rdf %s %s %s"%(compute_name,nbins,i,j))
                L.fix("%s all ave/time 1 1 1 c_%s[*] file %s_%s.dat mode vector" %(fix_name,compute_name,compute_name,file_time))
    
    L.run(0)
    L.close()
    

def read_gr(fil, p_types):
    """
    Gathers all the g(r) from the different files for a given timestep
    """
    
    results=[]
    file_time=int(filter(str.isdigit,fil))
    for i in p_types:
        for j in p_types:
            if j>=i:
                compute_name="gr_%s%s"%(i,j)
                file_name="%s_%s.dat" %(compute_name,file_time)
                Data=pd.read_csv(file_name,sep=" ",skiprows=4,dtype=np.float64,header=None).values
                os.remove(file_name)
                results.append(Data)
    
    return results



def run_one_time(fil,p_types):
    """
    Runs everything for a sigle time step
    args:
        fil: is the filename
        p_types is an array containing the species
            
    """
    call_pylammps(fil,p_types)
    results=read_gr(fil, p_types)
    
    return results
    
"""
*******************************************************************************
Main
*******************************************************************************
"""

imin=0 #Number of configurations to skip
nbins=100
rmax=2.5

p_types=[1,2,3]
    
input_files = glob.glob("*.gz")
input_files.sort(key=lambda f: int(filter(str.isdigit, f)))
times=cf.extract_digits(input_files)
num_conf=len(times)

# %%

t=time.time()


num_cores = multiprocessing.cpu_count()
results=Parallel(n_jobs=num_cores,verbose=10)(delayed(run_one_time)(fil,p_types) for fil in input_files)


g_r=np.average(results,axis=0)
print (time.time()-t)

            

# %% Plots

plt.plot(g_r[0][:,1],g_r[0][:,2],label="solvent-solvent")
plt.plot(g_r[1][:,1],g_r[1][:,2],label="solvent-solute")
plt.plot(g_r[2][:,1],g_r[2][:,2],label="solvent-poly")
plt.plot(g_r[3][:,1],g_r[3][:,2],label="solute-solute")
plt.plot(g_r[4][:,1],g_r[4][:,2],label="solute-poly")

plt.legend()

plt.savefig("g_r.pdf")





