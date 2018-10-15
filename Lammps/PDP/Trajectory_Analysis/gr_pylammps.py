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




cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

imin=0 #Number of configurations to skip
nbins=100
rmax=2.5





"""
*******************************************************************************
Main
*******************************************************************************
"""


t=time.time()
from lammps import IPyLammps

L=IPyLammps()
L.region("box block", 0, 20, 0, 1.9469999999999999e+01, 0, 30)
L.create_box(3, "box")
L.command("create_atoms 1 random 8292 123823 NULL")

L.command("read_dump dumpfile320000.gz 320000  x y z ix iy iz vx vy vz box yes add yes")
L.command("mass * 1.0")
L.pair_style("lj/cut", 2.5)
L.command("pair_coeff * * 1.0 1.0 2.5")

L.command("compute gr all rdf 100 1 1")
L.fix(" 2 all ave/time 1 1 1 c_gr[*] file gr.dat mode vector")
L.write_script("input.lmp")
L.run(0)

input_files = glob.glob("*1*.gz")
input_files.sort(key=lambda f: int(filter(str.isdigit, f)))
times=cf.extract_digits(input_files)
num_conf=len(times)


input_files.sort(key=lambda f: int(filter(str.isdigit, f)))
times=cf.extract_digits(input_files)
num_conf=len(times)

Data=pd.read_csv("gr.dat",sep=" ",skiprows=4,dtype=np.float64,header=None).values

plt.plot(Data[:,1],Data[:,2],label="lammps")
plt.legend()
plt.show()


print time.time()-t
#num_cores = multiprocessing.cpu_count()
#
#results=Parallel(n_jobs=num_cores,verbose=10)(delayed(compute_one_configuration)(fil) for fil in input_files)


#results=[]
#
#for fil in input_files:
#    
#    
#
#
#g_r=np.average(results,axis=0)



#g_poly_solvent=g_poly_solvent/num_conf
# 
#import matplotlib.pyplot as plt
#
#plt.plot(g_r[0][:,0],g_r[0][:,2],label="poly-solvent")
#plt.plot(g_r[1][:,0],g_r[1][:,2],label="poly-solute")
#plt.plot(g_r[2][:,0],g_r[2][:,2],label="poly-poly")
#
#plt.legend()
#
#


#Data=pd.read_csv("gr.dat",sep=" ",skiprows=4,dtype=np.float64,header=None).values
#
#plt.plot(Data[:,1],Data[:,2],label="lammps")
#plt.legend()
#plt.show()




