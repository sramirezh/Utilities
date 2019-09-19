#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 15:15:12 2018
This scripts reads the MSD generated from lammps and also computes it from the trajectories
to find the diffusion coefficient for both cases and validating the approach.
@author: sr802
"""


import numpy as np
import pandas as pd
import argparse
import os
import sys
import glob
from . import poly_analysis as pa
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

from joblib import Parallel, delayed
import multiprocessing

def lammps_MSD(delta_t):


    data=pd.read_csv("Parameters.dat",sep=" ",dtype=np.float64).values[:,:-1]

    times=data[:,0]-data[0,0]
    times=times*delta_t
    msd=data[:,-1]

    fig,ax=plt.subplots()

    ax.plot(times,msd,label="LAMMPS")

    out=np.polyfit(times,msd,1)

    ax.plot(times,out[0]*times,label="fit")

    D=out[0]/(2*3)

    print("The diffusion coefficient from Lammps MSD is %s"%D)
    return times,msd

def compute_one_time(pos_init,fil):
    """
    Computes the g(r) of the polymer particles with other particle types, so only 3
    three distri
    """

    #Try to change this for cf.read_data_file
    Data=pd.read_csv(fil,sep=" ",skiprows=9,dtype=np.float64,header=None).sort_values(0).values[:,:-1]
    pos=pa.real_position(Data,L) #Real positions of all the atoms

    delta_sqr_components=(pos-pos_init)**2
    delta_sqr=np.sum(delta_sqr_components,axis=1)

    msd_comp=np.average(delta_sqr_components,axis=0)
    msd=np.average(delta_sqr)

    return np.hstack([msd_comp,msd])

"""
*******************************************************************************
Main
*******************************************************************************
"""


delta_t= 0.005

#times_l,msd_l = lammps_MSD(delta_t)

input_files = glob.glob("conf/dump*")
input_files.sort(key=lambda f: int(list(filter(str.isdigit, f))))
times=cf.extract_digits(input_files)*delta_t
num_conf=len(times)
Box,L=pa.Box_limits(input_files[0])

Data_init=pd.read_csv(input_files[0],sep=" ",skiprows=9,dtype=np.float64,header=None).sort_values(0).values[:,:-1]
pos_init=pa.real_position(Data_init,L) #Real positions of all the atoms



num_cores = multiprocessing.cpu_count()
results=Parallel(n_jobs=num_cores,verbose=10)(delayed(compute_one_time)(pos_init,fil) for fil in input_files)


#SERIAL solution
#results=[]
#for i,fil in enumerate(input_files):
#    results.append(compute_one_time(pos_init,fil))

results=np.array(results)
#plt.plot(times_l,msd_l,label="lammps")
#plt.plot(times,results[:,3],label="My_algorithm")

out=np.polyfit(times,results[:,3],1)

D=out[0]/(2*3)

print("The diffusion coefficient from My calculations is %s"%D)
