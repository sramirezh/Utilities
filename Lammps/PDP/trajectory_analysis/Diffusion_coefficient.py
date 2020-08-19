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
import os
import sys
import glob
import poly_analysis as pa
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.General.thermo_analyser as ta


cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

from joblib import Parallel, delayed
import multiprocessing

def lammps_MSD(delta_t, logger):
    """
    Computes the diffusion coefficient from Lammps mean-squared displacement(MSD)
    
    Args: 
        delta_t the interval between measurements
    logger: an instance of the logger file to write the output
    
    Returns:
        times and msd 
    """

    ta.thermo_analyser("log.lammps")
    data=pd.read_csv("Parameters.dat",sep=" ",dtype=np.float64).values[:,:-1]

    times=data[:,0]-data[0,0]
    times=times*delta_t
    msd=data[:,-1]

    fig,ax=plt.subplots()

    ax.plot(times,msd,label="LAMMPS")

    out=np.polyfit(times,msd,1)

    ax.plot(times,out[0]*times,label="fit")

    D=out[0]/(2*3)

    logger.info("The diffusion coefficient from Lammps MSD is %s"%D)
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
plot_dir = "plots/diffusion_coefficient"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(), plot_dir)    


delta_t= 0.005

times_l, msd_l = lammps_MSD(delta_t, logger)

input_files = glob.glob("conf/dump*")
input_files.sort(key=lambda f: list(filter(str.isdigit, f)))
times = cf.extract_digits(input_files)[0] * delta_t
times = times - times[0]
num_conf = len(times)
Box,L = pa.Box_limits(input_files[0])

Data_init=pd.read_csv(input_files[0],sep=" ",skiprows=9,dtype=np.float64,header=None).sort_values(0).values[:,:-1]
pos_init=pa.real_position(Data_init,L) #Real positions of all the atoms



num_cores = multiprocessing.cpu_count()
results=Parallel(n_jobs=num_cores,verbose=10)(delayed(compute_one_time)(pos_init,fil) for fil in input_files)


#SERIAL solution
#results=[]
#for i,fil in enumerate(input_files):
#    results.append(compute_one_time(pos_init,fil))

results=np.array(results)

# Need to discard the short times
discard = int(0.5*len(times))
out = np.polyfit(times[discard:],results[discard:,3],1)

D = out[0]/(2*3)

logger.info("The diffusion coefficient from My calculations is %s"%D)


# =============================================================================
# Ploting the MSD
# =============================================================================

cf.set_plot_appearance()

# Ploting the comparison between BD-NEMD, FD-NEMD and Theory
fig, ax = plt.subplots()


ax.plot(times_l, msd_l,label="lammps")
ax.plot(times,results[:,3],label="My algorithm")
ax.plot(times[discard:], np.polyval(out,times[discard:]), label ="my fit")


#ax.errorbar(dcv[:,0], -dcv[:,1], yerr= dcv[:,2], label="BD-NEMD", ls = '--', fmt='o')
#ax.errorbar(fdmd_free[:,0],-fdmd_free[:,1],yerr= fdmd_free[:,2],label="FD-NEMD", fmt='o', ls = '--')
##ax.errorbar(fdmd_fixed[:,0],fdmd_fixed[:,1],yerr= fdmd_fixed[:,2],label="EFMD fixed", fmt='o')
##ax.scatter(theory[:,0], theory[:,1], label="Theory v0")
#ax.plot(theory[:,0], theory[:,2], label="Theory", marker = "s", ls = '--')
##ax.plot(theory_r[:,0], theory_r[:,2], label="Theory $r=R_h$")
#ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
#
#ax.set_xlabel(r'$\varepsilon_{cs} $')
#ax.set_ylabel(r'$v_c^x$')
#ymin,ymax=plt.ylim()
ax.set_ylim(0, None)
ax.set_xlim(0, None)
ax.legend(loc='upper left')
fig.tight_layout()
plt.savefig("%s/msd_comparison.pdf"%plot_dir)
logger.info("plotted %s/msd_comparison.pdf"%plot_dir)



