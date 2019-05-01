#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 15:15:12 2018
This scripts reads the pos.dat that contains the cm positions of a polymer and computes the diffusion coefficient 
@author: sr802
"""

from __future__ import division
import numpy as np
import pandas as pd
import os
import sys
import glob
import poly_analysis as pa
import matplotlib.pyplot as plt
import time
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

def compute_one_msd(pos,delta):
    """
    Computes the MSD for the positions every certain delta
    Args:
        pos: all the positions of the cm of the polymer
        delta: evert this number, we take the positions to compute the msd
    """

    pos=pos[::delta]
    delta_sqr_components=(pos-np.roll(pos,-1,axis=0))**2
    delta_sqr=np.sum(delta_sqr_components,axis=1)[:-1] #The last contribution is the last-the initial msd
    msd_comp=np.average(delta_sqr_components,axis=0)
    msd=np.average(delta_sqr)

    return np.hstack([msd_comp,msd])

"""
*******************************************************************************
Main
*******************************************************************************
"""
cf.set_plot_appearance()

delta_t= 0.005

print "Assuming that the time step is %s"%delta_t
#times_l,msd_l = lammps_MSD(delta_t)

input_file = "pos.dat"
Data=cf.read_data_file(input_file)
data1=Data.values
times=data1[:,0]*delta_t

pos=data1[:,1::]

msd=[]
t=[]
max_delta=int(len(times)*0.05) #Maximum delta of time to measure the MSD
propD=[]
for i in xrange(max_delta):
    msd_t=compute_one_msd(pos,i+1)[3]
    msd.append(msd_t)
    dt=times[i]
    t.append(dt)
    propD.append(msd_t/dt/(2*3))
    

out=np.polyfit(t,msd,1)
fit_line=np.polyval(out,t)
D=out[0]/(2*3)
print "The diffusion coefficient is %s"%D



plt.close('all')
fig1,(ax1,ax12)=plt.subplots(2,1)
ax1.plot(t,msd)
ax1.plot(t,fit_line,'--')
ax1.set_ylabel(r'$MSD$')
ax12.plot(t,propD)
ax12.axhline(y=D, xmin=0, xmax=1,ls='--',c='black')
ax12.set_xlabel(r'$\Delta t$')
ax12.set_ylabel(r'$D$')
plt.tight_layout()
plt.savefig("Diffusio_coefficient.pdf")
plt.show()

