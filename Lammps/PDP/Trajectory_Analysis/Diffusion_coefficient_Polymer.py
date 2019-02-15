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

from joblib import Parallel, delayed
import multiprocessing

def compute_one_delta(pos):
    """
    Computes the g(r) of the polymer particles with other particle types, so only 3
    three distri
    """


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
axis_font=24
tick_font=20
legend_font=18
xoffset=0.05
yoffset=0.8
error_cap=4

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
propD=[]
for i in xrange(3000):
    msd_t=compute_one_delta(pos[::(i+1)])[3]
    msd.append(msd_t)
    dt=times[i]
    t.append(dt)
    propD.append(msd_t/dt/(2*3))
    

out=np.polyfit(t,msd,1)
D=out[0]/(2*3)
print "The diffusion coefficient is %s"%D

plt.close('all')
fig1,(ax1,ax12)=plt.subplots(2,1)
ax1.plot(t,msd)
ax1.tick_params(labelsize=tick_font,direction='in',top=True, right=True)
ax1.set_ylabel(r'$MSD$',fontsize=axis_font)
ax12.plot(t,propD)
ax12.set_xlabel(r'$\Delta t$',fontsize=axis_font)
ax12.set_ylabel(r'$D$',fontsize=axis_font)
ax12.tick_params(labelsize=tick_font,direction='in',top=True, right=True)
plt.show()

