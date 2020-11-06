#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 15:07:55 2019

@author: sr802
"""

import numpy as np
import multiprocessing
import os
import sys
import time
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.Pore.qsub.simulation_results as sr
import Lammps.lammps_utilities as lu
import glob
from tqdm import tqdm  
from joblib import Parallel, delayed
import Others.Statistics.FastAverager as stat

from Others.Statistics.Functions import autocorrelation_error, blocking_error, blocking


from uncertainties import unumpy,ufloat


data=cf.read_data_file('vdata.dat')


num_cores = multiprocessing.cpu_count()

data1 = data.values
delta_t = 0.5
times=(data1[:,0]-data1[0,0])*delta_t
max_delta=int(len(times)*0.004) #Maximum delta of time to measure the correlation
var1 = data1[:,1]
var2 = data1[:,1]


#def compute_correlation_dt(var1,var2,delta):
#    """
#    This is a VERY GENERAL function
#    
#    *****
#    It is important to keep it outside the class as it is going to be evaluated in parallel and if it is a method of the class, it would require
#    to load the instance on each processor which could be very expensive
#    *****
#    
#    Computes the correlation for two variables for a given delta t
#    Args:
#        var1: time series for the first variable 
#        var2: time series for the first variable
#        Note that var1 and var2 have to have the same the same time series
#        delta: every this number of steps, we take the variable 2 to compute the correlation
#    """
#    cf.blockPrint()
#    if delta != 0:
#        var1 = var1[::delta]
#        var2 = var2[::delta]
#        cor = (var1*np.roll(var2,-1,axis=0))
#        cor = cor[:-1]#The last contribution is the last-the initial msd
#    else:
#        cor = var1*var2 
#
#    
#    average = stat.fast_averager(cor)[0]
#    
#    cf.enablePrint()
#    
#    return ufloat(average[1],average[2]) 

ti = time.time()
cor1 = var1*var2 

average1 = stat.fast_averager(cor1)[0]

print ("The final time without the multiplication is %f"%(time.time()-ti))



ti = time.time()
cor2 = var1*var2*4699**2 

average2 = stat.fast_averager(cor2)[0]

print ("The final time with the multiplication is %f"%(time.time()-ti))

#np.concatenate(var1,var1)
#
#var1 = data=np.reshape(var1,(len(var1),1))
#c,ct = autocorrelation_error(np.reshape(var1,(1,len(var1))),[np.average(var1)])

w = np.transpose(np.stack((cor1,cor1)))

c,ct = autocorrelation_error(w)


ti = time.time()
B = blocking(w)
print ("The final time of one blocking is is %f"%(time.time()-ti))


data2=cf.read_data_file('data1.txt')

u1= data2.values

final_error,Error = blocking_error(u1,True)



Correlation,time = autocorrelation_error(u1)
error_c = np.sqrt(Correlation[0,:]*2*time/(len(u1)+1))
print (error_c)