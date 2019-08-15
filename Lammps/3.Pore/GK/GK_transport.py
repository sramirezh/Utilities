#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:55:09 2019
This scripts computes the correlation between fluxes
@author: sr802
"""

from __future__ import division
import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
from scipy import optimize
import warnings
warnings.filterwarnings("ignore")
import argparse
from tqdm import tqdm  
from joblib import Parallel, delayed
import multiprocessing

try:
    from uncertainties import unumpy,ufloat
except ImportError as err2:
    print err2
import Others.Statistics.FastAverager as stat


def compute_correlation(var1,var2,delta):
    """
    This is a VERY GENERAL function
    
    Computes the MSD the correlation for two variables every certain delta t
    Args:
        var1: time series for the first variable 
        var2: time series for the first variable
        Note that var1 and var2 have to have the same the same time series
        delta: every this number of steps, we take the variable 2 to compute the correlation
    """
    cf.blockPrint()
    global average
    
    if delta!=0:
        var1=var1[::delta]
        var2=var2[::delta]
        correlation=(var1*np.roll(var2,-1,axis=0))
        correlation=correlation[:-1]#The last contribution is the last-the initial msd
    else:
        correlation=var1*var2
    
    average=stat.fast_averager(correlation)[0]
    cf.enablePrint()
        
    return ufloat(average[1],average[2]) 
    


class flux(object):
    """
    Flux is a vectorial or scalar entity
    """
    def __init__(self,components,times,name):
        """
        Args:
            Components: is a matrix( or vector) containing the time series of the components in [x,y,z]
            Name: Is the name that is going to appear in the plots in latex format, example "r'J_s-c_s^BQ'"
        """
        self.components = components
        self.dimension =np.shape(components)[1]
        self.name = name
        self.times = times
        
        
class correlation(object):
    def __init__(self,flux1,flux2, max_delta):
        self.flux1 = flux1
        self.flux2 = flux2
        
    def initial_check():
        """
        Checks if the time series are equal and the fluxes have the same number of components
        """
        
    

"""
Computes the diffusion coefficient
"""
#Input parameters
input_file='fluxes.dat'
delta_t=0.005

Data=cf.read_data_file(input_file)
data1=Data.values
times=(data1[:,0]-data1[0,0])*delta_t

total_flux=flux(data1[:,[1,3,5]],times,"r'Q'")
Solute_excess=flux(data1[:,[2,4,6]],times,"r'J_s-c_s^BQ'")





#var1=data1[:,1]
#var2=data1[:,2]
#
#t=[]
#max_delta=int(len(times)*0.004) #Maximum delta of time to measure the correlation
#
#
#num_cores = multiprocessing.cpu_count()
#
##correlation=compute_correlation(var1,var2,0)
#
#correlation=Parallel(n_jobs=num_cores)(delayed(compute_correlation)(var1,var1,i) for i in tqdm(xrange(max_delta)))
#
#correlation_initial=compute_correlation(var1,var1,0) #To normalise
#
#correlation_norm_x=np.array(correlation)/correlation_initial.nominal_value
#
#
#
#var1=data1[:,3]
#var2=data1[:,4]
#
#t=[]
#max_delta=int(len(times)*0.004) #Maximum delta of time to measure the correlation
#
#
#num_cores = multiprocessing.cpu_count()
#
##correlation=compute_correlation(var1,var2,0)
#
#correlation=Parallel(n_jobs=num_cores)(delayed(compute_correlation)(var1,var1,i) for i in tqdm(xrange(max_delta)))
#
#correlation_initial=compute_correlation(var1,var1,0) #To normalise
#
#correlation_norm_y=np.array(correlation)/correlation_initial.nominal_value
#
#
#var1=data1[:,5]
#var2=data1[:,6]
#
#t=[]
#max_delta=int(len(times)*0.004) #Maximum delta of time to measure the correlation
#
#
#num_cores = multiprocessing.cpu_count()
#
##correlation=compute_correlation(var1,var2,0)
#
#correlation=Parallel(n_jobs=num_cores)(delayed(compute_correlation)(var1,var1,i) for i in tqdm(xrange(max_delta)))
#
#correlation_initial=compute_correlation(var1,var1,0) #To normalise
#
#correlation_norm_z=np.array(correlation)/correlation_initial.nominal_value
#
#
##D_inst=[] #Array with the instantaneous diffusion coefficient
#
#Total=(correlation_norm_z+correlation_norm_y+correlation_norm_x)/3
#
#
#for i in xrange(max_delta):
#    dt=times[i]
#    t.append(dt)
#
## Put this inside the argparser
#try:
#    import matplotlib
#    matplotlib.use('agg')
#    import matplotlib.pyplot as plt
#except ImportError as err:
#    print err
#cf.set_plot_appearance()
#plt.close('all')
#fig,ax=plt.subplots()
#
##For x
#y=np.array([i.n for i in correlation_norm_x])
#y_error=np.array([i.s for i in correlation_norm_x])
#
#ax.plot(t,y)
#ax.fill_between(t, y-y_error, y+y_error ,alpha=0.4)
#
##For y
#y=np.array([i.n for i in correlation_norm_y])
#y_error=np.array([i.s for i in correlation_norm_y])
#
#ax.plot(t,y)
#ax.fill_between(t, y-y_error, y+y_error ,alpha=0.4)
#
##For z
#y=np.array([i.n for i in correlation_norm_z])
#y_error=np.array([i.s for i in correlation_norm_z])
#
#ax.plot(t,y)
#ax.fill_between(t, y-y_error, y+y_error ,alpha=0.4)
#
#
##For total
#y=np.array([i.n for i in Total])
#y_error=np.array([i.s for i in Total])
#
#ax.plot(t,y,color='w')
#ax.fill_between(t, y-y_error, y+y_error ,alpha=0.4,color='w')
#
#ax.set_xscale('log')
#
#plt.savefig("correlation.pdf")



#
#
#
##Writing arrays of averages and errors
#t=np.array(t)
#msd_error=unumpy.std_devs(msd)
#msd_average=unumpy.nominal_values(msd)
#
#
#D_inst_error=unumpy.std_devs(D_inst)
#D_inst_ave=unumpy.nominal_values(D_inst)
#pfinal,cov=fit_line(t,msd_average,msd_error,initial_index,final_ratio,step)
#D=pfinal[0]/(2*3)
#D_err=np.sqrt(cov[0][0])*D
#
#plot_diffusion(t,msd_average,msd_error,D_inst_ave,D_inst_error,pfinal,D)
#
#print "The diffusion coefficient is %s +/- %s"%(D,D_err)
#f=open("Diffusion.out",'w')
#f.write("The diffusion coefficient is %s +/- %s \n"%(D,D_err))
#f.close

