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
import cPickle as pickle

try:
    from uncertainties import unumpy,ufloat
except ImportError as err2:
    print err2
import Others.Statistics.FastAverager as stat


# Put this inside the argparser
try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
except ImportError as err:
    print err


def compute_correlation_dt(var1,var2,delta):
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
        
        #To reshape in order to slice easier in the correlation
        if self.dimension==1: 
            self.components = np.reshape(self.components,(len(self.components),1))
        
        
class correlation(object):
    def __init__(self,flux1,flux2, max_delta):
        self.flux1 = flux1
        self.flux2 = flux2
        self.max_delta = max_delta
        self.initial_check()
        self.cor = (self.dimension+1)*[0]  #list to store the correlations, the last is the total
        self.norm = (self.dimension+1)*[0]  #list to store the normalisation correlation with t=0 (var1(0) var2(0))
        
    def initial_check(self):
        """
        Checks if the time series are equal and the fluxes have the same number of components
        """
        self.dimension=self.flux1.dimension
        if self.flux1.dimension != self.flux2.dimension:
            print "The fluxes do not have the same dimension"
        else:
            self.dimension=self.flux1.dimension
            
        if self.flux1.times.all != self.flux2.times.all:
            print "The fluxes were not measured for the same times"
        else:
            self.times = self.flux1.times[:max_delta]
            
    def compute_correlation_dt(self,var1,var2,delta):
        """
        This is a VERY GENERAL function
        
        Computes the correlation for two variables for a given delta t
        Args:
            var1: time series for the first variable 
            var2: time series for the first variable
            Note that var1 and var2 have to have the same the same time series
            delta: every this number of steps, we take the variable 2 to compute the correlation
        """
        cf.blockPrint()
    
        
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
    
    def correlate_one_d(self,dim):
        """
        Performs a correlation between 1d components by evaluating the products at at different delta t.
        Args:
            dim is the component to be evaluated
        """
        num_cores = multiprocessing.cpu_count()
        var1=self.flux1.components[:,dim]
        var2=self.flux2.components[:,dim]
        max_delta = self.max_delta
        norm = self.compute_correlation_dt(var1,var2,0)
        self.norm[dim] = norm
        cor=Parallel(n_jobs=num_cores)(delayed(compute_correlation_dt)(var1,var2,i) for i in tqdm(xrange(max_delta)))
        self.cor[dim]=np.array(cor)/norm.nominal_value
        
    def evaluate(self):
        """
        Performs the correlations of the 1d components,calluing correlate_one_d, and adds them up to the total.
        
        """
        total = np.zeros(max_delta)
        for dim in xrange(self.dimension):
            self.correlate_one_d(dim)
            total = total + self.cor[dim]
        total = total/3
        self.cor[-1] = total
        self.norm[-1] = 1 #It must be normalised already
        
    def plot_individual(self,fig,ax,dim,alpha=0.4,every=1):
        """
        Args:
            ax axes object
            fig Figure 
            dim is the dimension, for example:in a 3D vector, 0-x, 1-y, 2-z and 3-total.
            alpha is the transparency of the filling
            every to not have so many points
            The axis label is given here but it could be renamed later
        """
        dic_label={0:'x',1:'y',2:'z',self.dimension:'Total'}
        y=np.array([i.n for i in self.cor[dim]])
        y_error=np.array([i.s for i in self.cor[dim]])
        ax.plot(self.times[::every],y[::every],label=dic_label[dim])
        ax.fill_between(self.times, y-y_error, y+y_error ,alpha=0.4)
        
        
        return fig,ax
    
    def plot_all(self,fig,ax,alpha=0.4):
        for dim in xrange(self.dimension+1):
            fig,ax = self.plot_individual(fig,ax,dim)
        
        ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
        return fig,ax
    
    def save(self,file_name):
        """
        Saves the instance 
        """
        afile = open(r'%s.pkl'%file_name, 'wb')
        pickle.dump(self, afile)
        afile.close()    
        
def load_instance(file_name):
    """
    Loads the data structure to be used later
    """
    file1 = open(file_name, 'rb')
    instance = pickle.load(file1)
    file1.close()
    
    return instance



    

##writing the structure
#afile = open(r'p.pkl', 'wb')
#pickle.dump(final_p, afile)
#afile.close()    
        


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
solute_excess=flux(data1[:,[2,4,6]],times,"r'J_s-c_s^BQ'")

max_delta=int(len(times)*0.004) #Maximum delta of time to measure the correlation

#c11=correlation(total_flux,total_flux,max_delta)

#c11.evaluate()
#c11.correlate_one_d(0)
c11 = load_instance("c11.pkl")
c22 = load_instance("c22.pkl")
c12 = load_instance("c12.pkl")
c21 = load_instance("c21.pkl")




plt.close('all')
cf.set_plot_appearance()




# =============================================================================
# Ploting all the correlations
# =============================================================================

# For c11
fig,ax=plt.subplots()

c11.plot_all(fig,ax)
ax.set_xscale('log')
plt.legend(loc='upper right')
plt.savefig("correlation11.pdf")

# For c12
fig,ax=plt.subplots()

c12.plot_all(fig,ax)
ax.set_xscale('log')
plt.legend(loc='upper right')
plt.savefig("correlation12.pdf")

# For c21
fig,ax=plt.subplots()

c21.plot_all(fig,ax)
ax.set_xscale('log')
plt.legend(loc='upper right')
plt.savefig("correlation21.pdf")

# For c22
fig,ax=plt.subplots()

c22.plot_all(fig,ax)
ax.set_xscale('log')
plt.legend(loc='upper right')
plt.savefig("correlation22.pdf")





#
#var1=data1[:,1]
#var2=data1[:,2]
#
#t=[]
#
#
#
#num_cores = multiprocessing.cpu_count()
#
##correlation=compute_correlation(var1,var2,0)
#
#correlation=Parallel(n_jobs=num_cores)(delayed(compute_correlation_dt)(var1,var1,i) for i in tqdm(xrange(max_delta)))
#
#correlation_initial=compute_correlation_dt(var1,var1,0) #To normalise
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

#
#cf.set_plot_appearance()
#plt.close('all')
#fig,ax=plt.subplots()

##For x
#y1=np.array([i.n for i in c11.cor[0]])
#y1_error=np.array([i.s for i in c11.cor[0]])
##
#ax.plot(c11.times,y1)
#ax.fill_between(c11.times, y1-y1_error, y1+y1_error ,alpha=0.4)
#
#
##for the old x
#
#y2=np.array([i.n for i in correlation_norm_x])
#y2_error=np.array([i.s for i in correlation_norm_x])
#
#ax.plot(t,y2)
#ax.fill_between(t, y2-y2_error, y2+y2_error ,alpha=0.4)

##For x
#y=np.array([i.n for i in c11.cor[1]])
#y_error=np.array([i.s for i in c11.cor[1]])
##
#ax.plot(c11.times,y)
#ax.fill_between(t, y-y_error, y+y_error ,alpha=0.4)
#
#
##For x
#y=np.array([i.n for i in c11.cor[2]])
#y_error=np.array([i.s for i in c11.cor[2]])
##
#ax.plot(c11.times,y)
#ax.fill_between(t, y-y_error, y+y_error ,alpha=0.4)
#
###For total
#y=np.array([i.n for i in c11.total_cor])
#y_error=np.array([i.s for i in c11.total_cor])
#
#ax.plot(t,y,color='w')
#ax.fill_between(t, y-y_error, y+y_error ,alpha=0.4,color='w')


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






