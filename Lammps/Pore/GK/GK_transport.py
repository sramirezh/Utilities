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
import Lammps.Pore.qsub.simulation_results as sr

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

cwd = os.getcwd() #current working directory


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
        ax.set_ylabel(r'$\langle %s(t)%s(0) \rangle$'%(self.flux1.name,self.flux2.name))
        ax.set_xlabel('time')
        return fig,ax
    
    def save(self,file_name):
        """
        Saves the instance 
        """
        afile = open(r'%s.pkl'%file_name, 'wb')
        pickle.dump(self, afile)
        afile.close()    
        

root_pattern="run*"
directory_pattern='[0-9]*'
parameter_id='number'

dictionary={'vx_Solv':r'$v^x_{f}$','vx_Solu':r'$v^x_{s}$','vx_Sol':r'$v^x_{sol}$'}


bundles_mu=sr.initialise_sim_bundles(root_pattern,parameter_id,directory_pattern,dictionary)
final_mu=sr.simulation_bundle(bundles_mu,parameter_id,3,cwd,dictionary=dictionary)






"""
Computes the diffusion coefficient
"""
#Input parameters
input_file='fluxes.dat'
delta_t=0.005

Data=cf.read_data_file(input_file)
data1=Data.values
times=(data1[:,0]-data1[0,0])*delta_t

total_flux=flux(data1[:,[1,3,5]],times,"Q")
solute_excess=flux(data1[:,[2,4,6]],times,"J_s-c_s^BQ")

max_delta=int(len(times)*0.004) #Maximum delta of time to measure the correlation


# =============================================================================
# Defining, evaluating and saving the correlations
# =============================================================================

def run():
    """
    Very specific function
    """
    c11=correlation(total_flux,total_flux,max_delta)
    c11.evaluate()
    c11.save('c11')
    
    c12=correlation(total_flux,solute_excess,max_delta)
    c12.evaluate()
    c12.save('c12')
    
    c21=correlation(solute_excess,total_flux,max_delta)
    c21.evaluate()
    c21.save('c21')
    
    c22=correlation(solute_excess,solute_excess,max_delta)
    c22.evaluate()
    c22.save('c22')


# =============================================================================
# Loading all
# =============================================================================
c11 = cf.load_instance("c11.pkl")
c22 = cf.load_instance("c22.pkl")
c12 = cf.load_instance("c12.pkl")
c21 = cf.load_instance("c21.pkl")




plt.close('all')
cf.set_plot_appearance()

# =============================================================================
# Ploting all the correlations
# =============================================================================

# For c11
fig,ax=plt.subplots()

fig,ax = c11.plot_all(fig,ax)
ax.set_xscale('log')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("correlation11.pdf")

# For c12
fig,ax=plt.subplots()

fig,ax = c12.plot_all(fig,ax)
ax.set_xscale('log')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("correlation12.pdf")

# For c21
fig,ax=plt.subplots()

fig,ax = c21.plot_all(fig,ax)
ax.set_xscale('log')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("correlation21.pdf")

# For c22
fig,ax=plt.subplots()

fig,ax = c22.plot_all(fig,ax)
ax.set_xscale('log')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("correlation22.pdf")




# =============================================================================
# Plot the cross coefficients
# =============================================================================


fig,ax = plt.subplots()

fig,ax = c12.plot_individual(fig,ax,3)
ax.lines[-1].set_label(r'$\langle (%s(t))(%s(0)) \rangle$'%(c12.flux1.name,c12.flux2.name)) #modifying the label of the last created line
fig,ax = c21.plot_individual(fig,ax,3)
ax.lines[-1].set_label(r'$\langle (%s(t))(%s(0)) \rangle$'%(c21.flux1.name,c21.flux2.name)) #modifying the label of the last created line
ax.set_xscale('log')

plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("crossed.pdf")




 Create something similar to the bundle but that just stores the correlations because the other things are very heavy. so the object is just a list of the total correlations
 and you can just average