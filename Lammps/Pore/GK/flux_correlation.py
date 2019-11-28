#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 17:22:14 2019
Contains the classes flux and correlation that help to compute correlations
@author: sr802
"""
import multiprocessing
import numpy as np
from tqdm import tqdm  
from joblib import Parallel, delayed
import pickle as pickle
import Lammps.core_functions as cf
import Others.Statistics.FastAverager as stat

try:
    from uncertainties import unumpy,ufloat
except ImportError as err2:
    print(err2)
    
# =============================================================================
# Class definition
# =============================================================================
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
    """
    TODO, here that the error is computed as a simple error 
    """
    def __init__(self,flux1,flux2, max_delta):
        self.flux1 = flux1
        self.flux2 = flux2
        self.max_delta = max_delta
        self.initial_check()
        dimension = self.dimension
        self.cor = (dimension+1)*[0]  #list to store the correlations, the last is the total
        self.norm = (dimension+1)*[0]  #list to store the normalisation correlation with t=0 (var1(0) var2(0))
        self.cor_norm = (dimension+1)*[0]  #list to store the normalised correlations, the last is the total
        self.dic_label={0:'x',1:'y',2:'z',self.dimension:'Total'}
        
    def initial_check(self):
        """
        Checks if the time series are equal and the fluxes have the same number of components
        """
        self.dimension=self.flux1.dimension
        if self.flux1.dimension != self.flux2.dimension:
            print("The fluxes do not have the same dimension")
        else:
            self.dimension=self.flux1.dimension
            
        if self.flux1.times.all != self.flux2.times.all:
            print("The fluxes were not measured for the same times")
        else:
            self.times = self.flux1.times[:self.max_delta]
            

    
    def correlate_one_d(self,dim):
        """
        Performs a correlation between 1d components by evaluating the products at at different delta t.
        Args:
            dim is the component to be evaluated
        """
        num_cores = multiprocessing.cpu_count()
        var1 = self.flux1.components[:,dim]
        var2 = self.flux2.components[:,dim]
        max_delta = self.max_delta
        cor = Parallel(n_jobs=num_cores)(delayed(compute_correlation_dt)(var1,var2,i) for i in tqdm(range(max_delta)))
        norm = cor[0].nominal_value
        self.norm[dim] = norm
        self.cor[dim] = np.array(cor)
        self.cor_norm[dim]=np.array(cor)/norm
        
    def evaluate(self):
        """
        Performs the correlations of the 1d components,calluing correlate_one_d, and adds them up to the total.
        
        """
        total = np.zeros(self.max_delta)
        for dim in range(self.dimension):
            self.correlate_one_d(dim)
            total = total + self.cor[dim]
        total = total/3
        self.cor[-1] = total
        self.norm[-1] = total[0].nominal_value
        self.cor_norm[-1]=total/total[0].nominal_value
        
        
    def plot_individual(self,fig,ax,dim=0,alpha=0.4,every=1,norm=True):
        """
        Args:
            ax axes object
            fig Figure 
            dim is the dimension, for example:in a 3D vector, 0-x, 1-y, 2-z and 3-total.
            alpha is the transparency of the filling
            every to not have so many points
            norm True if normalised
            The axis label is given here but it could be renamed later
        """
        if norm==True:
            cor=self.cor_norm[dim]
        else:
            cor=self.cor[dim]
            
        
        y=np.array([i.n for i in cor])
        y_error=np.array([i.s for i in cor])
        ax.plot(self.times[::every],y[::every],label=self.dic_label[dim])
        ax.fill_between(self.times, y-y_error, y+y_error ,alpha=0.4)
        
        
        return fig,ax
    
    def plot_all(self,fig,ax,alpha=0.4):
        for dim in range(self.dimension+1):
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
        


class bundle_correlation(correlation):
    def __init__(self,cor,times,flux1_name,flux2_name):
        self.dimension=1
        self.cor = [cor]
        self.times = times
        self.norm = cor[0].nominal_value
        self.cor_norm = [cor/self.norm]
        self.flux1_name = flux1_name
        self.flux2_name = flux2_name
        
        
    def plot(self,fig,ax,dim=0,alpha=0.4,every=1,ax_label = True,norm = True):
        """
        Args:
            ax axes object
            fig Figure 
            dim is the dimension, for example:in a 3D vector, 0-x, 1-y, 2-z and 3-total.
            alpha is the transparency of the filling
            every to not have so many points
            norm True if normalised
            The axis label is given here but it could be renamed later
        """
        if norm == True:
            cor = self.cor_norm[dim]
        else:
            cor = self.cor[dim]
            
        
        y=np.array([i.n for i in cor])
        print(y)
        y_error=np.array([i.s for i in cor])
        ax.plot(self.times[::every],y[::every])
        ax.fill_between(self.times, y-y_error, y+y_error ,alpha=0.4)
        
        # It mostly means that the plot will not be further used, for example not used to compare correlations
        if ax_label==True:
            ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
            ax.set_ylabel(r'$\langle %s(t)%s(0) \rangle$'%(self.flux1_name,self.flux2_name))
            ax.set_xlabel('time')

        return fig,ax
    
    def transport_coeff(self, T, pref, xmin, xmax):
        """
        T is the temperature
        pref is the prefactor, eg.system volume etc
        xmin
        """
    
        y = np.array([i.n for i in self.cor[0]])
        x = self.times
        I = cf.integrate(x,y,xmin,xmax)
        self.coeff= (pref/T)*I
        
        return self.coeff

        
def compute_correlation_dt(var1,var2,delta):
    """
    This is a VERY GENERAL function
    
    *****
    It is important to keep it outside the class as it is going to be evaluated in parallel and if it is a method of the class, it would require
    to load the instance on each processor which could be very expensive
    *****
    
    Computes the correlation for two variables for a given delta t
    Args:
        var1: time series for the first variable 
        var2: time series for the first variable
        Note that var1 and var2 have to have the same the same time series
        delta: every this number of steps, we take the variable 2 to compute the correlation
    """
    cf.blockPrint()
    if delta != 0:
        var1 = var1[::delta]
        var2 = var2[::delta]
        cor = (var1*np.roll(var2,-1,axis=0))
        cor = cor[:-1]#The last contribution is the last-the initial msd
    else:
        cor = var1*var2 
    
    average = stat.fast_averager(cor)[0]
    
    cf.enablePrint()
    
    return ufloat(average[1],average[2]) 