#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:14:47 2019
Estimates the velocity only using the concentration profile.
    
Define Viscosity
@author: sr802
"""

import numpy as np
import matplotlib.pyplot as plt


import os
import sys
import warnings
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
warnings.filterwarnings("ignore")
import Lammps.core_functions as cf
import pandas as pd

def plateau_finder(data,tol=0.0003):
    """
    Function that finds the plateaus of a distribution y, that has x spacing constant
    Args:
        data 1D data that has delta in the other dimension constant
        tol tolerance for the variance
    """
    from scipy.ndimage.filters import generic_filter
    filt_data = generic_filter(data, np.std, size=3)
    plat_index=np.where(filt_data<(np.min(filt_data)+tol))[0]
    
    plateaus=group_consecutive(plat_index)
    
    return plateaus


def group_consecutive(data):
    """
    Groups consecutive indexes in an array an return a list with all the groups
    Args:
        data 1D data array
    Returns:
        results a list with all the consecutive indexes grouped
    """
    from itertools import groupby
    from operator import itemgetter
    results=[]
    for k, g in groupby(enumerate(data), lambda i_x: i_x[0]-i_x[1]):
        results.append(list(map(itemgetter(1), g)))
        
         
    return results

def plot_plateau(a,data,plat_indexes,indexes_box,rh_origin):
    if rh_origin=='K':
        name='plateau_k.pdf'
    if rh_origin=='md':
        name='plateau_md.pdf'
        
    cs_bulk = data[plat_indexes[-1],3]
    fig,ax = plt.subplots()
    ax.plot(data[indexes_box,1],data[indexes_box,3])
    ax.plot(data[plat_indexes,1],data[plat_indexes,3],'--', label = "plateau" )
    ax.set_xlabel(r'$r[\sigma] $')
    ax.set_ylabel(r'$c_s^B [\sigma^{-3}]$')
    ymin,ymax=plt.ylim()
    ax.axhline(y=cs_bulk, xmin=0, xmax=1,ls='--',c='black')
    fig.tight_layout()
    #ax.set_xlim(6,box_size/2)
    #ax.set_ylim(0.37,0.38)
    ax.axvline(x=a, ymin=0, ymax=1,ls='--',c='black')
    ax.legend()
    plt.savefig(name)
    plt.show()
    


def velocity_colloid(a,T,eta,box_size,grad_mu,rh_origin='K',plot=True):
    data_limit=box_size/2 #Where the sampling volumes are outside the box
    
    #Uncomment to deal with Lammps data

    #data_file="prof_u.dat"
    
    #data=pd.read_csv(data_file,sep=" ",header=None,skiprows=4).dropna(axis=1,how='all').values
    
    # created with density distribution_from atom
    data =cf.read_data_file("prof_u.dat").values
    #Indexes inside the box, to avoid spherical shells outside the box
    indexes_box=np.where(data[:,1]<data_limit)[0] 
    

    plat_indexes=max(plateau_finder(data[indexes_box,3],tol = 0.003))

    indexes=np.arange(0,plat_indexes[-1])
    
    cs_bulk=data[indexes[-1],3]
    
    print ("The concentration in the bulk is %s"%cs_bulk)
    alpha=beta*grad_mu*cs_bulk
    
    
    
    c_excess=data[indexes,3]-cs_bulk
    y=(data[indexes,1]-a)
    
#    gamma=cf.integrate(y,c_excess,0,data_limit)
    
    integrand_k=c_excess/cs_bulk
    
    
    """There is a mistake as H << K"""
    K=cf.integrate(y,integrand_k,0,y[-1])
    
    integrand_1=integrand_k*y
    
    L=cf.integrate(y,integrand_1,0,y[-1])
    U_0=alpha/(beta*eta)*L
    
    integrand_2=0.5*integrand_k*y**2
    
    H=cf.integrate(y,integrand_2,0,data_limit)/L
    
    U_1=-U_0*(H+K)/a
    U=U_0+U_1
    if plot == True:
        plots(a,indexes,data,data_limit,y,c_excess,rh_origin)
        plot_plateau(a,data,plat_indexes,indexes_box,rh_origin)
        
    return K,L,H,U_0,U_1,U


def plots(a,indexes,data,data_limit,y,c_excess,rh_origin):
    
    if rh_origin=='K':
        name1='Solute_concentration_k.pdf'
        name2='Solute_excess_k.pdf'
    if rh_origin=='md':
        name1='Solute_concentration_md.pdf'
        name2='Solute_excess_md.pdf'
    
    """Plots"""
    plt.close('all')
    cf.set_plot_appearance()
    cs_bulk=data[indexes[-1],3]
    
    """Plot Solute concentration"""
    fig,ax=plt.subplots()
    ax.set_xlim(0,data_limit)
    ax.plot(data[indexes,1],data[indexes,3])
    ax.set_xlabel(r'$r[\sigma] $')
    ax.set_ylabel(r'$c_s^B [\sigma^{-3}]$')
    ymin,ymax=plt.ylim()
    ax.axhline(y=cs_bulk, xmin=0, xmax=1,ls='--',c='black')
    fig.tight_layout()
    ax.axvline(x=a, ymin=0, ymax=1,ls='--',c='black')
    plt.savefig(name1)
    
    
    """Plot Excess_solute"""
    fig,ax=plt.subplots()
    ax.plot(y,c_excess,label="Excess solute")
    ax.set_ylabel(r'$e^{-\beta\phi(y)}-1$')
    ax.set_xlabel(r'$y $')
    xmin,xmax=plt.xlim()
    ax.set_xlim(0,xmax)
    ax.axhline(y=0, xmin=xmin, xmax=xmax,ls='--',c='black')
    fig.tight_layout()
    plt.savefig(name2)
    
    plt.show()


print("Remember to define the viscosity, box_size,Rh...")



"""BE CAREFUL WITH THIS"""
#

"""Need to recompute"""


"""Colloid data"""
grad_mu = 0.6
T = 1
beta = 1/T
box_size = 20
a = 7.542431653 #Hydrodynamic radius
eta = 2.249983808
T = 1

results=velocity_colloid(a,T,eta, box_size,grad_mu,'md')
print('K,L,H,U_0,U_1,U')
print(results)

