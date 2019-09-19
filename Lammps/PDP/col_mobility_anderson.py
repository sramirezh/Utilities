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




def velocity_colloid(a,T,mu,box_size,alpha):

    data_limit=(box_size/2) #Where the sampling volumes are outside the box
    beta=1/T
    
    data_file="prof_u.dat"
    data=pd.read_csv(data_file,sep=" ",header=None,skiprows=4).dropna(axis=1,how='all').values
    
    indexes=np.where(data[:,1]<data_limit)[0] 
    cs_bulk=data[indexes[-1],3]
    
    
    
    c_excess=data[indexes,3]-cs_bulk
    y=(data[indexes,1]-a)
    
    gamma=cf.integrate(y,c_excess,0,data_limit)
    
    integrand_k=c_excess/cs_bulk
    
    
    """There is a mistake as H << K"""
    K=cf.integrate(y,integrand_k,0,y[-1])
    
    integrand_1=integrand_k*y
    
    
    U_0=alpha/(beta*mu)*cf.integrate(y,integrand_1,0,y[-1])
    
    integrand_2=0.5*integrand_k*y**2
    
    H=cf.integrate(y,integrand_2,0,data_limit)/cf.integrate(y,integrand_1,0,y[-1])
    
    U=U_0*(1-(H+K)/a)
    
    plots(indexes,data,data_limit,y,c_excess)
    
    return U_0,K,H,U


def plots(indexes,data,data_limit,y,c_excess):
    """Plots"""
    plt.close('all')
    cf.set_plot_appearance()
    
    """Plot Solute concentration"""
    fig,ax=plt.subplots()
    print(data_limit)
    ax.set_xlim(0,data_limit)
    ax.plot(data[indexes,1],data[indexes,3])
    ax.set_xlabel(r'$r[\sigma] $')
    ax.set_ylabel(r'$c_s^B [\sigma^{-3}]$')
    ymin,ymax=plt.ylim()
    fig.tight_layout()
    ax.axvline(x=a, ymin=0, ymax=1,ls='--',c='black')
    plt.savefig("Solute_concentration.pdf")
    
    
    """Plot Excess_solute"""
    fig,ax=plt.subplots()
    ax.plot(y,c_excess,label="Excess solute")
    ax.set_ylabel(r'$e^{-\beta\phi(y)}-1$')
    ax.set_xlabel(r'$y $')
    xmin,xmax=plt.xlim()
    ax.set_xlim(0,xmax)
    ax.axhline(y=0, xmin=xmin, xmax=xmax,ls='--',c='black')
    fig.tight_layout()
    plt.savefig("Solute_excess.pdf")
    
    plt.show()



print("Remember to define the viscosity, box_size,Rh...")



"""BE CAREFUL WITH THIS"""
#

"""Need to recompute"""


"""Colloid data"""
alpha=0.002
box_size=16
a=4.225 #Hydrodynamic radius
mu=1.2
T=0.845

U_0,K,H,U=velocity_colloid(a,T,mu, box_size,alpha)
print(U_0,K,H,U)

