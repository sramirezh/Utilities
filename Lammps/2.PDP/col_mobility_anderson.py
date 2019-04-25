#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:14:47 2019
Estimates the velocity only using the concentration profile.
@author: sr802
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


import os
import sys
import warnings
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
warnings.filterwarnings("ignore")
import Lammps.core_functions as cf
import pandas as pd

grad_mu=0.1

a=3.652945 #Hydrodynamic radius
T=1 #They do not mention it explicitely

"""Need to recompute"""
mu=2.3 


data_limit=9 #Where the sampling volumes are outside the box
beta=1

data_file="prof_u.dat"
data=pd.read_csv(data_file,sep=" ",header=None,skiprows=4).dropna(axis=1,how='all').values

indexes=np.where(data[:,1]<data_limit)[0] 
cs_bulk=data[indexes[-1],3]
alpha=grad_mu*cs_bulk  #Should be dependent on the 

c_excess=data[indexes,3]-cs_bulk
y=(data[indexes,1]-a)

gamma=cf.integrate(y,c_excess,0,data_limit)

integrand_k=c_excess/cs_bulk


"""There is a mistake as H << K"""
K=cf.integrate(y,integrand_k,0,y[-1])

integrand_1=integrand_k*y


U_0=alpha/(beta*mu)*cf.integrate(y,integrand_1,0,data_limit)

integrand_2=0.5*integrand_k*y**2

H=cf.integrate(y,integrand_2,0,data_limit)/cf.integrate(y,integrand_1,0,data_limit)

U=U_0*(1-(H+K)/a)

print U_0,U






"""Plots"""
plt.close('all')
cf.set_plot_appearance()

"""Plot Solute concentration"""
fig,ax=plt.subplots()

ax.set_xlim(0,data_limit)
ax.plot(data[:,1],data[:,3])
ax.set_xlabel(r'$r[\sigma] $')
ax.set_ylabel(r'$c_s^B [\sigma^{-3}]$')
ymin,ymax=plt.ylim()
ax.axvline(x=a, ymin=0, ymax=1,ls='--',c='black')
fig.tight_layout()
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



def r_g(N):
    """
    Estimation for the gyration radius dependence on the number of monomers
    from Dunweg et al"Corrections to scaling in the hydrodynamic properties of dilute polymer solutions"
    Args:
        N the number of monomers
    """
    return (0.2706*N**(1.1754)-0.32*N**(0.62))**0.5


def r_h(N):
    """
    Estimation for the Hydrodynamic radius dependence on the number of monomers
    from Dunweg et al"Corrections to scaling in the hydrodynamic properties of dilute polymer solutions"
    Args:
        N the number of monomers
    """
    return 1/(3.131*N**(-0.5877)-3.04/N)