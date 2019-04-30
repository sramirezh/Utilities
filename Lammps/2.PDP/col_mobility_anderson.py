#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:14:47 2019
Estimates the velocity only using the concentration profile.

To do:
    
    
Define Viscosity
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




def mobility_colloid(a,T,mu,box_size,alpha):

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


def mobility_polymer(a,T,mu,box_size,grad_mu,cut_off):
    data_limit=(box_size/2)*cut_off #Where the sampling volumes are outside the box
    beta=1/T
    
    data_file="prof_u.dat"
    data=pd.read_csv(data_file,sep=" ",header=None,skiprows=4).dropna(axis=1,how='all').values
    
    indexes=np.where(data[:,1]<data_limit)[0] 
    cs_bulk=data[indexes[-1],3]
    alpha=0.1  #Should be dependent on the 
    
    
    
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
    
    U_1=-U_0*(H+K)/a
    U=U_0+U_1
    
#    plots(indexes,data,data_limit,y,c_excess)
    return U_0,K,H,U_1,U


def plots(indexes,data,data_limit,y,c_excess):
    """Plots"""
    plt.close('all')
    cf.set_plot_appearance()
    
    """Plot Solute concentration"""
    fig,ax=plt.subplots()
    
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



print "Remember to define the viscosity, box_size,Rh..."



"""Polymer"""
grad_mu=0.1
box_size=19.4754
a=3.652945 #Hydrodynamic radius
T=1
mu=1.55639001606637 
cut_off=0.9
mobility_polymer(a,T,mu, box_size,grad_mu,cut_off)

rh_vect=np.linspace(2.5,7)

results=[]
for a in rh_vect:
    
    results.append(mobility_polymer(a,T,mu, box_size,grad_mu,cut_off))

results=np.array(results)


"""rh_analysis"""
plt.close('all')
fig,ax=plt.subplots()
ax.plot(rh_vect,results[:,0],label=r"$U_0$")
ax.plot(rh_vect,results[:,3],label=r"$U_1$")
ax.plot(rh_vect,results[:,4],label=r"$U=U_0+U_1$")
ax.set_xlabel(r'$R_h $')
xmin,xmax=plt.xlim()
ax.axhline(y=0, xmin=xmin, xmax=xmax,ls='--',c='black')
ax.set_xlim(xmin,xmax)

fig.tight_layout()
plt.legend(loc='upper rigth')
plt.savefig("Rh_vs_U.pdf")



fig,ax=plt.subplots()
ax.plot(rh_vect,results[:,1],label=r"$K$")
ax.plot(rh_vect,results[:,2],label=r"$H$")
ax.set_xlabel(r'$R_h $')
xmin,xmax=plt.xlim()
ax.axhline(y=0, xmin=xmin, xmax=xmax,ls='--',c='black')
ax.set_xlim(xmin,xmax)

fig.tight_layout()
plt.legend(loc='upper rigth')
plt.savefig("Rh_vs_KH.pdf")

plt.show()




"""BE CAREFUL WITH THIS"""
#

"""Need to recompute"""


"""Colloid data"""
alpha=0.002
box_size=14
a=2.225 #Hydrodynamic radius
mu=1.2
T=0.845

#U_0,K,H,U=mobility_colloid(a,T,mu, box_size,alpha)
#print U_0,K,H,U







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