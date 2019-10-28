#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:14:47 2019
Estimates the velocity only using the concentration profile.


The file prof_u.dat contains # Chunk Coord1 Ncount density/number

the file Parameters, has N R_H^K R_H^{md}
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
         



def plateau_finder(data,tol=0.0003):
    """
    Function that finds the plateaus of a distribution y, that has x spacing constant
    Args:
        data 1D data that has delta in the other dimension constant
        tol tolerance for the variance
    """
    from scipy.ndimage.filters import generic_filter
    tol=0.0003
    filt_data = generic_filter(data, np.std, size=3)
    plat_index=np.where(filt_data<(np.min(filt_data)+tol))[0]
    
    plateaus=group_consecutive(plat_index)
    
    return plateaus

def velocity_polymer(a,T,eta,box_size,grad_mu,rh_origin='K',plot=True):
    data_limit=box_size/2 #Where the sampling volumes are outside the box
    beta=1/T
    
    data_file="prof_u.dat"
    data=pd.read_csv(data_file,sep=" ",header=None,skiprows=4).dropna(axis=1,how='all').values
    
    #Indexes inside the box, to avoid spherical shells outside the box
    indexes_box=np.where(data[:,1]<data_limit)[0] 
    
    plat_indexes=max(plateau_finder(data[indexes_box,3]))

    indexes=np.arange(0,plat_indexes[-1])
    
    cs_bulk=data[indexes[-1],3]
    alpha=beta*grad_mu*cs_bulk
    
    
    
    c_excess=data[indexes,3]-cs_bulk
    y=(data[indexes,1]-a)
    
    gamma=cf.integrate(y,c_excess,0,data_limit)
    
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
    if plot==True:
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



def rh_analysis(rmin,rmax,T,eta, box_size, grad_mu):
    """
    Analyising the Change of the parameters while varyin the hydrodynamic radius between rmin and rmax
    """
    rh_vect=np.linspace(rmin,rmax)
    
    results=[]
    for a in rh_vect:
        
        results.append(velocity_polymer(a,T,eta, box_size,grad_mu,'K',plot=False))
    
    results=np.array(results)
    
    
    """rh_analysis"""
    plt.close('all')
    fig,ax=plt.subplots()
    ax.plot(rh_vect,results[:,0],label=r"$U_0$")
    ax.plot(rh_vect,results[:,3],label=r"$U_1$")
    ax.plot(rh_vect,results[:,4],label=r"$U=U_0+U_1$")
    ax.set_xlabel(r'$R_h $')
    ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
    ax.axvline(x=a_k, ymin=0, ymax=1,ls='--',c='black',label=r'$R_h^{K}$')
    ax.axvline(x=a_md, ymin=0, ymax=1,ls='--',c='blue',label=r'$R_h^{\lambda}$')
    fig.tight_layout()
    plt.legend(loc='upper rigth')
    plt.savefig("Rh_vs_U.pdf")
    
    
    
    fig,ax=plt.subplots()
    ax.plot(rh_vect,results[:,1],label=r"$K$")
    ax.plot(rh_vect,results[:,2],label=r"$H$")
    ax.set_xlabel(r'$R_h $')
    xmin,xmax=plt.xlim()
    ax.axvline(x=a_k, ymin=0, ymax=1,ls='--',c='black',label=r'$R_h^{K}$')
    ax.axvline(x=a_md, ymin=0, ymax=1,ls='--',c='blue',label=r'$R_h^{\lambda}$')
    ax.set_xlim(xmin,xmax)
    
    fig.tight_layout()
    plt.legend(loc='upper rigth')
    plt.savefig("Rh_vs_KH.pdf")
    
    plt.show()



def plot_plateau(a,data,plat_indexes,indexes_box,rh_origin):
    if rh_origin=='K':
        name='plateau_k.pdf'
    if rh_origin=='md':
        name='plateau_md.pdf'
        
    cs_bulk=data[plat_indexes[-1],3]
    fig,ax=plt.subplots()
    ax.plot(data[indexes_box,1],data[indexes_box,3])
    ax.plot(data[plat_indexes,1],data[plat_indexes,3],'--')
    ax.set_xlabel(r'$r[\sigma] $')
    ax.set_ylabel(r'$c_s^B [\sigma^{-3}]$')
    ymin,ymax=plt.ylim()
    ax.axhline(y=cs_bulk, xmin=0, xmax=1,ls='--',c='black')
    fig.tight_layout()
    #ax.set_xlim(6,box_size/2)
    #ax.set_ylim(0.37,0.38)
    ax.axvline(x=a, ymin=0, ymax=1,ls='--',c='black')

    plt.savefig(name)
    plt.show()
    
    
# =============================================================================
# MAIN
# =============================================================================

print("Remember to define the viscosity")

cf.set_plot_appearance()

"""Polymer"""
grad_mu=0.1
T=1
beta=1/T
eta=1.55639001606637 

#Automatically load the parameters
parameters=cf.read_data_file('../data.dat').values
N=int(cf.extract_digits(os.getcwd().split('/')[-1])[0]) #Number of monomers obtained from the path

print("Running for N=%d" %N)


index_1=np.where(parameters==N)[0][0] #gets the line with the parameters of this Number of monomers
box_size=parameters[index_1,1]

print('\nUsing the Rh estimation from Kirkwood')

a_k=parameters[index_1,2]

results=velocity_polymer(a_k,T,eta, box_size,grad_mu,rh_origin='K')
mobility=results[-1]/grad_mu
print('K,L,H,U_0,U_1,U,mobility')
print(results,mobility)

print('\nUsing the Rh estimation from the mobility')

D=parameters[index_1,3]
a_md=T/(6*np.pi*eta*D)

results=velocity_polymer(a_md,T,eta, box_size,grad_mu,'md')
mobility=results[-1]/grad_mu
print('K,L,H,U_0,U_1,U,mobility')
print(results,mobility)



"""
(UNCOMMENT) Rh dependence analyisis
"""
#rh_analysis(2.5,box_size/2-1,T,eta, box_size, grad_mu)






