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

a=4 #Radius of gyration
T=1 #They do not mention it explicitely
mu=2.3 #Need to recompute


data_limit=6 #Where the sampling volumes are outside the box
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
K=cf.integrate(y,integrand_k,0,data_limit)

integrand_1=integrand_k*y


U_0=alpha/(beta*mu)*cf.integrate(y,integrand_1,0,data_limit)

integrand_2=0.5*integrand_k*y**2

H=cf.integrate(y,integrand_2,0,data_limit)/cf.integrate(y,integrand_1,0,data_limit)

U=U_0*(1-(H+K)/a)


plt.close('all')
fig,ax=plt.subplots()

ax.set_xlim(0,10)
ax.plot(data[:,1],data[:,3],label="Solute concentration")

plt.legend()
plt.show()


fig,ax=plt.subplots()

ax.set_xlim(0,10)
ax.plot(y,c_excess,label="Excess_solute")
ax.set_xlabel(r'$y $')
plt.legend()
plt.show()


print U_0,U

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

#integrand=((2*a*T)/(9*mu)*((0.5/r_a)-(3/2*r_a)+(r_a)**2))*A1*d_phi_r

"""
Check that i am integrating over r and the integrand is a function of r_a
"""
#
#v_x=cf.integrate(r,integrand,rmin,rmax) 
#
#from scipy.interpolate import splev,splrep,splint,splder
#
#tck=splrep(r,phi)
#spline_x=np.linspace(rmin,rmax,2000)
#spline_y=splev(spline_x,tck)

#
#
#v_x_spline=cf.integrate(spline_x,spline_y,rmax,rmin) 
#print(v_x,v_x_spline)


#fig,ax=plt.subplots()
#ax.set_ylabel(r'$\phi(r) $')
#ax.set_xlabel(r'$r $')
#ax.plot(r,phi,'o',label="Points")
#ax.plot(spline_x,spline_y,label="Spline")
#plt.legend()
#
#
#tckder=splder(tck)
#
#der=splev(spline_x,tck, der=1)
#
#
#fig,ax=plt.subplots()
#ax.set_ylabel(r'$\frac{d\phi(r)}{dr} $')
#ax.set_xlabel(r'$r $')
#ax.plot(r,dome_ra,'o',label="Points")
#ax.plot(spline_x,der,label="Spline")
#plt.legend()
#plt.tight_layout()
