#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:14:47 2019
Estimates the velocity using Sharifi-Mood
a is the particle radius
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


alpha=0.0077
a=3.23
sigma=1
c23=3
c13=1
n=136
T=1 #They do not mention it explicitely
mu=2.3
delta=0.2
rcut=4.68
rmin=3.23
rmax=9
r=np.linspace(rmin,rmax,200)


r_a=r/a

def concentration(n,a): 
    """
    Density of particles in the colloidal particle
    """
    return n/(4./3*np.pi*a**3)

c=concentration(n,a)

""""OVerrided"""
c=0.73

omega=(1/(r_a))*((2/(r_a+1)**3)+(2/(r_a-1)**3)+(1/(r_a+1)**2)-(1/(r_a-1)**2))
L3=(4*np.pi*sigma**6*c/(12*T))*(c23-c13) 
L=(L3)**(1/3)
phi=-(L/(sigma*a))**3*omega



#expression without diffusion in the outer layer
def A1(r,rcut):
    
    r_a=r/a
    indexes=np.where(r<rcut)[0]
    A1=alpha*a*(np.exp(2*(L/sigma)**3/((r-a)**3))-1+r_a)
    A1[indexes]=A1(rcut,0)
    return A1

def A2(r):
    r_a=r/a
    return alpha*a*((3/2)*(np.exp(2*(L/sigma)**3/((r-a)**3))-1)+r_a+0.5*(r_a)**(-2))


d_ome_ra=((-1/(r_a)**2)*(2/(r_a+1)**3)+(2/(r_a-1)**3)+(1/(r_a+1)**2)-(1/(r_a-1)**2)+(1/r_a)*((-6/(r_a+1)**4)+(-6/(r_a-1)**4)+(-2/(r_a+1)**3)-(-2/(r_a-1)**3)))
d_ome_r=(1/a)*d_ome_ra

d_phi_r=-((L/(sigma*a))**3)*d_ome_r

plt.close('all')
fig,ax=plt.subplots()


ax.plot(r,A1(r,rcut),label="CSWOD")
ax.plot(r,A2(r),label="CSWD")
ax.set_xlim(4,9)
ax.set_ylim(0,0.11)
plt.legend()



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
