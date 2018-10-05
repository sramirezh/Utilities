#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 15:07:44 2018

@author: sr802
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

Rc=2.5

def LJ(r, epsilon, sigma, x, y):
    """
    computes the lennard-jones potential for different exponents.
    see the development in my notebook.
    Args:
        r       : the position at which the potential is computed.
        x       : the exponent of the repulsive part of the potential.
        y       : the exponent of the atractive part of the potential.
        epsilon : Energy of the minimum.
        sigma   : Radiuos at which the potential goes to zero.
    Returns:
        V Lennard Jones interaction energy
    """
    A=((x/y)**(x/(x-y))/((x/y)-1))
    
     
    V=A*epsilon*((sigma/r)**x-(sigma/r)**y) #-4*Epsilon*((Sigma/Rc)**12-(Sigma/Rc)**6)
    
    return V

xmin=1
r_vec=np.linspace(xmin,Rc)

LJ_12_6=LJ(r_vec,1,1,12,6)-LJ(Rc,1,1,12,6)
LJ_12_6_l=LJ(r_vec,10,1,12,6)-LJ(Rc,10,1,12,6)
LJ_12_10=LJ(r_vec,1,1,12,10)-LJ(Rc,1,1,12,10)

fig,ax=plt.subplots()
ax.plot(r_vec, LJ_12_6, label='$U^{LJ}_{12-6}$')
ax.plot(r_vec, LJ_12_6_l, label='$U^{LJ}_{12-6}Long$')
ax.plot(r_vec, LJ_12_10, label='$U^{LJ}_{12-10}$')
ax.axhline(y=0, xmin = xmin, xmax=Rc, ls=':',c='black')

plt.legend()