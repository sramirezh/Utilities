#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:53:27 2019
This is to plot the potential in table.dat that is in Sharifi's simulation and comparte with my prediction
@author: sr802
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path


from Lammps.PDP.Plots.LJ_plotter import LJ
import Lammps.core_functions as cf



def sharifi_LJ(r, epsilon, sigma, lam):
    """

    """
    V=4*epsilon*((sigma/r)**12-lam*(sigma/r)**6) #-4*Epsilon*((Sigma/Rc)**12-(Sigma/Rc)**6)

    return V

data = np.loadtxt("table.dat")

xmin = 1
Rc = 2.5
lam = 3 #CIJ from sharifi

# Getting the sharifi's datawith sigma and epsilon modified




# Getting the real sharifi expression

r_vec=np.linspace(xmin,Rc)
r_vec_s=np.linspace(0.01,Rc,300)
LJ_sharifi=sharifi_LJ(r_vec_s,1,1,3)-sharifi_LJ(Rc,1,1,3)

LJ_12_6=LJ(r_vec,1,1,12,6)-LJ(Rc,1,1,12,6)


plt.close('all')
fig,ax=plt.subplots()
ax.plot(r_vec, LJ_12_6, label='$U^{LJ}_{12-6}$')
ax.plot(r_vec_s[100:], LJ_sharifi[100:],'o',ms =2, label='LJ sharifi')
ax.plot(data[220:,1],data[220:,2], ls='--', label = 'Lammps')
ax.axhline(y=0, xmin = 0, xmax=1, ls=':',c='black')
plt.legend(loc = "bottom right")


plt.savefig("sharifi_potential.pdf")