#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 13:24:39 2019
This scripts plots the results inside 0.GCMC_single/Fig_5.8 
@author: sr802
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import warnings
warnings.filterwarnings("ignore")

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Others.JohnsonEOS.LJEOS as ljeos


simulation_data=cf.read_data_file("fig.dat")

data=simulation_data.values
cf.set_plot_appearance()


T=2.0 #Temperature


rho_vect=np.linspace(0.01,1)

p_vect=[]
mu_ex_vect=[]
for i in range(np.size(rho_vect)):
    p_vect.append(ljeos.pressure(rho_vect[i],T))
    mu_ex_vect.append(ljeos.mu_ex(rho_vect[i],T))



plt.close("all")
names=[r'$\mu^{ex} $',r'$P $']
fig1,ax1=plt.subplots()

ax1.plot(rho_vect,mu_ex_vect,label="EOS")
ax1.plot(data[:,1],data[:,2],".",label="Simulation")

ax1.set_xlabel(r'$\rho $')
ax1.set_ylabel(r'$\mu^{ex} $')
ax1.set_xlim(0,1)
ax1.set_ylim(-5,10)
plt.legend(loc='upper_left')
fig1.tight_layout()
fig1.savefig("mu_ex.pdf")


fig2,ax2=plt.subplots()
ax2.plot(rho_vect,p_vect,label="EOS")
ax2.plot(data[:,1],data[:,3],".",label="Simulation")
ax2.set_xlabel(r'$\rho $')
ax2.set_ylabel(r'$P $')
ax2.set_xlim(0,1)
ax2.set_ylim(-5,10)
plt.legend(loc='upper_left')
fig2.tight_layout()
fig2.savefig("pressure.pdf")

plt.show()