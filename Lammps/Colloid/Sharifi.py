#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 12:18:38 2019

Simple script to get Sharifi concentrations:
    
It is an ideal solution in the bulk

@author: sr802
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path

import Lammps.core_functions as cf

grad_c = np.array([0.0062, 0.0077, 0.0093, 0.0108])

n0 = 50
lx = 53.86
ly = 43.09
lz = 43.09

nbins = 8

delta_x = lx/nbins

vol_bin = delta_x*ly*lz

delta_N = grad_c * vol_bin*delta_x

fig,ax = plt.subplots()
cf.set_plot_appearance()

xpos = [i*delta_x for i in range(nbins+1)]


fig,ax = plt.subplots()

for grad in grad_c:
    delta_N = grad * vol_bin*delta_x
    ns = np.array([n0+i*delta_N for i in range(nbins+1)])
    ax.step(xpos,ns,label = r'$\nabla c_s = %s $'%grad, where = 'post')
    
    #Exact distribution
    ns_exact = [n0+x*grad*vol_bin for x in xpos]
    color=ax.lines[-1].get_color()
    
    ax.plot(xpos,ns_exact, color = color,linestyle='dashed' )

ax.set_xlim(0, lx)   
ymin,ymax=plt.ylim()
ax.set_ylim(0, 1.2*ymax)  
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$N_s(x)$")
ax.legend()
plt.tight_layout()
plt.savefig("Grad_c.pdf")

fig,ax = plt.subplots()

for grad in grad_c:
    delta_N = grad * vol_bin*delta_x
    cs = np.array([n0+i*delta_N for i in range(nbins+1)])/vol_bin
    ax.step(xpos,cs,label = r'$\nabla c_s = %s $'%grad, where = 'post')
    
    #Exact distribution
    cs_exact = np.array([n0+x*grad*vol_bin for x in xpos])/vol_bin
    color=ax.lines[-1].get_color()
    
    ax.plot(xpos,cs_exact, color = color,linestyle='dashed' )

ax.set_xlim(0, lx)   
ymin,ymax=plt.ylim()
ax.set_ylim(0, 1.2*ymax)  
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$c_s(x)$")
ax.legend()
plt.tight_layout()
plt.savefig("conc.pdf")


cs = ns/vol_bin



    

