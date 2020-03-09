#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:06:41 2020
Short script to see the different was to sample velocities for a given temperature T
@author: sr802
"""


import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell
Utilities_path=os.path.join(os.path.dirname(__file__), '../../')
sys.path.append(Utilities_path) #This falls into Utilities path
import Lammps.core_functions as cf


def box_muller(u1,u2):
  
  """
  Performs the box muller transformation
  """
  
  z1 = np.sqrt(-2*np.log(u1))*np.cos(2*np.pi*u2)
  z2 = np.sqrt(-2*np.log(u1))*np.sin(2*np.pi*u2)
  
  return z1,z2


def maxwell_dist(a,x):
    """
    My own maxwell distribution to see if everything is working
    it is as defined by https://mathworld.wolfram.com/MaxwellDistribution.html
    
    a is the scaling factor
    x are the velocities
    """
    
    dist = np.sqrt(2/np.pi)*x**2*np.exp((-x**2)/(2*a**2))/a**3
    
    return dist
    
    

# uniformly distributed values between 0 and 1
u1 = np.random.rand(10000)
u2 = np.random.rand(10000)
u3 = np.random.rand(10000)
u4 = np.random.rand(10000)
u5 = np.random.rand(10000)
u6 = np.random.rand(10000)




T = 1
m = 1 
scale = np.sqrt(T/m)


velocities =maxwell.rvs(size = 10000, scale = scale)


x = np.linspace(0,4)
pdist = maxwell.pdf(x)
mydist = maxwell_dist(scale,x)


# run the transformation
v1x,v2x = box_muller(u1,u2)
v1y,v2y = box_muller(u3,u4)
v1z,v2z = box_muller(u5,u6)


v1 =np.sqrt(v1x**2+v1y**2+v1z**2)


# Plotting
cf.set_plot_appearance()

plt.hist(scale*v1, density = True, label='Box Muller') 
plt.hist(velocities, density = True,label='scipy.maxwell')
plt.plot(x,mydist, label='Analytic function')
plt.plot(x,pdist,'.', label='scipy.maxwell')

plt.legend(loc = 'top right')

plt.show()