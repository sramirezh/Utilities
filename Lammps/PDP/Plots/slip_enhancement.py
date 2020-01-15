#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:02:17 2020


Fitting two scattered lines that are assumed to be just shifted by a constant value.

See some results in the logbook Colloid January 15 


@author: sr802
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf

from scipy import optimize




def fitfunc(p,x):
    """
    This gets the value of the theory and fit them to the simulation, with -p being the fitting parameter
    """
    index = np.where(x_vect == x)[0]
    y = theory[index]*p   # theory - fitting parameter that is a displacement

    return y
    
    
    
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err #To include the error in the least squares

theory = 10*np.random.rand(10) # Theory values, with no errors
x_vect = np.linspace(0,10,10)
noise = np.random.normal(0,1,10) # to be added to the simulation and to account for error, just for simplicity

displacement = 3
simulation =  (theory*displacement)+noise





pinit = 0 
out = optimize.leastsq(errfunc, pinit, args=(x_vect, simulation, noise), full_output=1)
#cov=out[1] #Covariance in the
pfinal = out[0] #fitting coefficients



#Plot 

cf.set_plot_appearance()

fig,ax = plt.subplots()



ax.plot(x_vect, theory, label='Theory')
ax.errorbar(x_vect, simulation ,yerr = noise , fmt='o', label='Simulations')
ax.plot(x_vect, theory*pfinal[0], label='Fitted theory')  
#ax.axvline(x=L[0]/2, ymin = 0, ymax=1, ls=':',c='black')
#ax.set_xlim(0, L[0])   
#ax.axvspan(L[0]/2-r_colloid,L[0]/2+r_colloid, alpha=0.5, color='green')
#ax.axvspan(6*lattice_constant,9*lattice_constant, alpha=0.5, color='blue')
#ax.axvspan(21*lattice_constant,24*lattice_constant, alpha=0.5, color='red')
#
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$v_x(x)$")
ax.legend()
plt.tight_layout()
plt.savefig("best_fit.pdf")


print ("The shift factor is %s"%pfinal)
