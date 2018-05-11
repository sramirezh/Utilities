#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 15:13:58 2018
Computes the solute adsorption in the radial direction for a monomer
reading the prof_u
@author: sr802
"""

import numpy as np
import pandas as pd


# =============================================================================
# Function Definition
# =============================================================================

def Integrate(x,y,xmin,xmax):
    """
    Integrate the data in x and y from xmin to xmax
    """
    MinIndex=np.min(np.where(x>=xmin))
    MaxIndex=np.max(np.where(x<=xmax))
    I=np.trapz(y[MinIndex:MaxIndex],x[MinIndex:MaxIndex])

    return I


def bulk_concentration(Positions,Densities,bulk_min):
    """
    Args:
        Position vector.
        Density vector, for each position.
        bulk_min, defined by taking a look at the density distribution.

    Returns:
        bulk_c the average concentration in the bulk
    """

    bulk_i=np.where(Positions>=bulk_min)[0][:-1] #The last to avoid taking the last density that explodes
    bulk_c=np.average(Densities[bulk_i])
    return bulk_c

file_name="prof_u.dat"
Data=pd.read_csv(file_name,sep=" ",dtype=np.float64,header=3).as_matrix() #3 Columns, Position, Count, Density
bulk_min=5
bc=bulk_concentration(Data[:,0],Data[:,2],bulk_min)

xmin=0.0
xmax=8.0
Integrand=4*np.pi*(Data[:,2]/bc-1.0)*Data[:,0]**2

SolAbso=Integrate(Data[:,0],Integrand,xmin,xmax)

print "Gamma %lf"%SolAbso

#import matplotlib.pyplot as plt
#plt.plot(Data[:50,0],Data[:50,2])
