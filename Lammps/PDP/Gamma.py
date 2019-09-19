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
Data=pd.read_csv(file_name,sep=" ",dtype=np.float64,header=3).values #3 Columns, Position, Count, Density
bulk_min=5
bc=bulk_concentration(Data[:,0],Data[:,2],bulk_min)

xmin=0.0
xmax=8.0




gamma_integrand=4*np.pi*(Data[:,2]/bc-1.0)*Data[:,0]**2
ls_integrand=4*np.pi*(Data[:,2]/bc-1.0)*Data[:,0]**3

gamma=Integrate(Data[:,0],gamma_integrand,xmin,xmax)
ls=(Integrate(Data[:,0],ls_integrand,xmin,xmax)/gamma)


print("Gamma %lf"%gamma)

print("Legth %lf"%ls)
print("Bulk Concentration %lf" %bc)

1
#Test

gamma_test_integrand=(Data[:,2]/bc-1.0)
gamma_test=Integrate(Data[:,0],gamma_test_integrand,xmin,xmax)
ls_integrand_test=(Data[:,2]/bc-1.0)*Data[:,0]
ls_test=(Integrate(Data[:,0],ls_integrand_test,xmin,xmax)/gamma_test)

print("Gamma_test %lf"%gamma_test)
print("ls_test %lf"%ls_test)

#Values for N_30 E_1.5_1.5

mobility=0.25428875 #0.1867325
rg=5.74288222199 #3.96062536265
k=1.5*rg**2


viscosity=2.061 #from Molecular Dynamics Study of the Lennard-Jones Fluid Viscosity: Application to Real Fluids
2
mobility_theo=bc*gamma*(ls+k**0.5)/viscosity
mob_test=-bc*gamma_test*(ls_test+k**0.5)/viscosity
print("mobility theo %lf %lf"%(mobility_theo,mob_test))



#import matplotlib.pyplot as plt
#plt.plot(Data[:50,0],Data[:50,2])
