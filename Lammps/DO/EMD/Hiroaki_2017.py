#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 20 13:42:54 2020
Script to  validate the results by Hiroaki et al 2017
@author: simon
"""

import numpy as np
import sys
import os
from shlex import split
import matplotlib.pyplot as plt
import linecache
import re
import density_analysis as da
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

# =============================================================================
# Main
# =============================================================================

def bulk_velocity(Cs_exc,z):
    """
    Computes the bulk velocity, should be 0.012
    """
    eta = 1.56
    grad_mu = - 0.125
    
    integrand = Cs_exc*z
    
    integral = cf.integrate(z,integrand,0,10)
    
    vx = -grad_mu/eta*integral
    
    return vx


# =============================================================================
# Integrating with other techniques 
# =============================================================================
    
import scipy.integrate as integrate




solute_hiroaki = da.DensityDistribution("Hiroaki.dat",'rBulk')
solute = da.DensityDistribution("Sproperties_short.dat",'rBulk')
solvent = da.DensityDistribution("Fproperties_short.dat",'rBulk')
fluid = da.DensityDistribution("properties_short.dat",'rBulk')

# Need the lower limit from all the solution
lower_limit = fluid.lower_limit

z_wall = 0

z_min = lower_limit

print("\nAssuming that the wall is at z = %s"%z_min)

solute.compute_all_properties(solute.lower_limit)
solute_hiroaki.compute_all_properties(solute_hiroaki.lower_limit)
solvent.compute_all_properties(solvent.lower_limit)

Kb = 1
T = 1
eta =  1.57
grad_mu_s = - 0.125

# For the solutes
grad_c_s = grad_mu_s * solute.rho_bulk
v_0_s = - (Kb* T)* solute.l*solute.k/eta
vx_s = v_0_s * grad_c_s


# For the solutes Hiroaki
grad_c_s_h = grad_mu_s * solute_hiroaki.rho_bulk
v_0_s_h = - (Kb* T)* solute_hiroaki.l*solute_hiroaki.k/eta
vx_s_h = v_0_s_h * grad_c_s_h

print("The velocity for hiroaki is %s"%vx_s_h)

# For the solvents
grad_mu_f = - (solute.rho_bulk / solvent.rho_bulk)* grad_mu_s
grad_c_f = grad_mu_f * solvent.rho_bulk
v_0_f = - (Kb* T)* solvent.l*solvent.k/eta
vx_f = v_0_f * grad_c_f



vx_total = vx_s + vx_f


print (vx_s)
print (vx_f)


print ("total is %s" % (vx_total))

vx_2 = bulk_velocity(solute.rho_exc, solute.positions)

print (vx_2)



# =============================================================================
# Computing the velocity using both species
# =============================================================================
print("The property distributions from the solutes are:" )
print(solute.property_dist)
#solute.plot_property_dist()


solute.get_k(0)
solute.plot_property_dist("integrand_k", 0,8)

#integrand = solute.exc_con * solute.positions



#def velocity(Cs_exc, zmax):
#    """
#    Computes the velocity for each position from the wall
#    
#    
#    
#    """
#    z = np.linspace(0,zmax)
#    
#    internal_values = [cf.integrate(z,Cs_exc,0,zlim) for zlim in z]
#    integral = cf.integrate(z,internal_values,0,zlim)
#    return internal_values
#
#


#    
#print ("\nThe bulk velocity is %f"%bulk_velocity(Cs_exc,z))
#print ("Gamma is %f\n K is %f\n L is %f\n H is %f \n"%(gamma,K,L,H))
#    
#
#
#
#
#
#
##Excess Solute concentration
#fig,ax=plt.subplots()
#ax.plot(SProperties[:,1],Cs_exc)
#ax.set_ylabel(r'$C_s(y)-C_s^B$')
#ax.set_xlabel(r'$y $')
#xmin,xmax=plt.xlim()
#ax.set_xlim(0,zmax)
#ax.axhline(y=0, xmin=xmin, xmax=xmax,ls='--',c='black')
#fig.tight_layout()
#plt.savefig("Excess.pdf")
#
## Density distribution
#
#fig2,ax2=plt.subplots()
#ax2.plot(FProperties[:,1],FProperties[:,4],label = 'Solvents', c ='b')
#ax2.plot(SProperties[:,1],SProperties[:,4],label = 'Solutes', c ='r')
#ax2.set_ylabel(r"$c[\sigma^{-3}]$", fontsize = 30)
#ax2.set_xlabel(r'$d[\sigma]$', fontsize = 30)
#ax2.set_xlim(zmin, zmax)
#ax2.set_ylim(0,ax2.get_ylim()[-1])
#ax2.set_xticks(np.arange(zmin, zmax, x_tick_distance))
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#ax2.axhline(y=Cs_B, xmin=0, xmax=1,ls=':',c='black')
#ax2.axhline(y=Cf_B, xmin=0, xmax=1,ls=':',c='black')
#plt.legend(loc = 'upper right', fontsize = 15)
#plt.tight_layout()
#fig2.savefig("density_dist.pdf")

