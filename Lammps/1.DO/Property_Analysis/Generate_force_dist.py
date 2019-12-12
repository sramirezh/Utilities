#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
New script that uses the density distributions to plot the force profiles and theoretical quantities as gamma, K,L, and V.
TODO this code assumes the surface at z = 0
File created copying what is in Property analysis but just to generate the force distribution
@author: sr802
"""

import numpy as np
import sys
import os
from shlex import split
import matplotlib.pyplot as plt
import linecache
import re
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.General.Log_Analysis.Thermo_Analyser as ta
import Lammps.lammps_utilities as lu


FProperties=cf.read_data_file("Fproperties_short.dat").values #Solvent properties
SProperties=cf.read_data_file("Sproperties_short.dat").values #Solute properties
AProperties=cf.read_data_file("properties_short.dat").values #All properties




# TODO this funtion could be replaced with the one that I used in the colloid/polymer theoretical mobility
def PosConstant(x,y,Tol):
    """
    Finds where a Force distribution type function no longer changes under the given Tolerance, analysing
    the points after the minimum and the maximum of the function. The chosen
    point corresponds to the minimum in the set of points


    :param x:
    :param y:
    :param Tol: Tolerance
    :return: The index of the position of the last point to be taken into account
    """


    BinSize = x[1] - x[0]
    der1=np.gradient(y,BinSize)
    IndexMax=np.argmax(y)
    IndexMin=np.argmin(y)
    ExtremeIndex=max(IndexMax,IndexMin)
    ExtremeIndex = max(IndexMax, IndexMin)
    Index_c = np.where(np.abs(der1) < Tol)[0]  # Indexes where the function is approximately constant
    LastIndex = np.where(Index_c > ExtremeIndex)[0][0]
    index = Index_c[LastIndex]
    return index


def Force(Ns,Nf,SProperties,FProperties,AProperties):
    """
    Force calculation
    :param Ns: Solute Number in the Bulk
    :param Nf: Solvent Number in the Bulk
    :param SProperties: Solute Properties
    :param FProperties: Solvent Properties
    :param AProperties: Fluid Properties
    :return:
    The Force distribution 

    """

    Fs = 1
    Ff = -Fs * Ns / Nf
    n,m = np.shape(FProperties)
    Force = np.zeros((n, 2))
    Force[:, 0] = FProperties[:, 1]

    for i in range(n):
        if AProperties[i, 4] == 0:
            Force[i, 1] = 0

        else:
            Force[i, 1] = (Ff * FProperties[i, 4] + Fs * SProperties[i, 4]) / AProperties[i, 4]

    np.savetxt("Force.dat", Force)
    print("\n*********************Getting the force profile*********************\n")
    print("The force on the Solutes is %f, on the Solvents %f" %(Fs,Ff))
    print("Created the force distribution File Force.dat ")

    return Force





# Todo put the read_box_limits and the read_bulk heigh as general functions that can read from the desired lin
    #the numbers of lines and extract the digits from there. 


# Obtaining the thermodynamic data


thermo_data = ta.thermo_analyser("log.lammps")
volume, limits  = lu.read_box_limits("log.lammps")


Ns = float(thermo_data.at['v_cBSolu','Average'])
Nf = float(thermo_data.at['v_cBSolv','Average'])




# =============================================================================
#  Getting the concentration in the bulk, using the results from log. lammps
# =============================================================================

#Getting the bulk limits
h_B, limits_B = lu.read_region_height("rBulk")

# Getting if its 2d or 3d
input_name = "input.lmp"
out,err = cf.bash_command("""grep -n "enforce2d" %s | awk -F":" '{print $1}' """%input_name)

if out:
    print("\nThis is a 2d Simulation\n")
    v_B=(volume*(h_B/(limits[1][1]-limits[0][1])))/0.5 #Added the 0.5 because of the 2D correction
    index  = 1
    
else:
    print("\nThis is a 3d Simulation\n")
    v_B=(volume*(h_B/(limits[1][2]-limits[0][2]))) 
    index = 2

Cs_B = Ns/v_B
Cf_B = Nf/v_B


Force = Force(Ns,Nf,SProperties,FProperties,AProperties)


# =============================================================================
# Finding the cutoff for the force
# =============================================================================

# To find the plateau after the maximum, as before there might be a plateau
imax =np.argmax(Force[:,1])

cut_off = imax+min(cf.plateau_finder(Force[imax::,1]))[-1]

zmax_force = SProperties[cut_off,1]






# =============================================================================
# Creating the force distribution
# =============================================================================
Zpos = Force[:cut_off+1, 0]
MuF = np.transpose(Force[:cut_off, 1])
print("The Force Cut-off is %f, this is where the region of applied forces finishes"%np.max(Zpos))
print("Creating the Files to iterate in Lammps")
print ("In the bulk, the solute concentration is %s, the solvent is %s"%(Cs_B,Cf_B))
np.savetxt("Zpos_iterate.dat", Zpos)
np.savetxt("Force_iterate.dat", MuF)



# =============================================================================
# Computing important quantities
# =============================================================================
# TODO make this a class 

z = SProperties[:,1]

Cs_exc=SProperties[:,4]-Cs_B

gamma=cf.integrate(z,Cs_exc,0,zmax_force)

integrand_k=Cs_exc/Cs_B

K=cf.integrate(z,integrand_k,0,zmax_force)

integrand_1=integrand_k*z

L=cf.integrate(z,integrand_1,0,zmax_force)

integrand_2=0.5*integrand_k*z**2

H=cf.integrate(z,integrand_2,0,zmax_force)/L


def velocity(Cs_exc, zmax):
    """
    Computes the velocity for each position from the wall
    
    
    
    """
    z = np.linspace(0,zmax)
    
    internal_values = [cf.integrate(z,Cs_exc,0,zlim) for zlim in z]
    integral = cf.integrate(z,internal_values,0,zlim)
    return internal_values


def bulk_velocity(Cs_exc,z):
    """
    Computes the bulk velocity, should be 0.012
    """
    eta = 1.56
    grad_mu = - 0.125
    
    integrand = Cs_exc*z
    
    integral = cf.integrate(z,integrand,0,25)
    
    vx = grad_mu/eta*integral
    
    return vx
    
print ("\nThe bulk velocity is %f"%bulk_velocity(Cs_exc,z))
print ("Gamma is %f\n K is %f\n L is %f\n H is %f \n"%(gamma,K,L,H))
    

# =============================================================================
# Plots
# =============================================================================

#Plots without shifting
zmin = 0
zmax = limits[1][index]
x_tick_distance = 5

#cf.set_plot_appearance()


# Force distribution
fig1,ax1=plt.subplots()
ax1.plot(Force[:,0],Force[:,1])
ax1.set_ylabel(r"$F$")
ax1.set_xlabel(r'$d[\sigma]$')
ax1.set_xlim(zmin, zmax)
ax1.set_xticks(np.arange(zmin, zmax, x_tick_distance))
ax1.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
ax1.axvline(x=Zpos[-1], ymin=0, ymax=1,ls='-.',c='black')
ax1.axvline(x=10, ymin=0, ymax=1,ls='-.',c='b')
ax1.axvline(x=20, ymin=0, ymax=1,ls='-.',c='b')
plt.tight_layout()
fig1.savefig("Force_dist.pdf")


#Excess Solute concentration
fig,ax=plt.subplots()
ax.plot(SProperties[:,1],Cs_exc)
ax.set_ylabel(r'$C_s(y)-C_s^B$')
ax.set_xlabel(r'$y $')
xmin,xmax=plt.xlim()
ax.set_xlim(0,zmax)
ax.axhline(y=0, xmin=xmin, xmax=xmax,ls='--',c='black')
fig.tight_layout()
plt.savefig("Excess.pdf")

# Density distribution

fig2,ax2=plt.subplots()
ax2.plot(FProperties[:,1],FProperties[:,4],label = 'Solvents', c ='b')
ax2.plot(SProperties[:,1],SProperties[:,4],label = 'Solutes', c ='r')
ax2.set_ylabel(r"$c[\sigma^{-3}]$")
ax2.set_xlabel(r'$d[\sigma]$')
ax2.set_xlim(zmin, zmax)
ax2.set_ylim(0,ax2.get_ylim()[-1])
ax2.set_xticks(np.arange(zmin, zmax, x_tick_distance))

ax2.axhline(y=Cs_B, xmin=0, xmax=1,ls=':',c='black')
ax2.axhline(y=Cf_B, xmin=0, xmax=1,ls=':',c='black')
plt.legend(loc = 'upper right')
plt.tight_layout()
fig2.savefig("density_dist.pdf")
