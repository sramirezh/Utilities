#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 09:58:32 2020
Script to analyse density distributions
@author: simon
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
import Lammps.General.Thermo_Analyser as ta
import Lammps.lammps_utilities as lu


class DensityDistribution(object):
    """
    TODO This could be the parent class of the timestep class in chunk utilities
    """

    def __init__(self, filename, log_file = []):
        """
        
        """
        self.filename = filename
        self._get_data()
        self.log_file = []
        self._initial_check()
        self._properties()
    
    def _initial_check(self):
        """
        Checks if the required files were given
        TODO check if the filename also exists
        """
        if len(self.log_file) !=0:
            self.sim = lu.Simulation(self.log_file)


    def _get_data(self):
        self.data_frame = cf.read_data_file(self.filename)
        
    def _properties(self):
        self.properties = self.data_frame.columns



class DDWithBulk(DensityDistribution):
    """
    Child including the simulation inside and a possible bulk region definition
    """
    
    def __init__(self, filename, bulk_name, log_file = []):
        """
        
        
        Attributes:
            properties: returns the name of the properties in the chunk
        """
        super().__init__(filename, log_file)
        self.bulk_name = bulk_name
        self._set_bulk()




    def _get_data(self):
        self.data_frame = cf.read_data_file(self.filename)
        
    
    def _set_bulk(self):
        """
        
        
        """
        self.h_b, self.limits_b = lu.read_region_height(self.bulk_name,geom_file = "log.lammps")
    
    def get_bulk_property(self, name):
        """
        name is the column name, eg. density/mass
        """
        positions = self.data_frame['Coord1'].values
        data = self.data_frame[name].values
        
        indexes = np.where((positions> self.limits_b[0]) & (positions<self.limits_b[1]))
        
        ave_property = np.average(data[indexes])
        
        return ave_property
        
    

solute = DDWithBulk("Sproperties_short.dat",'rBulk')


#
## =============================================================================
## Main
## =============================================================================
#
#sim = lu.Simulation("log.lammps")
#

#solvent = DensityDistribution("Fproperties_short.dat")
#fluid = DensityDistribution("properties_short.dat")
#
#ns = sim.thermo_ave.at['v_cBSolu','Average']
#nf = sim.thermo_ave.at['v_cBSolv','Average']
#
## =============================================================================
##  Getting the concentration in the bulk, using the results from log. lammps
## =============================================================================
#
#h_b, limits_b = lu.read_region_height("rBulk") # Getting bulk limits
#
#if sim.is_2d:
#    print("\nThis is a 2d Simulation\n")
#    
#    #Added the 0.5 because of the 2D correction
#    v_b = (sim.volume * (h_b / (sim.limits[1][1] - sim.limits[0][1]))) / 0.5 
#    index  = 1
#    
#else:
#    print("\nThis is a 3d Simulation\n")
#    v_b = (sim.volume * (h_b / (sim.limits[1][2] - sim.limits[0][2]))) 
#    index = 2
#
## TODO this can be computed from the distribution, knowing where is the bulk
## Concentration of species in the bulk
#cs_b = ns / v_b
#cf_b = nf / v_b
#
#
#
#
## =============================================================================
## Computing important quantities
## =============================================================================
## TODO make this a class 
#
#z = SProperties[:,1]
#
#Cs_exc=SProperties[:,4]-Cs_B
#
#gamma=cf.integrate(z,Cs_exc,0,zmax_force)
#
#integrand_k=Cs_exc/Cs_B
#
#K=cf.integrate(z,integrand_k,0,zmax_force)
#
#integrand_1=integrand_k*z
#
#L=cf.integrate(z,integrand_1,0,zmax_force)
#
#integrand_2=0.5*integrand_k*z**2
#
#H=cf.integrate(z,integrand_2,0,zmax_force)/L
#
#
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
#def bulk_velocity(Cs_exc,z):
#    """
#    Computes the bulk velocity, should be 0.012
#    """
#    eta = 1.56
#    grad_mu = - 0.125
#    
#    integrand = Cs_exc*z
#    
#    integral = cf.integrate(z,integrand,0,25)
#    
#    vx = grad_mu/eta*integral
#    
#    return vx
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
