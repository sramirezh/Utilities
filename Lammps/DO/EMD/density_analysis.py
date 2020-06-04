#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 09:58:32 2020
Utilities to analyse density distributions

See an example of the use in DO/EMD/analytic_resutls.py
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


class PropertyDistribution(object):
    """
    TODO This could be the parent class of the timestep class in chunk utilities
    """

    def __init__(self, filename, log_file = "log.lammps", directory = [] ):
        """
        
        """
        self.file_name = filename
        self.directory = directory
        self.log_file = log_file
        self._initial_check()
        self._get_data()
        self.positions = self.data_frame['Coord1'].values
    
    def _initial_check(self):
        """
        Checks if the required files were given
        TODO check if the filename also exists
        """
        if not self.directory:
            self.directory = os.getcwd()
        
        # Adding the path to each file
        self.log_file = "%s/%s"%(self.directory, self.log_file)
        self.sim = lu.Simulation(self.log_file)

        self.file_name = "%s/%s"%(self.directory, self.file_name)

    def _get_data(self):
        self.data_frame = cf.read_data_file(self.file_name)

    @property
    def property_dist(self):
        """
        Retuns the names of all the properties given by the columns in the df 
        """
        return self.data_frame.columns
    
    def plot_property_dist(self, property_name, xlims=[], ylims=[], ax=[], save=True):
        """
        Plots a property distribution, to get the name of the properties, run
        self.property_dist.

        Args:
            property_names: as in the columns of the data_frame (self.property_dist)
            xlims: an array with the limits in x
            ylims:an array with the limits in y
            ax: axis 
            save: True or False to save the figure
        
        Returns:
            fig, ax: containing all the plot information
        """
        cf.set_plot_appearance()
        if not ax:
            fig, ax = plt.subplots()
        else:
            fig = plt.figure()
        ax.plot(self.positions, self.data_frame[property_name].values)
        #ax.set_ylabel(r'$C_s(y)-C_s^B$')
        #ax.set_xlabel(r'$y $')
        #xmin,xmax=plt.xlim()
        if xlims:
            ax.set_xlim(xlims[0], xlims[1])
        
        if ylims:
            ax.set_ylim(ylims[0], ylims[1])

        #ax.axhline(y=0, xmin=xmin, xmax=xmax,ls='--',c='black')
        fig.tight_layout()
        property_name = property_name.replace('/', '_')
        if save == True:
            plt.savefig("%s.pdf" % (property_name))
        
        return ax



class DensityDistribution(PropertyDistribution):
    """
    Child including the simulation inside and a possible bulk region definition
    and computing all the theoretical properties as per Anderson1984

    TODO if no bulk given, infer from the density distribution
    """
    
    def __init__(self, filename, bulk_name, log_file = "log.lammps", directory = []):
        """
        
        
        Attributes:
            properties: returns the name of the properties in the chunk
        """
        super().__init__(filename, log_file, directory)
        self.bulk_name = bulk_name
        self._set_bulk()
        self.rho_dist = self.data_frame['density/mass'].values
        self.rho_bulk = self.get_bulk_property('density/mass')
        self._lower_limit()
        
    def _lower_limit(self):
        """
        Estimates the position of the last chunk without any particles,
        Getting this for the entire fluid, gives the lower limit of the
        integrals
        """
        index = np.min(np.where(self.rho_dist > 0))
        if index > 0: 
            index = index - 1
        self.lower_limit = self.positions[index]
        
    def _set_bulk(self):
        """
        
        
        """
        self.h_b, self.limits_b = lu.read_region_height(self.bulk_name, geom_file = self.log_file)
    
    def get_bulk_property(self, name):
        """
        name is the column name, eg. density/mass
        """
        bulk_property = self.get_property_ave(name, self.limits_b)
        
        return bulk_property


    def get_property_ave(self, name, limits):
        """
        Get the property average between the limits
        name is the column name, eg. density/mass
        """
        
        data = self.data_frame[name].values
        indexes = np.where((self.positions> limits[0]) & (self.positions<limits[1]))
        
        ave_property = np.average(data[indexes])
        
        return ave_property
    
    def _get_rho_exc(self):
        """
        Excess concentration
        """
        self.rho_exc = self.rho_dist - self.rho_bulk
        self.data_frame['rho_exc'] = self.rho_exc
        return self.rho_exc
    
    
    def get_gamma(self, lower_limit = 0, upper_limit = []):
        """
        lower_limit: For the integration, has to be from the wall.
        upper_limit: If not assumed, it will be the position of the 
        """
        
        if not upper_limit:
            upper_limit = self.limits_b[0]
            
        self._get_rho_exc()

        # Transforming the positions
        positions = self.positions - lower_limit
        gamma = cf.integrate(positions, self.rho_exc, 0, upper_limit)
        
        self.gamma = gamma
        
        return gamma
    
    def get_k(self, lower_limit = 0, upper_limit = []):
        """
        lower_limit: For the integration, has to be from the wall.
        upper_limit: If not assumed, it will be the position of the 
        """
        
        if not upper_limit:
            upper_limit = self.limits_b[0]
            
        self._get_rho_exc()
        # Transforming the positions
        positions = self.positions - lower_limit

        integrand = self.rho_exc / self.rho_bulk

        self.data_frame["integrand_k"] = integrand
        k = cf.integrate(positions, integrand, 0, upper_limit)
        
        self.k = k
        
        return k
    
    
    def get_l(self, lower_limit = 0, upper_limit = []):
        """
        Computes l* as described by Anderson1984
        !!!!! The positions are all shifted, the zero is defined as per 
        lower_limit which is the position of the first bin with particles!!!
        lower_limit: For the integration, has to be from the wall.
        upper_limit: If not assumed, it will be the position of the 
        """
        
        if not upper_limit:
            upper_limit = self.limits_b[0]
            
        self._get_rho_exc()
        self.get_k(lower_limit)
        
        # Transforming the positions
        positions = self.positions - lower_limit

        # The positions below 0 have to be set to zero
        # Notice that this is only relevant for the integrand_first distribution
        # But not for the integral as the integral starts at zero

        index_negative = np.where(positions < 0)
        positions[index_negative] = 0
        integrand = positions * (self.rho_exc / self.rho_bulk)

        self.data_frame["integrand_first"] = integrand
        self.first_moment = cf.integrate(positions, integrand, 0, upper_limit) 
        self.l = self.first_moment / self.k

        return self.l

    # TODO this could be added to the class
    def _inner_integral(self, low_limit, sim):
        """
        Computes the inner integral in equation (16) Anderson1984, from the lower limit
        to the beginning of the bulk as given by species.limits_b[0]
        Args:
            low_limit: The lower limit for the integral
            species: instance of the class da.DensityDistribution
        
        Returns:
            integral: the integral from the lower limit to the the start of the bulk 
        """
        integrand = self.data_frame['integrand_k'].values / sim.eta
        integral = cf.integrate(self.positions, integrand , low_limit, self.limits_b[0])
        
        return integral


    def vx(self, z, sim, grad_c, low_limit = []):
        """
        Computes the velocity in the x direction at the heigth z from the wall,
        as per equation (16) Anderson1984
        
        Args:
            z: Heigth from the wall in which the velocity is computed
            sim: Instance of the class SimulationEMD
            grad_c: Concentration gradient of the species.
            low_limit: Lower limit of the integral
            
        Returns:
            The velocity at the given position
            
            
        """
        if not low_limit:
            low_limit = self.lower_limit
        
        if self.data_frame['vx_integrand'].empty:
            self.data_frame['vx_integrand'] = [self._inner_integral(x,sim) for x in self.positions]
        integral = cf.integrate(self.positions, self.data_frame['vx_integrand'], low_limit, z) 
        vx_z = - (sim.T) * grad_c * integral 
        return    vx_z
    
    def vx_dist(self, sim, grad_c, low_limit = []):
        """
        Computes the velocity distribution in the x direction as a function
        of the position in z
            Args:
            z: Heigth from the wall in which the velocity is computed
            sim: Instance of the class SimulationEMD
            grad_c: Concentration gradient of the species.
            low_limit: Lower limit of the integral
        """
        self.data_frame['vx_integrand'] = [self._inner_integral(x,sim) for x in self.positions]
        self.data_frame['vx_z'] = [self.vx(z, sim, grad_c, low_limit) for z in self.positions]
        self.vx_bulk = self.data_frame['vx_z'].values[-1]
        return self.data_frame['vx_z'].copy()

    
    def compute_all_properties(self, lower_limit, upper_limit = []):
        
        if not upper_limit:
            upper_limit = self.limits_b[0]
        
        self.get_gamma(lower_limit, upper_limit) 
        self.get_k(lower_limit, upper_limit) 
        self.get_l(lower_limit, upper_limit) 
        
    def compute_ladm(self, sigma, lower_limit = 0, upper_limit = [], epsilon = 0.00001):
        """
        Computes the Local Average Density Model (LADM) as in Bitsanis1988
        but in the direction perpendicular to z.
        sigma is the difference between the integral limits in eq 2.2
        lower_limit: for the LADM application, in general starting from the wall, i.e =0
        upper_limit: for the LADM application, after this, the density gets 
        the values from self.rho_dis.
        epsilon: to make open intervals in the limits of the integration.
        """
        if not upper_limit:
            upper_limit = self.limits_b[0]

        indexes = cf.get_interval(self.positions, lower_limit, upper_limit)
        ladm = self.rho_dist.copy()

        for index in indexes:
            position = self.positions[index]  # point where ladm is computed
            left_limit = position - 0.5 * sigma + epsilon 
            right_limit = position + 0.5 * sigma - epsilon
            indexes_2 = cf.get_interval(self.positions, left_limit, right_limit )

            # Length of the integration interval
            delta_z = self.positions[indexes_2[-1]] - self.positions[indexes_2[-0]]
            ladm[index] = cf.integrate(self.positions, self.rho_dist, 
                                    left_limit, right_limit) / delta_z

        self.data_frame['ladm'] = ladm
        return ladm
            
    def print_summary(self):
        return 0
    
    
    
    
        
