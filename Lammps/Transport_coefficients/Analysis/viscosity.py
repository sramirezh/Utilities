#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:09:58 2020
Script to compute the viscosity based on the stress autocorrelation function, it uses two outputs from lammps:
    


Based on GK transport mu_mu
@author: simon
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import uncertainties as un
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
from scipy.spatial.distance import pdist,squareform
from tqdm import tqdm
from uncertainties import unumpy,ufloat
import Others.Statistics.FastAverager as stat
from scipy import optimize
import glob
from joblib import Parallel, delayed
import multiprocessing
import pandas as pd
import Lammps.lammps_utilities as lu



import Lammps.Pore.GK.flux_correlation as fc

cwd = os.getcwd() #current working directory




class SimulationViscosity(lu.SimulationType):

    def __init__(self, name, max_tau, tau_integration,
                 log_file, pressure_file):
        """
        Args:
            name: Identifier of the simulation
            time_step: Simulation time_step
            max_tau: maximum length of the correlation analysis in fs
            tau_integration: upper time limit of the correlation integration
            temp: system temperature in Kelvin
            log_file: lammps log file
            flux file: file containing all the fluxes
            
            
        Other atributes:
            sampling_interval
        """
        self.name = name
        self.time_step = []
        self.max_tau = max_tau
        self.tau_integration = tau_integration
        self.temp = []
        self.log_file = log_file
        self.pressure_file = pressure_file
#        self.flux_file_prefix = self.flux_file.split('.')[0]
        self._get_thermo()
        self._get_prefactor()
        
    def run_analysis(self):

        data = cf.read_data_file("pressures.dat")
        
        data1 = data.values
        times = (data1[:,0]-data1[0,0])*self.time_step
        
        self.sampling_interval = times[1]-times[0] # Assuming times are homogeneous
        
        self.max_delta = int(self.max_tau/self.sampling_interval)  # int(len(times)*0.1) #Maximum delta of time to measure the correlation
        
        components = ["v_pxy", "v_pxz", "v_pyz"]
        
        p_off_diag = data[components].values
        
        #pxy = np.reshape(np.)
        
        p_flux = fc.flux(p_off_diag ,times,"P")
        
        
        etha11 = fc.correlation(p_flux, p_flux, self.max_delta)
        etha11.evaluate_acf()
        etha11.save("etha11")
        
        return etha11
    
    def _get_thermo(self):
        """
        Gets the thermo data as an instance of lu.SImulation
        TODO add this to the parent class
        
        """
        
        self.thermo = lu.Simulation(self.log_file)
    
    
    def _get_prefactor(self):
        """
        Obtains the relevant parameters from the simulation input and 
        converts from Lammps real units to SI
        """
        
        # Conversion factors
        atm2pa = 101325.0
        fs2s = 10**-15
        ang2m = 10**-10  
        
        Kb = 1.380649*10**-23 #J/K
        
        # Getting data 
        self.temp = lu.read_value_from(self.log_file, 'myTemp')[0]
        volume = self.thermo.volume
        # converting into SI
        
        factor = volume/(Kb * self.temp)
        
        scale = ang2m**3 * atm2pa**2 * fs2s
    
        self.prefactor = scale * factor
        

# =============================================================================
# Main
# =============================================================================

# Defining all the types of simulations

# High sampling frequency
N2 = SimulationViscosity("N2", 1000000,1000000, "log.lammps", "pressures.dat") 
Octane = SimulationViscosity("Octane", 200000, 200000, "log.lammps", "pressures.dat") 

cwd = os.getcwd()
plot_dir = "plots/0.viscosity"

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   

sim = Octane.copy

sim.plot_dir = plot_dir
sim.time_step = lu.read_value_from("input.lmp", 'myStep')[0]
sim.print_params(logger)

# Loading the correlations if they were already computed, this saves almost 50 minutes in 16 cores (dexter)
if (len(glob.glob('etha11*')) == 1):
    
    logger.info ("\nThere is a pkl, we need to load etha11!\n")
    etha11 = cf.load_instance("etha11.pkl")
    
else:
    logger.info ("There is no pkl file, so we need to analyse")
    etha11 = sim.run_analysis()
    

# =============================================================================
# Ploting all the correlations
# =============================================================================

# Loading lammps, on the fly calculations
lammps_df = cf.read_data_file("profile.gk.2d")
lammps_df = lammps_df.iloc[-400:]

plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()

fig,ax = etha11.plot_individual(fig, ax, norm = False, dim = 3, every = 10)

#ax.set_xlim(1000,2000)
#ax.set_ylim(-0.0005,0.0005)
#ax.plot(lammps_df["TimeDelta"].values[::10] * sim.time_step, 
#        lammps_df["v_pxy*v_pxy"].values[::10],label =r"$P_{xy}^{Lammps}$")
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')
ax.set_xscale('log')
ax.set_xlabel(r'$\tau[fs]$')

plt.legend( loc = 'upper right', fontsize = 12, ncol = 2)
ax.set_ylabel(r'$\langle %s(\tau) %s(0)\rangle$'%(etha11.flux1.name,etha11.flux2.name))
fig.tight_layout()
fig.savefig("%s/correlation11.pdf"%sim.plot_dir)

# =============================================================================
# # Plot eta vs tau
# =============================================================================

fig,ax = plt.subplots()

tau_array = np.linspace(etha11.times[0], etha11.times[-1], 20)
# Need to add the integration time
tau_array = np.sort(np.append(tau_array, sim.tau_integration)) 

eta_interval = []

# Computes the integral per intervals such that there is no need to recompute 
# Parts of the integral
#TODO be sure that the integration is continuos beween the intervals
initial_t = 0
for t in tau_array:
    eta_interval.append(sim.prefactor * etha11.transport_coeff_comp(initial_t, t, -1))
    initial_t = t
    

# Neef to rewrite the first term, very quick fix
eta_interval[0] = ufloat(0,0)

eta = np.cumsum(eta_interval)
eta_error = np.array([el.s for el in eta])
eta_array = np.array([el.n for el in eta])
    
ax.plot(tau_array,eta_array)
ax.fill_between(tau_array, eta_array - eta_error, eta_array + eta_error, alpha=0.4)
ax.set_xlabel(r'$t^*$')
#ax.axhline(y = eta_tau_int.n , xmin=0, xmax=1,ls='--',c='black', label = r'$\eta^* = %2.4f$' %eta_tau_int.n)
#ax.axvline(x = tau_integration,ls='-.',c='black')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlim(0, etha11.times[-1])
ax.set_ylim(0, ymax)
ax.set_ylabel(r'$\eta^*$')
plt.tight_layout()
plt.savefig("%s/eta_vs_tau.pdf"%sim.plot_dir)
#
#
#
#
#
#integral = etha11.transport_coeff(0, etha11.times[-1]) # In atmospheres**2*fs
#
#eta = integral * sim.prefactor # in Pa s
#
## It is a unumpy array
#total_eta = sum(eta)/len(eta)
#logger.info("The 3 components of the viscosity are %s" %eta)
#logger.info("The average viscosity is %2.4e +/- %2.4e [Pa s]"%(total_eta.n, total_eta.s))
#
#integral_lammps = cf.integrate(lammps_df["TimeDelta"].values*sim.time_step,lammps_df["v_pxy*v_pxy"].values,0,20000)
#
#eta_lammps = integral_lammps * sim.prefactor
#logger.info("On the fly estimation is %s" %eta_lammps)


