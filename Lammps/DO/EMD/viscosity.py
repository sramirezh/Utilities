#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:09:58 2020
Script to compute the viscosity based on the stress autocorrelation function, 
it uses two outputs from lammps:

    pressure.dat that contains the off diagonal pressures "v_pxy", "v_pxz", "v_pyz"
    S0St.dat contains lammps on-the-fly stress autocorrelation function
    


Based on viscosity from transport_coefficient project
@author: simon
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import uncertainties as un
import glob
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.lammps_utilities as lu
import Lammps.Pore.EMD.flux_correlation as fc


def run_analysis(press_file, time_step, max_tau):
    """
    File containing the pressure
    """

    data = cf.read_data_file(press_file)
    
    data1 = data.values
    times = (data1[:,0]-data1[0,0]) * time_step
    
    
    
    sampling_interval = times[1]-times[0] # Assuming times are homogeneous
    
    
    
    max_delta = int(max_tau/sampling_interval)  # int(len(times)*0.1) #Maximum delta of time to measure the correlation
    
    components = ["v_pxy", "v_pxz", "v_pyz"]
    
    p_off_diag = data[components].values
    
    #pxy = np.reshape(np.)
    
    pxy_flux = fc.flux(p_off_diag ,times,"P")
    
    
    eta11 = fc.correlation(pxy_flux,pxy_flux, max_delta)
    eta11.evaluate_acf()
    eta11.save("eta11")
    
    return eta11

# =============================================================================
# Main
# =============================================================================
logger = cf.log(__file__, os.getcwd())

time_step = 0.005 # In LJ time
max_tau = 100 # in LJ time  Max tau to analyse
tau_integration = 10 # limit for the integration
Kb = 1  
Temperature = 1
log_file = "log.lammps"
press_file = "pressures.dat"
lammps_sacf = "S0St.dat"

logger.info("Using delta_t = %s " %time_step)
logger.info("Using max tau = %s for the analysis "%max_tau)
logger.info("Using integration tau = %s for the autocorrelation "%tau_integration)
logger.info("Getting the box size from %s" %log_file)

volume, limits = lu.read_box_limits("log.lammps", -1)

prefactor = volume / (Kb * Temperature)

# =============================================================================
# Computing the viscosity from the stress data
# =============================================================================

# Loading the correlations if they were already computed, this saves almost 50 minutes in 16 cores (dexter)
if (len(glob.glob('eta11*')) == 1):
    
    logger.info("\nThere is a pkl, we need to load eta11!\n")
    eta11 = cf.load_instance("eta11.pkl")
    
else:
    logger.info("There is no pkl file, so we need to analyse")
    eta11 = run_analysis(press_file, time_step, max_tau)
    
    
integral = eta11.transport_coeff(0, tau_integration) 
eta_tau_int = integral * prefactor 

logger.info("The viscosity is %2.5f+/- %2.5f in LJ units"%(eta_tau_int[-1].n, eta_tau_int[-1].s))

# =============================================================================
# Computing the viscosity from lammps on-the-fly correlation
# =============================================================================
lammps_df = cf.read_data_file("S0St.dat")
lammps_df = lammps_df.iloc[-400:]
integral_lammps1 = cf.integrate(lammps_df["TimeDelta"].values*time_step,lammps_df["v_pxy*v_pxy"].values,0,max_tau)
integral_lammps2 = cf.integrate(lammps_df["TimeDelta"].values*time_step,lammps_df["v_pxz*v_pxz"].values,0,max_tau)
integral_lammps3 = cf.integrate(lammps_df["TimeDelta"].values*time_step,lammps_df["v_pyz*v_pyz"].values,0,max_tau)

eta_lammps = (integral_lammps1 + integral_lammps2 + integral_lammps3) * prefactor/3

logger.info("The viscosity from LAMMPS is %2.5f in LJ units"%(eta_lammps))
 
# =============================================================================
# Ploting all the correlations
# =============================================================================

plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()

ax = eta11.plot_all(ax, norm = False)

#ax.set_xlim(1000,2000)
#ax.set_ylim(-0.0005,0.0005)
#ax.plot(lammps_df["TimeDelta"]*time_step, lammps_df["v_pxy*v_pxy"],label =r"$P_{xy}^{Lammps}$")

ax.axvline(x = tau_integration,ls='-.',c='black')
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xscale('log')
xmin,xmax = ax.get_xlim()
ax.set_xlim(xmin, max_tau)
ax.set_xlabel(r'$t^*$')

plt.legend([r'$xy$',r'$xz$',r'$yz$',"Total"], loc = 'upper right', fontsize = 12, ncol = 2)
ax.set_ylabel(r'$\langle \Pi(t) \Pi(0)\rangle$')
plt.tight_layout()
plt.savefig("correlation11.pdf")

# =============================================================================
# # Plot eta(average in all directions) vs tau
# =============================================================================

fig,ax = plt.subplots()

tau_array = np.linspace(eta11.times[1], eta11.times[-1])
# Need to add the integration time
tau_array = np.sort(np.append(tau_array, tau_integration)) 

eta_array = []
eta_error = []
plotting_component = -1

for t in tau_array:
    # Taking only the nominal value
    eta = prefactor * eta11.transport_coeff_comp(0, t,plotting_component)
    error = 0
    if isinstance(eta, un.UFloat):
        error = eta.s
        eta = eta.n
        
    eta_error.append(error)
    eta_array.append(eta)

eta_error = np.array(eta_error)
eta_array = np.array(eta_array)
    
ax.plot(tau_array,eta_array)
ax.fill_between(tau_array, eta_array - eta_error, eta_array + eta_error, alpha=0.4)
ax.set_xlabel(r'$t^*$')
ax.axhline(y = eta_tau_int[plotting_component].n , xmin=0, xmax=1,ls='--', 
           c='black', label = r'$\eta^* = %2.4f$' %eta_tau_int[plotting_component].n)
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlim(0, eta11.times[-1])
ax.set_ylim(0, ymax)
ax.set_ylabel(r'$\eta^*$')
plt.legend( loc = 'lower right', fontsize = 20)
plt.tight_layout()
plt.savefig("eta_vs_tau.pdf")
