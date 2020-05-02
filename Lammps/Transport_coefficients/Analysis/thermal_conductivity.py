#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:09:58 2020
Script to compute the thermal conductivity based on the sflux autocorrelation function

The flux is computed in lammps as:
    
compute      myKE all ke/atom
compute      myPE all pe/atom
compute      myStress all stress/atom NULL virial
compute      flux all heat/flux myKE myPE myStress
variable     Jx equal c_flux[1]/vol


Based on GK transport mu_mu
@author: simon
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.transformations as tr
import MDAnalysis.analysis.rdf as rdf
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



import Lammps.Pore.GK.flux_correlation as fc

cwd = os.getcwd() #current working directory

def run_analysis():

    data = cf.read_data_file("fluxes.dat")
    
    data1 = data.values
    times = (data1[:,0]-data1[0,0])*time_step
    
    
    
    sampling_interval = times[1]-times[0] # Assuming times are homogeneous
    
    
    
    max_delta = int(max_tau/sampling_interval)  # int(len(times)*0.1) #Maximum delta of time to measure the correlation
    
    components = ["v_Jx", "v_Jy", "v_Jz"]
    
    j_components = data[components].values
    
    #pxy = np.reshape(np.)
    
    j_heat = fc.flux(j_components ,times,"J")
    
    
    kappa = fc.correlation(j_heat, j_heat, max_delta)
    kappa.evaluate()
    kappa.save("kappa")
    
    return kappa







# =============================================================================
# Main
# =============================================================================



# =============================================================================
# # Input parameters
# =============================================================================

time_step = 10 # In femptoseconds
box_side =  400.57
Temperature =  273.15 #K 
max_tau = 1000000
lammps_d = 400 # See lammps code $p*$s

print ("\nUsing a box side of %s"%box_side)
print ("Using delta_t = %s fs" %time_step)
print ("Using max tau = %s fs"%max_tau)
print ("Temperature = %s fs" %Temperature)
volume = box_side**3 # Angs**3 

# Loading the correlations if they were already computed, this saves almost 50 minutes in 16 cores (dexter)
if (len(glob.glob('kappa*')) == 1):
    
    print ("\nThere is a pkl, we need to load kappa!\n")
    kappa = cf.load_instance("kappa.pkl")
    
else:
    print ("There is no pkl file, so we need to analyse")
    kappa = run_analysis()


# =============================================================================
# Ploting all the correlations
# =============================================================================


lammps_df = cf.read_data_file("J0Jt.dat")
lammps_df = lammps_df.iloc[-lammps_d:]


plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()

fig,ax = kappa.plot_all(fig, ax, norm = False)

#ax.set_xlim(1000,2000)
#ax.set_ylim(-0.0005,0.0005)
ax.plot(lammps_df["TimeDelta"]*time_step, lammps_df["c_flux[1]*c_flux[1]"],label ="JxLammps")
ax.plot(lammps_df["TimeDelta"]*time_step, lammps_df["c_flux[2]*c_flux[2]"],label ="JyLammps")
ax.plot(lammps_df["TimeDelta"]*time_step, lammps_df["c_flux[3]*c_flux[3]"],label ="JzLammps")


ax.set_xscale('log')
ax.set_xlabel(r'$\tau[fs]$')

plt.legend( loc = 'upper right', fontsize = 12, ncol = 2)
ax.set_ylabel(r'$\langle %s(\tau) %s(0)\rangle$'%(kappa.flux1.name,kappa.flux2.name))
plt.tight_layout()
plt.savefig("correlation11.pdf")





# converting into SI
Kb = 1.380649*10**-23 #J/K
kCal2J = 4186.0/(6.02214*10**23)
fs2s = 10**-15
ang2m = 10**-10  

prefactor = 1/(Kb * Temperature**2 * volume)


scale =  kCal2J**2/( fs2s*ang2m )

prefactor = scale * prefactor

thermal_cond = kappa.transport_coeff(prefactor, 0, kappa.times[-1]) # In W/m K



print ("The thermal conductivity is %2.4e"%thermal_cond)

integral_lammps = cf.integrate(lammps_df["TimeDelta"]*time_step,lammps_df["c_flux[1]*c_flux[1]"],0,20000)
#
cond_lammps = integral_lammps



# =============================================================================
# # Plot kappa vs tau
# =============================================================================

fig,ax = plt.subplots()

tau_array = np.arange(1,kappa.times[-1],100)

kappa_array = []


for t in tau_array:
    # We compute the transport coefficients from the total (i.e with the average of the 3 components)
    kappa_array.append(kappa.transport_coeff(prefactor, 0, t))
    
    
ax.plot(tau_array,kappa_array)

ax.set_xlabel(r'$\tau[fs]$')
ax.axhline(y = thermal_cond, xmin=0, xmax=1,ls='--',c='black', label = "%2.4e" %thermal_cond)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_ylabel(r'$\kappa[W/(m\cdot K)]$')
plt.legend( loc = 'lower right')
plt.tight_layout()
plt.savefig("eta_vs_tau.pdf")
