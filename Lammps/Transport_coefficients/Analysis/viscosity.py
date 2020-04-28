#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:09:58 2020
Script to compute the viscosity based on the stress autocorrelation function

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

data = cf.read_data_file("pressures.dat")


time_step = 10 # In femptoseconds
data1 = data.values
times = (data1[:,0]-data1[0,0])*time_step



sampling_interval = times[1]-times[0] # Assuming times are homogeneous
max_tau = 20000


max_delta = int(max_tau/sampling_interval)  # int(len(times)*0.1) #Maximum delta of time to measure the correlation


pxy = data['v_pxy'].values

#pxy = np.reshape(np.)

pxy_flux = fc.flux(pxy,times,"Pxy")


etha11 = fc.correlation(pxy_flux,pxy_flux, max_delta)
etha11.evaluate()


print ("Using delta_t = %s fs" %time_step)


# =============================================================================
# Ploting all the correlations
# =============================================================================


lammps_df = cf.read_data_file("profile.gk.2d")
lammps_df = lammps_df.iloc[-400:]



plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()

fig,ax = etha11.plot_individual(fig, ax, norm = False)

#ax.set_xlim(1000,2000)
#ax.set_ylim(-0.0005,0.0005)
ax.plot(lammps_df["TimeDelta"]*time_step, lammps_df["v_pxy*v_pxy"],label ="Lammps")


ax.set_xscale('log')
ax.set_xlabel(r'$\tau[fs]$')

plt.legend(["Post-processing","On-the-fly"], loc = 'upper right')
ax.set_ylabel(r'$\langle %s(\tau) %s(0)\rangle$'%(etha11.flux1.name,etha11.flux2.name))
plt.tight_layout()
plt.savefig("correlation11.pdf")



# Units

box_side = 400.57


Kb = 1.380649*10**-23 #J/K
Temperature = 273.15 #K 
volume = box_side**3 # Angs**3 

# converting into SI

atm2pa = 101325.0
fs2s = 10**-15
ang2m = 10**-10  

prefactor = volume/(Kb*Temperature)

scale = ang2m**3 * atm2pa**2 * fs2s

prefactor = scale* prefactor

integral = etha11.transport_coeff(1, 0, etha11.times[-1]) # In atmospheres**2*fs


eta = integral * prefactor # in Pa s


print ("The viscosity is %2.4e"%eta)

integral_lammps = cf.integrate(lammps_df["TimeDelta"]*time_step,lammps_df["v_pxy*v_pxy"],0,20000)

eta_lammps = integral_lammps * prefactor

