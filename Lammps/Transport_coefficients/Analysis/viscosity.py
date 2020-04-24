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


data1 = data.values
times = (data1[:,0]-data1[0,0])
max_delta = int(len(times)*0.004) #Maximum delta of time to measure the correlation


pxy = data['v_pxy'].values

#pxy = np.reshape(np.)

pxy_flux = fc.flux(pxy,times,"Pxy")


etha11 = fc.correlation(pxy_flux,pxy_flux, max_delta)
etha11.evaluate()


delta_t = 1#10*10 # Sampling every 10 steps each with 10 fs
print ("Using delta_t = %s fs" %delta_t)


# =============================================================================
# Ploting all the correlations
# =============================================================================


#lammps_df = cf.read_data_file("profile.gk.2d")
#lammps_df = lammps_df.iloc[-400:]



plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()

fig,ax = etha11.plot_individual(fig, ax, norm = False)

#ax.set_xlim(1000,2000)
#ax.set_ylim(-0.0005,0.0005)
#ax.plot(lammps_df["TimeDelta"], lammps_df["v_pxy*v_pxy"],label ="Lammps")


#ax.set_xscale('log')
plt.legend()
plt.tight_layout()
plt.savefig("correlation11.pdf")



# Units

Kb = 1.380649*10**-23 #J/K
Temperature = 273.15 #K 
time_step = 10 # Fs
volume = 200.285**3 # Angs**3 

# converting into SI

atm2pa = 101325.0
fs2s = 10**-15
ang2m = 10**-10  



prefactor = volume*ang2m**3/(Kb*Temperature)*time_step 



integral = etha11.transport_coeff(1, 0, etha11.times[-1]) # In atmospheres**2*fs

integral = atm2pa**2 * fs2s * integral

eta = integral*prefactor # in Pa s






