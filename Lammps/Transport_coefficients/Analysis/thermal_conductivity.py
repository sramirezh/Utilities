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

data = cf.read_data_file("fluxes.dat")




# =============================================================================
# # Input parameters
# =============================================================================

time_step = 10 # In femptoseconds
box_side = 400.57
Temperature = 273.15 #K 
print ("\nUsing a box side of %s"%box_side)
print ("Using delta_t = %s fs" %time_step)
print ("Temperature = %s fs" %Temperature)




volume = box_side**3 # Angs**3 








data1 = data.values
times = (data1[:,0]-data1[0,0])*time_step



sampling_interval = times[1]-times[0] # Assuming times are homogeneous
max_tau = 20000


max_delta = int(max_tau/sampling_interval)  # int(len(times)*0.1) #Maximum delta of time to measure the correlation


Jx = data['v_Jx'].values

#pxy = np.reshape(np.)

Jx = fc.flux(Jx,times,"J_x")


etha11 = fc.correlation(Jx, Jx, max_delta)
etha11.evaluate()





# =============================================================================
# Ploting all the correlations
# =============================================================================


lammps_df = cf.read_data_file("J0Jt.dat")
lammps_df = lammps_df.iloc[-400:]



plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()

fig,ax = etha11.plot_individual(fig, ax, norm = False)

#ax.set_xlim(1000,2000)
#ax.set_ylim(-0.0005,0.0005)
ax.plot(lammps_df["TimeDelta"]*time_step, lammps_df["c_flux[1]*c_flux[1]"]/volume**2,label ="Lammps")


ax.set_xscale('log')
ax.set_xlabel(r'$\tau[fs]$')

plt.legend(["Post-processing","On-the-fly"], loc = 'upper right')
ax.set_ylabel(r'$\langle %s(\tau) %s(0)\rangle$'%(etha11.flux1.name,etha11.flux2.name))
plt.tight_layout()
plt.savefig("correlation11.pdf")





# converting into SI
Kb = 1.380649*10**-23 #J/K
kCal2J = 4186.0/(6.02214*10**23)
fs2s = 10**-15
ang2m = 10**-10  

prefactor = volume/(Kb*Temperature**2)


scale =  kCal2J**2/( fs2s*ang2m )

prefactor = scale* prefactor

integral = etha11.transport_coeff(prefactor, 0, etha11.times[-1]) # In atmospheres**2*fs

thermal_cond = integral # in Pa s


print ("The thermal conductivity is %2.4e"%thermal_cond)

integral_lammps = cf.integrate(lammps_df["TimeDelta"]*time_step,lammps_df["c_flux[1]*c_flux[1]"],0,20000)
#
cond_lammps = integral_lammps

