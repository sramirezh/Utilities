#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:09:58 2020
Script to compute the thermal conductivity based on the sflux autocorrelation function

The flux is computed in lammps as (The results are in dexter Transport_coeff/1.Moltemplate/3.N2/T_273.15/LEE_KIM/Heat_cond_test)
    
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
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import glob
from statsmodels.tsa.stattools import acf
import Lammps.Pore.GK.flux_correlation as fc


cwd = os.getcwd() #current working directory

def run_analysis(sim):

    data = cf.read_data_file("fluxes.dat")
    
    data1 = data.values
    times = (data1[:,0]-data1[0,0])*sim.ts 
    
    sampling_interval = times[1]-times[0] # Assuming times are homogeneous
    
    max_delta = int(np.floor(sim.max_tau/sampling_interval))  # int(len(times)*0.1) #Maximum delta of time to measure the correlation
    
    components = ["v_Jx", "v_Jy", "v_Jz"]
    
    j_components = data[components].values
    
    #pxy = np.reshape(np.)
    
    j_heat = fc.flux(j_components ,times,"J")
    
    
    kappa = fc.correlation(j_heat, j_heat, max_delta)
    kappa.evaluate_acf() # Only because they are the same variables
    kappa.save("kappa")
    
    return kappa




class simulation(object):
    def __init__(self, name, ts, box_side, temp, max_tau, d):
        self.name = name
        self.ts = ts
        self.box_side = box_side
        self.temp = temp
        self.max_tau = max_tau
        self.d = d   # Lammps p*s sampling_interval
        
    def print_params(self):
        print ("\nUsing a box side of %s"%self.box_side)
        print ("Using delta_t = %s fs" %self.ts)
        print ("Using max tau = %s fs"%self.max_tau)
        print ("Temperature = %s fs" %self.temp)
        



# =============================================================================
# Main
# =============================================================================

example = simulation("argon", 4, 21.504, 70, 1990*4, 200 )
mine = simulation("N2", 10, 400.57, 273.15, 1000000, 400 )


# this is the only thing to define
sim = mine

print ("We are using the parameters from the %s simulation"%sim.name)
sim.print_params()

volume = sim.box_side**3 # Angs**3 

# Loading the correlations if they were already computed, this saves almost 50 minutes in 16 cores (dexter)
if (len(glob.glob('kappa.pkl')) == 1):
    
    print ("\nThere is a pkl, we need to load kappa!\n")
    kappa = cf.load_instance("kappa.pkl")
    
else:
    print ("There is no pkl file, so we need to analyse")
    kappa = run_analysis(sim)




# =============================================================================
# Getting the correlations with acf and error
# =============================================================================
acf_stat, condifence = acf(kappa.flux1.components[:,0],nlags = len(kappa.times)-1, alpha =.05, fft = True)
acf_stat = kappa.norm[0] * acf_stat
confidence = kappa.norm[0]*  condifence 
std_error = (acf_stat-confidence[:,0])/2 
# =============================================================================
# Ploting all the correlations
# =============================================================================


lammps_df = cf.read_data_file("J0Jt.dat")
lammps_df = lammps_df.iloc[-sim.d:]


plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()
#
fig,ax = kappa.plot_all(fig, ax, norm = False)


# Results from lammps
ax.plot(lammps_df["TimeDelta"]*sim.ts, lammps_df["c_flux[1]*c_flux[1]"],label ="JxLammps")
ax.plot(lammps_df["TimeDelta"]*sim.ts, lammps_df["c_flux[2]*c_flux[2]"],label ="JyLammps")
ax.plot(lammps_df["TimeDelta"]*sim.ts, lammps_df["c_flux[3]*c_flux[3]"],label ="JzLammps")

ax.plot(kappa.times, acf_stat,  label = r"$J_x^{acf}$",ls='--', c='black')
ax.fill_between(kappa.times,  acf_stat-std_error, acf_stat+std_error , alpha=0.4)

ax.plot()
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

prefactor = 1/(Kb * sim.temp**2 * volume)


scale =  kCal2J**2/( fs2s*ang2m )

prefactor = scale * prefactor

thermal_cond = prefactor * kappa.transport_coeff(1, 0, kappa.times[-1]) # In W/m K



print ("The thermal conductivity is %2.4e"%thermal_cond)

integral_lammps_x = cf.integrate(lammps_df["TimeDelta"]*sim.ts,lammps_df["c_flux[1]*c_flux[1]"],0, sim.max_tau)
integral_lammps_y = cf.integrate(lammps_df["TimeDelta"]*sim.ts,lammps_df["c_flux[2]*c_flux[2]"],0, sim.max_tau)
integral_lammps_z = cf.integrate(lammps_df["TimeDelta"]*sim.ts,lammps_df["c_flux[3]*c_flux[3]"],0, sim.max_tau)



cond_lammps = prefactor * (integral_lammps_x+ integral_lammps_y+ integral_lammps_z)/3

print ("The thermal conductivity on the flight is %2.4e"%cond_lammps)

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
plt.savefig("kappa_vs_tau.pdf")




# =============================================================================
# # Plot kappa vs tau (Lammps)
# =============================================================================

fig,ax = plt.subplots()

tau_array = np.arange(1,lammps_df["TimeDelta"].values[-1]*sim.ts,100)

correlations = (lammps_df["c_flux[1]*c_flux[1]"]+lammps_df["c_flux[2]*c_flux[2]"]+lammps_df["c_flux[3]*c_flux[3]"])/3


kappa_lammps = []
for  t in  tau_array:
    kappa_lammps.append(cf.integrate(lammps_df["TimeDelta"]*sim.ts, correlations, 0, t))
     

ax.plot(tau_array,kappa_lammps)

ax.set_xlabel(r'$\tau[fs]$')
ax.axhline(y = thermal_cond, xmin=0, xmax=1,ls='--',c='black', label = "%2.4e" %cond_lammps)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_ylabel(r'$\kappa[W/(m\cdot K)]$')
plt.legend( loc = 'lower right')
plt.tight_layout()
plt.savefig("kappa_vs_tau_lammps.pdf")





