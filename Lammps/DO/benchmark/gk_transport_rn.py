#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:40:14 2020
Script to compute the viscosity based on the stress autocorrelation function,
in the folders rn
 
it uses 1 output from lammps:

    fluxes.dat that contains TimeStep v_Qx v_Jsx_exc v_Qy v_Jsy_exc v_Qz v_Jsz_exc
    
and if new combinations of the fluxes are wanted, use:
    vdata.dat
    v_vx_Solv v_vx_Solu v_vx_Sol v_cSolu v_cSolv
    
Run inside 6.GJ_analysis
Based on viscosity from DO project
@author: simon
"""




import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.stats import sem
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.lammps_utilities as lu



import Lammps.Pore.GK.flux_correlation as fc



def run_analysis(flux_file, time_step, max_tau):
    """
    Specific analysis where all the fluxes are their components are defined
    
    """

    data = cf.read_data_file(flux_file)
    
    data1 = data.values
    times = (data1[:,0]-data1[0,0]) * time_step
    
    
    # Assuming times are homogeneous
    sampling_interval = times[1]-times[0] 

    # int(len(times)*0.1) #Maximum delta of time to measure the correlation
    max_delta = int(max_tau/sampling_interval)  
    
    # Name of the components
    Q_components = ["v_Qx", "v_Qy", "v_Qz"]
    Js_components = ["v_Jsx_exc", "v_Jsy_exc", "v_Jsz_exc" ]
    
    Q = data[Q_components].values
    Js = data[Js_components].values
    
    Q_flux = fc.flux(Q ,times, "Q")
    Js_flux = fc.flux(Js ,times, "J_s")
    
    
    etha11 = fc.correlation(Q_flux, Js_flux, max_delta)
    etha11.evaluate()
    etha11.save("etha11")
    
    return etha11

# =============================================================================
# Main
# =============================================================================
plot_dir = "plots"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   


time_step = 0.005 # In LJ time
max_tau = 1000 # in LJ time  Max tau to analyse
tau_integration = 125 # limit for the integration
Kb = 1  
Temperature = 1
log_file = "log.lammps"
flux_file = "fluxes.dat"


# Getting the sampling times
#spd = lu.read_value_from("log.lammps", "fluxes.dat")

logger.info("Using delta_t = %s " %time_step)
logger.info("Using max tau = %s for the analysis "%max_tau)
logger.info("Using integration tau = %s for the autocorrelation "%tau_integration)
logger.info("Getting the box size from %s" %log_file)



folders = glob.glob('r*')

# Array to store all the correlations
correlation_dist =[]
transport_coeff=[]
cwd = os.getcwd()

for i, folder in enumerate(folders):

    # Computing this only once, assuming all runs have the same
    os.chdir(folder)
    if i == 0:    
        volume, limits = lu.read_box_limits("log.lammps", -1)
        
        # TODO We need the volume where the fluxes were measured
#        volume_effective =
        prefactor = volume / (Kb * Temperature)

    
    # =============================================================================
    # Computing the viscosity from the stress data
    # =============================================================================
    
    # Loading the correlations if they were already computed, this saves almost 50 minutes in 16 cores (dexter)
    if (len(glob.glob('etha11*')) == 1):
        
        logger.info("\nThere is a pkl, we need to load etha11!\n")
        etha11 = cf.load_instance("etha11.pkl")
        
    else:
        logger.info("There is no pkl file, so we need to analyse")
        etha11 = run_analysis(flux_file, time_step, max_tau)
    
    
    transport_coeff.append( prefactor * etha11.transport_coeff(0, tau_integration) )
    
    # Getting all the correlation distributions only in x
    correlation_dist.append(etha11.cor[0])
    os.chdir(cwd)




# =============================================================================
# Ploting all the correlations only for the last etha
# =============================================================================

plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()

fig,ax = etha11.plot_all(fig, ax, norm = False)

#ax.set_xlim(1000,2000)
#ax.set_ylim(-0.0005,0.0005)
#ax.plot(lammps_df["TimeDelta"]*time_step, lammps_df["v_pxy*v_pxy"],label =r"$P_{xy}^{Lammps}$")

ax.axvline(x = tau_integration,ls='-.',c='black')
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xscale('log')
xmin,xmax = ax.get_xlim()
ax.set_xlim(xmin, max_tau)
ax.set_xlabel(r'$t^*$')

plt.legend([r'$xx$',r'$yy$',r'$zz$',"Total"], loc = 'upper right', fontsize = 12, ncol = 2)
ax.set_ylabel(r'$\langle %s(t^*) %s(0)\rangle$'%(etha11.flux1.name,etha11.flux2.name))
plt.tight_layout()
plt.savefig("%s/correlation11.pdf"%plot_dir)




# =============================================================================
# Some operations on the correlation distribution in the x-direction
# =============================================================================
correlation_dist = prefactor * np.array(correlation_dist)
corr_dist_ave = np.average(correlation_dist, axis = 0)
corr_dist_error = sem(correlation_dist, axis = 0)



# Now plotting
fig,ax = plt.subplots()

ax.plot(etha11.times, corr_dist_ave)
ax.fill_between(etha11.times, corr_dist_ave - corr_dist_error, corr_dist_ave + corr_dist_error, alpha=0.4)
ax.axvline(x = tau_integration,ls='-.',c='black')
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xscale('log')
xmin,xmax = ax.get_xlim()
ax.set_xlim(xmin, max_tau)
ax.set_xlabel(r'$t^*$')


for cor in correlation_dist:
    ax.plot(etha11.times, cor)

#plt.legend([r'$xx$',r'$yy$',r'$zz$',"Total"], loc = 'upper right', fontsize = 12, ncol = 2)
ax.set_ylabel(r'$\langle %s(t^*) %s(0)\rangle$'%(etha11.flux1.name,etha11.flux2.name))
plt.tight_layout()
plt.savefig("%s/correlation_ave_xx.pdf"%plot_dir)


# =============================================================================
# # Plot eta vs tau for the average correlation distribution in x
# =============================================================================

fig,ax = plt.subplots()

tau_array = np.linspace([1], etha11.times[-1])
# Need to add the integration time
tau_array = np.sort(np.append(tau_array, tau_integration)) 

eta_array = []
eta_error = []


for t in tau_array:
    # Taking only the nominal value
    eta = cf.integrate(etha11.times, corr_dist_ave, 0, t,)
    error = 0        
    eta_error.append(error)
    eta_array.append(eta)

eta_error = np.array(eta_error)
eta_array = np.array(eta_array)
    
ax.plot(tau_array,eta_array)
ax.fill_between(tau_array, eta_array - eta_error, eta_array + eta_error, alpha=0.4)
ax.set_xlabel(r'$t^*$')
#ax.axhline(y = eta_tau_int.n , xmin=0, xmax=1,ls='--',c='black', label = r'$\eta^* = %2.4f$' %eta_tau_int.n)
ax.axvline(x = tau_integration,ls='-.',c='black')
#xmin,xmax = ax.get_xlim()
#ymin,ymax = ax.get_ylim()
##ax.set_xlim(0, etha11.times[-1])
##ax.set_ylim(0, ymax)
ax.set_ylabel(r'$\eta^*$')
plt.legend( loc = 'lower right')
plt.tight_layout()
plt.savefig("%s/Mqs_vs_tau.pdf"%plot_dir)

