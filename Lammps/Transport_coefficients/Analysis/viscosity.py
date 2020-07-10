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
from uncertainties import ufloat
import glob
import time
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
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
        
        pressures = fc.flux(p_off_diag ,times,"P")
        
        # Creating the correlation object
        eta = fc.correlation(pressures, pressures, self.max_delta)
        eta.evaluate_acf()
        eta.save("eta")
        
        return eta
    
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
        
from scipy.signal import savgol_filter

def ave_differences(array, averaging_cycles):
    """
    Smoothens a function by taking the average between each point and the nearest neightbors
    Args:
        array to smothen
        averaging_cycles: Number of averaging iteractions.
        
    Returns the averaged array
    """
    corr_ave = np.copy(array)
    
    for j in range(averaging_cycles):
        averaging = []
        for i,point in enumerate(corr_ave):
            
            # For the extremes
            if i == 0: 
                average = (2*point+corr_ave[i+1])/3
            elif i == len(corr_ave)-1:
                average = (2*point+corr_ave[i-1])/3
            # For the middle points
            else:
                average = (2*point+corr_ave[i+1]+corr_ave[i-1])/4
            
            averaging.append(average)
            
        corr_ave = averaging 
    
    return corr_ave



def savgol_n(data, order, window, poly):
    """
    Performs Apply a Savitzky-Golay filter to an array several times
    Args:
        data: data to smooth
        order: Number of successive smoothening steps
        window: number of points to average
        poly: order of the polynomial to perform a least squares fitting.
        
    retunrs:
        the smoothen array
    """
    smooth =np.copy(data)
    for i in range(order):
        smooth = savgol_filter(smooth, window, poly)
    
    return smooth
# =============================================================================
# Main
# =============================================================================

# Defining all the types of simulations

# High sampling frequency
N2 = SimulationViscosity("N2", 1000000,1000000, "log.lammps", "pressures.dat") 
Octane = SimulationViscosity("Octane", 300000, 150000, "log.lammps", "pressures.dat") 

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
if (len(glob.glob('eta.pkl')) == 1):
    
    logger.info ("\nThere is a pkl, we need to load eta!\n")
    t0 = time.time()
    eta = cf.load_instance("eta.pkl")
    logger.info ("It took %2.4f sec to load the file"%(time.time()-t0))
    
else:
    logger.info ("There is no pkl file, so we need to analyse")
    eta = sim.run_analysis()
    

 

# =============================================================================
# Ploting all the correlations
# =============================================================================
t0 = time.time()
# Loading lammps, on the fly calculations
lammps_df = cf.read_data_file("profile.gk.2d")
lammps_df = lammps_df.iloc[-400:]

plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()

eta.plot_all(ax, norm = False)

#ax.set_xlim(1000,2000)
#ax.set_ylim(-0.0005,0.0005)
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')
ax.set_xscale('log')
ax.set_xlabel(r'$\tau[fs]$')
ax.set_ylabel(r'$\langle %s(\tau) %s(0)\rangle$'%(eta.flux1.name,eta.flux2.name))
ax.legend( loc = 'upper right', fontsize = 12)
fig.tight_layout()
plt.savefig("%s/correlation.pdf"%sim.plot_dir)
ax.get_legend().remove()
ax.plot(lammps_df["TimeDelta"].values[::10] * sim.time_step, 
        lammps_df["v_pxy*v_pxy"].values[::10],label =r"$P_{xy}^{Lammps}$")
ax.legend( loc = 'upper right', fontsize = 12, ncol = 2)
fig.tight_layout()
plt.savefig("%s/correlation_lammps.pdf"%sim.plot_dir)



logger.info("It took %2.4f sec sec to create the correlation plots"%(time.time()-t0))

# =============================================================================
# # Plot eta vs tau
# =============================================================================

fig,ax = plt.subplots()

tau_array = np.linspace(eta.times[0], eta.times[-1], 20)
# Need to add the integration time
tau_array = np.sort(np.append(tau_array, sim.tau_integration)) 

eta_interval = []

# Computes the integral per intervals such that there is no need to recompute 
# Parts of the integral
#TODO be sure that the integration is continuos beween the intervals
initial_t = 0
for t in tau_array:
    eta_interval.append(sim.prefactor * eta.transport_coeff_comp(initial_t, t, -1))
    initial_t = t
    

# Neef to rewrite the first term, very quick fix
eta_interval[0] = ufloat(0,0)

eta_tau = np.cumsum(eta_interval)
eta_error = np.array([el.s for el in eta_tau])
eta_array = np.array([el.n for el in eta_tau])


# Getting the data explicitly for tau_integration
index = np.where(tau_array==sim.tau_integration)[0][0]
    
ax.plot(tau_array,eta_array)
ax.fill_between(tau_array, eta_array - eta_error, eta_array + eta_error, alpha=0.4)
ax.set_xlabel(r'$t[fs]$')
ax.axhline(y = eta_array[index] , xmin=0, xmax=1,ls='--',c='black', label = r'$\eta^* = %2.4f$' %eta_array[index])
ax.axvline(x = sim.tau_integration,ls='-.',c='black')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlim(0, eta.times[-1])
ax.set_ylim(0, ymax)
ax.set_ylabel(r'$\eta[Pa\cdot s]$')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
plt.tight_layout()
plt.savefig("%s/eta_vs_tau.pdf"%sim.plot_dir)



# Getting the 3 components for the integration cut-off
eta_tau = sim.prefactor * eta.transport_coeff(0,sim.tau_integration)

eta_tau_average = sum(eta_tau)/len(eta_tau)
logger.info("The 3 components of the viscosity are %s" %eta_tau)
logger.info("The average viscosity is %2.4e +/- %2.4e [Pa s]"%(eta_tau_average.n, eta_tau_average.s))

integral_lammps = cf.integrate(lammps_df["TimeDelta"].values*sim.time_step,lammps_df["v_pxy*v_pxy"].values,0,20000)

eta_lammps = integral_lammps * sim.prefactor
logger.info("On the fly estimation is %s [the tau cut-off is %s fs]" %(eta_lammps,lammps_df["TimeDelta"].values[-1]*sim.time_step))


# =============================================================================
# Smoothening data
# =============================================================================

correlation = [e.n for e in eta.cor[-1]]
times = eta.times

# First try
initial = 1000 
window = 100 

times_smooth = times[initial::window] # center of the intervals
time_intervals =np.arange(initial-int(window/2),times[-1],100, dtype = int) 

#y_filter = savgol_filter(correlation,51, 3) 

# Smoothing by averaging over the window 
smoothed_data = np.zeros(len(times_smooth))
for i in range(len(times_smooth)):
    smoothed_data[i] = np.average(correlation[time_intervals[i]:time_intervals[i+1]])

# Second technique
averaging_cycles = 100
corr_ave = ave_differences(correlation[initial:], averaging_cycles)

#savgol_filter
savgol_order = 100
savgol_window = 101
poly = 1
y_savgol = savgol_n(correlation, savgol_order, 101, poly)

# =============================================================================
# Plots smoothed data
# =============================================================================

fig,ax = plt.subplots()

ax.plot(times_smooth, smoothed_data, label='windows w = %s'%window )
#ax.plot(times[initial:], corr_ave, label='smooth %s'%averaging_cycles )
ax.plot(times[initial:], y_savgol[initial:], label='savgol int = %s window =%s poly =%s'%(savgol_order,savgol_window,poly))
#ax.plot(times[initial:], y_conv, label='convolve w = %s'%window )

xmin = 1000
xmax = times[-1]
cf.plot_zoom(ax, [xmin, xmax])
#ax.set_yscale('log')

ax.set_xlabel(r'$\tau[fs]$')
ax.set_ylabel(r'$\langle \sigma(\tau) \sigma(0)\rangle$')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')
ax.legend( loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig("%s/correlation_smooth.pdf"%sim.plot_dir)


# =============================================================================
# Computations with smoothed data
# =============================================================================

# Defining the points where the fitting is going to be 
initial_fitting = 5000
example_cut = 19000 # Cut for the testing on my laptop
real_cut = 100000 # cut for the whole dataset
cut = example_cut

# Numerical integration until this point

# point where the analytical integration starts
index_r_initial = np.where(tau_array>initial_fitting)[0][0]
r_initial_fitting = tau_array[index_r_initial]


# Get a fitting of the exponential tail
transformed = np.log(y_savgol[initial_fitting:cut])
polynomial = np.polyfit(times[initial_fitting:cut],transformed,1)

times_eval = np.arange(r_initial_fitting,cut)
fit_exp = np.exp(np.polyval(polynomial,times_eval))


#estimating the value of the analytical part of the integral

eta_tail = sim.prefactor * (np.exp(polynomial[1])/polynomial[0]*(np.exp(polynomial[0]*times_eval[-1])-np.exp(polynomial[0]*times_eval[0])))
eta_hybrid = eta_tail+eta_array[index_r_initial]

logger.info("The analytical tail is %s"%eta_tail)
logger.info("The estimation of the viscosity including the analytical tail is %s"%eta_hybrid)
# =============================================================================
# Plotting the smoothed tail
# =============================================================================

# First getting the indexes where the time is in the interest interval
indexes_tau = np.where((times>=initial) & (times <= cut)) # Cutting there to not have the log negative

fig,ax = plt.subplots()

ax.plot(times[indexes_tau], y_savgol[indexes_tau], label='savgol int = %s window =%s poly =%s'%(savgol_order,savgol_window,poly))
ax.plot(times_eval, fit_exp)

ax.set_yscale('log')
ax.set_xlim(initial, cut)
ax.set_xlabel(r'$\tau[fs]$')
ax.set_ylabel(r'$\langle \sigma(\tau) \sigma(0)\rangle$')
#ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
#ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')
#ax.legend( loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig("%s/correlation_smooth_savgol.pdf"%sim.plot_dir)