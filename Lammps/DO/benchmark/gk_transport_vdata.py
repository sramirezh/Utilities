#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 17:40:14 2020
Script to compute the viscosity based on the stress autocorrelation function, 
it uses 1 output from lammps:

    vdata.dat
    
and if new combinations of the fluxes are wanted, use:
    vdata.dat
    v_vx_Solv v_vx_Solu v_vx_Sol v_cSolu v_cSolv
    
Run inside 6.GJ_analysis

Define which folder pattern to use

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



class SimulationGK(lu.SimulationType):
    def __init__(self, name, time_step, max_tau, tau_integration, temp,
                 log_file, flux_file):
        """
        Args:
            name: Identifier of the simulation
            time_step: Simulation time_step
            max_tau: maximum length of the correlation analysis
            tau_integration: upper time limit of the correlation integration
            temp: system temperature
            log_file: lammps log file
            flux file: file containing all the fluxes
            
            
        Other atributes:
            sampling_interval
        """
        self.name = name
        self.time_step = time_step
        self.max_tau = max_tau
        self.tau_integration = tau_integration
        self.temp = temp
        self.log_file = log_file
        self.flux_file = flux_file
        self.flux_file_prefix = self.flux_file.split('.')[0]
        
    def get_thermo(self):
        """
        Gets the thermo data as an instance of lu.SImulation
        
        """
        
        self.thermo = lu.Simulation(self.log_file)
        

    def run_analysis(self):
        """
        This asssumes that the cwd is where all the files are
        TODO: this could be a method of the SimulationGK class
        Specific analysis where all the fluxes are their components are defined
        
        """
    
        data = cf.read_data_file(self.flux_file)
        
        data1 = data.values
        times = (data1[:,0]-data1[0,0]) * self.time_step
        
        
        # Assuming times are homogeneous
        self.sampling_interval = times[1]-times[0] 
    
        # int(len(times)*0.1) #Maximum delta of time to measure the correlation
        self.max_delta = int(self.max_tau/self.sampling_interval)  
        
        # Name of the components
        Q_components = ["v_Qx", "v_Qy"]
        if self.name == "stat":
            Js_components = ["v_Jsx_exc_cons", "v_Jsy_exc_cons"]
        else:
            Js_components = ["v_Jsx_exc", "v_Jsy_exc"]
        
        Q = data[Q_components].values
        Js = data[Js_components].values
        
        Q_flux = fc.flux(Q ,times, "Q")
        Js_flux = fc.flux(Js ,times, "J_s'")
        
        #  Creating the correlation instance
        m_qs = fc.correlation(Q_flux, Js_flux, self.max_delta)
        m_qs.evaluate()
        m_qs.save("Mqs_%s"%sim.save_name)
        
        m_sq = fc.correlation( Js_flux,Q_flux, self.max_delta)
        m_sq.evaluate()
        m_sq.save("Msq_%s"%sim.save_name)
        
        return m_qs, m_sq
    
    @property
    def save_name(self):
        """
        Returns the sufix to save the output files, if the flux file has a number
        adds this to the original name of the simulation
        """
        self.number = cf.extract_digits(self.flux_file)
        
        # If there are no numbers, the name does not change
        if len(self.number) == 0:
            self.number =''
        else:
            self.number = int(self.number[0])
        
        save_name =  '%s%s'%(self.name,self.number)
        return save_name


def remove_pkl(root):
    
    """
    Removes all the pkl inside the directory recursively
    TODO generalise and put in cf
    """
    files = glob.glob('%s/*.pkl'%root,  recursive = True)
    for file in files:
        os.remove(file)



# =============================================================================
# Main
# =============================================================================
plot_dir = "plots/0.gk_transport"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   


# =============================================================================
# Simulation type definitions
# =============================================================================

# High sampling frequency
high = SimulationGK("high", 0.005, 1000, 150, 1, "log.lammps", "fluxes_high.dat") 
# Low sampling frequency but averaging
low = SimulationGK("low", 0.005, 1000, 150, 1, "log.lammps", "fluxes_low.dat") 
# Low sampling frequency but averaging
other = SimulationGK("other", 0.005, 1000, 150, 1, "log.lammps", "fluxes_other.dat") 
# Low sampling frequency but averaging with static prefactor
stat = SimulationGK("stat", 0.005, 1000, 150, 1, "log.lammps", "fluxes_stat_low.dat") 


#sim_array = [high, low, other]

# Define the type of simulation
sim = other.copy
sim.print_params(logger)

# Getting the sampling times
#spd = lu.read_value_from("log.lammps", "fluxes.dat")

folder_pattern = 'r*'

logger.info("\nUsing the folder pattern %s"%folder_pattern )

folders = glob.glob(folder_pattern)

logger.info("Analysing the results in %s"%folders )

# Array to store all the correlations
correlation_dist = []
correlation_dist1 = []
transport_coeff = []
transport_coeff1 = []
cwd = os.getcwd()



# Preparing the plots

plt.close('all')
cf.set_plot_appearance()

#For m_qs
fig, ax = plt.subplots()
#For m_sq
fig1, ax1 = plt.subplots()

for i, folder in enumerate(folders):
    Kb = 1
    logger.info("Analysing the folder %s"%folder)
    # Computing this only once, assuming all runs have the same
    os.chdir(folder)
    partial_runs = glob.glob("%s*"%sim.flux_file_prefix)
    for j, file in enumerate(partial_runs):
        if i == 0 and j == 0:            
            # Heigth of the volume where the fluxes are measured
            sim.get_thermo()
            lz_system = sim.thermo.limits[1,2] - sim.thermo.limits[0,2]
            lz_box = lu.read_value_from("in.geom", 'rSystem')[-1]
            
            # We need the volume where the fluxes were measured
            volume_eff = sim.thermo.volume * lz_box/lz_system
            
            prefactor = volume_eff / (Kb * sim.temp)
        
        # The name of the 
        sim.flux_file = file
        logger.info("\nAnalysing the file %s"%file)
        # =============================================================================
        # Computing the viscosity from the stress data
        # =============================================================================
        
        # Loading the correlations if they were already computed, this saves almost 50 minutes in 16 cores (dexter)
        if (len(glob.glob("Mqs_%s.pkl"%sim.save_name)) == 1):
            
            logger.info("\nThere is a pkl, we need to load Mqs_%s.pkl\n"%sim.save_name)
            m_qs = cf.load_instance("Mqs_%s.pkl"%sim.save_name)
            m_sq = cf.load_instance("Msq_%s.pkl"%sim.save_name)
            
        else:
            logger.info("There is no pkl file, so we need to analyse")
            m_qs, m_sq = sim.run_analysis()
        
        transport_coeff.append( prefactor * m_qs.transport_coeff(0, sim.tau_integration) )
        transport_coeff1.append( prefactor * m_sq.transport_coeff(0, sim.tau_integration) )
        
        # Getting all the correlation distributions in x and y
        correlation_dist.extend([m_qs.cor[0], m_qs.cor[1]])
        correlation_dist1.extend([m_sq.cor[0], m_sq.cor[1]])
        
        fig,ax = m_qs.plot_individual(fig, ax, dim = 0, norm = False)
        ax.lines[-1].set_label(r"$%s_x$"%sim.save_name)
        fig,ax = m_qs.plot_individual(fig, ax, dim = 1, norm = False)
        ax.lines[-1].set_label(r"$%s_y$"%sim.save_name)
        
        
        fig1,ax1 = m_sq.plot_individual(fig, ax, dim = 0, norm = False)
        ax1.lines[-1].set_label(r"$%s_x$"%sim.save_name)
        fig1,ax1 = m_sq.plot_individual(fig, ax, dim = 1, norm = False)
        ax1.lines[-1].set_label(r"$%s_y$"%sim.save_name)
        
    os.chdir(cwd)




# =============================================================================
# Ploting all the correlations only for the last etha
# =============================================================================

#ax.set_xlim(1000,2000)
#ax.set_ylim(-0.0005,0.0005)
#ax.plot(lammps_df["TimeDelta"]*time_step, lammps_df["v_pxy*v_pxy"],label =r"$P_{xy}^{Lammps}$")

ax.axvline(x = sim.tau_integration,ls='-.',c='black')
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xscale('log')
xmin,xmax = ax.get_xlim()
ax.set_xlim(xmin, sim.max_tau)
ax.set_xlabel(r'$t^*$')

plt.legend(loc = 'upper right', fontsize = 8, ncol = 4)
ax.set_ylabel(r'$\langle %s(t^*) %s(0)\rangle$'%(m_qs.flux1.name, m_qs.flux2.name))
plt.tight_layout()
plt.savefig("%s/c_qs_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))

ax1.axvline(x = sim.tau_integration,ls='-.',c='black')
ax1.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')

ax1.set_xscale('log')
xmin,xmax = ax1.get_xlim()
ax1.set_xlim(xmin, sim.max_tau)
ax1.set_xlabel(r'$t^*$')

plt.legend(loc = 'upper right', fontsize = 8, ncol = 4)
ax1.set_ylabel(r'$\langle %s(t^*) %s(0)\rangle$'%(m_sq.flux1.name, m_sq.flux2.name))
plt.tight_layout()
plt.savefig("%s/c_sq_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))
#

# Testing the different types of averaging


## Test for acf
#
#from  scipy.signal import correlate
#fig, ax = plt.subplots()
#
#
#Q_flux = fc.flux(m_qs.flux1.components[:,0],m_qs.flux1.times, "Q")
#m_qs = fc.correlation(Q_flux, Q_flux, sim.max_delta)
#
## Using my method
#m_qs.evaluate()
#m_qs.transport_coeff(0, sim.tau_integration)
#print("the transport using evaluate is %s"%m_qs.transport_coeff(0, sim.tau_integration))
#fig,ax = m_qs.plot_individual(fig, ax, dim = 0, norm = False)
#
## Using stat acf
#m_qs.evaluate_acf()
#m_qs.transport_coeff(0, sim.tau_integration)
#print("the transport using evaluate_acf is %s"%m_qs.transport_coeff(0, sim.tau_integration))
#fig,ax = m_qs.plot_individual(fig, ax, dim = 0, norm = False)
#
#
##m_qs.evaluate_ccf()
##m_qs.transport_coeff(0, sim.tau_integration)
##print("the transport using evaluate_ccf is %s"%m_qs.transport_coeff(0, sim.tau_integration))
#
## Using scipy
#correlation = correlate(m_qs.flux1.components[:,0],m_qs.flux1.components[:,0])
#amplitude = np.correlate(m_qs.flux1.components[:,0],m_qs.flux1.components[:,0])/len(m_qs.flux1.components[:,0])
#correlation = amplitude * correlation
#
#
#ax.axvline(x = sim.tau_integration,ls='-.',c='black')
#ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')
#
#ax.set_xscale('log')
#xmin,xmax = ax.get_xlim()
#ax.set_xlim(xmin, sim.max_tau)
#ax.set_xlabel(r'$t^*$')
#plt.tight_layout()
#plt.savefig("%s/test_acf_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))


# =============================================================================
# Some operations on the correlation distribution in the x-direction
# =============================================================================

correlation_dist = prefactor * np.array(correlation_dist)
corr_dist_ave = np.average(correlation_dist, axis = 0)
corr_dist_error = sem(correlation_dist, axis = 0)



# Now plotting
fig, ax = plt.subplots()

ax.plot(m_qs.times, corr_dist_ave)
ax.fill_between(m_qs.times, corr_dist_ave - corr_dist_error, corr_dist_ave + corr_dist_error, alpha=0.4)
ax.axvline(x = sim.tau_integration,ls='-.',c='black')
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xscale('log')
xmin,xmax = ax.get_xlim()
ax.set_xlim(xmin, sim.max_tau)
ax.set_xlabel(r'$t^*$')

#plt.legend([r'$xx$',r'$yy$',r'$zz$',"Total"], loc = 'upper right', fontsize = 12, ncol = 2)
ax.set_ylabel(r'$\langle %s(t^*) %s(0)\rangle$'%(m_qs.flux1.name,m_qs.flux2.name))
plt.tight_layout()
plt.savefig("%s/correlation_ave_qs_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))



correlation_dist1 = prefactor * np.array(correlation_dist1)
corr_dist_ave1 = np.average(correlation_dist1, axis = 0)
corr_dist_error1 = sem(correlation_dist1, axis = 0)



# Now plotting
fig,ax = plt.subplots()

ax.plot(m_sq.times, corr_dist_ave1)
ax.fill_between(m_sq.times, corr_dist_ave1 - corr_dist_error1, corr_dist_ave1 + corr_dist_error1, alpha=0.4)
ax.axvline(x = sim.tau_integration,ls='-.',c='black')
ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')

ax.set_xscale('log')
xmin,xmax = ax.get_xlim()
ax.set_xlim(xmin, sim.max_tau)
ax.set_xlabel(r'$t^*$')

#plt.legend([r'$xx$',r'$yy$',r'$zz$',"Total"], loc = 'upper right', fontsize = 12, ncol = 2)
ax.set_ylabel(r'$\langle %s(t^*) %s(0)\rangle$'%(m_sq.flux1.name,m_sq.flux2.name))
plt.tight_layout()
plt.savefig("%s/correlation_ave_sq_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))



# =============================================================================
# # Plot M_qs vs tau for the average correlation distribution in x
# =============================================================================

fig,ax = plt.subplots()

tau_array = np.linspace([1], m_qs.times[-1])
# Need to add the integration time
tau_array = np.sort(np.append(tau_array, sim.tau_integration)) 

m_qs_array = []


for t in tau_array:
    # Taking only the nominal value
    m_qs_inst = cf.integrate(m_qs.times, corr_dist_ave, 0, t,)
    m_qs_array.append(m_qs_inst)


m_qs_array = np.array(m_qs_array)
    
ax.plot(tau_array, m_qs_array)
#ax.fill_between(tau_array, eta_array - eta_error, eta_array + eta_error, alpha=0.4)
ax.set_xlabel(r'$t^*$')
#ax.axhline(y = eta_tau_int.n , xmin=0, xmax=1,ls='--',c='black', label = r'$\eta^* = %2.4f$' %eta_tau_int.n)
ax.axvline(x = sim.tau_integration,ls='-.',c='black')
#xmin,xmax = ax.get_xlim()
#ymin,ymax = ax.get_ylim()
##ax.set_xlim(0, etha11.times[-1])
##ax.set_ylim(0, ymax)
ax.set_ylabel(r'$M^{QS}$')
#plt.legend( loc = 'lower right')
plt.tight_layout()
plt.savefig("%s/Mqs_vs_tau_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))


# =============================================================================
# # Plot M_qs vs tau for all
# =============================================================================

fig,ax = plt.subplots()

tau_array = np.linspace([1], m_qs.times[-1], 20)
# Need to add the integration time
tau_array = np.sort(np.append(tau_array, sim.tau_integration)) 

m_qs_array = []
m_qs_error = []

for t in tau_array:
    m_qs_t = []
    for cor in correlation_dist:
        # Taking only the nominal value
        m_qs_t.append(cf.integrate(m_qs.times, cor, 0, t,))
    
    m_qs_array.append(np.average(m_qs_t))
    m_qs_error.append(sem(m_qs_t))


m_qs_array = np.array(m_qs_array)
m_qs_error = np.array(m_qs_error)
    
ax.errorbar(tau_array, m_qs_array, yerr = m_qs_error)
#ax.fill_between(tau_array, eta_array - eta_error, eta_array + eta_error, alpha=0.4)
ax.set_xlabel(r'$t^*$')
#ax.axhline(y = eta_tau_int.n , xmin=0, xmax=1,ls='--',c='black', label = r'$\eta^* = %2.4f$' %eta_tau_int.n)
ax.axvline(x = sim.tau_integration,ls='-.',c='black')
#xmin,xmax = ax.get_xlim()
#ymin,ymax = ax.get_ylim()
##ax.set_xlim(0, etha11.times[-1])
##ax.set_ylim(0, ymax)
ax.set_ylabel(r'$M^{QS}$')
#plt.legend( loc = 'lower right')
plt.tight_layout()
plt.savefig("%s/Mqs_vs_tau_error_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))





fig,ax = plt.subplots()

tau_array = np.linspace([1], m_qs.times[-1], 20)
# Need to add the integration time
tau_array = np.sort(np.append(tau_array, sim.tau_integration)) 

m_sq_array = []
m_sq_error = []

for t in tau_array:
    m_sq_t = []
    for cor in correlation_dist1:
        # Taking only the nominal value
        m_sq_t.append(cf.integrate(m_sq.times, cor, 0, t,))
    
    m_sq_array.append(np.average(m_sq_t))
    m_sq_error.append(sem(m_sq_t))


m_sq_array = np.array(m_sq_array)
m_sq_error = np.array(m_sq_error)
    
ax.errorbar(tau_array, m_sq_array, yerr = m_sq_error)
#ax.fill_between(tau_array, eta_array - eta_error, eta_array + eta_error, alpha=0.4)
ax.set_xlabel(r'$t^*$')
#ax.axhline(y = eta_tau_int.n , xmin=0, xmax=1,ls='--',c='black', label = r'$\eta^* = %2.4f$' %eta_tau_int.n)
ax.axvline(x = sim.tau_integration,ls='-.',c='black')
#xmin,xmax = ax.get_xlim()
#ymin,ymax = ax.get_ylim()
##ax.set_xlim(0, etha11.times[-1])
##ax.set_ylim(0, ymax)
ax.set_ylabel(r'$M^{SQ}$')
#plt.legend( loc = 'lower right')
plt.tight_layout()
plt.savefig("%s/Msq_vs_tau_error_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))
#
