#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 31 11:53:49 2020
Script to create a plot of all the results depending on the concentrations

The EMD analysis is base on analytic_results.py
@author: simon
"""


import glob
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.DO.EMD.density_analysis as da
import Lammps.lammps_utilities as lu


class SimulationEMD(lu.SimulationType):
    def __init__(self, name, T, eta, grad_mu_s):
        """
        name: simulation name, eg Yoshida.
        T: Reduced temperature in LJ units
        eta: viscosity of the bulk fluid
        grad_mu_s: gradient of chemical potential for the solutes
        """
        self.name = name
        self.T = T
        self.eta = eta
        self.grad_mu_s = grad_mu_s
        

class EMDResults(object):
    pass
        
# =============================================================================
# Main
# =============================================================================
logger = cf.log(__file__, os.getcwd())    

folders = glob.glob('x*')

# sorting the folders by the digit
a = cf.extract_digits(folders)
index_files = a[1]
folders = [folders[i] for i in index_files]

directory = "plots"
if not os.path.exists(directory):
    os.makedirs(directory)

# Removing the extremely low concentrations
folders = folders[2:]



# =============================================================================
# Preparing the distribution plots
# =============================================================================
cf.set_plot_appearance()

# Figure for the density distributions
fig1, ax1 = plt.subplots()
# Figure for the normalised density distributions
fig2, ax2 = plt.subplots()
# Figure for the integrand_k
fig3, ax3 = plt.subplots()

data_array = []
for folder in folders:
    x_0 = float(folder.split('_')[-1])
    path = '%s/3.Measuring/'%folder
    solute = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = path) 
    solute.compute_all_properties(solute.lower_limit)
    solvent = da.DensityDistribution("Fproperties_short.dat", "rBulk", directory = path) 
    solvent.compute_all_properties(solvent.lower_limit)
    
    # Getting properties from the bulk
    temp = solute.sim.thermo_ave.loc['Temp','Average']
    temp_error = solute.sim.thermo_ave.loc['Temp','Error_blocking']
    press = solute.sim.thermo_ave.loc['Pzz','Average']
    press_error = solute.sim.thermo_ave.loc['Pzz','Error_blocking']
    cbSolu = solute.sim.thermo_ave.loc['v_cBSolu','Average']
    cbSoluError = solute.sim.thermo_ave.loc['v_cBSolu','Error_blocking' ]
    
    data_array.append([x_0, temp, temp_error, press, press_error, solute.rho_bulk, solvent.rho_bulk,
                         solute.gamma, solute.k, solute.l, solute.k*solute.l,
                         solvent.gamma, solvent.k, solvent.l, solvent.k*solute.l,
                         cbSolu, cbSoluError ])
    
    
    # Plot how does the Number of solutes in the bulk, changes (To see if there is equilibration)
    fig, ax = plt.subplots()
    solute.sim.thermo_data.plot(x = 'Step', y = 'v_cBSolu', ax = ax, label = False)
    ax.set_xlabel("x label")
    ax.set_ylabel("y label")
    ax.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
    fig.savefig('%s/Nsolu_%s.pdf'%(directory, x_0))
    
    ax1.plot(solute.positions, solute.rho_dist, label = '%s'%x_0)
    ax2.plot(solute.positions, solute.rho_dist/max(solute.rho_dist), label = '%s'%x_0)
    ax3.plot(solute.positions, solute.data_frame['integrand_k'], label = '%s'%x_0)


# =============================================================================
# Distribution plots
# =============================================================================

ax1.set_ylabel(r'$c_s(z)$')
ax1.set_xlabel(r'$z$')
#ax1.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
ax1.legend(loc = 'upper right')
fig1.tight_layout()
fig1.savefig('density_distributions.pdf')


ax2.set_ylabel(r'$c_s(z)/c_s^{max}$')
ax2.set_xlabel(r'$z$')
#ax1.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
ax2.legend(loc = 'upper right')
fig2.tight_layout()
fig2.savefig('density_distributions_norm.pdf')



ax3.set_ylabel(r'$c_s(z)/c_s^B-1$')
ax3.set_xlabel(r'$z$')
#ax1.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
ax3.legend(loc = 'upper right')
fig3.tight_layout()
fig3.savefig('integrand_k.pdf')



# =============================================================================
# Preparing the data_file
# =============================================================================

columns = ['x0', 'temp', 'press', 'solute.rhobulk', 'solvent.rhobulk',
                         'solute.gamma', 'solute.k', 'solute.l','solute.kl',
                         'solvent.gamma', 'solvent.k', 'solvent.l', 'solvent.kl',
                         'cbSolu', 'cbSoluError']

df = pd.DataFrame(data_array, columns = columns)    
df = df.sort_values(['solute.rhobulk'], ascending = True)

# =============================================================================
# Plotting individual parameters
# =============================================================================

# Initial concentration in the box vs Number of solutes in the bulk after
# Equilibration
fig, ax = plt.subplots()
ax.errorbar(df['x0'], df['cbSolu'],yerr = df['cbSoluError'] )
ax.set_ylabel(r'$c_s^B$')
ax.set_xlabel(r'$x_0$')
ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('cbSolu.pdf')



# Initial concentration in the box vs equilibrated concentration in the bulk 
fig, ax = plt.subplots()
ax.scatter(df['x0'], df['solute.rhobulk'])
ax.plot(df['x0'], df['x0'], label = r'$x_0 = x_0$')
ax.set_ylabel(r'$c_s^i$')
ax.set_xlabel(r'$c_s^B$')
ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('solute_rho_bulk.pdf')


# Gamma vs concentration in the bulk
fig, ax = plt.subplots()
ax.scatter(df['solute.rhobulk'],df['solute.gamma'],label = 'Solutes')
ax.scatter(df['solute.rhobulk'],df['solvent.gamma'],label = 'Solvents')
ax.set_ylabel(r'$\Gamma $')
ax.set_xlabel(r'$c_s^B$')
ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('solute_gamma.pdf')

# L* vs concentration in the bulk
fig, ax = plt.subplots()
ax.scatter(df['solute.rhobulk'],df['solute.l'],label = 'Solutes')
ax.scatter(df['solute.rhobulk'],df['solvent.l'],label = 'Solvents')
ax.set_ylabel(r'$L^*$')
ax.set_xlabel(r'$c_s^B$')
ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('solute_l.pdf')


# K vs concentration in the bulk
fig, ax = plt.subplots()
ax.scatter(df['solute.rhobulk'],df['solute.k'],label = 'Solutes')
ax.scatter(df['solute.rhobulk'],df['solvent.k'],label = 'Solvents')
ax.set_xlabel(r'$c_s^B$')
ax.set_ylabel(r'$K$')
ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('solute_k.pdf')

# K L* vs concentration in the bulk
fig, ax = plt.subplots()
ax.scatter(df['solute.rhobulk'],df['solute.kl'],label = 'Solutes')
ax.scatter(df['solute.rhobulk'],df['solvent.kl'],label = 'Solvents')
ax.set_xlabel(r'$c_s^B$')
ax.set_ylabel(r'$KL^*$')
ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('solute_kl.pdf')



# For the solvents
#
#
#fig, ax = plt.subplots()
#df.plot(x='x0', y='solvent.rhobulk', kind = 'scatter', ax = ax)
#fig.tight_layout()
#fig.savefig('solvent_rho_bulk.pdf')
#
#fig, ax = plt.subplots()
#df.plot(x='solute.rhobulk',y = 'solvent.gamma', kind = 'scatter', ax = ax)
#fig.tight_layout()
#fig.savefig('solvent_gamma.pdf')
#
#fig, ax = plt.subplots()
#df.plot(x='solute.rhobulk',y = 'solvent.l', kind = 'scatter', ax = ax)
#fig.tight_layout()
#fig.savefig('solvent_l.pdf')
#
#fig, ax = plt.subplots()
#df.plot(x='solute.rhobulk',y = 'solvent.k', kind = 'scatter', ax = ax)
#fig.tight_layout()
#fig.savefig('solvent_k.pdf')
#
#fig, ax = plt.subplots()
#df.plot(x='solute.rhobulk',y = 'solvent.kl', kind = 'scatter', ax = ax)
#fig.tight_layout()
#fig.savefig('solvent_kl.pdf')

## Loading all the simulations 
#
## Write all the new types of simulations here
#Yoshida = SimulationEMD("Yoshida0125", 1, 1.57, -0.125)
#
## Define the type of simulation
#sim = Yoshida
#sim.print_params(logger)
#
#
#
##folder for the small bin simulations
#dir_f1= "4.Applying_force_0.125/1"
#dir_f2 = "4.Applying_force_0.063/1"
#dir_f3 = "4.Applying_force_0.025/1"
#dir_theo = "3.Measuring"
#
##logger.info("The data for the main data is from %s" %dir_normal_bin)
##logger.info("The data for the insert data is from %s" %dir_small_bin)
#
## Loading the NEMD data
#fluid_f1 = da.PropertyDistribution("properties_short.dat", directory = dir_f1) 
#fluid_f2 = da.PropertyDistribution("properties_short.dat", directory = dir_f2) 
#fluid_f3 = da.PropertyDistribution("properties_short.dat", directory = dir_f3) 
#
## Loading the EMD data
#solute = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = dir_theo) 
#solvent = da.DensityDistribution("Fproperties_short.dat", "rBulk", directory = dir_theo) 

### Loading the data for the bin 0.25 \sigma
##fluid_s = da.DensityDistribution("properties_short.dat", "rBulk", directory = dir_small_bin) 
#
## =============================================================================
## Computing the analytic properties
## =============================================================================
#solute.compute_all_properties(solute.lower_limit)
#solvent.compute_all_properties(solvent.lower_limit)
#
#
#grad_mu_s = [-0.125, -0.063, -0.025]
#velocity_t_dist = []
#velocity_s_dist = []
#
#for grad in grad_mu_s:
#    grad_c_s = grad * solute.rho_bulk
#    grad_mu_f = - (solute.rho_bulk / solvent.rho_bulk)* grad
#    grad_c_f = grad_mu_f * solvent.rho_bulk
#    
#    # Getting the contribution from both species
#    v_s = solute.vx_dist(sim, grad_c_s, solute.lower_limit) # zero for solutes
#    v_f =solvent.vx_dist(sim, grad_c_f, solvent.lower_limit) # zero for solvents
#    
#    v_total = v_s + v_f
#    velocity_s_dist.append(v_s)
#    velocity_t_dist.append(v_total)
#    
## =============================================================================
## Plotting the velocities
## =============================================================================
#cf.set_plot_appearance()
#
#plt.close("all")
#fig1, ax1 = plt.subplots()
#
#
#ax1.plot(fluid_f1.positions, fluid_f1.data_frame["vx"],marker = 'o', ls = '--', label = r'$\nabla \mu_s = -0.125$')
#ax1.plot(solute.positions,velocity_t_dist[0],c = ax1.lines[-1].get_color(), ls = '-.')
#ax1.plot(solute.positions,velocity_s_dist[0],c = ax1.lines[-1].get_color(), ls = ':')
#ax1.plot(fluid_f2.positions, fluid_f2.data_frame["vx"],marker = 'o', ls = '--', label = r'$\nabla \mu_s = -0.063$')
#ax1.plot(solute.positions,velocity_t_dist[1],c = ax1.lines[-1].get_color(), ls = '-.')
#ax1.plot(solute.positions,velocity_s_dist[1],c = ax1.lines[-1].get_color(), ls = ':')
#ax1.plot(fluid_f3.positions, fluid_f3.data_frame["vx"],marker = 'o', ls = '--', label = r'$\nabla \mu_s = -0.025$')
#ax1.plot(solute.positions,velocity_t_dist[2],c = ax1.lines[-1].get_color(), ls = '-.')
#ax1.plot(solute.positions,velocity_s_dist[2],c = ax1.lines[-1].get_color(), ls = ':')
#
#
##ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
#
#
#ax1.set_xlim(0, 30)
#ymin, ymax = ax1.get_ylim()
#ax1.set_ylim(0, 1.25* ymax)
#ax1.set_xlabel(r'$z[\sigma] $')
#ax1.set_ylabel(r'$v_x(z)$')
#ax1.legend(loc = 'upper right')
##
#### Adding the insert
##
##left, bottom, width, height = [0.55, 0.30, 0.4, 0.30]
##ax2 = fig1.add_axes([left, bottom, width, height])
##
##solute_s.plot_property_dist("vx", ax = ax2)
##solvent_s.plot_property_dist("vx", ax = ax2)
##fluid_s.plot_property_dist("vx", ax = ax2)
##
##
##xmin = 0
##xmax = 2
##ax2 = cf.plot_zoom(ax2, [xmin, xmax])
##ax2.tick_params(axis='both', which='major', labelsize = 12)
###ax2.set_ylabel(r'$V(r)$',fontsize =17, labelpad=-5)
###ax2.set_xlabel(r'$r$' ,fontsize =17, labelpad=-5)
##
##
#fig1.tight_layout()
#fig1.savefig('v_theo_sim.pdf')


