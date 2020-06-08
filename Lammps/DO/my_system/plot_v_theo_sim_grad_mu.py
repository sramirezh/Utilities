#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat June 1st 11:53:49 2020
velocity profile using NEMD and theory for simulations changing the concentration
at constan chemical potential gradient

The EMD analysis is base on analytic_results.py
@author: simon
"""


import glob
import os
import sys
import pandas as pd
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.DO.EMD.density_analysis as da
import Lammps.lammps_utilities as lu
from Lammps.DO.EMD.Meyer_2018 import eta_meyer

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
        

def plot_sim_theo(fluid, solute, solvent):
    """
    Plots the velocity distribution using theorerical predictions and simulations
    Args:
        fluid: is a da.PropertyDistribution instance
        solute: is a  da.DensityDistribution instance
        solvent: is a  da.DensityDistribution instance
    Note: for the solvent and solute instances, both the compute_all and vx_dist
    methods had to be run in advance
    Returns:
        ax,fig
    """
    fig1, ax1 = plt.subplots()
    
    v_total = solute.data_frame['vx_z']+solvent.data_frame['vx_z']
    ax1.plot(solute.positions, solute.data_frame['vx_z'], label = "Solutes")
    ax1.plot(solute.positions, solvent.data_frame['vx_z'],label = "Solvents")
    ax1.plot(solute.positions,v_total, label = "Total")
    ax1.plot(fluid.positions, fluid.data_frame["vx"],marker = 'o', ls = '--', label = "NEMD")

    
    ax1.set_xlim(0, 10)
    ymin, ymax = ax1.get_ylim()
    ax1.set_ylim(None, 1.25* ymax)
    ax1.set_xlabel(r'$z[\sigma] $')
    ax1.set_ylabel(r'$v_x(z)$')
    ax1.legend(loc = 'upper right', ncol =2)
    
    fig1.tight_layout()
    
    return fig1, ax1
        
# =============================================================================
# Main
# =============================================================================
plot_dir = "plots/1.theo_sim_predictions"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)    

mine = SimulationEMD("mine0125", 1, 1.57, -0.125)

# Define the type of simulation
sim = deepcopy(mine)
sim.print_params(logger)

# Getting all the folders with simulations applying forces inside 
folders = glob.glob("x*/4.Applying_force")
folders =  [f.split('/')[0] for f in folders]

# sorting the folders by the digit
a = cf.extract_digits(folders)
index_files = a[1]
folders = [folders[i] for i in index_files]


# Removing the extremely low concentrations
#folders = folders[2:]
#
#folders = ['x_0.05', 'x_0.3', 'x_0.4', 'x_0.5', 'x_0.9']
#folders = ['x_0.05','x_0.9']

# =============================================================================
# Preparing the distribution plots
# =============================================================================
plt.close('all')
cf.set_plot_appearance()

data_array = []
for folder in folders:
    
    logger.info(folder)
    x_0 = float(folder.split('_')[-1])
    path_theo = '%s/3.Measuring/'%folder
    path_nemd = '%s/4.Applying_force/'%folder
    
    fluid = da.DensityDistribution("properties_short.dat", "rBulk" , directory = path_nemd)
    fluid_theo = da.DensityDistribution("properties_short.dat", "rBulk" , directory = path_theo)
    solute = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = path_theo) 
    solvent = da.DensityDistribution("Fproperties_short.dat", "rBulk", directory = path_theo) 
    
    
#    # =============================================================================
#    # #changing the viscosity to the average local
#    # =============================================================================
#    logger.info("\n\n!!!!!!!!!!!!!!!!!!!!!!!Performing LADM!!!!!!!!!!!!\n\n\n\n")
#    rho_ladm = fluid_theo.compute_ladm(1)
#
#    viscosity_array = fluid_theo.rho_dist.copy()
#    viscosity_array = np.array([eta_meyer(rho, sim.T) for rho in rho_ladm])
##    viscosity_array[:]  = mine.eta
##    viscosity_array[6] = viscosity_array[5]
#    sim.eta = viscosity_array
#    fluid_theo.data_frame['eta_ladm'] = viscosity_array
    
#    logger.info("Changed the viscosity from %s to %s using Meyer et al for the average viscosity inint he diffusive layer"%(mine.eta,sim.eta))

    solute.compute_all_properties(solute.lower_limit)
    solvent.compute_all_properties(solvent.lower_limit)
    
    
    # Computing the chemical potential gradients for each case
    grad_c_s = sim.grad_mu_s * solute.rho_bulk
    grad_mu_f = - (solute.rho_bulk / solvent.rho_bulk)* sim.grad_mu_s
    logger.info("The chemical potential gradient for the solvent is  %s for x_0 = %s"%(grad_mu_f, x_0))
    grad_c_f = grad_mu_f * solvent.rho_bulk
    
    # Computing the theoretical velocities
    v_s = solute.vx_dist(sim, grad_c_s, solute.lower_limit) # zero for solutes
    v_f = solvent.vx_dist(sim, grad_c_f, solvent.lower_limit) # zero for solvents
    
    # Getting the average vx in the bulk region from simulations
    vx_sim = fluid.get_bulk_property('vx')
    # Plotting
    
    fig1, ax1 = plot_sim_theo(fluid, solute, solvent)
    fig1.savefig('%s/v_theo_sim_%s.pdf'%(plot_dir,x_0))
    

    
#    v_total = v_s + v_f
#    velocity_s_dist.append(v_s)
#    velocity_t_dist.append(v_total)
##    
    data_array.append([x_0, solute.rho_bulk, solvent.rho_bulk,
                       grad_c_s, grad_mu_f,grad_c_f, solute.gamma, solute.k, solute.l, solute.k*solute.l,
                         solvent.gamma, solvent.k, solvent.l, solvent.k*solute.l,
                         solute.vx_bulk, solvent.vx_bulk, vx_sim])
#    
#    
#    # Plot how does the Number of solutes in the bulk, changes (To see if there is equilibration)
#    fig, ax = plt.subplots()
#    solute.sim.thermo_data.plot(x = 'Step', y = 'v_cBSolu', ax = ax, label = False)
#    ax.set_xlabel("x label")
#    ax.set_ylabel("y label")
#    ax.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
#    fig.savefig('%s/Nsolu_%s.pdf'%(directory, x_0))
#    
#    ax1.plot(solute.positions, solute.rho_dist, label = '%s'%x_0)
#    ax2.plot(solute.positions, solute.rho_dist/max(solute.rho_dist), label = '%s'%x_0)
#    ax3.plot(solute.positions, solute.data_frame['integrand_k'], label = '%s'%x_0)


# =============================================================================
# Distribution plots
# =============================================================================
#
#ax1.set_ylabel(r'$c_s(z)$')
#ax1.set_xlabel(r'$z$')
#
#ax1.set_xlim(0, 30)
#ax1.set_ylim(0, None)
#ax1.legend(loc = 'upper right')
#fig1.tight_layout()
#fig1.savefig('density_distributions.pdf')
#
#
#ax2.set_ylabel(r'$c_s(z)/c_s^{max}$')
#ax2.set_xlabel(r'$z$')
#ax2.set_xlim(0, 30)
#ax2.set_ylim(0, None)
##ax1.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
#ax2.legend(loc = 'upper right')
#fig2.tight_layout()
#fig2.savefig('density_distributions_norm.pdf')
#
#
#
#ax3.set_ylabel(r'$c_s(z)/c_s^B-1$')
#ax3.set_xlabel(r'$z$')
#ax3.set_xlim(0, 30)
#
##ax1.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
#ax3.legend(loc = 'upper right')
#fig3.tight_layout()
#fig3.savefig('integrand_k.pdf')
#
#
#
# =============================================================================
# Preparing the data_file
# =============================================================================

columns = ['x_0', 'solute.rho_bulk', 'solvent.rho_bulk',
                       'grad_c_s', 'grad_mu_f','grad_c_f','solute.gamma', 'solute.k', 'solute.l','solute.kl',
                         'solvent.gamma', 'solvent.k', 'solvent.l', 'solvent.kl',
                         'solute.vx_bulk', 'solvent.vx_bulk','vx_sim']

df = pd.DataFrame(data_array, columns = columns)    
df = df.sort_values(['solute.rho_bulk'], ascending = True)
#
## =============================================================================
## Plotting individual parameters
## =============================================================================


# Bulk velocity contribution vs x_0
fig, ax = plt.subplots()
ax.scatter(df['x_0'], df['solute.vx_bulk'], label = "Solute")
ax.scatter(df['x_0'], df['solvent.vx_bulk'], label = "Solvent")
ax.scatter(df['x_0'], df['solvent.vx_bulk']+df['solute.vx_bulk'], label = "Fluid")
ax.scatter(df['x_0'], df['vx_sim'], label = "NEMD")
ax.set_ylabel(r'$v_x^B$')
ax.set_xlabel(r'$x_0$')
ymin, ymax = ax.get_ylim()
ax.set_ylim(0, 1.25 * ymax )
ax.set_xlim(0, None)
ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('%s/vel_bulk_vs_x0.pdf'%(plot_dir))

# Bulk velocity contribution vs cs^B
fig, ax = plt.subplots()
ax.scatter(df['solute.rho_bulk']/df['solvent.rho_bulk'], df['solute.vx_bulk'], label = "Solute")
ax.scatter(df['solute.rho_bulk']/df['solvent.rho_bulk'], df['solvent.vx_bulk'], label = "Solvent")
ax.scatter(df['solute.rho_bulk']/df['solvent.rho_bulk'], df['solvent.vx_bulk']+df['solute.vx_bulk'], label = "Fluid")
ax.scatter(df['solute.rho_bulk']/df['solvent.rho_bulk'], df['vx_sim'], label = "NEMD")
ax.set_ylabel(r'$v_x^B$')
ax.set_xlabel(r'$c_s^B/c_f^B$')
ymin, ymax = ax.get_ylim()
ax.set_ylim(0, 1.25 * ymax )
#ax.set_xlim(0, None)
ax.set_xscale('log')
ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('%s/vel_bulk_vs_c_ratio.pdf'%(plot_dir))


# Chemical potential solvent contribution vs c^B ration
fig, ax = plt.subplots()
ax.scatter(df['solute.rho_bulk']/df['solvent.rho_bulk'], df['grad_mu_f'])
ax.set_ylabel(r'$\nabla \mu_f$')
ax.set_xlabel(r'$c_s^B/c_f^B$')
ymin, ymax = ax.get_ylim()
#ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('%s/grad_mu_f.pdf'%(plot_dir))




# Testing LADM

#df = solute.data_frame
#df['eta'] = fluid_theo.data_frame['eta_ladm']
#df['rho_fluid'] = fluid_theo.data_frame['density/mass']
#df['rho_fluid_lamd'] = fluid_theo.data_frame['ladm']
#
#
#
#v_s = solute.vx_dist(sim, grad_c_s, solute.lower_limit) # zero for solutes
#
#fig, ax = plt.subplots()
#ax.plot(solute.positions, solute.data_frame['integrand_k'], label = 'integrand k', marker ='o', markersize=2)
#ax.plot(fluid_theo.positions, fluid_theo.data_frame['ladm'], label = 'ladm', marker ='o', markersize=2)
#ax.plot(fluid_theo.positions, fluid_theo.data_frame['density/mass'], label = 'density', marker ='o', markersize=2)
#ax.plot(fluid_theo.positions, fluid_theo.data_frame['eta_ladm'], label = 'eta ladm', marker ='o', markersize=2)
#ax.plot(solute.positions, solute.data_frame['vx_integrand'], label = 'integrand vx', marker ='o', markersize=2)
#ax.plot(solute.positions, solute.data_frame['vx_z'], label = 'vx', marker ='o', markersize=2)
#
#ax.set_xlabel(r'$z$')
#ax = cf.plot_zoom(ax,[0,4])
#ax.legend(loc = 'upper right')
#ax.axhline(y = 0, ls='--',c='black')
#ax.axvline(x = solute.lower_limit, ls='--',c='black')
#fig.tight_layout()
#fig.savefig('%s/test_integrand.pdf'%(plot_dir))

