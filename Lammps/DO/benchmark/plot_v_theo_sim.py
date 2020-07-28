#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 11:53:49 2020
Script to create the plot of the velocity profile using NEMD and theory
Run inside Benchmark/x_0.2

The EMD analysis is base on analytic_results.py
@author: simon
"""

import os
import sys
import matplotlib.pyplot as plt
from copy import deepcopy
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
        
# =============================================================================
# Main
# =============================================================================
plot_dir = "plots/1.theo_sim_predictions"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   

# Write all the new types of simulations here
Yoshida = SimulationEMD("Yoshida0125", 1, 1.57, -0.125)

# Define the type of simulation
sim = deepcopy(Yoshida)
sim.print_params(logger)



#folder for the small bin simulations
dir_f1= "4.Applying_force_0.125/1"
dir_f2 = "4.Applying_force_0.063/1"
dir_f3 = "4.Applying_force_0.025/1"
dir_theo = "3.Measuring"
dir_benchmark= "4.Applying_force_hiroaki"

#logger.info("The data for the main data is from %s" %dir_normal_bin)
#logger.info("The data for the insert data is from %s" %dir_small_bin)

# Loading the NEMD data
fluid_f1 = da.PropertyDistribution("properties_short.dat", directory = dir_f1) 
fluid_f2 = da.PropertyDistribution("properties_short.dat", directory = dir_f2) 
fluid_f3 = da.PropertyDistribution("properties_short.dat", directory = dir_f3) 
fluid_f4 = da.PropertyDistribution("properties_short.dat", directory = dir_benchmark) 

# Loading the EMD data
solute = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = dir_theo) 
solvent = da.DensityDistribution("Fproperties_short.dat", "rBulk", directory = dir_theo) 
solution = da.DensityDistribution("properties_short.dat", "rBulk", directory = dir_theo) 

## Loading the data for the bin 0.25 \sigma
#fluid_s = da.DensityDistribution("properties_short.dat", "rBulk", directory = dir_small_bin) 


# =============================================================================
# #changing the viscosity to the average local
# =============================================================================

rho_ladm = solution.compute_ladm(1)

viscosity_array = solution.rho_dist.copy()
viscosity_array = [eta_meyer(rho, 1) for rho in rho_ladm]
sim.eta = viscosity_array

#average_density = solution.get_property_ave('density/mass',[solution.lower_limit, solution.limits_b[0]])
#sim.eta = eta_meyer(average_density, sim.T)
#
#logger.info("Changed the viscosity from %s to %s using Meyer et al for the average viscosity int he diffusive layer"%(Yoshida.eta,sim.eta))

# =============================================================================
# Computing the analytic properties
# =============================================================================
solute.compute_all_properties(solute.lower_limit)
solvent.compute_all_properties(solvent.lower_limit)

grad_mu_s = [-0.125, -0.063, -0.025]
velocity_t_dist = []
velocity_s_dist = []

for grad in grad_mu_s:
    grad_c_s = grad * solute.rho_bulk
    grad_mu_f = - (solute.rho_bulk / solvent.rho_bulk)* grad
    grad_c_f = grad_mu_f * solvent.rho_bulk
    
    # Getting the contribution from both species
    v_s = solute.vx_dist(sim, grad_c_s, solute.lower_limit) # zero for solutes
    v_f =solvent.vx_dist(sim, grad_c_f, solvent.lower_limit) # zero for solvents
    
    v_total = v_s + v_f
    velocity_s_dist.append(v_s)
    velocity_t_dist.append(v_total)
    
    
    
    
    
# =============================================================================
# Plotting 
# =============================================================================
cf.set_plot_appearance()
plt.close("all")
# =============================================================================
# Plotting Density distribution
# =============================================================================
fig1, ax1 = plt.subplots()

ax1.plot(solution.positions, solution.data_frame["density/mass"], label = "Equilibrium")
ax1.plot(fluid_f1.positions, fluid_f1.data_frame["density/mass"], label = r'$\nabla \mu_s = -0.125$')
ax1.plot(fluid_f2.positions, fluid_f2.data_frame["density/mass"], label = r'$\nabla \mu_s = -0.063$')
ax1.plot(fluid_f3.positions, fluid_f3.data_frame["density/mass"], label = r'$\nabla \mu_s = -0.025$')

ax1.set_xlim(0, 30)
ymin, ymax = ax1.get_ylim()
ax1.set_ylim(0, 1.25* ymax)
ax1.set_xlabel(r'$z[\sigma] $')
ax1.set_ylabel(r'$c(z)$')
ax1.legend(loc = 'upper right')

fig1.tight_layout()
fig1.savefig('%s/densities.pdf'%plot_dir)


# =============================================================================
# Plotting the velocities
# =============================================================================
cf.set_plot_appearance()


fig1, ax1 = plt.subplots()


ax1.plot(fluid_f1.positions, fluid_f1.data_frame["vx"],marker = 'o', ls = '--', label = r'$\nabla \mu_s = -0.125$')
ax1.plot(solute.positions,velocity_t_dist[0],c = ax1.lines[-1].get_color(), ls = '-.')
ax1.plot(solute.positions,velocity_s_dist[0],c = ax1.lines[-1].get_color(), ls = ':')
ax1.plot(fluid_f2.positions, fluid_f2.data_frame["vx"],marker = 'o', ls = '--', label = r'$\nabla \mu_s = -0.063$')
ax1.plot(solute.positions,velocity_t_dist[1],c = ax1.lines[-1].get_color(), ls = '-.')
ax1.plot(solute.positions,velocity_s_dist[1],c = ax1.lines[-1].get_color(), ls = ':')
ax1.plot(fluid_f3.positions, fluid_f3.data_frame["vx"],marker = 'o', ls = '--', label = r'$\nabla \mu_s = -0.025$')
ax1.plot(solute.positions,velocity_t_dist[2],c = ax1.lines[-1].get_color(), ls = '-.')
ax1.plot(solute.positions,velocity_s_dist[2],c = ax1.lines[-1].get_color(), ls = ':')


#ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')


ax1.set_xlim(0, 30)
ymin, ymax = ax1.get_ylim()
ax1.set_ylim(0, 1.25* ymax)
ax1.set_xlabel(r'$z[\sigma] $')
ax1.set_ylabel(r'$v_x(z)$')
ax1.legend(loc = 'upper right')
#
### Adding the insert
#
#left, bottom, width, height = [0.55, 0.30, 0.4, 0.30]
#ax2 = fig1.add_axes([left, bottom, width, height])
#
#solute_s.plot_property_dist("vx", ax = ax2)
#solvent_s.plot_property_dist("vx", ax = ax2)
#fluid_s.plot_property_dist("vx", ax = ax2)
#
#
#xmin = 0
#xmax = 2
#ax2 = cf.plot_zoom(ax2, [xmin, xmax])
#ax2.tick_params(axis='both', which='major', labelsize = 12)
##ax2.set_ylabel(r'$V(r)$',fontsize =17, labelpad=-5)
##ax2.set_xlabel(r'$r$' ,fontsize =17, labelpad=-5)
#
#
fig1.tight_layout()
fig1.savefig('%s/v_theo_sim.pdf'%plot_dir)

# Plotting With Hiroaki's results

fig1, ax1 = plt.subplots()


ax1.plot(fluid_f1.positions, fluid_f1.data_frame["vx"],marker = 'o', ls = '--', label = r'$F_s = -\nabla \mu_s = 0.125$')
ax1.plot(solute.positions,velocity_t_dist[0], ls = '-.', label = "Theory solutes + solvents")
ax1.plot(solute.positions,velocity_s_dist[0], ls = ':', label = "Theory only solutes")
ax1.plot(fluid_f4.positions, fluid_f4.data_frame["vx"],marker = 'o', ls = '--', label = r'$F_s = -\frac{N^B-N_s^B}{N^B}\nabla \mu_s$')
ax1.set_xlim(0, 30)
ymin, ymax = ax1.get_ylim()
ax1.set_ylim(0, 1.25* ymax)
ax1.set_xlabel(r'$z[\sigma] $')
ax1.set_ylabel(r'$v_x(z)$')
ax1.legend(loc = 'lower right')
fig1.tight_layout()
fig1.savefig('%s/theo_sim_bench.pdf'%plot_dir)


# =============================================================================
# Plotting LADM
# =============================================================================
logger.info("\n\n!!!!!!!!!!!!!!!!!!!!!!!Performing LADM!!!!!!!!!!!!\n\n\n\n")
# Densities
fig, ax = plt.subplots()
ax.plot(solution.positions, solution.rho_dist, label = 'Measured')
ax.plot(solution.positions, solution.data_frame['ladm'], label = 'LADM')
ax.set_ylabel(r'$c(z)$')
ax.set_xlabel(r'$z$')
ax.legend(loc = 'upper right')
ax.set_ylim(0, None)
ax.set_xlim(0, 8)
fig.tight_layout()
fig.savefig('%s/ladm.pdf'%(plot_dir))

# Viscosity
fig, ax = plt.subplots()
ax.plot(solution.positions, sim.eta)
ax.axhline(y=Yoshida.eta, xmin=0, xmax=1,ls='--',c='black')
ax.set_ylabel(r'$\eta(z)$')
ax.set_xlabel(r'$z$')
ax.set_ylim(0, None)
ax.set_xlim(0, 8)
fig.tight_layout()
fig.savefig('%s/eta.pdf'%(plot_dir))




## =============================================================================
## # Testing LADM
## =============================================================================
#
#v_s = solute.vx_dist(sim, grad_c_s, solute.lower_limit) # zero for solutes
#
#fig, ax = plt.subplots()
#ax.plot(solute.positions, solute.data_frame['integrand_k'], label = 'integrand k', marker ='o', markersize=2)
#ax.plot(solution.positions, solution.data_frame['ladm'], label = 'ladm', marker ='o', markersize=2)
#ax.plot(solution.positions, solution.data_frame['density/mass'], label = 'density', marker ='o', markersize=2)
#ax.plot(solution.positions, solution.data_frame['eta_ladm'], label = 'eta ladm', marker ='o', markersize=2)
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



