#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 20 13:42:54 2020
Script to get the theoretical estimation of the diffusio-osmotic velocity
Needs to be run inside Benchmark/3.Measurements
@author: simon
"""

import sys
import os
import matplotlib.pyplot as plt
import density_analysis as da
from copy import deepcopy
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
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

# Reading all the final averaged chunks
solute = da.DensityDistribution("Sproperties_short.dat",'rBulk')
solvent = da.DensityDistribution("Fproperties_short.dat",'rBulk')
fluid = da.DensityDistribution("properties_short.dat",'rBulk')

# Loading the data from the simualtion
sim_log = lu.Simulation("log.lammps")
z_wall = sim_log.limits[0, 2]

# Define the type of simulation
sim = deepcopy(Yoshida)
sim.print_params(logger)

# =============================================================================
# Computing the analytic properties
# =============================================================================
solute.compute_all_properties(solute.lower_limit)
solvent.compute_all_properties(solvent.lower_limit)

# =============================================================================
#  LADM
# =============================================================================
logger.info("\n\n!!!!!!!!!!!!!!!!!!!!!!!Performing LADM!!!!!!!!!!!!\n\n\n\n")
rho_ladm = fluid.compute_ladm(1)
viscosity_array = fluid.rho_dist.copy()
viscosity_array = [eta_meyer(rho, 1) for rho in rho_ladm]
fluid.data_frame['eta_ladm'] = viscosity_array
sim.eta = viscosity_array

# =============================================================================
# Getting the velocities as in Anderson1984
# =============================================================================
# For the solutes
grad_c_s = sim.grad_mu_s * solute.rho_bulk
v_0_s = - (sim.T)* solute.l*solute.k/Yoshida.eta
vx_s = v_0_s * grad_c_s

# For the solvents
grad_mu_f = - (solute.rho_bulk / solvent.rho_bulk)* sim.grad_mu_s
grad_c_f = grad_mu_f * solvent.rho_bulk
v_0_f = - (sim.T)* solvent.l*solvent.k/Yoshida.eta
vx_f = v_0_f * grad_c_f

vx_total = vx_s + vx_f

logger.info("The contribution to the velocity due to the solutes is %s"%vx_s)
logger.info("The contribution to the velocity due to the solvents is %s"%vx_f)
logger.info("The predicted DO velocity is  %s" % (vx_total)) 

# =============================================================================
# Getting the velocity distributions
# =============================================================================

# For the solutes assuming 3 different zeros

vx_z_0 = solute.vx_dist(sim, grad_c_s, 0.000000001) # wall is zero
vx_z_f = solute.vx_dist(sim, grad_c_s, solvent.lower_limit) #zero for solvents
solute.vx_dist(sim, grad_c_s, solute.lower_limit) # zero for solutes

# Plot comparison
cf.set_plot_appearance()
fig, ax = plt.subplots()

solute.plot_property_dist("vx_z", ax = ax)
ax.plot(solute.positions, vx_z_f)
ax.plot(solute.positions, vx_z_0)

ax.set_xlim(0, 20)
ax.set_ylim(0, None)
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$v_x(z)$')

fig.tight_layout()
ax.legend([r'$z_0 = %s$'%solute.lower_limit, r'$z_0 = %s$'%solvent.lower_limit, r'$z_0 = 0$'], loc = 'lower right')

fig.savefig('%s/vprofile_theo_solutes.pdf'%(plot_dir))

# For the solvents


vx_z_0 = solvent.vx_dist(sim, grad_c_f, 0.000000001) # wall is zero
vx_z_s = solvent.vx_dist(sim, grad_c_f, solute.lower_limit) #zero for solutes
solvent.vx_dist(sim, grad_c_f, solvent.lower_limit) # zero for solvents

# Plot comparison
fig, ax = plt.subplots()

ax.plot(solvent.positions, vx_z_s)
solvent.plot_property_dist("vx_z", ax = ax)
ax.plot(solvent.positions, vx_z_0)

ax.set_xlim(0, 20)
ax.set_ylim(0, None)
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$v_x(z)$')

fig.tight_layout()
ax.legend([r'$z_0 = %s$'%solute.lower_limit, r'$z_0 = %s$'%solvent.lower_limit, r'$z_0 = 0$'], loc = 'lower right')

fig.savefig('%s/vprofile_theo_solvents.pdf'%(plot_dir))


# =============================================================================
# # Plot total velocity
# =============================================================================
fig, ax = plt.subplots()
solute.plot_property_dist("vx_z", ax = ax)
solvent.plot_property_dist("vx_z", ax = ax)

total_velocity = solute.data_frame['vx_z'] + solvent.data_frame['vx_z']

ax.plot(solvent.positions, total_velocity )

ax.set_xlim(0, 20)
ax.set_ylim(0, None)
ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$v_x(z)$')

ax.legend(["Solute", "Solvent", "Fluid"], loc = 'upper right')

fig.tight_layout()
fig.savefig('%s/vprofile_total.pdf'%(plot_dir))

# =============================================================================
# Density, Gamma, Integrand K, Integrand L
# =============================================================================
fig,(ax1,ax2,ax3)=plt.subplots(3,1,sharex='col')

solute.plot_property_dist("density/mass", ax = ax1)
solute.plot_property_dist("rho_exc", ax = ax2)
solute.plot_property_dist("integrand_k", ax = ax3)

plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0.1)

fig.savefig("%s/distributions.pdf"%(plot_dir), transparent=True)

# =============================================================================
# Showing L integrand to ephasise the effect of z0
# =============================================================================

fig, ax = plt.subplots()

# Redefines the computation with the lower limit at solute's lower limit
solute.get_l(solute.lower_limit)
solute.plot_property_dist("integrand_first", ax = ax)

solute.get_l(solvent.lower_limit)
solute.plot_property_dist("integrand_first", ax = ax)

solute.get_l(0)
solute.plot_property_dist("integrand_first", ax = ax)


ax.set_xlabel(r'$z[\sigma] $')
ax.set_ylabel(r'$I_L(z)$')


ax.legend([r'$z_0 = %s$'%solute.lower_limit, r'$z_0 = %s$'%solvent.lower_limit, r'$z_0 = 0$'], loc = 'upper right')

ax = cf.plot_zoom(ax, [0,8])
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
fig.tight_layout()
fig.savefig('%s/first_moment_integrand.pdf'%(plot_dir))



# =============================================================================
# # Testing LADM
# =============================================================================

v_s = solute.vx_dist(sim, grad_c_s, solute.lower_limit) # zero for solutes

fig, ax = plt.subplots()
ax.plot(solute.positions, solute.data_frame['integrand_k'], label = 'integrand k', marker ='o', markersize=2)
ax.plot(fluid.positions, fluid.data_frame['ladm'], label = 'ladm', marker ='o', markersize=2)
ax.plot(fluid.positions, fluid.data_frame['density/mass'], label = 'density', marker ='o', markersize=2)
ax.plot(fluid.positions, fluid.data_frame['eta_ladm'], label = 'eta ladm', marker ='o', markersize=2)
ax.plot(solute.positions, solute.data_frame['vx_integrand'], label = 'integrand vx', marker ='o', markersize=2)
ax.plot(solute.positions, solute.data_frame['vx_z'], label = 'vx', marker ='o', markersize=2)

ax.set_xlabel(r'$z$')
ax = cf.plot_zoom(ax,[0,4])
ax.legend(loc = 'upper right')
ax.axhline(y = 0, ls='--',c='black')
ax.axvline(x = solute.lower_limit, ls='--',c='black')
fig.tight_layout()
fig.savefig('%s/test_integrand.pdf'%(plot_dir))



