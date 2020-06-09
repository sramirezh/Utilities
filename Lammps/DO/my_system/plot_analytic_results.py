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
plot_dir = "plots/0.concentration_dependence_grad_mu"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)    

folders = glob.glob('x*')

# sorting the folders by the digit
a = cf.extract_digits(folders)
index_files = a[1]
folders = [folders[i] for i in index_files]

folders = folders[3:]

# =============================================================================
# Preparing the distribution plots
# =============================================================================
cf.set_plot_appearance()

# Figure for the solute density distributions
fig1, ax1 = plt.subplots()
# Figure for the normalised density distributions
fig2, ax2 = plt.subplots()
# Figure for the integrand_k
fig3, ax3 = plt.subplots()
# Figure for the solvent density distributions
fig4, ax4 = plt.subplots()
# Figure for the solute excess distributions
fig5, ax5 = plt.subplots()
# Figure for the solvent excess distributions
fig6, ax6 = plt.subplots()
# Figure for the solution concentration
fig7, ax7 = plt.subplots()
# Figure for the solution concentration LADM only
fig8, ax8 = plt.subplots()
# Figure for the solution viscosity LADM only
fig9, ax9 = plt.subplots()




data_array = []
for folder in folders:
    x_0 = float(folder.split('_')[-1])
    logger.info('\nAnalysing inside %s'%x_0)
    path = '%s/3.Measuring/'%folder
    
    solute = da.DensityDistribution("Sproperties_short.dat", "rBulk", directory = path) 
    solute.compute_all_properties(solute.lower_limit)
    
    solvent = da.DensityDistribution("Fproperties_short.dat", "rBulk", directory = path) 
    solvent.compute_all_properties(solvent.lower_limit)
    
    solution = da.DensityDistribution("properties_short.dat", "rBulk", directory = path) 
    
    # LADM
    logger.info("Performing LADM")
    rho_ladm = solution.compute_ladm(1)
    viscosity_array = solution.rho_dist.copy()
    viscosity_array = [eta_meyer(rho, 1) for rho in rho_ladm]

    
    # Getting properties from the bulk
    temp = solute.sim.thermo_ave.loc['Temp','Average']
    temp_error = solute.sim.thermo_ave.loc['Temp','Error_blocking']
    press = solute.sim.thermo_ave.loc['Pzz','Average']
    press_error = solute.sim.thermo_ave.loc['Pzz','Error_blocking']
    cbSolu = solute.sim.thermo_ave.loc['v_cBSolu','Average']
    cbSoluError = solute.sim.thermo_ave.loc['v_cBSolu','Error_blocking' ]
    
    data_array.append([x_0, temp, temp_error, press, press_error, solute.rho_bulk, solvent.rho_bulk,
                         solute.gamma, solute.k, solute.l, solute.k*solute.l,
                         solvent.gamma, solvent.k, solvent.l, solvent.k*solvent.l,
                         cbSolu, cbSoluError ])
    
    
    # Plot how does the Number of solutes in the bulk, changes (To see if there is equilibration)
    fig, ax = plt.subplots()
    ax.plot(solute.sim.thermo_data['Step'], solute.sim.thermo_data['v_cBSolu'])
#    solute.sim.thermo_data.plot(x = 'Step', y = 'v_cBSolu', ax = ax, label = False)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$N_s^B$")
    ax.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
    fig.tight_layout()
    fig.savefig('%s/Nsolu_%s.pdf'%(plot_dir, x_0))
    
    ax1.plot(solute.positions, solute.rho_dist, label = '%s'%x_0)
    ax2.plot(solute.positions, solute.rho_dist/max(solute.rho_dist), label = '%s'%x_0)
    ax3.plot(solute.positions, solute.data_frame['integrand_k'], label = '%s'%x_0)
    ax4.plot(solvent.positions, solvent.rho_dist, label = '%s'%x_0)
    ax5.plot(solute.positions, solute.rho_exc, marker ='o', markersize=2, label = r'$\Gamma = %1.2f, x_0 =%s$'%(solute.gamma,x_0))
    ax6.plot(solvent.positions, solvent.rho_exc, marker ='o', markersize=2, label = r'$\Gamma = %1.2f, x_0 =%s$'%(solvent.gamma,x_0))
    ax7.plot(solution.positions, solution.rho_dist, marker ='o', markersize=2,lw = 0.5,  label = '%s'%x_0)
    ax7.plot(solution.positions, solution.data_frame['ladm'],lw = 0.5, color = ax7.lines[-1].get_color())
    ax8.plot(solution.positions, solution.data_frame['ladm'], marker ='o', markersize=2,lw = 0.5,  label = '%s'%x_0)
    ax9.plot(solution.positions, viscosity_array, marker ='o', markersize=2,lw = 0.5,  label = '%s'%x_0 )

# =============================================================================
# Distribution plots
# =============================================================================

ax1.set_ylabel(r'$c_s(z)$')
ax1.set_xlabel(r'$z$')
ax1.set_xlim(0, 30)
ax1.set_ylim(0, None)
ax1.legend(loc = 'upper right')
fig1.tight_layout()
fig1.savefig('%s/solute_density_distributions.pdf'%(plot_dir))


ax2.set_ylabel(r'$c_s(z)/c_s^{max}$')
ax2.set_xlabel(r'$z$')
ax2.set_xlim(0, 30)
ax2.set_ylim(0, None)
#ax1.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
ax2.legend(loc = 'upper right')
fig2.tight_layout()
fig2.savefig('%s/density_distributions_norm.pdf'%(plot_dir))



ax3.set_ylabel(r'$c_s(z)/c_s^B-1$')
ax3.set_xlabel(r'$z$')
ax3.set_xlim(0, 30)
#ax1.axhline(y=solute.sim.thermo_ave.loc['v_cBSolu','Average'], xmin=0, xmax=1,ls='--',c='black')
ax3.legend(loc = 'upper right')
fig3.tight_layout()
fig3.savefig('%s/integrand_k.pdf'%(plot_dir))


ax4.set_ylabel(r'$c_f(z)$')
ax4.set_xlabel(r'$z$')
ax4.set_xlim(0, 30)
ax4.set_ylim(0, None)
ax4.legend(loc = 'upper right')
fig4.tight_layout()
fig4.savefig('%s/solvent_density_distributions.pdf'%(plot_dir))



ax5.set_ylabel(r'$c_s(z)-c_s^B$')
ax5.set_xlabel(r'$z$')
ax5.set_xlim(0, 30)
#ax5.set_ylim(0, None)
ax5.legend(loc = 'upper right')
ax5 = cf.plot_zoom(ax5,[0,8])
ax5.axvline(x=solute.lower_limit,ls='-.',c='black')
fig5.tight_layout()
fig5.savefig('%s/solute_excess.pdf'%(plot_dir))


ax6.set_ylabel(r'$c_f(z)-c_f^B$')
ax6.set_xlabel(r'$z$')
#ax6.set_xlim(0, 30)
ax6 = cf.plot_zoom(ax6,[0,8])
#ax5.set_ylim(0, None)
ax6.legend(loc = 'upper right')
ax6.axvline(x=solvent.lower_limit,ls='-.',c='black')
fig6.tight_layout()
fig6.savefig('%s/solvent_excess.pdf'%(plot_dir))


ax7.set_ylabel(r'$c(z)$')
ax7.set_xlabel(r'$z$')
ax7 = cf.plot_zoom(ax7,[0,8])
ax7.legend(loc = 'upper right')
ax7.axvline(solution.lower_limit,ls='--',c='black')
ax7.axvline(solution.lower_limit+1,ls='--',c='black')
fig7.tight_layout()
fig7.savefig('%s/solution_density_distributions.pdf'%(plot_dir))


ax8.set_ylabel(r'$c(z)$')
ax8.set_xlabel(r'$z$')
ax8 = cf.plot_zoom(ax8,[0,8])
ax8.legend(loc = 'upper right')
ax8.axvline(solution.lower_limit,ls='--',c='black')
ax8.axvline(solution.lower_limit+1,ls='--',c='black')
fig8.tight_layout()
fig8.savefig('%s/solution_density_distributions_LADM.pdf'%(plot_dir))


ax9.set_ylabel(r'$\eta(z)$')
ax9.set_xlabel(r'$z$')
ax9 = cf.plot_zoom(ax9,[0,8])
ax9.legend(loc = 'upper right')
ax9.axvline(solution.lower_limit,ls='--',c='black')
ax9.axvline(solution.lower_limit+1,ls='--',c='black')
fig9.tight_layout()
fig9.savefig('%s/eta_LADM.pdf'%(plot_dir))



# =============================================================================
# Preparing the data_file
# =============================================================================

columns = ['x0', 'temp', 'temp_error', 'press','press_error', 'solute.rhobulk', 'solvent.rhobulk',
                         'solute.gamma', 'solute.k', 'solute.l','solute.kl',
                         'solvent.gamma', 'solvent.k', 'solvent.l', 'solvent.kl',
                         'cbSolu', 'cbSoluError']

df = pd.DataFrame(data_array, columns = columns)    
df = df.sort_values(['solute.rhobulk'], ascending = True)

# =============================================================================
# Plotting individual parameters
# =============================================================================

# Temperatures
fig, ax = plt.subplots()
ax.errorbar(df['solute.rhobulk'], df['temp'],yerr = df['temp_error'] )
ax.set_xlabel(r'$c_s^B$')
ax.set_ylabel(r'$T$')
ax.axhline(y=1, xmin=0, xmax=1,ls='--',c='black')
fig.tight_layout()
fig.savefig('%s/temperatures.pdf'%(plot_dir))

# Pressures
fig, ax = plt.subplots()
ax.errorbar(df['solute.rhobulk'], df['press'],yerr = df['press_error'] )
ax.set_xlabel(r'$c_s^B$')
ax.set_ylabel(r'$P$')
ax.axhline(y=1, xmin=0, xmax=1,ls='--',c='black')
fig.tight_layout()
fig.savefig('%s/pressures.pdf'%(plot_dir))


# Initial concentration in the box vs Number of solutes in the bulk after
# Equilibration
fig, ax = plt.subplots()
ax.errorbar(df['x0'], df['cbSolu'],yerr = df['cbSoluError'] )
ax.set_ylabel(r'$c_s^B$')
ax.set_xlabel(r'$x_0$')
fig.tight_layout()
fig.savefig('%s/cbSolu.pdf'%(plot_dir))

# Initial concentration in the box vs equilibrated concentration in the bulk 
fig, ax = plt.subplots()
ax.scatter(df['x0'], df['solute.rhobulk'])
ax.plot(df['x0'], df['x0'], label = r'$x_0 = x_0$')
ax.set_ylabel(r'$c_s^i$')
ax.set_xlabel(r'$c_s^B$')
ax.legend(loc = 'upper right')
fig.tight_layout()
fig.savefig('%s/solute_rho_bulk.pdf'%(plot_dir))

# Gamma vs concentration in the bulk
fig, ax = plt.subplots()
ax.scatter(df['solute.rhobulk'],df['solute.gamma'],label = 'Solutes')
ax.scatter(df['solute.rhobulk'],df['solvent.gamma'],label = 'Solvents')
ax.set_ylabel(r'$\Gamma $')
ax.set_xlabel(r'$c_s^B$')
ax.legend(loc = 'upper right')
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
fig.tight_layout()
fig.savefig('%s/gamma.pdf'%(plot_dir))

# L* vs concentration in the bulk
fig, ax = plt.subplots()
ax.scatter(df['solute.rhobulk'],df['solute.l'],label = 'Solutes')
ax.scatter(df['solute.rhobulk'],df['solvent.l'],label = 'Solvents')
ax.set_ylabel(r'$L^*$')
ax.set_xlabel(r'$c_s^B$')
ax.legend(loc = 'upper right')
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
fig.tight_layout()
fig.savefig('%s/l.pdf'%(plot_dir))


# K vs concentration in the bulk
fig, ax = plt.subplots()
ax.scatter(df['solute.rhobulk'],df['solute.k'],label = 'Solutes')
ax.scatter(df['solute.rhobulk'],df['solvent.k'],label = 'Solvents')
ax.set_xlabel(r'$c_s^B$')
ax.set_ylabel(r'$K$')
ax.legend(loc = 'upper right')
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
fig.tight_layout()
fig.savefig('%s/k.pdf'%(plot_dir))

# K L* vs concentration in the bulk
fig, ax = plt.subplots()
ax.scatter(df['solute.rhobulk'],df['solute.kl'],label = 'Solutes')
ax.scatter(df['solute.rhobulk'],df['solvent.kl'],label = 'Solvents')
ax.set_xlabel(r'$c_s^B$')
ax.set_ylabel(r'$KL^*$')
ax.legend(loc = 'upper right')
ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')
fig.tight_layout()
fig.savefig('%s/kl.pdf'%(plot_dir))



## Trying LADM

#fig, ax = plt.subplots()
#ax.plot(solution.positions, solution.rho_dist, label = 'Measured')
#ax.plot(solution.positions, ladm, label = 'LADM')
#ax.set_ylabel(r'$c(z)$')
#ax.set_xlabel(r'$z$')
#ax.legend(loc = 'upper right')
#ax.set_ylim(0, None)
#ax.set_xlim(0, 8)
#fig.tight_layout()
#fig.savefig('%s/ladm.pdf'%(plot_dir))
#
#
#indexes = cf.get_interval(solution.positions, solution.lower_limit, 8)
#viscosity_array = [eta_meyer(rho, 1) for rho in ladm[indexes]]
#
#eta = 1.57
#
#fig, ax = plt.subplots()
#ax.plot(solution.positions[indexes], viscosity_array)
#ax.set_ylabel(r'$\eta(z)$')
#ax.set_xlabel(r'$z$')
#ax.set_ylim(0, None)
#ax.set_xlim(0, 8)
#fig.tight_layout()
#fig.savefig('%s/eta.pdf'%(plot_dir))
#
#
## Now plotting the viscosity




