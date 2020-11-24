# -*- coding: utf-8 -*-
"""
Spyder Editor
Script to estimate the viscosity from lammps


"""

import numpy as np
import pandas as pd
import os
import sys
from copy import deepcopy
from io import StringIO
from scipy.stats import sem
from scipy import optimize
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.General.chunk_utilities as cu
import Lammps.core_functions as cf
import Others.Statistics.FastAverager as stat
import Lammps.DO.EMD.density_analysis as da
import Lammps.lammps_utilities as lu


def err_func(p, x, y, err):
    return (y - model(p, x)) / err

def model_disp(p,x):
    return p[0] + p[1] * (x-p[2])**2
    
def model(p,x):
    return p[0] + p[1] * x**2

def residual(p,x,y):
    return y - model(p,x)


# =============================================================================
# main
# =============================================================================
plot_dir = "plots/2.Poiseuille"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   


# =============================================================================
#  Problem parameters and simple calculations
# =============================================================================

grad_p_array = ['0.002', '0.001', '0.00063','0.00050', '0.00025' ]

# Preparing the plots
plt.close('all')
cf.set_plot_appearance()
fig1, ax1 = plt.subplots() # Plotting velocity distribution

viscosity_array = []
for grad_p in grad_p_array:
    
    dir_flow = "5.Applying_force_p_%s"%(grad_p)
    
    # if there is subfolder
    if os.path.exists("%s/1"%dir_flow):
        dir_flow = "%s/1"%dir_flow
        
    logger.info("\nWorking inside %s"%dir_flow)
    fluid = da.DensityDistribution("properties_short.dat",'rBulk', directory = dir_flow)
    
    force = 1/fluid.rho_bulk * float(grad_p)
    
    # Getting the bulk indexes
    indexes_bulk = cf.get_interval(fluid.positions, 10, 25)
    
    # Shifting all the positions, just for the fitting, so the zero is at the top
    # of the box
    fluid.data_frame['pos_shift'] = fluid.positions - fluid.positions[-1]
    
    
# =============================================================================
#     # Fitting
# =============================================================================
    pinit = [1.0, -1.0]
    out = optimize.leastsq(residual, pinit,
                       args=(fluid.data_frame['pos_shift'][indexes_bulk],
                             fluid.data_frame['vx'][indexes_bulk]), full_output=1)
    pfinal = out[0]
    viscosity = -((fluid.rho_bulk*force)/(2*pfinal[1]))
    
    viscosity_array.append(viscosity)
    
    vx_predicted = model(pfinal, fluid.data_frame['pos_shift'][indexes_bulk]) 
    
# =============================================================================
#     # Plotting
# =============================================================================
    ax1.plot(fluid.positions, fluid.data_frame['vx'], label = r'$F^e = %s$'%grad_p, ls = '--')
    ax1.plot(fluid.positions[indexes_bulk], vx_predicted, c = ax1.lines[-1].get_color(), ls = '-')
 
ax1.set_xlabel(r'$z[\sigma] $')
ax1.set_ylabel(r'$v_x(z)$')
ax1.legend(loc = 'upper left')
ax1.legend(loc = 'upper left')
ax1.set_xlim(0, 30)
ax1.set_ylim(0, None)
fig1.tight_layout()
fig1.savefig("%s/fitting.pdf"%plot_dir)     

eta_ave = np.average(viscosity_array)
eta_sem = sem(viscosity_array)

logger.info("The viscosity is %s +/- %s"%(eta_ave, eta_sem))

