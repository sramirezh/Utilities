#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 18:03:55 2020
This set of functions computes the diffusion coefficient
@author: simon
"""

import sys
import os
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf


fitfunc = lambda p, x: p[0] * x + p[1] #Fitting to a line
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / (err+10**-8)

def fit_line(x,y,yerr, initial_index = 50):
    """
    Performs a least square fitting to a line from data including error
    
    Args:
        initial_index to avoid fitting the intial behaviour [Default=50]
        final_ration percentage of the data to be taken into account [Default=0.5]
        step take data every this step [Default=10]
        
    Return:
        pfinal coefficients of the linear fit
        cov covariance matrix of the fit
    """
    pinit=[1,-1]
    out = optimize.leastsq(errfunc, pinit, args=(x[initial_index:],y[initial_index:],yerr[initial_index:]), full_output=1)
    pfinal = out[0] #fitting coefficients
    cov=out[1] #Covariance
    
    return pfinal,cov    

def lammps_MSD(delta_t, data):
    """
    delta_t from the simulations
    data is a pandas data frame which contains in the first column the timestep
    and the second the msd
    """
    data = data.values
    times = data[:,0]-data[0,0]
    times = times*delta_t
    msd = data[:,1]

    fig,ax = plt.subplots()

    ax.plot(times,msd,label="LAMMPS")

    out = np.polyfit(times,msd,1)

    ax.plot(times,out[0]*times,label="fit")
    
    ax.legend()
    ax.set_xlabel(r'$\Delta t(fs)$')
    ax.set_ylabel(r'$MSD[{\AA}^2]$')
    plt.savefig("msd.pdf")
    D = out[0]/(2*3)
    error = out[1]/(2*3)

    print("The diffusion coefficient from Lammps MSD is %s +/- %s"%(D,error))
    
    return times,msd

def plot_diffusion(t, msd_average, msd_error, D_inst_ave, D_inst_error, pfinal, D, initial_index, dim):
    """
    
    """
    cf.set_plot_appearance()
    plt.close('all')
    fig1,(ax1,ax12)=plt.subplots(2,1, sharex='col')
    ax1.plot(t, msd_average)
    ax1.fill_between(t, msd_average-msd_error, msd_average+msd_error ,alpha=0.4)
    ax1.plot(np.unique(t),fitfunc(pfinal,np.unique(t)),linestyle='--',c='black')
    ax1.plot(t[initial_index], msd_average[initial_index], marker = 'o')
    ax1.set_ylabel(r'$MSD [{\AA}^2]$')
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax12.plot(t,D_inst_ave)
    ax12.fill_between(t, D_inst_ave-D_inst_error, D_inst_ave+D_inst_error ,alpha=0.4)
    ax12.axhline(y = D, xmin=0, xmax=1,ls='--',c='black', label =r'$D = %2.3f$'%D )
    ax12.set_xlabel(r'$\Delta t[fs]$')
    ax12.set_ylabel(r'$D[{\AA}^2/fs]$')
    ax12.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.legend(loc = "lower right", fontsize = 10)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0.1)
    plt.savefig("Diffusio_coefficient%s.pdf"%dim)
    


def compute_one_msd(pos_init, pos_final):
    """
    Computes the msd between two positons
    Returns the 3 components and the total
    """

    delta_sqr_components = (pos_final-pos_init)**2
    msd = np.average(delta_sqr_components,axis=0)
    msd = np.append(msd,np.sum(msd))

    return msd
    
def one_delta_t(delta, positions, max_delta):
    """
    returns an array with all the msd for a given delta_t
    TODO [Note that it assumes that centroids_traj and max_delta are loaded in memory, so it can be accesed by all the threads, probably not the most efficeint way of doing]
    
    Args:
        delta is the delta in sampling times that is going to be analysed
    """
    num_cores = multiprocessing.cpu_count()    
    pos = zip(positions[:max_delta],positions[delta:max_delta+delta])
    msd_array_t = Parallel(n_jobs = num_cores)(delayed(compute_one_msd)(*p) for p in pos) # * unzips 
    return np.array(msd_array_t)



def one_delta_t_parallel(delta, centroids_traj, max_delta):
    """
    returns an array with all the msd for a given delta_t
    TODO [Note that it assumes that centroids_traj and max_delta are loaded in memory, so it can be accesed by all the threads, probably not the most efficeint way of doing]
    
    Args:
        delta is the delta in sampling times that is going to be analysed
    """
    msd_array_t = []
    for j in range(max_delta):        
        msd_array_t.append(compute_one_msd(centroids_traj[j,:,:],centroids_traj[j+delta,:,:]))
    return msd_array_t

def msd(positions, max_delta):
    dim = np.shape(positions)
    msd_array = []
    for i in tqdm(range(1,max_delta)):
        msd_array.append(one_delta_t(i, positions, max_delta))
    # The first one is zero in all dimensions
    msd_array.insert(0,np.zeros((max_delta,dim[-1]+1))) 
    cf.save_instance(msd_array,"msd_array")
    return msd_array
    