#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 11:26:33 2020
Computes the diffusion coefficient of the molecules using MD analysis to read the configuration outputs

I will use copy some of the functions from PDP/trajectory_analysis/Diffusion_coefficient
@author: simon
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.transformations as tr
import MDAnalysis.analysis.rdf as rdf
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
from scipy.spatial.distance import pdist,squareform
from tqdm import tqdm
from uncertainties import unumpy,ufloat
import Others.Statistics.FastAverager as stat
from scipy import optimize
import glob
from joblib import Parallel, delayed
import multiprocessing

def compute_one_msd(pos_init,pos_final):
    """
    Computes the msd between two positons
    """

    delta_sqr_components = (pos_final-pos_init)**2
    delta_sqr = np.sum(delta_sqr_components,axis=1)

    #msd_comp = np.average(delta_sqr_components,axis=0)
    msd = np.average(delta_sqr)

    return msd

def lammps_MSD(delta_t, data):
    """
    delta_t from the simulations
    data is a pandas data frame which contains in the first column the timestep and the second the msd
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


def plot_diffusion(t, msd_average, msd_error,D_inst_ave, D_inst_error, pfinal, D, initial_index):
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
    plt.savefig("Diffusio_coefficient.pdf")
    
    
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
    
    
    
    


def compute_centroids():
    """
    Gets all the positions and images and saves an object with the positions of the centroids
    """
    u = mda.Universe("system.data", "dcd_nvt.dcd")  # The universe for all the atoms
    v = mda.Universe("system.data","per_image.dat", format = "LAMMPSDUMP" ) # Reading all the periodic images
    
    
    n_molecules = u.atoms.n_residues
    time_steps = u.trajectory.n_frames
    centroids_traj = np.empty((time_steps, n_molecules, 3 ))
    
    
    # =============================================================================
    #  wrapping the coordinates
    # =============================================================================
    ag = u.atoms
    transform = tr.wrap(ag)
    u.trajectory.add_transformations(transform)
    
    
    centroids_traj = np.empty((time_steps, n_molecules, 3 )) # Will contain all the centroid including the effect of the image
    
    for i,ts in enumerate(tqdm(u.trajectory, file = sys.stdout)):
        
        v.trajectory[i]    
        
        # Converting into real units including the images
        u.atoms.positions = u.atoms.positions+v.atoms.positions 
        
        centroids_traj[i, :,:] = u.atoms.residues.atoms.centroid(compound = 'residues')
        
    # universe with the real position of all the centroids 
    # Notice that this can be packed again to be inside the box with u_new.atoms.pack_into_box(box = u.dimensions)
    
#    u_new = mda.Universe.empty(n_molecules, trajectory = True)
#    u_new.load_new(centroids_traj)
    
    cf.save_instance(centroids_traj,"centroids_traj")
    
    return centroids_traj, time_steps


# =============================================================================
# Main
# =============================================================================
    


if os.path.exists("centroids_traj.pkl"):
    print ("Reading 'centroids_traj.pkl'")
    centroids_traj = cf.load_instance("centroids_traj.pkl")
    time_steps = len(centroids_traj)
else:
    centroids_traj, time_steps = compute_centroids()




max_delta = int(time_steps*0.5) #Maximum delta of time to measure the MSD as per Keffer2001

delta_t_arr =[]
msd_array = []

# Sampling interval I am printing configurations every 100 time steps and time step equal 10 fs
mult_t = 100*10



for i in range(max_delta):
    print("analysing delta = %s" %i)
    msd_array_t = []
    delta_t_arr.append(i*mult_t)
    
    for j in range(max_delta):        
        msd_array_t.append(compute_one_msd(centroids_traj[j,:,:],centroids_traj[j+i,:,:]))
        
    msd_array.append(msd_array_t)
    

# Getting the error
    

ave_msd =[]
for el in msd_array:
    
    ave = (stat.fast_averager(np.array(el)))[0]
    ave_msd.append(ufloat(ave[1],ave[3])) #Average and blocking error


  
    

D_inst=[0] #Array with the instantaneous diffusion coefficient
for i in range(1,max_delta):
    dt = delta_t_arr[i]
    D_inst.append(ave_msd[i]/dt/(2*3))



#Writing arrays of averages and errors
t = np.array(delta_t_arr)
msd_error = unumpy.std_devs(ave_msd)
msd_average = unumpy.nominal_values(ave_msd)


D_inst_error = unumpy.std_devs(D_inst)
D_inst_ave = unumpy.nominal_values(D_inst)



pfinal,cov = fit_line(t,msd_average,msd_error, initial_index = 500)

D = pfinal[0]/(2*3)

D_err=np.sqrt(cov[0][0])*D

plot_diffusion(t,msd_average,msd_error,D_inst_ave,D_inst_error,pfinal, D, initial_index = 500 )

print("\nThe diffusion coefficient is %s +/- %s"%(D,D_err))
f = open("Diffusion.out",'w')
f.write("The diffusion coefficient is %s +/- %s \n"%(D,D_err))

f.close


## =============================================================================
## From lammps chunk/msd, which only has one origin
## =============================================================================
#
#cf.set_plot_appearance()
#
delta_t = 10 # fs
#
#print ("Delta t in the simulations is %s"%delta_t)
data_lammps = cf.read_data_file('diffusion_data.dat')


times_l,msd_l = lammps_MSD(delta_t,data_lammps)



fig,ax = plt.subplots()

ax.plot(times_l,msd_l,label="Single Origin")
ax.plot(t, msd_average, label = "Multiple Origin",ls='--')
ax.legend()
ax.set_xlabel(r'$\Delta t(fs)$')
ax.set_ylabel(r'$MSD[{\AA}^2]$')
plt.tight_layout()
plt.savefig("MSD_Lammps_comparison.pdf", transparent = True)



