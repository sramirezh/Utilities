#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 11:26:33 2020
Computes the diffusion coefficient of the molecules using MD analysis to read 
the configuration outputs

I will use copy some of the functions from:
PDP/trajectory_analysis/Diffusion_coefficient
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
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm
from uncertainties import unumpy, ufloat
import Others.Statistics.FastAverager as stat
import diffusion_coefficient_utils as dcu
import glob


    
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
    
    np.save("centroids_traj",centroids_traj)
    
    return centroids_traj, time_steps


class simulation(object):
    def __init__(self, name, ts, d, initial_index):
        self.name = name
        self.ts = ts
        self.d = d   # sampling_interval myDump
        
    def print_params(self):
        print (" Using the parameters from %s"%self.name)
        print ("\nUsing a sampling threshold (myDump) of  %s"%self.d)
        print ("Using delta_t = %s fs" %self.ts)


# =============================================================================
# Main
# =============================================================================
    
octane = simulation("octane", 1, 100, 5000 )
nitrogen = simulation("N2", 10, 100, 500 )

# this is the only thing to define
sim = octane

sim.print_params()


# Getting the centroids
if os.path.exists("centroids_traj.npy"):
    print ("Reading 'centroids_traj.npy'")
    centroids_traj = np.load("centroids_traj.npy")
    time_steps = len(centroids_traj)
else:
    centroids_traj, time_steps = compute_centroids()




max_delta = int(time_steps*0.5) #Maximum delta of time to measure the MSD as per Keffer2001
mult_t = sim.d*sim.ts
delta_t_arr = np.arange(max_delta)*mult_t



# Computing the MSD array
if os.path.exists("msd_array.pkl"):
    print ("Reading msd data")
    msd_array = cf.load_instance("msd_array.pkl")
else:
    print ("Computing the msd array")
    msd_array = dcu.msd(centroids_traj, max_delta)


    



# Computing the average msd for each tau
if os.path.exists("ave_msd.pkl"):
    print ("Reading ave msd")
    ave_msd = cf.load_instance("ave_msd.pkl")
else:
    print ("Computing the average msd")
    ave_msd =[]
    for el in msd_array:
        ave = (stat.fast_averager(np.array(el)))[0]
        
        ######hdre make it for each direction
        ave_msd.append(ufloat(ave[1],ave[3])) #Average and blocking error 
    
    cf.save_instance(ave_msd,"ave_msd")


  
    

D_inst = [0] #Array with the instantaneous diffusion coefficient
for i in range(1,max_delta):
    dt = delta_t_arr[i]
    D_inst.append(ave_msd[i]/dt/(2*3))



#Writing arrays of averages and errors
t = np.array(delta_t_arr)
msd_error = unumpy.std_devs(ave_msd)
msd_average = unumpy.nominal_values(ave_msd)


D_inst_error = unumpy.std_devs(D_inst)
D_inst_ave = unumpy.nominal_values(D_inst)

# TODO This is a rough estimate, check that the blue point in the plot is correct
initial_index = int(len(t)*0.5)

pfinal,cov = dcu.fit_line(t,msd_average,msd_error, initial_index  = initial_index)

D = pfinal[0]/(2*3)

D_err=np.sqrt(cov[0][0])*D

dcu.plot_diffusion(t,msd_average,msd_error,D_inst_ave,D_inst_error,pfinal, D, initial_index  = initial_index )



print("\nThe diffusion coefficient is %s +/- %s [Angstrom**2/femptoseconds]"%(D,D_err))
f = open("Diffusion.out",'w')
f.write("The diffusion coefficient is %s +/- %s [Angstrom**2/femptoseconds]\n"%(D,D_err))

f.close


## =============================================================================
## From lammps chunk/msd, which only has one origin
## =============================================================================
#
#cf.set_plot_appearance()
#
delta_t = sim.ts # fs
#
#print ("Delta t in the simulations is %s"%delta_t)
#data_lammps = cf.read_data_file('diffusion_data.dat')


#times_l,msd_l = lammps_MSD(delta_t,data_lammps)



fig,ax = plt.subplots()

#ax.plot(times_l,msd_l,label="Single Origin")
ax.plot(t, msd_average, label = "Multiple Origin",ls='--')
ax.legend()
ax.set_xlabel(r'$\Delta t(fs)$')
ax.set_ylabel(r'$MSD[{\AA}^2]$')
plt.tight_layout()
plt.savefig("MSD_Lammps_comparison.pdf", transparent = True)



