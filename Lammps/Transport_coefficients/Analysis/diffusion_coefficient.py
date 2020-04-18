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

# =============================================================================
# Main
# =============================================================================
u = mda.Universe("system.data", "dcd_nvt.dcd")  # The universe for all the atoms
v = mda.Universe("system.data","per_image.dat", format = "LAMMPSDUMP" ) # Reading all the periodic images

scale = v.dimensions[0]

n_molecules = u.atoms.n_residues
molecules = u.atoms.residues # I can select them here and then iterate later through the trajectory

time_steps = u.trajectory.n_frames

centroids_traj = np.empty((time_steps, n_molecules, 3 ))

distances = []

# =============================================================================
#  Unwtrapping the coordinates
# =============================================================================
ag = u.atoms
transform = tr.wrap(ag)
u.trajectory.add_transformations(transform)


centroids_traj = np.empty((time_steps, n_molecules, 3 )) # Will contain all the centroid including the effect of the image

for i,ts in enumerate(tqdm(u.trajectory, file = sys.stdout)):
    
    v.trajectory[i]    
    # Converting into real units
    u.atoms.positions = u.atoms.positions+v.atoms.positions 
    
    centroids_traj[i, :,:] = u.atoms.residues.atoms.centroid(compound = 'residues')
    
# universe with the real position of all the centroids 
# Notice that this can be packed again to be inside the box with u_new.atoms.pack_into_box(box = u.dimensions)

u_new = mda.Universe.empty(n_molecules, trajectory = True)
u_new.load_new(centroids_traj)



max_delta = int(time_steps*0.20) #Maximum delta of time to measure the MSD

delta_t_arr =[]
msd_array = []

# Sampling interval I am printing configurations every 100 time steps and time step equal 10 fs
mult_t = 100*10

for i in range(max_delta):
    msd_array_t = []
    delta_t_arr.append(i*mult_t)
    for j in range(len(centroids_traj[::i+1])):
        print (i,j)
        
        msd_array_t.append(compute_one_msd(centroids_traj[j,:,:],centroids_traj[j+i,:,:]))
        
    msd_array.append(msd_array_t)
    

ave_msd =[]
for el in msd_array:
    
    ave_msd.append(np.average(np.array(el)))


# =============================================================================
# From lammps chunk/msd, which only has one origin
# =============================================================================

cf.set_plot_appearance()

delta_t = 10 # fs

print ("Delta t in the simulations is %s"%delta_t)
data_lammps = cf.read_data_file('diffusion_data.dat')


times_l,msd_l = lammps_MSD(delta_t,data_lammps)



fig,ax = plt.subplots()

ax.plot(times_l,msd_l,label="LAMMPS")
ax.plot(delta_t_arr, ave_msd, label = "Mine",ls='--')
ax.legend()
ax.set_xlabel(r'$\Delta t(fs)$')
ax.set_ylabel(r'$MSD[{\AA}^2]$')
plt.tight_layout()
plt.savefig("MSD_Lammps_comparison.pdf", transparent = True)



