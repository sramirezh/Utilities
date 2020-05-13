#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 07:54:52 2020
Old paralellisation of the code

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
from tqdm import tqdm
import diffusion_coefficient_utils as dcu
from uncertainties import unumpy
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../'))  
import Lammps.core_functions as cf


def compute_centroids():
    """
    Gets all the positions and images for the particles, then computes the 
    position of the centroid of all particles.

    Returns:
    centroids_traj: an object with the positions of the centroids of
    the particles
    time_steps: the time steps from the trajectory
    """
    # The universe for all the atoms
    u = mda.Universe("system.data", "dcd_nvt.dcd") 
    # Reading all the periodic images
    v = mda.Universe("system.data", "per_image.dat", format="LAMMPSDUMP") 
    n_molecules = u.atoms.n_residues
    time_steps = u.trajectory.n_frames

    # Will contain all the centroid including the effect of the image
    centroids_traj = np.empty((time_steps, n_molecules, 3))

    # =============================================================================
    #  wrapping the coordinates
    # =============================================================================
    ag = u.atoms
    transform = tr.wrap(ag)
    u.trajectory.add_transformations(transform)

    for i, ts in enumerate(tqdm(u.trajectory, file=sys.stdout)):
        
        v.trajectory[i]    
        # Converting into real units including the images
        u.atoms.positions = u.atoms.positions + v.atoms.positions 
        
        centroids_traj[i, :, :] = u.atoms.residues.atoms.centroid(compound='residues')
    
    np.save("centroids_traj", centroids_traj)
    
    return centroids_traj, time_steps


class simulation(object):
    """Class to define the type of simulation analysed, for example N2 or
    Octane
    """
    def __init__(self, name, ts, d, initial_index):
        """
        Args:
            name: to identify the simulation type, eg. N2
            ts: Time step in [fs]
            d: Lammps sampling interval, how many time steps does lammps print
            the trajectory.
            initial_i: after this index, it could be said that the MSD has a
            long time behavior, it appears as a blue dot in the plots.
        """
        self.name = name
        self.ts = ts  # Time in femptoseconds
        self.d = d   # sampling_interval myDump
        self.initial_i = initial_index  # Initial index, for long time behavior
        
    def print_params(self):
        print(" Using the parameters from %s"%self.name)
        print("\nUsing a sampling threshold (myDump) of  %s"%self.d)
        print("Using delta_t = %s fs" %self.ts)


# =============================================================================
# Main
# =============================================================================
    
octane = simulation("octane", 1, 100, 3000)
nitrogen = simulation("N2", 10, 100, 500)

# this is the only thing to define
sim = octane

sim.print_params()


# =============================================================================
# # Getting the centroids
# =============================================================================
if os.path.exists("centroids_traj.npy"):
    print("Reading 'centroids_traj.npy'")
    centroids_traj = np.load("centroids_traj.npy")        
    time_steps = len(centroids_traj)
else:
    centroids_traj, time_steps = compute_centroids()

dimensions = np.shape(centroids_traj)[-1]

#Maximum delta of time to measure the MSD as per Keffer2001
max_delta = int(time_steps * 0.5) 
mult_t = sim.d * sim.ts
delta_t_arr = np.arange(max_delta) * mult_t

# =============================================================================
#  Computing the MSD array
# =============================================================================
if os.path.exists("msd_array.npy"):
    print("Reading msd data")
    msd_array = np.load("msd_array.npy")
else:
    print("Computing the msd array")
    msd_array = dcu.msd_np_parallel(centroids_traj, max_delta)

# =============================================================================
# # Computing the average msd for each tau
# =============================================================================
if os.path.exists("ave_msd.pkl"):
    print("Reading ave msd")
    ave_msd = cf.load_instance("ave_msd.pkl")
else:
    print("Computing the average msd")
    ave_msd = dcu.ave_serial_no_autocorr(msd_array)

# =============================================================================
# Final computations
# =============================================================================
D_inst = [ave_msd[0]]  # Array with the instantaneous diffusion coefficient
for i in range(1, max_delta):
    dt = delta_t_arr[i]
    D_inst.append(ave_msd[i] / dt / 2)

D_inst = np.array(D_inst)

D_inst[:, -1] = D_inst[:, -1] / 3  # To account for the 3 dimension

t = np.array(delta_t_arr)

f = open("Diffusion.out", 'w')
for dim in range(dimensions + 1):
    #Writing arrays of averages and errors
    
    # Mean square displacements
    msd_error = unumpy.std_devs(ave_msd[:, dim])
    msd_average = unumpy.nominal_values(ave_msd[:, dim])
    
    # Instantaneous diffusion coefficients
    D_inst_error = unumpy.std_devs(D_inst[:, dim])
    D_inst_ave = unumpy.nominal_values(D_inst[:, dim])
    
    pfinal, cov = dcu.fit_line(t, msd_average, msd_error, initial_index=sim.initial_i)
    
    if dim == dimensions: 
        D = pfinal[0] / (2 * 3)  # For the total D, accounting for x,y,z
    else:
        D = pfinal[0] / 2
    
    D_err = np.sqrt(cov[0][0]) * D
    
    dcu.plot_diffusion(t, msd_average, msd_error, D_inst_ave, D_inst_error, 
                       pfinal, D, sim.initial_i, dim)
    
    print("\nThe diffusion coefficient is %s +/- %s [Angstrom**2/femptoseconds]"%(D, D_err))
    
    f.write("The diffusion coefficient is %s +/- %s [Angstrom**2/femptoseconds]\n"%(D, D_err))
    
f.close

# =============================================================================
# From lammps chunk/msd, which only has one origin
# =============================================================================
#
#cf.set_plot_appearance()
#
delta_t = sim.ts 
#
#print ("Delta t in the simulations is %s"%delta_t)
#times_l,msd_l = lammps_MSD(delta_t,data_lammps)
fig, ax = plt.subplots()
#ax.plot(times_l,msd_l,label="Single Origin")
ax.plot(t, msd_average, label="Multiple Origin", ls='--')
ax.legend()
ax.set_xlabel(r'$\Delta t(fs)$')
ax.set_ylabel(r'$MSD[{\AA}^2]$')
plt.tight_layout()
plt.savefig("MSD_Lammps_comparison.pdf", transparent=True)
