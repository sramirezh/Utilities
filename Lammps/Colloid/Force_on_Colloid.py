#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:55:46 2020
measurement of the force on the colloid 
@author: sr802
"""




import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
from joblib import Parallel, delayed
import multiprocessing

from Lammps.PDP.Plots.LJ_plotter import LJ
import Lammps.core_functions as cf
import Others.Statistics.FastAverager as stat
#import Lammps.PDP.trajectory_analysis.poly_analysis as pa

from ovito.modifiers import *
import ovito.io as ov

#
#import sys
#import time
#from PyQt5.QtCore import *
#from PyQt5.QtGui import *
#from PyQt5.QtWidgets import *
#
#app = QApplication(sys.argv)
#label = QLabel("message")
#label.setWindowFlags(Qt.SplashScreen)
#label.show()
#QTimer.singleShot(6000, app.quit)
#app.exec_()


def fdivr(r, epsilon, sigma):
    """
    See Frenekl Pag 68-69
    computes the lennard-jones force divided by r (f/r)
    Args:
        r       : the position at which the potential is computed..
        epsilon : Energy of the minimum.
        sigma   : Radius at which the potential goes to zero.
    Returns:
        F Lennard Jones force
    """
    r2 = sigma/r**2
    r6 = r2**3
    virial = 48*epsilon*r6*r2*(r6 - 0.5)

    return virial


def spherical_coordinates(X):
    """
    Converts from cartesian to spherical coordinates (With the elevation angle defined from the XY plane)
    Args:
        Vector with positions [xi,yi,zi]
    Return:
        Vector [r,tetha,phi] 
        tetha and phi in radians
    """
    n,m=np.shape(X)
    new_X=np.zeros((n,3))
    xy=X[:,0]**2+X[:,1]**2
    new_X[:,0]=np.sqrt(xy+X[:,2]**2)
    new_X[:,1]=np.arctan2(X[:,2],X[:,0]) #for elevation angle defined from the xy plane, from the x-positive axis
    new_X[:,2]=np.arctan2(X[:,1],X[:,0]) # phi defined from the x axis 

    return new_X

r_colloid  = 3.23

node = ov.import_file("all.atom", multiple_frames = True)


#Getting the properties of the box

#Assuming that the box does not change size

data = node.compute(0)
list(data.particle_properties.keys())

box = node.compute(0).cell # getting the properties of the box
L = np.diag(box.matrix) 
center = L/2.0
r_cut = 2.5
r_max = r_colloid+r_cut


# LJ parameters
epsilon1 = 1
epsilon2 = 1.5
sigma2 = r_colloid
sigma1 = r_colloid

print ("Beware that epsilon is %s"%epsilon2)


discard = 0

n_frames = node.source.num_frames
f_x_time = []
solu_count_time = []


# To use either total number or percentage to discard
if discard< 1:
    discard = discard*n_frames 

for frame in range(n_frames):
    if frame >= discard:
        print ("Analysing frame %s of %s)"%(frame,n_frames))
        data = node.compute(frame)
        pos = data.particle_properties.position.array
        types = data.particle_properties.particle_type.array
        relative_pos = pos-center
        
        
        r = np.linalg.norm(relative_pos,axis = 1 )
        
        
        #Limiting by radius
        ind_r = np.where(r < r_max)[0]
        relative_pos = relative_pos[ind_r]
        types = types[ind_r]
        pos_sphe = spherical_coordinates(relative_pos)
        
        
        
        ind_solv = np.where(types == 1)
        ind_solu = np.where(types == 2)
        
        
        
        n_theta = 15 # Number of partitions in theta defined as measured from xy plane
        theta = np.linspace(0,np.pi,180/n_theta)
        theta_low = theta[:-1]
        theta_high = theta[1:]
        
        
        indexes = []

        
        

        
        #Angular distribution of partilces and forces
        
        force = [] #
        theta = []
        solute = [] 
        
        for k,tl in enumerate(theta_low):
            index = (np.where((np.abs(pos_sphe[:,1]) > tl) & (np.abs(pos_sphe[:,1]) <= theta_high[k]))[0])
            indexes.append(index)
            f_x = 0 
            theta.append((tl+theta_high[k])/2)
            solu_count = 0
            for i in index:
                if types[i] == 1:
                    r = pos_sphe[i,0]
                    dx = relative_pos[i,0]
                    f_x += - fdivr(r,epsilon1,sigma1)*dx  #The minus is because is the force that the colloid feels, not the particle
                if types[i] == 2:
                    r = pos_sphe[i,0]
                    dx = relative_pos[i,0]
                    f_x += - fdivr(r,epsilon2,sigma2)*dx 
                    solu_count += 1
            solute.append(solu_count)
            force.append(f_x)
        f_x_time.append(force)
        print (np.sum(force))
        solu_count_time.append(solute)
        
            
            
theta = np.array(theta)*180/np.pi
f_dist = np.average(f_x_time,axis = 0)
solute_dist = np.average(solu_count_time,axis = 0)
ftotal = np.sum(f_dist)

total_force_time = np.sum(f_x_time, axis =1)


print("\nThe average total force in x is %s\n"%np.average(total_force_time))
print (stat.fast_averager(total_force_time))

#
## Figures
#fig,ax = plt.subplots()
#
#
#
#ax.plot(theta,f_x_time[0], label='Force', ls = '--')
##ax.plot(x, rho_total_c, label='Total')  
#ax.axhline(y=0, xmin = 0, xmax=1, ls='--',c='black')
#ax.axvline(x=180, ymin = 0, ymax=1, ls='--',c='black')
#
#ax.set_xlabel(r"$x$")
#ax.set_ylabel(r"$fx$")
#ax.legend()
#plt.tight_layout()
#plt.savefig("force_angular.pdf")



# =============================================================================
# Plotting just to check
# =============================================================================
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#plt.close('all')
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#for indexes_positive in index:
#    ax.scatter(relative_pos[indexes_positive,0], relative_pos[indexes_positive,1], relative_pos[indexes_positive,2])
#
## draw sphere
#u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
#x = r_colloid*np.cos(u)*np.sin(v)
#y = r_colloid*np.sin(u)*np.sin(v)
#z = r_colloid*np.cos(v)
#ax.plot_wireframe(x, y, z, color="r")
#
#plt.show()
