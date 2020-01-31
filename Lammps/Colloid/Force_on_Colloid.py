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

#import Lammps.PDP.trajectory_analysis.poly_analysis as pa

from ovito.modifiers import *
import ovito.io as ov


import sys
import time
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

app = QApplication(sys.argv)
label = QLabel("message")
label.setWindowFlags(Qt.SplashScreen)
label.show()
QTimer.singleShot(6000, app.quit)
app.exec_()


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

node = ov.import_file("equil.dat", multiple_frames = False)


#Getting the properties of the box

#Assuming that the box does not change size

data = node.compute(0)
list(data.particle_properties.keys())

box = node.compute(0).cell # getting the properties of the box
L = np.diag(box.matrix) 
center = L/2.0


pos = data.particle_properties.position.array
types = data.particle_properties.particle_type.array


r_cut = 2.5

relative_pos = pos-center


r = np.linalg.norm(relative_pos,axis = 1 )


#Limiting by radius
ind_r = np.where(r < r_colloid+r_cut)[0]


relative_pos = relative_pos[ind_r]

r_sphe = spherical_coordinates(relative_pos)


indexes_positive = np.where(np.abs(r_sphe[:,1]) < 1)[0]


plt.scatter(relative_pos[indexes_positive,0],relative_pos[indexes_positive,2])



from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(relative_pos[indexes_positive,0], relative_pos[indexes_positive,1], relative_pos[indexes_positive,2])

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = r_colloid*np.cos(u)*np.sin(v)
y = r_colloid*np.sin(u)*np.sin(v)
z = r_colloid*np.cos(v)
ax.plot_wireframe(x, y, z, color="r")


