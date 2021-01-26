#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 17:21:33 2019
Plot the streamlines and the solute concentration distribution from 2d distributions of velocities obtained in cylindrical coordinates.


compute cc1 all chunk/atom bin/cylinder x center 0.5 ${yo} ${zo} 0.0 10.0 20

Cylinder with radius 10 and 20 circle bins in that direction

in the x direction, the partitions are 0.5

bin/cylinder args = dim origin delta c1 c2 rmin rmax ncbin
    dim = x or y or z = axis of cylinder axis
    origin = lower or center or upper or coordinate value (distance units)
    delta = thickness of spatial bins in dim (distance units)
    c1,c2 = coords of cylinder axis in other 2 dimensions (distance units)
    crmin,crmax = bin from cylinder radius rmin to rmax (distance units)
    ncbin = # of concentric circle bins between rmin and rmax
    

yo=Ly/2
zo=Lz/2
compute cc1 all chunk/atom bin/cylinder x center 0.5 ${yo} ${zo} 0.0 10.0 20
The first dimension is along the cylinder axis, the second dimension is radially outward from the cylinder axis. 

in the output 

The first dimension is along the cylinder axis, the second dimension is radially outward from the cylinder axis. 

The created bins (and hence the chunk IDs) are numbered consecutively from 1 to the number of bins = Nchunk
For bin/cylinder, the numbering varies most rapidly in the dimension along the cylinder axis and most slowly in the radial direction
@author: simon
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from joblib import Parallel, delayed
import multiprocessing

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Others.Statistics.FastAverager as stat


from ovito.modifiers import *
import ovito.io as ov

def semi_circle(origin,r):
    tetha=np.linspace(0,np.pi,1000)
    x=r*np.cos(tetha)+origin[0]
    y=r*np.sin(tetha)+origin[1]
    
    return x,y


def data_contour(x,y,z):
    """
    Retourn data ready for a contourf plot
    """
    n_x=len(np.unique(x))
    n_y=len(np.unique(y))
    
    x=np.reshape(x,(n_y,n_x))
    y=np.reshape(y,(n_y,n_x))
    z=np.reshape(z,(n_y,n_x))
    
    return x,y,z

data = np.loadtxt("prof2d_vel.dat",skiprows=4)

data_rho = np.loadtxt("prof2d_con.dat",skiprows=4)

#data_f = np.loadtxt("prof2d_force.dat", skiprows = 4)

R_h = 3.23 #4.67520164
print((("Remember to define the Hydrodynamic radius, at the moment it is %s")%R_h))


    


#Reading data and removing anomalous
r = data[:,1]
x = data[:,2]

n_x = len(np.unique(x))
n_r = len(np.unique(r))


#Deleting the last 
r = r[:-n_x]
x = x[:-n_x]
density = data_rho[:-n_x,4]
vx = data[:-n_x,5]
vr = data[:-n_x,6]

#fx = data_f[:-n_x,5]
#fr = data_f[:-n_x,6]

xmesh,rmesh,density = data_contour(x,r,density)



# Getting the box properties using Ovito
node = ov.import_file("equil.dat", multiple_frames = False)
box = node.compute(0).cell # getting the properties of the box
L = np.diag(box.matrix) 
center = L/2.0


circle = semi_circle([center[0],0],R_h)

cf.set_plot_appearance()

# =============================================================================
# Plotting
# =============================================================================

"""
Velocity field and density contour
Here the velocities do  have the polymer velocity substracted, as it was done in LAMMPS in the in.velprof
"""
delta_r=0.2
fig,(cax,ax)=plt.subplots(nrows=2,gridspec_kw={"height_ratios":[0.05, 1]})
plt.close('all')
ax.axes.set_aspect('equal')
cntr1=ax.contourf(xmesh,rmesh,density,alpha=0.8,cmap="RdBu_r") #cnap also could be set
cbar=fig.colorbar(cntr1, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=10) 
cax.set_xlabel(r'$c_s$')
cax.xaxis.set_label_position('top') 
ax.quiver(x,r,vx,vr)
ax.plot(circle[0],circle[1],color='black')
ax.set_ylabel(r'$r=\sqrt{y^2+z^2}$')
ax.set_xlabel(r'$x$')
ax.set_ylim(np.min(r)-delta_r,np.max(r)+delta_r)
fig.tight_layout()
fig.savefig('vfield.pdf', transparent = True)



# """
# Force field and density contour
# """
# delta_r=0.2
# fig,(cax,ax)=plt.subplots(nrows=2,gridspec_kw={"height_ratios":[0.05, 1]})
# plt.close('all')
# ax.axes.set_aspect('equal')
# cntr1=ax.contourf(xmesh,rmesh,density,alpha=0.8,cmap="RdBu_r") #cnap also could be set
# cbar=fig.colorbar(cntr1, cax = cax, orientation='horizontal')
# cbar.ax.tick_params(labelsize=10) 
# cax.set_xlabel(r'$c_s$')
# cax.xaxis.set_label_position('top') 
# ax.quiver(x,r,fx,fr,scale =200)
# ax.plot(circle[0],circle[1],color='black')
# ax.set_ylabel(r'$r=\sqrt{y^2+z^2}$')
# ax.set_xlabel(r'$x$')
# ax.set_ylim(np.min(r)-delta_r,np.max(r)+delta_r)
# fig.tight_layout()
# fig.savefig('force_field.pdf')



#"""
##average velocity 
#"""
#data2=cf.read_data_file('vel_sol.dat').get_values()
##Substracting the velocity of the polymer
#for i in range(1,4):
#    data2[:,i]=data2[:,i]-data2[:,3]
#
#fig,ax=plt.subplots()
#
#
#averages=stat.fast_averager(data2[:,1:3],min_limit=0.3)
#every=200
#ax.plot(data2[::every,0],data2[::every,1],label='Bulk',c='r')
#ax.plot(data2[::every,0],data2[::every,2],label='Inside',c='b')
#ax.axhline(y=averages[0][1], xmin=0, xmax=1,ls=':',c='black')
#ax.axhline(y=averages[1][1], xmin=0, xmax=1,ls=':',c='white')
#
#ax.set_ylabel(r'$v_x-v_x^{cm}$')
#ax.set_xlabel(r'$step $')
#plt.legend()
#fig.tight_layout()
#fig.savefig('vx.pdf')
#
#fig.show()
#vel_inside=stat.fast_averager(data2[:,2]-data2[:,1],min_limit=0)[0]
#print(vel_inside)



"""
#Only Fluid inside polymer
Substracting the average bulk velocity from all the others
"""


delta_r=0.2
fig,(cax,ax)=plt.subplots(nrows=2,gridspec_kw={"height_ratios":[0.05, 1]})
plt.close('all')
ax.axes.set_aspect('equal')
cntr1=ax.contourf(xmesh,rmesh,density,alpha=0.8,cmap="RdBu_r") #cnap also could be set
cbar=fig.colorbar(cntr1, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=15) 
cax.set_xlabel(r'$c_s$')
cax.xaxis.set_label_position('top') 
ax.arrow( 10, 0.42, -3.0,0, fc="k", ec="k",head_width=0.6, head_length=1,width=0.3 )
ax.plot(circle[0],circle[1],color='black')
ax.set_ylabel(r'$r=\sqrt{y^2+z^2}$')
ax.set_xlabel(r'$x$')
ax.set_ylim(np.min(r)-delta_r,np.max(r)+delta_r)
fig.tight_layout()
fig.savefig('vel_pol.pdf', transparent = True)



"""
Vx contour
The velocity of the fluid without substracting the bulk velocity, in the x direction
"""
xmesh,rmesh,vx_contour=data_contour(x,r,vx)
delta_r=0.2
fig,(cax,ax)=plt.subplots(nrows=2, gridspec_kw={"height_ratios":[0.05, 1]})
plt.close('all')
ax.axes.set_aspect('equal')
cntr1=ax.contourf(xmesh,rmesh,vx_contour,alpha=0.8,cmap="RdBu_r") #cnap also could be set
cbar=fig.colorbar(cntr1, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=10) 
cax.set_xlabel(r'$v_x$')
cax.xaxis.set_label_position('top') 
ax.plot(circle[0],circle[1],color='black')
ax.set_ylabel(r'$r=\sqrt{y^2+z^2}$')
ax.set_xlabel(r'$x$')
ax.set_ylim(np.min(r)-delta_r,np.max(r)+delta_r)
fig.tight_layout()
fig.savefig('vx_contour.pdf', transparent = True)


"""
V contour
The velocity of the fluid without substracting the bulk velocity
"""
v=np.sqrt(vx**2+vr**2)
xmesh,rmesh,v_contour=data_contour(x,r,v)
delta_r=0.2
fig,(cax,ax)=plt.subplots(nrows=2, gridspec_kw={"height_ratios":[0.05, 1]})
plt.close('all')
ax.axes.set_aspect('equal')
cntr1=ax.contourf(xmesh,rmesh,v_contour,alpha=0.8,cmap="RdBu_r") #cnap also could be set
cbar=fig.colorbar(cntr1, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=10) 
cax.set_xlabel(r'$|v|$')
cax.xaxis.set_label_position('top') 
ax.plot(circle[0],circle[1],color='black')
ax.set_ylabel(r'$r=\sqrt{y^2+z^2}$')
ax.set_xlabel(r'$x$')
ax.set_ylim(np.min(r)-delta_r,np.max(r)+delta_r)
fig.tight_layout()
fig.savefig('v_contour.pdf', transparent = True)




# """
# F Contour
# The velocity of the fluid without substracting the bulk velocity
# """
# f = fx
# xmesh,rmesh,f_contour=data_contour(x,r,f)
# delta_r=0.2
# fig,(cax,ax)=plt.subplots(nrows=2, gridspec_kw={"height_ratios":[0.05, 1]})
# plt.close('all')
# ax.axes.set_aspect('equal')
# cntr1=ax.contourf(xmesh,rmesh,f_contour,alpha=0.8,cmap="RdBu_r") #cnap also could be set
# cbar=fig.colorbar(cntr1, cax=cax, orientation='horizontal')
# cbar.ax.tick_params(labelsize=10) 
# cax.set_xlabel(r'$f$')
# cax.xaxis.set_label_position('top') 
# ax.plot(circle[0],circle[1],color='black')
# ax.set_ylabel(r'$r=\sqrt{y^2+z^2}$')
# ax.set_xlabel(r'$x$')
# ax.set_ylim(np.min(r)-delta_r,np.max(r)+delta_r)
# fig.tight_layout()
# fig.savefig('f_contour.pdf')



# """
# V contour
# The velocity of the fluid without substracting the bulk velocity
# """
# v=np.sqrt(vx**2+vr**2)
# xmesh,rmesh,v_contour=data_contour(x,r,v)
# delta_r=0.2
# fig,(cax,ax)=plt.subplots(nrows=2, gridspec_kw={"height_ratios":[0.05, 1]})
# plt.close('all')
# ax.axes.set_aspect('equal')
# cntr1=ax.contourf(xmesh,rmesh,v_contour,alpha=0.8,cmap="RdBu_r") #cnap also could be set
# cbar=fig.colorbar(cntr1, cax=cax, orientation='horizontal')
# cbar.ax.tick_params(labelsize=10) 
# cax.set_xlabel(r'$|v|$')
# cax.xaxis.set_label_position('top') 
# ax.plot(circle[0],circle[1],color='black')
# ax.set_ylabel(r'$r=\sqrt{y^2+z^2}$')
# ax.set_xlabel(r'$x$')
# ax.set_ylim(np.min(r)-delta_r,np.max(r)+delta_r)
# fig.tight_layout()
# fig.savefig('v_contour.pdf', transparent = True)








##"""
##Analysis of the discarding percentage
##"""
#percentages=np.linspace(0,0.9,10)
#num_cores = multiprocessing.cpu_count()
#vel_array=Parallel(n_jobs=num_cores,verbose=10)(delayed(stat.fast_averager)(data2[:,2]-data2[:,1],min_limit=perc) for perc in percentages)
#
#
##Because the list has one additional []
#vel_matrix=[]
#for el in vel_array:
#   vel_matrix.append(el[0])
#    
#
#vel_matrix=np.array(vel_matrix)
#fig,ax=plt.subplots()
#ax.errorbar(percentages,vel_matrix[:,1], yerr=vel_matrix[:,2])
#ax.set_ylabel(r'$v_{x}^{\prime in}-v_{x}^{\prime Bulk}$')
#ax.set_xlabel("Discarded percentage")
#ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
#fig.tight_layout()
#fig.savefig('error_evolution.pdf')
