#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 17:21:33 2019
Streamlines plot 


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
data=np.loadtxt("prof2d_vel.dat",skiprows=4)
data_rho=np.loadtxt("prof2d_con.dat",skiprows=4)
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf



R_h=4.48
def semi_circle(origin,r):
    tetha=np.linspace(0,np.pi,1000)
    x=r*np.cos(tetha)+origin[0]
    y=r*np.sin(tetha)+origin[1]
    
    return x,y
    
    
finish=760

r=data[:finish,1]
x=data[:finish,2]
density=data_rho[:finish,4]
vx=data[:finish,5]
vr=data[:finish,6]


n_x=len(np.unique(x))
n_r=len(np.unique(r))

X,Y=np.meshgrid(x,r)
Z=np.reshape(density,(40,19))

#xreshaped=np.transpose(np.reshape(x,(20,40)))
#rreshaped=np.transpose(np.reshape(r,(20,40)))
#vxreshaped=np.transpose(np.reshape(vx,(20,40)))
#vrreshaped=np.transpose(np.reshape(vr,(20,40)))

xmesh,ymesh=np.mgrid[np.min(x):np.max(x):np.complex(0,n_x), np.min(r):np.max(r):np.complex(0,n_r)]
xmesh=np.reshape(x,(19,40))
ymesh=np.reshape(r,(19,40))
vxmesh=np.reshape(density,(19,40))
#
#u = -1 - xmesh**2 + ymesh
#v = 1 + xmesh - ymesh**2


#Reshaping the arrays

#x=np.reshape(x,(n_x,n_y))
#y=np.reshape(y,(n_x,n_y))
#u=np.reshape(u,(n_x,n_y))
#v=np.reshape(v,(n_x,n_y))

circle=semi_circle([10,0],R_h)

cf.set_plot_appearance()

#Plotting
fig,ax=plt.subplots()
plt.close('all')
#ax.streamplot(ymesh,xmesh, v, u)

#ax.streamplot (xreshaped,rreshaped,vxreshaped,vrreshaped)

ax.axes.set_aspect('equal')
cntr1=ax.contourf(xmesh,ymesh,vxmesh,alpha=0.5,cmap="RdBu_r")

fig.colorbar(cntr1, ax=ax)
ax.quiver(x,r,vx,vr)
ax.plot(circle[0],circle[1])
ax.set_ylim(0,9)
#ax.contourf(xreshaped,rreshaped,vxreshaped)
fig.savefig('vfield.pdf')



"""
#average velocity
"""
data2=cf.read_data_file('vel_sol.dat').get_values()
#Substracting the velocity of the polymer
for i in xrange(1,4):
    data2[:,i]=data2[:,i]-data2[:,3]

fig,ax=plt.subplots()

every=100
ax.plot(data2[::every,0],data2[::every,1],label='Solution')
ax.plot(data2[::every,0],data2[::every,2],label='Inside')
ax.plot(data2[::every,0],data2[::every,3],label='Polymer')
ax.set_ylabel(r'$v_x $')
ax.set_xlabel(r'$step $')
plt.legend()
plt.tight_layout()

fig.show()
