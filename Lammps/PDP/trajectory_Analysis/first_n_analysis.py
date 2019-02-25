#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 16:23:44 2019

Reads the xyz taken from VMD and prepares a trajectory file that centers the cm of the polymer and leaves the new trajectory ready for a video
also analyses the distances between the nearest neighbors to the polymer monomers


I could do a direct analysis of the g(r) without taking into account the 
@author: sr802
"""
from __future__ import division
import numpy as np
import pandas as pd
import argparse
import linecache
import os
import sys
from scipy.spatial.distance import pdist,squareform
import re
import glob 
import warnings
import itertools




warnings.filterwarnings("ignore")

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf


try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
#    from matplotlib.backends.backend_pdf import PdfPages
except ImportError as err:
    print err


utilities_path=str(os.path.join(os.path.dirname(__file__), '../../../') )
"""
*******************************************************************************
Functions
*******************************************************************************
"""

def read_box_limits(log_name):
    """
    Reads the box limits from log.lammps
    ONLY required for .xyz not for .dump
    Args:
        None: log_name name of the log file
    returns:
        volume
        limits

    """
    out,err=cf.bash_command("""grep -n "orthogonal box" %s | awk -F":" '{print $1}' """%log_name)
    line=int(out.split()[0])
    limits=linecache.getline(log_name, line)
    limits=re.findall(r"[-+]?\d*\.?\d+", limits)
    limits=np.array(np.reshape(limits,(2,3)),dtype=float) #To have the box as with columns [x_i_min,x_i_max]
    volume=(limits[1,0]-limits[0,0])*(limits[1,1]-limits[0,1])*(limits[1,2]-limits[0,2])

    return volume,limits



def xyz_splitter(file_name):
    
    cf.bash_command("""%s/Lammps/Trajectory_Analysis/Trajectory_Splitter.sh -i %s""" %(utilities_path,file_name))
    
    return



def com_pbc(positions,box_length):
    """
    Computes the center of mass for xyz coordinates in a box with periodic boundary conditions,
    using the procedure in https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
    
    Args:
        positions a matrix with the [X Y Z] coordinates
        box_length a vector with Lx Ly Lz
    """
    
    com=np.zeros(3)
    for i in xrange(3):
        
        theta=positions[:,i]/box_length[i]*2*np.pi
        epsilon=np.cos(theta)
        sigma=np.sin(theta)
        theta_new=np.arctan2(-np.average(sigma),-np.average(epsilon))+np.pi
        com[i]=box_length[i]*theta_new/(2*np.pi)
    
    return com



def shift_coordinates(positions,box_length,shift):
    """
    Shifting coordinates in PBC
    
    args:
        positions a matrix with the [X Y Z] coordinates
        box_length a vector with Lx Ly Lz
        shift a vector with the shifs in each direction
    Returns:
        new_pos shifted coordinates
    """
    n_particles,m=np.shape(positions)
    new_pos=np.zeros((n_particles,3))
    
    for dim in xrange(3):
        for i in xrange(n_particles):
            if shift[dim]<0:
                new_pos[i,dim]=(positions[i,dim]-shift[dim])-np.rint(0.5*(positions[i,dim]-shift[dim])/box_length[dim])*box_length[dim]
            else:
                new_pos[i,dim]=(positions[i,dim]-shift[dim])-np.rint(0.5*(positions[i,dim]-shift[dim]-box_length[dim])/box_length[dim])*box_length[dim]
                
    return new_pos


def write_trajectory(file_name,frame,data):
    """
    Creates a xyz trajectory file with variable number of particles per time step
    Args:
        file_name name of the trajectory file
        frame number of the frame or timestep.
        data matrix containing [type X Y Z]
    """
    n_atoms,m=np.shape(data)
    f=open(file_name,'ab') 
    f.write("%s \n" %n_atoms)
    f.write("frame:%s \n" %frame)
    np.savetxt(f,data,fmt='%d %f %f %f')
    f.close()
    


def nnb_one_file(data):
    """
    Evaluates what is the nearest particle to each monomer per type,
    also calls the gr
    types as in the trajectory file from LAMMMPS [1-solvent, 2-solute, 3-polymer]
    
    Args:
        data: Contains the [Type X Y Z]
        
    Returns:
        nearest_solvent a list with the distances of the nearest solvents to monomers 
        nearest_solute a list with the distances of the nearest solutes to monomers
    """

    n,m=data.shape
    #Particle indexes
    p_types=np.array([1,2,3])
    p_types_frame=np.unique(data[:,0]).astype(int)
    particles=[np.where(data[:,0]==j)[0] for j in p_types]
    pos=data[:,1::]

    dist=squareform(pdist(pos))
    np.fill_diagonal(dist, 1000) #To avoid self contributions
    
    nearest_solute=[]
    nearest_solvent=[]
    if 1 in p_types_frame:
        nearest_solvent=computation_first_n(particles,p_types,dist,3,1)
    if 2 in p_types_frame:
        nearest_solute=computation_first_n(particles,p_types,dist,3,2)
        
    
    return nearest_solvent,nearest_solute


def computation_first_n(particles,p_types,dist,i,j):
    """
    
    computes the distance of the closest particle of type j around type i
    
    Returns:
        An array with all the positions of the nearest particles
    
    """
    i=np.where(p_types == i)[0][0]
    j=np.where(p_types == j)[0][0]
    
    i_axis0=[]
    i_axis1=[]
    
    for k in xrange(len(p_types)):
        if k!=i:
            i_axis0.append(particles[k])
        if k!=j:
            i_axis1.append(particles[k])
            
    dist = np.delete(dist,np.hstack(i_axis0), axis=0)
    dist = np.delete(dist,np.hstack(i_axis1), axis=1)
    
    
    ind_min=np.argmin(dist,axis=1)
    
    distances=[]
    for l,m in enumerate(ind_min):
        index=np.argmin(dist[:,m])
        if l==index:
            distances.append(dist[index,m])
    
    return distances



def computation_gr(data):
    """
    Computes the histogram for the g(r) for all the particles
    types as in the trajectory file from LAMMMPS [1-solvent, 2-solute, 3-polymer]
    
    Args:
        data: Contains the [Type X Y Z]
        
    Returns:
        nearest_solvent a list with the distances of the nearest solvents to monomers 
        nearest_solute a list with the distances of the nearest solutes to monomers
    """

    n,m=data.shape
    #Particle indexes
    p_types=np.array([1,2,3])
    particles=[np.where(data[:,0]==j)[0] for j in p_types]
    pos=data[:,1::]

    dist=squareform(pdist(pos))
    np.fill_diagonal(dist, 1000) #To avoid self contributions

        
    #Values to compare with vmd
    delta_r=0.1
    rmax=10
    nbins=int(rmax/delta_r)
    gr=g_r(particles,p_types,dist,2,3,nbins,rmax)
    
    return gr



def g_r(particles,p_types,dist,i,j,nbins, rmax):
    """
    computes the gr for particles type j around type i
    
    args:
        i,j particles types as in the trajectory file from LAMMMPS [1-solvent, 2-solute, 3-polymer]
        
    Returns:
        bin_count contains
        -the position of the center of the bin.
        -The number of j particles in the shell. 
        -The density of particles j in the bin
    """
    i=np.where(p_types == i)[0][0]
    j=np.where(p_types == j)[0][0]
    
    
    if len(p_types)>1:
        #indexes to delete if there is more than one type of particles
        i_axis0=[]
        i_axis1=[]
        for k in xrange(len(p_types)):
            if k!=i:
                i_axis0.append(particles[k])
            if k!=j:
                i_axis1.append(particles[k])
        dist = np.delete(dist,np.hstack(i_axis0), axis=0)
        dist = np.delete(dist,np.hstack(i_axis1), axis=1)
    
    
    
    bin_count = np.zeros((nbins,3))
    #bin_ends = -rmax*np.cos(np.linspace(np.pi/2,np.pi,num=nbins+1))
    bin_ends=np.linspace(0,rmax,num=nbins+1)

    for i in xrange(nbins):
        bin_count[i,0]=0.5*(bin_ends[i+1]+bin_ends[i]) #Count position in the middle of the bin only needed in the first
        rmax_bin=bin_ends[i+1]  
        indexes=np.where(dist<=rmax_bin)
        dist[indexes]=1000
        bin_count[i,1]=len(indexes[0])/len(particles[j])
    print sum(bin_count[:,1])
        
    

    delta_r=rmax/(nbins)       
    vol_total=(4/3*np.pi*rmax**3)
    n_total=np.sum(bin_count[:,1])
    
    return bin_count





###############################################################################
# Main    
###############################################################################

parser = argparse.ArgumentParser(description='This script evaluates the trajectory file',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file_name', metavar='InputFile',help='Input filename',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-log', help='lammps logfile name to read the box limits', default="log.lammps", type=str)
parser.add_argument('-split', help='True if trajectory file need to be splitter', default=False, type=bool)

args = parser.parse_args()


if args.split==True:
    xyz_splitter(args.file_name)
else:
    print "The Trajectory file was not splitted"


trajectory_name='new_trajectory.xyz'

open(trajectory_name, 'w').close()


files = glob.glob("*.cxyz")
files.sort(key=lambda f: int(filter(str.isdigit, f)))


nearest_solvent=[]
nearest_solute=[]
h_r=[] #as defined by allen and Tindsley

for file_name in files:
    print file_name
    data=pd.read_csv(file_name,sep=" ",dtype=np.float64,skiprows=2,header=None).values

    volume,limits=read_box_limits(args.log)
    box_length=limits[1,:]-limits[0,:]

    poly_coords=data[np.where(data[:,0]==3)[0],:][:,1::]
    all_coords=data[:,1::]
    
    #Getting the image 

    #Cm in x with respect to the first monomer
    n_monomers,m=np.shape(poly_coords)
    dist=np.zeros(n_monomers)
    
    #Main algorithm
    com_vect=com_pbc(poly_coords,box_length)
    
    shift=com_vect-box_length*0.5
    
    new_pos=shift_coordinates(all_coords,box_length,shift)
    data[:,1::]=new_pos
    write_trajectory(trajectory_name,cf.extract_digits(file_name)[0],data)
    
    #Building the Histogram
    
    solvent,solute=nnb_one_file(data)
    nearest_solvent.append(solvent)
    nearest_solute.append(solute)
    
    
    #Building the g(r)
    #Make this self contained
    
    result_gr=computation_gr(data)
    h_r.append(result_gr)


# =============================================================================
# Normalizing the g(r)
# Using the procedure in Allen Ch 8.2
# =============================================================================

n_average=np.sum(np.average(h_r,axis=0),axis=0) #Average number of particles as the trajectory has it variable
gr=h_r/1
    

nearest_solvent=list(itertools.chain(*nearest_solvent))    
nearest_solute=list(itertools.chain(*nearest_solute))  



# Creating the histogram
plt.close('all')
plt.hist(nearest_solvent,bins='auto')
plt.hist(nearest_solute,bins='auto')
plt.figure()
plt.plot(g_r[:,0],g_r[:,2],label="poly-solute")
plt.show()





#file_name="14.cxyz"
#data=pd.read_csv(file_name,sep=" ",dtype=np.float64,skiprows=2,header=None).values
#
#
#volume,limits=read_box_limits(args.log)
#box_length=limits[1,:]-limits[0,:]
#
#poly_coords=data[np.where(data[:,0]==3)[0],:][:,1::]
#all_coords=data[:,1::]
##Getting the image 
#
##Cm in x with respect to the first monomer
#n_monomers,m=np.shape(poly_coords)
#dist=np.zeros(n_monomers)
#
##Main algorithm
#com_vect=com_pbc(poly_coords,box_length)
#
#shift=com_vect-box_length*0.5
#
#new_pos=shift_coordinates(poly_coords,box_length,shift)
















#
#

#plt.plot(poly_coords[:,0],poly_coords[:,1],'o')
#plt.plot(new_pos[:,0],new_pos[:,1],'or')
##plt.xlim(limits[:,0])
##plt.ylim(limits[:,1])
#plt.show()

    
    

