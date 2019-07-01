#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 15:37:04 2018
This scripts analyses the pair correlation function from the *.gz files
@author: simon
"""


from __future__ import division
import numpy as np
import pandas as pd
import argparse
import os
import sys
import glob
import poly_analysis as pa
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf


from scipy.spatial.distance import pdist,squareform
import time

from joblib import Parallel, delayed
import multiprocessing

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script


def computation_gr(particles,p_types,dist,i,j,nbins, rmax):
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
    bin_ends = -rmax*np.cos(np.linspace(np.pi/2,np.pi,num=nbins+1))

    vol_old=0
    for i in xrange(nbins):
        bin_count[i,0]=0.5*(bin_ends[i+1]+bin_ends[i]) #Count position in the middle of the bin only needed in the first
        rmax_bin=bin_ends[i+1]
        indexes=np.where(dist<=rmax_bin)
        dist[indexes]=1000
        bin_count[i,1]=len(indexes[0])/len(particles[j])
        print len(particles[j])
        vol_new=4/3*np.pi*rmax_bin**3
        bin_count[i,2]=bin_count[i,1]/(vol_new-vol_old)

    rho_ave=256/6.71838**3 #np.sum(bin_count[:,1])/(4/3*np.pi*rmax**3)

    print rho_ave

    bin_count[:,2]=bin_count[:,2]/rho_ave**2  #g(r)=rho(r)/rho_ave

    return bin_count


def compute_one_configuration_not_mine(fil):
    """
    Computes the g(r) of the polymer particles with other particle types, so only 3
    three distri
    """


    Data=pd.read_csv(fil,sep=" ",skiprows=9,dtype=np.float64,header=None).values[:,:-1]
    n,m=Data.shape
    pos=pa.real_position(Data,L) #Real positions of all the atoms

    #Particle indexes
    p_types=np.unique(Data[:,1]).astype(int)
    particles=[np.where(Data[:,1]==j)[0] for j in p_types]

    dist=squareform(pdist(pos))
    np.fill_diagonal(dist, 1000) #To avoid self contributions
    print np.shape(dist)

    results=np.zeros((3,nbins,3)) #Factorial is the number of gr that we have

    cont=0

    for m in p_types:
            results[cont]=computation(particles,p_types,dist,1,m,nbins, rmax)
            cont+=1
    return results


def compute_one_configuration(fil):
    """
    Computes the g(r) of the polymer particles with other particle types, so only 3
    three distri
    """


    Data=pd.read_csv(fil,sep=" ",skiprows=9,dtype=np.float64,header=None).values[:,:-1]
    n,m=Data.shape
    pos=pa.real_position(Data,L) #Real positions of all the atoms

    #Particle indexes
    p_types=np.unique(Data[:,1]).astype(int)
    particles=[np.where(Data[:,1]==j)[0] for j in p_types]

    dist=squareform(pdist(pos))
    np.fill_diagonal(dist, 1000) #To avoid self contributions
    print np.shape(dist)

    results=np.zeros((3,nbins,3)) #Factorial is the number of gr that we have

    cont=0

    for m in p_types:
            results[cont]=computation(particles,p_types,dist,1,m,nbins, rmax)
            cont+=1
    return results


"""
*******************************************************************************
Main
*******************************************************************************
"""

imin=0 #Number of configurations to skip
nbins=100
rmax=2.5

#input_files = glob.glob("*1*.gz")
#input_files.sort(key=lambda f: int(filter(str.isdigit, f)))
#times=cf.extract_digits(input_files)
#num_conf=len(times)
input_files = glob.glob("config*.dat")

input_files.sort(key=lambda f: int(filter(str.isdigit, f)))
times=cf.extract_digits(input_files)
num_conf=len(times)

Box,L=pa.Box_limits(input_files[0])






#num_cores = multiprocessing.cpu_count()
#
#results=Parallel(n_jobs=num_cores,verbose=10)(delayed(compute_one_configuration)(fil) for fil in input_files)


results=[]

for fil in input_files:
    results.append(compute_one_configuration(fil))



g_r=np.average(results,axis=0)



#g_poly_solvent=g_poly_solvent/num_conf
#
import matplotlib.pyplot as plt

plt.plot(g_r[0][:,0],g_r[0][:,2],label="poly-solvent")
#plt.plot(g_r[1][:,0],g_r[1][:,2],label="poly-solute")
#plt.plot(g_r[2][:,0],g_r[2][:,2],label="poly-poly")
#
#plt.legend()
#
#


Data=pd.read_csv("gr.dat",sep=" ",skiprows=4,dtype=np.float64,header=None).values

plt.plot(Data[:,1],Data[:,2],label="lammps")
plt.legend()
plt.show()
