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

imin=0 #Number of configurations to skip
nbins=30
rmax=3


def computation(particles,p_types,dist,i,j,nbins, rmax):
    """
    computes the gr for particles type j around type i
    
    args:
        i,j particles types as in the trajectory file from LAMMMPS [1-solvent, 2-solute, 3-polymer]
        
    Returns:
        bin_count:
    """
    i=np.where(p_types == i)[0][0]
    j=np.where(p_types == j)[0][0]
    
    #indexes to delete
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
        vol_new=4/3*np.pi*rmax_bin**3
        bin_count[i,2]=bin_count[i,1]/(vol_new-vol_old)
        
    return bin_count




def compute_one_configuration(fil):
    

    Data=pd.read_csv(fil,sep=" ",skiprows=9,dtype=np.float64,header=None).values[:,:-1]
    n,m=Data.shape
    pos=pa.real_position(Data,L) #Real positions of all the atoms

    #Particle indexes
    p_types=np.unique(Data[:,1]).astype(int)
    particles=[np.where(Data[:,1]==j)[0] for j in p_types]


    dist=squareform(pdist(pos))
    
    
    g_poly_solvent=computation(particles,p_types, dist,3,1,nbins,rmax)
    return g_poly_solvent


"""
*******************************************************************************
Main
*******************************************************************************
"""

input_files = glob.glob("*.gz")

input_files.sort(key=lambda f: int(filter(str.isdigit, f)))

Box,L=pa.Box_limits(input_files[0])
times=cf.extract_digits(input_files)
num_conf=len(times)


#initialising

g_poly_solvent=np.zeros((nbins,3))

#I can do a function here and a loop inside tha calls computation for all the indexes that I want

num_cores = multiprocessing.cpu_count()

results=Parallel(n_jobs=num_cores,verbose=10)(delayed(compute_one_configuration)(fil) for fil in input_files)


#g_poly_solvent=g_poly_solvent/num_conf
# 
#import matplotlib.pyplot as plt
#
#plt.plot(g_poly_solvent[:,0],g_poly_solvent[:,2])
#plt.show()




