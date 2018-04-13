#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 17:55:31 2018

Analysing the polymer dynamics, this file needs to be run in the directory of the trajectory file

Args:
    The poly.atom file should have the next structure in the first 8 columns,
    id type x y z ix iy iz
    Notice that this file is not sorted so always we need to see the indexes of the particles.
    
    
@author: simon
"""
from __future__ import division
import numpy as np 
from subprocess import Popen,PIPE
from shlex import split
import pandas as pd
import argparse
import re
import linecache
import matplotlib.pyplot as plt
import os

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg
        


        

        

parser = argparse.ArgumentParser(description='This script evaluates the trajectory file of a polymer')
parser.add_argument('FileName', metavar='InputFile',help='Input filename',type=lambda x: is_valid_file(parser, x))


args = parser.parse_args()
InputFile=args.FileName






#InputFile="1.5_f0.02.atom"

#Uncomment to split the trajectory again
InputFile=cwd+'/'+InputFile
p1 = Popen(split('bash '+dir_path+'/Trajectory_poly.sh -i '+InputFile),stdout=PIPE)
out,err=p1.communicate()

#
def Box_limits():
    """
    read the box limits
    Args:
        
    Returns:
        limits, an array with the box limits in every dimension.
        L, box size per dimension
    """
    command=str('grep -n -m1 "BOUNDS" '+InputFile)
    p2=Popen(split(command),stdout=PIPE)
    out2,err2=p2.communicate()
    LineNumber=int(out2.split(":")[0])
    
    #a=file.readlines()[LineNumber,LineNumber+3]
    limits=[]
    for i in range(LineNumber,LineNumber+3):
        limits.append(linecache.getline(InputFile, i+1).strip('\n').split())
    limits=np.array(limits,dtype='double')
    L=limits[:,1]-limits[:,0]
    
    return limits,L

cm=lambda x: np.average(x,axis=0)

def real_position(Data):
    """
    return the real positions of the particles after taking into account the box image of the atom.
    Args:
        Data should have the next structure in the first 8 columns,
        id type x y z ix iy iz, where ix represents the image box of the particle as described in
        Lammps dum
    Returns:
        Array which has only 3 columns:  x y z
    """
    n,m=np.shape(Data)
    array=np.zeros((n,3))
    for i in xrange(3):
        array[:,i]=Data[:,2+i]+Data[:,5+i]*L[i]
        
    return array

def spherical_coordinates(X):
    """
    Converts from cartesian to spherical coordinates
    Args:
        Vector with positions [xi,yi,zi]
    Return:
        Vector [r,tetha,phi]
    """
    n,m=np.shape(X)
    new_X=np.zeros((n,3))
    xy=X[:,0]**2+X[:,1]**2
    new_X[:,0]=np.sqrt(xy+X[:,2]**2)
    new_X[:,1]=np.arctan2(np.sqrt(xy),X[:,2])
    new_X[:,2]=np.arctan2(X[:,1],X[:,0])
    
    return new_X

def relative_position(pos):
    """
    Computes the positions with respect to the center of mass of the system
    Args:
        pos is a vector with [xi,yi,zi]
    """
    v_cm=cm(pos)
    pos_rel=np.zeros((n,3))
    for i in range(3):
        pos_rel[:,i]=pos[:,i]-v_cm[i]
    return pos_rel

def radial_distribution(nbins,rmax,pos_sphere):
    """
    Computes the radial distribution of a centered point distribution in spherical coordinates
    args:
        nbins: Number of bins in the radial direction
        rmax: Maximum radious of analysis
        pos_sphere: vector with r,tetha, fhi
    Returns:
        bin_count: Array with the positions of the centers of the bins and the count number. 
    """
    delta=rmax/nbins
    bin_count=np.zeros((nbins,2))
    bin_count[:,0]=np.linspace(delta,rmax,num=nbins)-delta/2.
    for i in xrange(nbins):
        if np.size(pos_sphere)==0:break #To stop counting after the 
        rmax_bin=delta*(i+1)
        indexes=np.where(pos_sphere[:,0]<=rmax_bin)
        bin_count[i,1]=np.size(indexes[0])
        pos_sphere=np.delete(pos_sphere,indexes,axis=0)
        
    return bin_count
    

#Reading the initial data
Box,L=Box_limits()
Times=pd.read_csv("Times.dat").as_matrix()
x=np.size(Times)

rel_tail=[]
rel_head=[]

nbins=30
rmax=15

av_rd_positive=np.zeros(nbins)
av_rd_negative=np.zeros(nbins)

for k in xrange(x): #Runs over the sampled times.
    print("Reading configuration %d of %d" %(k,x-1))
    File_Name=str(int(Times[k]))+".cxyz"
    # As there is a space after the las column, pandas read it as a column of nan, then we need to avoid it
    Data=pd.read_csv(File_Name,sep=" ",dtype=np.float64,header=None).as_matrix()[:,:-1]  
    n,m=Data.shape
    pos=real_position(Data) #Real positions of all the atoms
    pos_relative=relative_position(pos) #
    #Evaluating the positions of the head and the tail respect to v_cm
    i_head=np.where((Data[:,0]==1))[0][0]
    i_tail=np.where((Data[:,0]==n))[0][0] 
    rel_tail.append(pos_relative[i_tail,0])
    rel_head.append(pos_relative[i_head,0])
    
    "Getting the points in front and in the back"
    #First i get the indexes of the points in front and then I delete those indexes to get the points in the back
    pos_sphere=spherical_coordinates(pos_relative) 
    i_front=np.where(np.abs(pos_sphere[:,2])<np.pi/2.)[0]
    pos_semi_positive=pos_sphere[i_front,:]
    pos_semi_negative=np.delete(pos_sphere,i_front,axis=0)
    
    """Computing the polymeric distribution"""
    rd_positive=radial_distribution(nbins,rmax,pos_semi_positive)
    rd_negative=radial_distribution(nbins,rmax,pos_semi_negative)
    av_rd_positive+=rd_positive[:,1]
    av_rd_negative+=rd_negative[:,1]
    rd_negative[:,0]=rd_negative[:,0]*-1
    
av_rd_positive=av_rd_positive/x
av_rd_negative=av_rd_negative/x


#To plot the radial distribution
plt.figure()
plt.plot(rd_positive[:,0],av_rd_positive,'b')  
plt.plot(rd_negative[:,0],av_rd_negative,'r')
plt.grid()
plt.savefig('radial_distribution.png')
plt.close()
#To plot the alternation of tail and head in the x axis with respect to    
plt.figure() 
plt.plot(Times[:,0],rel_tail)
plt.plot(Times[:,0],rel_head,'r')
plt.grid()
plt.savefig('tip_behaviour.png')
plt.close()

#plt.plot(pos[:,2],pos[:,3],'o')
#plt.plot(v_cm[0],v_cm[1],'*')


















    

            


    





    
 
    


 


 

        
            
    
    
    
    