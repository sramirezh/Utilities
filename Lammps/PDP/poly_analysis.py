#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 17:55:31 2018

Analysing the polymer dynamics

Args:
    The poly.atom file should have the next structure in the first 8 columns,
    id type x y z ix iy iz
@author: simon
"""

import numpy as np 
from subprocess import Popen,PIPE
from shlex import split
import pandas as pd
#import argparse
import re
import linecache
import matplotlib.pyplot as plt


#parser = argparse.ArgumentParser(description='This script evaluates the trajectory file of a polymer')
#parser.add_argument('FileName', metavar='InputFile', type=str,help='Input filename')

#parser.add_argument('--min', help='Number of timesteps to be discarded', default=1000, type=int)


#args = parser.parse_args()
#InputFile=args.FileName

InputFile="1.5_f0.02.atom"
#Uncomment to split the trajectory again
#p1 = Popen(split('bash Trajectory_poly.sh -i'+InputFile),stdout=PIPE)
#out,err=p1.communicate()

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
        Array which has only 5 columns: id type x y z
    """
    array=Data[:,:5]
    for i in xrange(3):
        array[:,2+i]=Data[:,2+i]+Data[:,5+i]*L[i]
        
    return array

    
    


Box,L=Box_limits()



Times=pd.read_csv("Times.dat").as_matrix()
x=np.size(Times)

rel_tail=[]
rel_head=[]
for k in xrange(x): #Runs over the sampled times.
    print("Reading configuration %d of %d" %(k,x-1))
    File_Name=str(int(Times[k]))+".cxyz"
    # As there is a space after the las column, pandas read it as a column of nan, then we need to avoid it
    Data=pd.read_csv(File_Name,sep=" ",dtype=np.float64,header=None).as_matrix()[:,:-1]  
    n,m=Data.shape
    pos=real_position(Data) #Real positions of all the atoms
    v_cm=cm(pos[:,2::])
    
    #Evaluating the positions of the head and the tail respect to v_cm
    i_head=np.where((pos[:,0]==1))[0][0]
    i_tail=np.where((pos[:,0]==n))[0][0] 
    rel_tail.append(pos[i_tail,2]-v_cm[0])
    rel_head.append(pos[i_head,2]-v_cm[0])
    
    


#plt.plot(pos[:,2],pos[:,3],'o')
#plt.plot(v_cm[0],v_cm[1],'*')



    
    
plt.plot(Times[:,0],rel_tail)
#plt.plot(Times[:,0],rel_head,'r')

 


 

        
            
    
    
    
    