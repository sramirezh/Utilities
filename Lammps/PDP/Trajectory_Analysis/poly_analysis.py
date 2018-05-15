#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 17:55:31 2018

Analysing the polymer dynamics, this file needs to be run in the directory of the trajectory file

Args:
    The poly.atom file should have the next structure in the first 8 columns,
    id type x y z ix iy iz
    Notice that this file is not sorted so always we need to see the indexes of the particles.

Returns:
    Writes files
    rdist_negative.dat and rdist_positive.dat containing the polymer density distributions from the center
    of mass.
    tip_behaviour.dat contains the head and tail behaviour with respect to the center of mass.



@author: simon
"""
from __future__ import division
import numpy as np
from subprocess import Popen,PIPE
from shlex import split
import pandas as pd
import argparse
import linecache
import os
import sys


sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
from Lammps.linux import bash_command



try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
except ImportError as err:
    print err

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


parser = argparse.ArgumentParser(description='This script evaluates the trajectory file of a polymer')
parser.add_argument('FileName', metavar='InputFile',help='Input filename',type=lambda x: is_valid_file(parser, x))
parser.add_argument('--split', help='True if trajectory file need to be splitter', default=False, type=bool)
parser.add_argument('--Nbins', help='Number of bins for the chunk calculations', default=60, type=int)

args = parser.parse_args()
InputFile=args.FileName


#InputFile="1.5_f0.02.atom"

#Uncomment to split the trajectory again
if args.split==True:
    InputFile=InputFile
    p1 = Popen(split('bash '+dir_path+'/Trajectory_poly.sh -i '+InputFile),stdout=PIPE)
    out,err=p1.communicate()
else:
    print "The Trajectory file was not splitted"

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

def number_of_monomers():
    """
    reads the number of particles
    Args:

    Returns:
        The number of monomers

    """
    out,err=bash_command('grep -n -m1 "NUMBER OF ATOMS" '+InputFile)
    line_number=int(out.split(":")[0])
    number_part=int(linecache.getline(InputFile, line_number+1))


    return number_part

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
        bin_count: Array with the positions of the centers of the bins, the count number and the density
    """
    delta=rmax/nbins
    bin_count=np.zeros((nbins,3))
    bin_count[:,0]=np.linspace(delta,rmax,num=nbins)-delta/2.
    vol_old=0
    for i in xrange(nbins):
        if np.size(pos_sphere)==0:break #To stop counting after there are no more particles
        rmax_bin=delta*(i+1)
        indexes=np.where(pos_sphere[:,0]<=rmax_bin)
        count=np.size(indexes[0])
        bin_count[i,1]=count
        vol_new=4/3*np.pi*rmax_bin**3
        bin_count[i,2]=count/(vol_new-vol_old)
        vol_old=vol_new
        pos_sphere=np.delete(pos_sphere,indexes,axis=0) #Deletes the particles taken into account

    return bin_count

def gyration_radious_squared(pos):
    """
    Computes the gyration radious, assuming all the particles have the same mass=1
    Args:
        Pos, vector with the positions with measured from the cm
    Returns: The scalar gyration radious squared
    """
    gr2=np.sum(np.average(np.square(pos),axis=0))


    return gr2







#Reading the initial data
Box,L=Box_limits()
times=pd.read_csv("Times.dat",header=None).as_matrix()
x=np.size(times)

rel_tail=[]
rel_head=[]

nbins=args.Nbins
rmax=number_of_monomers()/2 #Assumes the maximum radius is number_particles/2

av_rd_positive=np.zeros((nbins,2))
av_rd_negative=np.zeros((nbins,2))

r_gyration_2=[]
for k in xrange(x): #Runs over the sampled times.
    print("Reading configuration %d of %d" %(k,x-1))
    File_Name=str(int(times[k]))+".cxyz"
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
    i_front=np.where(np.abs(pos_sphere[:,2])<=np.pi/2.)[0]
    pos_semi_positive=pos_sphere[i_front,:]
    pos_semi_negative=np.delete(pos_sphere,i_front,axis=0)

    """Computing the polymeric distribution"""
    rd_positive=radial_distribution(nbins,rmax,pos_semi_positive)
    rd_negative=radial_distribution(nbins,rmax,pos_semi_negative)
    av_rd_positive+=rd_positive[:,1:]
    av_rd_negative+=rd_negative[:,1:]
    rd_negative[:,0]=rd_negative[:,0]*-1

    """Other properties"""
    r_gyration_2.append(gyration_radious_squared(pos_relative))





"""
###############################################################################
Other properties
###############################################################################
"""
r_gyration=np.sqrt(r_gyration_2)


av_rd_positive=av_rd_positive/x
av_rd_negative=av_rd_negative/x
rd_positive[:,1:]=av_rd_positive
rd_negative[:,1:]=av_rd_negative

#Need to multiply by two the density as I counted only on the semisphere.
rd_positive[:,2]*=2
rd_negative[:,2]*=2

""" Radial distribution output"""

np.savetxt('rdist_positive.dat',rd_positive)
np.savetxt('rdist_negative.dat',rd_negative)

plt.figure()
plt.plot(rd_positive[:,0],rd_positive[:,2],'b')
plt.plot(rd_negative[:,0],rd_negative[:,2],'r')
plt.grid()
plt.xlabel("$r-r_{cm}$")
plt.ylabel("Concentration")
plt.savefig('radial_distribution.pdf')
plt.close()

#"""Fourier Transformation"""
#
## Number of samplepoints
#N = x
## sample spacing
#T = Times[1]-Times[0]
#
#yf_tail = scipy.fftpack.fft(rel_tail)
#xf_tail = np.linspace(0.0, 1.0/(2.0*T), N/2)
#
#yf_head = scipy.fftpack.fft(rel_head)
#xf_head = np.linspace(0.0, 1.0/(2.0*T), N/2)




"""alternation of tail and head in the x axis with respect to"""
plt.figure()


plt.subplot(311)
plt.plot(times[:,0],rel_tail)
plt.grid()

plt.subplot(312)
plt.plot(times[:,0],rel_head,'r')
plt.ylabel("$x-x_{cm}$")
plt.grid()

plt.subplot(313)
plt.plot(times[:,0],rel_head,'r')
plt.plot(times[:,0],rel_tail)
plt.grid()


#plt.subplot(223)
#plt.plot(xf_tail, yf_tail[:N//2])
#plt.grid()
#
#plt.subplot(224)
#plt.plot(xf_head, yf_head[:N//2])
#plt.grid()

plt.xlabel("Time")
plt.savefig('tip_behaviour.pdf')
plt.close()

tip_behaviour=np.transpose(np.vstack([times[:,0],rel_tail,rel_head]))
np.savetxt('tip_behaviour.dat',tip_behaviour,header="time_step tail_position head_position")




#"""
#Testing the radial distribution function
#"""

#def spherical_integral(Data):
#    delta=Data[1,0]-Data[0,0]
#    r_max=Data[:,0]+delta/2
#    n,m=np.shape(Data)
#    vol_old=0
#    part_bin=[]
#    for i in xrange(n):
#        vol_new=4./3.*r_max[i]**3*np.pi
#        part_bin.append((vol_new-vol_old)*Data[i,2])
#        vol_old=vol_new
#    return part_bin
#
#
#distrib_pos=np.loadtxt("rdist_positive.dat")
#distrib_negative=np.loadtxt("rdist_negative.dat")
#
#part_bin=np.array(spherical_integral(distrib_pos))/2
#print np.sum(part_bin),np.sum(distrib_pos[:,1])
#part_bin_n=np.array(spherical_integral(distrib_negative))/2
#print np.sum(part_bin_n),np.sum(distrib_negative[:,1])
