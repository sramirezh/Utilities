#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 16:23:44 2019

Reads the xyz taken from VMD and prepares a trajectory file that centers the cm of the polymer and leaves the new trajectory ready for a video
also analyses the distances between the nearest neighbors to the polymer monomers
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

def xyz_to_lammps(file_xyz,log_name):
    """
    This could be made of a function that reads the coordinates and then writes as lammps
    """
    
    return

   


def xyz_splitter(file_name):
    
    cf.bash_command("""%s/Lammps/Trajectory_Analysis/Trajectory_Splitter.sh -i %s""" %(utilities_path,file_name))
    
    return






parser = argparse.ArgumentParser(description='This script evaluates the trajectory file',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file_name', metavar='InputFile',help='Input filename',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-log', help='lammps logfile name to read the box limits', default="log.lammps", type=str)
parser.add_argument('-split', help='True if trajectory file need to be splitter', default=False, type=bool)

args = parser.parse_args()


if args.split==True:
    xyz_splitter(args.file_name)
else:
    print "The Trajectory file was not splitted"


files = glob.glob("*.cxyz")
files.sort(key=lambda f: int(filter(str.isdigit, f)))
file_name=files[3]
data=pd.read_csv(file_name,sep=" ",dtype=np.float64,skiprows=2,header=None).values

volume,limits=read_box_limits(args.log)
box_length=limits[1,:]-limits[0,:]

poly_coords=data[np.where(data[:,0]==3)[0],:]

#Getting the image 

#Cm in x with respect to the first monomer
n_monomers,m=np.shape(poly_coords)
dist=np.zeros(n_monomers)

n_bins=5
partition=np.linspace(limits[0,0],limits[1,0],n_bins+1)
xpos=poly_coords[:,1]
occupancy=np.zeros(n_bins)

for i in xrange(n_bins):
    indexes=np.where(xpos<partition[i+1])[0]
    occupancy[i]=len(indexes)
    xpos=np.delete(xpos,indexes)

if occupancy[0]!=0:
    ind=np.where(occupancy==0)[0]
    center_of_void=0


theta=pos_x/
    
    
    

