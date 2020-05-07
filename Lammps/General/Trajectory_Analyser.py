"""
This script analyzes the trajectory files splitted with Trajectory_Splitter.sh.
The trajectory can have a variable Number of particles

The analysis is performed for a given volume otherwise in the simulation box

It generates two files "Sconcentration.dat" and "FConcentration.dat" that has average the concentration for the Solutes and Solvents

It prints the force factor for pressure driven simulations and also other geometrical parameters.


In My particle definition
1=Solvent
2=Solutes
3=Lower Solid Wall
4=Upper Solid Wall

"""

import numpy as np
import argparse
import os
import sys
import linecache
import re
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.lammps_utilities as lu
"""
###############################################################################
Function definitions
###############################################################################
"""

def check_analysis_box(cv_limits,simulation_limits):
    """
    checks if the analysis box is consistent with the simulation box
    Args:
        cv_limits control box limits
        simulation_limits simulation box limits
    Returns:
        cv_limit after checking that its valid
    """
    vol_args=np.size(cv_limits)
    #Checking if the input is valid
    if vol_args!=6 and (not args.Volume)==False: #If the argument number is not 6 or 0
        sys.exit("There is a mistake with the analysis volume, 6 inputs required")
    elif (not args.Volume): #If no arguments, sets the control volume as the simulation box.
        print("\nNo volume given, assumend the entire box")
        cv_limits=simulation_limits
    else: #If 6 arguments given, check if the volume is within the simulation box.
        cv_limits=np.transpose(np.reshape(np.array(cv_limits),(3,2)))
        for i in range(3):
            for j in range(2):
                #if (cv_limits[j,i]<simulation_limits[0,i] or cv_limits[j,i]>simulation_limits[1,i]):
                if (simulation_limits[0,i]<=cv_limits[j,i] and simulation_limits[1,i]>=cv_limits[j,i]):
                    continue
                else:
                    sys.exit("There is an error with the limit of the analysis box in the %d-th direction"%i)
                    break

    return cv_limits

def read_times():
    """
    Read the Times.dat file created with the trajectory splitter
    Retuns:
        x number of timesteps
        Times an array with all the timesteps

    """
    times=np.sort(pd.read_csv("Times.dat",header=None).values,axis=0)
    x=np.size(times)
    times=np.reshape(times,(x,1))
    return x,times


def one_dim_slicer(data,limits,index):
    """
    returns the particles that have the positions in the given direction within the cv_limits
    Args:
        data has the paricle type, posx,poy, posz
        limits: array containing min and max in the direction of analysis [min,max]
        index indicates the index for the given direction
    """
    indexes=np.where(np.logical_and(data[:,index]>=limits[0],data[:,index]<=limits[1]))[0]
    s_data=data[indexes,:]
    return s_data

def particles_cv(data,cv_limits,box_limits):
    """
    returns the particles inside the cv
    Args:
        data has the paricle type, posx,poy, posz
        cv_limits: array containing min and max each direction for the control box
        box_limits: array containing min and max each direction for the simulatio box
    """
    if np.array_equal(cv_limits,box_limits)==False:
        for i in range(3):
            if np.array_equal(cv_limits[:,i],box_limits[:,i])==True:
                continue
            else:
                data=one_dim_slicer(data,cv_limits[:,i],i+1)
    return data

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


"""
###############################################################################
Argument control
###############################################################################
"""

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script


parser = argparse.ArgumentParser(description='This script evaluates the trajectory file')
parser.add_argument('FileName', metavar='InputFile',help='Input filename',type=lambda x: is_valid_file(parser, x))
parser.add_argument('--Volume', help='Insert the limits of the Analysis box, if not the total volume is assumed', nargs='+', type=float)
parser.add_argument('--bin', help='the bin size for the analysis', default=0.05, type=float)
parser.add_argument('--log', help='lammps logfile name', default="log.lammps", type=str)
parser.add_argument('--split', help='True if trajectory file need to be splitter', default=False, type=bool)

args = parser.parse_args()
bin_size=args.bin
cv_limits=args.Volume
log_name=args.log
InputFile=args.FileName

if args.split==True:
    bash_command("""%s/Trajectory_Splitter.sh -i %s"""  %(dir_path,InputFile))
else:
    print("The Trajectory file was not splitted")




"""
###############################################################################
Control Volume, or volume where the analysis is going to be performed
###############################################################################
"""

#For the simulation box
box_volume,box_limits=lu.read_box_limits(log_name)

#For the analysis box
cv_limits=check_analysis_box(cv_limits,box_limits)
cv_length=np.diff(cv_limits,axis=0)[0]
cv_volume=np.prod(cv_length)
number_bins=int(cv_length[0]/bin_size)

n_tsteps,times=read_times()

bin_limits=np.linspace(cv_limits[0,0],cv_limits[1,0],number_bins+1)

print("Delta is %lf"%bin_size)

bin_centers=bin_limits[:-1]+0.5*bin_size




"""
###############################################################################
Computing the averages and other parameters
###############################################################################
"""

Ns=np.zeros(number_bins)
Nf=np.zeros(number_bins)

for k in range(n_tsteps): #Runs over the sampled times.
    print(("Reading configuration %d of %d" %(k,n_tsteps-1)))
    file_name=str(int(times[k]))+".cxyz"
    data=pd.read_csv(file_name,sep=" ",dtype=np.float64,skiprows=2,header=None).values

    """
    Getting the position of the surface
    """
    if k==0:
        if np.max(data[:,0])<=2:
            print("There is no solid surface")
        else:
            print("Analysing the solid surface")
            lu.solid_surface(data,3)

    data=particles_cv(data,cv_limits,box_limits)
    n,m=data.shape

    for i in range(n):
        if data[i,0]==1:
            errorv=np.minimum(int(np.floor((data[i,1]-cv_limits[0,0])/bin_size)),number_bins)
            Nf[np.minimum(int(np.floor((data[i,1]-cv_limits[0,0])/bin_size)),number_bins-1)]+=1 #The -xmin is to avoid negative indexes

        if data[i,0]==2:
            Ns[np.minimum(int(np.floor((data[i,1]-cv_limits[0,0])/bin_size)),number_bins-1)]+=1
        ifailure=i+1

np.array_equal(cv_limits,box_limits)


Chunk_Volume=bin_size*cv_length[1]*cv_length[2]
num_fluid=int(np.sum(Nf)/n_tsteps)
Ns=Ns/n_tsteps/Chunk_Volume
Nf=Nf/n_tsteps/Chunk_Volume



print("The volume of the control volume is %lf, Number of fluid particles is %d and the Fp is %lf" %(cv_volume,num_fluid,cv_volume/num_fluid))

#Creating the output file
Ns=np.column_stack((bin_centers,Ns))
Nf=np.column_stack((bin_centers,Nf))
header=str(np.resize(np.transpose(cv_limits),6))
np.savetxt("SConcentration.dat",Ns,header=header)
np.savetxt("FConcentration.dat",Nf,header=header)
