"""
This script analyzes the trajectory files splitted with Trajectory_Splitter.sh.
The trajectory can have a variable Number of particles

It generates two files "Sconcentration.dat" and "FConcentration.dat" that has average the concentration for the Solutes and Solvents

It prints the force factor for pressure driven simulations and also other geometrical parameters.


In My particle definition
1=Solvent
2=Solutes
3=Lower Solid Wall
4=Upper Solid Wall



"""
from __future__ import division
import numpy as np
import argparse
import os
import sys
import linecache
import re

sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
from Lammps.linux import bash_command

"""
Function definitions
"""

def solid_surface(atom_type):
    """
    Computes the limits of the solid surface
    Args:
        atom_type atom type in the trajectory file of the solid surface.
    Returns:
        Characteristics of the solid surface
        Writes a file Zshift.dat with the maximun height to be used by other codes. 
    """
    
    #Getting the maximum position of the surface.
    Maxz=-100
    Minz=10000
    Nfluid=0


    for i in xrange(n):
        if Data[i,0]==3: #3 is for solid surface, 2 for solutes, 1 for solvents.
            Maxz=max(Maxz,Data[i,atom_type])
            Minz=min(Minz,Data[i,atom_type])
        else:
            Nfluid+=1

    print "The maximum height of the solid surface is %lf" %Maxz
    print "The minimum height of the solid surface is %lf" %Minz
    print "The height of the solid surface is %lf" %(Maxz-Minz)

    #Writing the Zshift
    f=open("Zshift.dat",'w')
    f.writelines("%lf \n" %Maxz)
    f.close
    

def read_box_limits(log_name):
    """
    Reads the box limits from log.lammps
    Args:
        None: log_name name of the log file
    returns:
        volume
        limits 
    
    """
    out,err=bash_command("""grep -n "Created orthogonal" %s | awk -F":" '{print $1}' """%log_name)
    line=int(out.split()[0])
    limits=linecache.getline(log_name, line)
    limits=re.findall(r"[-+]?\d*\.?\d+", limits)
    limits=np.array(np.reshape(limits,(2,3)),dtype=float) #To have the box as with columns [x_i_min,x_i_max]
    volume=(limits[1,0]-limits[0,0])*(limits[1,1]-limits[0,1])*(limits[1,2]-limits[0,2])

    return volume,limits


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
        sys.exit("There is a mistake with the analysis volume, 6 required")
    elif (not args.Volume): #If no arguments, sets the control volume as the simulation box.
        print "\nNo volume given, assumend the entire box"
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
    

"""
Argument control
"""
parser = argparse.ArgumentParser(description='This script evaluates the trajectory file of a polymer')
#parser.add_argument('FileName', metavar='InputFile',help='Input filename',type=lambda x: is_valid_file(parser, x))
parser.add_argument('--Volume', help='Insert the limits of the Analysis box, if not the total volume is assumed', nargs='+', type=float)
parser.add_argument('--bin', help='the bin size for the analysis', default=0.05, type=float)
parser.add_argument('--log', help='lammps logfile name', default="log.lammps", type=str)

args = parser.parse_args()
binS=args.bin
cv_limits=args.Volume
log_name=args.log


"""
Analysis
"""

box_volume,box_limits=read_box_limits(log_name)
cv_limits=check_analysis_box(cv_limits,box_limits)

box_length=np.diff(cv_limits,axis=0)[0]

number_bins=int(box_length[0]/binS)

#Reading the times to make it easier to read the file by chunks
Times=np.loadtxt("Times.dat",dtype=int)
x=np.size(Times)
Times=np.reshape(Times,(x,1))

#Getting the shape of the data array
File_Name=str(int(Times[0]))+".cxyz"
Data=np.genfromtxt(File_Name,skip_header=2)

Xarray=np.linspace(xmin,xmax,Nbins)
delta=binS

print "Delta is %lf"%delta

CenterPos=Xarray[:-1]+0.5*binS
l=CenterPos.size
Ns=np.zeros(l)
Nf=np.zeros(l)
"""
Computing the averages and other parameters
"""

for k in xrange(x): #Runs over the sampled times.
    print("Reading configuration %d of %d" %(k,x-1))
    File_Name=str(int(Times[k]))+".cxyz"
    Data=np.genfromtxt(File_Name,skip_header=2)
    n,m=Data.shape

    """
    Getting the position of the surface
    """
    if k==0:
        if np.max(Data[:,0])<=2:
            print "There is no solid surface"
        else:
            print "Analysing the solid surface"
            solid_surface(3)

    for i in xrange(n):
        if Data[i,0]==1:
            Nf[np.minimum(int(np.floor((Data[i,1]-xmin)/delta)),l-1)]+=1 #The -xmin is to avoid negative indexes

        if Data[i,0]==2:
            Ns[np.minimum(int(np.floor((Data[i,1]-xmin)/delta)),l-1)]+=1


Ly=20 #30
Lz=20 #37.02016

Chunk_Volume=delta*Ly*Lz
Ns=Ns/x/Chunk_Volume
Nf=Nf/x/Chunk_Volume

#Getting the force factor for pressure driven calculations
Volume=Ly*Lz*L
Nfluid=(np.sum(Ns)+np.sum(Nf))*Chunk_Volume
Fp=Volume/(Nfluid)

print "The vol is %lf, Number of fluid particles is %lf and the Fp is %lf" %(Volume,Nfluid,Fp)

#Creating the output file
Ns=np.column_stack((CenterPos,Ns))
Nf=np.column_stack((CenterPos,Nf))

np.savetxt("SConcentration.dat",Ns)
np.savetxt("FConcentration.dat",Nf)
#
#
##For Testing Porpuses
#import matplotlib.pyplot as plt
#plt.plot(CenterPos,Ns[:,1],'*')
