"""
This script analyzes the trajectory files splitted with Trajectory_Splitter.sh.
The trajectory can have a variable Number of particles


It generates two files "Sconcentration.dat" and "FConcentration.dat" that has average the concentration for the Solutes and Solvents.

"""
from __future__ import division
import numpy as np




print "Remember to define the Volume and the length in X so the bins are the same as in Lammps"

xmin=-30
xmax=30
L=xmax-xmin
binS=0.05
Nbins=int(L/0.05)

print Nbins
#Reading the times to make it easier to read the file by chunks
Times=np.loadtxt("Times.dat",dtype=int)
x=np.size(Times)

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
    print File_Name
    Data=np.genfromtxt(File_Name,skip_header=2)
    n,m=Data.shape

    #Checking if the solutes go through the surface,
    for i in xrange(n):
        if Data[i,0]==1:
            Nf[np.minimum(int(np.floor((Data[i,1]-xmin)/delta)),l-1)]+=1 #The -xmin is to avoid negative indexes

        if Data[i,0]==2:
            Ns[np.minimum(int(np.floor((Data[i,1]-xmin)/delta)),l-1)]+=1


Ly=20 #30
Lz=20 #37.02016

Volume=delta*Ly*Lz
Ns=Ns/x/Volume
Nf=Nf/x/Volume
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
