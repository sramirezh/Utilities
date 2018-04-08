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

def SolidSurface():
    #Getting the maximum position of the surface.
    Maxz=-100
    Minz=10000
    Maxf=0
    Nfluid=0


    for i in xrange(n):
        if Data[i,0]==3: #3 is for solid surface, 2 for solutes, 1 for solvents.
            Maxz=max(Maxz,Data[i,3])
            Minz=min(Minz,Data[i,3])
        else:
            Nfluid+=1

    print "The maximum height of the solid surface is %lf" %Maxz
    print "The minimum height of the solid surface is %lf" %Minz
    print "The height of the solid surface is %lf" %(Maxz-Minz)

    #Writing the Zshift
    f=open("Zshift.dat",'w')
    f.writelines("%lf \n" %Maxz)
    f.close



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
Times=np.reshape(Times,(x,1))

#Getting the shape of the data array
File_Name=str(int(Times[0]))+".cxyz"
Data=np.genfromtxt(File_Name,skip_header=0)

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
    Data=np.genfromtxt(File_Name,skip_header=0)
    n,m=Data.shape

    """
    Getting the position of the surface
    """
    if k==0:
        if np.max(Data[:,0])<=2:
            print "There is no solid surface"
        else:
            print "Analysing the solid surface"
            SolidSurface()

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
