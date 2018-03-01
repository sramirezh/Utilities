"""
This script reads a final configuration from Lammps and performs widom insertion method of two different species 
in two different reservoirs to compute the chemical potential.
"""
from __future__ import division
import numpy as np


def PotentialEnergy(i):
    """
    Evaluates the potential energy of the particle i interacting with all the other particles.
    """
    V=0
    for j in xrange(i+1,n):
        
        dx=Data[i,1]-Data[j,1] 
        dy=Data[i,2]-Data[j,2]
        dz=Data[i,3]-Data[j,3]
        
        
        
        #Apply Boundary conditions
        dx-=Lx*np.rint(dx/Lx)
        dy-=Lx*np.rint(dx/Lx)
        dz-=Lx*np.rint(dx/Lx)
        
        R2=dx**2+dy**2+dz**2
        if R2<=Rc2:
            V+=LJ(R2,Sigma,Epsilon)-LJcut

    return V             
    
    
    
def LJ(r2,Sigma,Epsilon):
    s2r=Sigma**2/r2
    
    V=4*Epsilon*((s2r)**6-(s2r)**3.)
    
    return V

def InsertTrial():
    Vinsertion=0
    np.random.rand(2)
    
    return Vinsertion



Data=np.loadtxt("trajectory.xyz", skiprows=2)
n,m=np.shape(Data)

#Box limits
Box=[-30,30,-10,10,-10,10]
Lx=Box[1]-Box[0]
Lz=Box[3]-Box[2]
Ly=Box[5]-Box[4]

#Input data
Sigma=1.0
Rc=(2)**(1/6)
Rc2=Rc**2
Epsilon=1.0
Temperature=1
LJcut=LJ(Rc**2,Sigma,Epsilon)


Pe=0
for i in xrange(n):
    print i
    Pe+=PotentialEnergy(i)
print Pe  
    






