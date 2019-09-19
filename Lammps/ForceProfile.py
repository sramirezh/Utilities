"""
Script to produce the average force per particle in function of the distance from the wall using Yawei's approach
It gives Fig 3.c
"""

import numpy as np
import matplotlib.pyplot as plt

"""Input parameters"""
rhoS=0.02
rhoF=0.74
T=0.846091996
GrhoS=1

"""Loading the properties"""
F=np.loadtxt("A_dir/full-mean") #Solvent properties
S=np.loadtxt("B_dir/full-mean") #Solute properties
All=np.loadtxt("stress_dir/full-mean") #All properties

"""Computing the forces"""
fS=-T/rhoS*GrhoS
fF=-fS*rhoS/rhoF
              
n,m=np.shape(F)
Force=np.zeros((n,2))
Force[:,0]=F[:,1]
for i in range(n):
    if All[i,5]==0:
        Force[i,1]=0
    else:
        Force[i,1]=(fF*F[i,4]+fS*S[i,4])/All[i,5]
    
np.savetxt("Force.dat",Force)    
