"""
Script to produce the average force per particle in function of the distance from the wall with Hiroaki approach and Yawei approach

Many parameters have to be set.

Generates YForce.dat and HForce.dat
"""

import numpy as np
import matplotlib.pyplot as plt

"""Input parameters"""
#Taking them from the log.out
print "\nBe sure that the input parameters are set properly\n"

#Hiroaki data
Ns=391.576673
Nb=2126.836227 #Number of fluid particles in the bulk
GradChem=-0.125

#Yawei data
rhoS=0.135973
rhoF=0.602 #Estimation from the plot 
T=1

"""Loading the properties"""
f=np.loadtxt("LAverages.dat") #Solvent properties
s=np.loadtxt("SAverages.dat") #Solute properties
All=np.loadtxt("AAverages.dat") #All properties

              
"""Computing Yawei forces"""
FsY=-GradChem
FfY=-FsY*rhoS/rhoF

"""Computing Hiroaki Forces"""
#Do not include the GradChem
FsH=(Nb-Ns)/Nb
FfH=-FsH*(Ns/(Nb-Ns))
              
n,m=np.shape(f)
YForce=np.zeros((n,2))
YForce[:,0]=f[:,1]
HForce=np.zeros((n,2))
HForce[:,0]=f[:,1]
for i in xrange(n):
    if All[i,4]==0:
        YForce[i,1]=0
        HForce[i,1]=0
        
    else:
        YForce[i,1]=(FfY*f[i,4]+FsY*s[i,4])/All[i,4]
        HForce[i,1]=(FfH*f[i,4]+FsH*s[i,4])/All[i,4]
    
plt.plot(YForce[:,0],YForce[:,1])
plt.plot(HForce[:,0],HForce[:,1])
np.savetxt("YForce.dat",YForce)
np.savetxt("HForce.dat",HForce)    


#==============================================================================
# Generatig the Files to be iterated in Lammps
#==============================================================================

zBulk=8
Index=np.where(HForce[:,0]<zBulk)
BinSize=HForce[1,0]-HForce[0,0]

Zpos=np.append(HForce[Index,0],HForce[Index[0][-1],0]+BinSize)
HF=np.transpose(HForce[Index,1])

np.savetxt("Zpos_iterate.dat",Zpos)
np.savetxt("YForce_iterate.dat",HF)





