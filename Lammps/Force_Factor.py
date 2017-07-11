"""
This script analyzes chunk-averaged data generated by LAMMPS
The chunk files should be SPlitted before using Concentration.sh
In order to perform the Concentration analysis, the chunk extensions, should be as follow:

.chunks for solute properties. 
.chunk for solute+solvent properties.

It prints the value of the force factor

"""

import numpy as np


#Reading the times to make it easier to read the file by chunks
Times=np.loadtxt("Times.dat",dtype=int)
x=np.size(Times)

#Getting the shape of the data array 
File_Name=str(int(Times[0]))+".chunk"
FirstChunk=np.loadtxt(File_Name,skiprows=1)
n,m=np.shape(FirstChunk) 

Averages=np.zeros((n,2))

print "\nRemember to set the zdrift properly\n"

#zshifts 
#0.5, 0.8   is -2.118560 
#1.0, 1.0   is -2.149910 
#1.5, 1.5   is -1.714010 

#SC
#0.5, 0.8   is -2.109250
#1.0, 1.0   is -2.089910
#1.5, 1.5   is -1.918170

#1 Layer
#0.5, 0.8   is -1.064190
#1.0, 1.0   is -1.100260
#1.5, 1.5   is -0.833002


zshift=-1.714010 
Averages[:,0]=FirstChunk[:,1]-zshift
"""
Computing the averages and other parameters
"""
BulkSolutes=[]
BulkTotal=[]

#These parameters define the bulk region
zmin=15
zmax=25

for k in xrange(x): #Runs over the sampled times.
   # print("Reading configuration %d of %d" %(k,x-1))
    SFile_Name=str(int(Times[k]))+".chunks"
    File_Name=str(int(Times[k]))+".chunk"
    
    Chunk_Results=np.loadtxt(File_Name,skiprows=1) 
    SChunk_Results=np.loadtxt(SFile_Name,skiprows=1)
    BSvector=[]
    BTvector=[]
    for i in xrange(n):          
        #Computing the bulk concentration as defined by Bocquet.
        if Averages[i,0]>=zmin and Averages[i,0]<=zmax:
            BSvector.append(SChunk_Results[i,2])
            BTvector.append(Chunk_Results[i,2])
            
    BulkSolutes.append(np.sum(BSvector))
    BulkTotal.append(np.sum(BTvector))
Ns=np.average(BulkSolutes)
Nb=np.average(BulkTotal)
Ffactor=Nb/(Nb-Ns)

print "The average number of solutes in the bulk is Ns %f and total particles in the bulk Nb %f" %(Ns,Nb)
print "The force factor is %f" %Ffactor
                       
#"""
#Creating the output file
#"""
#np.savetxt("SConcentration.dat",Concentration)
#
#
#"""
#Uncomment for testing
#"""
# 
#import matplotlib.pyplot as plt
#plt.figure(1)
#plt.plot(BulkSolutes)
#plt.plot(BulkTotal)

#x=np.linspace(min(Concentration[:,0]),max(Concentration[:,0]))
#y=np.zeros(len(x))
#y[:]=BulkC
#plt.plot(x,y)
#plt.xlim([0,25])


#plt.figure(2)
#
#plt.plot(Concentration[:,0],Integrand)
#x=np.linspace(min(Concentration[:,0]),max(Concentration[:,0]))
#y=np.zeros(len(x))
#plt.plot(x,y)
#plt.xlim([0,8])