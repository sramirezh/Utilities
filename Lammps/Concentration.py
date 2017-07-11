"""
This script analyzes chunk-averaged data generated by LAMMPS
The chunk files should be SPlitted before using Chunk_Splitter.sh
In order to perform the Concentration analysis, the chunk extensions, should be as follow:

.chunk for all atoms properties.
.chunks for solute properties. 

It generates a file "Sconcentration.dat" that has the averages of all data.

"""

import numpy as np

#Reading the times to make it easier to read the file by chunks
Times=np.loadtxt("Times.dat",dtype=int)
x=np.size(Times)

#Getting the shape of the data array 
File_Name=str(int(Times[0]))+".chunk"
n,m=np.shape(np.loadtxt(File_Name,skiprows=1)) 


"""
Computing the averages and other parameters
"""
Concentration=np.zeros((n,2))
for k in xrange(x): #Runs over the sampled times.
   # print("Reading configuration %d of %d" %(k,x-1))
    File_Name=str(int(Times[k]))+".chunk"
    SFile_Name=str(int(Times[k]))+".chunks"
    Chunk_Results=np.loadtxt(File_Name,skiprows=1) 
    SChunk_Results=np.loadtxt(SFile_Name,skiprows=1) 
    
    for i in xrange(n):  
        if Chunk_Results[i,2]==0: continue
        Concentration[i,1]=Concentration[i,1]+SChunk_Results[i,2]/Chunk_Results[i,2]
        
Concentration[:,0]=Chunk_Results[:,1]
Concentration[:,1]=Concentration[:,1]/(k+1)
"""
Creating the output file
"""
np.savetxt("SConcentration.dat",Concentration)


##For Testing Porpuses
#import matplotlib.pyplot as plt
#plt.plot(Concentration[:,0],Concentration[:,1],'*')
#y=np.linspace(-0.1,1.1)
#x=np.zeros(len(y))
#x[:]=-2.11856
#plt.plot(x,y)