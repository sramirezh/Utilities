"""
This script computes the solute concentration=Nsolutes/(Nsolute+Nsolvent)
Before you need to run Concentration.sh
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

for k in xrange(1): #Runs over the sampled times.
   # print("Reading configuration %d of %d" %(k,x-1))
    File_Name=str(int(Times[k]))+".chunk"
    SFile_Name=str(int(Times[k]))+".chunks"
    Chunk_Results=np.loadtxt(File_Name,skiprows=1)  
    SChunk_Results=np.loadtxt(SFile_Name,skiprows=1)  

Nsolutes=np.sum(SChunk_Results[:,2])
Ntotal=np.sum(Chunk_Results[:,2])

print ("The solute concentration is %f"%(Nsolutes/Ntotal))