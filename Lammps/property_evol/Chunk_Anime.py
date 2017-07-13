"""
This script analyzes chunk-averaged data generated by LAMMPS
The chunk files should be SPlitted before using Chunk_Splitter.sh

It generates a file "Averages.dat" that has the averages of all data and then if there are stress per atom, computes the 
pressure.

The positions of the stresses have to be given explicitly as CSV, and the index counting starts at zero in chunk, for example:
    ['Chunk', 'Coord1', 'Ncount', 'vx', 'vz', 'density/mass', 'c_Stress[3]'], c_Stress[3] is the 6th entry. 
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import time

print "\nRemember to define the surface shift and run only Chunk_Splitter!!!!\n"

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

"""
Functions
"""

def Parameter_Finder(List, String):
    """
    Finds a string on a List and returns the position on the list
    """
    cont=0
    indexes=[]
    for s in List:
        if String in s: 
            indexes.append(cont) 
        cont+=1
    return indexes 

"""
Reading the necessary data
"""

InputArg=sys.argv #Aditional arguments

#Confirmation details
if len(InputArg)==1:
    pressFlag=0
    print "\nNo stress calculations explicitly required, searching for the word 'stress' in any of the Parameters...\n"
else:
    pressFlag=1
    Stresses=map(int,sys.argv[1].split(","))	
    print "The positions of the stress inputs are:\n" 
    for entry in Stresses:
	print entry 

#Reading the header
f=open("header", 'r')
Parameters=f.readlines()[2].split() 
f.close()
Nparam=len(Parameters)-1 
ExcludeP=2 #Excluded parameters, like the position of the chunk
Parameters.remove("#")
print "The parameters are:\n"
for f in Parameters:
    print f

#Checking if there are stress calculations
IsPress=np.zeros(Nparam)
if pressFlag==0: 
    Stresses=Parameter_Finder(Parameters,"tress")
    if len(Stresses)>0:
        pressFlag=1
        print "Found %d possible stress parameter in the input:" %(len(Stresses))
        for elem in Stresses:
            print Parameters[elem]
IsPress[Stresses]=1    
index=Parameter_Finder(Parameters,"density/mass")[0]



#Reading the times to make it easier to read the file by chunks
Times=np.loadtxt("Times.dat",dtype=int)
x=np.size(Times)

#Getting the shape of the data array 
File_Name=str(int(Times[0]))+".chunk"
n,m=np.shape(np.loadtxt(File_Name,skiprows=1)) 
   

"""
Computing the averages and other parameters
"""

Averages=np.zeros((n,Nparam))
IsPressure=0
for k in xrange(x): #Runs over the sampled times.
   # print("Reading configuration %d of %d" %(k,x-1))
    File_Name=str(int(Times[k]))+".chunk"
    Chunk_Results=np.loadtxt(File_Name,skiprows=1) 
    plt.plot(Chunk_Results[:,1],Chunk_Results[:,3])
    plt.ylim(-0.3,0.1)
    plt.xlim(zshift,25)
    plt.pause(0.01)
    plt.clf()
    for l in xrange(ExcludeP,Nparam): #Runs over the parameter
        if IsPress[l]==1:
            Averages[:,l]=Averages[:,l]-Chunk_Results[:,l]*Chunk_Results[:,index] #Stress per atom*mass/density
        else:
            Averages[:,l]=Averages[:,l]+Chunk_Results[:,l]
        
    #Averages[:,1]=Averages[:,1]+Chunk_Results[:,5]*Chunk_Results[:,6]  #To set the Pyy    

#Building the final array
for i in xrange(Nparam):
    if i<ExcludeP:
        Averages[:,i]=Chunk_Results[:,i]
    else:
        Averages[:,i]=Averages[:,i]/(k+1) 



Averages[:,1]=Averages[:,1]-zshift #Adding the shift due the surface
"""
Creating the output file
"""
np.savetxt("Averages.dat",Averages)


#For Testing Porpuses
#import matplotlib.pyplot as plt
#
#plt.plot(Averages[:,1],Averages[:,3])
#plt.xlim([0,25])
