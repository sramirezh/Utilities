#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
script to be used after running Gatherer.sh


@author: sr802
"""
import numpy as np
import matplotlib.pyplot as plt

Times=np.loadtxt("../times.dat",dtype=int)
x=np.size(Times)

#Opening the first array
data=np.loadtxt("full-mean."+str(int(Times[0])))

for k in xrange(1,x): #Runs over the sampled times.
    print("Reading configuration %d of %d" %(k,x-1))
    #Reading the Results 
    name="full-mean."+str(int(Times[k]))
    array=np.loadtxt(name)
    data[:,2::]=data[:,2::]+array[:,2::]
    
#Normalizing
data[:,2::]=data[:,2::]/x
np.savetxt("Averages.dat",data)

#plt.plot(data[:,1],data[:,3])
#plt.plot(array[:,1],array[:,3])


    