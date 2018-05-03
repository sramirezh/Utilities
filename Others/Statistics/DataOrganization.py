# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 14:52:26 2015

@author: simon
"""
import numpy as np

data1=np.loadtxt("data1.txt",usecols=[1,2])
data2=np.loadtxt("data2.txt",usecols=[1,2])
data3=np.loadtxt("data3.txt")
data=np.concatenate((data1,data2,data3),axis=1)
data=np.delete(data,4,1)
np.savetxt("data",data)

variance=np.var(data,axis=0)
Error=np.sqrt(variance/25000)
print(Error)
