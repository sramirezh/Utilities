#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 17:05:14 2018
Trying to paralelise things in python

"""

from joblib import Parallel, delayed
import multiprocessing
import itertools


# what are your inputs, and what operation do you want to 
# perform on each input. For example...
inputs = range(100) 
jnputs =range(100)

def processInput(i):
    return i**2
 
num_cores = multiprocessing.cpu_count()

pairs=[]
for i in inputs:
    pairs.append([i,jnputs[i]])
results=Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)



#Trying to build an entrance pair

a = ['foo', 'bar', 'baz']
b = ['x', 'y', 'z', 'w']
input_param=[]
constant=3
for r in itertools.product(a, b): input_param.append([r[0],r[1],constant])