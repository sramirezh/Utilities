# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 01:46:04 2015

@author: simon
"""

import numpy as np
error=np.loadtxt("jkerror")
CorError=np.loadtxt("corError")
RelErr=np.zeros((5,5))
for i in range(5):
    for j in range(5):
        if j!=i:
            RelErr[i,j]=np.abs(error[i,j]-CorError[i,j])/(0.5*(error[i,j]+CorError[i,j]))*100.0