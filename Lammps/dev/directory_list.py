#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 18:56:26 2018
Read the files in a directory
@author: sr802
"""

import os
import re
import glob
import numpy as np

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

#onlyfiles = [f for f in os.listdir(cwd) if os.isfile(os.join(cwd, f))]

numbers=[]
os.chdir(cwd)
files=glob.glob('*.gz')
for f in files:
    numbers.append(re.findall(r"[-+]?\d*\.?\d+",f))
    
numbers=np.sort(np.array(numbers,dtype=float).reshape((len(numbers))))


