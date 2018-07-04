#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 17:02:55 2018
This script supposes that the header is commented with "#" and that the last line of the header contines the data variable names.
@author: sr802
"""

from __future__ import division
import numpy as np
import pandas as pd
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf


input_file="prof_sol.dat.tmp"

# This is to remove all the headers and just start reading in the columns if there are

header_lines=0
last_pound_pos=-1
with open(input_file, 'r') as data_file:
    while(data_file.read(1)=='#'):
        last_pound_pos = data_file.tell()
        header=data_file.readline()
        header_lines+=1

    #Read the next lines
    data_1=data_file.readline().split()
    data_2=data_file.readline().split()

    data_file.seek(last_pound_pos+1) #Goes back to the last line of the header
    if len(data_1)!=len(data_2): #If there is a line containing the number of particles, 
        data_file.readline()
    data_file.readline()
    data=pd.read_csv(data_file,sep=" ",header=None).dropna(axis=1,how='all') 

if header_lines>0:
    data.columns=header.split()
    
        
        
    
    
#   data_file.seek(last_pound_pos+1) #The 1 is added to avoid reading a non-existing column in the names as there lines are "# " Seek goes to the byte of the file
#
#    data_file.seek(last_pound_pos+1)
#    data=pd.read_csv(data_file,sep=" ",header=0)