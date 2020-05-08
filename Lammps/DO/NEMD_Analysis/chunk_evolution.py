#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 17:27:17 2020
This script checks the chunk evolution
@author: simon
"""

import numpy as np
import pandas as pd
import argparse
import os
import sys
from io import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.General.chunk_utilities as cu
import Lammps.core_functions as cf
import matplotlib.pyplot as plt




# =============================================================================
# Reading data
# =============================================================================
results = cu.chunk_reader("properties.dat") 


# Creating all the frames

#f = open(results.filename)

series = []

byte_pos = results.offsets

#for i, start in enumerate(byte_pos[:-1]):
#    print ("Reading frame %s/%s"%(i,results.n_frames))
#    f.seek(start)
#    step, n_chunks, _ = f.readline().split() # Timestep Number-of-chunks Total-count
#    start = f.tell()
#    stuff = f.read(byte_pos[i+1] - start)
#    data = pd.read_csv(StringIO(stuff),sep=" ",header=None).dropna(axis=1,how='all')
#    data.columns = results.header
#    series.append(cu.timestep(step, n_chunks, data))
#    
#    data['vx'] get this one and append


# Get the first frame

