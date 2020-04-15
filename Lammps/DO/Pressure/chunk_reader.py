#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:01:03 2020
based on based.py, LAMMPS.py and LAMMPSparser.py from mdanalysis
@author: simon
"""
import numpy as np
import os
import six
import pandas as pd
from io import StringIO

# TODO add the value of the timestep
class timestep(object):
    def __init__(self, frame, n_chunks, data):
        """
        data is a dataframe
        """
        self.frame = frame
        self._n_chunks = n_chunks
        self.data = data
        self.header = []
                
    def __len__(self):
        return self.n_chunks

    def __iter__(self):
        """Iterate over chunks

        ``for chunk in timestep``

            iterate of the chunks
        """
        for i in range(self.n_chunks):
            yield self[i]
    

class chunk_reader(object):
    
    def __init__(self, filename):
        self.filename = filename
#        self._read_first_frame()     #An underscore (_) at the beginning is used to denote private variables in Python
#        self._reopen()
        self._get_initial()
        self._get_position_frames()
        self.frames = []
#        self._read_next_timestep()
#        root, ext = os.path.splitext(self.filename)
#        self._cache = {}
        
        
#    def _reopen(self):
#        self.close()
#        self._file = open(self.filename)
##        self.ts = timestep(self.n_chunks) # Creates an element of the class timestep
#        self.ts.frame = -1
    
    def close(self):
        if hasattr(self, '_file'):
            self._file.close()
    
    def _get_initial(self):
        """
        reads the header and estimates the number of chunks
        """
        header_lines = 0
        last_pound_pos=-1
        f = open( self.filename)
        while(f.read(1)=='#'):
              last_pound_pos = f.tell()
              header_lines += 1
              self.header = f.readline().split()
        self._first_line = f.tell()      
        self.n_chunks = int(f.readline().split()[1])  # Timestep Number-of-chunks Total-count

        f.close()
        
    
    def _get_position_frames(self):
        """
        Gets the starting position of each frame (Offsets in bytes), that can be used with seek
        Assumes that the number of chunks is constant
        """
        
        f = open( self.filename)
        f.seek(self._first_line) # start reading after header
        lines_per_frame = self.n_chunks + 1
        offsets = []
        counter = 0
        
        line = True
        while line:
            if not counter % lines_per_frame:
                offsets.append(int(f.tell()))
            line = f.readline()
            counter += 1
        array_offsets = offsets[:-1]  # last is EOF
        self.offsets = np.array(array_offsets, dtype = int)
        
        f.close()
        
        #TODO use something as in MDANALYSIS "read_frame", or read next frame, 
        #ir iter, such that when I am iterating I do not open and close the file
        
#    def get_all_data(self):
#        
#        self.frames.append(timestep())
#        f = open( self.filename)
#        f.seek(self._first_line) # start reading after header
#        lines_per_frame = self.n_chunks + 1
#        offsets = []
#        counter = 0
#        
#        line = True
#        while line:
#            if not counter % lines_per_frame:
#                offsets.append(int(f.tell()))
#            line = f.readline()
#            counter += 1
#        array_offsets = offsets[:-1]  # last is EOF
#        self._offsets = np.array(array_offsets, dtype = int)
#        
#        f.close()
        
        
# One way of using it
#results = chunk_reader("velocity_all.dat") 
#        
#    
## Creating all the frames
#
#f = open( results.filename)
#
#series = []
#
#byte_pos = results.offsets
#
#for i, start in enumerate(byte_pos[:-1]):
#
#    f.seek(start)
#    step, n_chunks, _ = f.readline().split() # Timestep Number-of-chunks Total-count
#    start = f.tell()
#    stuff = f.read(byte_pos[i+1] - start)
#    data = pd.read_csv(StringIO(stuff),sep=" ",header=None).dropna(axis=1,how='all')
#    data.columns = results.header
#    series.append(timestep(step, n_chunks, data))


