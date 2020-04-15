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

# TODO add the value of the timestep
class timestep(object):
    def __init__(self, n_chunks):
        self.frame = -1
        self._n_chunks = n_chunks
        self.data = []
                
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
        self._get_initial()
#        self._read_first_frame()     #An underscore (_) at the beginning is used to denote private variables in Python
        self._reopen()
#        self._read_next_timestep()
#        root, ext = os.path.splitext(self.filename)
#        self._cache = {}
        
    
    def _get_initial(self):
        """
        reads the header and estimates the number of chunks
        """
        header_lines = 0
        last_pound_pos=-1
        with open("velocity_all_steps.dat") as f:
            while(f.read(1)=='#'):
                  last_pound_pos = f.tell()
                  header_lines += 1
                  self.header = f.readline()
            self._first_line = f.tell()      
            self.n_chunks = int(f.readline().split()[1])  # Timestep Number-of-chunks Total-count

        f.close()
    
    def next(self):
        """Forward one step to next frame."""
        try:
            ts = self._read_next_timestep()
        except (EOFError, IOError):
            self.rewind()
            six.raise_from(StopIteration, None)
        else:
            for auxname in self.aux_list:
                ts = self._auxs[auxname].update_ts(ts)

            ts = self._apply_transformations(ts)

        return ts

    def __next__(self):
        """Forward one step to next frame when using the `next` builtin."""
        return self.next()

    def rewind(self):
        """Position at beginning of trajectory"""
        self._reopen()
        self.next()
        
    def _reopen(self):
        self.close()
        self._file = open(self.filename)
        self.ts = timestep(self.n_chunks) # Creates an element of the class timestep
        self.ts.frame = -1
        
    def close(self):
        if hasattr(self, '_file'):
            self._file.close()
            
    def _read_frame(self, frame):
        self._file.seek(self._offsets[frame])
        self.ts.frame = frame - 1  # gets +1'd in next

        return self._read_next_timestep()
    
    def _read_next_timestep(self):
        f = self._file
        ts = self.ts
        ts.frame += 1
        if ts.frame >= len(self):
            raise EOFError
            
        data = f.readline() # ITEM TIMESTEP
        ts.data = data
#        step_num = int(f.readline())
#        ts.data['step'] = step_num
#
#        f.readline() # ITEM NUMBER OF ATOMS
#        n_atoms = int(f.readline())
#        if n_atoms != self.n_atoms:
#            raise ValueError("Number of atoms in trajectory changed "
#                             "this is not suported in MDAnalysis")
#
#        triclinic = len(f.readline().split()) == 9  # ITEM BOX BOUNDS
#        if triclinic:
#            xlo, xhi, xy = map(float, f.readline().split())
#            ylo, yhi, xz = map(float, f.readline().split())
#            zlo, zhi, yz = map(float, f.readline().split())
#
#            box = np.zeros((3, 3), dtype=np.float64)
#            box[0] = xhi - xlo, 0.0, 0.0
#            box[1] = xy, yhi - ylo, 0.0
#            box[2] = xz, yz, zhi - zlo
#
#            xlen, ylen, zlen, alpha, beta, gamma = mdamath.triclinic_box(*box)
#        else:
#            xlo, xhi = map(float, f.readline().split())
#            ylo, yhi = map(float, f.readline().split())
#            zlo, zhi = map(float, f.readline().split())
#            xlen = xhi - xlo
#            ylen = yhi - ylo
#            zlen = zhi - zlo
#            alpha = beta = gamma = 90.
#        ts.dimensions = xlen, ylen, zlen, alpha, beta, gamma
#
#        indices = np.zeros(self.n_atoms, dtype=int)
#
#        f.readline()  # ITEM ATOMS etc
#        for i in range(self.n_atoms):
#            idx, _, xs, ys, zs = f.readline().split()
#
#            indices[i] = idx
#            ts.positions[i] = xs, ys, zs
#
#        order = np.argsort(indices)
#        ts.positions = ts.positions[order]
#        # by default coordinates are given in scaled format, undo that
#        ts.positions = distances.transform_StoR(ts.positions, ts.dimensions)

        return data

        


results = chunk_reader("velocity_all_steps.dat")





#lines_per_frame = n_chunks + 1
#offsets = []
#counter = 0
#def n_frames(self):
#with open("velocity_all_steps.dat") as f:
#    line = True
#    
#    while line:
#        if not counter % lines_per_frame:
#            offsets.append(int(f.tell()))
#        line = f.readline()
#        counter += 1
#array_offsets = offsets[:-1]  # last is EOF
#array_offsets = np.array(array_offsets, dtype = int)
#
#f = open("velocity_all_steps.dat")
#
#f.seek(array_offsets[1])
