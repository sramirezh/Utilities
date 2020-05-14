#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 09:58:32 2020
Script to analyse density distributions
@author: simon
"""

import sys
import os
from shlex import split
import matplotlib.pyplot as plt
import linecache
import re
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.General.Thermo_Analyser as ta
import Lammps.lammps_utilities as lu



class DensityDistribution(object):
    
    def __init__(self, filename, log_file):
        self.filename = filename
        self._get_data()
        self.sim = lu.Simulation(log_file)

    def _get_data(self):
        self.data_frame = cf.read_data_file(self.filename)
        


    


