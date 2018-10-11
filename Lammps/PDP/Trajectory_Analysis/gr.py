#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 15:37:04 2018
This scripts analyses the pair correlation function from the *.gz files
@author: simon
"""


from __future__ import division
import numpy as np
import pandas as pd
import argparse
import linecache
import os
import sys
import poly_analysis as pa
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf


cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script



