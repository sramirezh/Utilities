#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:09:58 2020
Script to compute the viscosity based on the stress autocorrelation functio
@author: simon
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.transformations as tr
import MDAnalysis.analysis.rdf as rdf
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
from scipy.spatial.distance import pdist,squareform
from tqdm import tqdm
from uncertainties import unumpy,ufloat
import Others.Statistics.FastAverager as stat
from scipy import optimize
import glob
from joblib import Parallel, delayed
import multiprocessing



data = cf.read_data_file("pressures.dat")