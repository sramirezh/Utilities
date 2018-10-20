#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 12:56:45 2018
This files creates the table for the frenkel potential to be used in Lammps
@author: simon
"""

from __future__ import division
import argparse
import pandas as pd
import numpy as np
import warnings
import sys 
import os
import glob
import bisect 

warnings.filterwarnings("ignore")


sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path


