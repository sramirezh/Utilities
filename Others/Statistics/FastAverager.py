#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This scripts reads a Lammps output file and averages all the columns

Args:
    Input filen name
Returns:
    A message with the averages.
    A File called statistics.dat

@author: simon
"""


import pandas as pd
import numpy as np
import argparse
from scipy import stats
import os
import sys
from Functions import Autocorrelation 

sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf

"""
*******************************************************************************
Functions
*******************************************************************************
"""

def exclude_parameters(names, p_to_exclude):
    """
    return an array containig the indexes of the parameters to exclude
    Args:
        names is the list of strings containing the names of the parameters.
        p_to_exclude is a list of strings containing the parameters to exclude
    """
    i_exclude=list(map(lambda x: cf.parameter_finder(names,x),exclude))
    i_exclude=np.concatenate(np.array(i_exclude))
    
    return np.int(i_exclude)
    
    

"""
*******************************************************************************
Main
*******************************************************************************
"""



parser = argparse.ArgumentParser(description='This script evaluates the average of a quantity')
parser.add_argument('FileName', metavar='InputFile', type=str,
                    help='Input filename')

parser.add_argument('--min', help='Number of timesteps to be discarded', default=1000, type=int)


args = parser.parse_args()
min_limit=args.min
InputFile=args.FileName




# This is to remove all the headers and just start reading in the columns if there are
with open(InputFile, 'r') as data_file:
    while(data_file.read(1)=='#'):
        last_pound_pos = data_file.tell()
        data_file.readline()
    data_file.seek(last_pound_pos+1) #The one is added to avoid reading a non-existing column in the names as there lines are "# "
    data=pd.read_csv(data_file,sep=" ",header=0)
    

names= list(data.columns.values)
data1=data.as_matrix()[min_limit::]


#Excluding some data that does not need to be analysed
exclude=["time", "Chunk", "Coord1"]
i_delete=exclude_parameters(names, exclude)
data_to_analyse=np.delete(data1,i_delete,axis=1)
names_to_analyse=np.delete(names,i_delete)
size=len(names_to_analyse)
averages=np.average(data_to_analyse,axis=0)


#Autocorrelation analysis
n,m=np.shape(data_to_analyse)
Correlation,time=Autocorrelation(data_to_analyse, averages)
error=np.sqrt(Correlation[0,:]*2*time/(n+1))


print "The averages are:\n"
file=open("statistics.dat",'w')
for i in xrange(size):
    print "%s = %lf %lf"%(names_to_analyse[i],averages[i],error[i])
    file.write("%s = %lf %lf\n"%(names_to_analyse[i],averages[i],error[i]))
file.close()

print "\nCreated a file statistics.dat with the averages"
