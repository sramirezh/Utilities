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
import numpy as np
import argparse
import os
import sys
from Functions import autocorrelation_error, blocking_error
from scipy import stats

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

parser.add_argument('--min', help='Number of samples to be discarded', default=1000, type=int)


args = parser.parse_args()
min_limit=args.min
input_file=args.FileName


data=cf.read_data_file(input_file)


names= list(data.columns.values)
data1=data.values[min_limit::]


#Excluding some data that does not need to be analysed
exclude=["time", "Chunk", "Coord1"]
if isinstance(names[0],basestring)==True:
    i_delete=exclude_parameters(names, exclude)
    data_to_analyse=np.delete(data1,i_delete,axis=1)
    names_to_analyse=np.delete(names,i_delete)
else: #If the data does not have header
    data_to_analyse=data1
    names_to_analyse=names


size=len(names_to_analyse)
averages=np.average(data_to_analyse,axis=0)


"""
*******************************************************************************
Error Analysis
*******************************************************************************
"""

#Autocorrelation Analysis

n,m=np.shape(data_to_analyse)
Correlation,time=autocorrelation_error(data_to_analyse, averages)
error_c=np.sqrt(Correlation[0,:]*2*time/(n+1))


#Simple error
error_s=stats.sem(data_to_analyse)

#Blocking analysis
error_b=blocking_error(data_to_analyse)

print "The Results are:\n"
print "Property    Average    Error_autocorrelation    Error_blocking    Error_simple"
file=open("statistics.dat",'w')
file.write("#print Property    Average    Error_autocorrelation    Error_blocking    Error_simple\n")
for i in xrange(size):
    print "%s = %lf %lf %lf %lf"%(names_to_analyse[i],averages[i],error_c[i], error_b[i], error_s[i])
    file.write("%s = %lf %lf %lf %lf\n"%(names_to_analyse[i],averages[i],error_c[i], error_b[i], error_s[i]))
file.close()

print "\nCreated a file statistics.dat with the averages"
