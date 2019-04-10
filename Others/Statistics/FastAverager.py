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
    i_exclude=list(map(lambda x: cf.parameter_finder(names,x),p_to_exclude))
    i_exclude=np.concatenate(np.array(i_exclude))
    i_exclude=np.unique(np.array(i_exclude,dtype=int))

    return i_exclude



"""
*******************************************************************************
Main
*******************************************************************************
"""
def read_from_file(input_file):
    """
    Reads the input file and returns the data to be evaluates
    Args:
        input_file that contains the data
        returns the data
    """
    data=cf.read_data_file(input_file)
    names= list(data.columns.values)
    data1=data.values

    return data1,names

def fast_averager(input,min_limit=0,output_file="statistics.dat"):
    """
    Function to call from another python script
    Args:
        input: could be a file or data(np.array)
        min_limit Number of samples to be discarded (default 0)
        output file name of the file with the analysis output (default statistics.dat)
    Returns:
        output_array containing the [Name_of variable, average, Error_autocorrelation, Error_blocking, Error_simple,variance simple]
    """
    try:
        if os.path.exists(input):
            data,names=read_from_file(input)
            calculations(data,min_limit,output_file,names)
    except:
        data=input

    output_array=calculations(data,min_limit,output_file)

    return output_array




def calculations(data, min_limit,output_file, names=None):
    """
    This script evaluates the average of a quantity

    Args:
        names list containing the column names
        data is a numpy array containing the data
        min_limit Number of samples to be discarded

    It creates a file statistics.dat with the averages, containing the valiable name,
    the average, the Error_autocorrelation,  the  Error_blocking  and the  Error_simple
    """

    #If the data does not have header
    if names is None:

        # If it is a 1-D array it will need to be reshaped
        if len(np.shape(data))==1:
            data=np.reshape(data,(len(data),1))

        data1=data[min_limit::]
        data_to_analyse=data1
        n,size=np.shape(data_to_analyse)
        names_to_analyse=np.arange(0,size)

    else:
        #Excluding some data that does not need to be analysed
        exclude=["time", "Chunk", "Coord1","step"]
        if isinstance(names[0],basestring)==True:
            try:
                i_delete=exclude_parameters(names, exclude)
                data_to_analyse=np.delete(data1,i_delete,axis=1)
                print "skipped the next parameters from the analysis:"
                for j,i in enumerate(i_delete):
                    print "%s. %s"%(j,names[i])
                names_to_analyse=np.delete(names,i_delete)
            except:
                print "No columns to skip"
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
    variance_s=stats.variance(data_to_analyse)

    #Blocking analysis
    error_b=blocking_error(data_to_analyse)


    print "The Results are:\n"
    print "Property    Average    Error_autocorrelation    Error_blocking    Error_simple"
    file=open(output_file,'w')
    file.write("#print Property    Average    Error_autocorrelation    Error_blocking    Error_simple\n")
    for i in xrange(size):
        print "%s = %lf %lf %lf %lf"%(names_to_analyse[i],averages[i],error_c[i], error_b[i], error_s[i], variance_s[i])
        file.write("%s = %lf %lf %lf %lf\n"%(names_to_analyse[i],averages[i],error_c[i], error_b[i], error_s[i],variance_s[i] ))
    file.close()

    output_array=[]
    for i in xrange(size):
        output_array.append([names_to_analyse[i],averages[i],error_c[i], error_b[i], error_s[i],variance_s[i]])

    print "\nCreated a file %s with the averages"%output_file

    return output_array


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='This script evaluates the average of a quantity',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('filename', metavar='InputFile', type=str,
                        help='Input filename')

    parser.add_argument('--min', help='Number of samples to be discarded', default=0, type=int)
    parser.add_argument('--output',help='Name of the output file',default="statistics.dat",type=str)
    args = parser.parse_args()
    min_limit=args.min
    input_file=args.filename
    output_file=args.output

    fast_averager(input_file,min_limit,output_file)
