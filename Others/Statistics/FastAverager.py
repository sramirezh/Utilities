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

Names= list(data.columns.values)
size=len(Names)
data1=data.as_matrix()
Averages=np.average(data1[min_limit::],axis=0)

print "The averages are:\n"
for i in xrange(size):
    print "%s = %lf"%(Names[i],Averages[i])

np.savetxt('statistics.dat',Averages)
print "\nCreated a file statistics.dat with the averages"
