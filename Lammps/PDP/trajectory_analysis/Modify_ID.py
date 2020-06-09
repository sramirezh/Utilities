#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 16:56:21 2018


This scripts reads an atom file and replaces the position


@author: simon
"""
import re
import numpy as np 
from subprocess import Popen,PIPE
from shlex import split
import pandas as pd
#import argparse
import re
import linecache
import matplotlib.pyplot as plt


#parser = argparse.ArgumentParser(description='This script evaluates the trajectory file of a polymer')
#parser.add_argument('FileName', metavar='InputFile', type=str,help='Input filename')

#parser.add_argument('--min', help='Number of timesteps to be discarded', default=1000, type=int)


#args = parser.parse_args()
#InputFile=args.FileName


file_name="poly.atom"

#Reading the file
file=open(file_name,'r')
content=file.read()
file.close()

print() 
initial=1
final=int(linecache.getline(file_name, 4))
linecache.clearcache()


#Finding the first element of the polymer
#Finds a beginning of a line with 30 a space and 3, replaces with 30 1
initial_search=r'^'+re.escape(str(initial))+r'\s3'
initial_replace=re.escape(str(initial))+r' 1'
intermediate = re.sub(initial_search, initial_replace, content, flags = re.M) 

#Finding the last element of the polymer 
final_search=r'^'+re.escape(str(final))+r'\s3'
final_replace=re.escape(str(final))+r' 2'
output = re.sub(final_search, final_replace, intermediate, flags = re.M)

#Writing the output
file2=open("pmodif.atom",'w')
file2.write(output)
file2.close()
