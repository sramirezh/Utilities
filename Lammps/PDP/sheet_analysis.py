#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This scripts reads the polymerDp spreadsheet
Args:
    Input filen name
Returns:


@author: simon
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
#import argparse
#
#
#
#parser = argparse.ArgumentParser(description='This script evaluates the average of a quantity')
#parser.add_argument('FileName', metavar='InputFile', type=str,
#                    help='Input filename')
#
#parser.add_argument('--min', help='Number of timesteps to be discarded', default=1000, type=int)
#
#
#args = parser.parse_args()
#min_limit=args.min
#InputFile=args.FileName


input_file="Results.dat"

# This is to remove all the headers and just start reading in the columns if there are
with open(input_file, 'r') as data_file:
    while(data_file.read(1)=='#'):
        last_pound_pos = data_file.tell()
        data_file.readline()
    data_file.seek(last_pound_pos+1) #The one is added to avoid reading a non-existing column in the names as there lines are "# "
    data=pd.read_csv(data_file,sep="\t",header=0)
    

#data=np.genfromtxt(input_file)
Names= list(data.columns.values)
size=len(Names)
data1=data.as_matrix()


datasol=pd.read_csv("Statistic_summary.dat",sep="\t",header=0)

"""
*******************************************************************************
CLASS DEFINITION
*******************************************************************************
"""
class LJInteraction(object):
    """
    Every pair of LJ interactions has its own
    """
    
    def __init__(self,lj_parameters):
        self.epsilon=float(lj_parameters[0])
        self.sigma=float(lj_parameters[1])
        self.forces=[]
        self.properties=[]
        self.propertie_names=[]
        
    def addforce(self,new_force):
        force=float(new_force.strip('\n/dDP'))/1000
        self.forces.append(force)
        
    def addproperties(self,properties):
        values=[]
        names=[]
        for element in properties:
            name,value=element.strip("\n").split("=")
            values.append(float(value))
            names.append(name)
        self.properties.append(values)
        self.propertie_names.append(names)
    
"""
*******************************************************************************
Functions
*******************************************************************************
"""

def buil_data():
    """
    Function to initialise all the elements of the class
    """

    interactions=[]
    with open("Statistic_summary.dat", 'r') as f:
      lines = f.readlines()
    
    i=0
    nproperties=11 #The number of properties per force in the input file (TO IMPROVE)
    inter=[]
    count=0
    
    while i<len(lines):
        if re.search("\AE_*",lines[i] ): #Finding the LJ parameters
            interactions.append(LJInteraction(re.findall(r"[-+]?\d*\.?\d+", lines[i])))
            force=[]
            print lines[i]
            i=i+1
            count+=1
        else: 
            if re.search("\AdDP*",lines[i] ):
                interactions[-1].addforce(lines[i])
                i+=1
                properties=[]
                for j in xrange(nproperties):
                    properties.append(lines[i])
                    i+=1
                interactions[-1].addproperties(properties)
                    
            i+=1
        
    
"""
*******************************************************************************
CLASS DEFINITION
*******************************************************************************
"""
        



#
#for line in file:
#    if re.search("\AE_*", line):
#        counter+=1
#        print line
#        interactions.append(LJInteraction(re.findall(r"[-+]?\d*\.?\d+", line)))
#        for line in file:
#            if re.search("dDP*", line):
#                print subline
        
        

        
        

    







#
#"""
################################################################################
#Starting the plot
################################################################################
#"""
#
#plt.close("all")
#
#
#ax=data.plot.scatter(x="Delta_Cs",y="Mobility/rg")
#ax.set_xlabel("Delta $C_s$ [$1/\sigma^3$]",fontsize=16)
#ax.set_ylabel("$\mu/R_g$",fontsize=16)
#plt.grid()
#plt.savefig("MobilityRg_Delta_Cs.pdf")
#
#
#
#ax=data.plot.scatter(x="Delta_Cs",y="Av_Mobility")
#ax.set_xlabel("Delta $C_s$ [$1/\sigma^3$]",fontsize=16)
#ax.set_ylabel("$\mu$",fontsize=16)
#plt.grid()
#plt.savefig("Mobility_Delta_Cs.pdf")






