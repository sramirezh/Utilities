#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 17:37:13 2018
s
@author: sr802
"""



import os
import glob
import numpy as np
import sys
import time

Utilities_path=os.path.join(os.path.dirname(__file__), '../../../')
sys.path.append(Utilities_path) #This falls into Utilities path
import Lammps.core_functions as cf

s=0#Initial step
d=10000 #Interval
n=100000#FinalTimestep
dmin=0 #Samples to be discarded from vdata.dat

os.chdir('/home/sr802/Delete/compute_statistics')
cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

#onlyfiles = [f for f in os.listdir(cwd) if os.isfile(os.join(cwd, f))]


def make_dp_poly():
    
    path_dp_poly=Utilities_path+"Lammps/PDP/f_programs/"
    out, error = cf.bash_command("""make -C %s -f Makefile""" %path_dp_poly)
    print out, error
    path_dp_poly+="dp_poly"
    return path_dp_poly



numbers=[]
os.chdir(cwd)
directories=glob.glob('E_*')
f=open(cwd+"/Statistics_summary.dat",'w')

path_dp_poly=make_dp_poly()
path_to_averager=Utilities_path+"Others/Statistics/FastAverager.py"

for directory in directories:
    
    f.write( "############################################################################\n")
    f.write(directory+"\n")        
    
    os.chdir(directory)
    print os.getcwd()
    print directory
    forces=glob.glob("dDP*")
    for force in forces:
        print force
        f.write("\n"+force+"\n")
        os.chdir(force)
        files=glob.glob('conf/*.gz')
        #n=int(cf.extract_digits(files)[-1])
        
        #Results from dp_poly
        out,error=cf.bash_command("%s -s %s -d %s -n %s" %(path_dp_poly,s,d,n))
        with open("average_info.dat",'r') as ave_info:
            f.writelines(ave_info.readlines())
            
        #Results from Fast averager    
        out, error = cf.bash_command("""python %s vdata.dat --min %s""" %(path_to_averager,dmin))
        with open("statistics.dat",'r') as ave_info:
            f.writelines(ave_info.readlines()[1:])
        
        os.chdir("../")

    os.chdir(cwd)
        
        
        
        
        

    
    
    
    
            
    
    


f.close()
