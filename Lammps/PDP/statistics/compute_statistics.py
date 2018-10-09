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

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

s=0#Initial step
d=10000 #Interval
n=100000#FinalTimestep
dmin=0 #Samples to be discarded from vdata.dat

os.chdir('/home/sr802/Delete/compute_statistics')
cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

#onlyfiles = [f for f in os.listdir(cwd) if os.isfile(os.join(cwd, f))]




numbers=[]
os.chdir(cwd)
directories=glob.glob('E_*')
f=open(cwd+"/Statistics_summary.dat",'w')
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
        out,error=cf.bash_command("../../programs/dp_poly -s %s -d %s -n %s" %(s,d,n))
        with open("average_info.dat",'r') as ave_info:
            f.writelines(ave_info.readlines())
        cf.bash_command("""python ~/Utilities/Others/Statistics/FastAverager.py vdata.dat --min %s""" %min)
        with open("statistics.dat",'r') as ave_info:
            f.writelines(ave_info.readlines()[1:])
        
        os.chdir("../")

    os.chdir(cwd)
        
        
        
        
        

    
    
    
    
            
    
    


f.close()
