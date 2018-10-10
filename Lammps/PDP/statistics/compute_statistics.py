#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 17:37:13 2018
s
@author: sr802
"""
import os
import glob
import sys
from joblib import Parallel, delayed
import multiprocessing
import time

Utilities_path=os.path.join(os.path.dirname(__file__), '../../../')
sys.path.append(Utilities_path) #This falls into Utilities path
import Lammps.core_functions as cf
import Others.Statistics.FastAverager as stat

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

def make_dp_poly():
    
    path_dp_poly=Utilities_path+"Lammps/PDP/f_programs/"
    out, error = cf.bash_command("""make -C %s -f Makefile""" %path_dp_poly)
    path_dp_poly+="dp_poly"
    return path_dp_poly

def gather_statistics(directories):
    """
    Gathers all the information from a previous statistics analyisis

    creates the file "Statistics_summary.dat"
    """
    
    os.chdir(cwd)
    f=open(cwd+"/Statistic_summary.dat",'w')
    
    for directory in directories:
        
        f.write( "############################################################################\n")
        f.write(directory+"\n")        
        
        os.chdir(directory)
        print os.getcwd()
        print directory
        forces=glob.glob("dDP*")
        forces.sort(key=lambda f: int(filter(str.isdigit, f)))
        for force in forces:
            print force
            f.write("\n"+force+"\n")
            os.chdir(force)
            
            #Results from dp_poly
            with open("average_info.dat",'r') as ave_info:
                f.writelines(ave_info.readlines())
                
            #Results from Fast averager    
            with open("statistics.dat",'r') as ave_info:
                f.writelines(ave_info.readlines()[1:])
            
            os.chdir("../")
    
        os.chdir(cwd)
            
    f.close()
    return
    

def run_analysis(interaction,force,path_dp_poly,s,d,n,dmin):
    
    initial_directory=os.getcwd()
    
    os.chdir("%s/%s"%(interaction,force))
    files=glob.glob('conf/*.gz')
    
    n=int(cf.extract_digits(files)[-1])

    #Results from dp_poly
    out,error=cf.bash_command("%s -s %s -d %s -n %s" %(path_dp_poly,s,d,n))        
    #Results from Fast averager    
    stat.fast_averager("vdata.dat",dmin)
    os.chdir(initial_directory)
    return

def compute_statistics(directories,s, d, n, dmin):
    
    """
    runs the statistics anaylisis with dp_poly and fast_averager
    
    Args:
    interactions: folders where this algoritm is going to run
    s: initial step
    d: interval
    n: final time step  [This is always overrided]
    dmin: samples to be discarded from vdata
    
    creates the file "Statistics_summary.dat"
    """    
    
    os.chdir(cwd)

    path_dp_poly=make_dp_poly()
    
    num_cores = multiprocessing.cpu_count()
    
    parameters=[]
    for directory in directories:
        os.chdir(directory)
        forces=glob.glob("dDP*")
        for force in forces:
            parameters.append([directory,force])
        os.chdir("../")
    
    t=time.time()
    
    #Parallel analysis
    Parallel(n_jobs=num_cores,verbose=10)(delayed(run_analysis)(param[0], param[1] ,path_dp_poly,s,d,n,dmin) for param in parameters)
    

    print "This is the time in paralell %f" %(time.time()-t)
    
    gather_statistics(directories)
    
    return
    