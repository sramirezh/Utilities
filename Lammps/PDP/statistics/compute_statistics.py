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
import Lammps.PDP.trajectory_analysis.poly_analysis as ta

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

def make_dp_poly():
    """
    Compiles dp_poly with fortran
    """
    
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
    
    
    directories=filter_directories(directories)
    
    f=open(cwd+"/Statistic_summary.dat",'w')
    
    for directory in directories:
        os.chdir(cwd)
        f.write( "############################################################################\n")
        f.write(directory+"\n")        
        
        os.chdir(directory)
        print(os.getcwd())
        print(directory)
        forces=glob.glob("dDP*")
        forces.sort(key=lambda f: int(list(filter(str.isdigit, f))))
        for force in forces:
            print(force)
            f.write("\n"+force+"\n")
            os.chdir(force)
                
            #Results from Fast averager    
            with open("statistics.dat",'r') as ave_info:
                f.writelines(ave_info.readlines()[1:])
            ave_info.close()
            
            #Results from Trajectory analysis
            
            with open("stat_strajectory.dat",'r') as ave_info:
                f.writelines(ave_info.readlines()[1:])
            ave_info.close()
            
            os.chdir("../")   
        os.chdir(cwd)
            
    f.close()
    return

def check_terminated_simulation(force):
    """
    First checks if the simulation started
    then checks if it finished or the last step
    Returns a counter that is 0 if the simulation chrashed before starting.
    """
    counter=1
    os.chdir("%s"%(force))
    if os.path.isfile("vdata.dat")==False:
        print("The simulation crashed before starting")
        counter=0
    else:
        tail=cf.bash_command("""tail -1 vdata.dat""")
        if "Total wall" in tail:
            print("This simulation terminated")
        else:
            last_step=cf.extract_digits(tail[0])[0]
            print("This simulation stoped at %s" %last_step)
    os.chdir("../")
    return counter
    
def filter_directories(directories):
    """
    checks if all the simulations inside all the directories finished, if not, deletes the directory from the analysis
    """
    os.chdir(cwd)
    
    for directory in directories:
        print(directory)
        os.chdir(directory)
        forces=glob.glob("dDP*")
        forces_finished=1
        for force in forces:
            print(force)
            forces_finished*=check_terminated_simulation(force)
        if forces_finished==0:
            directories.remove(directory)
        os.chdir("../")
        
    
    
    
    return directories
    

def run_analysis(interaction,force,dmin):
    
    initial_directory=os.getcwd()
    
    os.chdir("%s/%s"%(interaction,force))
    
    #Results from Fast averager   
    stat.fast_averager("vdata.dat",dmin, "statistics.dat")
    
    #Results from Poly Analysis
    
    cf.blockPrint()
    ta.poly_analysis("poly.atom",True)
    cf.enablePrint()
    
    stat.fast_averager("radius.dat",dmin, "stat_strajectory.dat" )
    
    #To delete all trajectory chunks
    #The bash command did not work 
    for f in glob.glob("*.cxyz"):
        os.remove(f)
    
    os.chdir(initial_directory)
    return



def compute_statistics(directories, dmin):
    
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

    directories=filter_directories(directories)
    
    
    parameters=[]
    for directory in directories:
        os.chdir(directory)
        forces=glob.glob("dDP*")
        for force in forces:
            parameters.append([directory,force])
        os.chdir("../")
    
    t=time.time()
    num_cores = multiprocessing.cpu_count()
    #Parallel analysis
    Parallel(n_jobs=num_cores,verbose=10)(delayed(run_analysis)(param[0], param[1],dmin) for param in parameters)
    

    print("This is the time in paralell %f" %(time.time()-t))
    
    gather_statistics(directories)
    
    return
    
