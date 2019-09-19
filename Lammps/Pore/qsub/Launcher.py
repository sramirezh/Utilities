#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 09:24:03 2019
"Creates replicas of simulations starting from configurations during the equilibration"
@author: sr802
"""

import glob
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import shutil
from .simulation_utilities import simulation

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script





    

#Main
        
#Getting the path to all the restart files
files=glob.glob('./Try/particle/*')
template=cwd+'/Template'
home=cwd+'/mu_force'

times=cf.extract_digits(files)
times=[str(int(time)) for time in times]
#Takign the last N configurations

n_conf=4

conf_times=times[-n_conf:]
files_analysis=cf.parameter_finder(files,conf_times)

shutil.rmtree(home,ignore_errors=True)

for i in files_analysis:
      
    #The extraction of the parameters for the simulation comes here
    
    path="%s/%s"%(home,times[i])
    time=int(cf.extract_digits(files[i])[-1])
    name=str(time)
    restart=files[i]
    
    #Creating the simulation instance 
    
    sim=simulation(home,template,name,restart)
    sim.create_folder()
    sim.create_qsub('test',1,1,1,'input.lmp',)
    #Mofications to the files here 
    
    file_name="input.lmp"
    file_path=sim.folder+'/'+file_name
    value_modify=sim.initial_conf.split('/')[-1]
    cf.modify_file(file_path,'read_restart','read_restart\t%s\n'%value_modify)
    
    
