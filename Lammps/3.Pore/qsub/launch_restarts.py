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
from simulation_utilities import simulation
import argparse

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script



# =============================================================================
# Main
# =============================================================================
def main(name,root,template,conf_folder,n_conf,epsilon,force,run):       
    #Getting the path to all the restart files
    files=glob.glob('%s/*'%conf_folder)
    
    home=root+'/'+name+'_%s'%force
    
    times=cf.extract_digits(files)
    times=[str(int(time[-1])) for time in times]
    #Takign the last N configurations
    
    
    conf_times=times[-n_conf:]
    files_analysis=cf.parameter_finder(files,conf_times)
    
    shutil.rmtree(home,ignore_errors=True)
    
    for i in files_analysis:
          
    # =============================================================================
    #     The extraction of the parameters for the simulation comes here
    # =============================================================================
    
        time=int(cf.extract_digits(files[i])[-1])
        name=str(time)
        restart=files[i]
        
    # =============================================================================
    #     Creating the simulation instance 
    # =============================================================================
        
        sim=simulation(home,template,name,restart)
        sim.create_folder()
        sim.create_qsub('short',1,16,24,'input.lmp')
    # =============================================================================
    #     #Mofications to the files here (THIS IS SPECIFIC)
    # =============================================================================
        
        file_name="input.lmp"
        file_path=sim.folder+'/'+file_name
        value_modify=sim.initial_conf.split('/')[-1]
        cf.modify_file(file_path,'read_restart','read_restart\t%s\n'%value_modify)
        
        value_modify=force
        cf.modify_file(file_path,'force','variable\tforce equal %s\n'%value_modify)
        
        
        file_name="in.interaction"
        file_path=sim.folder+'/'+file_name
        value_modify=epsilon
        cf.modify_file(file_path,'2 3','pair_coeff\t2 3 %s 1.0\n'%value_modify)
        
    # =============================================================================
    #     Running the simulation
    # =============================================================================
        if run==True:
            sim.run_simulation()
	        os.chdir(cwd)
            
if __name__ == "__main__":
    """
    THIS IS VERY SPECIFIC
    The arguments of this depend on the application
    """
    parser = argparse.ArgumentParser(description='Launch simulations from restart',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-name', metavar='name',help='Name of the folder to keep all the simulations',default='mu_force')
    parser.add_argument('-conf_folder', metavar='path_restart',help='Directory of the restart files',default='./Try/particle/')
    parser.add_argument('-template', metavar='path_template',help='Directory to take as template',default=cwd+'/Template')
    parser.add_argument('-root', metavar='root directory',help='Directory to create the folder for the simulations',default=cwd)
    parser.add_argument('-n_conf',metavar='n conf',help='number of configurations starting from the last',default=5,type=int)
    parser.add_argument('-epsilon',metavar='epsilon',help='monomer solute interaction',default=3.0,type=float)
    parser.add_argument('-force',metavar='force',help='Force on the solutes',default=0.01,type=float)
    parser.add_argument('-run',metavar='run',help='Define if run simulations or not. If not, just creates the folder structure',default=False,type=cf.str2bool)
    
    args = parser.parse_args()

    
    main(args.name,args.root,args.template,args.conf_folder,args.n_conf,args.epsilon,args.force,args.run)
