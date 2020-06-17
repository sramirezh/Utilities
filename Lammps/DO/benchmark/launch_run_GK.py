#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Jun 17 09:24:03 2020
Creates replicas of simulations of GK based on a template
@author: sr802
"""

import glob
import sys
import os
from random import seed
from random import randint
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import shutil
from Lammps.simulation_utilities import simulation_launcher
import argparse

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

# Initiating the random number generator
seed(1)
# =============================================================================
# Main
# =============================================================================
def main(name, root, template, n_run, run):
    """
    Args:
        name:Name of the folder to keep all the simulations
        root:Directory to create the folder for the simulations
        template: Directory to take as template
        n_run: number of copies to run
        identifier: Numerical identifier
        
    """       
    home = root+'/'+name
    shutil.rmtree(home, ignore_errors = True)
    
    for i in range(n_run):
          
    # =============================================================================
    #     The extraction of the parameters for the simulation comes here
    # =============================================================================
    
        name= 'r%s'%i
        
    # =============================================================================
    #     Creating the simulation instance 
    # =============================================================================
        
        sim = simulation_launcher(home, template, name)
        sim.create_folder(keep_qsub = True)

    # =============================================================================
    #     #Mofications to the files here (THIS IS SPECIFIC)
    # =============================================================================
        
        file_name = "input.lmp"
        file_path = sim.folder+'/'+file_name
        value_modify = randint(0, 1000000)
        cf.modify_file(file_path, 'vel_seed', 'variable\tvel_seed equal %s \n'%value_modify)
        
        
        file_name = "run.qsub"
        file_path = sim.folder+'/'+file_name
        cf.modify_file(file_path, '#PBS -N', '#PBS -N r%s\n'%i)
        
    # =============================================================================
    #     Running the simulation
    # =============================================================================
        if run == True:
            sim.run_simulation()
        os.chdir(cwd)
            
if __name__ == "__main__":
    """
    THIS IS VERY SPECIFIC
    The arguments of this depend on the application
    """
    parser = argparse.ArgumentParser(description='Launch copy of template',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-name', metavar='name',help='Name of the folder to keep all the simulations',default='test')
    parser.add_argument('-template', metavar='path_template',help='Directory to take as template',default=cwd+'/Template')
    parser.add_argument('-root', metavar='root directory',help='Directory to create the folder for the simulations',default=cwd)
    parser.add_argument('-n_run',metavar='n conf',help='number of copies to run',default=5,type=int)
    parser.add_argument('-run',metavar='run',help='Define if run simulations or not. If not, just creates the folder structure',default=False,type=cf.str2bool)
    
    args = parser.parse_args()

    
    main(args.name, args.root, args.template, args.n_run, args.run)
