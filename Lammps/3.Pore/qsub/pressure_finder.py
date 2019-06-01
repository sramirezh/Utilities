#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 09:24:03 2019
"Creates replicas of simulations starting at different chemical potentials to try to find the mu for the desired pressure
@author: sr802
"""

import glob
import sys
import os
import numpy as np
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
def main(name,root,template,n_simulations,limits):   
    print ("Remember that in the template the input file is named input.lmp")   

    
    home=root+'/'+name
    
    mu=np.linspace(limits[0],limits[1],20)
    
    
    shutil.rmtree(home,ignore_errors=True)
    
    for i in mu:
        
    # =============================================================================
    #     Creating the simulation instance 
    # =============================================================================
        folder_name='mu_%s'%i
        sim=simulation(home,template,folder_name)
        sim.create_folder()
        sim.create_qsub('test',1,1,1,'input.lmp')
        
    # =============================================================================
    #     #Mofications to the files here 
    # =============================================================================1
        
        file_name="input.lmp"
        file_path=sim.folder+'/'+file_name
        cf.modify_file(file_path,'mu1','variable\tmu1 equal %s\n'%i)
        cf.modify_file(file_path,'mu2','variable\tmu2 equal %s\n'%i)
        cf.modify_file(file_path,'Temp','variable\tTemp equal 1.0\n'%i)
        
        #running the statistics analyisis
        file_name='run.qsub'
        file_path=sim.folder+'/'+file_name
        lines_to_add='python ~/Utilities/Lammps/0.General/Log_Analysis/Thermo_Analyser.py log.lammps --min 0.3'
        cf.modify_file(file_path,'echo',lines_to_add,n_ocurrence=-1)
    
    # =============================================================================
    #     Running the simulation
    # =============================================================================
        sim.run_simulation()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Launch simulations from restart',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-name_folder', metavar='name',help='Name of the folder to keep all the simulations',default='mu_variation')
    parser.add_argument('-limits', metavar='limits',help='minimum and maximum mu to evaluate',default=[-3 ,3],nargs='+',type=float)
    parser.add_argument('-template', metavar='path_template',help='Directory to take as template',default=cwd+'/Template')
    parser.add_argument('-root', metavar='root directory',help='Directory to create the folder for the simulations',default=cwd)
    parser.add_argument('-n_simulations',metavar='n conf',help='Number of configurations starting from the last',default=5,type=int)
    parser.add_argument('-prefix',metavar='prefix ',help='Prefix for each folder',default='mu_',type=str)
    args = parser.parse_args()

    
    main(args.name_folder,args.root,args.template,args.n_simulations,args.limits)
