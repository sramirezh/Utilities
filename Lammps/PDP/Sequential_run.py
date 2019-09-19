#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
File created to run all the sequential equilibrations
@author: simon
"""


import os
import sys
from shutil import copyfile

sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf

cwd=os.getcwd()
path=os.path.dirname(__file__)

#evaluate stopping condition
epsilon_stop=12
epsilon_current=cf.extract_digits(cwd.split("/")[-2])[0]
if epsilon_current>epsilon_stop:
    sys.exit("Reached the maximum epsilon")
else:
    os.chdir("../../")
    out,err=cf.bash_command("""bash %s/bash_scripts/new_interaction.sh %s 1.0"""%(path,epsilon_current+0.5))
    os.chdir("E_%s_S_1.0/Equilibration"%(epsilon_current+0.5))
    copyfile(cwd+"/final_conf.dat", os.getcwd()+"/initial_conf.dat")
    out,err=cf.bash_command("""qsub run.qsub""")
    os.chdir(cwd)
