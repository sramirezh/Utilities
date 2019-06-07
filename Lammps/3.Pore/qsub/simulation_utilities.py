#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 17:48:19 2019

This is a modification of the class created by Raman Ganti
"""
from __future__ import division
from __future__ import print_function

import os
import glob
import shutil
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf


def sec_to_pbs_time(seconds, nodays=False):
    """
    clean solution from
    https://stackoverflow.com/questions/21323692/convert-seconds-to-weeks-days-hours-minutes-seconds-in-python
    "The idea behind this is that the number of seconds in your answer is the remainder after dividing 
    them in minutes; minutes are the remainder of dividing all minutes into hours etc... This version 
    is better because you can easily adjust it to months, years, etc..."
    """
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    if nodays:
        return "{0:02d}:{1:02d}:{2:02d}".format(int(hours), int(minutes), int(seconds))
    else:
        days, hours = divmod(hours, 24)
        return "{0:02d}:{1:02d}:{2:02d}:{3:02d}".format(int(days), int(hours), int(minutes), int(seconds))
    
    
    
class simulation(object):
    def __init__(self,home,template,name,initial_conf=None):
        """

        home is where the folder for all the simulation is going to be created
        template points to the folder where the template is hosted
        name ex. DP_0020, the name of the folder that contains everything
        initial_conf is the initial configuration file if it is required
        """
        self.initial_conf=initial_conf
        self.home=home
        self.template=template
        self.name=name

    def clean_template(self,keep_qsub=False):
        """
        Clean the template files except for in* or optionally the qsub file
        """
        useful_files=glob.glob(self.template+'/in*')
        if keep_qsub==True:
            useful_files.extend(glob.glob(self.template+'/*.qsub'))
        all_files=glob.glob(self.template+'/*')
        remove_files=[f for f in all_files if f not in useful_files]
        for fil in remove_files:
            os.remove(fil)
        
    
    def create_folder(self,keep_qsub=False):
        """
        Copies the template and adds the input configuration to the folder
        """
        
        #Deletes innecesary files and folders in template
        self.clean_template(keep_qsub)
        
        #Copy the template
        self.folder=self.home+'/%s'%self.name

        #Copy the template folder
        shutil.copytree(self.template,self.folder)
        
        #Copy the initial configuration
        if self.initial_conf!=None:
            shutil.copy(self.initial_conf,self.folder)

        
    def create_qsub(self,queue_type, nodes, cores, wall_time, lammps_script,lammps_version='/home/sr802/Programs/lammps-12Dec18/src/lmp_dexter'):
        
        self.qsub=BuildPBSScript(queue_type,nodes,cores,wall_time,"input.lmp")
        self.qsub.writePBSscript(self.folder+'/run.qsub',self.name)
        
        
    def run_simulation(self):
        os.chdir(self.folder)
        cf.bash_command("""qsub run.qsub""")
        return 0



    
class BuildPBSScript(object):
    """
    This is a completely general class, it does not and should not know anything about the naming conventions
    and structures of any particular library. If any such knowledge is necessary please make a derived class
    and overload the member functions 
    *queue_type [string]
    *nodes [int]
    *core [int] processors per node
    *walltime [hours]
    *command [string]: command line to execute e.g. python parallel_tempering.py args
    *outdir is the directory where to redirect the standard output
    """

    def __init__(self, queue_type, nodes, cores, wall_time, lammps_input, output_dir=None, no_days=False, spatial=False,
                 out_dir=None,lmp_version='/home/sr802/Programs/lammps-12Dec18/src/lmp_dexter'):
        self.qtype = queue_type
        self.nodes = nodes
        self.cores = cores
        self.s_wtime = wall_time * 60 * 60  # convert hours to seconds
        self.dhms_wtime = sec_to_pbs_time(self.s_wtime, nodays=no_days)  # DD:HH:MM:SS time
        self.lammps_input = lammps_input
        self.pbs_ready = False
        if out_dir and not os.path.isabs(out_dir):
            out_dir = os.path.abspath(out_dir)
        self.out_dir = out_dir
        self.current_dir = os.getcwd()
        self.output_dir = output_dir
        self.lammps_input=lammps_input
        self.lmp_version=lmp_version


    def writePBSscript(self, fname, job_name, spatial=False):
        """
        *fname [string]: name of the pbs bash script where to write
        *job_name [string]: name of the pbs job
        """
        print("writing PBS file")
#        if ".sh" not in fname:
#            fname += ".sh"
        f = open(fname, 'w')
        f.write('#PBS -N {0} \n'.format(job_name))
        f.write('#PBS -q {0} \n'.format(self.qtype))
        f.write('#PBS -l nodes={0}:ppn={1} \n'.format(self.nodes, self.cores))
        f.write('#PBS -l walltime={0} \n'.format(self.dhms_wtime))
        f.write('#PBS -o output.pbs \n')  # this directive places all the outputs in that file
        f.write('#PBS -e error.pbs \n')  # this directive places all the outputs in that file
        f.write('cd $PBS_O_WORKDIR')
        if self.out_dir!=None:
            f.write('#PBS -o {0} \n'.format(self.out_dir))
            f.write('\n')
            f.write('cd {0}\n'.format(self.output_dir))
        f.write('\n')
        f.write('echo Starting job $PBS_JOBID \n')
        f.write('echo\n')
        f.write('echo PBS assigned me this node: \n')
        f.write('cat $PBS_NODEFILE \n')
        f.write('echo \n')
        f.write('echo \"Running ${job_name}\" \n')
        f.write('echo \n')
        # f.write('mkdir {0}\n'.format(self.output_dir))
        # f.write('cp {0} {1}\n'.format(self.lammps_script, self.output_dir))
        # f.write('cd {0}\n'.format(self.output_dir))
        # f.write('cp {0} {1} \n'.format(self.lammps_script, self.output_dir))
        
        f.write('mpirun -np {0} {1} -in {2}\n'.format(self.nodes * self.cores,self.lmp_version,
                                                                                        self.lammps_input))

        f.write('echo \n')
        f.write('echo \"Job finished. PBS details are:\" \n')
        f.write('echo \n')
        f.write('qstat -f ${PBS_JOBID} \n')
        f.write('echo \n')
        f.write('echo Finished at \`date\` \n')
        f.close()
        self.pbs_ready = True


