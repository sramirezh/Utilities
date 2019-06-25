#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 17:35:55 2019

Gathers the flow from diffusio-osmotic simulations

Missing:
    Create class for the force...

@author: sr802
"""

import os
import sys
import glob
Utilities_path=os.path.join(os.path.dirname(__file__), '../../../')
sys.path.append(Utilities_path) #This falls into Utilities path
import Lammps.core_functions as cf
import Others.Statistics.FastAverager as stat
import Lammps.General.Log_Analysis.Thermo_Analyser as thermo
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import copy
import tqdm

cwd = os.getcwd() #current working directory



def check_terminated_simulation(folder_name):
    """
    GENERIC FUNCTION
    First checks if the simulation started
    then checks if it finished or the last step
    Returns a counter that is 0 if the simulation chrashed before starting.
    """
    counter=1
    cwd=os.getcwd()
    os.chdir(folder_name)
    if os.path.isfile("log.lammps")==False:
        print "The simulation crashed before starting"
        counter=0
    else:
        tail,error=cf.bash_command("""tail -2 log.lammps""")
        if "Total wall" in tail:
            print"This simulation terminated"
        else:
            last_step=cf.extract_digits(tail)
            print "This simulation stoped at %s" %last_step[0]
    os.chdir(cwd)
    return counter


def check_terminated_by_file(file_name):
    """
    GENERIC FUNCTION
    Checks for the existance of the file
    """
    counter=1
    if not os.path.exists(file_name):
        print("The file %s does not exist!" % file_name)
        counter=0
    
    return counter
    

def filter_directories(directories,key_file):
    """
    THIS IS A VERY SPECIFIC FUNCTION
    checks if all the simulations inside all the directories finished, if not, deletes the directory from the analysis
    Args:
        key_file: is the file that is checked inside the directories, if it is not there the directory is delated.
    """
    os.chdir(cwd)
    directories=copy.copy(directories)
    for directory in directories:
        print '\n%s' %directory
        finished=1
        finished*=check_terminated_simulation(directory)
        finished*=check_terminated_by_file(directory+'/'+key_file)
        if finished==0:
            directories.remove(directory)
        
    return directories


def extract_property_file(property_name,file_name):
    """
    Return the property value from a file
    """
    f=open(file_name,'r')
    lines=f.readlines()
    line_number=line_number=cf.parameter_finder(lines,property_name)[0]
    value=cf.extract_digits(lines[line_number])
    f.close()
    return value


def run_stat_file(directories,file_analysis,dmin,output_name):
    """
    THIS IS A GENERAL FUNCTION (Kind of replaces compute statistics from PDP)
    runs the statistical analysis on a given file for a list of folders
    Args:
        directories: list of directories to run the analysis
        file_analysis: The file to be analysed with the fast_averager
        dmin: the minimum of steps to be discarded as defined in fast_averager
        output_name for each file analysed
    """
    cwd=os.getcwd()
    for folder in directories:
        os.chdir(folder)
        stat.fast_averager(file_analysis,dmin,output_name)
        os.chdir(cwd)
        
def run_thermo_directories(directories,log_name,dmin):
    """ 
    THIS IS A VERY GENERAL FUNCTION (Kind of replaces compute statistics from PDP)
    Runs the thermo analysis inside the specific directories
    
    Args:
        directories: list of directories to run the analysis
        dmin: the minimum of steps to be discarded as defined in fast_averager
    """
    cwd=os.getcwd()
    for folder in directories:
        os.chdir(folder)
        thermo.thermo_analyser(log_name,dmin)
        os.chdir(cwd)
        
        
def gather_statistics(directories,folder_name,files=["statistics.dat","thermo.dat"]):
    """
    Creates statistics summary
    
    Gathers all the information from a previous statistics analyisis

    creates the file "Statistics_summary.dat"
    """
    os.chdir(cwd)
    print ("\nCreating the Statistic_summary.dat\n")
    f=open(cwd+"/Statistic_summary.dat",'w')
    array=[]
    for directory in directories:
        
        dir_array=[]
        dir_array.append([folder_name,os.path.split(directory)[-1]])
        
        print os.path.split(directory)[-1]
        os.chdir(cwd)
        f.write( "############################################################################\n")
        f.write(directory+"\n")        
        
        os.chdir(directory)
        
        p_array=[]

        for fil in files:
            with open(fil,'r') as ave_info:
                lines=ave_info.readlines()[1:]
                for line in lines:
                    line=line.strip()
                    p_name=line.split('=')[0][2::] #Removing the v_ or c_ from lammps
                    p_value=np.array(line.split('=')[1].split(),dtype=float) #Value generated by fast averager
                    p_array.append([p_name,p_value])
            
                f.writelines(ave_info.readlines()[1:])
            ave_info.close()
        dir_array.append(p_array)
        array.append(dir_array)
        os.chdir(cwd)
    f.close()
    return array


def parameter_finder(List, String):
    """
    Finds a string on a List and returns the position on the list
    """
    cont=0
    indexes=[]
    for s in List:
        if String in s:
            indexes.append(cont)
        cont+=1
    length=len(indexes)
    if length>1: print "There were several ocurrences"
    if length==0: print "No ocurrences found"
    return indexes


"""
*******************************************************************************
CLASS DEFINITION
*******************************************************************************
"""
class time(object):
    """
    Every pair of LJ interactions has its own
    """
    total = 0 #This is a class atribute

    @staticmethod  #Does not have self as is intended to be used by the class and not an object.
    def number_interactions():
        "Return the number of instances created in the class"
        return time.total


    def __init__(self,time):
        self.time=float(time)
        self.properties=[]
        self.property_names=[]
        time.total += 1

    def __str__(self):
        self.name=("time=%s" %(self.time))
        return self.name

    def __lt__(self,other):
        """
        In order to sort the times by their value
        """
        return self.time < other.time


#    def addforce(self,new_force):
#        force=float(new_force.strip('\n/dDP'))/1000
#        self.forces.append(force)

    def addproperties(self,properties):
        values=[]
        names=[]
        for element in properties:
            name,value=element.strip("\n").split("=")
            name=name.replace(" ","_")
            if 'NaN' in value:
                value=float('nan')
            elif len(value)>1:
                value=value.split()

            values.append(np.double(value))
            names.append(name)
        self.properties.append(values)
        self.property_names.append(names)


#    def compute_mobility(self):
#        self.mobility=[]
#        self.mob_rg=[] #Mobility over Rg
#        count=0
#        index_vx=parameter_finder(self.property_names[count],"vx_relative")[0]
#        index_rg=parameter_finder(self.property_names[count],"r_gyration")[0]
#        for force in self.forces:
#            if force!=0:
#                velocity=ufloat(self.properties[count][index_vx][0],self.properties[count][index_vx][1])
#                rg=ufloat(self.properties[count][index_rg][0],self.properties[count][index_rg][1])
#                mobility=-velocity/force
#                self.mobility.append(mobility)
#                if rg.n==0:
#                   self.mob_rg.append(10**8) #To avoid division by 0
#                else:
#                    self.mob_rg.append(mobility/rg)
#            count+=1

    def get_property(self,name):
        """
        function to get the specific property
        """
        index=parameter_finder(self.property_names[0],name)[0]
        count=0
        prop=[]
        for force in self.forces:
            prop_f=self.properties[count][index]
            prop.append(prop_f)
            count+=1

        return prop



    



        
# =============================================================================
# main 
# =============================================================================
    
directories=glob.glob('mu_force_0.1/*')



# =============================================================================
# Checking if the simulation and the required files are on each folder
# =============================================================================

print "\nChecking if the simulations finished with vdata\n"
dir_fin=filter_directories(directories,"vdata.dat")

print "\nChecking if the simulations finished with statistics\n"
dir_stat=filter_directories(dir_fin,"statistics.dat")

print "\nChecking if the simulations finished with thermo\n"
dir_thermo=filter_directories(dir_fin,"thermo.dat")



# =============================================================================
# Running the necessary analysis
# =============================================================================

#directories to run statistics
dir_run_vdata=[x for x in dir_fin if x  not in dir_stat]
run_stat_file(dir_run_vdata,"vdata.dat",0.3,"statistics.dat")



#Directories to run thermo analysis
dir_run_thermo=[x for x in dir_fin if x  not in dir_thermo]
run_thermo_directories(dir_run_thermo,"log.lammps",0.3)


array=gather_statistics(dir_fin,'Time')

#cf.set_plot_apparence()






#class force(object):
#    
#    def 


#pressure=[]
#mu=[]
#for d in directories:
#    f=d+'/thermo.dat'
#    pressure.append(extract_property_file('Press',f))
#    mu.append(cf.extract_digits(d)[0])
#    
#pressure=np.array(pressure)
#cf.set_plot_appearance()
#
#
##Fitting to a curve
#
#popt, pcov=curve_fit(fitfunc, pressure[:,0], mu, sigma=pressure[:,1] ,p0=[0]*3)
#
#
#mu_fit=np.polyval(popt,np.sort(pressure[:,0]))
#
#p_target=1
#mu_desired=np.polyval(popt,p_target)
#
#plt.close('all')
#fig,ax=plt.subplots()
#ax.errorbar(mu,pressure[:,0],yerr=pressure[:,1],fmt='o')
#ax.plot(mu_fit,np.sort(pressure[:,0]))
#ax.axhline(y=1,c='black',ls=':')
#ax.plot()
#ax.set_xlabel(r'$\mu$')
#ax.set_ylabel(r'$P$')
#
#fig.tight_layout()
#
#fig.savefig("mu_vs_p.pdf")
#
#print ('\nThe desired chemical potential is %s' %mu_desired)
#
##Now need to fit and get the mu, print it in the plot








