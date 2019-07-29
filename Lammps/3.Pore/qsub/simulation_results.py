#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 19:08:44 2019
Here I will include all the functions an classes necessary to gather and analyse simulations
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
import copy
import re
from scipy import optimize

try:
    from uncertainties import ufloat
except ImportError as err2:
    print err2

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

    for folder in directories:
        os.chdir(folder)
        thermo.thermo_analyser(log_name,dmin)
        os.chdir(cwd)
        
        

def initialise_sim_bundles(root_patter,directory_pattern,dictionary):
    
    #Needed parameters
    roots=glob.glob(root_pattern)
    mu=cf.extract_digits(roots, sort=False)
    
    
    bundles=[]

    for i,root in enumerate(roots):
        directories=glob.glob('%s/%s'%(root,directory_pattern))
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
        
        
        
        # Creating the statistic summary (MAYBE Get rid of this)
        array=gather_statistics(dir_fin,'Time',root)
    
    
        #Building the simulations
        times=construct_simulations(directories)
    
        #Creating the directory for the plots
        
        directory="%s/plots/"%root
        if not os.path.exists(directory):
            os.mkdir(directory)
    
        
        #Creating the bundle
        bundles.append(simulation_bundle(times,"mu",mu[i],root,dictionary=dictionary))
    
        #Plot for all the properties
    #    for prop in bundles[-1].simulations[-1].property_names:
    #        print "\ncreating the plot of %s"%prop
    #        if prop!="time":
    #            bundles[-1].plot_property(prop)
        bundles[-1].plot_all_properties()
        
        return bundles
        
def gather_statistics(directories,folder_name,root,files=["statistics.dat","thermo.dat"]):
    """
    Creates statistics summary
    
    Gathers all the information from a previous statistics analyisis

    creates the file "Statistics_summary.dat"
    
    Args:
        directories list of directories to look for properties
        folder_name name of the parameter that defines the subdirectories, ex time, temperature.
        root is where the output is going to be saved
    """
    os.chdir(cwd)
    print ("\nCreating the Statistic_summary.dat\n")
    f=open(cwd+"/%s/Statistic_summary%s.dat"%(root,root),'w')
    array=[]
    for directory in directories:
        
        dir_array=[]
        dir_array.append([folder_name,os.path.split(directory)[-1]])
        
        print os.path.split(directory)[-1]
        os.chdir(cwd)
        f.write( "############################################################################\n")
        f.write(directory.split('/')[-1]+"\n")        
        
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
#                    f.write("%s\t"%p_name+'\t'.join(map(str, p_value))+'\n')
                    f.write(line+'\n')
            ave_info.close()
        dir_array.append(p_array)
        array.append(dir_array)
        os.chdir(cwd)
    f.close()
    return array



"""
*******************************************************************************
CLASS DEFINITION
*******************************************************************************
"""
class simulation(object):
    """
    Defines a generic class of simulation that has a characteristic parameter ex= time, Temperature, etc and then fills the instance with properties
    """
    total = 0 #This is a class atribute

    @staticmethod  #Does not have self as is intended to be used by the class and not an object.
    def number_simulations():
        "Return the number of instances created in the class"
        return simulation.total


    def __init__(self,parameter_id,parameter_value):
        """
        Args:
            parameter_id= string with the name of the characteristic parameter of the simulation.
            parameter_value= value of the characteristic parameter
        """
        self.param_value=float(parameter_value)
        self.param_id=parameter_id
        simulation.total += 1
        self.properties=[]
        self.property_names=[]

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

    def add_properties(self,properties):
        """
        Args:properties is an array with elements [property_name, property_value], where property value could be an array
        """
        for element in properties:
            name,value=element
            if 'NaN' in value:
                value=float('nan')
            self.add_property(name,value)
        self.property_names.insert(0,self.param_id)
        self.properties.insert(0,self.param_value)
        
    def add_property(self,prop_name,prop_value):
        """
        adding just one property this could be used to add number of particles, etc
        """
        self.properties.append(prop_value)
        self.property_names.append(prop_name)

    def get_property(self,name,exact=False):
        """
        function to get the specific property
        Args:
            name: a string containing the name that wants to be found
            exact: True if the exact name is required
        """
        index=cf.parameter_finder(self.property_names,name,exact=exact)
        prop_names=[]
        prop_values=[]
        for i in index:
            prop_names.append(self.property_names[i])
            prop_values.append(self.properties[i])
        return prop_names,prop_values


fitfunc1 = lambda p, x: p * x  #Fitting to a line that goes through the origin
errfunc1 = lambda p, x, y, err: (y - fitfunc1(p, x)) / err #To include the error in the least squares

"""
*******************************************************************************
Class Inheritage creating the superclass
*******************************************************************************
"""

class simulation_bundle(simulation):
    def __init__(self,simulations,parameter_id,parameter_value,root,dictionary=None):
        """
        
        root is the directory where the plot folder and the statistic summary is going to be created
        """
        self.simulations=simulations
        self.param_value=float(parameter_value)
        self.param_id=parameter_id
        self.properties=[]
        self.property_names=[]
        self.root=root
        #Ask Shaltiel if this is oK?
        self.add_properties()
        self.dictionary=dictionary
        
    def add_properties(self):
        
        self.property_names.extend(self.simulations[-1].property_names)  #is this very Ugly?????
        
        #Now I will add the average with error
        
        for prop in self.property_names:
            array=[]
            for sim in self.simulations:
                value=sim.get_property(prop,exact=True)[1][0]
                if np.size(value)>1:  #Has error
                    array.append(ufloat(value[0],value[1])) 
                else:
                    array.append(ufloat(value,0))
            ave=sum(array)/len(array)
            self.properties.append([ave.n,ave.s])
        
        #Adding the property
        self.property_names.insert(0,self.param_id)
        self.properties.insert(0,self.param_value)
            
            
    def plot_property(self,p_name,plot_name=None,x_name=None,y_name=None,fit=False):
        """
        Could be extended to plot agains any parameter
        """
        
        directory="%s/plots_%s/"%(self.root,self.param_id)
        if not os.path.exists(directory):
            os.mkdir(directory)
    
        cf.set_plot_appearance()
        #Create array time vs vx_Sol
        x=[]
        y=[]
        y_error=[]
        
        for sim in self.simulations:
            x.append(sim.param_value)
            values=sim.get_property(p_name,exact=True)[1][0]
            if np.size(values)>1: #Has error
                y.append(values[0])
                y_error.append(values[1])
            else:
                y.append(values)
                
        
        fig,ax=plt.subplots()
        ax.errorbar(x,y,yerr=y_error,xerr=None,fmt='o')
        
        #Checking if there are axis names
        if x_name==None:
            ax.set_xlabel(self.property_names[1]) #The zeroth-property is the param_id
        else:
            ax.set_xlabel(x_name)
        if y_name==None:

            if p_name in self.dictionary.keys():
                y_name=self.dictionary[p_name]
            else:
                y_name=re.sub('_',' ',p_name)
            ax.set_ylabel(y_name)
        else:
            ax.set_ylabel(y_name)
        
        # Adding a fitting line
        if fit==True and np.size(values)>1: #Has error
            #Same procedure as in statistic_parameters.py
            pinit=[1.0]
            print x,y,y_error
            out = optimize.leastsq(errfunc1, pinit, args=(x, y, y_error), full_output=1)
            pfinal = out[0] #fitting coefficients
            x=np.insert(x,0,0)
            ax.plot(np.unique(x),fitfunc1(pfinal,np.unique(x)),linestyle='--')
        else:
            #Just plot the average
            ax.axhline(y=self.get_property(p_name,exact=True)[1][0][0],c='black',ls=':')
            
        
        plt.tight_layout()
        fig.savefig("%s/plots_%s/%s.pdf"%(self.root,self.param_id,p_name), transparent=True)
        
    
    def plot_all_properties(self):
        for i,prop in enumerate(self.simulations[-1].property_names):
            
            if prop!="time" and i>0:
                print "\ncreating the plot of %s"%prop
                self.plot_property(prop)

    

def read_properties(directory,files):
    """
    """
    os.chdir(directory)
    global p_array
    #Reading the properties
    p_array=[]
    for fil in files:
        with open(fil,'r') as ave_info:
            lines=ave_info.readlines()[1:]
            for line in lines:
                line=line.strip()
                p_name=line.split('=')[0][2::] #Removing the v_ or c_ from lammps
                p_name=p_name.strip()
                p_value=np.array(line.split('=')[1].split(),dtype=float) #Value generated by fast averager
                p_array.append([p_name,p_value])
        ave_info.close()
    
    return p_array
    


def construct_simulations(directories,files=["statistics.dat","thermo.dat"]):
    """

    """
    global properties
    times=[]
    os.chdir(cwd)
    
    print ("\nConstructing the simulation instances\n")
    
    
    for directory in directories:
        
        t=float(os.path.split(directory)[-1])
        times.append(simulation("time",t))
        print os.path.split(directory)[-1]
        
        properties=read_properties(directory,files)
        times[-1].add_properties(properties)
        os.chdir(cwd)
    return times
