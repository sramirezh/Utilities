#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 17:35:55 2019
Basic algorithm to get the results, this can be improved generalizing compute statistics and statistic_parameters
@author: sr802
"""

import os
import sys
import glob
Utilities_path=os.path.join(os.path.dirname(__file__), '../../../')
sys.path.append(Utilities_path) #This falls into Utilities path
import Lammps.core_functions as cf
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import Others.Statistics as stat
import copy

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
        print '\n %s' %directory
        finished=1
        finished*=check_terminated_simulation(directory)
        finished*=check_terminated_by_file(directory+'/'+key_file)
        if finished==0:
            directories.remove(directory)
        
    return directories

    

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
        print os.getcwd()
        print directory
        forces=glob.glob("dDP*")
        forces.sort(key=lambda f: int(filter(str.isdigit, f)))
        for force in forces:
            print force
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


fitfunc = lambda  x,*p: p[0] * x**2 + p[1]*x + p[2] #Fitting to parabola
# =============================================================================
# main 
# =============================================================================
    
directories=glob.glob('mu_variation/mu_*')


directories=filter_directories(directories)



pressure=[]
mu=[]
for d in directories:
    f=d+'/thermo.dat'
    pressure.append(extract_property_file('Press',f))
    mu.append(cf.extract_digits(d)[0])
    
pressure=np.array(pressure)
cf.set_plot_appearance()


#Fitting to a curve

popt, pcov=curve_fit(fitfunc, pressure[:,0], mu, sigma=pressure[:,1] ,p0=[0]*3)


mu_fit=np.polyval(popt,np.sort(pressure[:,0]))

p_target=1
mu_desired=np.polyval(popt,p_target)

plt.close('all')
fig,ax=plt.subplots()
ax.errorbar(mu,pressure[:,0],yerr=pressure[:,1],fmt='o')
ax.plot(mu_fit,np.sort(pressure[:,0]))
ax.axhline(y=1,c='black',ls=':')
ax.plot()
ax.set_xlabel(r'$\mu$')
ax.set_ylabel(r'$P$')

fig.tight_layout()

fig.savefig("mu_vs_p.pdf")

print ('\nThe desired chemical potential is %s' %mu_desired)

#Now need to fit and get the mu, print it in the plot








