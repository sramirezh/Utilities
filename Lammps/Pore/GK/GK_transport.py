#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:55:09 2019
This scripts computes the correlation between fluxes
TODO optimise and simplify
@author: sr802
"""


import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.Pore.qsub.simulation_results as sr
import glob

import flux_correlation as fc 


try:
    from uncertainties import unumpy,ufloat
except ImportError as err2:
    print(err2)


import matplotlib.pyplot as plt


cwd = os.getcwd() #current working directory

def run_correlation_analysis(folder,input_file,delta_t, save = "True"):
    """
    Very specific function
    
    Computes the correlations using the general class correlation
    Args:
        folder is where the input_file is located
        input_file is the name of the file containing the fluxes
        save True if you want to save the correlation instances as pkl
    """
    
    cwd = os.getcwd()
    os.chdir(folder)
    Data=cf.read_data_file(input_file)
    
    
    data1=Data.values
    times=(data1[:,0]-data1[0,0])*delta_t
    max_delta=int(len(times)*0.004) #Maximum delta of time to measure the correlation
    
    # Definitions of the fluxes
    total_flux=fc.flux(data1[:,[1,3,5]],times,"Q")
    solute_excess=fc.flux(data1[:,[2,4,6]],times,"J_s-c_s^BQ")
    
    
    
    c11=fc.correlation(total_flux,total_flux,max_delta)
    c11.evaluate()

    
    c12=fc.correlation(total_flux,solute_excess,max_delta)
    c12.evaluate()

    
    c21=fc.correlation(solute_excess,total_flux,max_delta)
    c21.evaluate()

    
    c22=fc.correlation(solute_excess,solute_excess,max_delta)
    c22.evaluate()
    
    
    if save=="True":
        c11.save('c11')
        c12.save('c12')
        c21.save('c21')
        c22.save('c22')
    
    
    os.chdir(cwd)
    
    return [c11,c12,c21,c22]


def load_correlations(folder):
    cwd = os.getcwd()
    os.chdir(folder)
    
    c11 = cf.load_instance("c11.pkl")
    c22 = cf.load_instance("c22.pkl")
    c12 = cf.load_instance("c12.pkl")
    c21 = cf.load_instance("c21.pkl")
    
    os.chdir(cwd)
    
    return [c11,c12,c21,c22]
    

# Input parameters

delta_t=0.005
root = "."
directory_pattern='[0-9]*'
input_file='fluxes.dat'



# =============================================================================
#  Checking the directories
# =============================================================================

prep_results = sr.check_n_analyse(root,directory_pattern)
prep_results.check_finished("fluxes.dat")
prep_results.check_stat("c11.pkl") #Assuming that if the first correlation is not ready, need to compute all

finished_directories=prep_results.dir_fin
unfinished_correlation=prep_results.dir_stat

print("\nFinished directories %s\n"%finished_directories)
print("\nUnfinished correlations in %s\n"%unfinished_correlation) 


# Array with the total correlations
c11_array=[]
c12_array=[]
c21_array=[]
c22_array=[]

# perform the analysis if there are no instances saved in the main folder
if len(glob.glob('c1*'))<1:
    for folder in finished_directories:
        #When Correlation still need to be run 
        #TODO change this to just gathering the results
        print("Now in folder %s\n"%folder)
        
        if folder in unfinished_correlation[0]:
            print("Running the correlation analysis in %s\n"%folder)
            c11,c12,c21,c22=run_correlation_analysis(folder,input_file, delta_t)
        
        else:
            print("Loading the results in %s\n"%folder)
            c11,c12,c21,c22=load_correlations(folder)
            
        c11_array.append(c11.cor[-1])
        c12_array.append(c12.cor[-1])
        c21_array.append(c21.cor[-1])
        c22_array.append(c22.cor[-1])
        
        #Extracting the times from one of the correlation
        times = c11.times
    
    
    
    c11_total=np.sum(np.array(c11_array),axis=0)/len(c11_array)
    c12_total=np.sum(np.array(c12_array),axis=0)/len(c12_array)
    c21_total=np.sum(np.array(c21_array),axis=0)/len(c21_array)
    c22_total=np.sum(np.array(c22_array),axis=0)/len(c22_array)
        
    
    
    
    #Creating the bundle instances 
    c11=fc.bundle_correlation(c11_total,times,"Q","Q")
    c12=fc.bundle_correlation(c12_total,times,"Q","J_s-c_s^BQ")
    c21=fc.bundle_correlation(c21_total,times,"J_s-c_s^BQ","Q")
    c22=fc.bundle_correlation(c22_total,times,"J_s-c_s^BQ","J_s-c_s^BQ")
    
    #TODO this could be done iterating over all the instances
    c11.save('c11')
    c12.save('c12')
    c21.save('c21')
    c22.save('c22')

else:
    print ("loading c11\n")
    c11 = cf.load_instance("c11.pkl")
    print ("loading c22\n")
    c22 = cf.load_instance("c22.pkl")
    print ("loading c12\n")
    c12 = cf.load_instance("c12.pkl")
    print ("loading c21\n")
    c21 = cf.load_instance("c21.pkl")
    


# =============================================================================
# Ploting all the correlations
# =============================================================================


plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax=plt.subplots()

fig,ax = c11.plot(fig,ax)
ax.set_xscale('log')
ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
plt.tight_layout()
plt.savefig("correlation11.pdf")


#For c12
fig,ax=plt.subplots()

fig,ax = c12.plot(fig,ax)
ax.set_xscale('log')
ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
plt.tight_layout()
plt.savefig("correlation12.pdf")



#For c21
fig,ax=plt.subplots()

fig,ax = c21.plot(fig,ax)
ax.set_xscale('log')
ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
plt.tight_layout()
plt.savefig("correlation21.pdf")



#For c22
fig,ax=plt.subplots()

fig,ax = c22.plot(fig,ax)
ax.set_xscale('log')
ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
plt.tight_layout()
plt.savefig("correlation22.pdf")







# =============================================================================
# Plot the cross coefficients
# =============================================================================

xmax = 10 # Integration limitc
fig,ax = plt.subplots()

fig,ax = c12.plot(fig,ax,ax_label=False)
ax.lines[-1].set_label(r'$\langle (%s(t))(%s(0)) \rangle$'%(c12.flux1_name,c12.flux2_name)) #modifying the label of the last created line
fig,ax = c21.plot(fig,ax,ax_label=False)
ax.lines[-1].set_label(r'$\langle (%s(t))(%s(0)) \rangle$'%(c21.flux1_name,c21.flux2_name)) #modifying the label of the last created line
ax.set_xscale('log')
ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
ax.axvline(x = xmax, c='black')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("crossed.pdf")




# =============================================================================
# Performing the integration
# =============================================================================

V = 20**3 # Simulation box volume

#Todo, this could be added to each integral

print("The c11 is %s\n" %c11.transport_coeff(1,V,0,xmax))
print("The c12 is %s\n" %c12.transport_coeff(1,V,0,xmax))
print("The c21 is %s\n" %c21.transport_coeff(1,V,0,xmax))
print("The c22 is %s\n" %c22.transport_coeff(1,V,0,xmax))

f=open("GK.out",'w')
f.write("The c11 is %s\n" %c11.transport_coeff(1,V,0,xmax))
f.write("The c12 is %s\n" %c12.transport_coeff(1,V,0,xmax))
f.write("The c21 is %s\n" %c21.transport_coeff(1,V,0,xmax))
f.write("The c22 is %s\n" %c22.transport_coeff(1,V,0,xmax))
f.close




