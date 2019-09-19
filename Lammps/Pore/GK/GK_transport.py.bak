#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:55:09 2019
This scripts computes the correlation between fluxes
TODO optimise and simplify
@author: sr802
"""

from __future__ import division
import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.Pore.qsub.simulation_results as sr
import glob

from flux_correlation import flux, correlation


try:
    from uncertainties import unumpy,ufloat
except ImportError as err2:
    print err2


# Put this inside the argparser
try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
except ImportError as err:
    print err

cwd = os.getcwd() #current working directory


# =============================================================================
# Class definition
# =============================================================================
# A child class that i need just to plot the correlations
class bundle_correlation(correlation):
    def __init__(self,cor,times,flux1_name,flux2_name):
        self.dimension=1
        self.cor = [cor]
        self.times = times
        self.norm = cor[0].nominal_value
        self.cor_norm = [cor/self.norm]
        self.flux1_name = flux1_name
        self.flux2_name = flux2_name
        
        
    def plot(self,fig,ax,dim=0,alpha=0.4,every=1,ax_label=True,norm=True):
        """
        Args:
            ax axes object
            fig Figure 
            dim is the dimension, for example:in a 3D vector, 0-x, 1-y, 2-z and 3-total.
            alpha is the transparency of the filling
            every to not have so many points
            norm True if normalised
            The axis label is given here but it could be renamed later
        """
        if norm==True:
            cor=self.cor_norm[dim]
        else:
            cor=self.cor[dim]
            
        
        y=np.array([i.n for i in cor])
        y_error=np.array([i.s for i in cor])
        ax.plot(self.times[::every],y[::every])
        ax.fill_between(self.times, y-y_error, y+y_error ,alpha=0.4)
        
        # It mostly means that the plot will not be further used, for example not used to compare correlations
        if ax_label==True:
            ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
            ax.set_ylabel(r'$\langle %s(t)%s(0) \rangle$'%(self.flux1_name,self.flux2_name))
            ax.set_xlabel('time')

        return fig,ax
    
    def transport_coeff(self, T, V, xmin, xmax):
        """
        T is the temperature
        V is the system volume
        xmin
        """
    
        y = np.array([i.n for i in self.cor[0]])
        x = self.times
        I = cf.integrate(x,y,xmin,xmax)
        self.coeff= V/T*I
        
        return self.coeff
    




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
    total_flux=flux(data1[:,[1,3,5]],times,"Q")
    solute_excess=flux(data1[:,[2,4,6]],times,"J_s-c_s^BQ")
    
    
    
    c11=correlation(total_flux,total_flux,max_delta)
    c11.evaluate()

    
    c12=correlation(total_flux,solute_excess,max_delta)
    c12.evaluate()

    
    c21=correlation(solute_excess,total_flux,max_delta)
    c21.evaluate()

    
    c22=correlation(solute_excess,solute_excess,max_delta)
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

print "\nFinished directories %s\n"%finished_directories
print "\nUnfinished correlations in %s\n"%unfinished_correlation 


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
        print "Now in folder %s\n"%folder
        
        if folder in unfinished_correlation[0]:
            print "Running the correlation analysis in %s\n"%folder
            c11,c12,c21,c22=run_correlation_analysis(folder,input_file, delta_t)
        
        else:
            print "Loading the results in %s\n"%folder
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
    c11=bundle_correlation(c11_total,times,"Q","Q")
    c12=bundle_correlation(c12_total,times,"Q","J_s-c_s^BQ")
    c21=bundle_correlation(c21_total,times,"J_s-c_s^BQ","Q")
    c22=bundle_correlation(c22_total,times,"J_s-c_s^BQ","J_s-c_s^BQ")
    
    #TODO this could be done iterating over all the instances
    c11.save('c11')
    c12.save('c12')
    c21.save('c21')
    c22.save('c22')

else:
    c11 = cf.load_instance("c11.pkl")
    c22 = cf.load_instance("c22.pkl")
    c12 = cf.load_instance("c12.pkl")
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

print "The c11 is %s\n" %c11.transport_coeff(1,V,0,xmax)
print "The c12 is %s\n" %c12.transport_coeff(1,V,0,xmax)
print "The c21 is %s\n" %c21.transport_coeff(1,V,0,xmax)
print "The c22 is %s\n" %c22.transport_coeff(1,V,0,xmax)




