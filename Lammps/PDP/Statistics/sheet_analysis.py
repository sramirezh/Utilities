#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This scripts reads the polymerDp spreadsheet
Args:
    Input filen name
Returns:


@author: simon
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re


"""
*******************************************************************************
Functions
*******************************************************************************
"""

def build_data():
    """
    Function to initialise all the elements of the class
    """

    interactions=[]
    with open("Statistic_summary.dat", 'r') as f:
      lines = f.readlines()
    
    i=0
    nproperties=11 #The number of properties per force in the input file (TO IMPROVE)
    count=0
    
    while i<len(lines):
        if re.search("\AE_*",lines[i] ): #Finding the LJ parameters
            interactions.append(LJInteraction(re.findall(r"[-+]?\d*\.?\d+", lines[i])))
            print "Analysing %s"%lines[i]
            i=i+1
            count+=1
        else: 
            if re.search("\AdDP*",lines[i] ):
                interactions[-1].addforce(lines[i])
                i+=1
                properties=[]
                for j in xrange(nproperties):
                    properties.append(lines[i])
                    i+=1
                interactions[-1].addproperties(properties)
                    
            i+=1
        
    return interactions


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


def plot_force_individuals(interactions):
    """
    Plots the parameters from statistic_summary for each force, for each interaction
    """
    
    directory="plots/individual"
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    n_properties=len(interactions[0].properties[0]) #Number of properties
    
    for property_index in xrange(n_properties):
        plt.figure()
        for ljpair in interactions:
            n=0
            re2=[]
            force_list=[]
            for force in ljpair.forces:
                re2.append(ljpair.properties[n][property_index]) #end to end radious squared
                force_list.append(force)
                n+=1
            plt.plot(force_list,re2,label='$\epsilon$=%s $\sigma$=%s '%(ljpair.epsilon,ljpair.sigma))
            #plt.legend("" %(ljpair.epsilon,ljpair.sigma))
        name=interactions[-1].property_names[0][property_index].replace(" ","_")
        plt.title(name)
        plt.legend()
        plt.savefig("plots/individual/%s.pdf"%name)
        plt.close()



"""
*******************************************************************************
CLASS DEFINITION
*******************************************************************************
"""
class LJInteraction(object):
    """
    Every pair of LJ interactions has its own
    """
    
    def __init__(self,lj_parameters):
        self.epsilon=float(lj_parameters[0])
        self.sigma=float(lj_parameters[1])
        self.forces=[]
        self.properties=[]
        self.property_names=[]

    def addforce(self,new_force):
        force=float(new_force.strip('\n/dDP'))/1000
        self.forces.append(force)
        
    def addproperties(self,properties):
        values=[]
        names=[]
        for element in properties:
            name,value=element.strip("\n").split("=")
            name=name.replace(" ","_")
            values.append(float(value))
            names.append(name)
        self.properties.append(values)
        self.property_names.append(names)
            
        
    def compute_mobility(self):
        self.mobility=[]
        self.mob_rg=[] #Mobility over Rg
        count=0
        index_vx=parameter_finder(self.property_names[count],"vx_relative")[0]
        index_rg=parameter_finder(self.property_names[count],"rg_ave")[0]
        for force in self.forces:
            if force!=0:
                velocity=self.properties[count][index_vx]
                rg=self.properties[count][index_rg]
                mobility=velocity/force
                self.mobility.append(mobility)
                self.mob_rg.append(mobility/rg)
            count+=1
            
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
            
        
        




#import argparse
#
#
#
#parser = argparse.ArgumentParser(description='This script evaluates the average of a quantity')
#parser.add_argument('FileName', metavar='InputFile', type=str,
#                    help='Input filename')
#
#parser.add_argument('--min', help='Number of timesteps to be discarded', default=1000, type=int)
#
#
#args = parser.parse_args()
#min_limit=args.min
#InputFile=args.FileName
       
    
"""
*******************************************************************************
Main program
*******************************************************************************
"""
     
interactions=build_data()
#plot_force_individuals(interactions)




"""
*******************************************************************************
Building the averaged data, excluding force equal=0
*******************************************************************************
"""


ave_data=[]
for interaction in interactions:    
    name='E_%sS_%s '%(interaction.epsilon,interaction.sigma)
    interaction.compute_mobility()
    
    ave_mobility=np.average(interaction.mobility)
    ave_concentration_rg=np.average(interaction.get_property("concentration")[1:]) #Solute concentration inside rg,In order to exlude f=0
    ave_concentration_bulk=np.average(interaction.get_property("conc_bulk")[1:]) #Solute concentration inn the bulk
    delta_cs=ave_concentration_rg-ave_concentration_bulk
    ave_rg=np.average(interaction.get_property("rg_ave")[1:]) #Average Rg 
    mobility_rg=ave_mobility/ave_rg #Mobility divided by Rg
    
    
    data_interaction=[name,ave_mobility, ave_concentration_rg, ave_concentration_bulk, delta_cs, ave_rg, mobility_rg ]
    ave_data.append(data_interaction)

header_data="lj_interaction,ave_mobility, ave_concentration_rg, ave_concentration_bulk, delta_cs, ave_rg, mobility_rg"
ave_data=np.array(ave_data)
pd_data=pd.DataFrame(ave_data,columns=['LJ_interaction','ave_mobility', 'ave_concentration_rg','ave_concentration_bulk','delta_cs','ave_rg','mobility_rg'])    

pd_data.to_csv("Results.dat",sep=' ',index=False)




#"""
################################################################################
#Starting the plot
################################################################################
#"""
#

#ax=pd_data.plot.scatter(x="delta_cs",y="mobility_rg")
#ax.set_xlabel("Delta $C_s$ [$1/\sigma^3$]",fontsize=16)
#ax.set_ylabel("$\mu/R_g$",fontsize=16)
#plt.grid()
#plt.savefig("MobilityRg_Delta_Cs.pdf")


#directory="plots/all"
#if not os.path.exists(directory):
#    os.makedirs(directory)
#
#
#fig,ax=plt.subplots()
#ax.scatter(ave_data[:,4],ave_data[:,1])
#
#for i, txt in enumerate(ave_data[:,0]):
#    ax.annotate(txt, (ave_data[i,4],ave_data[i,1]))
#ax.set_xlabel("$\Delta c_s$ [$1/\sigma^3$] Solutes inside and outside",fontsize=16)
#ax.grid()
#ax.set_ylabel("b [t/m]",fontsize=16)
#fig.savefig("plots/all/Mobility_Delta_Cs.pdf")
#
#plt.close()






