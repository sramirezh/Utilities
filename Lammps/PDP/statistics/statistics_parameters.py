#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script gets the results created by dp_poly and the averages of vdata.dat
and computes relevant quantities and generates plots, It has to be run inside every N_X

Args:
    Input filen name
Returns:


@author: simon
"""

import os
import sys
import pandas as pd
import numpy as np
import re
import argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
from Lammps.linux import bash_command


try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
#    from matplotlib.backends.backend_pdf import PdfPages
except ImportError as err:
    print err

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
            print "\nReading data from  %s"%lines[i]
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
    print "\nGenerating Plots..."
    directory="plots/individual"
    if not os.path.exists(directory):
        os.makedirs(directory)

    n_properties=len(interactions[0].properties[0]) #Number of properties

    for property_index in xrange(n_properties):

        prop_name=interactions[-1].property_names[0][property_index] #Crude property name

        if "Time" in prop_name: continue #To avoid plotting the timestep

        plt.figure()

        for ljpair in interactions:
            n=0
            re2=[]
            force_list=[]
            for force in ljpair.forces:
                re2.append(ljpair.properties[n][property_index]) #end to end radious squared
                force_list.append(force)
                n+=1
            plt.plot(force_list,re2,'-o',label='$\epsilon$=%s $\sigma$=%s '%(ljpair.epsilon,ljpair.sigma))
            #plt.legend("" %(ljpair.epsilon,ljpair.sigma))

        file_name=re.sub('^_|^v_|^c_',"",prop_name).strip('_')
        name=re.sub('_',' ',file_name)

        print "\nplotting the %s" %name

        plt.title(name)
        file_name=name.replace(" ","_")
        plt.legend()
        plt.grid()
        plt.xlabel("Force")
        plt.savefig("plots/individual/%s.pdf"%file_name)
        plt.close()
    print "\nGenerated plots for the individual properties vs forces, find them in '%s' " %directory


def is_valid_file(parser,arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!, see source options" % arg)



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




"""
###############################################################################
Argument Parser
###############################################################################
"""

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script gets the results created by dp_poly and the averages of vdata.dat and computes relevant quantities and generates plots, It has to be run inside every N_X",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--source',choices=['READ','RUN','GATHER'], default="READ",help='Decides if the if the file Statistics_summary.dat needs to be READ, RUN, GATHERED  ')
args = parser.parse_args()
source=args.source

if source=="RUN":
    print "\nRunning the statistics analysis"
    bash_command("""bash %s/compute_statistics.sh"""%dir_path)
elif source=="GATHER":
    print "\nGathering the statistics analysis results"
    bash_command("""bash %s/gather_statistics.sh"""%dir_path)
else:
    is_valid_file(parser,"Statistic_summary.dat")










"""
*******************************************************************************
Main program
*******************************************************************************
"""

print "\nAnalizing the results"
interactions=build_data()
plot_force_individuals(interactions)




"""
*******************************************************************************
Building the averaged data, excluding force equal=0
*******************************************************************************
"""


ave_data=[]
for interaction in interactions:
    name='E_%s_S_%s '%(interaction.epsilon,interaction.sigma)
    interaction.compute_mobility()

    ave_mobility=np.average(interaction.mobility)
    ave_concentration_rg=np.average(interaction.get_property("concentration")[1:]) #Solute concentration inside rg,In order to exlude f=0
    ave_concentration_bulk=np.average(interaction.get_property("conc_bulk")[1:]) #Solute concentration in the bulk
    delta_cs=ave_concentration_rg-ave_concentration_bulk
    ave_rg=np.average(interaction.get_property("rg_ave")[1:]) #Average Rg
    mobility_rg=ave_mobility/ave_rg #Mobility divided by Rg


    data_interaction=[name,ave_mobility, ave_concentration_rg, ave_concentration_bulk, delta_cs, ave_rg, mobility_rg ]
    ave_data.append(data_interaction)

header_data="lj_interaction,ave_mobility, ave_concentration_rg, ave_concentration_bulk, delta_cs, ave_rg, mobility_rg"
ave_data=np.array(ave_data)
pd_data=pd.DataFrame(ave_data,columns=['LJ_interaction','ave_mobility', 'ave_concentration_rg','ave_concentration_bulk','delta_cs','ave_rg','mobility_rg'])


pd_data.to_csv("Results.dat",sep=' ',index=False)

"""
###############################################################################
Starting the plot
###############################################################################
"""


directory="plots/all"
if not os.path.exists(directory):
    os.makedirs(directory)

"""
Mobility vs Delta Cs
"""
fig,ax=plt.subplots()


ax.scatter(ave_data[:,4],ave_data[:,1])
x=np.array(ave_data[:,4]).astype(np.float)
y=np.array(ave_data[:,1]).astype(np.float)


for i in xrange(len(interactions)):
    txt="%.2lf,%.2lf"%(interactions[i].epsilon,interactions[i].sigma)
    ax.annotate(txt, (x[i]-0.006,y[i]+0.002))
ax.set_xlabel(r'$\Delta c_s [1/\sigma^3] $',fontsize=16)
#ax.grid()
ax.set_ylabel(r'$b [\tau/m]$',fontsize=16)

ax.tick_params(labelsize=14)

#ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))

ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
ax.axvline(x=0, ymin=0, ymax=1,ls=':',c='black')
plt.tight_layout()

fig.savefig("plots/all/Mobility_Delta_Cs.pdf")
plt.close()

"""
Mobility/Rg vs Delta Cs
"""

fig,ax=plt.subplots()
ax.scatter(ave_data[:,4],ave_data[:,-1])


x=np.array(ave_data[:,4]).astype(np.float)
y=np.array(ave_data[:,-1]).astype(np.float)

for i in xrange(len(interactions)):
    txt="%.2lf,%.2lf"%(interactions[i].epsilon,interactions[i].sigma)
    ax.annotate(txt, (x[i]-0.006,y[i]+0.002))
ax.set_xlabel(r'$\Delta c_s [1/\sigma^3]$',fontsize=16)
#ax.grid()
ax.set_ylabel(r'$ b/R_g [\tau/m\sigma]$',fontsize=16)

ax.tick_params(labelsize=14)

#Using np.unique(x) instead of x handles the case where x isn't sorted or has duplicate values.
#ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))

ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
ax.axvline(x=0, ymin=0, ymax=1,ls=':',c='black')
plt.tight_layout()

fig.savefig("plots/all/Mobility_rg_Delta_Cs.pdf")
plt.close()




pd_data.to_csv("plots/all/Results.dat",sep=' ',index=False)

print "\nGenerated average results Results.dat and plots in '%s'"%directory
