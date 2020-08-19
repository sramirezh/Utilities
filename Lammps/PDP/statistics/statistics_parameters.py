    #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script gets the results created by dp_poly and the averages of vdata.dat
and computes relevant quantities and generates plots, It has to be run inside every N_X

The parameters that have error measurements are taken as ufloats
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
from scipy import optimize
import glob
import warnings
import seaborn as sns
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import compute_statistics as cstat

try:
    from uncertainties import ufloat,unumpy
except ImportError as err2:
    print(err2)


"""
*******************************************************************************
Functions
*******************************************************************************
"""

def number_properties(lines):
    """
    Gets the number of properties based on the pattern in Statistics_summary.dat
    Args:
        lines are the list of lines in the file.
    Returns:
        nproperties are the number of properties analysed

    """
    indexes=cf.parameter_finder(lines, "dDP")
    nproperties=indexes[1]-indexes[0]-2
    return nproperties

def build_data():
    """
    Function to initialise all the elements of the class
    """

    interactions=[]
    with open("Statistic_summary.dat", 'r') as f:
      lines = f.readlines()


    #Finding the number of properties
    indexes=cf.parameter_finder(lines,"dDP")
    nproperties=indexes[1]-indexes[0]-2
    print("\nFound %d properties in Statistic_summary.dat\n" %nproperties)

    i=0
    count=0

    while i<len(lines):
        if re.search("E_*",lines[i] ): #Finding the LJ parameters
            interactions.append(LJInteraction(cf.extract_digits(lines[i])[-2::]))
            print("\nReading data from  %s"%lines[i])
            i=i+1
            count+=1
        else:
            print(lines[i])
            if re.search("\AdDP*",lines[i] ):
                interactions[-1].addforce(lines[i])
                i+=1
                properties=[]
                for j in range(nproperties):
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
    if length>1: print("There were several ocurrences")
    if length==0: print("No ocurrences found")
    return indexes


def ncols(nparameters, row_per_column):
    """
    Returns the ideal number of columns for the desired number of rows per column
    """
    ncols=nparameters/row_per_column

    return ncols

# define our (line) fitting function



def plot_force_individuals(interactions):
    """
    Plots the parameters from statistic_summary for each force, for each interaction
    """


    #General plot parameters
    axis_font=24
    tick_font=20
    legend_font=18
    xoffset=0.1
    yoffset=0.1
    error_cap=4


    colors=['r','b','k','g', sns.xkcd_rgb['tangerine'], sns.xkcd_rgb['dark red']]
    if len(interactions)>len(colors):
        palette=sns.color_palette(palette= "bright" , n_colors = len(interactions))

    #This Dict is going to be compared with the variable file_name
    dic_yaxis={'conc_bulk':r'$c_s^B [\sigma^{-3}]$','vx_poly':r'$v_p^x[\sigma/\tau]$','r_gyration':r'$R_g [\sigma]$','rRg2':r'$R_{g}^2 [\sigma^2]$'}
    dic_fit={'vx_poly':1}
    print("\nGenerating Plots...")
    directory="plots/individual"

    if not os.path.exists(directory):
        os.makedirs(directory)

    n_properties=len(interactions[0].properties[0]) #Number of properties

    has_error=np.ones((n_properties), dtype=bool)

    for property_index in range(n_properties):
        if np.size(interactions[0].properties[0][property_index])==1:
            has_error[property_index]=False

        prop_name=interactions[-1].property_names[0][property_index] #Crude property name
        if "Time" in prop_name: continue #To avoid plotting the timestep
        file_name=re.sub('^_|^v_|^c_',"",prop_name).strip('_')
        name=re.sub('_',' ',file_name)
        print("\nplotting the %s" %name)


        interactions.sort() #sorts using the method __lt__

        fig,ax=plt.subplots()
        for i,ljpair in enumerate(interactions):
            yvalue=np.empty(0)
            yerror=np.empty(0)
            force_list=[]
            for n,force in enumerate(ljpair.forces):

                if has_error[property_index]==True:
                    yvalue=np.append(yvalue,ljpair.properties[n][property_index][0]) #Be careful np.append inserts in reverse order compared to .append
                    yerror=np.append(yerror,ljpair.properties[n][property_index][1]) #Taking the autocorrelation error
                else:
                    yerror=None
                    yvalue=np.append(yvalue,ljpair.properties[n][property_index])
                force_list.append(force)


            #Defining the first colors from array and the rest by random numbers
            if i<len(colors):color=colors[i]
            else:
                color=palette[i]


            plt.errorbar(force_list,yvalue,yerr=yerror,xerr=None,fmt='o',label=r'$\varepsilon_{ms}$=%s $\sigma_{ms}$=%s '%(ljpair.epsilon,ljpair.sigma),
                         color=color,capsize=error_cap)

            """
            Linear fit
            """

            fitfunc1 = lambda p, x: p[0] * x  #Fitting to a line that goes through the origin
            fitfunc2 = lambda p, x: p[0] * x + p[1] #Fitting to a line
            errfunc1 = lambda p, x, y, err: (y - fitfunc1(p, x)) / err #To include the error in the least squares
            errfunc2 = lambda p, x, y, err: (y - fitfunc2(p, x)) / err

            try:
                x=np.array(force_list)
                y=yvalue
                #The velocity is the only one that is forced to go through the origin
                if file_name=="vx_poly":
                    pinit=[1.0]
                    out = optimize.leastsq(errfunc1, pinit, args=(x, y, yerror), full_output=1)
                    pfinal = out[0] #fitting coefficients
                    x=np.insert(x,0,0)
                    if ljpair.epsilon==1.0 and ljpair.sigma==1.0:
                        ax.plot(np.unique(x),np.zeros(len(np.unique(x))),color=color,linestyle='--')
                    else:
                        ax.plot(np.unique(x),fitfunc1(pfinal,np.unique(x)),color=color,linestyle='--')
                else:
                    # This is for not vx_poly"
                    pinit=[1.0,-1.0]
                    out = optimize.leastsq(errfunc2, pinit, args=(x, y, yerror), full_output=1)
                    pfinal = out[0] #fitting coefficients
                    ax.plot(np.unique(x),fitfunc2(pfinal,np.unique(x)),color=color,linestyle='--')
                cov=out[1] #Covariance

                print("for Epsilon=%s and Sigma =%s The slope is %f error is %f" %(ljpair.epsilon,ljpair.sigma,pfinal[0],np.sqrt(cov[0][0])))

            except:
                pass

        file_name=name.replace(" ","_")


        """Legend"""
        ncols=int(np.ceil(len(interactions)/4))
#        plt.legend(fontsize=legend_font,loc='upper left',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
#           frameon=True, fancybox=False, edgecolor='k')
        plt.legend(fontsize=legend_font/ncols,loc='upper left',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
           frameon=True, fancybox=False, edgecolor='k',ncol=ncols)



        """Axis"""
        try:
            ylabel=dic_yaxis[file_name]
            ax.set_ylabel(ylabel,fontsize=axis_font)
        except:
            ylabel=file_name

        ax.set_xlabel(r'$F_{s}^{\mu}=-\nabla_x \mu_s [\varepsilon/\sigma]$',fontsize=axis_font)
        ax.tick_params(labelsize=tick_font,direction='in',top=True, right=True)
        ylabel=file_name

        ymin,ymax=plt.ylim()
        deltay=ymax-ymin
        ax.set_ylim(ymin-deltay*yoffset,ymax+deltay*0.45)

        xmin,xmax=plt.xlim()
        deltax=xmax-xmin
        ax.set_xlim(0,0.12)



        plt.xticks(np.arange(0.02,0.12,0.02))
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)


        """Lines"""
        if ymin*ymax<0:
            ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')

        """General"""

        plt.grid(False)
        try:
            plt.rcParams["mathtext.fontset"] = "cm"
            plt.rcParams["text.usetex"] =True
        except:
            pass
        plt.tight_layout()
        plt.savefig("plots/individual/%s.pdf"%file_name,transparent=True)
        plt.close()
    print("\nGenerated plots for the individual properties vs forces, find them in '%s' " %directory)


def is_valid_file(parser,arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!, see source options" % arg)


def compute_statistics_param(dpolymin):
    """
    Gets the parameters for the dpolymin that should be used to call the compute_statistics.sh
    It reads them from the firs input.lmp that finds.

    Parameters
    ----------
    dpolymin :  int
                Number of samples to be discarded in DPpoly

    Returns
    -------
    s : int
        Minimum time step for the dp_poly analysis
    d : int
        time step between samples
    n : int
        Total number of steps to analyse

    """
    tfile,err=cf.bash_command("""find . -name "*.lmp" -path "*/dDP*" -print -quit""")#Assuming all the input files have the same parameters.
    print(tfile)
    out,err=cf.bash_command("""grep -m 1 "myDump equal" %s"""%tfile)
    d=int(cf.extract_digits(out)[0])

    try: #For the old version of the simulations
        out2,err=cf.bash_command("""grep -m 1 "myStepsEach equal"  %s"""%tfile)
        n1=int(cf.extract_digits(out2)[0])
        out3,err=cf.bash_command("""grep -m 1 "myLoop loop"  %s"""%tfile)
        n2=int(cf.extract_digits(out3)[0])
        n = n1 * n2
    except: #New version without loop
        out4,err = cf.bash_command("""grep -m 1 "myRun equal" %s"""%tfile)
        n = n1=int(cf.extract_digits(out4)[0])

    s=d*dpolymin
    return [s,d,n]




def plot_parameter_vs_epsilon(y_label,name_key,file_name):
    """
    Plots any parameter in ave_data as a function of the epsilon
    Args:
        y_label is a latex type expression for the y axis label for example r'$\Gamma_{ps} [\tau/m]$'
        name_key is a key string that could be found in the column names to identify the column number.
        file name is the name of the including the extension.

    """

    axis_font=24
    tick_font=20
    xoffset=0.16
    yoffset=0.1
    error_cap=4
    ave_data_index=cf.parameter_finder(column_names,name_key)[0]-1
    epsilon_vect=[]
    for i in range(len(interactions)):
        epsilon_vect.append(interactions[i].epsilon)

    fig,ax=plt.subplots()
    ax.errorbar(epsilon_vect,ave_data[:,ave_data_index],yerr=ave_data[:,ave_data_index+1],fmt='o',capsize=error_cap,color='b')

#    x=np.array(epsilon_vect)
#    y=np.array(ave_data[:,0])


    """Axis"""
    ax.set_xlabel(r'$\varepsilon_{ms} $',fontsize=axis_font)
    ax.grid(False)
    ax.set_ylabel(y_label,fontsize=axis_font)
    ax.tick_params(labelsize=tick_font,direction='in',top=True, right=True)

    ymin,ymax=plt.ylim()
    deltay=ymax-ymin
    ax.set_ylim(ymin-deltay*yoffset,ymax+deltay*yoffset)

    xmin,xmax=plt.xlim()
    #ax.set_xlim(0,xmax+deltax*xoffset)
    ax.set_xlim(0,xmax+1)
    plt.xticks(np.arange(0,xmax+1,1))




    """Lines"""
    ax.axhline(y=0, xmin=0, xmax=1,ls='--',c='black')

    """General"""
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["text.usetex"] =True
    plt.tight_layout()
    fig.savefig("plots/all/%s"%file_name, transparent=True)
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
    total = 0 #This is a class atribute

    @staticmethod  #Does not have self as is intended to be used by the class and not an object.
    def number_interactions():
        "Return the number of instances created in the class"
        return LJInteraction.total


    def __init__(self,lj_parameters):
        self.epsilon=float(lj_parameters[0])
        self.sigma=float(lj_parameters[1])
        self.forces=[]
        self.properties=[]
        self.property_names=[]
        LJInteraction.total += 1

    def __str__(self):
        self.name=("E=%s S=%s" %(self.epsilon,self.sigma))
        return self.name

    def __lt__(self,other):
        """
        In order to sort the interactions by the epsilon
        """
        return self.epsilon < other.epsilon


    def addforce(self,new_force):
        force=float(new_force.strip('\n/dDP'))/1000
        self.forces.append(force)

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


    def compute_mobility(self):
        self.mobility=[]
        self.mob_rg=[] #Mobility over Rg
        count=0
        index_vx=parameter_finder(self.property_names[count],"vx_relative")[0]
        index_rg=parameter_finder(self.property_names[count],"r_gyration")[0]
        for force in self.forces:
            if force!=0:
                velocity=ufloat(self.properties[count][index_vx][0],self.properties[count][index_vx][1])
                rg=ufloat(self.properties[count][index_rg][0],self.properties[count][index_rg][1])
                mobility=-velocity/force
                self.mobility.append(mobility)
                if rg.n==0:
                   self.mob_rg.append(10**8) #To avoid division by 0
                else:
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
    parser = argparse.ArgumentParser(description="This script gets the results created by dp_poly and the averages of vdata.dat, computes relevant quantities and generates plots, It has to be run inside every N_X",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-folder_names',help='Folders where the analysis is going to run i.e E_*',default=glob.glob('E_*'),nargs='+')
    parser.add_argument('-s','--source',choices=['read','run','gather'], default="read",help='Decides if the if the file Statistics_summary.dat needs to be read, run, gather  ')
    parser.add_argument('--vdatamin', help='Number of samples to be discarded in vdata.dat', default=1000, type=int)
    parser.add_argument('--dpolymin', help='Number of samples to be discarded in DPpoly', default=100, type=int)
    args = parser.parse_args()
    directories=args.folder_names
    source=args.source

if source == "read":
    is_valid_file(parser,"Statistic_summary.dat")


if source=="run":
    print("\nRunning the statistics analysis, using the following parameters")
    dppoly_params=compute_statistics_param(args.dpolymin)
    print("Initial dp_poly step=%d"%dppoly_params[0])
    print("Interval dp_poly=%d"%dppoly_params[1])
    print("Final dp_poly step=%d"%dppoly_params[2])
    print("vdata discarded steps =%d"%args.vdatamin)
    print(" ")
    cstat.compute_statistics(directories,args.vdatamin)


elif source=="gather":
    print("\nGathering the statistics analysis results")
    cstat.gather_statistics(directories)


"""
*******************************************************************************
Main program
*******************************************************************************
"""
print("\nAnalizing the results")
interactions=build_data()
plot_force_individuals(interactions)

"""
*******************************************************************************
Building the averaged data, excluding force equal=0
*******************************************************************************
"""

def average_uncertainties(A):
    """
    Returns the average of an array with uncertainties
    Args:
        A is a nx2 array where the first colum is the mean, the second contains the std deviation

    Returns:
        ave is the average as a ufloat
    """
    array=unumpy.uarray(A[:,0],A[:,1])
    ave=sum(array)/len(array)

    return ave



ave_data=[]
for interaction in interactions:
    name='E_%s_S_%s '%(interaction.epsilon,interaction.sigma)
    interaction.compute_mobility()

    mobilities=np.array((interaction.mobility))
    ave_mobility=sum(mobilities)/len(mobilities)

    #For non uarrays
    ave_concentration_bulk=average_uncertainties(np.array(interaction.get_property("conc_bulk"))[:,:2]) #Solute concentration in the bulk
    ave_rg=average_uncertainties(np.array(interaction.get_property("r_gyration"))[:,:2])
    ave_rhyd=average_uncertainties(np.array(interaction.get_property("r_hyd"))[:,:2])

    #For uarrays
    mobility_rg=sum(interaction.mob_rg)/len(interaction.mob_rg) #Using the properties of uarrays

    data_interaction=[name,ave_mobility.n,ave_mobility.s, ave_concentration_bulk.n, ave_concentration_bulk.s,ave_rg.n,ave_rg.s, ave_rhyd.n,ave_rhyd.s ] #the n=nominal, s standard deviation from ufloat
    column_names=['LJ_interaction','ave_mobility','mobility_error','ave_concentration_bulk', 'concentration_bulk_error','ave_rg','rg_error','ave_rhyd','rhyd_error']
    ave_data.append(data_interaction)

ave_data=np.array(ave_data)
pd_data=pd.DataFrame(ave_data,columns=column_names)


pd_data.to_csv("Results.dat",sep=' ',index=False)


"""
###############################################################################
Starting the plot
###############################################################################
"""


directory="plots/all"
if not os.path.exists(directory):
    os.makedirs(directory)


ave_data=np.array(ave_data[:,1::],dtype=float) #Avoiding the first column which contains the interactions.


#Mobility vs epsilon
plot_parameter_vs_epsilon(r'$\Gamma_{ps} [\tau/m]$','mobility','Mobility_vs_epsilon.pdf')

#Rg vs epsilon ms
plot_parameter_vs_epsilon(r'$R_g [\sigma]$','rg','R_g_vs_epsilon.pdf')

#Rhyd vs epsilon ms
plot_parameter_vs_epsilon(r'$R_{hyd} [\sigma]$','rhyd','R_hyd_vs_epsilon.pdf')

print("\nGenerated average results Results.dat and plots in '%s'"%directory)
