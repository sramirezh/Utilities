#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 12:00:01 2018
Script to plot all the results from different monomer number,
in the options you can decide to plot all the dat files or just the ones you decide

Still have to add the interactions by hand in the interactions array
@author: sr802
"""

import argparse
import pandas as pd
import numpy as np
import warnings
import sys
import os
import glob
import bisect
from scipy import optimize

warnings.filterwarnings("ignore")


sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

import matplotlib.pyplot as plt


"""
*******************************************************************************
Functions
*******************************************************************************
"""



def plot_results(all_data, interactions, theory = False, out_results = False):
    """
    return the plot, and if asked the an object with all the result.
    
    Args:
    
    Returns:
        results is an object that has for each epsilon, sigma those values in 
        the first entry and in the second, the lenght, mobility and the error
    """
    cf.set_plot_appearance()
    fig,ax=plt.subplots()
    
    results = []
    
    for j,interaction in enumerate(interactions):
        mobility=[]
        error_mobility=[]
        one_set = []
        for i,ave_data in enumerate(all_data):
            mobility.append(ave_data[j,0])
            error_mobility.append(ave_data[j,1])
        
        if theory == True:
            mobility=np.abs(mobility)
        
        
        if j<len(colors):color = colors[j]
        else: color=np.random.rand(3)
    
#        x=np.array(lengths).astype(np.float)
#        y=np.array(mobility).astype(np.float)
#
#    
#        yerror=error_mobility
        #pinit = [1.0,-1.0]
        #out = optimize.leastsq(errfunc, pinit, args=(x, y, yerror), full_output=1)
        #cov=out[1] #Covariance in the
        #pfinal = out[0] #fitting coefficients
        #logger.info "for %s The slope is %f error is %f" %(interaction,pfinal,np.sqrt(cov))
    
        epsilon=float(cf.extract_digits(interaction)[0])
        sigma=float(cf.extract_digits(interaction)[1])
    
#        #Plotting the fitting curve
#        if epsilon==1.0 and sigma==1.0:
#            ax.plot(np.unique(x),np.zeros(len(np.unique(x))),color=colors[j],linestyle='--')
#        else:
#            ax.plot(np.unique(x),fitfunc(pfinal,np.unique(x)),color=colors[j],linestyle='--')
#
#        color=ax.lines[-1].get_color() #Color of the last line ploted, it takes each point in error bar a a different line
        
        
        #To use the results later
        one_set.append([epsilon,sigma])
        one_set.append(np.column_stack((lengths, mobility, error_mobility)))
        
        results.append(one_set)
        ax.errorbar(lengths,mobility,yerr=error_mobility,label=interaction, color=color, fmt='o')
    
        """
        Prints the fitting factors and their errors
        """
    
        logger.info("For epsilon=%s and sigma=%s" %(epsilon,sigma))
#        logger.info("The slope is %f and the error is %f" %(pfinal[1],np.sqrt(cov[1,1])))
#        logger.info("The intercept is %f and the error is %f" %(pfinal[0],np.sqrt(cov[0,0])))
        
    if out_results == True:
        return fig, ax, results
    else:
        return fig, ax
    
        




"""
*******************************************************************************
Main
*******************************************************************************
"""
parser = argparse.ArgumentParser(description='This script reads the Results.dat ' \
                                 'from several simualtions and plots the mobility'\
                                 'vs N, as long as they are named "*Results.dat"',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-file_name', metavar='InputFile',help='Input filename',nargs='+',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-theory', metavar='theory',help='Theoretical filename, containing, N U_0 U',type=bool,default=False)
args = parser.parse_args()
files=args.file_name

# Creating the logger
logger = cf.log(__file__, os.getcwd())   

logger.info("Using the following arguments for the paser")
logger.info(args)



if files==None:
    dat_files=glob.glob('*Results.dat') #Set this as the option for no input.
else:
    dat_files=files



#Sorting with the number of polymers
all_data=[] #Every position keeps the data for one number of particles
data_pd=[]
lengths=[]

for f in dat_files:
    """
    Inserts the data from the files in order based on the length of the polymer.
    """
    logger.info("reading file %s \n"%f)
    length=int(cf.extract_digits(f)[0])
    position=bisect.bisect(lengths,length)
    lengths.insert(position,length)
    data=pd.read_csv(f,sep=" ").values
    all_data.insert(position,np.array(data[:,1:],dtype=float))
    data_pd.insert(position,data)


colors=['r','b','k','g']
"""
###############################################################################
Starting the plot
###############################################################################
"""
plt.close('all')

interactions = [ r'$\varepsilon_{ms}=0.5\, \sigma_{ms}=1.0 $',r'$\varepsilon_{ms}=1.5 \, \sigma_{ms}=1.0 $']
#interactions=[r'$\epsilon_{ms}=0.5\, \sigma_{ms}=1.0 $',r'$\epsilon_{ms}=1.0 \,\sigma_{ms}=1.0 $',r'$\epsilon_{ms}=1.5 \, \sigma_{ms}=1.0 $']

fig,ax = plot_results(all_data, interactions, theory = False, out_results =False)
fig2,ax2,results = plot_results(all_data, interactions, theory = True, out_results =True) #Including theoretical results



# =============================================================================
# Adding the theoretical results
# =============================================================================
if args.theory==True:
    L=ax2.legend()
    
    theory_ads=cf.read_data_file("Theoretical_1.5.dat").values
    theory_dep=cf.read_data_file("Theoretical_0.5.dat").values
#    ax2.scatter(theory_ads[:,0],theory_ads[:,1],marker='v',color="blue",label=r'$R_h^K$')
    ax2.scatter(theory_dep[:,0],np.abs(theory_dep[:,2]),marker='x',color="red",label=' '+r'$\epsilon=0.5$')
    ax2.scatter(theory_ads[:,0],theory_ads[:,2],marker='x',color="blue",label=r'$R_h^{\lambda}$')
#    ax2.scatter(theory_dep[:,0],np.abs(theory_dep[:,1]),marker='v',color="red",label=r'$R_h^K$')
    
    
    




"""Axis"""

xoffset=0.05
yoffset=0.8
# =============================================================================
# For the simulation results
# =============================================================================
ax.set_xlabel(r'$N_m $')
ax.grid(False)
ax.set_ylabel(r'$M_{ps} [\tau/m]$')

ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
#ax.axvline(x=0, ymin=0, ymax=1,ls=':',c='black')
xmin,xmax=plt.xlim()
deltax=xmax-xmin

ax.set_xticks(np.arange(len(lengths)+1)*30)
ax.set_xlim(0,70)


ymin,ymax=ax.get_ylim()
deltay=ymax-ymin


ax.set_ylim(ymin,ymax+deltay*yoffset)
ax.set_yticks(np.arange(-0.2,0.4,0.1))
#ax.set_xlim(xmin-deltax*xoffset,xmax+deltax*xoffset)

ax.legend(loc='upper left',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
           frameon=True, fancybox=False, edgecolor='k',ncol=1)

# =============================================================================
# For the comparison
# =============================================================================
xoffset=0.05
y2offset=1

ax2.set_xlabel(r'$N_m $')
ax2.grid(False)
ax2.set_ylabel(r'$\ln |M_{ps}|$')

ax2.set_yscale('log')
ax2.set_ylim(0,1)
ax2.set_xlim(0,70)

ymin,ymax=ax2.get_ylim()
deltay=ymax-ymin

ax2.set_ylim(ymin,ymax+deltay*y2offset)




#The original legend
#ax2.legend(loc='upper left',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],frameon=True, fancybox=False, edgecolor='k',ncol=2)

#Rewriting the legend
ax2.legend(['Theo '+r'$\varepsilon=0.5$','Theo '+r'$\varepsilon=1.5$','Sim '+r'$\varepsilon=0.5$','Sim '+r'$\varepsilon=1.5$'],loc='upper left',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
           frameon=True, fancybox=False, edgecolor='k',ncol=2)



"""General"""
fig.tight_layout()
fig.savefig("Mobility_N.pdf",transparent=True)
fig2.tight_layout()
fig2.savefig("Theoretical_comparison.pdf",transparent=True)
plt.close()




# =============================================================================
# Results including the slip correction
# =============================================================================


def fitfunc(p,x):
    """
    This gets the value of the theory and fit them to the simulation, with -p being the fitting parameter
    """
    index = np.where(x_vect == x)[0]
    y = theory_vect[index]*p   # theory - fitting parameter that is a displacement

    return y
    
    
    
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err #To include the error in the least squares


theory = [theory_dep,theory_ads]


if args.theory==True:
    
    theory_ind = 2 # 1 for Kirkwood, 2 lambda
    logger.info("Defined the theory ind as %s,1 for Kirkwood, 2 lambda"%theory_ind)
    
    fig3, ax3 = plt.subplots()

    color = ["blue","red"]
    color2 = ["black","green"]
    for i,interaction in enumerate(results):
        logger.info(interaction,color[i])
    
        ax3.errorbar(interaction[1][:,0],interaction[1][:,1],yerr=interaction[1][:,2],label=r'Simulation $\varepsilon_{ms}=%s$'%interaction[0][0], color=color[i], fmt='o')
        
        

        

        
        # Insert for the theoretical results
        left, bottom, width, height = [0.55, 0.25, 0.4, 0.25]
        ax4 = fig3.add_axes([left, bottom, width, height])
        ax4.plot(theory[i][:,0],np.abs(theory[i][:,theory_ind]),marker = 'v', ls = '--' , color=color[i],label=r'$R_h^K$')
        ymin, ymax = plt.ylim()
        ax4.set_ylim(0, ymax*1.1)
        ax4.tick_params(axis='both', which='major', labelsize=14)
    

        
        
        
        
        #Finding the slip length
        # naming things, reducing to the simulations that have equivalent in the theory
        indexes = [i for i, e in enumerate(interaction[1][:,0]) if e in theory[1][:,0]]
        
        
        x_vect = theory[i][:,0]
        theory_vect = np.abs(theory[i][:,theory_ind])
        simulation_vect = interaction[1][indexes,1]
        error_vect = interaction[1][indexes,2]
        
        pinit = 0
        out = optimize.leastsq(errfunc, pinit, args=(x_vect, simulation_vect, error_vect), full_output=1)
        #cov=out[1] #Covariance in the
        pfinal = out[0] #fitting coefficients
        
        logger.info (pfinal[0])
        
        ax3.plot(theory[i][:,0],pfinal[0]*np.abs(theory[i][:,theory_ind]),marker = 's', ls = '--', color=color[i],
                 label=r'$1+b/L= %1.2f$'%pfinal[0])
    ymin,ymax=ax3.get_ylim()
    ax3.set_ylim(0,ymax*1.05)
    ax3.legend(loc='upper left',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
           frameon=True, fancybox=False, edgecolor='k',ncol=2, fontsize = 10)
    
    ax3.set_xlabel(r'$N_m $')
    ax3.grid(False)
    ax3.set_ylabel(r'$|M_{ps}|$')
    fig3.tight_layout()
    
    
    fig3.savefig("Theoretical_comparison_improved.pdf",transparent=True)


#    L=ax2.legend()
#    
#    theory_ads=cf.read_data_file("Theoretical_1.5.dat").values
#    theory_dep=cf.read_data_file("Theoretical_0.5.dat").values
##    ax2.scatter(theory_ads[:,0],theory_ads[:,1],marker='v',color="blue",label=r'$R_h^K$')
#    ax2.scatter(theory_dep[:,0],np.abs(theory_dep[:,2]),marker='x',color="red",label=' '+r'$\epsilon=0.5$')
#    ax2.scatter(theory_ads[:,0],theory_ads[:,2],marker='x',color="blue",label=r'$R_h^{\lambda}$')
##    ax2.scatter(theory_dep[:,0],np.abs(theory_dep[:,1]),marker='v',color="red",label=r'$R_h^K$')



#"""
#Mobility vs N for same interaction
#"""
#fig,ax=plt.subplots()
#n_interactions=len(data)
#for i in xrange(n_interactions):
#    label=data_pd[0][i][0]
#    y=[]
#    x=names
#    for j in xrange(len(x)):
#        y.append(data_pd[j][i][1])
#    ax.scatter(x,y,label=label)
#    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
#    i=i+1
#
#ax.set_xlabel(r'$N $',fontsize=axis_font)
##ax.grid()
#ax.set_ylabel(r'$ b/R_g [\tau/m\sigma]$',fontsize=axis_font)
#
#ax.tick_params(labelsize=tick_font)
#ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
#ax.axvline(x=0, ymin=0, ymax=1,ls=':',c='black')
#plt.tight_layout()
#ax.legend(fontsize=legend_font)
#fig.savefig("Mobility_vs_N.pdf",transparent=True)
#plt.close()
