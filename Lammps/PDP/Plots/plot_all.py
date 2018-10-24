#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 12:00:01 2018
Script to plot all the results from different polymer number,
in the options you can decide to plot all the dat files or just the ones you decide
@author: sr802
"""
from __future__ import division
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


try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
except ImportError as err:
    print err

"""
*******************************************************************************
Functions
*******************************************************************************
"""

"""
Linear fit
"""

fitfunc = lambda p, x: p[1] * x+p[0]
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err #To include the error in the least squares


"""
*******************************************************************************
Main
*******************************************************************************
"""
parser = argparse.ArgumentParser(description='This script reads the Results.dat ' \
                                 'from several simualtions and plots the mobility'\
                                 'vs N, as long as they are named "*Results.dat"',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-file_name', metavar='InputFile',help='Input filename',nargs='+',type=lambda x: cf.is_valid_file(parser, x))
args = parser.parse_args()
files=args.file_name


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
    print "reading file %s \n"%f
    length=int(cf.extract_digits(f)[0])
    position=bisect.bisect(lengths,length)
    lengths.insert(position,length)
    data=pd.read_csv(f,sep=" ").as_matrix()
    all_data.insert(position,np.array(data[:,1:],dtype=float))
    data_pd.insert(position,data)


colors=['r','b','k','g']
"""
###############################################################################
Starting the plot
###############################################################################
"""

axis_font=24
tick_font=20
legend_font=18
xoffset=0.05
yoffset=0.8
error_cap=4

"""
Mobility vs N
"""
fig,ax=plt.subplots()
interactions=[r'$\epsilon_{ms}=0.5\, \sigma_{ms}=1.0 $',r'$\epsilon_{ms}=1.0 \,\sigma_{ms}=1.0 $',r'$\epsilon_{ms}=1.5 \, \sigma_{ms}=1.0 $']

for j,interaction in enumerate(interactions):
    mobility=[]
    error_mobility=[]
    for i,ave_data in enumerate(all_data):
        mobility.append(ave_data[j,0])
        error_mobility.append(ave_data[j,1])

    if j<len(colors):color=colors[j]
    else: color=np.random.rand(3)

    x=np.array(lengths).astype(np.float)
    y=np.array(mobility).astype(np.float)
    yerror=error_mobility
    pinit = [1.0,-1.0]
    out = optimize.leastsq(errfunc, pinit, args=(x, y, yerror), full_output=1)
    cov=out[1] #Covariance in the 
    pfinal = out[0] #fitting coefficients
    #print "for %s The slope is %f error is %f" %(interaction,pfinal,np.sqrt(cov))
    
    epsilon=float(cf.extract_digits(interaction)[0])
    sigma=float(cf.extract_digits(interaction)[1])
    
    if epsilon==1.0 and sigma==1.0:
        ax.plot(np.unique(x),np.zeros(len(np.unique(x))),color=colors[j],linestyle='--')
    else:    
        ax.plot(np.unique(x),fitfunc(pfinal,np.unique(x)),color=colors[j],linestyle='--')
    color=ax.lines[-1].get_color() #Color of the last line ploted, it takes each point in error bar a a different line
    ax.errorbar(lengths,mobility,yerr=error_mobility,label=interaction, color=color, fmt='o',capsize=error_cap)
    
    """
    Printing the fitting factors and their errors
    """
    
    print "For epsilon=%s and sigma=%s" %(epsilon,sigma)
    print "The slope is %f and the error is %f" %(pfinal[1],np.sqrt(cov[1,1]))
    print "The intercept is %f and the error is %f" %(pfinal[0],np.sqrt(cov[0,0]))
    
    




"""Axis"""
ax.set_xlabel(r'$N_m $',fontsize=axis_font)
ax.grid(False)
ax.set_ylabel(r'$\Gamma_{ps} [\tau/m]$',fontsize=axis_font)
ax.tick_params(labelsize=tick_font,direction='in',top=True, right=True)



ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
#ax.axvline(x=0, ymin=0, ymax=1,ls=':',c='black')



xmin,xmax=plt.xlim()
deltax=xmax-xmin

plt.xticks(np.arange(len(lengths)+1)*30)
ax.set_xlim(20,130)


ymin,ymax=plt.ylim()
deltay=ymax-ymin

ax.set_ylim(ymin,ymax+deltay*yoffset)
plt.yticks(np.arange(-0.1,0.2,0.1))
#ax.set_xlim(xmin-deltax*xoffset,xmax+deltax*xoffset)

"""Legend"""
plt.legend(fontsize=legend_font,loc='upper left',labelspacing=0.5,borderpad=0.4,scatteryoffsets=[0.6],
           frameon=True, fancybox=False, edgecolor='k')


"""General"""
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["text.usetex"] =True
plt.tight_layout()
fig.savefig("Mobility_N.pdf",transparent=False)
plt.close()



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
