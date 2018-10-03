#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 12:00:01 2018
Script to plot all the results from different polymer number
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


all_data=[] #Every position keeps the data for one number of particles
data_pd=[]

for file in dat_files:
    print "reading file %s \n"%file 
    data=pd.read_csv(file,sep=" ").as_matrix()
    all_data.append(np.array(data[:,1:],dtype=float))
    data_pd.append(data)
    
    
    
#data=pd.read_csv("N_30_Results.dat",sep=" ").as_matrix()
#data1=pd.read_csv("N_60_Results.dat",sep=" ").as_matrix()
#data2=pd.read_csv("N_90_Results.dat",sep=" ").as_matrix()
#data3=pd.read_csv("N_120_Results.dat",sep=" ").as_matrix()
##all_data.append(np.array(data0[:,1:],dtype=float))
#all_data.append(np.array(data[:,1:],dtype=float))
#all_data.append(np.array(data1[:,1:],dtype=float))
#all_data.append(np.array(data2[:,1:],dtype=float))
#all_data.append(np.array(data3[:,1:],dtype=float))
#
##data_pd.append(data0)
#data_pd.append(data)
#data_pd.append(data1)
#data_pd.append(data2)
#data_pd.append(data3)





#names=[0,30,60,90,120]
names=[30,60,90,120]
colors=['r','b','k']
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
    print j
    mobility=[]
    error_mobility=[]
    for i,ave_data in enumerate(all_data):
        mobility.append(ave_data[j,0])
        error_mobility.append(ave_data[j,1])
    #    x=np.array(ave_data[3,4]).astype(np.float)
    #    y=np.array(ave_data[3,0]).astype(np.float)
    #    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
    x=np.array(names).astype(np.float)
    y=np.array(mobility).astype(np.float)
    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),color=colors[j],linestyle='--')
    color=ax.lines[-1].get_color() #Color of the last line ploted, it takes each point in error bar a a different line
    ax.errorbar(names,mobility,yerr=error_mobility,label=interaction, color=color, fmt='o',capsize=error_cap)




"""Axis"""
ax.set_xlabel(r'$N_m $',fontsize=axis_font)
ax.grid(False)
ax.set_ylabel(r'$\Gamma_{ps} [\tau/m]$',fontsize=axis_font)
ax.tick_params(labelsize=tick_font,direction='in',top=True, right=True)



ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
#ax.axvline(x=0, ymin=0, ymax=1,ls=':',c='black')



xmin,xmax=plt.xlim()
deltax=xmax-xmin

plt.xticks(np.arange(len(names)+1)*30)
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

