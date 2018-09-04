#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 12:00:01 2018
Script to plot all the results, it is really rough yet, so need to be improved
@author: sr802
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



all_data=[] #Every position keeps the data for one number of particles
data_pd=[]
data=pd.read_csv("N_30_Results.dat",sep=" ").as_matrix()
data1=pd.read_csv("N_60_Results.dat",sep=" ").as_matrix()
data2=pd.read_csv("N_90_Results.dat",sep=" ").as_matrix()
data3=pd.read_csv("N_120_Results.dat",sep=" ").as_matrix()
all_data.append(np.array(data[:,1:],dtype=float))
all_data.append(np.array(data1[:,1:],dtype=float))
all_data.append(np.array(data2[:,1:],dtype=float))
all_data.append(np.array(data3[:,1:],dtype=float))

data_pd.append(data)
data_pd.append(data1)
data_pd.append(data2)
data_pd.append(data3)





names=[30,60,90,120]
"""
###############################################################################
Starting the plot
###############################################################################
"""
axis_font=24
tick_font=20
legend_font=18

"""
Mobility vs Delta Cs
"""
fig,ax=plt.subplots()
i=0
for ave_data in all_data:
    label=str("N=%d"%names[i])
    ax.scatter(ave_data[:,3],ave_data[:,0],label=label)
    x=np.array(ave_data[:,3]).astype(np.float)
    y=np.array(ave_data[:,0]).astype(np.float)
    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
    i=i+1

ax.set_xlabel(r'$\Delta c_s [1/\sigma^3] $',fontsize=axis_font)
#ax.grid()
ax.set_ylabel(r'$b [\tau/m]$',fontsize=axis_font)

ax.tick_params(labelsize=tick_font)
ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
ax.axvline(x=0, ymin=0, ymax=1,ls=':',c='black')

ax.legend(fontsize=legend_font)


plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["text.usetex"] =True
plt.tight_layout()
fig.savefig("Mobility_Delta_Cs_all.pdf")
plt.close()




"""
Mobility_rg vs Delta Cs
"""
fig,ax=plt.subplots()
i=0
for ave_data in all_data:
    label=str("N=%d"%names[i])
    ax.scatter(ave_data[:,3],ave_data[:,-1],label=label)
    x=np.array(ave_data[:,3]).astype(np.float)
    y=np.array(ave_data[:,-1]).astype(np.float)
    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
    i=i+1

ax.set_xlabel(r'$\Delta c_s [1/\sigma^3] $',fontsize=axis_font)
#ax.grid()
ax.set_ylabel(r'$ b/R_g [\tau/m\sigma]$',fontsize=axis_font)

ax.tick_params(labelsize=tick_font)
ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
ax.axvline(x=0, ymin=0, ymax=1,ls=':',c='black')
plt.tight_layout()
ax.legend(fontsize=legend_font)
fig.savefig("Mobility_rg_Delta_Cs_all.pdf")
plt.close()



#"""
#Mobility_rg vs N for same interaction
#"""
#fig,ax=plt.subplots()
#n_interactions=len(data)
#for i in xrange(n_interactions):
#    label=data_pd[0][i][0]
#    y=[]
#    x=names
#    for j in xrange(len(x)):
#        y.append(data_pd[j][i][-1])
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
#fig.savefig("Movility_rg_vs_N.pdf")
#plt.close()
#
#
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
#fig.savefig("Movilityvs_N.pdf")
#plt.close()






#"""
#Mobility/Rg vs Delta Cs
#"""
#
#fig,ax=plt.subplots()
#ax.scatter(ave_data[:,4],ave_data[:,-1])
#
#
#x=np.array(ave_data[:,4]).astype(np.float)
#y=np.array(ave_data[:,-1]).astype(np.float)
#
#for i in xrange(len(interactions)):
#    txt="%.2lf,%.2lf"%(interactions[i].epsilon,interactions[i].sigma)
#    ax.annotate(txt, (x[i]-0.006,y[i]+0.002))
#ax.set_xlabel(r'$\Delta c_s [1/\sigma^3]$',fontsize=16)
##ax.grid()
#ax.set_ylabel(r'$ b/R_g [\tau/m\sigma]$',fontsize=16)
#
#ax.tick_params(labelsize=14)
#
##Using np.unique(x) instead of x handles the case where x isn't sorted or has duplicate values.
#ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
#
#ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
#ax.axvline(x=0, ymin=0, ymax=1,ls=':',c='black')
#plt.tight_layout()
#
#fig.savefig("plots/all/Mobility_rg_Delta_Cs.pdf")
#plt.close()
