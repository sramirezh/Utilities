#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 15:15:12 2018
This scripts reads the pos.dat that contains the cm positions of a polymer and computes the diffusion coefficient
@author: sr802
"""

from __future__ import division
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
from scipy import optimize

try:
    from uncertainties import unumpy,ufloat
except ImportError as err2:
    print err2

import Others.Statistics.FastAverager as stat
cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script

def compute_one_msd(pos,delta):
    """
    Computes the MSD for the positions every certain delta
    Args:
        pos: all the positions of the cm of the polymer
        delta: evert this number, we take the positions to compute the msd
    """
    global msd_comp,msd_comp_u,msd_u
    pos=pos[::delta]
    delta_sqr_components=(pos-np.roll(pos,-1,axis=0))**2
    delta_sqr_components=delta_sqr_components[:-1]#The last contribution is the last-the initial msd
    msd_comp=np.array(stat.fast_averager(delta_sqr_components))
    msd_comp=unumpy.uarray(msd_comp[:,1],msd_comp[:,2]) #average and autocorrelation error
    msd=msd_comp[0]+msd_comp[1]+msd_comp[2]

    return np.hstack([msd_comp,msd])

"""
*******************************************************************************
Main
*******************************************************************************
"""
cf.set_plot_appearance()

delta_t= 0.005

print "Assuming that the time step is %s"%delta_t
#times_l,msd_l = lammps_MSD(delta_t)

input_file = "pos.dat"
Data=cf.read_data_file(input_file)
data1=Data.values
times=data1[:,0]*delta_t

pos=data1[:,1::]

msd=[]
t=[]
max_delta=int(len(times)*0.05) #Maximum delta of time to measure the MSD

D_inst=[] #Array with the instantaneous diffusion coefficient
for i in xrange(max_delta):
    msd_t=compute_one_msd(pos,i+1)[3]
    msd.append(msd_t)
    dt=times[i]
    t.append(dt)
    D_inst.append(msd_t/dt/(2*3))



#Writing arrays of averages and errors
t=np.array(t)
msd_error=unumpy.std_devs(msd)
msd_average=unumpy.nominal_values(msd)


D_inst_error=unumpy.std_devs(D_inst)
D_inst_ave=unumpy.nominal_values(D_inst)





#New fit
fitfunc = lambda p, x: p[0] * x + p[1] #Fitting to a line
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

pinit=[1,-1]

out = optimize.leastsq(errfunc, pinit, args=(t[::10],msd_average[::10],msd_error[::10]), full_output=1)

pfinal = out[0] #fitting coefficients
cov=out[1] #Covariance





D=pfinal[0]/(2*3)
D_err=np.sqrt(cov[0][0])*D
print "The diffusion coefficient is %s +/- %s"%(D,D_err)



plt.close('all')
fig1,(ax1,ax12)=plt.subplots(2,1)
ax1.plot(t,msd_average)
ax1.fill_between(t, msd_average-msd_error, msd_average+msd_error ,alpha=0.3)
ax1.plot(np.unique(t),fitfunc(pfinal,np.unique(t)),linestyle='--',c='black')
ax1.set_ylabel(r'$MSD$')
ax12.plot(t,D_inst_ave)
ax12.fill_between(t, D_inst_ave-D_inst_error, D_inst_ave+D_inst_error ,alpha=0.3)
ax12.axhline(y=D, xmin=0, xmax=1,ls='--',c='black')
ax12.set_xlabel(r'$\Delta t$')
ax12.set_ylabel(r'$D$')
plt.tight_layout()
plt.savefig("Diffusio_coefficient.pdf")
plt.show()
