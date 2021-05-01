#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 15:37:04 2018
This scripts analyses the pair correlation function from the *.gz files loading
them in pylammps
@author: simon
"""



import numpy as np
import pandas as pd
import argparse
import os
import sys
import pickle as pickle
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

import time

from joblib import Parallel, delayed
import multiprocessing


try:
    from lammps import IPyLammps
except ImportError as err:
    print(err)

cwd = os.getcwd() #current working directory
dir_path = os.path.dirname(os.path.realpath(__file__))#Path of this python script


"""
*******************************************************************************
Functions
*******************************************************************************
"""

def call_pylammps(fil,p_types):
    """
    Calls pylammps to read the configuration and compute the g(r)
    fil is the name of the dump file
    p_types is an array with the particle types.
    """

    file_time=int(list(filter(str.isdigit,fil)))

    L=IPyLammps()
    L.region("box block", 0, 20, 0, 1.9469999999999999e+01, 0, 30)
    L.create_box(len(p_types), "box")

    L.command("read_dump %s %s x y z ix iy iz vx vy vz box yes add yes"%(fil,file_time))
    L.command("mass * 1.0")
    L.pair_style("lj/cut", 2.5)
    L.command("pair_coeff * * 1.0 1.0 2.5")


    for i in p_types:
        for j in p_types:
            if j>=i:
                compute_name="gr_%s%s"%(i,j)
                fix_name="f_%s" %compute_name
                L.command("compute %s all rdf %s %s %s"%(compute_name,nbins,i,j))
                L.fix("%s all ave/time 1 1 1 c_%s[*] file %s_%s.dat mode vector" %(fix_name,compute_name,compute_name,file_time))

    L.run(0)
    L.close()


def read_gr(fil, p_types):
    """
    Gathers all the g(r) from the different files for a given timestep
    """

    results=[]
    file_time=int(list(filter(str.isdigit,fil)))
    for i in p_types:
        for j in p_types:
            if j>=i:
                compute_name="gr_%s%s"%(i,j)
                file_name="%s_%s.dat" %(compute_name,file_time)
                Data=pd.read_csv(file_name,sep=" ",skiprows=4,dtype=np.float64,header=None).values
                os.remove(file_name)
                results.append(Data)

    return results



def run_one_time(fil,p_types):
    """
    Runs everything for a sigle time step
    args:
        fil: is the filename
        p_types is an array containing the species

    """
    call_pylammps(fil,p_types)
    results=read_gr(fil, p_types)

    return results


def run_analysis(input_files,p_types):

    input_files.sort(key=lambda f: int(list(filter(str.isdigit, f))))
    t=time.time()


    num_cores = multiprocessing.cpu_count()
    results=Parallel(n_jobs=num_cores,verbose=10)(delayed(run_one_time)(fil,p_types) for fil in input_files)


    g_r=np.average(results,axis=0)
    print("The total time was %s seconds, the g(r) was saved as gr.pkl "%(time.time()-t))

    afile = open(r'gr.pkl', 'wb')
    pickle.dump(g_r, afile)
    afile.close()

    return g_r

def load_gr(file_name):
    """
    Loads the data structure to be used later
    """
    file1 = open(file_name[0], 'rb')
    g_r = pickle.load(file1)
    file1.close()

    return g_r


def plot_gr(g_r):
    """
    ###############################################################################
    Starting the plot
    ###############################################################################
    """
    colors=['r','b','k','g','c']
    axis_font=24
    tick_font=20
    legend_font=18


    fig,ax=plt.subplots()
    names=[r"$g_{ff}$",r"$g_{fs}$",r"$g_{mf}$",r"$g_{ss}$",r"$g_{ms}$"]
    markers=['.','^','v','s','D']
    markevery=[5,6,5,6,5]

    for i in range(5):
        if i<len(colors):color=colors[i]
        else: color=np.random.rand(3)

        ax.plot(g_r[i][:,1],g_r[i][:,2],label=names[i],color=color,marker=markers[i],markersize=3,markevery=markevery[i])

    rmax=np.max(g_r[:][:,1])

    """Axis"""
    ax.set_xlabel(r'$r[\sigma] $',fontsize=axis_font)
    ax.grid(False)
    ax.set_ylabel(r'$g(r)$',fontsize=axis_font)
    ax.tick_params(labelsize=tick_font,direction='in',top=True, right=True)



    #ax.axhline(y=1, xmin=0, xmax=1,ls='--',c='black')
    ax.axvline(x=2**(1/6), ymin=0, ymax=1,ls=':',c='black')



    xmin,xmax=plt.xlim()

    #plt.xticks(np.arange(len(lengths)+1)*30)
    ax.set_xlim(0,rmax)


    ymin,ymax=plt.ylim()

    #ax.set_ylim(0,ymax)
    #plt.yticks(np.arange(-0.1,0.2,0.1))
    #ax.set_xlim(xmin-deltax*xoffset,xmax+deltax*xoffset)

    """Legend"""
    plt.legend(fontsize=legend_font,loc='upper right',labelspacing=0.1,borderpad=0.4,scatteryoffsets=[0.6],
               frameon=True, fancybox=False, edgecolor='k')


    """General"""
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["text.usetex"] =True
    plt.tight_layout()
    fig.savefig("g_r.pdf",transparent=False)
    plt.close()

"""
*******************************************************************************
Main
*******************************************************************************
"""
parser = argparse.ArgumentParser(description='This script reads trajectory files and computes the pair distribution function',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-mode',help='Indicate if you want to run the analyisis or just read a previous computed g(r)',choices=['read','run'],default="read")
parser.add_argument('-input_files', metavar='InputFile',help='Input filename',nargs='+',type=lambda x: cf.is_valid_file(parser, x))
parser.add_argument('-nbins', help='Number of bins in the radial direction', default=100, type=float)
parser.add_argument('-nmin', help='Number of timesteps to be discarded', default=500, type=int)


args = parser.parse_args()
input_files=args.input_files


logger = cf.log(__file__, os.getcwd())  
logger.info("Using the following arguments for the paser")
logger.info(args)


if args.mode=="run":

    imin=args.nmin
    nbins=args.nbins
    rmax=2.5
    input_files = input_files[imin:]
    p_types=[1,2,3]
    g_r=run_analysis(input_files, p_types)

else:
#    try:
    g_r = load_gr(input_files)
#    except:
#        sys.exit("There is no file to load, run the analysis first!")

plot_gr(g_r)
