"""
This script is going to reproduce yawei's figure 3.c with the distribution force,
It needs to be run inside a folder with raman style (A_dir,B_dir, Stress_dir)
"""
import numpy as np


"""Input parameters"""
#This can be obtained reading the files in the folder
rhos=0.02
rhof=0.74
Temp=0.846091996
Grads=0.0025


"""Reading the parameter files"""
s=np.loadtxt("B_dir/full-mean")
f=np.loadtxt("A_dir/full-mean")
All=np.loadtxt("stress_dir/full-mean")
n,m=np.shape(f)

#Computing the forces
Fs=-Temp/rhos*Grads
Ff=-Fs*rhos/rhof

Force=np.zeros((n,2))
Force[:,0]=f[:,1]

for i in xrange(n):
    if All[i,5]==0:
        Force[i,1]=0
    else:
        Force[i,1]=(Fs*s[i,4]+Ff*f[i,4])/All[i,5]
    
np.savetxt("Force_dist.dat",Force)

