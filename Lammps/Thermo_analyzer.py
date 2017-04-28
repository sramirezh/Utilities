"""
This analyzer loads the dump file containing the instantaneous values printed by the thermostyle in LAMMPS
It is required that the log.lammps file was run before to generate the parameters.dat file
"""
import numpy as np
import matplotlib.pyplot as plt



Nsampling=1000 #This is the number of interations, it took the system to stabilize (It is guessed)
Param=np.loadtxt("parameters.dat",skiprows=1)
#plt.plot(Param[:,6])

#plt.figure(2)
#plt.plot(Param[:,5])

#plt.figure(2)
#plt.plot(Param[:,1])

P_average=np.average(Param[Nsampling::,3])
H_average=np.average(Param[Nsampling::,2]) 
print P_average,H_average
