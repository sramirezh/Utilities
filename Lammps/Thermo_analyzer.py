"""
This analyzer loads the dump file containing the instantaneous values printed by the thermostyle in LAMMPS
It is required that the log.lammps file was run before to generate the parameters.dat file
"""
import numpy as np
import matplotlib.pyplot as plt


#Next two lines are for testing 
#PATH="/home/sr802/Dropbox/PhD/Cambridge/Academical/3.Simulation/0.Lammps/2.Solid_fluid/2D"
#TotalPath=PATH+"/parameters.dat"

TotalPath="parameters.dat"

#Display the header file
f=open(TotalPath, 'r')
FirstLine=f.readline().split()
f.close()

Param=np.loadtxt(TotalPath,skiprows=1)
n,m=np.shape(Param)
Nsampling=int(n*0.5) #This is the number of interations, it took the system to stabilize (It is guessed)
#plt.plot(Param[:,6])

#plt.figure(2)
#plt.plot(Param[:,5])

#plt.figure(2)
#plt.plot(Param[:,1])

for i in xrange(m):
    print "The Average %s is %lf " %(FirstLine[i],np.average(Param[Nsampling::,i]))

#plt.plot(Param[:,0],Param[:,5])
