"""
This script analyzes a chunk of the trajectory file

the trajectory.xyz has to be splitted before with 1Trajectory_Splitter.sh 

It prints the force factor for pressure driven simulations and also other geometrical parameters. 

"""
import numpy as np

#Getting the shape of the data array 
Data=np.loadtxt("1Trajectory.xyz",skiprows=2)
n,m=np.shape(Data) 


#Getting the maximum position of the surface. 
Maxz=-100
Minz=10000
Maxf=0
Nfluid=0


for i in xrange(n):
    if Data[i,0]==3: #3 is for solid surface, 2 for solutes, 1 for solvents.
        Maxz=max(Maxz,Data[i,3])
        Minz=min(Minz,Data[i,3])
    else:
        Nfluid+=1
        
print "The maximum height of the solid surface is %lf" %Maxz
print "The minimum height of the solid surface is %lf" %Minz
print "The height of the solid surface is %lf" %(Maxz-Minz)

#Writing the Zshift
f=open("Zshift.dat",'w')
f.writelines("%lf \n" %Maxz)
f.close        

#Computing the density

#Getting the volume
f=open("Parameters.dat")
FirstLine=f.readline().split()
SecondLine=f.readline().split()
f.close()

Volume=np.float(SecondLine[FirstLine.index("Volume")])
Ff=Volume/Nfluid

print "The vol is %lf" %Volume
print "The Number of fluid particles is %lf" %Nfluid
print "The Ff is %lf" %Ff
