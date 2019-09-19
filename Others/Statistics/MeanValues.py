import numpy as np
from .Functions import statistics, Autocorrelation 
import matplotlib.pyplot as plt



data=np.loadtxt("data")
n,m=np.shape(data)
Results=statistics(data)
Correlation,time=Autocorrelation(data, Results[0:m])
Error=np.sqrt(Correlation[0,:]*2*time/(n+1))
Variance=Error**2*n #variance of individual data U1, U2, U3...
Average=Results[0:m] #Averages of measurements

"""Independent error Analysis"""
#Working with Rij=Ui/Uj
IndError=np.zeros((m,m))
for i in range(m):
    for j in range(m):
        if i!=j:
            IndError[j,i]=np.sqrt((1.0/n)*((Average[i]/Average[j])**2)*(Variance[j]/Average[j]**2+Variance[i]/Average[i]**2))

"""Correct error Analysis"""
NewVariable=np.zeros((n,m))
CorError=np.zeros((m,m))
###Indexes are like in the pdf, Ui/Uj
for j in range(m):
    for i in range(m):
        if i!=j:
            NewVariable[:,i]=data[:,i]/Average[i]-data[:,j]/Average[j]
    Results=statistics(NewVariable)
    
    #Parameter calculation for the new variable using Autocorrelation Analysis.
    Correlation,time=Autocorrelation(NewVariable, Results[0:m])
    Error=np.sqrt(Correlation[0,:]*2*time/(n+1))
    Variance=Error**2*n 
    NAverage=Results[0:m]
    cont=0
    for i in range(m):
        if i!=j:
            CorError[j,i]=np.sqrt((1.0/n)*((Average[i]/Average[j])**2)*Variance[i])

###Print results ready to be put in a table in Latex
np.savetxt("IndError",IndError,fmt='%6.5f',delimiter=" & ")            
np.savetxt("CorError",CorError,fmt='%6.5f',delimiter=" & ")


RelErr=np.zeros((5,5))
for i in range(5):
    for j in range(5):
        RelErr[i,j]=np.abs(IndError[i,j]-CorError[i,j])/(0.5*(IndError[i,j]+CorError[i,j]))*100.0
