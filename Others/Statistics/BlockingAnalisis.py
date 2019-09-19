import numpy as np
import matplotlib.pyplot as plt
from .Functions import statistics, blocking
    
data=np.loadtxt("data")
n,m=np.shape(data)
MaxIter=np.round(np.log(n)/np.log(2))
Results=np.zeros((MaxIter,3*m))
Error=np.zeros((MaxIter,m))


for j in range(m):
    i=0
    diff=1
    data1=data
    while diff>0 and i<MaxIter:    ###Always ascending results that's why it is required diff>0
        Results[i,:]=statistics(data1)
        data1=blocking(data1)
        Error[i,j]=np.sqrt((2.0**i/n)*Results[i,2*m+j])
        if i>0:
            diff=Error[i,j]-Error[i-1,j]
        i+=1

for i in range(m):
    plt.plot(Error[:np.max(np.nonzero(Error[:,i])),i],'o-',label='$U_{%s}$' %i )
    plt.plot()
plt.legend(loc=2)
plt.xlabel("Blocking step",fontsize=14)
plt.ylabel("$\sigma$",fontsize=20)
plt.show()
Nonblock=statistics(data)
print("Blocking analysis")
for i in range(m):
    print((Error[np.max(np.nonzero(Error[:,i]))-1,i]))
print(("Non Blocking Results",Nonblock[m:2*m]))