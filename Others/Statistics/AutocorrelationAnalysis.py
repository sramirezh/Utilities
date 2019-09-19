import numpy as np
from .Functions import statistics, Autocorrelation 
import matplotlib.pyplot as plt



data=np.loadtxt("data")
n,m=np.shape(data)
Results=statistics(data)
Correlation,time=Autocorrelation(data, Results[0:m])
Error=np.sqrt(Correlation[0,:]*2*time/(n+1))
print(Error)
for i in range(m):
    plt.plot(Correlation[:np.max(np.nonzero(Correlation[:,i])),i],'o-',label='$U_{%s}$' %i )
    plt.plot()
plt.legend()
plt.xlabel("n",fontsize=14)
plt.ylabel("C(n)",fontsize=14)
axes = plt.gca()
axes.set_ylim([-1,8])
plt.show()
    
    
        