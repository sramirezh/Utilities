'''
This scripts determines the shear viscosity from pressure driven simulations
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splev,splrep,splint


# =============================================================================
# input parameters 
# =============================================================================
Data=np.loadtxt("AAverages.dat")
Coeff=np.polyfit(Data[:,1],Data[:,3],2)
PGrad=-0.0005


x=np.linspace(np.min(Data[:,1]),np.max(Data[:,1]),500)
y=np.polyval(Coeff,x)

plt.plot(x,y)
plt.plot(Data[:,1],Data[:,3],'.')

Viscosity=PGrad/(2*Coeff[0])

print Viscosity