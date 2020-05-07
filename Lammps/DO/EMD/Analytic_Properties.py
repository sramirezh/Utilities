"""
This OLD (for FYR) script analyzes the solute concentraion file, either Raman style or my style.

It computes the diffusio-osmotic mobility and the velocity in the bulk. 


"""

import numpy as np
from scipy.interpolate import splev,splrep,splint
import matplotlib.pyplot as plt

# =============================================================================
# Function Definition
# =============================================================================


def Integrate(x,y,xmin,xmax):
    """
    Integrate the data in x and y from xmin to xmax
    """
    MinIndex=np.min(np.where(x>=xmin))
    MaxIndex=np.max(np.where(x<=xmax))
    I=np.trapz(y[MinIndex:MaxIndex],x[MinIndex:MaxIndex])
    
    return I


# =============================================================================
# Reading data
# =============================================================================

print (" \n Remember to choose which style are you going to use, and comment the other\n")
Concentration=np.loadtxt("Concentration.dat")
Concentration=np.loadtxt("Averages.dat", usecols=[1,4]) #To use with results with Raman Style
n,m=np.shape(Concentration)

# =============================================================================
# Computing the bulk Concentration
# =============================================================================

BulkMin=15
BulkMax=25
BulkC=0
cont=0
for i in range(n):
    if (Concentration[i,0]>=BulkMin and Concentration[i,0]<BulkMax):
        cont=cont+1
        BulkC=BulkC+Concentration[i,1]
BulkC=BulkC/cont



# =============================================================================
# Gamma Calculations
# =============================================================================


#Integration limits as defined by Bocquet

xmin=0.0
xmax=8.0


Integrand=Concentration[:,1]/BulkC-1.0
                      
SolAbso=Integrate(Concentration[:,0],Integrand,xmin,xmax)


"""
Spline Calculations 
"""
tck=splrep(Concentration[:,0],Integrand)
tck0=splrep(Concentration[:,0],Concentration[:,1])
SplineX=np.linspace(xmin,xmax,400)
SplineY=splev(SplineX,tck)

Gamma=splint(xmin,xmax,tck) #Solute Absorption


# =============================================================================
# Diffusio-Osmotic Movility K_DO
# =============================================================================

etha=2.60 #Viscosity

IntegrandKdo=(Concentration[:,1]/BulkC-1.0)*Concentration[:,0]

Kdo=Integrate(Concentration[:,0],IntegrandKdo,xmin,xmax)/etha



# =============================================================================
# Bulk Velocity
# =============================================================================
IntegrandVb=(BulkC-Concentration[:,1])*Concentration[:,0]
Vb=Integrate(Concentration[:,0],IntegrandVb,xmin,xmax)*0.063

print(Vb)


#
##Generating Spline Concentration
#TSplineX=np.linspace(0,25,1500)
#TSplineY=splev(TSplineX,tck0)
#
#print "The average concentration is %f" %(Average_c)
#print "The bulk concentration is %f" %(BulkC)
#print "The solute adsorption is %f" %SolAbso
#print "The solute adsorption using splines is %f" %Gamma
#                       
#"""
#Creating the output file
#"""
#np.savetxt("Concentration.dat",Concentration)
#np.savetxt("Spline_Concentration.dat", np.transpose([TSplineX,TSplineY]))
#np.savetxt("Integrand.dat",np.transpose([Concentration[:,0],Integrand]))
#np.savetxt("Spline_Integrand.dat",np.transpose([SplineX,SplineY]))
#
#print "Generated Concentration.dat Spline_Concentration.dat Integrand.dat Spline_Integrand.dat"
#"""
#Uncomment for testing
#"""
# 
##import matplotlib.pyplot as plt
##plt.figure(1)
##plt.plot(Concentration[:,0],Concentration[:,1],'*')
##
##x=np.linspace(min(Concentration[:,0]),max(Concentration[:,0]))
##y=np.zeros(len(x))
##y[:]=BulkC
##plt.plot(x,y)
##plt.plot(TSplineX,TSplineY)
##plt.xlim([0,25])
##
##
##
##plt.figure(2)
##
##plt.plot(Concentration[:,0],Integrand,'*')
##x=np.linspace(min(Concentration[:,0]),max(Concentration[:,0]))
##y=np.zeros(len(x))
##plt.plot(x,y)
##plt.plot(SplineX,SplineY)
##plt.xlim([0,8])
#
#
