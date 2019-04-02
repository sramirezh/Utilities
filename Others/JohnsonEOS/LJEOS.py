"""
This algorithm is based on the LJ EOS by Johnson et al. 
The Lennard-Jones equation of state revisited.

Reproduces the plot 5.8 in Daan's Book
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import warnings
warnings.filterwarnings("ignore")

sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf


X=np.loadtxt("LJEOS.param")
gamma=3.0 #Defined by historical reasons

def a(T):
    """
    Return a vector of the parameters a for a given 
    temperature T
    """
    A=np.zeros(8)
    A[0] = X[0]*T + X[1]*T**0.5 + X[2] + X[3]/T + X[4]/T**2
    A[1] = X[5]*T + X[6]+ X[7]/T + X[8]/T**2
    A[2] = X[9]*T + X[10] + X[11]/T
    A[3] = X[12]
    A[4] = X[13]/T + X[14]/T**2
    A[5] = X[15]/T
    A[6] = X[16]/T + X[17]/T**2
    A[7] = X[18]/T**2
    return A

def b(T):
    """
    Return a vector of the parameters b for a given 
    temperature T
    """
    B=np.zeros(6)
    B[0] = X[19]/T**2 + X[20]/T**3
    B[1] = X[21]/T**2 + X[22]/T**4
    B[2] = X[23]/T**2 + X[24]/T**3
    B[3] = X[25]/T**2 + X[26]/T**4
    B[4] = X[27]/T**2 + X[28]/T**3
    B[5] = X[29]/T**2 + X[30]/T**3 + X[31]/T**4
    return B
     

def funcF(rho):
    """
    Returns the parameter F
    """    
    F=np.exp(-gamma*rho**2)
    return F

def funcG(rho):
    """
    Return a vector of the parameters G for a given 
    density rho
    """
    F=funcF(rho)
    G=np.zeros(6)
    G[0] = (1-F)/(2.*gamma)
    G[1] = -(F*rho**2-2*G[0]) / (2.*gamma)
    G[2] = -(F*rho**4-4*G[1]) / (2.*gamma)
    G[3] = -(F*rho**6-6*G[2]) / (2.*gamma)
    G[4] = -(F*rho**8-8*G[3]) / (2.*gamma)
    G[5] = -(F*rho**10-10*G[4]) / (2.*gamma)
    return G




def pressure(rho,T):
    """
    Pressure as given by equation 7
    """
    sum1=0
    sum2=0
    A=a(T)
    B=b(T)
    F=funcF(rho)
    for i in xrange(8):
        sum1=sum1+A[i]*rho**((i+1)+1)
        if i<6:
            sum2=sum2+B[i]*rho**(2*(i+1)+1)
    
    P=rho*T+sum1+F*sum2 #Pressure 
    return P



def free_energy(rho,T):
    """
    Residual Helmholtz free energy (excess) of the fluid equation 5
    
    """
    sum1=0
    sum2=0
    A=a(T)
    B=b(T)
    G=funcG(rho)
    for i in xrange(8):
        sum1=sum1+A[i]*rho**((i+1))/(i+1)
        if i<6:
            sum2=sum2+B[i]*G[i]
    A=sum1+sum2 #Pressure 
    return A

def mu_ex(rho,T):
    """
    Residual chemical potential (excess) given by equation 7
    """
    mu=free_energy(rho,T)+pressure(rho,T)/rho-T
    
    return mu
    
def mu_id(rho,T):
    return T*np.log(rho)


def main():
    """
    Beginning of calculations 
    """



cf.set_plot_appearance()    
T=2.0 #Temperature
rho=0.8 #Density 

P=pressure(rho,T)

rho_vect=np.linspace(0.01,1)

p_vect=[]
mu_ex_vect=[]
for i in xrange(np.size(rho_vect)):
    p_vect.append(pressure(rho_vect[i],T))
    mu_ex_vect.append(mu_ex(rho_vect[i],T))


names=[r'$\mu^{ex} $',r'$P $']
fig,ax=plt.subplots()
ax.plot(rho_vect,mu_ex_vect)
ax.plot(rho_vect,p_vect)
ax.set_xlabel(r'$\rho $')
ax.set_xlim(0,1)
ax.set_ylim(-5,10)
plt.legend(names,loc='upper_left')
fig.tight_layout()

if   __name__=="__main__":
    main()



     
     