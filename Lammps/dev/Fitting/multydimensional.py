#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:10:12 2019

@author: simon
"""
import numpy as np
from scipy.optimize import curve_fit

#Lets create the known polynomial

def p_known(x,y):
    z=3*x*y**2+4*y*x**2+5*x**2
    return z

powerlaw = lambda x, amp, index: amp * (x**index)



def arbitrary_poly(point, *params):
    """
    Creates a polymer with sqrt(Number of params), where params are the coefficients
    Args:
        params: coefficients
        point: contains the the two independent variables x and y
        
        p00+p01*y+...+p11*x*y
    
    """
    x=point[0]
    y=point[1]
    dim=int(np.sqrt(np.size(params)))
    params=np.reshape(params,(dim,dim))
    poly=0
    for i in xrange(dim):
        for j in xrange(dim):
            poly+=params[i,j]*x**i*y**j
            
    return poly


n_points=1000

x=np.linspace(0,10,n_points)
y=np.linspace(0,10,n_points)
variables=np.stack((x,y),axis=0)
z=p_known(x,y)

zerr= np.random.rand(n_points)                     # simulated errors (10%)
d=2
popt, pcov = curve_fit(arbitrary_poly, variables[:,:], z, sigma=zerr, p0=[1]*(d+1)**2)

popt_matrix=np.reshape(popt,((d+1),(d+1)))


#Analizing the estimation

z_predict=[]
for i in xrange(n_points):
    z_predict.append(arbitrary_poly(variables[:,i],popt))

diff=z-z_predict