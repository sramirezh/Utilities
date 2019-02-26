#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 15:07:44 2018

@author: sr802
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate

Rc=2.5

def LJ(r, epsilon, sigma, x, y):
    """
    computes the lennard-jones potential for different exponents.
    see the development in my notebook.
    Args:
        r       : the position at which the potential is computed.
        x       : the exponent of the repulsive part of the potential.
        y       : the exponent of the atractive part of the potential.
        epsilon : Energy of the minimum.
        sigma   : Radiuos at which the potential goes to zero.
    Returns:
        V Lennard Jones interaction energy
    """
    A=((x/y)**(x/(x-y))/((x/y)-1))


    V=A*epsilon*((sigma/r)**x-(sigma/r)**y) #-4*Epsilon*((Sigma/Rc)**12-(Sigma/Rc)**6)

    return V

def frenkel(r, epsilon, sigma, Rc, n):
    """
    Computes the frenkel potential, that naturally goes to zero at Rc
    """

    alpha=2*n*Rc**2*((1+2*n)/(2*n*(Rc**2-1)))**(2*n+1)
    phi = alpha*epsilon*((sigma/r)**2-1)*((Rc/r)**2-1)**(2*n)

    return phi


def force_frenkel(r, epsilon, sigma, Rc, n):
    """
    computes the force given by the frenkel potential
    """

    alpha=2*n*Rc**2*((1+2*n)/(2*n*(Rc**2-1)))**(2*n+1)

    Rc_term=(Rc/r)**2-1
    sigma_term=(sigma/r)**2-1
    first_term=-alpha*epsilon
    second_term=(-2*sigma**2/(r**3))*Rc_term**(2*n)
    third_term=-(Rc**2*4*n*Rc_term**(2*n-1)*sigma_term/r**3)

#    print "The first term is %f, the second term is %f and the third is %f"%(first_term, second_term, third_term)

    force=first_term*(second_term+third_term)

    return force






if __name__=="__main__":
    xmin=1
    r_vec=np.linspace(xmin,Rc)


    r_spline=np.linspace(0.8,Rc+0.1)


    #Frenkel force numerical fitting with splines

    tck = interpolate.splrep(r_spline, frenkel(r_spline,2.5,1,1.6,4), s=0)
    phi_interpol=interpolate.splev(r_vec, tck, der=0)
    force_interpol=-1*interpolate.splev(r_vec, tck, der=1)

    LJ_12_6=LJ(r_vec,1,1,12,6)-LJ(Rc,1,1,12,6)
    LJ_12_6_l=LJ(r_vec,10,1,12,6)-LJ(Rc,10,1,12,6)
    LJ_12_10=LJ(r_vec,1,1,12,10)-LJ(Rc,1,1,12,10)
    PhiF_n_4=frenkel(r_vec,2.5,1,1.6,4)
    force_n_4=force_frenkel(r_vec,1,1,Rc,4)


    plt.close('all')
    fig,ax=plt.subplots()
    ax.plot(r_vec, LJ_12_6, label='$U^{LJ}_{12-6}$')
    ax.plot(r_vec, PhiF_n_4, label='$Frenkel \, n=4$')
    ax.plot(r_vec, phi_interpol,'*',label='Interpolates',)
    ax.axhline(y=0, xmin = xmin, xmax=Rc, ls=':',c='black')
    plt.legend()


    fig2,ax2=plt.subplots()
    ax2.plot(r_vec, force_n_4, label='$ForceFrenkel \, n=4$')
    ax2.plot(r_vec, force_interpol,'*', label='$Force interpol \, n=4$')

    plt.legend()

    #Generating the table for lammps
