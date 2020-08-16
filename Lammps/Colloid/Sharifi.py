#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 12:18:38 2019

Simple script to get Sharifi concentrations:
    
It is an ideal solution in the bulk

@author: sr802
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

import Lammps.core_functions as cf



def concentration_distribution(distance,c_gradient, vol_bin, n0 = 50):
    concentration = (n0+distance*grad*vol_bin)/vol_bin
    
    
    return np.array(concentration)


# =============================================================================
# Main
# =============================================================================
plot_dir = "plots/1.Sharifi"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)    


grad_c = np.array([0.0062, 0.0077, 0.0093, 0.0108])

n0 = 50
lx = 53.86
ly = 43.09
lz = 43.09
r_colloid = 3.23
T = 1
nbins = 8

delta_x = lx/nbins

vol_bin = delta_x*ly*lz

delta_N = grad_c * vol_bin*delta_x

cf.set_plot_appearance()

xpos = [i*delta_x for i in range(nbins+1)]


fig,ax = plt.subplots()

for grad in grad_c:
    delta_N = grad * vol_bin*delta_x
    ns = np.array([n0+i*delta_N for i in range(nbins+1)])
    ax.step(xpos,ns,label = r'$\nabla c_s = %s $'%grad, where = 'post')
    
    #Exact distribution
    ns_exact = [n0+x*grad*vol_bin for x in xpos]
    color=ax.lines[-1].get_color()
    
    ax.plot(xpos,ns_exact, color = color,linestyle='dashed' )

ax.set_xlim(0, lx)   
ymin,ymax=plt.ylim()
ax.set_ylim(0, 1.2*ymax)  
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$N_s(x)$")
ax.legend()
plt.tight_layout()
plt.savefig("%s/Grad_c.pdf"%plot_dir)
logger.info("%s/Grad_c.pdf"%plot_dir)

fig,ax = plt.subplots()
fig2,ax2 = plt.subplots()

grad_mu = []
m_cs = []
cs_middle = []
grad_mu_middle = []

dist = np.linspace(0,lx)

for grad in grad_c:
    delta_N = grad * vol_bin*delta_x
    cs = np.array([n0+i*delta_N for i in range(nbins+1)])/vol_bin
    ax.step(xpos,cs,label = r'$\nabla c_s = %s $'%grad, where = 'post')
    m_cs.append(cs)
    #Exact distribution
    cs_exact = np.array([n0+x*grad*vol_bin for x in xpos])/vol_bin
    cs_middle.append((n0+(lx/2)*grad*vol_bin)/vol_bin)
    grad_mu_middle.append(grad/((n0+(lx/2)*grad*vol_bin)/vol_bin))
    
    color=ax.lines[-1].get_color()
    
    ax.plot(xpos,cs_exact, color = color,linestyle='dashed' )
    
    con_dist = concentration_distribution(dist,grad,vol_bin)
    
    
    grad_mu = T* grad/con_dist 
    
    ax2.plot(dist,grad_mu,label = r'$\nabla c_s = %s $'%grad )
    
    
grad_mu = np.array(grad_mu)
m_cs = np.array(m_cs)



# Concentration distribution
ax.set_xlim(0, lx)   
ymin,ymax=plt.ylim()
#ax.set_ylim(0, 1.05*ymax)  
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$c_s(x)$")
ax.axvspan(lx/2-r_colloid,lx/2+r_colloid, alpha=0.5, color='red')
ax.legend()
plt.tight_layout()
plt.savefig("%s/conc.pdf"%plot_dir)
logger.info("plotted %s/conc.pdf"%plot_dir)

#Chemical potential gradient
ax2.set_xlim(22,32)
ax2.set_ylim(0.03,0.042)
ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"$\nabla_{\mu}(x)$")
ax2.axvspan(lx/2-r_colloid,lx/2+r_colloid, alpha=0.5, color='red')
ax2.legend()
plt.tight_layout()
plt.savefig("%s/grad_mu_dist.pdf"%plot_dir)
logger.info("Plotted %s/grad_mu_dist.pdf"%plot_dir)


# Sharifi Results
vx = [0.0028,0.0032, 0.004,  0.0046]




    

