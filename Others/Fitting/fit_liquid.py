#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:10:12 2019

@author: simon
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from cycler import cycler
import os
import shutil
import sys
import fitting_functions as ff
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf

# =============================================================================
# MAIN
# =============================================================================


    
#def main(rho_ref,beta_ref,deg_x,deg_y,p_ref,e_ref,file_p,file_e):
#    '''
#    x is the density
#    y is beta
#    z is the property
#    '''
    


rho_ref = 0 #0.32
beta_ref = 0  #0 1/1.04
deg_x = [2,8] 
deg_y = [0,6]
file_p='pressure/pressure_liquid_rc2.0.txt' # Temp rho press Sigmapress
file_e='energy/energy_liquid_rc2.0.txt'
p_ref = 0
e_ref = 0



# =============================================================================
#     Here is where I define the polynomial for the pressure
# =============================================================================
#    print "Basic information about the polynomial for the Pressure"
poly_p = ff.polynomial(deg_x,deg_y,[1,-1],[1],[1,0],[1,0])



# =============================================================================
#     Here is where I define the polynomial for the pressure
# =============================================================================
#    print "Basic information about the polynomial for the Energy"
poly_e  = ff.polynomial(deg_x,deg_y,[1],[1,0],[1,0],[1,-1])



print("rho_ref = %s"%rho_ref)
print("beta_ref = %s"%beta_ref)

# =============================================================================
#     #Reading the data for the pressure
# =============================================================================
x_p, y_p, z_p, zerr_p = ff.read_data(file_p, p_ref)
x_p = x_p - rho_ref
y_p = y_p - beta_ref



# =============================================================================
#     #Reading the data for the energy
# =============================================================================
x_e, y_e, z_e, zerr_e = ff.read_data(file_e, e_ref)
x_e = x_e - rho_ref
y_e = y_e - beta_ref




#Simultaneous Fitting
popt, pcov, variables = ff.fit_two_poly([[x_p, y_p, z_p, zerr_p], poly_p],[[x_e,y_e,z_e,zerr_e],poly_e])


plt.close('all')


#
## =============================================================================
##     #Testing for Pressure
## =============================================================================
#print ('\nTesting the pressure results\n')
#variables_p=np.stack((x_p,y_p),axis=0)
#er_results_p=test_prediction(popt,variables_p,z_p,poly_p)
#
#outputs(popt,pcov,er_results_p, deg_x, deg_y,'p')
#
#fig1 = plt.figure()
#ax = fig1.add_subplot(111, projection='3d')
#ax.scatter(x_p, y_p, z_p, zdir='z',marker='.',label="Simulation",color='r')
#
#
#
#
#
##Creating the surface
#x,y=np.meshgrid(np.linspace(np.min(x_p),np.max(x_p),20),np.linspace(np.min(y_p),np.max(y_p),20))
#
#
#
#z=np.asarray(x)
#variables_p=np.stack((x.flatten(),y.flatten()),axis=0)
#er_results=test_prediction(popt,variables_p,z.flatten(),poly_p)
#z_mesh=np.reshape(er_results[:,1],np.shape(x))
#ax.plot_wireframe(x,y,z_mesh,color='b')
#
#ax.set_xlabel(r'$\rho-\rho^*$')
#ax.set_ylabel(r'$\beta-\beta^*$')
#ax.set_zlabel(r'$\Delta \beta P_{exc}$')
#fig1.legend()
#fig1.show()
#fig1.savefig("3Dplot_p.pdf")


# =============================================================================
#     #Testing for Energy
# =============================================================================
print ('\nTesting the energy results\n')
variables_e=np.stack((x_e,y_e),axis=0)
er_results_e = ff.test_prediction(popt,variables_e,z_e,poly_e, e_ref)
#    
ff.outputs(popt,pcov,er_results_e, deg_x, deg_y,'e')
#    
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(x_e, y_e, z_e, zdir='z',marker='.',label="Simulation",color='r')
ax2.scatter(x_e, y_e, er_results_e[:,1], zdir='z',label="Fitting",color='black')
#    
#Creating the surface
#x,y=np.meshgrid(np.linspace(np.min(x_e),np.max(x_e),20),np.linspace(np.min(y_e),np.max(y_e),20))




#z=np.asarray(x)
#variables_e=np.stack((x.flatten(),y.flatten()),axis=0)
#er_results=test_prediction(popt,variables_e,z.flatten(),poly_e)
#z_mesh=np.reshape(er_results[:,1],np.shape(x))
#ax2.plot_wireframe(x,y,z_mesh,color='b')
ax2.set_xlabel(r'$\rho-\rho^*$')
ax2.set_ylabel(r'$\beta-\beta^*$')
ax2.set_zlabel(r'$\Delta \rho \ e_{exc}$')

fig2.legend()
fig2.show()
fig2.savefig("3Dplot_e.pdf")





## =============================================================================
## #    Single fitting
## =============================================================================
#
#set_plot_appearance()
#popt1,pcov1,variables1=fit_poly(x_p,y_p,z_p,zerr_p,poly_p)
#single_results_p=test_prediction(popt1,variables1,z_p,poly_p)
#
#
#popt2,pcov2,variables2=fit_poly(x_e,y_e,z_e,zerr_e,poly_e)
#single_results_e=test_prediction(popt2,variables2,z_e,poly_e)
#
#
#def plot_slices():
#    # Plot slices in beta
#    
#    dir_name="slices_beta_p"
#    shutil.rmtree(dir_name,ignore_errors=True) 
#    
#    os.mkdir(dir_name)
#    slices=np.unique(y_p)
#    for sli in slices:
#        
#        fig3=plt.figure()
#        ax3=fig3.add_subplot(111)
#        error_cap=4
#        indexes=np.where(y_p==sli)
#        ax3.errorbar(x_p[indexes],z_p[indexes],yerr=zerr_p[indexes],fmt='o',capsize=error_cap,label='beta=%s'%sli)
#        ax3.plot(x_p[indexes],er_results_p[indexes[0],1],label='Fit 2')
#        ax3.plot(x_p[indexes],single_results_p[indexes[0],1],label='Fit 1')
#        fig3.legend(loc='upper right')
#        fig3.tight_layout()
#        fig3.savefig("%s/slice_%s.pdf"%(dir_name,sli))
#        
#        
#    dir_name="slices_beta_e"
#    shutil.rmtree(dir_name,ignore_errors=True) 
#    os.mkdir(dir_name)
#    
#    slices=np.unique(y_e)
#    for sli in slices:
#        
#        fig3=plt.figure()
#        ax3=fig3.add_subplot(111)
#        error_cap=4
#        indexes=np.where(y_e==sli)
#        ax3.errorbar(x_e[indexes],z_e[indexes],yerr=zerr_e[indexes],fmt='o',capsize=error_cap,label='beta=%s'%sli)
#    
#        ax3.plot(x_e[indexes],er_results_e[indexes[0],1],label='Fit 2')
#        ax3.plot(x_e[indexes],single_results_e[indexes[0],1],label='Fit 1')
#        fig3.legend(loc='upper right')
#        fig3.tight_layout()
#        fig3.savefig("%s/slice_%s.pdf"%(dir_name,sli))


#if __name__ == "__main__":
#    parser = argparse.ArgumentParser(description='This script fits a 2 variable data to a poly of deg n',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#    parser.add_argument('-file_p', metavar='file_p',help='Input filename with pressure results',default='LP')
#    parser.add_argument('-file_e', metavar='file_e',help='Input filename with energy results',default='LE')
#    parser.add_argument('-beta_ref', metavar='beta reference',help='reference beta',default=1/1.1,type=float)
#    parser.add_argument('-rho_ref',metavar='rho reference',help='reference rho',default=0,type=float )
#    parser.add_argument('-P_ref', metavar='Reference Pressure',help='reference pressure',default=0,type=float)
#    parser.add_argument('-E_ref', metavar='Reference Energy',help='reference energy',default=0,type=float)
#    parser.add_argument('-deg_x',metavar='poly degree in rho',help='Degree of the poly, or min->max degree',nargs='+', default=[0,5],type=int)
#    parser.add_argument('-deg_y',metavar='poly degree in beta',help='Degree of the poly',nargs='+',default=[5],type=int)
#    parser.add_argument('-exc_x',metavar='Exponents to exclude in x',help='Exponents to exclude in x as a list i.e 0 3 4',nargs='+',type=int)
#    parser.add_argument('-exc_y',metavar='Exponents to exclude in y',help='Exponents to exclude in y as a list i.e 0 3 4',nargs='+',type=int)
#    
#    
#    args = parser.parse_args()
#    
#    
#    main(args.rho_ref,args.beta_ref,args.deg_x,args.deg_y,args.P_ref,args.E_ref,args.file_p,args.file_e)
#    
    
