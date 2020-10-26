#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:10:12 2019

@author: simon
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import sys
import fitting_functions as ff
sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf

# =============================================================================
# MAIN
# =============================================================================
plot_dir = "plots/simulataneous_fit"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
    
logger = cf.log(__file__, os.getcwd(),plot_dir)    

rho_ref = 0 #0.32
beta_ref = 0  #0 1/1.04
deg_x = [2,8] 
deg_y = [0,6]
file_p='pressure/pressure_liquid_rc2.0.txt' # Temp rho press Sigmapress
file_e='energy/energy_liquid_rc2.0.txt'
p_ref = 0
e_ref = 0

logger.info("rho_ref = %s"%rho_ref)
logger.info("beta_ref = %s"%beta_ref)

# =============================================================================
#     Here is where I define the polynomials
# =============================================================================

# Pressure
poly_p = ff.Polynomial(deg_x,deg_y,[1,-1],[1],[1,0],[1,0])
# Energy
poly_e  = ff.Polynomial(deg_x,deg_y,[1],[1,0],[1,0],[1,-1])

# =============================================================================
#     Reading the data
# =============================================================================

# Pressure
x_p, y_p, z_p, zerr_p = ff.read_data(file_p, p_ref)
x_p = x_p - rho_ref
y_p = y_p - beta_ref

# Energy
x_e, y_e, z_e, zerr_e = ff.read_data(file_e, e_ref)
x_e = x_e - rho_ref
y_e = y_e - beta_ref

# =============================================================================
#     Simultaneous fitting
# =============================================================================

popt, pcov, variables = ff.fit_general([[x_p, y_p, z_p, zerr_p], poly_p],[[x_e,y_e,z_e,zerr_e],poly_e])


plt.close('all')

# =============================================================================
#     Testing Results
# =============================================================================
print ('\nTesting the energy results\n')

# Computing the predictions
variables_e = np.stack((x_e, y_e),axis=0)
z_e_predicted = poly_e.evaluate(variables_e, popt)

# Computing the error
er_results_e = ff.test_prediction(z_e, z_e_predicted, e_ref)

# Writing the output data
ff.outputs(popt,pcov, er_results_e, deg_x, deg_y,'e')

# Plotting the results
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(x_e, y_e, z_e, zdir = '-z',marker ='.', label="Simulation", color ='r')
ax1.scatter(x_e, y_e, z_e_predicted, zdir='-z', marker='.', label="Fitting", color = 'b')
ax1.set_xlabel(r'$\rho$', labelpad = 5)
ax1.set_ylabel(r'$\beta$', labelpad = 5)
ax1.set_zlabel(r'$\rho \ e_{\text {exc}}$', labelpad = 5)
fig1.subplots_adjust(left=0, right=1, bottom=0, top=1) 
ax1.legend(loc=(0.25,0.65),frameon=0)
fig1.show()

fig1.savefig("%s/3Dplot_e.pdf"%plot_dir)
logger.info("plotted %s/3Dplot_e.pdf"%plot_dir)


print ('\nTesting the pressure results\n')

# Computing the predictions
variables_p = np.stack((x_e,y_e),axis=0)
z_p_predicted = poly_p.evaluate(variables_p, popt)

# Computing the error
er_results_p = ff.test_prediction(z_p, z_p_predicted, p_ref)

# Writing the output data
ff.outputs(popt, pcov, er_results_p, deg_x, deg_y, 'p')

# Plotting the results
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(y_p, x_p, z_p, zdir='z',marker='.',label="Simulation", color='r')
ax2.scatter(y_p, x_p, z_p_predicted, zdir='z', marker='.',label="Fitting", color='b')
ax2.set_ylabel(r'$\rho$', labelpad = 5)
ax2.set_xlabel(r'$\beta$', labelpad = 5)
ax2.set_zlabel(r'$ \beta\ P_{\text {exc}}$', labelpad = 5)

# To reduce the empty space around the figure
fig2.subplots_adjust(left=0, right=1, bottom=0, top=1) 
ax2.legend(loc=(0.55,0.65),frameon=0)
fig2.show()
fig2.savefig("%s/3Dplot_p.pdf"%plot_dir)
logger.info("plotted %s/3Dplot_p.pdf"%plot_dir)

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
    

