#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:10:12 2019

@author: simon
"""
import numpy as np
import argparse
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from cycler import cycler
import os
import shutil

def set_plot_appearance():

    """
    Defines the appearence of the plot

    use rcParams.keys() to see the available parameters
    """
    axis_font=24
    tick_font=20
    legend_font=18


    #plt.rcParams['lines.linewidth'] = 1.5
    #plt.rcParams['lines.markeredgewidth'] = 1.0
    #plt.rcParams['lines.markersize'] = 2.5

    # Colors
    plt.rcParams['axes.prop_cycle'] = cycler(color='rbkgcmy')


    # Fonts and symbols
    #plt.rcParams['font.family'] = 'serif'
    #plt.rcParams['font.serif'] = 'Times New Roman'
    #plt.rcParams['font.weight'] = 'normal'
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams['text.usetex'] = True
    #plt.rcParams['mathtext.rm'] = 'serif'
    #plt.rcParams['mathtext.it'] = 'serif:italic'
    #plt.rcParams['mathtext.fontset'] = 'stix'


    # Axes
    plt.rcParams['axes.edgecolor'] = (0.0, 0.0, 0.0)
    #plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['axes.spines.right'] = True
    plt.rcParams['axes.spines.top'] = True
    plt.rcParams['axes.labelsize'] = axis_font
    plt.rcParams['axes.grid'] = False

    # Ticks
    plt.rcParams['xtick.direction'] = "in"
    plt.rcParams['ytick.direction'] = "in"
    plt.rcParams['ytick.labelsize'] = tick_font
    plt.rcParams['xtick.labelsize'] = tick_font

    # Legend

    plt.rcParams['legend.fontsize'] = legend_font
    plt.rcParams['legend.loc'] ='upper left'
    plt.rcParams['legend.labelspacing'] = 0.5
    plt.rcParams['legend.borderpad'] =0.4
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.fancybox'] = False
    plt.rcParams['legend.edgecolor'] = 'k'
    # Fonts and symbols


    # Errorbar plots
    plt.rcParams['errorbar.capsize'] = 4
    
# =============================================================================
# Class polynomial
# =============================================================================


class polynomial(object):
    """
    Class to identify a polynomial 
    
    p(x,y)=\sum_{i,j}^{n,m} [ fn(i) fm(j) c_{i,j} x^{f_expx(i)} y^{f_expy(i)}
    
    with f(n),f(m) additional funcitons that are independen of the fitting coefficient
    """
    
    def __init__(self, n,m,fn,fm,f_expx,f_expy,exc_n=[],exc_m=[]):
        """
        
        Args:
            n is the exponent in x
            m is the exponent in y
            exc_n are the excluded exponents in x
            exc_m are the excluded exponents in y
            fn array describing a poly to describe f(n) that is independent of c_{nm}
            fm array describing a poly to describe f(m) that is independen independent of c_{nm}
            
            f_expx array describing the function that generates the exponent of x
            f_expy array describing the function that generates the exponent of y
            
            fn,fm,f_an,f_am are going to be evaluated with np.polyval, so see the syntaxis
        
        """
        self.exponents=[get_list(n),get_list(m)]
        self.excluded=[get_list(exc_n),get_list(exc_m)]
        self.get_limits()
        self.get_exponents()
        self.get_dimensions()
        self.func_coeff=[fn,fm]
        self.func_exp=[f_expx,f_expy]
        
        

    
        
    def get_limits(self):
        limits=[]
        for exp in self.exponents:
            if len(exp)==1:
                exp=[0,exp[0]]
            limits.append(exp)
        self.limits=limits
        
        
    
    def get_exponents(self):
        """
        creates the list of exponents in the i-th dimension, exluding the given exponents
        
        """
        for i,lim in enumerate(self.limits):
            exponents=np.arange(lim[0],lim[-1]+1)
            if len(self.excluded[i])!=0:
                exponents=delete_values(exponents,self.excluded[i])
            self.exponents[i]=exponents
        return self.exponents
    
    def get_dimensions(self):
        self.get_exponents()
        dimensions=[]
        for exp in self.exponents:
            dimensions.append(len(exp))
        self.dim=dimensions
        
    
    def print_function(self):
        """
        Creates a polymer
        pijx^iy^j with i of i=0,...,n, j=0,...,m
        
        Args:
            params: coefficients
            data: contains the the two independent variables x and y
            
            pijx^iy^j
        
        """
        ndim,mdim=self.dim
        for i,n in enumerate(self.exponents[0]):
            for j,m in enumerate(self.exponents[1]):
                
                #Getting the n,m dependent coefficients and exponents
                coeff_n=coeff(self.func_coeff[0],n)
                coeff_m=coeff(self.func_coeff[1],m)
                x_exp=coeff(self.func_exp[0],n)
                y_exp=coeff(self.func_exp[1],m)
                print '%s  %s c_{%s %s} x^{%s} y^{%s} +'%(coeff_n,coeff_m,n,m,x_exp,y_exp)
                
#    def print_initial_msg(self):
#        print "The exponents in x are:%s" %self.exponents[0]
#        print "The exponents in y are:%s\n" %self.exponents[1]
     
        
def coeff(function,point):

    coefficient=np.polyval(function,point)
    
    return coefficient


def get_list(a):
    """
    Converts to a list if it is not
    """
    if isinstance(a,list):
        return a
    else:
        return [a]
    
    
def delete_values(vector,delete_val):
    """
    Delete values from a vector 
    
    args:
        vector with all values from which you are going to delete
        delete_val values to delete from value
    """
    indexes=[]
    
    for val in delete_val:
        indexes.append(np.where(vector==val)[0])
    
    vector=np.delete(vector,np.array(indexes))
    
    return vector

def arbitrary_poly(data, *params):
    """
    Creates a polymer
    p(x,y)=\sum_{i,j}^{n,m} [ fn(i) fm(j) c_{i,j} x^{f_expx(i)} y^{f_expy(i)}
    
    Args:
        params: all the fitting parameters
        data: contains the the two independent variables x and y and an instance of the polynomial class containing all the information of it

    """
    
    points=data[0]
    x=points[0]
    y=points[1]
    poly=data[1]
    ndim,mdim=poly.dim
    params=np.reshape(params,(ndim,mdim))
    function=0
    
#    print 'Inside arbitraty poly %s %s'%(np.shape(x),np.shape(y))
    
    for i,n in enumerate(poly.exponents[0]):
        for j,m in enumerate(poly.exponents[1]):
            
            #Getting the n,m dependent coefficients and exponents
            coeff_n=coeff(poly.func_coeff[0],n)
            coeff_m=coeff(poly.func_coeff[1],m)
            if coeff_m==0 or coeff_n==0:
                function+=0
            else: 
                x_exp=coeff(poly.func_exp[0],n)
                y_exp=coeff(poly.func_exp[1],m)
                function+=params[i,j]*coeff_n*coeff_m*x**(x_exp)*y**(y_exp)
    return function

def poly_sum(data,*params):
    """
    Adds two arbitrary polynomials"""
    points=data[0]
    poly=data[1]
    res1=arbitrary_poly([points,poly[0]],params)
    res2=arbitrary_poly([points,poly[1]],params)
    return res1+res2

def two_poly(data,*params):
    
    """
    data contains[data1+data2,poly1+poly2] + means appended
    """
    data_1=data[:,:len(x_p)]  #Correct because x_e is global
    poly_1=poly_p
    data_2=data[:,len(x_p):]
    poly_2=poly_e
    z1=poly_sum([data_1,poly_1],params)
    z2=arbitrary_poly([data_2,poly_2],params)
    
    return np.append(z1,z2)


def fit_two_poly(data_1,data_2):
    """
    data_i contains the [[x,y,z,zerr],poly]
    """
    
    #Organising the dat in a suitable way
    x=np.append(data_1[0][0],data_2[0][0])
   
    
    y=np.append(data_1[0][1],data_2[0][1])
    z=np.append(data_1[0][2],data_2[0][2])
    

    zerr=np.append(data_1[0][3],data_2[0][3])
    poly=[data_1[1],data_2[1]]
    variables=np.stack((x,y))
    
    ndim,mdim=poly[1].dim
    
    popt, pcov = curve_fit(two_poly,variables[:,:],z,sigma=zerr,p0=[0]*ndim*mdim)
    popt_matrix=np.reshape(popt,(ndim,mdim))
    
    
    return popt_matrix,pcov,variables

def fit_poly(x,y,z,zerr,poly):
    """
    fits a polynomial using least square fitting
    the coefficients are pijx^i*y^j
    Args:
        x array containing the first variable 
        y array containing the second variable
        z array with the dependent variable
        zerr sigma on the dependent variable
        d degree of the polynomial
    
    Returns:
        popt: Optimal fitting coefficients.
        pcov: Covariant matrix of the fitting.
    """
    ndim,mdim=poly.dim
    variables=np.stack((x,y),axis=0)
    print np.shape(variables)
    popt, pcov = curve_fit(arbitrary_poly, [variables[:,:],poly], z, sigma=zerr,p0=[0]*ndim*mdim)
    popt_matrix=np.reshape(popt,(ndim,mdim))
    return popt_matrix,pcov,variables

def outputs(popt_matrix,pcov,e_results,n,m,name):
    
    print "\nCreated coefficients.dat containing all the fitting coefficients"
    np.savetxt('coefficients.dat', popt_matrix)
    print "\nCreated covariant.dat with the covariant matrix of the fitting"
    np.savetxt('covariant.dat', pcov)
    print "\nCreated error.dat containing the relative error between the property and the prediction given by the fitting evaluated at the same input points"
    header='%s\t%s_predicted\trelative_error'%(name,name)
    np.savetxt('error_%s.dat'%name,e_results, header=header,fmt='%12.8f')


def test_prediction(popt,variables,z,poly):
    """
    
    Args:
        popt is the fitting coefficient matrix
    
    Returns:
        Results with [z, z_predict, error]
    """
    m,n_point=np.shape(variables)
    z_predict=[]
    popt=np.reshape(popt,(np.size(popt)))
    
    for i in xrange(n_point):
        if np.size(poly)==1:
            z_predict.append(arbitrary_poly([variables[:,i],poly],popt))
        else:
            z_predict.append(poly_sum([variables[:,i],poly],popt))
        
#        print z_predict,z[i]
    z_predict=np.array(z_predict)

    error=np.abs(z-z_predict/z)*100
    
    results=np.transpose(np.vstack((z,z_predict,error)))
    
    return results
    

# =============================================================================
# Things not belonging to the class (This can be adapted to the problem)
# =============================================================================
    

def p_known(x,y):
    """
    known function for testing porpouses
    """
    z=3*x*y**2+4*y*x**2+5*x**2
    return z


def build_example(n_points=1000):
    """
    creates an example data file with the structure of the input
    """

    x=np.linspace(1,3,n_points)
    y=np.linspace(1,3,n_points)
    z=p_known(x,y)
    zerr= np.random.rand(n_points)
    
    data=np.column_stack([x,y,z,zerr])
    
    header='# density Temperature property sigma_property'
    np.savetxt('input_example.dat',data, header=header)
    
    
def build_example_grid(n_points=20):
    x=np.linspace(1,3,n_points)
    y=np.linspace(1,3,n_points)
    x,y=np.meshgrid(x,y)
    
    z=p_known(x,y)
    zerr= np.random.rand(*np.shape(x))
    
    
    data=np.column_stack([x.flatten(),y.flatten(),z.flatten(),zerr.flatten()])
    
    header='# density Temperature property sigma_property'
    np.savetxt('input_example_grid.dat',data, header=header)
    
def read_data(fname,prop_ref):
    """
    Reads the data from the input file that has [Temperature rho Property errorProperty]
    Args:
        fname: name of the file containing the data, run build example to see the structure 
        of the input_example.dat
        rho_ref: reference density
        beta_ref: reference beta
        prop_ref: Property at the reference point
        
    Returns:
        variables
    """
    data=np.loadtxt(fname)
    
    rho=data[:,1]
    temperature=data[:,0]
    prop=data[:,2]-prop_ref
    sigma_prop=data[:,3]
    beta=(1/temperature)
    
    return rho,beta,prop,sigma_prop


    

# =============================================================================
# MAIN
# =============================================================================


    
#def main(rho_ref,beta_ref,deg_x,deg_y,p_ref,e_ref,file_p,file_e):
#    '''
#    x is the density
#    y is beta
#    z is the property
#    '''
    


rho_ref=1.2
beta_ref=1/1.1
deg_x=6
deg_y=6
file_p='SP'
file_e='SE'
p_ref=17.67
e_ref=-4.66

# =============================================================================
#     Here is where I define the polynomial for the pressure
# =============================================================================
#    print "Basic information about the polynomial for the Pressure"
poly_p=[]
poly_p.append(polynomial(deg_x,deg_y,[1,-1],[1],[1,0],[1,0]))
poly_p.append(polynomial(deg_x,deg_y,[1,0],[1],[1,-1],[1,0])) 
#    poly_p.print_initial_msg()


# =============================================================================
#     Here is where I define the polynomial for the pressure
# =============================================================================
#    print "Basic information about the polynomial for the Energy"
poly_e=polynomial(deg_x,deg_y,[1],[1,0],[1,0],[1,-1])
#    poly_e.print_initial_msg()


print "rho_ref = %s"%rho_ref
print "beta_ref = %s"%beta_ref

# =============================================================================
#     #Reading the data for the pressure
# =============================================================================
x_p,y_p,z_p,zerr_p=read_data(file_p,p_ref)
x_p=x_p-rho_ref
y_p=y_p-beta_ref



# =============================================================================
#     #Reading the data for the energy
# =============================================================================
x_e,y_e,z_e,zerr_e=read_data(file_e,e_ref)
x_e=x_e-rho_ref
y_e=y_e-beta_ref


#Simultaneous Fitting
popt, pcov,variables=fit_two_poly([[x_p,y_p,z_p,zerr_p],poly_p],[[x_e,y_e,z_e,zerr_e],poly_e])


plt.close('all')



# =============================================================================
#     #Testing for Pressure
# =============================================================================
print ('\nTesting the pressure results\n')
variables_p=np.stack((x_p,y_p),axis=0)
er_results_p=test_prediction(popt,variables_p,z_p,poly_p)

outputs(popt,pcov,er_results_p, deg_x, deg_y,'p')

fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')
ax.scatter(x_p, y_p, z_p, zdir='z',marker='.',label="Simulation",color='r')





#Creating the surface
x,y=np.meshgrid(np.linspace(np.min(x_p),np.max(x_p),20),np.linspace(np.min(y_p),np.max(y_p),20))



z=np.asarray(x)
variables_p=np.stack((x.flatten(),y.flatten()),axis=0)
er_results=test_prediction(popt,variables_p,z.flatten(),poly_p)
z_mesh=np.reshape(er_results[:,1],np.shape(x))
ax.plot_wireframe(x,y,z_mesh,color='b')

ax.set_xlabel(r'$\rho-\rho^*$')
ax.set_ylabel(r'$\beta-\beta^*$')
ax.set_zlabel(r'$\Delta \beta P_{exc}$')
fig1.legend()
fig1.show()
fig1.savefig("3Dplot_p.pdf")


# =============================================================================
#     #Testing for Energy
# =============================================================================
print ('\nTesting the energy results\n')
variables_e=np.stack((x_e,y_e),axis=0)
er_results_e=test_prediction(popt,variables_e,z_e,poly_e)
#    
outputs(popt,pcov,er_results_e, deg_x, deg_y,'e')
#    
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(x_e, y_e, z_e, zdir='z',marker='.',label="Simulation",color='r')
#    ax.scatter(x, y, z_predict, zdir='z',label="Fitting",color='black')
#    
#Creating the surface
x,y=np.meshgrid(np.linspace(np.min(x_e),np.max(x_e),20),np.linspace(np.min(y_e),np.max(y_e),20))



z=np.asarray(x)
variables_e=np.stack((x.flatten(),y.flatten()),axis=0)
er_results=test_prediction(popt,variables_e,z.flatten(),poly_e)
z_mesh=np.reshape(er_results[:,1],np.shape(x))
ax2.plot_wireframe(x,y,z_mesh,color='b')
ax2.set_xlabel(r'$\rho-\rho^*$')
ax2.set_ylabel(r'$\beta-\beta^*$')
ax2.set_zlabel(r'$\Delta \rho \ e_{exc}$')

fig2.legend()
fig2.show()
fig2.savefig("3Dplot_e.pdf")





# =============================================================================
# #    Fitting only Pressure
# =============================================================================
#set_plot_appearance()
#popt1,pcov1,variables1=fit_poly(x_p,y_p,z_p,zerr_p,poly_p)
#single_results_p=test_prediction(popt1,variables1,z_p,poly_p)
#
#
#popt2,pcov2,variables2=fit_poly(x_e,y_e,z_e,zerr_e,poly_e)
#single_results_e=test_prediction(popt2,variables2,z_e,poly_e)



def plot_slices():
    # Plot slices in beta
    
    dir_name="slices_beta_p"
    shutil.rmtree(dir_name,ignore_errors=True) 
    
    os.mkdir(dir_name)
    slices=np.unique(y_p)
    for sli in slices:
        
        fig3=plt.figure()
        ax3=fig3.add_subplot(111)
        error_cap=4
        indexes=np.where(y_p==sli)
        ax3.errorbar(x_p[indexes],z_p[indexes],yerr=zerr_p[indexes],fmt='o',capsize=error_cap,label='beta=%s'%sli)
        ax3.plot(x_p[indexes],er_results_p[indexes[0],1],label='Fit 2')
    #    ax3.plot(x_p[indexes],single_results_p[indexes[0],1],label='Fit 1')
        fig3.legend(loc='upper right')
        fig3.tight_layout()
        fig3.savefig("%s/slice_%s.pdf"%(dir_name,sli))
        
        
    dir_name="slices_beta_e"
    shutil.rmtree(dir_name,ignore_errors=True) 
    os.mkdir(dir_name)
    
    slices=np.unique(y_e)
    for sli in slices:
        
        fig3=plt.figure()
        ax3=fig3.add_subplot(111)
        error_cap=4
        indexes=np.where(y_e==sli)
        ax3.errorbar(x_e[indexes],z_e[indexes],yerr=zerr_e[indexes],fmt='o',capsize=error_cap,label='beta=%s'%sli)
    
        ax3.plot(x_e[indexes],er_results_e[indexes[0],1],label='Fit 2')
    #    ax3.plot(x_e[indexes],single_results_e[indexes[0],1],label='Fit 1')
        fig3.legend(loc='upper right')
        fig3.tight_layout()
        fig3.savefig("%s/slice_%s.pdf"%(dir_name,sli))


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
    

