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
    for i,n in enumerate(poly.exponents[0]):
        for j,m in enumerate(poly.exponents[1]):
            
            #Getting the n,m dependent coefficients and exponents
            coeff_n=coeff(poly.func_coeff[0],n)
            coeff_m=coeff(poly.func_coeff[1],m)
            x_exp=coeff(poly.func_exp[0],n)
            y_exp=coeff(poly.func_coeff[1],m)
            function+=params[i,j]*coeff_n*coeff_m*x**(x_exp)*y**(y_exp)
    return function

def test_prediction(popt,variables,z,poly):
    """
    Returns the points generated by the predicted function
    Args:
        popt is the fitting coefficient matrix
    """
    m,n_point=np.shape(variables)
    z_predict=[]
    popt=np.reshape(popt,(np.size(popt)))
    for i in xrange(n_point):
        z_predict.append(arbitrary_poly([variables[:,i],poly],popt))
    z_predict=np.array(z_predict)

    error=np.abs(z-z_predict)/z
    
    return z_predict,error


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
    popt, pcov = curve_fit(arbitrary_poly, [variables[:,:],poly], z, sigma=zerr,p0=[1]*ndim*mdim)
    popt_matrix=np.reshape(popt,(ndim,mdim))
    return popt_matrix,pcov,variables

def outputs(popt_matrix,pcov,error,n,m):
    
    print "\nCreated coefficients.dat containing all the fitting coefficients"
    np.savetxt('coefficients.dat', popt_matrix)
    print "\nCreated covariant.dat with the covariant matrix of the fitting"
    np.savetxt('covariant.dat', pcov)
    print "\nCreated error.dat containing the relative error between the property and the prediction given by the fitting evaluated at the same input points"
    np.savetxt('error.dat',error)


# =============================================================================
# things not belonging to the class (This can be adapted to the problem)
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
    print np.shape(x),np.shape(y),np.shape(z),np.shape(zerr)
    
    
    data=np.column_stack([x.flatten(),y.flatten(),z.flatten(),zerr.flatten()])
    
    header='# density Temperature property sigma_property'
    np.savetxt('input_example_grid.dat',data, header=header)
    
def read_data(fname,rho_ref,beta_ref):
    """
    Reads the data from the input file
    Args:
        fname: name of the file containing the data, run build example to see the structure 
        of the input_example.dat
        rho_ref: reference density
        beta_ref
        
    Returns:
        variables
    """
    global rho, beta
    data=np.loadtxt(fname)
    
    
    rho=data[:,1]-rho_ref
    temperature=data[:,0]
    prope=data[:,2]
    sigma_prope=data[:,3]
    beta=1/temperature-beta_ref
    
    return rho,beta, prope, sigma_prope





        

    


    
    

# =============================================================================
# MAIN
# =============================================================================

def print_initial_msg(poly):
    print '\nRunning the script assuming:'
    print "The exponents in rho are:%s" %poly.exponents[0]
    print "The exponents in beta are:%s" %poly.exponents[1]
    
def main(input_file,rho_ref,beta_ref,deg_x,deg_y):
    
    global x,y,z,z_mesh,popt,poly_p,poly_e
    
    poly_p=polynomial(deg_x,deg_y,[1,-1],[1],[1,0],[1,0])
    poly_e=polynomial(deg_x,deg_y,[1],[1,0],[1,0],[1,-1])
    print_initial_msg(poly_p)
    
    print_initial_msg(poly_e)
    
    print "rho_ref = %s"%rho_ref
    print "beta_ref = %s"%beta_ref
    
    
#    x,y,z, zerr=read_data(
#            input_file,rho_ref,beta_ref)
#    popt, pcov,variables=fit_poly(x,y,z,zerr,poly)
#    
#    z_predict,error=test_prediction(popt,variables,z,poly)
#    
#    outputs(popt,pcov,error, deg_x, deg_y)
#    
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.scatter(x, y, z, zdir='z',marker='.',label="Simulation",color='r')
##    ax.scatter(x, y, z_predict, zdir='z',label="Fitting",color='black')
#    
#    #Creating the surface
#    x,y=np.meshgrid(np.linspace(np.min(x),np.max(x),20),np.linspace(np.min(y),np.max(y),20))
#    z=np.asarray(x)
#    variables=np.stack((x.flatten(),y.flatten()),axis=0)
#    z_mesh,error=test_prediction(popt,variables,z.flatten(),poly)
#    z_mesh=np.reshape(z_mesh,np.shape(x))
#    ax.plot_wireframe(x,y,z_mesh,color='b')
#    ax.set_xlabel("rho-rho*")
#    ax.set_ylabel("beta-beta*")
#    
#    fig.legend()
#    fig.savefig("3Dplot.pdf")
    
    




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script fits a 2 variable data to a poly of deg n',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file_name', metavar='input file',help='Input filename', type=str)
    parser.add_argument('-beta_ref', metavar='beta reference',help='reference beta',default=0,type=float)
    parser.add_argument('-rho_ref',metavar='rho reference',help='reference rho',default=0,type=float )
    parser.add_argument('-deg_x',metavar='poly degree in rho',help='Degree of the poly, or min->max degree',nargs='+', default=[0,3],type=int)
    parser.add_argument('-deg_y',metavar='poly degree in beta',help='Degree of the poly',nargs='+',default=[2],type=int)
    parser.add_argument('-exc_x',metavar='Exponents to exclude in x',help='Exponents to exclude in x as a list i.e 0 3 4',nargs='+',type=int)
    parser.add_argument('-exc_y',metavar='Exponents to exclude in y',help='Exponents to exclude in y as a list i.e 0 3 4',nargs='+',type=int)
    
    args = parser.parse_args()
    
    
    main(args.file_name,args.rho_ref,args.beta_ref,args.deg_x,args.deg_y)
    
    
    





##Diagonal matrix test
#
#n=4 
#vec=np.random.rand(4)
#mat=np.diag(vec)
#inv=np.linalg.inv(mat)
#
#mat_a=np.random.rand(4,4)
#diag=np.diagonal(mat_a)

