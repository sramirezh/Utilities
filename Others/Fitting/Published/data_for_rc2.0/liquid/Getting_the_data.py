#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 12:50:23 2020
This script is to predict the energy data from the fitted coefficients in the SI for 
the GLJ paper
@author: simon
"""
import numpy as np
import pandas as pd
import multid_fit_liquid_jup as mf
import matplotlib.pyplot as plt


coefficient_file = 'energy/fit_coefficients_energy_pressure_liquid_rc2.0.txt'
values_file = 'energy/energy_liquid_rc2.0.txt'



    
def read_file(file_name):
    """
    Reads a file assuming that the the first line is the header
    
    Args:
        file_name: Name of the file
        
    Returns:
        header: the header as an array of strings
        data: is a pandas dataframe
    """
    
    with open(file_name, 'r') as data_file:
        header = data_file.readline().split()
        data = pd.read_csv(data_file,delim_whitespace=True, header=None).dropna(axis=1, how='all')
        data =data.to_numpy()
    return header, data
    

def read_coeff_file(input_file):
    """
    Reads the fitting coefficients and the exponents for the polymer fit
    
    """
    
    header, data = read_file(input_file)
    m_values = [int(el.split('=')[-1]) for el in header]
    n_values = [int(el.split('=')[-1]) for el in data[:,0]]
        
    coefficients = data[:,1:].astype('float')
        
    
    return m_values, n_values, coefficients
    




def poly_coeff(func, point):
    """
    Returns the function f1(n), f2(m), f3(n) and f4(m) that appear in my notebook
    
    function is the list of exponents as described in polyval
    point is the value of n or m
    
    """

    f = np.polyval(func, point)
    
    return f


def arbitrary_poly(data, *params):
    """
    evaluates the polynomial
    p(x,y)=\sum_{i,j}^{n,m} [ fn(i) fm(j) c_{i,j} x^{f_expx(i)} y^{f_expy(i)}
    
    at at the point x, y 
    
    Args:
        params: all the fitting coefficients a_nm
        data: contains the the two independent variables x and y and an instance of the polynomial class containing all the information of it

    """
    
    points = data[0]
    x = points[0]
    y = points[1]
    poly = data[1]
    ndim,mdim = poly.dim
    params = np.reshape(params,(ndim,mdim))
    f_eval = 0
    
    
#    print 'Inside arbitraty poly %s %s'%(np.shape(x),np.shape(y))
    
    for i,n in enumerate(poly.exponents[0]):
        for j,m in enumerate(poly.exponents[1]):
            
            #Getting the n,m dependent coefficients and exponents
            f1_n = poly_coeff(poly.func_coeff[0],n)
            f2_m = poly_coeff(poly.func_coeff[1],m)
            if f1_n  == 0 or f2_m  == 0:
                pass
            else: 
                f3_n = poly_coeff(poly.func_exp[0],n)
                f4_m = poly_coeff(poly.func_exp[1],m)
                f_eval += params[i,j] * f1_n * f2_m * x ** (f3_n) * y ** (f4_m)
                
        
    return f_eval

def test_prediction(popt,variables,z,poly):
    """
    
    Args:
        popt is the fitting coefficient matrix
        variables: should contain the variables x = temperature and y =density in an array 2xNumber of points
                It should be TRANSPOSED
    
    Returns:
        Results with [z, z_predict, error]
    """
    m,n_point=np.shape(variables)
    z_predict=[]
    popt=np.reshape(popt,(np.size(popt)))
    for i in range(n_point):
        z_predict.append(arbitrary_poly([variables[:,i],poly],popt))
#        print z_predict,z[i]
    z_predict=np.array(z_predict)

    error = np.abs((z - z_predict) / z) * 100
    
    results=np.transpose(np.vstack((z, z_predict, error)))
    
    return results




# Reading the coefficients and exponents
m_values, n_values, coefficients = read_coeff_file(coefficient_file)

# Reading the energy
header, data = read_file(values_file)



poly_e = mf.polynomial(n_values,m_values,[1],[1,0],[1,0],[1,-1])




x_e = data[:,0]
y_e = data[:,1]
z_e = data[:,2]


# Getting the volume of the system

n_part = 1000
vol = n_part/y_e[0]



variables = np.copy(data[:,0:2])
# Getting beta
variables[:,0] = 1/variables[:,0]



first_point = arbitrary_poly([variables[1,:],poly_e],coefficients)

er_results_e = test_prediction(coefficients, np.transpose(data[:,0:2]), data[:,2], poly_e)

#mf.set_plot_appearance()
#
#fig2 = plt.figure()
#ax2 = fig2.add_subplot(111, projection='3d')
#ax2.scatter(x_e, y_e, z_e, zdir='z',marker='.',label="Simulation",color='r')
#ax2.set_xlabel(r'$T$')
#ax2.set_ylabel(r'$\rho$')
#ax2.set_zlabel(r'$E$')



# Testing the polynomial by hand



