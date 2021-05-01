# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 17:22:14 2019
Contains the classes flux and correlation that help to compute correlations
@author: sr802
"""
import multiprocessing
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed
import pickle as pickle
from scipy.stats import sem
import uncertainties as un
from uncertainties import unumpy  # Must be imported like this
from statsmodels.tsa.stattools import acf, ccf
import time as t
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf

# =============================================================================
# Class definition
# =============================================================================

print(os.getcwd())
class flux(object):
    """
    Flux is a vectorial or scalar entity
    """

    def __init__(self, components, times, name):
        """
        Args:
            Components: is a matrix( or vector) containing the time series of 
            the components in [x,y,z]
            Name: Is the name that is going to appear in the plots in latex 
            format, example "r'J_s-c_s^BQ'"
        """

        self.components = components
        self.name = name
        self.times = times
        self._reshape()

    def _reshape(self):
        """
        If it is a single column, it will reshape it into a matrix with one 
        column
        """

        shape = np.shape(self.components)

        if len(shape) == 1:
            self.components = np.reshape(self.components, 
                                         (len(self.components), 1))
            self.dimension = 1

        else:
            self.dimension = shape[1]
        

class correlation(object):
    """
    TODO, I could compute all the correlations for the given delta just with
    one np.cov

    Attributes:
        flux1: Instances from the class flux
        flux2: Instances from the class flux
        max_delta: is the maximum tau to compute correlations <v(\tau) v(0)>
        cor: contains the correlations as an array, the components
            are x,y,z, and total in 3D
        norm: normalisation factor for each dimension  t=0 (var1(0) var2(0))
        cor_norm: list to store the normalised correlations, the last is the
            total
        half_max: half of the time steps, in order to sample all the correlation
                  tau with the same number of samples, as for the 
                  diffusion coefficient.
                such that <v(0)v(tau)> + <v(1)v(tau+1)> + ...
                <v(half_max)v(tau+half_max)>

                where half_max is (max time step)/2 
    """
    def __init__(self, flux1, flux2, max_delta):
        """
        Args:
            flux1 and flux2 are instances of the class flux
            max_delta is the index of the maximum \tau for the analysis
            
        """
        self.flux1 = flux1
        self.flux2 = flux2
        self.max_delta = max_delta   
        self.initial_check()
        dimension = self.dimension
        self.cor = (dimension + 1) * [0]  
        self.norm = (dimension + 1) * [0]  
        self.cor_norm = (dimension + 1) * [0]  
        self.dic_label = {0: r'$x$', 1: r'$y$', 2: r'$z$', self.dimension: 'Total'}
        
    def initial_check(self):
        """
        Checks if the time series are equal and the fluxes have the same
        number of components
        """
        self.dimension = self.flux1.dimension
        if self.flux1.dimension != self.flux2.dimension:
            print("The fluxes do not have the same dimension")
        else:
            self.dimension = self.flux1.dimension
            
        if self.flux1.times.all != self.flux2.times.all:
            print("The fluxes were not measured for the same times")
        else:
            self.times = self.flux1.times[:self.max_delta]
        
        self.half_max = int(len(self.flux1.times)*0.5)

            
    def correlate_one_d(self, dim):
        """Performs a correlation between 1d components by evaluating the
        products at at different delta t.
        Args:
            dim is the component to be evaluated
        """
        num_cores = multiprocessing.cpu_count()
        var1 = self.flux1.components[:, dim]
        var2 = self.flux2.components[:, dim]

        # cor = []
        # for i in range(self.max_delta):
        #     cor.append(compute_correlation_dt(var1, var2, i))


        cor = Parallel(n_jobs=num_cores)(delayed(compute_correlation_dt)
                                         (var1, var2, i) for i in tqdm(range(self.max_delta)))
        norm = cor[0]
        self.norm[dim] = norm
        self.cor[dim] = np.array(cor)
        self.cor_norm[dim] = np.array(cor) / norm
   
    def evaluate(self):
        """
        DO NOT USE, THE CORRELATE USING COV does not work! see my logbook on Pore
        or my autocorrelation_and_GK jupyter
        Performs the correlations of the 1d components,calling
        correlate_one_d, and adds them up to the total.
        
        """
        total = np.zeros(self.max_delta)
        for dim in range(self.dimension):
            self.correlate_one_d(dim)
            total = total + self.cor[dim]
        total = total / 3
        self.cor[-1] = total
        self.norm[-1] = total[0]
        self.cor_norm[-1] = total / total[0]

    def evaluate_acf(self):
        """
        TODO: Generate an error if the fluxes are not the same
        Performs the Autocorrelation of the 1d components,using 
        statsmodels.tsa.stattools.acf
        See my ipython about autocorrelation and GK

        """
        total = np.zeros(self.max_delta)
        for dim in range(self.dimension):
            print("started the analysis for %s"%dim)
            t0 = t.time()
            x = self.flux1.components[:, dim]
            acf_array, confidence = acf(x, nlags = len(self.times[:self.max_delta])-1, fft = True, alpha = 0.05)
            # see autocorrelation_and_gk.ipynb for explanation
            std_error = (acf_array - confidence[:, 0]) / 2  
            amplitude = np.correlate(x, x) / len(x)
            self.norm[dim] = amplitude
            cor_norm = unumpy.uarray(acf_array, std_error)   
            self.cor_norm[dim] = cor_norm
            self.cor[dim] = amplitude * cor_norm

            total = total + self.cor[dim]
            print(t.time() - t0)
        
        total = total / 3
        self.cor[-1] = total
        self.norm[-1] = total[0]
        self.cor_norm[-1] = total / total[0]

    def evaluate_ccf(self):
        """
        TODO it is very slow. It should be improved following the source of the
        acf(statsmodels.tsa.stattools) which evaluates the acovf, tha uses
        fft to compute the correlation.
        TODO: Include the error
        Performs the correlation of the 1d components,using 
        statsmodels.tsa.stattools.ccf

        ----------------------------------------------
        Notes: See my ipython about autocorrelation and GK
        
        """
        total = np.zeros(self.max_delta)
        for dim in range(self.dimension):
            print ("started the analysis for %s"%dim)
            t0 = t.time()
            x = self.flux1.components[:, dim]
            y = self.flux2.components[:, dim]
            ccf_array = ccf(x, y)[:len(self.times)]
            amplitude = np.correlate(x,y)/len(x)
            self.norm[dim] = amplitude
            self.cor[dim] = amplitude * ccf_array
            self.cor_norm[dim] = ccf_array

            total = total + self.cor[dim]
            print (t.time()-t0)
        
        total = total / 3
        self.cor[-1] = total
        self.norm[-1] = total[0]
        self.cor_norm[-1] = total / total[0]
        
    def plot_individual(self, ax, dim=0, alpha=0.4, every=1, norm=True):
        """
        Plots the correlation for the given dimension
        Args:
            ax axes object
            fig Figure
            dim is the dimension, for example:in a 3D vector, 0-x, 1-y, 2-z
            and 3-total.
            alpha is the transparency of the filling
            every to not have so many points
            norm True if normalised
            The axis label is given here but it could be renamed later
        """
        if norm == True:
            cor = self.cor_norm[dim]
        else:
            cor = self.cor[dim]
            
        # checking the first entry to see if it is a ufloat   
        if isinstance(cor[0], un.UFloat):
    
            y = np.array([i.n for i in cor])
            y_error = np.array([i.s for i in cor])
            
            # Getting time, y, error, with the frequency
            times = self.times[::every]
            y = y[::every]
            y_error = y_error[::every]
            ax.plot(times, y, label=self.dic_label[dim])
            ax.fill_between(times, y - y_error, y + y_error, alpha=0.4)
        
        else:
            ax.plot(self.times[::every], cor[::every], 
                    label=self.dic_label[dim])
            
        return ax
    
    def plot_all(self, ax, alpha=0.4, norm=True):
        """
        Plots the correlation in all dimensions
        """
        for dim in range(self.dimension + 1):
           ax = self.plot_individual( ax, dim, norm=norm)
        
        ax.axhline(y=0, xmin=0, xmax=1, ls=':', c='black')
        ax.set_ylabel(r'$\langle %s(t)%s(0) \rangle$'%(self.flux1.name, 
                                                       self.flux2.name))
        ax.set_xlabel('time')
        return ax
    
    def save(self, file_name):
        """
        Saves the instance
        """
        afile = open(r'%s.pkl'%file_name, 'wb')
        pickle.dump(self, afile)
        afile.close()

    def transport_coeff_comp(self, xmin, xmax, component):
        """
        Returns the integral of the correlation in one component
        Args:
            xmin: lower limit of the integral
            xmax: upper limit of the integral
            component: could be 0,1,2,-1, ie. x,y,z,total
        Returns:
            integral: of the correlation
        """
        
        integral = cf.integrate(self.times, self.cor[component], xmin, xmax)

        return integral

    def transport_coeff(self, xmin, xmax):
        """
        computes the self.coeff, transport coefficients in all directions without 
        multiplying for any prefactor, so, only the integration of the correlation

        Args:
        ----------
            xmin: lower limit of the integral
            xmax: upper limit of the integral
        Returns:
        ----------
            Nothing
        """
        coeff = []
        if self.dimension > 1:
            for dim in range(self.dimension+1):
                integral = self.transport_coeff_comp(xmin, xmax, dim)
                coeff.append(integral)
            self.coeff = np.array(coeff)
        else:
            self.coeff = cf.integrate(self.times, self.cor, xmin, xmax )

        return self.coeff


class bundle_correlation(correlation):
    """
    The bundle is made of an array of correlations (correlation.cor), it is not
    a bundle of correlation objects as this might become very expensive.

    """
    def __init__(self, array, times, flux1_name, flux2_name):
        """
        Args:
            array: contanining all the correlations 
            times: times at which the correlations where measured. Couls be 
                    taken from one of the correlation instances as 
                    instance.times
            flux1_name: Name of the flux 1, also can be taken from the instance
            flux2_name: Name of the flux 2, also can be taken from the instance
        
        
        Atributes:
            cor: contains the correlations as an array, the components
            are x,y,z, and total in 3D
            norm: normalisation factor for each dimension  t=0 (var1(0) var2(0))
            cor_norm: list to store the normalised correlations, the last is the
            total
        """
        self.arrays = np.array(array)
        self.flux1_name = flux1_name
        self.flux2_name = flux2_name
        self.times = times
        self.initial_computes()
        self.dic_label = {0: r'$x$', 1: r'$y$', 2: r'$z$', -1: 'Total', self.dimension:'Total'}
        
        # self.norm = max(self.cor).nominal_value
        # self.cor_norm = [self.cor / self.norm]
        
        
    def initial_computes(self):
        """
        Gets the shape of the correlation array and calls the averaging, 
        Args:
            self

        Returns
        -------
        None.

        """
        
        self.n_runs, self.dimension, self.n_points = np.shape(self.arrays)
        # Initilising the lists
        self.norm = (self.dimension ) * [0]  
        self.cor_norm = (self.dimension ) * [0]  
        self.averaging()
        
        for dim in range(self.dimension):
            normalisation = self.cor[dim][0]
            self.norm[dim] = normalisation
            self.cor_norm[dim] = self.cor[dim]/normalisation

    def averaging(self):
        """
        np.shape(self.arrays) =[#runs, Dimensions(+average),points 
        Takes the average of the array of correlations.
        Sometimes the individual correlations have errors, therefore
        each point has (average + std)

        At the moment it strips the error from uncertainties and analyses
        everypoint as independent.

        Returns
        -------
        None.

        """
        ddof = 1
        array = self.arrays
        
        # To avoid errors if there is only one run
        if len(array) == 0:
            ddof = 0
        
        #Striping the uncertainties from individual runs
        if isinstance(array[0][0][0], un.UFloat):
            array = unumpy.nominal_values(array)
        
        
        self.cor = unumpy.uarray(np.average(array, axis=0), sem(array,
                                axis=0, ddof=ddof))
           
    def plot(self, fig, ax, dim = 0, alpha=0.4, every=1, ax_label=True,
             norm=True):
        """
        Args:
            ax: axes object
            fig: Figure
            dim: is the dimension,
            for example:in a 3D vector, 0-x, 1-y, 2-z and 3-total.
            alpha: is the transparency of the filling
            every: to not have so many points
            norm: True if normalised
            ax_label: The axis label is given here but it could be renamed
                      later

        Returns:
            fig, ax
        """
        if norm == True:
            cor = self.cor_norm[dim]
        else:
            cor = self.cor[dim]

        y = np.array([i.n for i in cor])
        y_error = np.array([i.s for i in cor])
        ax.plot(self.times[::every], y[::every])
        ax.fill_between(self.times, y - y_error, y + y_error, alpha=0.4)

        # It mostly means that the plot will not be further used, for example 
        # not used to compare correlations
        if ax_label == True:
            ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')
            ax.set_ylabel(r'$\langle %s(t)%s(0) \rangle$'%(self.flux1_name, 
                                                           self.flux2_name))
            ax.set_xlabel(r'$t$ ')

        return fig, ax

    # def plot_bundle(self, fig, ax, dim):
    #     """
    #     Plots all the correlations from the individual runs in the given direction

    #     Returns:
    #         fig, ax
    #     """

    #     for cor in self.arrays:
    #         ax.plot(self.times, cor[dim])
    #     ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')
    #     return fig, ax
        


    def plot_bundle(self, ax, dim, alpha=0.4, every = 1):
        """
        Plots all the correlations from the individual runs in the given direction
        
        TODO: Add the normalised option
        Args:
            ax axes object
            fig Figure
            dim is the dimension, for example:in a 3D vector, 0-x, 1-y, 2-z
            and 3-total.
            alpha is the transparency of the filling
            every to not have so many points
            The axis label is given here but it could be renamed later
        """
        
        
        for array in self.arrays:
            
            cor = array[dim]
            # checking the first entry to see if it is a ufloat   
            if isinstance(cor[0], un.UFloat):
        
                y = np.array([i.n for i in cor])
                y_error = np.array([i.s for i in cor])
                
                # Getting time, y, error, with the frequency
                times = self.times[::every]
                y = y[::every]
                y_error = y_error[::every]
                ax.plot(times, y, label = self.dic_label[dim])
                ax.fill_between(times, y - y_error, y + y_error, alpha=0.4)
            
            else:
                ax.plot(self.times[::every], cor[::every], 
                        label = self.dic_label[dim])
            
        return ax




# def average_unumpy_tensor(array):
#     """
#     Performs the average of unumpy values 
#     Parameters
#     ----------
#     array : np.array
#         it is an array that has [n,m,l] dimensions

#     Returns
#     -------
#     None.

#     """
    

def compute_correlation_dt(var1, var2, delta):
    """
    This is a VERY GENERAL function
    
    *****
    It is important to keep it outside the class as it is going to be 
    evaluated in parallel and if it is a method of the class, it would require
    to load the instance on each processor which could be very expensive
    *****
    
    Computes the correlation for two variables for a given delta t
    Args:
        var1: time series for the first variable 
        var2: time series for the first variable
        Note that var1 and var2 have to have the same time series
        delta: every this number of steps, we take the variable 2 to compute 
        the correlation
    """
    # cf.blockPrint()
    # Covariance need to be reduced as it returns a matrix that is why the 0, 1 
    # component is taken

    # Defining the half of the total of samples as in diffusion coefficient
    half_max = int(len(var1) * 0.5)

    if delta != 0:
        cor = np.cov(var1[:half_max], var2[delta:half_max + delta])[0][1]
    else:
        cor = np.cov(var1, var2)[0][1] 

    # cf.enablePrint()
    
    return cor
