#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:40:14 2020
Script to compute the Transport coefficients from DO simulations
it uses 1 output from lammps:

    fluxes.dat that contains TimeStep v_Qx v_Jsx_exc v_Qy v_Jsy_exc v_Qz v_Jsz_exc
    
and if new combinations of the fluxes are wanted, use run_analysis_vdata:
    vdata.dat 
    v_vx_Solv v_vx_Solu v_vx_Sol v_cSolu v_cSolv
    
Run inside 6.GJ_analysis

Define which folder pattern to use in 

folder_pattern eg "r*'


Run it as: gk_transport.py r1 -t "low"

Based on viscosity from DO project
@author: simon
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.stats import sem
from collections import defaultdict
import weakref
import argparse
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.lammps_utilities as lu
import Lammps.Pore.GK.flux_correlation as fc

class SimulationGK(lu.SimulationType):
    __refs__ = defaultdict(list)
    def __init__(self, name, max_tau, tau_integration, temp,
                 log_file, flux_file):
        """
        Args:
            name: Identifier of the simulation
            time_step: Simulation time_step
            max_tau: maximum length of the correlation analysis
            tau_integration: upper time limit of the correlation integration
            temp: system temperature
            log_file: lammps log file
            flux file: file containing all the fluxes
            
            
        Other atributes:
            sampling_interval
        """
        self.__refs__[self.__class__].append(weakref.ref(self))
        self.name = name
        self.time_step = []
        self.max_tau = max_tau
        self.tau_integration = tau_integration
        self.temp = temp
        self.log_file = log_file
        self.flux_file = flux_file
        self.flux_file_prefix = self.flux_file.split('.')[0]
        
    def get_thermo(self):
        """
        Gets the thermo data as an instance of lu.SImulation
        
        """
        
        self.thermo = lu.Simulation(self.log_file)
        

    def run_analysis(self):
        """
        Specific method for this problem
        This asssumes that the cwd is where all the files are
        Specific analysis where all the fluxes are their components are defined
        
        """
    
        data = cf.read_data_file(self.flux_file)
        
        data1 = data.values
        times = (data1[:,0]-data1[0,0]) * self.time_step
        
        
        # Assuming times are homogeneous
        self.sampling_interval = times[1]-times[0] 
    
        # int(len(times)*0.1) #Maximum delta of time to measure the correlation
        self.max_delta = int(self.max_tau/self.sampling_interval)  
        
        # Name of the components
        Q_components = ["v_Qx", "v_Qy"]
        if self.name == "stat":
            Js_components = ["v_Jsx_exc_cons", "v_Jsy_exc_cons"]
        else:
            Js_components = ["v_Jsx_exc", "v_Jsy_exc"]
        
        Q = data[Q_components].values
        Js = data[Js_components].values
        
        Q_flux = fc.flux(Q ,times, "Q")
        Js_flux = fc.flux(Js ,times, "J_s'")
        
        #  Creating the correlation instance
        m_qs = fc.correlation(Q_flux, Js_flux, self.max_delta)
        m_qs.evaluate()
        m_qs.save("Mqs_%s"%self.save_name)
        
        m_sq = fc.correlation( Js_flux,Q_flux, self.max_delta)
        m_sq.evaluate()
        m_sq.save("Msq_%s"%self.save_name)
        
        return m_qs, m_sq
    
    def run_analysis_vdata(self, vdata_file, vol_system, vol_bulk):
        """
        Specific method for this 
        This asssumes that the cwd is where all the files are
        Specific analysis where all the fluxes are their components are defined
        from a vdata.dat file
        
        args:
            vdata_file: name of the vdata file to analyse
            vol_system: is the omega, all the volume but far from the reflective
            wall
            
        """
        
        data =  cf.read_data_file(vdata_file)
        times = data['TimeStep']
        times = (times-times[0])* self.time_step
        
        # Assuming times are homogeneous
        self.sampling_interval = times[1]-times[0] 
    
        # int(len(times)*0.1) #Maximum delta of time to measure the correlation
        self.max_delta = int(self.max_tau/self.sampling_interval) 
        
        
        # Creating the fluxes 
        
        data['Qx'] = data['v_vx_Sol']
        data['Qy'] = data['v_vy_Sol']
        
        cs = data['v_cSolu']/vol_system # Concentration of solutes in vol
        
        Js_x = cs * data['v_vx_Solu']
        Js_y = cs * data['v_vy_Solu']
        
        pref = data['v_cBSolu']/(data['v_cBSolu']+data['v_cBSolv'])*data['v_rho_ave_sys']
        
        data['Jsx_exc'] = Js_x - pref * data['Qx']
        data['Jsy_exc'] = Js_y - pref * data['Qy']
        
        
        Q_components = ["Qx", "Qy"]
        Js_components = ["Jsx_exc", "Jsy_exc"]
        
        Q = data[Q_components].values
        Js = data[Js_components].values
        
        Q_flux = fc.flux(Q ,times, "Q")
        Js_flux = fc.flux(Js ,times, "J_s'")
        
        #  Creating the correlation instance
        m_qs = fc.correlation(Q_flux, Js_flux, self.max_delta)
        m_qs.evaluate()
        m_qs.save("Mqs_vdata%s"%self.save_name)
        
        m_sq = fc.correlation( Js_flux,Q_flux, self.max_delta)
        m_sq.evaluate()
        m_sq.save("Msq_vdata%s"%self.save_name)
        
        return m_qs, m_sq
        
    @property
    def save_name(self):
        """
        Returns the sufix to save the output files, if the flux file has a number
        adds this to the original name of the simulation
        """
        self.number = cf.extract_digits(self.flux_file)
        
        # If there are no numbers, the name does not change
        if len(self.number) == 0:
            self.number =''
        else:
            self.number = int(self.number[0])
        
        save_name =  '%s%s'%(self.name,self.number)
        return save_name
    
    @classmethod
    def get_instances(cls):
        """
        Gets all the instances created of this class
        https://stackoverflow.com/questions/328851/printing-all-instances-of-a-class
        """
        for inst_ref in cls.__refs__[cls]:
            inst = inst_ref()
            if inst is not None:
                yield inst


def remove_pkl(root):
    
    """
    Removes all the pkl inside the directory recursively
    TODO generalise and put in cf
    """
    files = glob.glob('%s/*.pkl'%root,  recursive = True)
    for file in files:
        os.remove(file)


def plot_coeff_vs_tau(transport, sim, correlation_dist):
    """
    Plots the evolution of the transport coefficient
    TODO this can be added to the correlation class
    
    """
    fig, ax = plt.subplots()
    
    times = transport.times
    tau_array = np.linspace([1], times[-1], 20)
    # Need to add the integration time
    tau_array = np.sort(np.append(tau_array, sim.tau_integration)) 
    
    coeff_array = []
    coeff_error = []
    
    for t in tau_array:
        coeff_t = [] 
        for cor in correlation_dist:
            # Taking only the nominal value
            coeff_t.append(cf.integrate(times, cor, 0, t,))
        
        coeff_array.append(np.average(coeff_t))
        coeff_error.append(sem(coeff_t))
    
    
    coeff_array = np.array(coeff_array)
    coeff_error = np.array(coeff_error)
    
    ax.plot(tau_array, coeff_array)
    ax.fill_between(tau_array, coeff_array - coeff_error, coeff_array + coeff_error, alpha=0.4)
    ax.set_xlabel(r'$t^*$')
    ax.axvline(x = sim.tau_integration,ls='-.',c='black')
    
    # Generating a table
    
    logger.info("Created a table with the coefficient vs tau in %s"%sim.plot_dir)
    header = "tau\tcoefficient\terror"
    
    data = np.transpose(np.stack((tau_array,coeff_array,coeff_error)))
    name = "Gamma_%s_%s"%(transport.flux1.name, transport.flux2.name)
    np.savetxt("%s/%s_%s.dat"%(sim.plot_dir, name, sim.name), data, header = header)
    
    return fig, ax


def main(pattern, sim_type):
    """
    Args:
        pattern: Folder pattern where all the simulations are
        sim_type: among all the simulation type definitions see here in the code.
    """
    
    for instance in SimulationGK.get_instances():
        if instance.name == sim_type:
            sim = instance.copy
    #sim = low.copy
    folder_pattern = pattern
    
    sim.plot_dir = plot_dir
    
    logger.info("\nUsing the folder pattern %s"%folder_pattern )
    folders = glob.glob(folder_pattern)
    logger.info("Analysing the results in %s"%folders )
    
    # Array to store all the correlations
    #m_qs
    correlation_dist = []
    transport_coeff = []
    #m_sq
    correlation_dist1 = []
    transport_coeff1 = []
    
    for i, folder in enumerate(folders):
        Kb = 1
        logger.info("Analysing the folder %s"%folder)
        # Computing this only once, assuming all runs have the same
        os.chdir(folder)
        partial_runs = glob.glob("%s*"%sim.flux_file_prefix)
        for j, file in enumerate(partial_runs):
            if i == 0 and j == 0: 
                
                # Getting some necessary parameters
                sim.get_thermo()
                lz_box = sim.thermo.limits[1,2] - sim.thermo.limits[0,2]
                lz_system = lu.read_value_from("in.geom", 'rSystem')[-1]
                lz_bulk = lu.read_value_from("in.geom", 'rBulk')
                lz_bulk = lz_bulk[-1] - lz_bulk[0]
                sim.time_step = lu.read_value_from("input.lmp", 'timestep')[0]
                sim.print_params(logger)
                
                # We need the volume where the fluxes were measured
                volume_system = sim.thermo.volume * lz_system/lz_box
                volume_bulk = sim.thermo.volume * lz_bulk/lz_box
                prefactor = volume_system / (Kb * sim.temp)
                
            # The name of the 
            sim.flux_file = file
            logger.info("\nAnalysing the file %s"%file)
            # =============================================================================
            # Computing the viscosity from the stress data
            # =============================================================================
            
            # Loading the correlations if they were already computed, this saves almost 50 minutes in 16 cores (dexter)
            if (len(glob.glob("Mqs_%s.pkl"%sim.save_name)) == 1):
                
                logger.info("\nThere is a pkl, we need to load Mqs_%s.pkl\n"%sim.save_name)
                m_qs = cf.load_instance("Mqs_%s.pkl"%sim.save_name)
                m_sq = cf.load_instance("Msq_%s.pkl"%sim.save_name)
                
            else:
                logger.info("There is no pkl file, so we need to analyse")
                m_qs, m_sq = sim.run_analysis()
            
            transport_coeff.append( prefactor * m_qs.transport_coeff(0, sim.tau_integration) )
            transport_coeff1.append( prefactor * m_sq.transport_coeff(0, sim.tau_integration) )
            
            # Getting all the correlation distributions in x and y
            correlation_dist.extend([m_qs.cor[0], m_qs.cor[1]])
            correlation_dist1.extend([m_sq.cor[0], m_sq.cor[1]])
            
    # =============================================================================
    #         test with other vdata
    # =============================================================================
    #        if (len(glob.glob("Mqs_vdata%s.pkl"%sim.save_name)) == 1):
    #    
    #            logger.info("\nThere is a pkl, we need to load Mqs_%s.pkl\n"%sim.save_name)
    #            m_qs_vdata = cf.load_instance("Mqs_vdata%s.pkl"%sim.save_name)
    #            m_sq_vdata = cf.load_instance("Msq_vdata%s.pkl"%sim.save_name)
    #            
    #        else:
    #            logger.info("There is no pkl file, so we need to analyse")
    #            m_qs_vdata, m_sq_vdata = sim.run_analysis_vdata("vdata_low.dat", volume_system, volume_bulk)
    #            
    #            
    #        fig,ax = m_qs_vdata.plot_individual(fig, ax, dim = 0, norm = False)
    #        ax.lines[-1].set_label(r"$vdata_x$")
    #        fig,ax = m_qs_vdata.plot_individual(fig, ax, dim = 1, norm = False)
    #        ax.lines[-1].set_label(r"$vdata_y$")
    #        
    #        
    #        fig1,ax1 = m_sq_vdata.plot_individual(fig1, ax1, dim = 0, norm = False)
    #        ax1.lines[-1].set_label(r"$vdata_x$")
    #        fig1,ax1 = m_sq_vdata.plot_individual(fig1, ax1, dim = 1, norm = False)
    #        ax1.lines[-1].set_label(r"$vdata_y$")
    #        
    #        transport_coeff.append( prefactor * m_qs_vdata.transport_coeff(0, sim.tau_integration) )
    #        transport_coeff1.append( prefactor * m_sq_vdata.transport_coeff(0, sim.tau_integration) )
        os.chdir(cwd)
    
    # =============================================================================
    # Creating the correlation bundles
    # =============================================================================
    correlation_dist = prefactor * np.array(correlation_dist)
    m_qs_bundle = fc.bundle_correlation(correlation_dist, m_qs.times, m_qs.flux1.name, m_qs.flux2.name)
    
    correlation_dist1 = prefactor * np.array(correlation_dist1)
    m_sq_bundle = fc.bundle_correlation(correlation_dist1, m_sq.times, m_sq.flux1.name, m_sq.flux2.name)
    
    # =============================================================================
    # Ploting all the correlations both in x and y
    # =============================================================================
    
    # Preparing the plots
    plt.close('all')
    cf.set_plot_appearance()
    
    fig, ax = plt.subplots()
    fig, ax = m_qs_bundle.plot_bundle(fig, ax)
    ax.axvline(x = sim.tau_integration,ls='-.',c='black')
    ax.set_xscale('log')
    xmin,xmax = ax.get_xlim()
    ax.set_xlim(xmin, sim.max_tau)
    ax.set_xlabel(r'$t^*$')
    fig.tight_layout()
    fig.savefig("%s/c_qs_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))
    
    fig, ax = plt.subplots()
    fig, ax = m_sq_bundle.plot_bundle(fig, ax)
    ax.axvline(x = sim.tau_integration,ls='-.',c='black')
    ax.set_xscale('log')
    xmin,xmax = ax.get_xlim()
    ax.set_xlim(xmin, sim.max_tau)
    ax.set_xlabel(r'$t^*$')
    fig.tight_layout()
    fig.savefig("%s/c_sq_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))
    
    # =============================================================================
    # Ploting the average correlation distribution for each coefficient
    # =============================================================================
    fig, ax = plt.subplots()
    fig, ax = m_qs_bundle.plot(fig, ax)
    ax.axvline(x = sim.tau_integration,ls='-.',c='black')
    ax.set_xscale('log')
    xmin,xmax = ax.get_xlim()
    ax.set_xlim(xmin, sim.max_tau)
    ax.set_xlabel(r'$t^*$')
    fig.tight_layout()
    fig.savefig("%s/correlation_ave_qs_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))
    
    # Adding the other
    fig, ax = m_sq_bundle.plot(fig, ax)
    ax.lines[-2].set_label(r'$I_{sq}$')
    ax.lines[0].set_label(r'$I_{qs}$')
    ax.set_ylabel(None)
    ax.set_xlabel(r'$t^*$')
    ax.legend(loc = 'upper right')
    fig.savefig("%s/correlation_ave_comparison_log_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))
    ax.set_xscale("linear")
    fig.savefig("%s/correlation_ave_comparison_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))
    
    fig, ax = plt.subplots()
    fig, ax = m_sq_bundle.plot(fig, ax)
    ax.axvline(x = sim.tau_integration,ls='-.',c='black')
    ax.set_xscale('log')
    xmin,xmax = ax.get_xlim()
    ax.set_xlim(xmin, sim.max_tau)
    ax.set_xlabel(r'$t^*$')
    fig.tight_layout()
    fig.savefig("%s/correlation_ave_sq_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))
    
    # =============================================================================
    # # Plot M_qs vs tau for all
    # =============================================================================
    fig, ax = plot_coeff_vs_tau(m_sq, sim, correlation_dist)
    ax.set_ylabel(r'$|Gamma^{QS}$')
    ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')
    fig.tight_layout()
    plt.savefig("%s/Mqs_vs_tau_error_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))
    
    fig, ax = plot_coeff_vs_tau(m_qs, sim, correlation_dist1)
    ax.set_ylabel(r'$|Gamma^{SQ}$')
    ax.axhline(y = 0, xmin=0, xmax=1,ls='--',c='black')
    fig.tight_layout()
    fig.savefig("%s/Msq_vs_tau_error_%s_%s.pdf"%(plot_dir, sim.name, folder_pattern))





# =============================================================================
# Main
# =============================================================================

cwd = os.getcwd()
plot_dir = "plots/0.gk_transport"

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   

# =============================================================================
# Simulation type definitions
# =============================================================================

# High sampling frequency
high = SimulationGK("high", 500, 200, 1, "log.lammps", "fluxes_high.dat") 
# Low sampling frequency but averaging
low = SimulationGK("low", 2000, 150, 1, "log.lammps", "fluxes_low.dat") 
# Low sampling frequency no averaging
other = SimulationGK("other", 2000, 150, 1, "log.lammps", "fluxes_other.dat") 
# Low sampling frequency but averaging with static prefactor
stat = SimulationGK("stat", 2000, 150, 1, "log.lammps", "fluxes_stat_low.dat") 

# Define the type of 

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Script to compute the Transport coefficients from DO simulations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('pattern', metavar='pattern',help = 'file patter, eg. "r*"', type=str)
    parser.add_argument('-type', metavar='type',help = 'simulation type. eg. low, stat, high', type=str)


    args = parser.parse_args()
    main(args.pattern, args.type)




