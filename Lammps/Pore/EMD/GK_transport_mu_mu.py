#!/usr/bin/env python3
"""
Created on Wed Aug 14 14:55:09 2019
This scripts computes the correlation between fluxes
TODO optimise and simplify
@author: sr802
"""


import numpy as np
import os
import sys
import glob
from uncertainties import unumpy,ufloat
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.Pore.qsub.simulation_results as sr
import Lammps.lammps_utilities as lu
import flux_correlation as fc 




class AnalysisGKPore(lu.SimulationType):
    
    def __init__(self, name, max_tau, tau_integration, temp,
                 root, dir_pattern, input_file ):
        """
        Args:
            name: Identifier of the simulation
            max_tau: maximum length of the correlation analysis in LJ time
            tau_integration: upper time limit of the correlation integration
            temp: system temperature in Kelvin
            root: Where all the simulations folders are in
            dir_pattern: the pattern of the simulations files
            input_file: file containing the fluxes
            
            
        Other atributes:
            sampling_interval
        """
        self.name = name
        self.max_tau = max_tau
        self.tau_integration = tau_integration
        self.temp = temp
        self.root = root
        self.dir_pattern = dir_pattern
        self.input_file = input_file
        
        
    def check_terminated(self, target_file):
        """
        Checks if the simulations inside the folders terminated and also if they
        had been analysed before by checking if one of the pkl is inside

        Parameters
        ----------
        target_file : name of the file that is going to be searched on every folder
        to check if the simulation finished

        Returns
        -------
        None.

        """
        
        prep_results = sr.check_n_analyse(self.root, self.dir_pattern)
        prep_results.check_finished(sim.input_file)
        prep_results.check_stat(target_file)
        self.finished_directories = prep_results.dir_fin
        self.unfinished_correlation = prep_results.dir_stat
        



def run_correlation_analysis(folder, input_file, save = "True"):
    
    """
    ******************Very specific function*************************
    
    Here all the fluxes have to be defined and the correlations
    
    Computes the correlations using the general class correlation
    *****************************************************************
    
    Args:
        
        input_file is the name of the file containing the fluxes
        save True if you want to save the correlation instances as pkl
    """
    
    cwd = os.getcwd()
    os.chdir(folder)
    data = cf.read_data_file(input_file)
    
    
    data1 = data.values
    delta_t = lu.read_value_from("log.lammps", 'timestep')[1]
    times = (data1[:,0]-data1[0,0])*delta_t
    max_delta = int(len(times)*0.005) #Maximum delta of time to measure the correlation
    
    
    
# =============================================================================
#      Preprocessing the data to generate the fluxes
# =============================================================================
    # Definitions of the fluxess
    # TODO the indexes of the fluxes could be found with core functions
    
    n_solv = lu.read_value_from("log.lammps", "atoms in group gSolv")
    n_solu = lu.read_value_from("log.lammps", "atoms in group gSolu")
    
    j_s_components = ["v_vx_Solu", "v_vy_Solu", "v_vz_Solu"]
    j_f_components = ["v_vx_Solv", "v_vy_Solv", "v_vz_Solv"]
    
    j_s = n_solu * data[j_s_components].values
    j_f = n_solv * data[j_f_components].values
    
    solute_flux = fc.flux(j_s, times, "J_s")
    solvent_flux = fc.flux(j_f, times, "J_f")
    
    c11 = fc.correlation(solute_flux, solute_flux, max_delta)
    c11.evaluate_acf()

    c12 = fc.correlation(solute_flux,solvent_flux,max_delta)
    c12.evaluate()

    c21 = fc.correlation(solvent_flux,solute_flux,max_delta)
    c21.evaluate()

    c22 = fc.correlation(solvent_flux,solvent_flux,max_delta)
    c22.evaluate_acf()
    
    
    # # TEST TO SEE WHY ACF is good but ccf is not
    # c11_ccf = fc.correlation(solute_flux, solute_flux, max_delta)
    # c11_ccf.evaluate_ccf()
    
    
    if save == "True":
        c11.save('c11')
        c12.save('c12')
        c21.save('c21')
        c22.save('c22')
    
    
    os.chdir(cwd)
    
    # return [c11,c12]
    return [c11, c12, c21, c22]


def load_correlations(folder):
    cwd = os.getcwd()
    os.chdir(folder)
    
    c11 = cf.load_instance("c11.pkl")
    c22 = cf.load_instance("c22.pkl")
    c12 = cf.load_instance("c12.pkl")
    c21 = cf.load_instance("c21.pkl")
    
    os.chdir(cwd)
    
    return [c11,c12,c21,c22]
    
    # return [c11, c12]
  

def remove_pkl(root):
    
    """
    Removes all the pkl inside the directory recursively
    TODO generalise and put in cf
    """
    files = glob.glob('%s/**/*.pkl'%root,  recursive = True)
    for file in files:
        os.remove(file)




# =============================================================================
# MAIN
# =============================================================================


# Defining all the types of simulations

run1 = AnalysisGKPore("Run1", 100, 10, 1, '.', '[0-9]*', 'vdata.dat')

cwd = os.getcwd()
plot_dir = "plots"

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   

# Defining the type of simulation to analyse
sim = run1.copy
sim.plot_dir = plot_dir
sim.print_params(logger)

# Checking if all the simulations finished
sim.check_terminated('c11.pkl')
logger.info("\nFinished directories %s\n"%sim.finished_directories)
logger.info("\nUnfinished correlations in %s\n"%sim.unfinished_correlation) 

# Array with the total correlations
c11_array=[]
c12_array=[]
c21_array=[]
c22_array=[]

# perform the analysis if there are no instances saved in the main folder
if (len(glob.glob('c1*'))<1):
    for folder in sim.finished_directories:
        #When Correlation still need to be run 
        #TODO change this to just gathering the results
        logger.info("Now in folder %s\n"%folder)
        
        if folder in sim.unfinished_correlation[0]:
            logger.info("Running the correlation analysis in %s\n"%folder)
            c11, c12, c21, c22 = run_correlation_analysis(folder, sim.input_file)
            # c11, c12 = run_correlation_analysis(folder, input_file)
        
        else:
            logger.info("Loading the results in %s\n"%folder)
            c11, c12, c21, c22 = load_correlations(folder)
            
            # c11, c12 = load_correlations(folder)
            # 
            
        # Appending only the average of x,y,z
        c11_array.append(c11.cor)
        c12_array.append(c12.cor)
        c21_array.append(c21.cor)
        c22_array.append(c22.cor)
    
    logger.info(np.shape(c12_array))
    logger.info(np.shape(np.ravel(c12_array)))
    
    #Extracting the times from one of the correlation
    times = c12.times
    
    #Creating the bundle instances 
    c11 = fc.bundle_correlation(c11_array, times, "J_s", "J_s")
    c12 = fc.bundle_correlation(c12_array, times, "J_s", "J_f")
    c21 = fc.bundle_correlation(c21_array, times, "J_f", "J_s")
    c22 = fc.bundle_correlation(c22_array, times, "J_f", "J_f")
    
    #TODO this could be done iterating over all the instances
    c11.save('c11')
    c12.save('c12')
    c21.save('c21')
    c22.save('c22')

else:
    logger.info("loading c11\n")
    c11 = cf.load_instance("c11.pkl")
    logger.info("loading c22\n")
    c22 = cf.load_instance("c22.pkl")
    logger.info("loading c12\n")
    c12 = cf.load_instance("c12.pkl")
    logger.info("loading c21\n")
    c21 = cf.load_instance("c21.pkl")
    


# =============================================================================
# Ploting all the correlations
# =============================================================================


plt.close('all')
cf.set_plot_appearance()

#For c11
fig,ax = plt.subplots()

fig,ax = c11.plot(fig,ax)
ax.set_xscale('log')
plt.tight_layout()
plt.savefig("%s/correlation11.pdf"%plot_dir)


#For c12
fig,ax=plt.subplots()

fig,ax = c12.plot(fig,ax)
ax.set_xscale('log')
plt.tight_layout()
plt.savefig("correlation12.pdf")



#For c21
fig,ax=plt.subplots()

fig,ax = c21.plot(fig,ax)
ax.set_xscale('log')
plt.tight_layout()
plt.savefig("correlation21.pdf")



#For c22
fig,ax = plt.subplots()

fig,ax = c22.plot(fig,ax)
ax.set_xscale('log')
plt.tight_layout()
plt.savefig("correlation22.pdf")


# =============================================================================
# Plot the cross coefficients
# =============================================================================

# xmax = 2000 # Integration limitc
# fig,ax = plt.subplots()

# fig,ax = c12.plot(fig,ax,ax_label=False)
# ax.lines[-1].set_label(r'$\langle (%s(t))(%s(0)) \rangle$'%(c12.flux1_name,c12.flux2_name)) #modifying the label of the last created line
# fig,ax = c21.plot(fig,ax,ax_label=False)
# ax.lines[-1].set_label(r'$\langle (%s(t))(%s(0)) \rangle$'%(c21.flux1_name,c21.flux2_name)) #modifying the label of the last created line
# ax.set_xscale('log')
# ax.axhline(y=0, xmin=0, xmax=1,ls=':',c='black')
# ax.axvline(x = xmax, c='black')
# plt.legend(loc='upper right')
# plt.tight_layout()
# plt.savefig("crossed.pdf")



# =============================================================================
# Performing the integration
# =============================================================================

V = 20**3 # Simulation box volume
T = 1.0 # Temperature
pref = 1/(V*T)
xmax = 10
#Todo, this could be added to each integral


logger.info("The c11 is %s\n" %(pref * c11.transport_coeff(0, xmax)))
logger.info("The c12 is %s\n" %(pref*c12.transport_coeff(0, xmax)))
logger.info("The c21 is %s\n" %(pref*c21.transport_coeff(0, xmax)))
logger.info("The c22 is %s\n" %(pref*c22.transport_coeff(0, xmax)))



# =============================================================================
# Checking the plotting methods
# =============================================================================

#For c11
# fig,ax = plt.subplots()

# ax = c11.plot_bundle(ax, dim = -1)
# ax.set_xscale('log')
# plt.tight_layout()
# plt.savefig("%s/correlation11_bundle.pdf"%plot_dir)



# =============================================================================
# # Plot c11 vs tau
# =============================================================================
xmax_plot = 100
fig,ax = plt.subplots()

tau_array = np.linspace(c11.times[0], c11.times[-1], 100)
# Need to add the integration time
tau_integration = 10
tau_array = np.sort(np.append(tau_array, tau_integration)) 

c_interval = []

# Computes the integral per intervals such that there is no need to recompute 
# Parts of the integral
#TODO be sure that the integration is continuos beween the intervals
initial_t = 0
for t in tau_array:
    c_interval.append(pref * c11.transport_coeff_comp(initial_t, t, -1))
    initial_t = t
    

# Neef to rewrite the first term, very quick fix
c_interval[0] = ufloat(0,0)

c_tau = np.cumsum(c_interval)
c_error = np.array([el.s for el in c_tau])
c_array = np.array([el.n for el in c_tau])


# Getting the data explicitly for tau_integration
index = np.where(tau_array==tau_integration)[0][0]
    
ax.plot(tau_array,c_array)
ax.fill_between(tau_array, c_array - c_error, c_array + c_error, alpha=0.4)
ax.set_xlabel(r'$t$')
ax.axhline(y = c_array[index] , xmin=0, xmax=1,ls='--',c='black', label = r'$\eta^* = %2.4f$' %c_array[index])
ax.axvline(x = tau_integration,ls='-.',c='black')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlim(0, xmax_plot)
ax.set_ylim(0, ymax)
ax.set_ylabel(r'$c_{11}$')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
plt.tight_layout()
plt.savefig("%s/c11.pdf"%plot_dir)


# =============================================================================
# # Plot c12 vs tau
# =============================================================================
xmax_plot = 100
fig,ax = plt.subplots()

tau_array = np.linspace(c12.times[0], c12.times[-1], 100)
# Need to add the integration time
tau_integration = 10
tau_array = np.sort(np.append(tau_array, tau_integration)) 

c_interval = []

# Computes the integral per intervals such that there is no need to recompute 
# Parts of the integral
#TODO be sure that the integration is continuos beween the intervals
initial_t = 0
for t in tau_array:
    c_interval.append(pref * c12.transport_coeff_comp(initial_t, t, -1))
    initial_t = t
    

# Neef to rewrite the first term, very quick fix
c_interval[0] = ufloat(0,0)

c_tau = np.cumsum(c_interval)
c_error = np.array([el.s for el in c_tau])
c_array = np.array([el.n for el in c_tau])


# Getting the data explicitly for tau_integration
index = np.where(tau_array==tau_integration)[0][0]
    
ax.plot(tau_array,c_array)
ax.fill_between(tau_array, c_array - c_error, c_array + c_error, alpha=0.4)
ax.set_xlabel(r'$t$')
ax.axhline(y = c_array[index] , xmin=0, xmax=1,ls='--',c='black', label = r'$\eta^* = %2.4f$' %c_array[index])
ax.axvline(x = tau_integration,ls='-.',c='black')
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
ax.set_xlim(0, xmax_plot)
ax.set_ylim(0, ymax)
ax.set_ylabel(r'$c_{12}$')
ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
plt.tight_layout()
plt.savefig("%s/c12_vs_tau.pdf"%plot_dir)