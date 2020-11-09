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


print(os.getcwd())


def run_correlation_analysis(folder, input_file, save = "True"):
    
    """
    ******************Very specific function*************************
    
    Here all the fluxes have to be defined and the correlations
    
    Computes the correlations using the general class correlation
    *****************************************************************
    
    Args:
        folder is where the input_file is located
        input_file is the name of the file containing the fluxes
        save True if you want to save the correlation instances as pkl
    """
    
    cwd = os.getcwd()
    os.chdir(folder)
    data = cf.read_data_file(input_file)
    
    
    data1 = data.values
    delta_t = lu.read_value_from("log.lammps", 'timestep')[1]
    times = (data1[:,0]-data1[0,0])*delta_t
    max_delta = int(len(times)*0.5) #Maximum delta of time to measure the correlation
    
    
    
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

    # c12 = fc.correlation(solute_flux,solvent_flux,max_delta)
    # c12.evaluate_ccf()

    # c21 = fc.correlation(solvent_flux,solute_flux,max_delta)
    # c21.evaluate_ccf()

    # c22 = fc.correlation(solvent_flux,solvent_flux,max_delta)
    # c22.evaluate_acf()
    
    
    # # TEST TO SEE WHY ACF is good but ccf is not
    # c11_ccf = fc.correlation(solute_flux, solute_flux, max_delta)
    # c11_ccf.evaluate_ccf()
    
    
    if save == "True":
        c11.save('c11')
        # c12.save('c12')
        # c21.save('c21')
        # c22.save('c22')
    
    
    os.chdir(cwd)
    
    return c11
    # return [c11, c12, c21, c22]


def load_correlations(folder):
    cwd = os.getcwd()
    os.chdir(folder)
    
    c11 = cf.load_instance("c11.pkl")
    # c22 = cf.load_instance("c22.pkl")
    # c12 = cf.load_instance("c12.pkl")
    # c21 = cf.load_instance("c21.pkl")
    
    os.chdir(cwd)
    
    # return [c11,c12,c21,c22]
    
    return c11
  

def remove_pkl(root):
    
    """
    Removes all the pkl inside the directory recursively
    TODO generalise and put in cf
    """
    files = glob.glob('%s/**/*.pkl'%root,  recursive = True)
    for file in files:
        os.remove(file)



'''****************************************************************************
# =============================================================================
# MAIN
# =============================================================================
****************************************************************************'''

cwd = os.getcwd()
plot_dir = "plots"

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
logger = cf.log(__file__, os.getcwd(),plot_dir)   

# =============================================================================
# Input parameters (SPECIFIC TO THE PROBLEM)
# =============================================================================

root = "."
directory_pattern = '[0-9]*'
input_file = 'vdata.dat'



#     # Testing for one folder


# cwd = os.getcwd()
# for folder in glob.glob(directory_pattern):
    
#     os.chdir(folder)
    
#     data = cf.read_data_file(input_file)
    
    
#     data1 = data.values
#     delta_t = lu.read_value_from("log.lammps", 'timestep')[1]
#     times = (data1[:,0]-data1[0,0])*delta_t
#     max_delta = int(len(times)*0.004) #Maximum delta of time to measure the correlation
    
    
    
#     n_solv = lu.read_value_from("log.lammps", "atoms in group gSolv")
#     n_solu = lu.read_value_from("log.lammps", "atoms in group gSolu")
    
    
#     j_s_components = ["v_vx_Solu", "v_vy_Solu", "v_vz_Solu"]
#     j_f_components = ["v_vx_Solv", "v_vy_Solv", "v_vz_Solv"]
    
#     j_s = n_solu * data[j_s_components].values
#     j_f = n_solv * data[j_f_components].values
    
    
#     solute_flux = fc.flux(j_s, times, "J_s")
#     solvent_flux = fc.flux(j_f, times, "J_f")
    

#     # Evaluation using my function
#     c1 = fc.correlation(solute_flux, solute_flux, max_delta)
#     c1.evaluate_acf()
#     print(c1.transport_coeff(0,10))


#     # Evaluation using my acf
#     c2 = fc.correlation(solute_flux, solute_flux, max_delta)
#     c2.evaluate_ccf()
#     print(c2.transport_coeff(0,10))
#     os.chdir(cwd)



    


# For testing ONLY
#remove_pkl(root)

# =============================================================================
#  Checking the directories (SPECIFIC TO THE PROBLEM)
# =============================================================================

prep_results = sr.check_n_analyse(root,directory_pattern)
prep_results.check_finished("vdata.dat")
# Checks inside each directory
prep_results.check_stat("c11.pkl") #Assuming that if the first correlation is not ready, need to compute all

finished_directories=prep_results.dir_fin
unfinished_correlation=prep_results.dir_stat

print("\nFinished directories %s\n"%finished_directories)
print("\nUnfinished correlations in %s\n"%unfinished_correlation) 


# Array with the total correlations
c11_array=[]
# c12_array=[]
# c21_array=[]
# c22_array=[]

# perform the analysis if there are no instances saved in the main folder
if (len(glob.glob('c1*'))<1):
    for folder in finished_directories:
        #When Correlation still need to be run 
        #TODO change this to just gathering the results
        print("Now in folder %s\n"%folder)
        
        if folder in unfinished_correlation[0]:
            print("Running the correlation analysis in %s\n"%folder)
            # c11, c12, c21, c22 = run_correlation_analysis(folder, input_file)
            c11 = run_correlation_analysis(folder, input_file)
        
        else:
            print("Loading the results in %s\n"%folder)
            # c11, c12, c21, c22 = load_correlations(folder)
            
            c11 = load_correlations(folder)
            
            
        # Appending only the average of x,y,z
        c11_array.append(c11.cor)
        # c12_array.append(c12.cor)
        # c21_array.append(c21.cor)
        # c22_array.append(c22.cor)
    
    print(np.shape(c11_array))
    print(np.shape(np.ravel(c11_array)))
    #Extracting the times from one of the correlation
    times = c11.times
    
    #Creating the bundle instances 
    c11 = fc.bundle_correlation(c11_array, times, "J_s", "J_s")
    # c12 = fc.bundle_correlation(c12_array, times, "J_s", "J_f")
    # c21 = fc.bundle_correlation(c21_array, times, "J_f", "J_s")
    # c22 = fc.bundle_correlation(c22_array, times, "J_f", "J_f")
    
    #TODO this could be done iterating over all the instances
    c11.save('c11')
    # c12.save('c12')
    # c21.save('c21')
    # c22.save('c22')

else:
    print ("loading c11\n")
    c11 = cf.load_instance("c11.pkl")
    # print ("loading c22\n")
    # c22 = cf.load_instance("c22.pkl")
    # print ("loading c12\n")
    # c12 = cf.load_instance("c12.pkl")
    # print ("loading c21\n")
    # c21 = cf.load_instance("c21.pkl")
    


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


# #For c12
# fig,ax=plt.subplots()

# fig,ax = c12.plot(fig,ax)
# ax.set_xscale('log')
# plt.tight_layout()
# plt.savefig("correlation12.pdf")



# #For c21
# fig,ax=plt.subplots()

# fig,ax = c21.plot(fig,ax)
# ax.set_xscale('log')
# plt.tight_layout()
# plt.savefig("correlation21.pdf")



# #For c22
# fig,ax = plt.subplots()

# fig,ax = c22.plot(fig,ax)
# ax.set_xscale('log')
# plt.tight_layout()
# plt.savefig("correlation22.pdf")


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
# logger.info("The c12 is %s\n" %(pref*c12.transport_coeff(0, xmax)))
# logger.info("The c21 is %s\n" %(pref*c21.transport_coeff(0, xmax)))
# logger.info("The c22 is %s\n" %(pref*c22.transport_coeff(0, xmax)))



# =============================================================================
# Checking the plotting methods
# =============================================================================

#For c11
fig,ax = plt.subplots()

ax = c11.plot_bundle(ax, dim = -1)
ax.set_xscale('log')
plt.tight_layout()
plt.savefig("%s/correlation11_bundle.pdf"%plot_dir)



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
plt.savefig("%s/c11_vs_tau.pdf"%plot_dir)

