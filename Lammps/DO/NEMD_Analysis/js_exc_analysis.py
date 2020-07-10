#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 17:54:21 2020
This script is to analyse in more detail what happens with the excess solute flux
In pressure driven simulations
@author: simon
"""

import os
import sys
Utilities_path=os.path.join(os.path.dirname(__file__), '../../../')
sys.path.append(Utilities_path) #This falls into Utilities path
import Lammps.core_functions as cf
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import Lammps.Pore.qsub.simulation_results as sr
from uncertainties import ufloat,unumpy
import glob
import re
import Lammps.lammps_utilities as lu
import Lammps.DO.EMD.density_analysis as da



# Reading necesary data
gradP = 0.00025

df = cf.read_data_file('vdata.dat')
folder = './'
box_volume, box_limits = lu.read_box_limits('%s/log.lammps'%folder)
bulk_heigth  = lu.read_region_height('rBulk','%s/in.geom'%folder)
#This is a region defined as the bulk+interface
sys_heigth  = lu.read_region_height('rSystem','%s/in.geom'%folder) 

vol_bulk = (box_limits[1][0]-box_limits[0][0])*(box_limits[1][1]-box_limits[0][1])*bulk_heigth[0]
vol_sys = (box_limits[1][0]-box_limits[0][0])*(box_limits[1][1]-box_limits[0][1])*(sys_heigth[0])


# =============================================================================
# Nomenclature for the averages
# # creating timewise post-processed variables (tp)
# # f on the fly
# # p post processed, only taking into account the averages
# Nomenclature for the excess flow
# the new way of defining things by 3 letters, ooo, means JS is for omega, prefactor omega
# and Q omega
# =============================================================================

# First just Js
df['js_omega_tp'] = (df['v_cSolu']/vol_sys) * df['v_vx_Solu']
js_omega_tp = df['js_omega_tp'].mean()
js_omega_f = df['v_Js'].mean()
js_omega_p = (df['v_cSolu'].mean()/vol_sys) * df['v_vx_Solu'].mean()

# Now the excess Js - c_s Q

df['js_exc_omega_tp'] = (df['v_cSolu']/vol_sys) * (df['v_vx_Solu']- df['v_vx_Sol'])
js_exc_omega_tp = df['js_exc_omega_tp'].mean()
js_exc_omega_f = df['v_Js_exc'].mean()
js_exc_omega_p = (df['v_cSolu'].mean()/vol_sys) *( df['v_vx_Solu'].mean()-df['v_vx_Sol'].mean())


# Now the excess Js - c_s^B Q

df['js_exc_obo_tp'] = (df['v_cSolu']/vol_sys) * df['v_vx_Solu']- (df['v_cBSolu']/vol_bulk)*df['v_vx_Sol']
js_exc_obo_tp = df['js_exc_obo_tp'].mean()

# Now the excess Js - phi*c Q
df['js_exc_hy_tp'] = (df['v_cSolu']/vol_sys) * df['v_vx_Solu']- (df['v_cBSolu']/(df['v_cBSolu']+df['v_cBSolv']))*((df['v_cSolu']+df['v_cSolv'])/vol_sys)*df['v_vx_Sol']
js_exc_hy_tp = df['js_exc_hy_tp'].mean()
js_exc_hy_p = (df['v_cSolu'].mean()/vol_sys) * df['v_vx_Solu'].mean()- (df['v_cBSolu'].mean()/(df['v_cBSolu'].mean()+df['v_cBSolv'].mean()))*((df['v_cSolu'].mean()+df['v_cSolv'].mean())/vol_sys)*df['v_vx_Sol'].mean()


# Now Js - (cs/cf) Jf
df['js_exc_jf_tp'] = (df['v_cSolu']/vol_sys) * (df['v_vx_Solu']- df['v_vx_Solv'])
js_exc_jf_tp = df['js_exc_jf_tp'].mean()


# Now Js - (cs^B/cf^B) Jf
df['js_exc_jf_obo_tp'] = (df['v_cSolu']/vol_sys) * df['v_vx_Solu']- (df['v_cBSolu']/df['v_cBSolv'])*(df['v_cSolv']/vol_sys)*df['v_vx_Solv']
js_exc_jf_obo_tp = df['js_exc_jf_obo_tp'].mean()


# Now the excess Js^B - c_s^B Q^B
df['js_exc_bbb_tp'] = (df['v_cBSolu']/vol_bulk) * df['v_vxB_Solu']- (df['v_cBSolu']/vol_bulk)*df['v_vxB_Sol']
js_exc_bbb_tp = df['js_exc_bbb_tp'].mean()


# =============================================================================
# Other comparisons
# =============================================================================


#  cQ (Omega) =? cQ (Bulk)
df['c_Q_o'] =  ((df['v_cSolu'].mean()+df['v_cSolv'].mean())/vol_sys)*df['v_vx_Sol'].mean()
df['c_Q_b'] = ((df['v_cBSolu'].mean()+df['v_cBSolv'].mean())/vol_bulk)*df['v_vxB_Sol'].mean()

cq_o_tp = df['c_Q_o'].mean()
cq_b_tp = df['c_Q_b'].mean()