#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 18:44:17 2019
File created copying what is in Property analysis but just to generate the force distribution
@author: sr802
"""

import numpy as np
import sys
import os
from shlex import split
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../')) #This falls into Utilities path
import Lammps.core_functions as cf
import Lammps.General.Log_Analysis.Thermo_Analyser as ta



FProperties=cf.read_data_file("Fproperties_short.dat") #Solvent properties
SProperties=cf.read_data_file("Sproperties_short.dat") #Solute properties
AProperties=cf.read_data_file("properties_short.dat") #All properties


# TODO this funtion could be replaced with the one that I used in the colloid/polymer theoretical mobility
def PosConstant(x,y,Tol):
    """
    Finds where a Force distribution type function no longer changes under the given Tolerance, analysing
    the points after the minimum and the maximum of the function. The chosen
    point corresponds to the minimum in the set of points


    :param x:
    :param y:
    :param Tol: Tolerance
    :return: The index of the position of the last point to be taken into account
    """


    BinSize = x[1] - x[0]
    der1=np.gradient(y,BinSize)
    IndexMax=np.argmax(y)
    IndexMin=np.argmin(y)
    ExtremeIndex=max(IndexMax,IndexMin)
    ExtremeIndex = max(IndexMax, IndexMin)
    Index_c = np.where(np.abs(der1) < Tol)[0]  # Indexes where the function is approximately constant
    LastIndex = np.where(Index_c > ExtremeIndex)[0][0]
    index = Index_c[LastIndex]
    return index


def Force(Ns,Nf,SProperties,FProperties,AProperties):
    """
    Force calculation
    :param Cs: Solute Concentration in the Bulk
    :param Cf: Solvent Concentration in the Bulk
    :param SProperties: Solute Properties
    :param FProperties: Solvent Properties
    :param AProperties: Fluid Properties
    :return:
    The Force distribution only up to IntUp, after that it is assumed that the force is zero

    """
    IntUp = 10
    Fs = 1
    Ff = -Fs * Ns / Nf
    Indexes=np.where(SProperties[:, 1]<=IntUp)[0]
    n=len(Indexes)
    Force = np.zeros((n, 2))
    Force[:, 0] = FProperties[Indexes, 1]

    for i in range(n):
        if AProperties[i, 4] == 0:
            Force[i, 1] = 0

        else:
            Force[i, 1] = (Ff * FProperties[i, 4] + Fs * SProperties[i, 4]) / AProperties[i, 4]

    np.savetxt("Force.dat", Force)
    print("\n*********************Running Yawei Calculations*********************\n")
    print("The force on the Solutes is %f, on the Solvents %f" %(Fs,Ff))
    print("Created the force distribution File Force.dat ")

    return Force



# Obtaining the thermodynamic data
    
thermo_data = ta.thermo_analyser("log.lammps")
Ns = float(thermo_data.at['v_cBSolu','Average'])
Nf = float(thermo_data.at['v_cBSolv','Average'])

Force = Force(Ns,Nf,SProperties,FProperties,AProperties)

Zshift = 0


Cut_off=PosConstant(Force[:,0],Force[:,1],0.001)+1

Zpos =Force[:Cut_off+1, 0]+Zshift
MuF = np.transpose(Force[:Cut_off, 1])
print("The Force Cut-off is %f, this is where the region of applied forces finishes"%np.max(Zpos))
print("Creating the Files to iterate in Lammps")
np.savetxt("Zpos_iterate.dat", Zpos)
np.savetxt("Force_iterate.dat", MuF)