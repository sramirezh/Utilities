"""
This script reads a final configuration from Lammps and performs widom insertion method of two different species 
in two different reservoirs to compute the chemical potential.
"""

import numpy as np

Data=np.loadtxt("trajectory.xyz", skiprows=2)

