#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:04:24 2020

@author: simon
"""

import MDAnalysis as mda


import MDAnalysis.analysis.rdf as rdf
import matplotlib.pyplot as plt

u = mda.Universe("system.data", "dcd_nvt.dcd", format="LAMMPS") 

n_molecules = u.atoms.n_residues
molecules = u.atoms.residues


for ts in u.trajectory[:2]:
    for mol in molecules:
        centroid = mol.atoms.centroid()
        print(centroid)

    
#    # analyze frame
#    if take_this_frame == True:
#    with mda.Writer('frame.data') as W:
#            W.write(u.atoms)
#         break

#u= mda.Universe("1.cxyz",format='xyz')
#
#monomers=u.select_atoms("type 3")
#solutes= u.select_atoms("type 2")
#
#g_r=rdf.InterRDF(monomers,solutes)
#g_r.run()
#
#plt.plot(g_r.bins, g_r.rdf)
#plt.show()

    
#g_r = rdf.InterRDF()