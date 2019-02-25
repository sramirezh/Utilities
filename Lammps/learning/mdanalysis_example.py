#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 16:33:02 2019
This file shows how to load a Universe and analyise the g(r)
@author: simon
"""
import MDAnalysis as mda
import MDAnalysis.analysis.rdf as rdf
import matplotlib.pyplot as plt


u= mda.Universe("1.cxyz",format='xyz')

monomers=u.select_atoms("type 3")
solutes= u.select_atoms("type 2")

g_r=rdf.InterRDF(monomers,solutes)
g_r.run()

plt.plot(g_r.bins, g_r.rdf)
plt.show()
