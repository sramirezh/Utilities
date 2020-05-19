#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 12:46:06 2020

Script to get the reduced densities and temperature for the paper

Borgelt, P., Hoheisel, C. and Stell, G. (1990) ‘Exact molecular dynamics and 
kinetic theory results for thermal transport coefficients of the Lennard-Jones 
argon fluid in a wide region of states’, Physical Review A, 42(2), pp. 789–794.

@author: simon
"""

ekb = 119.8 #  Kelvin
sigma = 3.405 # Angstrom
mass = 39.95  # A.U


rho = 1.414  # g/ cm -3


# conversion factors
g2part = (1/mass)*6.022*10**23 # from grams to particle number
cm2A = 10**8

rho_red = rho * g2part / (cm2A**3) 


