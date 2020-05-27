#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 12:24:05 2020
Todd 1995 precalculations for the low density
@author: simon
"""

rho_wall = 1.2549
lz_wall = 1.5625
lx = 5.2479
n_wall =  18*3
lz =  26.126

rho_wall_computed = n_wall/(lx*lx*lz_wall) 

a = (4/rho_wall)**(1/3)

d100 = a
d111 = 3**0.5*a/3
d110 = 2**0.5*a/2



# Supposing that we have lz_wall as the lattice parameter, then


a2 = lz_wall

rho_2 = 4/a2**3



# Getting the miller indexes 

d_spacing = lx 