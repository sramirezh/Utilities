#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 19:17:37 2019
https://stackoverflow.com/questions/22566284/matplotlib-how-to-plot-images-instead-of-points
https://stackoverflow.com/questions/31401812/matplotlib-rotate-image-file-by-x-degrees
@author: sr802
"""
import matplotlib.pyplot as plt
from scipy import ndimage
import sys
import os
import glob
import pandas as pd
name="cow.png"
image=plt.imread(name)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib import animation

sys.path.append(os.path.join(os.path.dirname(__file__), '../../')) #This falls into Utilities path
import Lammps.core_functions as cf

import warnings
warnings.filterwarnings("ignore")


def generate_rotated_cows():
    import pygame
    
    pygame.init()
    screen = pygame.display.set_mode([400, 400])
    pygame.display.set_caption('Rotating image example')
    clock = pygame.time.Clock()
    
    img = pygame.image.load('cow.png').convert_alpha()
    
    img_rect = img.get_rect(center = screen.get_rect().center)
    degree = 0
    
    while degree < 360:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                done = True
    
        # rotate image
        rot_img = pygame.transform.rotate(img, degree)
        img_rect = rot_img.get_rect(center = img_rect.center)
        imagen=pygame.image.save(rot_img, "rotated_%s.png"%degree)
        # copy image to screen
        screen.fill((0, 0, 0))
        screen.blit(rot_img, img_rect)
        pygame.display.flip()
    
        clock.tick(60)
        degree += 1
    
    pygame.quit()

def imscatter(x, y, image, ax=None, zoom=1):
    """
    Puts the image in the desired coordinates
    
    """
    if ax is None:
        ax = plt.gca()
    try:
        image = plt.imread("rotated_339.png")
    except TypeError:
        # Likely already an array...
        pass
    im = OffsetImage(image, zoom=zoom)
    degree=0
    artists = []
    for x0, y0 in zip(x, y):
        print x0,y0
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
#    ax.update_datalim(np.column_stack([x, y]))
    return artists

l_box=2.6726124191242437e+01

files = glob.glob("*.cxyz")
files.sort(key=lambda f: int(filter(str.isdigit, f)))

file_name="20.cxyz"
data=cf.read_data_file(file_name).values

X=data[:,2]
Y=data[:,3]



x = X
y = Y
image_path = name
fig, ax = plt.subplots()
ax.set_xlim(0,l_box)
ax.set_ylim(0,l_box)
imscatter(x, y, image_path, zoom=0.3, ax=ax)
plt.savefig('destination_path.png', format='png', dpi=1000)
plt.show()



