#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 19:17:37 2019
https://stackoverflow.com/questions/22566284/matplotlib-how-to-plot-images-instead-of-points
https://stackoverflow.com/questions/31401812/matplotlib-rotate-image-file-by-x-degrees

The only thing I need it the cow file and the trajectory files.
1. with the cow I generate the rotated views
2. with the trajectories I generate the frames
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
        pygame.image.save(rot_img, "rotated_%s.png"%degree)
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
    artists = []
    for x0, y0 in zip(x, y):
        angle=np.random.randint(0,359)
        name="rotated_%s.png"%angle
        image=plt.imread(name)
        im = OffsetImage(image, zoom=zoom)
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
#    ax.update_datalim(np.column_stack([x, y]))
    return artists



def cow_scatter(cows,ax=None,zoom=1):
    artists = []
    for cow in cows:
        image=plt.imread("cows/rotated_%s.png"%cow.angle)
        im = OffsetImage(image, zoom=zoom)
        ab = AnnotationBbox(im, (cow.x, cow.y), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
    return artists
    
    


class cow(object):
    def __init__(self,identity,x,y):
        """
        Args:
            identity in the lammps file
            angle rotation angle
            omega angular velocity
            position
        """
        self.identity=identity
        self.angle=np.random.randint(0,359)
        self.omega=np.random.randint(-5,5)
        self.x=x
        self.y=y
        
    def update_position(self,new_x,new_y):
        self.x=new_x
        self.y=new_y
    
    def update_angle(self):
        angle=self.angle+self.omega
        if angle<0:
            angle=359-angle
        angle=angle%359
        self.angle=angle
        
        

l_box=2.6726124191242437e+01

files = glob.glob("Trajectories/*.cxyz")
files.sort(key=lambda f: int(list(filter(str.isdigit, f))))

#intialise cow array




for j,file_name in enumerate(files):
    name=file_name.split('/')[1]
    time=name.split('.')[0]
    data=cf.read_data_file(file_name).values


    #Initialise cow array
    if j==0:
        cows=[]
        n,m=np.shape(data)
        for i in range(n):
            cows.append(cow(data[i,0],data[i,2],data[i,3]))
    #Update cows
    else:
        for i in range(n):
            cows[i].update_position(data[i,2],data[i,3])
            cows[i].update_angle()
            
    #Setting up the figure
    fig, ax = plt.subplots()
    ax.set_xlim(0,l_box)
    ax.set_ylim(0,l_box)
    #ax.axis("off")
    
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(3)


    #imscatter(x, y, image_path, zoom=0.3, ax=ax)
    fig_name="frames/frame_%s.png"%time
    cow_scatter(cows,zoom=0.3,ax=ax)
    plt.savefig(fig_name, format='png', dpi=300)




