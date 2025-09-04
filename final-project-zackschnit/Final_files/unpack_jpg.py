#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 18:23:57 2024

@author: zackaryschnitzer
"""

from PIL import Image
import numpy as np

def read_file(fileName, nx, ny):
    #these are defaults which I want everything to be a multiple of
    hTunnel = ny
    wTunnel = nx
    
    
    im = Image.open(fileName)
    im.load()
    
    width, height = im.size
    
    #Essentially will let me crop image to be the right aspect ratio. And I only want whole number increments
    mult = np.minimum(int(width/wTunnel), int(height/hTunnel))
    
    #all here to define the slice from the middle
    wantedHeight = hTunnel * mult
    wantedWidth = wTunnel * mult
    
    diffH = (height - wantedHeight)/2
    diffW = (width - wantedWidth)/2
    
    #actually get the pixels from the image
    pixels = np.array(im)
    
    #cropping out designated middle segment (using int to both convert to int and act as np.floor())
    pixels = pixels[int(diffH):int(height - diffH),int(diffW):int(width - diffW)]
    
    #"black and white" (1 if white, 0 if any other color)
    bw_image = np.zeros((wantedHeight, wantedWidth))
    for i in range(len(pixels)):
        for j in range(len(pixels[0])):
            if (pixels[i][j] == 255).all():
                bw_image[i][j] = 1
                
    #now to compressing this image, we will make use of the multiplier from earlier
    compressed_bw = np.zeros((ny,nx))
    
    for i in range(ny):
        for j in range(nx):
            miniArr = bw_image[mult * i:mult * (i+1), mult * j:mult * (j+1)]
            
            #can just count the nonzero elements since anything nonzero is just one
            countOnes = np.count_nonzero(miniArr)
            
            #if more than half of the spots in this chunk are white, then so should the compressed pixel
            if countOnes > (mult**2) / 2:
                compressed_bw[i][j] = 1

    return compressed_bw

