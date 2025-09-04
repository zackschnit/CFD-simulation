#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 21:36:50 2024

@author: zackaryschnitzer
"""
import numpy as np
import copy
import matplotlib.pyplot as plt

class Grid:
    
    def __init__(self, xsize, ysize, overrelaxation, dx, density):
        self.xsize = xsize
        self.ysize = ysize
        self.overrelaxation = overrelaxation
        self.dx = dx
        self.density = density
        #   u is horizontal velocity, v is vertical velocity on a staggered grid
        #   Staggered grid requires one more value to get both left and right of each cell
        #   for [i,j] in the colocated grid, [i,j] is the left / top line of the cell
        #   
        #   ___
        #   |• 
        #   
        self.u = np.zeros(((ysize + 1), (xsize+1)))
        self.v = np.zeros(((ysize + 1), (xsize+1)))
        
        #   shifting to colocated grid. Values are for center of each cell
        self.p = np.zeros((ysize,xsize))    #pressure
        self.smokeD = np.zeros((ysize,xsize))   #smoke density
        self.solid = np.zeros((ysize,xsize))   
        self.solid.fill(1)
        #is cell solid (0 -> yes), (1 -> no)
        #   Top, and bottom rows should always be solid. This is a tank essentially.
        #   Solid will be replaced with it's own object down the road to account for more complicated solids
    
    
    #change velocities in accordance to an acceleration field (just g for now)
    def integrateGravity(self, dt, gravity):
        for r in range(self.ysize-1):
            for c in range(self.xsize):
                #if the cell is solid, or if the cell below is solid don't deal with gravity
                if (self.solid[r][c] != 0 and self.solid[r+1][c] != 0):
                    self.v[r][c] += gravity * dt
                    
    def fixIncompressibility(self, numItts, dt):
        
        for i in range(numItts):
            
            #only care about the middle bits of the wind tunnel so skip over one and end one short
            for r in range(1, self.ysize - 1):
                for c in range(1, self.xsize - 1):
                    
                    if (self.solid[r][c] == 0):
                        continue
                    
                    solidLeft = self.solid[r][c-1]
                    solidRight = self.solid[r][c+1]
                    solidTop = self.solid[r-1][c]
                    solidBot = self.solid[r+1][c]
                    
                    solidSum = solidLeft + solidRight + solidTop + solidBot
                    
                    if solidSum == 0:
                        continue
                    
                    uLeft = self.u[r][c]
                    uRight = self.u[r][c+1]
                    vTop = self.v[r][c]
                    vBot = self.v[r+1][c]
                    
                    #alternating minus signs to record OUTWARD flow
                    divergence = -uLeft + uRight - vBot + vTop
                    
                    #   To make the divergence zero we add back a little bit to each side
                    #   Overrelaxation term makes convergence much quicker
                    correction = self.overrelaxation * (divergence / solidSum)
                    
                    
                    #the actual correction for velocities (left, right, top, bottom)
                    self.u[r][c] += correction * solidLeft
                    self.u[r][c+1] -= correction * solidRight
                    self.v[r][c] += correction * solidTop
                    self.v[r+1][c] -= correction * solidBot
                    
                    #accounting for changes in pressure
                    self.p[r][c] += correction * (self.density * self.dx / dt)
                    
                    


nothing = np.zeros((100,100))
for i in range(100):
    for j in range(100):
        nothing[i][j] = np.sin(((i**2+j**2)**0.5)/10)
        
g = Grid(100,100,1,0.1,1)
g.u = copy.copy(nothing)
g.v = copy.copy(nothing)

plt.imshow(g.p)
plt.pause(1)


g.fixIncompressibility(40,0.01)


    

        
        
        