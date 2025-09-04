#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:41:15 2024

@author: zackaryschnitzer
"""

import numpy as np
from numpy import linalg as linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as math
import scipy.signal as signal
import time
import copy

stTime = time.time()
c = 1
dx = 0.1
dt = 0.01

Lx = 25
Ly = 100

piL = math.pi
mat = np.zeros((Lx,Ly))
boundaries = np.zeros((Lx,Ly))
for i in range(1, Lx-1):
    for j in range(1, Ly-1):
        #if ((i-Lx/2)**2+(j-Ly/2)**2)<25:
        #    boundaries[i][j] = 1
        if (j<5):
            mat[i][j] = 10 * np.sin(piL / 5 * j)
        if j<20 and (i<11 or i>14):
            boundaries[i][j] = 1


#stores the previous two sheets for easy access. index 0 is 2 back, index 1 is 1 back
lastTwo = [copy.copy(mat), copy.copy(mat)]

def boundaryOrNo(boundary, value):
    if boundary == 1:
        return 0
    return value

def approxSolveSheet(previousSheet, previousprevious):
    result = np.zeros((Lx,Ly))
    for i in range(1, Lx-1):
        for j in range(1,Ly-1):
            previousCenter = boundaryOrNo(boundaries[i][j], previousprevious[i][j])
            center = boundaryOrNo(boundaries[i][j], previousSheet[i][j])
            jminus = boundaryOrNo(boundaries[i][j-1], previousSheet[i][j-1])
            jplus = boundaryOrNo(boundaries[i][j+1], previousSheet[i][j+1])
            iminus = boundaryOrNo(boundaries[i-1][j], previousSheet[i-1][j])
            iplus = boundaryOrNo(boundaries[i+1][j], previousSheet[i+1][j])
            
            #given from arranging finite difference equations
            result[i][j] = 2 * center - previousCenter + c**2 * dt**2 / dx**2 * (iminus + iplus + jminus + jplus - 4 * center)
    
    return result
    
#will store all of the data
waveLists = copy.copy(lastTwo)


itts = 0
#will itterate forwards in time and add the sheets to waveLists
while (itts<3000):    
    solution = approxSolveSheet(lastTwo[1], lastTwo[0])
    waveLists.append(solution)
    
    lastTwo[0], lastTwo[1] = copy.copy(lastTwo[1]), copy.copy(solution)
    waveLists.append(solution)
    itts += 1
    

waveLists = np.array(waveLists)

"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.array(range(Ly))
y = np.array(range(Lx))
X,Y = np.meshgrid(x,y)


#plot the first frame
plot = ax.plot_surface(X,Y, waveLists[0],cmap=plt.cm.viridis, vmin = -4, vmax = 4)
ax.axes.set_zlim3d(bottom=-4, top=4) 
ax.set_aspect('equal')

#call this method in FuncAnimation to get animation
def update(frame):
    ax.clear()
    plot = ax.plot_surface(X,Y, waveLists[frame],cmap=plt.cm.viridis, vmin = -4, vmax = 4)
    ax.axes.set_zlim3d(bottom=-4, top=4) 
    ax.set_aspect('equal')
    return plot

ani = animation.FuncAnimation(fig, update, interval=1, frames = np.linspace(0,itts+2,100, dtype = int))
plt.show()

"""


print(f"Time Elapsed: {round(time.time()-stTime, 4)} seconds")



#ax.plot_surface(X,Y,waveLists[0])
