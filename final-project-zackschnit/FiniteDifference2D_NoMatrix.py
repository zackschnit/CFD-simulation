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


c = 1
dx = 0.01
dt = 0.001

L = 50

piL = math.pi / L
mat = np.zeros((L,L))
for i in range(1, L-1):
    for j in range(1, L-1):
        mat[i][j] = (i + j) / 25 - 5


#stores the previous two sheets for easy access. index 0 is 2 back, index 1 is 1 back
lastTwo = [copy.copy(mat), copy.copy(mat)]


def empty2dList(x, y):
    result = []
    for i in range(y):
        l = []
        for j in range(x):
            l.append("empty")
        result.append(l)
    
    return result


def approxSolveSheet(guessVals, itts):
    #this should be -a but doing that gives a wildly wrong solution
    a = 2 * c * dt**2 - dx
    b = c * dt**2
    
    result = empty2dList(L, L)
    
    for i in range(L):
        for j in range(L):
            if i == 0 or j == 0:
                 result[j][i] = guessVals[j][i]   
                 
            else: 
                iminus = result[j][i-1]
                jminus = result[j-1][i]

                
                #if there isn't a value on this itteration use the guess value instead
                if iminus == "empty":
                    iminus = guessVals[j][i-1]

                if jminus == "empty":
                    jminus = guessVals[j-1][i]

                
                #given from arranging finite difference equations
                result[j][i] = (b * (iminus + jminus) + dx * (lastTwo[1][j][i] - 2 * lastTwo[0][j][i])) / a
       
    itts -= 1
    
    #keep itterating (recursively), using the result as the new guess values, for a set # of itts.
    if itts == 0:
        return result
    
    return approxSolveSheet(result, itts)
    
#will store all of the data
waveLists = copy.copy(lastTwo)


itts = 0
#will itterate forwards in time and add the sheets to waveLists
while (itts<500):    
    solution = approxSolveSheet(lastTwo[1],10)
    waveLists.append(solution)
    
    lastTwo[0], lastTwo[1] = copy.copy(lastTwo[1]), copy.copy(solution)
    waveLists.append(solution)
    itts += 1
    

waveLists = np.array(waveLists)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.array(range(L))
y = np.array(range(L))
X,Y = np.meshgrid(x,y)

#plot the first frame
plot = ax.plot_surface(X,Y, waveLists[0],cmap=plt.cm.viridis, vmin = -4, vmax = 4)
ax.axes.set_zlim3d(bottom=-4, top=4) 

#call this method in FuncAnimation to get animation
def update(frame):
    ax.clear()
    plot = ax.plot_surface(X,Y, waveLists[frame],cmap=plt.cm.viridis, vmin = -4, vmax = 4)
    ax.axes.set_zlim3d(bottom=-4, top=4) 
    return plot

ani = animation.FuncAnimation(fig, update, interval=1, frames = itts+2)
plt.show()



#ax.plot_surface(X,Y,waveLists[0])
