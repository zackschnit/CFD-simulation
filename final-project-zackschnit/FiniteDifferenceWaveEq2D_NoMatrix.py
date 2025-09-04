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



start_time = time.time()
c = 1
dx = 0.001
dt = 0.0001

L = 50

piL = math.pi / L
mat = np.zeros((L,L))
for i in range(1, L-1):
    for j in range(1, L-1):
        mat[i][j] = 2 * np.sin(piL * 3 * i) * np.sin(piL * 2 * j) + 2 * np.sin(piL * 2 * i) * np.sin(piL * 3 * j)


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
    a = -4 * c * dt**2 - dx**2
    b = c * dt**2
    
    for i in range(1, L-1):
        for j in range(1,L-1):
            iminus = guessVals[j][i-1]
            iplus = guessVals[j][i+1]
            jminus = guessVals[j-1][i]
            jplus = guessVals[j+1][i]
            
            #given from arranging finite difference equations
            guessVals[j][i] = (b * (iminus + iplus + jminus + jplus) - dx**2 * (lastTwo[1][j][i] - 2 * lastTwo[0][j][i])) / a
   
    itts -= 1
    
    #keep itterating (recursively), using the result as the new guess values, for a set # of itts.
    if itts == 0:
        return guessVals
    
    return approxSolveSheet(guessVals, itts)
    
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

ani = animation.FuncAnimation(fig, update, interval=1, frames = np.linspace(0,itts+2,200, dtype = int))
plt.show()

ani.save("FiniteDifferenceWaveEq_2D_Diverging.gif")

print(f"Time Elapsed: {round(time.time()-start_time,4)} seconds")

#ax.plot_surface(X,Y,waveLists[0])
