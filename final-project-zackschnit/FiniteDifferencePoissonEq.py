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

dx = 0.1

Lx = 50
Ly = 50

piL = math.pi
initial = np.zeros((Lx,Ly))
RHS = np.zeros((Lx,Ly))

RHS[int(Lx/4)][int(Ly/4)] = 1
RHS[int(3*Lx/4)][int(3*Ly/4)] = -1





def approxSolveSheet(previousSheet, itts):
    if itts == 0:
        return previousSheet
    
    result = np.zeros((Lx,Ly))
    for i in range(1, Lx-1):
        for j in range(1,Ly-1):

            RHSval = RHS[i][j]
            jminus = previousSheet[i][j-1]
            jplus = previousSheet[i][j+1]
            iminus = previousSheet[i-1][j]
            iplus = previousSheet[i+1][j]
            
            #given from arranging finite difference equations
            result[i][j] = (iminus + iplus + jminus + jplus) / 4 - RHSval / (4 * dx**2)
    
    return approxSolveSheet(result, itts-1)
    



solution = approxSolveSheet(initial, 100)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.array(range(Ly))
y = np.array(range(Lx))
X,Y = np.meshgrid(x,y)

#plot the first frame
plot = ax.plot_surface(X,Y, solution,cmap=plt.cm.viridis, vmin = -4, vmax = 4)
ax.axes.set_zlim3d(bottom=-4, top=4) 



plt.show()
print(time.time()-stTime)



#ax.plot_surface(X,Y,waveLists[0])
