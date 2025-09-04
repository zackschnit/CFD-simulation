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

start_time = time.time()

c = 10
dx = 0.01
dt = 0.001


#2 back initial conditions (0 -> n = -1, 1 -> n = 0)
initial0 = -(4*signal.square(np.array(range(500))*math.pi/200 + 1) - 4)
initial1 = -(4*signal.square(np.array(range(500))*math.pi/200 + 1) - 4)
L = len(initial0)

#create N x N matrix filled with zeroes to start. to be filled with finite-difference eqs
FDmatrix = np.zeros((L,L))

#first entry is just 1 meaning first grid pt is just equal to boundary condition
FDmatrix[0][0] = 1
FDmatrix[L-1][L-1] = 1

#populate the matrix with the correct equations
for i in range(1,L-1):
     FDmatrix[i][i-1] = -c * (dt**2)
     FDmatrix[i][i] = dx**2 + 2 * c * (dt**2)
     FDmatrix[i][i+1] = -c * (dt**2)

#the matrix stays the same so we only need to invert once
FDinv = linalg.inv(FDmatrix)



waveLists = [initial0, initial1]


itts = 0
while (itts<3000):    
    #what everything equals to start including boundary conditions of 0
    colVector = np.zeros((L,1))
    colVector[0][0] = 0
    colVector[L-1][0] = 0
    
    #populate the col vector with the correct equations
    for i in range(1,L-1):
        colVector[i][0] = (dx**2) * (2 * initial1[i] - initial0[i])
        
    
    solution = np.rot90(np.matmul(FDinv, colVector))[0]
    waveLists.append(solution)
    
    initial0, initial1 = initial1, solution
    itts += 1
    

waveLists = np.array(waveLists)

"""
fig, ax = plt.subplots()
ax.set_ylim([-10,10])

line, = ax.plot(waveLists[0])

def update(frame):
    line.set_ydata(waveLists[frame])

ani = animation.FuncAnimation(fig, update, interval=1, frames = np.linspace(0,itts,200, dtype = int))
plt.show()
"""

print(f"Time Elapsed: {round(time.time()-start_time,4)} seconds")

