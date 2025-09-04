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
dx = 0.01
dt = 0.001

L = 50

piL = math.pi / L
mat = np.zeros((L,L))
for i in range(1, L-1):
    for j in range(1, L-1):
        mat[i][j] = 2 * np.sin(piL * 1 * i) * np.sin(piL * 2 * j) + 1 * np.sin(piL * 2 * i) * np.sin(piL * 2 * j)


#stores the previous two sheets for easy access. index 0 is 2 back, index 1 is 1 back
lastTwo = [copy.copy(mat), copy.copy(mat)]




#this matrix stays the same for the wave equation as long as the box is square
def genMatrix():
    #create N x N matrix filled with zeroes to start. to be filled with finite-difference eqs
    fdMat = np.zeros((L**2,L**2))
    a = -3 * c * dt**2 - dx**2
    b = 0.5 * c * dt**2
    d = 0.25 * c * dt**2
    
    #populate the matrix
    cells = 0
    while cells < L:
        #the top left cell of the current along the diagonal
        startIndex = L * cells
        
        #case where it is identity matrix:       [I]
        if cells == 0 or cells == L-1:
            for i in range(L):
                fdMat[startIndex + i][startIndex + i] = 1
        
        #case where it is BAB:                   [BAB]
        else:
            #deal with the two B matricies first:       [[B]A[B]] (NOTE I CHANGED THE KERNEL TO THE BETTER ONE WHICH IS WHY THERE IS A d TERM HERE)
            for i in range(1,L-1):
                fdMat[startIndex + i][startIndex - L + i - 1] = d
                fdMat[startIndex + i][startIndex - L + i] = b
                fdMat[startIndex + i][startIndex - L + i + 1] = d
                
                fdMat[startIndex + i][startIndex + L + i - 1] = d
                fdMat[startIndex + i][startIndex + L + i] = b
                fdMat[startIndex + i][startIndex + L + i + 1] = d
                
            #deal with the central A matrix:            [B[A]B]
            fdMat[startIndex][startIndex] = 1
            fdMat[startIndex + L - 1][startIndex + L - 1] = 1
            for i in range(1, L-1):
                fdMat[startIndex + i][startIndex - 1 + i] = b
                fdMat[startIndex + i][startIndex + i] = a
                fdMat[startIndex + i][startIndex + 1 + i] = b
        
        
        cells += 1
        
    
    return fdMat



fdMat = genMatrix()
fdInv = linalg.inv(fdMat)

for i in range(len(fdInv)):
    for j in range(len(fdInv[0])):
        fdInv[i][j] = round(fdInv[i][j], 2)



def genColVec():
    #make empty array to start
    colVec = np.zeros((L**2,1))
    #essentially collapse the 2d array to a 1d array. Carry same info in different format
    for i in range(1, L-1):
        for j in range(1, L-1):
            colVec[i*L + j] = [dx**2 * (lastTwo[0][i][j] - 2 * lastTwo[1][i][j])]
    return colVec



def solveSheet():
    colVec = genColVec()
    
    solution = np.matmul(fdInv, colVec)
    
    #put the matrix in a more comfortable form
    solution = np.rot90(solution)[0]

    result = np.zeros((L,L))
    for i in range(L):
        for j in range(L):
            #pull out one row at a time and add it to the result array
            result[i][j] = solution[i*L + j]
    return result

    
#will store all of the sata
waveLists = copy.copy(lastTwo)


itts = 0
#will itterate forwards in time and add the sheets to waveLists
while (itts<1000):    
    solution = solveSheet()
    waveLists.append(solution)
    
    lastTwo[0], lastTwo[1] = copy.copy(lastTwo[1]), copy.copy(solution)
    waveLists.append(solution)
    itts += 1
    

waveLists = np.array(waveLists)

"""
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

ani = animation.FuncAnimation(fig, update, interval=1, frames = np.linspace(0,itts+2,100, dtype = int))
plt.show()
"""


print(f"Time Elapsed: {round(time.time()-stTime,4)} seconds")
#ax.plot_surface(X,Y,waveLists[0])
