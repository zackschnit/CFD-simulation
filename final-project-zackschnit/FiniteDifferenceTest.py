#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:41:15 2024

@author: zackaryschnitzer
"""

import numpy as np
from numpy import linalg as linalg
import matplotlib.pyplot as plt
import time as time

start_time = time.time()

def genMatrix(xmin, xmax, N, boundary):
    dx = (xmax-xmin)/N
    #create N x N matrix filled with zeroes to start. to be filled with finite-difference eqs
    FDmatrix = np.zeros((N,N))
    
    #what everything equals to start
    initialVector = np.zeros((N,1))
    initialVector[0][0] = boundary
    
    #first entry is just 1 meaning first grid pt is just equal to boundary condition
    FDmatrix[0][0] = 1
    
    #populate the matrix with the correct equations
    for i in range(1,N):
        FDmatrix[i][i-1] = -1
        FDmatrix[i][i] = 1 + dx
    
    solution = linalg.solve(FDmatrix, initialVector)
    return solution
    
sol = genMatrix(1,2,100,1)
#exp = np.exp(-np.array(range(100))/100)
#plt.plot(exp)
plt.plot(sol)

print(time.time()-start_time)