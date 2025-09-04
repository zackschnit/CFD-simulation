#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:41:15 2024

@author: zackaryschnitzer
"""

import numpy as np
from numpy import linalg as linalg
import matplotlib.pyplot as plt
import time


start_time = time.time()

def solveSystem(xmin, xmax, N, guessVals):
    dx = (xmax-xmin)/N
    #create N x N matrix filled with zeroes to start. to be filled with finite-difference eqs (LHS)
    FDmatrix = np.zeros((N,N))
    
    #create column vector for (RHS) 
    colVector = np.zeros((N,1))
    colVector[0][0] = guessVals[0] #boundary condition

    #first entry is just 1 meaning first grid pt is just equal to boundary condition
    FDmatrix[0][0] = 1
    
    for i in range(1,N):
        #populate the matrix with the correct equations
        FDmatrix[i][i-1] = -1
        FDmatrix[i][i] = 1 + 2 * dx * guessVals[i]
        #populate LHS
        colVector[i][0] = dx * (guessVals[i])**2
    
    solution = linalg.solve(FDmatrix, colVector)
    return solution
    
initialGuess = np.zeros(100)
initialGuess[0] = 1
sol = solveSystem(1,2,100,initialGuess)
#plt.plot(sol)
for i in range(10):
    plt.plot(sol, label = f"itt #{i}")
    sol = solveSystem(1,2,100,sol)

print(f"Time elapsed {round(time.time() - start_time,4)} seconds")

plt.legend()









