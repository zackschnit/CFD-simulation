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


#   
#----------------------------------------------------------------
#
#   For N = 10000 this code takes 0.10 seconds to run
#   The old code takes 50 seconds to do the same task
#   This is way more efficient with next to no difference 
#   in accuracy at this level. It's a no brainer to use this
#
#----------------------------------------------------------------  
#   


start_time = time.time()

def approxSolveSystem(xmin, xmax, N, itts, guessVals):
    dx = (xmax-xmin)/N
    #create an array to store values on this itteration
    result = N*["empty"]
    result[0] = guessVals[0] #boundary condition
    
    for i in range(1,N):
        if result[i-1] != "empty":
            result[i] = (result[i-1] + dx * guessVals[i]**2) / (1 + 2 * dx * guessVals[i])
            
        else:
            result[i] = (guessVals[i-1] + dx * guessVals[i]**2) / (1 + 2 * dx * guessVals[i])
    
    itts -= 1
    
    if itts == 0:
        return result
    
    return approxSolveSystem(xmin, xmax, N, itts, result)
    
initialGuess = np.zeros(10000)
initialGuess[0] = 1
sol = approxSolveSystem(1,2,10000,10,initialGuess)

print(time.time() - start_time)

plt.plot(sol)
plt.legend()










