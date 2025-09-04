#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 12:15:32 2024

@author: zackaryschnitzer
"""
import numpy as np
import scipy.optimize as sy


def windAndForce(u, v, p, obj, dx):
    windSpeed = []
    force = []
    for i in range(len(u)):
        mag = np.hypot(u[i], v[i])
        speed = mag[::,0].mean()
        windSpeed.append(speed)
        f, area = horizontalForceOnObject(p[i], obj, speed, dx)
        force.append(f)
        
    return windSpeed, force, area



def horizontalForceOnObject(p, obj, speed, dx):
    force = dx * (sum(p[15:25,34]) - sum(p[15:25,45]))
    area = 0
    
    for i in range(len(p)):
        objInRow = False
        for j in range(len(p[0])):
            if obj[i][j] == 0:
                objInRow = True
                
                if obj[i][j-1] == 1:
                    force += dx * p[i][j-1]
                
                if obj[i][j+1] == 1:
                    force -= dx * p[i][j+1]
        if objInRow:
            area += dx
    
    return force, area



def cleanUpLists(mainList, percent, *args):
    length = len(mainList)
        
    otherLists = []
    
    for lists in args:
        #make sure that it is actually a list type object first
        try: 
            len(lists)
        
        except:
            #do nothing but still let me know
            print("Error: tried using an object that wasn't a list")
        
        else:
            #only consider if it is the same length
            if len(lists) == length:
                otherLists.append(lists)
    
    cutoff = length - 1
    
    for i in range(length-1):
        if mainList[i] != 0:
            per = (abs(mainList[i] - mainList[i+1]) / mainList[i]) * 100
        
            #if the percent change is less than the cutoff we have found the place we will start
            if per < percent:
                cutoff = i
                break
            

        
    
    mainList = mainList[cutoff:length - 1]
    
    for i in range(len(otherLists)):
        otherLists[i] = otherLists[i][cutoff:length - 1]
    
    return cutoff, mainList, otherLists



def dragCalculations(windSpeed, force, rho, area):
    cutOff, windSpeed, others = cleanUpLists(windSpeed, 0.5, force)
    force = np.array(others[0])
    windSpeed = np.array(windSpeed)
    
    def forceFit(v, a, D, c):
        return (1/2) * rho * area * D * v**2 + a * v + c
    
    try:
        params, covar = sy.curve_fit(forceFit, windSpeed, force)
    except:
        return [-1, -1, -1]
    
    #params[0] is viscous drag coefficient
    #params[1] is pressure drag coefficient
    return params

        
    


