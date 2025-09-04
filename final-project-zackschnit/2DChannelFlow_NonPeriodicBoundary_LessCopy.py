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

#length in x and y
Lx = 2
Ly = 2


#num of points in x and y and t
nx = 41
ny = 41
nt = 500
numitts = 50

#used to correspond points in an array to points in space
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X,Y = np.meshgrid(x,y)

#density and viscosity
rho = 1
nu = 0.1

#other parameters which lead to convergence
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
dt = 0.001

#empty arrays (initial conditions and boundary conditions for the edges only)
u_initial = np.zeros((ny,nx))
v_initial = np.zeros((ny,nx))
p_initial = np.zeros((ny,nx))
#b is a helper array for calculating p
b_initial = np.zeros((ny,nx))

#wind into and out of tunnel at same speed
u_initial[1:-1,0] = 1
u_initial[1:-1,len(u_initial[0])-1] = 1



def helperB(u, v, p):
    bNew = np.empty_like(p)
    
    for i in range(1, len(u[0])-1):
        for j in range(1, len(u)-1):
            uiplus = u[j][i+1]
            uiminus = u[j][i-1]
            ujplus = u[j+1][i]
            ujminus = u[j-1][i]
            viplus = v[j][i+1]
            viminus = v[j][i-1]
            vjplus = v[j+1][i]
            vjminus = v[j-1][i]
            
            #equation for helper term for pressure Poisson eq
            bNew[j][i] = (1 / dt) * ((uiplus - uiminus) / (2*dx) + (vjplus - vjminus) / (2*dy)) - ((uiplus - uiminus) / (2*dx))**2 - 2*((ujplus - ujminus) / (2*dy))*((viplus-viminus) / (2*dx)) - ((vjplus - vjminus) / (2*dy))**2 
    
    return bNew
            


def poissonPressure(u, v, p):
    pNew = np.empty_like(p)
    pNew = p.copy()
    b = helperB(u, v, p)
    #itterate to stabilize solution
    for itts in range(numitts):
        
        for i in range(1, len(u[0])-1):
            for j in range(1, len(u)-1):
                piplus = pNew[j][i+1]
                piminus = pNew[j][i-1]
                pjplus = pNew[j+1][i]
                pjminus = pNew[j-1][i]
                bCenter = b[j][i]
                
                #pressure Poisson eq in form for finite difference method
                pNew[j][i] = ((piplus + piminus)*dy**2 + (pjplus + pjminus)*dx**2) / (2*(dx**2 + dy**2)) - (rho * dx**2 * dy**2) / (2*(dx**2 + dy**2)) * bCenter
                
    return pNew


def flow(u, v, p):
    uList = [u.copy()]
    vList = [v.copy()]
    pList = [p.copy()]
    
    for itts in range(nt):
        un = u.copy()
        vn = v.copy()
        
        p = poissonPressure(u, v, p)
        
        #this puts p "one back". Essentially ensures that each frame is physically correct before calculating the next
        if itts == 0:
            pList[0] = p
            
        else:
            pList.append(p)
        
        for i in range(1, len(u[0])-1):
            for j in range(1, len(u)-1):
                uhere = u[j][i]
                uiplus = u[j][i+1]
                uiminus = u[j][i-1]
                ujplus = u[j+1][i]
                ujminus = u[j-1][i]
                
                vhere = v[j][i]
                viplus = v[j][i+1]
                viminus = v[j][i-1]
                vjplus = v[j+1][i]
                vjminus = v[j-1][i]
                
                phere = p[j][i]
                piplus = p[j][i+1]
                piminus = p[j][i-1]
                pjplus = p[j+1][i]
                pjminus = p[j-1][i]
                
                un[j][i] = uhere - uhere*(dt/dx)*(uhere - uiminus) - vhere*(dt/dy)*(uhere - ujminus) - (dt / (rho * 2 * dx))*(piplus - piminus) + nu * ((dt / dx**2)*(uiplus - 2*uhere + uiminus) + (dt / dy**2)*(ujplus - 2*uhere + ujminus))
                vn[j][i] = vhere - uhere*(dt/dx)*(vhere - viminus) - vhere*(dt/dy)*(vhere - vjminus) - (dt / (rho * 2 * dy))*(pjplus - pjminus) + nu * ((dt / dx**2)*(viplus - 2*vhere + viminus) + (dt / dy**2)*(vjplus - 2*vhere + vjminus))
                
        uList.append(un)
        vList.append(vn)
        
        u = un
        v = vn

    pList.append(poissonPressure(u, v, p))
    
    #convert and return them as arrays now because we are done with the appending
    return uList, vList, pList
                

#call the function
uSol, vSol, pSol = flow(u_initial, v_initial, p_initial)

uSol = np.array(uSol)
vSol = np.array(vSol)
pSol = np.array(pSol)

    


fig = plt.figure()

#not in 3d yet!
#ax = fig.add_subplot(111, projection='3d')




#plot the first frame
plot = plt.contourf(X, Y, pSol[0], alpha=0.5, cmap=plt.cm.viridis)
plt.colorbar()
plt.contour(X, Y, pSol[0], cmap=plt.cm.viridis)
#plt.quiver(X[::2, ::2], Y[::2, ::2], uSol[0, ::2, ::2], vSol[0, ::2, ::2])
plt.streamplot(X, Y, uSol[0], vSol[0])

#call this method in FuncAnimation to get animation
def update(frame):
    fig.clear()
    plot = plt.contourf(X, Y, pSol[frame], alpha=0.5, cmap=plt.cm.viridis)
    plt.colorbar()
    plt.contour(X, Y, pSol[frame], cmap=plt.cm.viridis)
    plt.quiver(X[::2, ::2], Y[::2, ::2], uSol[frame, ::2, ::2], vSol[frame, ::2, ::2])
    #plt.streamplot(X, Y, uSol[frame], vSol[frame])
    return plot

ani = animation.FuncAnimation(fig, update, interval=1, frames = len(uSol))
plt.show()
print(time.time()-stTime)



#ax.plot_surface(X,Y,waveLists[0])
