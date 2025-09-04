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
import matplotlib.colors as mcolors
import math as math
import scipy.signal as signal
import time
import copy
import scipy.optimize as sy

stTime = time.time()
#
#
# -------------------------------------------------------------------------------------
#
#                           A WHOLE BUNCH OF PARAMETERS
#
# -------------------------------------------------------------------------------------
#
#
#length in x and y
Lx = 4
Ly = 2

#num of points in x and y and t
nx = int(20.5 * Lx)
ny = int(20.5 * Ly)

nt = 100
numitts = 50

#used to correspond points in an array to points in space
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X,Y = np.meshgrid(x,y)

#density and viscosity
rho = 1
nu = 0.00002

#other parameters which lead to convergence
dx = Lx / nx
dy = Ly / ny
dt = 0.001

#empty arrays (initial conditions and boundary conditions for the edges only)
u_initial = np.zeros((ny,nx))
v_initial = np.zeros((ny,nx))
p_initial = np.zeros((ny,nx))

#force which is 1 everywhere (may have to change when objects start being added to the channel)
F = 1*np.ones((ny,nx))

#
#
# -------------------------------------------------------------------------------------
#
#                DEFINING A FUNCTION WHICH CAN CHANGE THESE PARAMETERS
#           (I know that this is essentially a class, and I'm writing the constructor,
#      but I've already written the entire code and can't be bothered to add .self to everything)
#
# -------------------------------------------------------------------------------------
#
#

def changeParameters(**kwargs):
    #again I know this is bad code but I will mess it all up if I have to add .self in hundreds of spots
    global Lx, Ly, nt, rho, nu, nx, ny, dx, dy, dt, x, y, X, Y, u_initial, v_initial, p_initial, F

    if "xLength" in kwargs.keys():
        Lx = kwargs["xLength"]
        
    if "yLength" in kwargs.keys(): 
        Ly = kwargs["yLength"]
    
    if "timeItterations" in kwargs.keys():
        nt = kwargs["timeItterations"]
    
    if "density" in kwargs.keys():
        rho = kwargs["density"]
    
    if "viscocity" in kwargs.keys():
        nu = kwargs["viscosity"]
    
    if "force" in kwargs.keys():
        F = kwargs["force"]*np.ones((ny,nx))
        
    if "timeStep" in kwargs.keys():
        dt = kwargs["timeStep"]
    
    #num of points in x and y and t
    nx = int(20.5 * Lx)
    ny = int(20.5 * Ly)
    
    #other parameters which lead to convergence
    dx = Lx / nx
    dy = Ly / ny
    
    #used to correspond points in an array to points in space
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    X,Y = np.meshgrid(x,y)
    
    #empty arrays (initial conditions and boundary conditions for the edges only)
    u_initial = np.zeros((ny,nx))
    v_initial = np.zeros((ny,nx))
    p_initial = np.zeros((ny,nx))
    
        

#if there is an object, make it zero, else it is one
obj = np.ones((ny,nx))

#little square in the middle
obj[15:25, 35:45] = 0

def horizontalForceOnSquare(p, speed):
    force = dx * (sum(p[15:25,34]) - sum(p[15:25,45]))
    
    return force


#for a circle
#for i in range(len(obj[0])):
    #for j in range(len(obj)):
        #if ((i-40)**2+(j-20)**2)<25:
            #obj[j][i] = 0



#
#
# -------------------------------------------------------------------------------------
#
#
#                           SOLVING THE DIFFERENTIAL EQUATIONS
#
#
# -------------------------------------------------------------------------------------
#
#

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
        #will set velocity to zero if there is an object in the way
        #This is the implementation of no slip boundary condition
        u *= obj
        v *= obj
        
        un = u.copy()
        vn = v.copy()
        
        p = poissonPressure(u, v, p)
        
        #this puts p "one back". Essentially ensures that each frame is physically correct before calculating the next
        if itts == 0:
            pList[0] = p.copy()
            
        else:
            pList.append(p.copy())
        
        #will account for periodic boundary with if statements, so span the whole area for i
        for i in range(0, len(u[0])):
            for j in range(1, len(u)-1):
                uhere = u[j][i]
                ujplus = u[j+1][i]
                ujminus = u[j-1][i]
                
                vhere = v[j][i]
                vjplus = v[j+1][i]
                vjminus = v[j-1][i]
                
                pjplus = p[j+1][i]
                pjminus = p[j-1][i]
                
                #periodic boundary condition
                if i == 0:
                    #stuff that changes
                    uiminus = u[j][len(u[0])-1]
                    viminus = v[j][len(u[0])-1]
                    piminus = p[j][len(u[0])-1]
                    
                    #normal
                    uiplus = u[j][i+1]
                    viplus = v[j][i+1]
                    piplus = p[j][i+1]
                
                elif i == len(u[0])-1:
                    #stuff that changes
                    uiplus = u[j][0]
                    viplus = v[j][0]
                    piplus = p[j][0]
                    
                    #normal
                    uiminus = u[j][i-1]
                    viminus = v[j][i-1]
                    piminus = p[j][i-1]
                    
                else:
                    uiplus = u[j][i+1]
                    uiminus = u[j][i-1]
                    viplus = v[j][i+1]
                    viminus = v[j][i-1]
                    piplus = p[j][i+1]
                    piminus = p[j][i-1]
                    
                
                fhere = F[j][i]
                
                #Navier-Stokes equation in form for finite difference
                un[j][i] = uhere - uhere*(dt/dx)*(uhere - uiminus) - vhere*(dt/dy)*(uhere - ujminus) - (dt / (rho * 2 * dx))*(piplus - piminus) + nu * ((dt / dx**2)*(uiplus - 2*uhere + uiminus) + (dt / dy**2)*(ujplus - 2*uhere + ujminus)) + dt*fhere
                vn[j][i] = vhere - uhere*(dt/dx)*(vhere - viminus) - vhere*(dt/dy)*(vhere - vjminus) - (dt / (rho * 2 * dy))*(pjplus - pjminus) + nu * ((dt / dx**2)*(viplus - 2*vhere + viminus) + (dt / dy**2)*(vjplus - 2*vhere + vjminus))
                
        uList.append(un)
        vList.append(vn)
        print(f"Itteration #: {itts}")
        u = un.copy()
        v = vn.copy()

    pList.append(poissonPressure(u, v, p))
    
    #convert and return them as arrays now because we are done with the appending
    return np.array(uList), np.array(vList), np.array(pList)




def windAndDrag(u, v, p):
    windSpeed = []
    force = []
    for i in range(len(u)):
        mag = np.hypot(u[i], v[i])
    
        speed = mag[::,0].mean()
        
        windSpeed.append(speed)
        
        force.append(horizontalForceOnSquare(p[i], speed))
        
    return windSpeed, force
    


#
#
# -------------------------------------------------------------------------------------
#
#
#                              PLOTTING THE FLUID FLOW
#
#
# -------------------------------------------------------------------------------------
#
#

    

def plotFlow(u, v, p, **kwargs):
    mask = True
    streams = False
    vectors = True
    pressure = True
    drag_plot = False
    
    if "mask" in kwargs.keys():
        mask = kwargs["mask"]
        
    if "streams" in kwargs.keys():
        streams = kwargs["streams"]
    
    if "vectors" in kwargs.keys():
        vectors = kwargs["vectors"]
    
    if "pressure" in kwargs.keys():
        pressure = kwargs["pressure"]
        
    if "drag_plot" in kwargs.keys():
        drag_plot = kwargs["drag_plot"]
        
    
    
    
    if drag_plot:
        fig = plt.figure(figsize = (12,8))
        
    else:
        fig = plt.figure(figsize = (12,4))
        
    plt.suptitle("Flow", fontweight = "bold")
    
    #a bunch of things useful to make things look nice
    transparent = (1, 1, 1, 0)
    object_cmap, norm = mcolors.from_levels_and_colors([0, 0.5, 1], ['black', transparent])
    a = 4
    b = int(a/2)
    
    if mask:
        maskArr = 1 - obj
        maskArr = maskArr.astype(bool)
        
        for i in range(len(u)):
            u[i][maskArr == True] = np.nan
            v[i][maskArr == True] = np.nan
            p[i][maskArr == True] = np.nan
    
    windSpeedList, forceList = windAndDrag(u, v, p)
    windSpeedList = np.array(windSpeedList)
    forceList = np.array(forceList)
    

    #call this method in FuncAnimation to get animation
    def update(frame):
        fig.clear()
        
        #magnitude of the velocity
        mag = np.hypot(u[frame], v[frame])
        
        
        plt.suptitle(f"Wind Speed: {round(windSpeedList[frame],2)} m/s", fontweight = "bold")
        
        tunnelPlot = 111
        if drag_plot:
            tunnelPlot = 211
            plt.subplot(212)
            plt.title("Drag Force vs. Wind Speed")
            plt.xlabel("Wind Speed (m/s)")
            plt.ylabel("Drag Force (N)")
            plt.plot(windSpeedList[:frame], forceList[:frame])
                
        
        if pressure:
            plt.subplot(tunnelPlot)
            plt.imshow(p[frame], alpha=0.5, cmap="viridis", extent = (np.min(x), np.max(x), np.min(y), np.max(y)))
            cbarp = plt.colorbar(location = "right")
            cbarp.ax.set_title("Pressure (Pa)")
        
        if vectors:
            plt.quiver(X[b::a, b::a], Y[b::a, b::a], u[frame, b::a, b::a], v[frame, b::a, b::a], mag[b::a, b::a], cmap = "plasma")
            cbarv = plt.colorbar(location = "left")
            cbarv.ax.set_title("Flow Speed (m/s)")
        
        if streams:
            plt.streamplot(X[b::a, b::a], Y[b::a, b::a], u[frame, b::a, b::a], v[frame, b::a, b::a], color = mag[b::a, b::a], cmap = "plasma")
            cbarv = plt.colorbar(location = "left")
            cbarv.ax.set_title("Flow Speed (m/s)")
            
        
        #black if object there
        plt.imshow(obj, cmap=object_cmap, extent = (np.min(x), np.max(x), np.min(y), np.max(y)))
        plt.gca().set_aspect('equal')
        
        
        
    
        return fig
    
    ani = animation.FuncAnimation(fig, update, interval=1, frames = len(uSol))
    
    return ani

#call the function
changeParameters(timeItterations = 150, force = 100000, viscosity = 0.00002, timeStep = 0.000001)
uSol, vSol, pSol = flow(u_initial, v_initial, p_initial)

anim = plotFlow(uSol, vSol, pSol, drag_plot = True)

print(time.time()-stTime)

