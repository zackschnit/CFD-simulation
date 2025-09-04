#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 09:46:07 2024

@author: zackaryschnitzer
"""

import turtle 
import random 
  
# global colors 
col = ['red', 'yellow', 'green', 'blue', 
       'white', 'black', 'orange', 'pink'] 
  
# method to call on screen click 
def draw():

    global sc
    global points
    points = []
    
    
    def fxn(x,y):


        print(x)
        print(y)
        if x < -10 and y < -10:
            sc.exitonclick()
            return points
        
        points.append((x,y))
        turtle.goto(x,y)

    sc.onclick(fxn)
    
    
    

# set screen 
sc = turtle.Screen() 
sc.setup(400, 300) 

# call method on screen click 
print(draw())