#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 20:51:13 2024

@author: zackaryschnitzer
"""

#   This file will simply be the class definition for Eulerian tiles
#
#   Each tile will contain all of the information needed to describe
#   the gas in a small grid block.
#
#   This is just the tiles. NOT the whole grid
#
#   NOTE: this uses a staggered grid to store velocity values for the sides,
#   but also information about the center of the cell. General convention is that
#   the values for the sides are left and top parts of the cell.
#   
#   
#   
#
    
class eulerTile:
    
    
    def __init__(self, pressure, velocity, smokeDensity):
        """
        Parameters
        ----------
        pressure : float
            pressure.
        velocity : tuple
            (u,v) component of velocity of left, and top edge resp.
        smokeDensity : float
            density of smoke in sim.

        Returns
        -------
        None.

        """
        self.pressure = pressure
        self.u, self.v = velocity
        self.smokeDensity = smokeDensity


        
        