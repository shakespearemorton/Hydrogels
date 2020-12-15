#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 12:35:54 2020

@author: WilliamMorton
"""
import numpy as np
from movement import *
import cProfile
  
# --------- VARIABLES ---------
mass            = 1
lo              = 1                  #original length
lref            = 1                #
ko              = 1
n               = 0.2                 #viscosity
Fb              = 0.3               #breaking force
k               = ko*(lref/lo)
Fr              = 2e-5               #force rate
xl              = 21
yl              = 11
dt              = .001
tmax            = dt*1000
t               = 0
crack           = 5
  
# --------- CREATE SYSTEM ---------
pos,neigh       = squareNetwork(xl,yl,lo)
fixed           = fixbot(pos)
top             = topper(pos)
bonds           = pos[list(neigh[:,0])]-pos[list(neigh[:,1])]
w0              = np.sqrt(inner1d(bonds,bonds))
neigh,w0        = cracked(neigh,crack,w0)
prevpos         = pos
damp            = np.zeros(pos.shape)

# --------- RUN SIMULATION ---------
while t < tmax:
    pos,neigh,F,w0,crack = move2(pos, fixed, neigh, w0, k, crack, Fb, Fr, top, t, damp, n, dt)
    damp = pos-prevpos
    prevpos = pos
    if np.max(F) >= Fb:
            breaker(pos,fixed,Fprev,Fb)
            break
    Fprev = F
    printer(pos,neigh,t) #Prints an output at each timestep for viewing in Ovito as a LAMMPS bonds script
    t += dt