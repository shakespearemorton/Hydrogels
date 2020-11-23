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

pos,neigh       = squareNetwork(xl,yl,lo)
pos             = np.insert(pos, 2, 0, axis=1)
fixed           = fixbot(pos)
top             = topper(pos)
neigh           = cracked(neigh,crack)
bonds           = pos[list(neigh[:,0])]-pos[list(neigh[:,1])]
w0              = np.sqrt(inner1d(bonds,bonds))

while t < tmax:
    for i in top:
        pos[i,1] += 0.1
    pos,neigh = move(pos, fixed, neigh, w0, k, Fb)
    printer(pos,neigh,t)
    t += dt