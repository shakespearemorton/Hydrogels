#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:26:39 2020

@author: WilliamMorton
"""

import numpy as np
from movement import *
  
# --------- VARIABLES ---------
mass = 1
lo = 1                  #original length
lref = 1                #
ko = 1
n = 0.2                 #viscosity
Fb = 0.3               #breaking force
k = ko*(lref/lo)
Fr = 2e-5               #force rate
xl = 5
yl = 4
dt = 1
tmax = dt*11
t = 0
crack = 2

# ----------- EMPTY INPUT ARRAYS ----------
initpos,neigh = squareNetwork(xl,yl,lo)
#for i in range(crack):
#    neigh[i]=[]
#    for j in range(len(neigh)):
#        if i in neigh[j]:
#            neigh[j].remove(i)
v = np.zeros((len(initpos),2))
top=[]
for i in range(len(initpos)):
    current = initpos[i]
    if current[1] == yl-1:
        top.append(i)
bot=[]
for i in range(len(initpos)):
    current = initpos[i]
    if current[1] == 0:
        bot.append(i)
prevpos=np.zeros(len(initpos))
pos = initpos
prevacc = np.zeros((len(initpos),2))
t=0
f = np.zeros((sum( [ len(listElem) for listElem in neigh]),2))

# ---------- RUN ELASTIC STRAIN ----------
while t < tmax:
    f , acc = force(pos,neigh,k,lo,top,dt,t,Fr,n,prevpos,mass,f)
    acc = 0.5*(acc + prevacc) 
    prevacc = acc
    disp = v*dt+0.5*dt**2*acc
    v = (disp/dt)-0.5*dt**2*acc
    prevpos = pos 
    pos += disp
    for i in bot:
        pos[i]=initpos[i]
    printer(pos,f,neigh,t)
    neigh,f = breaker(f,Fb,neigh,pos)
    plt.scatter(pos[:,0],pos[:,1])
    
    t += dt
