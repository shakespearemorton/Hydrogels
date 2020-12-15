#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:12:55 2020

@author: WilliamMorton
"""
import numpy as np
from matplotlib import pylab as plt
from numpy.core.umath_tests import inner1d

# ---------- Develop square 2D network and bonds ----------
def squareNetwork(xl,yl,lo):
    pos = []
    neigh=[]
    for i in range(yl):
        for j in range(xl):
            pos.append([j*lo,i*lo])
    xl = xl - 1
    yl = yl - 1
    for k in pos:
        i,j = k
        if i == 0:
            if j == 0:
                neigh.append([pos.index([i,j+1]),pos.index([i+1,j])])
            elif j == yl:
                neigh.append([pos.index([i,j-1]),pos.index([i+1,j])])
            else:
                neigh.append([pos.index([i,j-1]),pos.index([i,j+1]),pos.index([i+1,j])])
        elif i == xl:
            if j == 0:
                neigh.append([pos.index([i-1,j]),pos.index([i,j+1])])
            elif j == yl:
                neigh.append([pos.index([i-1,j]),pos.index([i,j-1])])
            else:
                neigh.append([pos.index([i,j-1]),pos.index([i,j+1]),pos.index([i-1,j])])
        elif j == yl:
            neigh.append([pos.index([i,j-1]),pos.index([i+1,j]),pos.index([i-1,j])])
        elif j == 0:
            neigh.append([pos.index([i,j+1]),pos.index([i+1,j]),pos.index([i-1,j])])
        else:
            neigh.append([pos.index([i,j-1]),pos.index([i+1,j]),pos.index([i-1,j]),pos.index([i,j+1])])
    new = []
    for i in range(len(neigh)):
        for j in neigh[i]:
            current = [j,i] in new
            if  current == False:
                new.append([i,j])
    neigh = np.array(new)
    pos = np.asarray(pos).astype(float)
    return (pos,neigh)

# ---------- Finds nodes at the bottom of the array, which can then be fixed (0) or free(1) ----------                     
def fixbot(pos):
    bot=[]
    for i in range(len(pos)):
        current = pos[i]
        if current[1] == 0:
            bot.append(i)
    fixed = np.ones(len(pos))
    for i in bot:
        fixed[i]=0
    return(fixed)

# ---------- Finds nodes at the top of the array, which are then pulled up ----------   
def topper(pos):
    top=[]
    for i in range(len(pos)):
        current = pos[i]
        if current[1] == np.max(pos[:,1]):
            top.append(i)
    return(top)
    
# ---------- Calculates force on each node ----------
def move(pos, fixed, neigh, w0, k, crack, Fb, Fr, top, t):
    F = np.zeros(pos.shape)
    cycles=1000
    precision=0.001
    for i in range(cycles):
        F = np.zeros(pos.shape)
        F[top,1]+=Fr*t
        bonds = pos[neigh[:,1]]-pos[neigh[:,0]] 
        bonds = np.asarray(bonds).astype(float)
        w1 = np.sqrt(inner1d(bonds,bonds)) 
        bonds/=w1[:,None] 
        f = (w1 - w0)
        f = f[:,None] * bonds
        np.add.at(F, neigh[:,0], f)
        np.subtract.at(F, neigh[:,1], f)
        pos += k * F * 1 * fixed[:,None]
        if np.amax(F) < precision:
            break
    return (pos,neigh,F,w0,crack)


# ---------- Calculates force on each node ----------
def move2(pos, fixed, neigh, w0, k, crack, Fb, Fr, top, t, damp, n, dt):
    F = np.zeros(pos.shape)
    F[top,1]+=Fr*t
    bonds = pos[neigh[:,1]]-pos[neigh[:,0]] 
    bonds = np.asarray(bonds).astype(float)
    w1 = np.sqrt(inner1d(bonds,bonds)) 
    bonds/=w1[:,None] 
    f = (w1 - w0)
    f = f[:,None] * bonds * (n*damp/dt)
    np.add.at(F, neigh[:,0], f)
    np.subtract.at(F, neigh[:,1], f)
    pos += k * F * 1 * fixed[:,None]
    return (pos,neigh,F,w0,crack)

# ---------- Initiates a crack ----------   
def cracked(neigh,crack,w0):
    b = []
    for i in range(len(neigh)):
        if neigh[i,0]+1 <= crack or neigh[i,1]+1 <= crack:
            b.append(i)
    neigh = np.delete(neigh,b,0)
    w0 = np.delete(w0,b,0)
    return(neigh,w0)

# ---------- Output for when break occurs ----------   
def breaker(pos,fixed,F,Fb):
    test=[]
    for i in range(len(pos)):
        if fixed[i] == 0:
            x = F[i].sum()/Fb
            if x != 0:
                test.append(x)
    d = np.arange(1,len(test)+1)
    plt.loglog(d,test)
    return()

# ---------- Prints an output at each timestep for viewing in Ovito as a LAMMPS bonds script ----------
def printer(pos,neigh,t):
    with open( 'pos'+repr(t)+'.txt', 'w' ) as g:
        g.write("#Elastic Network\n#second line will be skipped\n\n" )
        g.write( "{0:.0f} atoms\n".format( len(pos)) )
        g.write( "{0:.0f} bonds\n".format( len(neigh)) )
        g.write("0 angles\n\n" )
        g.write("1 atom types\n1 bond types\n0 angle types\n-50 50 xlo xhi\n-50 50 ylo yhi\n-100 100 zlo zhi\n\nMasses\n\n1 1.0\n\nAtoms\n\n")

        for t in range(len(pos)):
            g.write( "{0:.0f} 1 1 {1:.0f}. {2:.0f}. {3:.0f}.\n".format( t+1, pos[ t, 0 ], pos[ t, 1 ],0) )

        g.write ("\nBonds\n\n")
        for k in range(len(neigh)):
            g.write ("{0:.0f} {1:.0f} {2:.0f} {3:.0f} \n".format( k+1, 1, neigh[k,0]+1, neigh[k,1]+1))
        g.write("\nAngles")
        g.close()
    return()