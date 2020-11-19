#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:12:55 2020

@author: WilliamMorton
"""
import numpy as np
from matplotlib import pylab as plt

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
 
    return (pos,neigh)

# ---------- Calculates radial distance between two points ----------            
def radDistance( x1, x2 ):
  x1 = np.array( x1 )
  x2 = np.array( x2 )
  r2 = ( ( x1 - x2 )**2 ).sum() 
  return (r2)            
    
# ---------- Calculates force on each node ----------
#this definitely isn't working yet, but as of now it pulls the top level up
def force(pos,neigh,k,lo,top,dt,t,Fr,n,prevpos,mass,f):
    acc = np.zeros((len(pos),2))
    k=0
    for i in range(len(pos)):
        for j in neigh[i]:
            l = radDistance(pos[i],pos[j])
            #hook = k*(l-lo)*((np.array(pos[j])-np.array(pos[i]))/l)
        if i in top:
            f[k,1] += Fr*(t+1)
            
        acc[i]+=f[k]/mass
        k+=1
    return(f,acc)


# ---------- If the force between the bonds is greater than the breaking force -> delete bond ----------
#if the bond is broken, the force goes to 0
def breaker(f,Fb,neigh,pos):
    k=0
    m=0
    for i in range(len(pos)):
        for j in neigh[i]:
            if f[k].sum() > Fb:
                neigh[i].remove(j)
                f[k] = 0
            k+=1
    return(neigh,f)


# ---------- Prints an output at each timestep for viewing in Ovito as a LAMMPS bonds script ----------
def printer(pos,f,neigh,t):
    with open( 'pos'+repr(t)+'.txt', 'w' ) as g:
        g.write("#Elastic Network\n#second line will be skipped\n\n" )
        g.write( "{0:.0f} atoms\n".format( len(pos)) )
        g.write( "{0:.0f} bonds\n".format( len(f)) )
        g.write("0 angles\n\n" )
        g.write("1 atom types\n1 bond types\n0 angle types\n-50 50 xlo xhi\n-50 50 ylo yhi\n-100 100 zlo zhi\n\nMasses\n\n1 1.0\n\nAtoms\n\n")

        for t in range(len(pos)):
            g.write( "{0:.0f} 1 1 {1:.0f}. {2:.0f}. {3:.0f}.\n".format( t+1, pos[ t, 0 ], pos[ t, 1 ],0) )

        g.write ("\nBonds\n\n")
        i=0
        for k in range(len(pos)):
            for j in neigh[k]:
                g.write ("{0:.0f} {1:.0f} {2:.0f} {3:.0f} \n".format( i+1, 1, k+1, j+1))
                i+=1
        g.write("\nAngles")
        g.close()
    return()