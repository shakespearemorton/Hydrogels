#Imports

from __future__ import division
import random
import scipy.special as sp
import numpy as np
import math

#Variables
w=0
l = 9    # degree
m = 4    # order
N = 1000  # number of particles
f = 1 #scaling (mesh size)
s = 0 #average radius
rfactor = 2.5
difference = 0.5
sigma = 1

#potential: Y = 1-r
limit1=rfactor - (difference/2)
limit2=rfactor + (difference/2)
ymax = (limit2-limit1)/2
rcut = 1-ymax
rdel = 0.3*sigma
rverlet = rcut + rdel


#Start this shit
#g= open("mid.txt","w+")
g= open("vmd.txt","w+")
g.write(repr(N*2) + '\n')
g.write('nanoparticle' + '\n')

#Place Nanoparticles

atomID = np.zeros((N+1,5))
n=0

while n<=N:
    r=random.uniform(-1,1)
    r2=random.uniform(-1,1)
    if r**2 + r2**2 < 1:
        x=2*r*math.sqrt(1-(r**2+r2**2))
        y=2*r2*math.sqrt(1-(r**2+r2**2))
        z=1-(2*(r**2+r2**2))
        x= x / math.sqrt( x**2 + y**2 + z**2 )
        y= y / math.sqrt( x**2 + y**2 + z**2 )
        z= z / math.sqrt( x**2 + y**2 + z**2 )
        t= np.arctan(y/x)
        #if t < 0:
            #t = abs(t)
        p= np.arccos(z/math.sqrt(x**2+y**2+z**2))
       # if p < 0:
           # p = abs(p)
        #p, t = np.mgrid[0:2*np.pi:300j, 0:np.pi:150j]
        R = sp.sph_harm(m, l, p, t).real
        X = x*(rfactor + R)
        Y = y*(rfactor + R)
        Z = z*(rfactor + R)
  
        atomID[n][0] = n
        atomID[n][1] = X
        atomID[n][2] = Y
        atomID[n][3] = Z
        g.write( 'Si ' + repr(X) + ' ' + repr(Y) + ' ' + repr(Z) + '\n')
        n+=1
    else:
        pass

#Meet the Neighbors
        
Neighbors = np.zeros((N+1,N))
i=0
j=0
k=2
r_i = 0
while i<=N:
    while j<=N:
        if j == i:
            j += 1
        if j < N+1:
            r = math.sqrt((atomID[i][1]-atomID[j][1])**2 + (atomID[i][2]-atomID[j][2])**2 + (atomID[i][3]-atomID[j][3])**2 )
            r_i += r
            Neighbors [i][1] = i
            if r < rcut:
                Neighbors [i][k] = j
                k +=1
        j+=1
    Neighbors [i][0] = k-2
    k=2
    j=0
    i+=1
r_i = r_i/N**2
    
#Honey I Minimized the Energy
deltat = 0.1
tmax = deltat*500000000
time=0
while time <=tmax:
#while w < 10:
    i = random.randint(0,N)
    k=2
    disp=0
    dispx = 0
    dispy = 0
    dispz = 0
    totpot1 = 0
    totpot2 = 0
    while k <= Neighbors [i][0]+2:
        j = int(Neighbors [i][k])
        if i == j:
            pass
        else:
            r = math.sqrt((atomID[i][1]-atomID[j][1])**2 + (atomID[i][2]-atomID[j][2])**2 + (atomID[i][3]-atomID[j][3])**2 )
            dispx += ((atomID[i][1]-atomID[j][1])/r)
            dispy += ((atomID[i][2]-atomID[j][2])/r)
            dispz += ((atomID[i][3]-atomID[j][3])/r)
            totpot1 += (1 - r)
        k+=1
    temp_posx =atomID [i][1] + dispx
    temp_posy =atomID [i][2] + dispy
    temp_posz =atomID [i][3] + dispz
    k=2
    while k <= Neighbors [i][0]+2:
        j = int(Neighbors [i][k])
        r = math.sqrt(((atomID[i][1]+dispx)-atomID[j][1])**2 + ((atomID[i][2]+dispy)-atomID[j][2])**2 + ((atomID[i][3]+dispz)-atomID[j][3])**2 )
        totpot2 += (1 - r)
        k+=1
    if totpot2 < totpot1:
        t= np.arctan(temp_posy/temp_posx)
        p= np.arccos(temp_posz/math.sqrt(temp_posx**2+temp_posy**2+ temp_posz**2))
        R = sp.sph_harm(m, l, p, t).real
        test1 = limit1+R
        test2 = limit2+R
        r = math.sqrt ((atomID[i][1]+dispx)**2 + (atomID[i][2]+dispy)**2 + (atomID[i][3]+dispz)**2)
        if r < test1:
            pass
        elif test2 < r:
            pass
        else:
            w +=1
            disp = math.sqrt(dispx**2 + dispy**2 + dispz**2 )
            atomID [i][4] += disp
            atomID [i][1] =atomID [i][1] + dispx
            atomID [i][2] =atomID [i][2] + dispy
            atomID [i][3] =atomID [i][3] + dispz
    if atomID[i][4] == disp:
        pass
    elif atomID[i][4] == 0:
        pass
    else:
        i=0
        j=1
        k=2
        while i<=N:
            while j<=N:
                if j == i:
                    j += 1
                if j < N+1:
                    r = math.sqrt((atomID[i][1]-atomID[j][1])**2 + (atomID[i][2]-atomID[j][2])**2 + (atomID[i][3]-atomID[j][3])**2 )
                    Neighbors [i][1] = i
                    if r < rcut:
                        Neighbors [i][k] = j
                        k +=1
                j+=1
            atomID[i][4] = 0
            Neighbors [i][0] = k-2
            k=2
            j=0
            i+=1
        
    time += deltat
#Meet the Neighbors
i=0
j=0
r_f = 0
while i<=N:
    while j<=N:
        if j == i:
            j += 1
        if j < N+1:
            r = math.sqrt((atomID[i][1]-atomID[j][1])**2 + (atomID[i][2]-atomID[j][2])**2 + (atomID[i][3]-atomID[j][3])**2 )
            r_f += r
        j+=1
    j=0
    i+=1
r_f = r_f/N**2
    

n=0       
while n <N:
        X=atomID[n][1]
        Y=atomID[n][2]
        Z=atomID[n][3]
        g.write( 'Au ' + repr(X) + ' ' + repr(Y) + ' ' + repr(Z) + '\n')
        n +=1
g.write('Average initial spacing ' + repr(r_i) +' Average final spacing ' + repr(r_f) + ' in ' + repr(w) + 'steps')        
g.close()
    


    
    
