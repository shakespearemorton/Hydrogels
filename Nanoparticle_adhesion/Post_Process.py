from __future__ import division
import random
import scipy.special as sp
import numpy as np
import math
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate
import os.path

def loadnpdata(typer):
    file =('dump.system.0')
    tester = os.path.exists(str(file))
    y=0
    if tester == False:
        pass
    else:
        f=open(file,"r")
        lines=f.readlines()

        poly = np.zeros((48500,3))
        nano = np.zeros((len(lines)-48509,3))
        bridger = np.zeros(48500)
        i=0
        j=0
        for x in lines:
            if y < 9:
                pass
            else:
                test = int((x.split(' ')[1]))
                if test == 1:
                    poly[i][0] = (x.split(' ')[2])
                    poly[i][1] = (x.split(' ')[3])
                    poly[i][2] = (x.split(' ')[4])
                    temp = int((x.split(' ')[0]))
                    bridger[i] = typer[temp-1][1]
                    i+=1
                elif test ==2:
                    nano[j][0] = (x.split(' ')[2])
                    nano[j][1] = (x.split(' ')[3])
                    nano[j][2] = (x.split(' ')[4])
                    j+=1
            y+=1
        f.close()
        return poly,nano,bridger

def distance2( x1, x2 ):
  x1 = np.array( x1 )
  x2 = np.array( x2 )
  r2 = ( ( x1 - x2 )**2 ).sum() 
  return r2    

    
def proxi( nano,poly,sigma,bridger ):
  neigh = [ [] for i in range( len(poly) ) ]
  num_attached = 0
  tot_contact = 0
  Janet = np.zeros((785,3))
  bridge1 = []
  bridge3 =[]
  for i in range( len(poly) ):
    attached = 0
    if poly[i][2] < 16:
        if poly[i][2]>-16:
            if poly[i][0]<16:
                if poly[i][0] > -16:
                    if poly[i][1]<16:
                        if poly[i][1]>-16:
                            for j in range( len(nano) ):
                              where = distance2( poly[ i, : ], nano[ j, : ] )
                              if where > sigma:
                                  if where < 2.5:
                                      attached = 1
                                      neigh[ i ].append( j )  
                                #neigh[ j ].append( i )  
                            if attached == 1:
                                Janet[num_attached,:] = poly[i,:]
                                num_attached +=1
                                tot_contact += len(neigh[i])
                                if bridger[i] == 1:
                                    bridge1.append(tot_contact)
                                elif bridger[i] == 3:
                                    bridge3.append(tot_contact)
  return neigh, num_attached, tot_contact, Janet, bridge1, bridge3

def writeVMD(pos,i,nParticles):
    with open( 'VMD.txt', 'a' ) as g:
        g.write(repr(nParticles)+'\n\n')
        t=0
        for t in range(nParticles):
            g.write( "{0:.5f} {1:.5f} {2:.5f} \n".format( pos[ t, 0 ], pos[ t, 1 ], pos[ t, 2 ] ) )
            
    return pos


def loadgel(gel):
    f=open("gel.txt",'r')
    lines=f.readlines()
    typer = np.zeros((48500,2))
    i=0
    y=0
    for x in lines:
        if y < 21:
            pass
        elif y > 48520:
            pass
        else:
            test = int((x.split(' ')[5]))
            if test > 0:
                typer[i][0] = (x.split(' ')[0])
                typer[i][1] = 1
                i+=1
            elif test < 0:
                typer[i][0] = (x.split(' ')[0])
                typer[i][1] = 3
                i+=1
        y+=1
    f.close()
    return typer

w=0
m=1
y=0
ss=0

Summary = np.zeros(1,7)
file =('sscurve_system.txt')
tester = os.path.exists(str(file))
if tester == False:
    pass
else:
    f=open(file,"r")
    lines=f.readlines()
    info = np.zeros((len(lines),2))
    for x in lines:
        if y ==0:
            pass
        else:
            info[y][0] = (x.split(' ')[0])
            info[y][1] = (x.split(' ')[1])
        y+=1
    f.close()
    
        
    
    n=2
    z=2000
    m=0
    sume = 0
    count = 0
    deltax = (info[3][0]-info[2][0])/2
    y = len(info)
    Summary[w][0]=ss
    #Summary[w][1]=info[:,m+1].max()
    marker = 0
    while z<=(len(info)+1):
        #if info[z][1] == Summary[w][1]:
        #    marker = 1
        #if marker == 1: 
        if info[z][1] < 0:
            Summary[w][2] = info[z][0]
            zfinal = z
            z = len(info)+1
        if z == len(info)-1:
            Summary[w][2] = 666
            z+=2
        z+=1
    if Summary[w][2] == 666:
        Summary[w][3] = 666
    data = np.zeros((zfinal+1,2))
    for i in range(zfinal):
        data[i][0]=info[i][0]
        data[i][1]=info[i][1]
    Summary[w][1]=data[:,1].max()
    Summary[w][3]=np.trapz(data[:,1],data[:,0])

w=0
m=1
y=0
gel = str('gel')
typer = loadgel(gel)
Touchy = np.zeros((1,3))
file =('dump.system.0')
tester = os.path.exists(str(file))
if tester == False:
    pass
else:
    poly,nano,bridger= loadnpdata(ss,typer)
    sigma = 1
    neigh,num_attached,tot_contact,Janet, bridge1,bridge3 = proxi(nano,poly,sigma,bridger)
    if num_attached > 0:
        effective = (tot_contact/num_attached)*5
    else:
        effective = 5
    Summary[0][4]=num_attached
    Summary[0][5]=effective
    B = num_attached / ((abs(len(bridge1)-len(bridge3)))+1)
    Summary[0][6]=B

Janet = writeVMD(Janet,0,len(Janet))
nano = writeVMD(nano,1,len(nano))


mat = np.matrix(Summary)
with open('_Summary.txt','wb') as f:
for line in mat:
    np.savetxt(f, line, fmt='%.2f')