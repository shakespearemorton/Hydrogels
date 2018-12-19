import numpy as np
import math

#inputs
xlo = -50
xhi = 50
ylo = -50
yhi = 50
zlo = -100
zhi = 100
zgap = 20
spacing = 10
layers = 2
total = 0
atoms = 3876582
bonds = 160292
angles = 12341334
m1 = 1
m2 = 0.1
type_atoms = 2
type_bonds = 2
type_angles = 2
type_masses = 2
N=100
n=0
r=0.8
bond = 1
atomID = np.zeros((86862,4))
bondID = np.zeros((10*644120,2))

#prepare document
#g= open("lammps_nanoparticles.txt","w+")
#g.write('#Gels with nanoparticles \n')
#g.write('#Second line will be skipped \n \n')
#g.write(repr(atoms) +' atoms\n' )
#g.write(repr(bonds) +' bonds\n' )
#g.write(repr(angles) +' angles\n' )
#g.write(repr(type_atoms) +' atom types\n' )
#g.write(repr(type_bonds) +' bond types\n' )
#g.write(repr(type_angles) +' angle types\n \n' )
#g.write(repr(xlo) +' ' +repr(xhi) +' xlo xhi\n' )
#g.write(repr(ylo) +' ' +repr(yhi) +' ylo yhi\n' )
#g.write(repr(zlo) +' ' +repr(zhi) +' zlo zhi\n \n' )
#g.write('Masses \n \n')
#g.write('1 ' +repr(m1) +'\n')
#g.write('2 ' +repr(m2) +'\n\n')
#g.write('Atoms \n')


#polymer placement
xstart = xlo
xfinish = xhi
ystart = ylo
yfinish = yhi
zstart = zlo
zfinish = zhi

xstep = xstart
ystep = ystart
zstep = zstart

while zstep <= zfinish:
    while ystep <= yfinish:
        while xstep <= xfinish:
            if abs(xstep) == xhi or abs(ystep) == yhi or abs(zstep) == zhi:
                pass
            else:
                line=repr(total) + ' 1 '+ '1 ' +repr(xstep) +' ' +repr(ystep) +' ' +repr(zstep) +'\n'
                atomID[total][0] = total
                atomID[total][1] = xstep
                atomID[total][2] = ystep
                atomID[total][3] = zstep
                #g.write(line)
                total += 1
            xstep += 1
        ystep += spacing
        xstep = xstart
    zstep += spacing
    ystep = ystart    
    if zstep == 0:
        zstep = 0+zgap/2
    else:
        pass
    
xstep = xstart
ystep = ystart
zstep = zstart
zfinish = zhi

while zstep <= zfinish:
    while xstep <= xfinish:
        while ystep <= yfinish:
            if abs(xstep) == xhi or abs(ystep) == yhi or abs(zstep) == zhi:
                pass
            else:
                line=repr(total) + ' 1 '+ '1 ' +repr(xstep) +' ' +repr(ystep) +' ' +repr(zstep) +'\n'
                atomID[total][0] = total
                atomID[total][1] = xstep
                atomID[total][2] = ystep
                atomID[total][3] = zstep
                #g.write(line)
                total += 1
            ystep += 1
        xstep += spacing
        ystep = ystart
    zstep += spacing
    xstep = xstart
    if zstep == 0:
        zstep= 0+zgap/2
    else:
        pass
  
xstep = xstart
ystep = ystart
zstep = zstart
zfinish = zhi
    
while xstep <= xfinish:
    while ystep <= yfinish:
        while zstep <= zfinish:
            if abs(xstep) == xhi or abs(ystep) == yhi or abs(zstep) == zhi:
                pass
            else:
                line=repr(total) + ' 1 '+ '1 ' +repr(xstep) +' ' +repr(ystep) +' ' +repr(zstep) +'\n'
                atomID[total][0] = total
                atomID[total][1] = xstep
                atomID[total][2] = ystep
                atomID[total][3] = zstep
                # g.write(line)
                total += 1
            zstep += 1
            if zstep == 0-zgap/2 + 1:
                zstep = 0+zgap/2
            else:
                pass
        ystep += spacing
        zstep = zstart
    xstep += spacing
    ystep = ystart


#nanoparticle parameters
f=5
theta_scale=10
tau_scale=15
N=100
n=0
r=0.8
total_use = total
i=0
#nanopartice placement
xstart = xlo + 0.5
xfinish = xhi - 0.5
ystart = ylo + 0.5
yfinish = yhi - 0.5
zstart = 0 - zgap/4
zfinish = 0 + zgap/4
xstep = xstart
ystep = ystart
zstep = zstart
while zstep <= zfinish:
    while ystep <= yfinish:
        while xstep <= xfinish:
            #center of the NP
            line = repr (total) + ' 2' + ' 2 ' + repr(xstep) + ' ' + repr(ystep) + ' ' + repr(zstep)
            atomID[total][0] = total
            atomID[total][1] = xstep
            atomID[total][2] = ystep
            atomID[total][3] = zstep
            i = total
            total += 1
           # g.write(line)
            while n<=N:
                c1 = np.random.random(1)
                c2 = np.random.random(1)
                theta = math.asin(c1)
                tau = math.asin(c2)
                R=r+((1/f)*math.sin(theta*theta_scale)*math.cos(tau*tau_scale))
                xpart = xstep + R*math.sin(tau)*math.cos(theta)
                ypart = ystep + R*math.sin(tau)*math.sin(theta)
                zpart = zstep + R*math.cos(tau)
                line=repr(total) + ' 2 '+ '2 ' +repr(xpart) +' ' +repr(ypart) +' ' +repr(zpart) +'\n'
                atomID[total][0] = total
                atomID[total][1] = xpart
                atomID[total][2] = ypart
                atomID[total][3] = zpart
                bondID[bond-1][0] = i
                #g.write(line)
                bondID[bond-1][1] = total
                total += 1
                bond +=1
                n +=1
            xstep += spacing
            n=0
        ystep += spacing
        xstep = xstart
    zstep += zgap/2
    ystep = ystart
#g.write('\nBonds \n')
print ('done with placement')

#gel bond calculations
i=0
j=0
while i < total_use:
    while j < total_use:
        x1 = (atomID[i][1])
        y1 = (atomID[i][2])
        z1 = (atomID[i][3])
        x2 = (atomID[j][1])
        y2 = (atomID[j][2])
        z2 = (atomID[j][3])
        test = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        if test == spacing:
            line = repr(bond) + ' 1 ' + repr(i) + ' ' + repr(j) +'\n'
            bondID[bond-1][0]=i
            bondID[bond-1][1]=j
            #g.write(line)
            bond +=1
        j +=1
    i += 1
    j=0

print('done with bonds')

#gel angle 1 calculation
i=0
j=0
k=0
angle=1
angleID = np.zeros((3*total_use,4))
while k < total_use:
    while i < total_use:
        while j < total_use:
            x1 = (atomID[i][1])
            y1 = (atomID[i][2])
            z1 = (atomID[i][3])
            x2 = (atomID[j][1])
            y2 = (atomID[j][2])
            z2 = (atomID[j][3])
            x3 = (atomID[k][1])
            y3 = (atomID[k][2])
            z3 = (atomID[k][3])
            test = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
            test2 = math.sqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
            if test+test2 == 2*spacing:
                test_a = math.sqrt((x1-x3)**2+(y1-y3)**2+(z1-z3)**2)
                if test_a == 2*spacing:
                    a=2
                else:
                    a=1
                if atomID[i][0] == atomID[j][0] or atomID[i][0] == atomID[k][0] or atomID[j][0] == atomID[k][0]:
                    pass
                else:
                    line = repr(angle) + repr(a) + ' ' + repr(i) + ' ' + repr(j) + ' ' +repr(k) + '\n'
                    angleID[angle-1][0]=a
                    angleID[angle-1][1]=i
                    angleID[angle-1][2]=j
                    angleID[angle-1][3]=k
                    #g.write(line)
                    angle +=1
            j +=1
        i += 1
        j=0
    k +=1
    i=0

print('done with angles')

    
    

#g.close()
print ('done')
