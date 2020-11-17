import scipy.special as sp
import numpy as np
import random
import scipy

def nlist( pos, nParticles, sigma2 ):
  neigh = [ [] for i in range( nParticles ) ]
  for i in range( nParticles ):
    for j in range( i+1, nParticles ):
      if distance2( pos[ i, : ], pos[ j, : ] ) < 2 * sigma2:
        neigh[ i ].append( j )  
        neigh[ j ].append( i )  
  return neigh

def distance2( x1, x2 ):
  x1 = np.array( x1 )
  x2 = np.array( x2 )
  r2 = ( ( x1 - x2 )**2 ).sum() 
  return r2

def rescale( r1, m, l, scaling, rfactor ):
  mod = np.sqrt( ( r1**2 ).sum() )
  teta = np.arctan( r1[ 1 ]  / r1[ 0 ] )
  if r1[ 0 ] > 0:
    phi = np.arccos( r1[ 2 ] / mod )
  else:
    phi = 2.0 * np.pi - np.arccos( r1[ 2 ] / mod ) 
  rvec = surface( m, l, phi, teta, scaling, rfactor )

  return rvec

def surface( m, l, p, t, scaling, rfactor ):
  r = scaling * ( rfactor + sp.sph_harm( m, l, p, t ).real)
  x = r * np.sin( p ) * np.cos( t )
  y = r * np.sin( p ) * np.sin( t )
  z = r * np.cos( p ) 
  return np.array( [ x, y, z ] ) 

def placeRandom( m, l, nParticles, scaling, rfactor ):
  pos = np.zeros( ( nParticles, 3 ) )
  for nn in range( nParticles ):
    check = True
    while check:
      i = random.uniform( -1, 1 )
      j = random.uniform( -1, 1 )
      if i**2 + j**2 <= 1.0:
        x = 2 * i * np.sqrt( 1.0 - ( i**2 + j**2 ) )
        y = 2 * j * np.sqrt( 1.0 - ( i**2 + j**2 ) )
        z = 1.0 - 2.0 * ( i**2 + j**2 ) 
        check = False
        r = np.sqrt( x**2 + y**2 + z**2 )
        x /= r
        y /= r
        z /= r
        teta = np.arctan( y / x )
        if x > 0:
          phi = np.arccos( z )
        else:
          phi = 2.0 * np.pi - np.arccos( z ) 
        rvec = surface( m, l, phi, teta, scaling, rfactor )
    pos[ nn, : ] = rvec
  return pos

def pairEnergy( rsq, eps, sigma2 ):
  r6 = ( sigma2 /rsq )**3
  r12 = r6**2 
  ene = 4.0 * eps * ( r12 - r6 )
  return ene

def totalEnergy( pos, nParticles, neigh, eps, sigma2, m, l ):
  totEne = 0.0
  for i in range( nParticles ):
    r1 = pos[ i, : ]
    for j in neigh [ i ]: 
      r2 = pos[ j, : ]
      d2 = distance2( r1, r2 )
      if d2 < 2**( 1.0 / 3.0 ) * sigma2:
        totEne += pairEnergy( d2, eps, sigma2 )
  return totEne

def setSystem( var ): 
  nParticles = var[ 'nParticles' ] 
  eps = var[ 'eps' ] 
  sigma2 = var[ 'sigma' ]**2
  m = var[ 'm' ]
  l = var[ 'l' ]
  temp = var[ 'temp' ]
  nSweeps = var[ 'nSweeps' ]
  maxdr = var[ 'maxdr' ]
  rfactor = var[ 'rfactor' ]
  projection = var[ 'Projection' ]
  if projection == 1:
    p_temp = scipy.linspace(0,6.28,200)
    t_temp = scipy.linspace(0,3.14,200)
    [p_temp,t_temp]=scipy.meshgrid(p_temp,t_temp)
    Y_temp = ( rfactor + sp.sph_harm( m, l, p_temp, t_temp ).real)
    maxi = np.amax(Y_temp)
    scaling = 10 / maxi
  elif projection == 0:   
    scaling = var[ 'scaling' ]
  SA = surfaceArea(l,m,scaling,rfactor)
  nParticles = int(SA) * nParticles
  #Place particles initially randomly
  pos = placeRandom( m, l, nParticles, scaling, rfactor )
  neigh = nlist( pos, nParticles, sigma2 )
  totEne = totalEnergy( pos, nParticles, neigh, eps, sigma2, m, l )
  # Try to minimise the energy. Particles must however stay on the surface
  # as well
  accepted = 0
  for nsw in range( nSweeps ):
    dr = np.zeros( 3 )
    dr[ 0 ] = maxdr * random.uniform( -1, 1 )
    dr[ 1 ] = maxdr * random.uniform( -1, 1 )
    dr[ 2 ] = maxdr * random.uniform( -1, 1 )
    atom = int( np.random.rand() * nParticles )
    ene = 0.0
    #current energy
    for i in neigh[ atom ]:
      rsq = distance2( pos[ atom, : ], pos[ i, : ] )
      ene += pairEnergy( rsq, eps, sigma2 )
    #new energy
    oldPos = pos.copy()
    pos[ atom, : ] += dr
    pos[ atom, : ] = rescale( pos[ atom, : ], m, l, scaling, rfactor ) 
    ene2 = 0.0
    for i in neigh[ atom ]: 
      rsq = distance2( pos[ atom, : ], pos[ i, : ] )
      ene2 += pairEnergy( rsq, eps, sigma2 )
    delta = ene2 - ene
    accept = True
    if delta > 0 :
      rand = np.random.rand( )
      if rand > np.exp( -delta / temp ):
        accept = False

    if accept:
      accepted += 1.0
      totEne += delta
    else:
      pos = oldPos

    nCheck = nParticles 

    if ( nsw % nCheck ) == 0:
      neigh = nlist( pos, nParticles, sigma2 )
      #print( "Step number", nsw )
      #print( "% accepted", accepted / nCheck )
      #print( "totEne", totEne )
      accepted = 0
#      with open( 'pos.xyz', 'a' ) as myF:
#        myF.write( "{0} \n".format( nParticles ) ) 
#        myF.write( "\n" ) 
#        for i in range( nParticles ):
#          myF.write( "C {0:.5f} {1:.5f} {2:.5f} \n".format( pos[ i, 0 ], pos[ i, 1 ], pos[ i, 2 ] ) ) 
      
  return pos

def writeLAMMPS(pos,nParticles):
    with open( 'system.txt', 'a' ) as g:
        g.write("#NP 10*10*20\n#second line will be skipped\n\n" )
        g.write( "{0:.0f}     atoms\n".format( nParticles+1) )
        g.write( "{0:.0f}     bonds\n".format( nParticles) )
        g.write("0     angles\n\n" )
        g.write("1 atom types\n1 bond types\n0 angle types\n\n-50 50 xlo xhi\n-50 50 ylo yhi\n-100 100 zlo zhi\n\nMasses\n\n2 1.0\n\n Atoms\n\n")
        part = 48501
        t=0
        for t in range (nParticles):
            g.write( "{0:.0f} 2 2 {1:.5f}. {2:.5f}. {3:.5f}. \n".format( part+t, pos[ t, 0 ], pos[ t, 1 ], pos[ t, 2 ] ) )
        t=0
        
        g.write ("{0:.0f} 2 2 {1:.5f}. {2:.5f}. {3:.5f}. \n".format( part+(nParticles), 0, 0, 0 ))
        g.write ("\nBonds\n\n")
        n1=0
        n=0
        while n1<nParticles:
            g.write ("{0:.0f} {1:.0f} {2:.0f} {3:.0f} \n".format( n+51801, 2, 48501+n1, part+(nParticles)))
            n+=1
            n1+=1
        g.write("\nAngles")
        g.close()
    return pos

def writeVMD(pos,nParticles):
    with open( 'VMD.txt', 'a' ) as g:
        g.write(repr(nParticles)+'\n\n')
        t=0
        for t in range(nParticles):
            g.write( "{0:.5f} {1:.5f} {2:.5f} \n".format( pos[ t, 0 ], pos[ t, 1 ], pos[ t, 2 ] ) )
            
    return pos
        
    
def rough( var ):
    m = var[ 'm' ]
    l = var[ 'l' ]
    scale = var[ 'scaling' ]
    rfactor = var[ 'rfactor' ]
    tot = 5000
    delta = tot
    p=0
    t=0
    rough = np.zeros((tot**2,4))
    rough2 = np.zeros((tot,tot))
    delta1 = 3.14/tot
    n=0
    m1=0
    m2=0
    w=0
    R1=0
    wave=[]
    while m2<tot:
        while m1<tot:
            R = scale*(rfactor+sp.sph_harm(m, l, p, t).real)
            if m1>0 and m1<tot-1:
                R1 = scale*(rfactor+sp.sph_harm(m, l, p, t+delta1).real)
                R2 = scale*(rfactor+sp.sph_harm(m, l, p, t-delta1).real)
                if R1<R and R2<R:
                    wave.append(t)
                    w+=1
            rough[n][0]=n
            rough[n][1]=R
            rough[n][2]=p
            rough[n][3]=t
            rough2[m1][m2]=R
            t+=delta
            n+=1
            m1+=1
        m1=0
        m2+=1
        t=-3.14/2
        p+=delta
    x=int(len(wave)/2)
    wavelength = np.zeros((x,1))
    n=0
    while n<x:
        wavelength[n] = np.abs(wave[n+1]-wave[n])
        n+=1
    n=1
    h=[]
    counter=0
    while n < tot**2:
        R = rough[n-1][1]
        R1 = rough[n][1]
        if n < (tot**2 - 1):
            R2 = rough[n+1][1]
            if R1 > R and R1 > R2:
                maxy=R1
                counter+=1
            if R1 < R and R1 < R2:
                miny=R1
                counter+=1
            if counter == 2:
                h.append(maxy-miny)
                counter = 1
        n+=1
    
    avg_wave=np.average(wavelength)
    avg_height=np.average(h)
    max_wave=np.max(wavelength)
    min_wave=np.min(wavelength)
    depth=np.max(rough2)-np.min(rough2)
    return(avg_wave, avg_height, max_wave, min_wave, depth)


def rescaler(pos,nParticles):
    nanosphere = np.zeros((nParticles,3))
    for i in range(len(pos)):
            nanosphere[i][0] = np.sqrt((pos[i][0]**2)+(pos[i][1]**2)+(pos[i][2]**2))
            
            if pos[i][0] > 0:
                nanosphere[i][1] = np.arccos(pos[i][2])
            else:
                nanosphere[i][1] = 2 * np.pi - np.arccos(pos[i][2])
            nanosphere[i][2] = np.arctan(pos[i][1]/pos[i][0])
    maxi = max(nanosphere[:,0])
    scale = maxi/10
    nanosphere[:,0] = nanosphere[:,0]/scale
    for i in range(len(pos)):
        pos[i][0] = nanosphere[i][0]*np.sin(nanosphere[i][1])*np.cos(nanosphere[i][2])
        pos[i][1] = nanosphere[i][0]*np.sin(nanosphere[i][1])*np.sin(nanosphere[i][2])
        pos[i][2] = nanosphere[i][0]*np.cos(nanosphere[i][1])
    return (pos)

def surfaceArea(l,m,scaling,rfactor):
    p_temp = scipy.linspace(0,6.28,200)
    t_temp = scipy.linspace(0,3.14,200)
    [p_temp,t_temp]=scipy.meshgrid(p_temp,t_temp)
    Y_temp = scaling*( rfactor + sp.sph_harm( m, l, p_temp, t_temp ).real)
    x = Y_temp*np.sin(t_temp)*np.cos(p_temp)
    y = Y_temp*np.sin(t_temp)*np.sin(p_temp)
    z = Y_temp*np.sin(t_temp)*np.cos(p_temp)
    w,n = z.shape
    area = 0
    v0 = np.zeros(3)
    v1 = np.zeros(3)
    v2 = np.zeros(3)
    v3 = np.zeros(3)
    for i in range(w-1):
        for j in range(n-1):
            v0[0] = x[i][j]
            v0[1] = y[i][j]
            v0[2] = z[i][j]
            v1[0] = x[i][j+1]
            v1[1] = y[i][j+1]
            v1[2] = z[i][j+1]
            v2[0] = x[i+1][j]
            v2[1] = y[i+1][j]
            v2[2] = z[i+1][j]
            v3[0] = x[i+1][j+1]
            v3[1] = y[i+1][j+1]
            v3[2] = z[i+1][j+1]
            a = v1 - v0;
            b = v2 - v0;
            c = v3 - v0;
            A =  0.5*((np.sum(np.cross(a,c)**2))**0.5+(np.sum(np.cross(b,c)**2))**0.5)
            area+=A
    print(area)
    return(area)