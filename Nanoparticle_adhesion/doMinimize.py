from minimize2 import *

var = dict()
var[ 'eps' ] = 1.0
var[ 'sigma' ] = 1.0 
var[ 'm' ] = REPLACEM
var[ 'l' ] = REPLACEL
var[ 'temp' ] = 0.10
var[ 'nSweeps' ] = 10**(3)
var[ 'maxdr' ] = 0.1
var[ 'scaling' ] = REPLACES
var[ 'rfactor' ] = REPLACER
var[ 'nParticles' ] = 1 #Scale so that nparticles = 1 * SurfaceArea
var[ 'Projection'] = REPLACEP # if Projection area = 0 it is not conserved. If Projection area = 1, it is scaled so the maximum radius is 10 (i.e. scaling inputed by user is redundant)
pos = setSystem( var )
nParticles = var[ 'nParticles' ]
#pos = writeVMD(pos,i,nParticles)
pos = writeLAMMPS(pos,nParticles)