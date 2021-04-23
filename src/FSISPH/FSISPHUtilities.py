from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()


for dim in dims:
    exec("""
def computeSurfaceNormals(dataBase, W):
    dataBase.updateConnectivityMap()
    cm = dataBase.connectivityMap()
    rho = dataBase.fluidMassDensity()
    pos = dataBase.fluidPosition()
    mass = dataBase.fluidMass()
    H = dataBase.fluidHfield()
    computeSurfaceNormals%(dim)sd(cm, W, pos, mass, H, rho)
    
""" % {"dim" : dim})
