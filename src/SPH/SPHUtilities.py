from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Update the mass densities of the FluidNodeLists in a DataBase with the sum 
# definition
#-------------------------------------------------------------------------------
def sumSPHMassDensityImpl(dataBase, W, method):
    dataBase.updateConnectivityMap()
    cm = dataBase.connectivityMap()
    rho = dataBase.fluidMassDensity()
    pos = dataBase.fluidPosition()
    mass = dataBase.fluidMass()
    H = dataBase.fluidHfield()
    method(cm, W, pos, mass, H, rho)
    return

for dim in dims:
    exec("""
def sumSPHMassDensity%(dim)sd(dataBase, W):
    sumSPHMassDensityImpl(dataBase, W, computeSPHSumMassDensity%(dim)sd)
""" % {"dim" : dim})
