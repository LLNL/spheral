from SpheralModules.Spheral.SPHSpace import *

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

def sumSPHMassDensity1d(dataBase, W):
    sumSPHMassDensityImpl(dataBase, W, computeSPHSumMassDensity1d)

def sumSPHMassDensity2d(dataBase, W):
    sumSPHMassDensityImpl(dataBase, W, computeSPHSumMassDensity2d)

def sumSPHMassDensity3d(dataBase, W):
    sumSPHMassDensityImpl(dataBase, W, computeSPHSumMassDensity3d)
