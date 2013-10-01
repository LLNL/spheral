from Spheral import *
from SpheralTestUtilities import *

def TestNodes(nodesName):
    nodes = eval(nodesName, globals())

    output(nodesName + ".numNodes")
    output(nodesName + ".numInternalNodes")
    output(nodesName + ".numGhostNodes")
    output(nodesName + ".numFields")
    output(nodesName + ".valid")
    output(nodesName + ".eos")

    output(nodesName + ".mass")
    output("len(" + nodesName + ".mass)")
    output(nodesName + ".haveField(" + nodesName + ".mass)")
    output(nodesName + ".mass[:]")

    output(nodesName + ".positions")
    output("len(" + nodesName + ".positions)")
    output(nodesName + ".haveField(" + nodesName + ".positions)")
    output(nodesName + ".positions[:]")

    output(nodesName + ".massDensity")
    output("len(" + nodesName + ".massDensity)")
    output(nodesName + ".haveField(" + nodesName + ".massDensity)")
    output(nodesName + ".massDensity[:]")

    output(nodesName + ".specificThermalEnergy")
    output("len(" + nodesName + ".specificThermalEnergy)")
    output(nodesName + ".haveField(" + nodesName + ".specificThermalEnergy)")
    output(nodesName + ".specificThermalEnergy[:]")

    output(nodesName + ".velocity")
    output("len(" + nodesName + ".velocity)")
    output(nodesName + ".haveField(" + nodesName + ".velocity)")
    output(nodesName + ".velocity[:]")

    output(nodesName + ".weight")
    output("len(" + nodesName + ".weight)")
    output(nodesName + ".haveField(" + nodesName + ".weight)")
    output(nodesName + ".weight[:]")

    output(nodesName + ".Hfield")
    output("len(" + nodesName + ".Hfield)")
    output(nodesName + ".haveField(" + nodesName + ".Hfield)")
    output(nodesName + ".Hfield[:]")

    output(nodesName + ".pressure")
    output("len(" + nodesName + ".pressure)")
    output(nodesName + ".haveField(" + nodesName + ".pressure)")
    output(nodesName + ".pressure[:]")

    global pmom
    pmom = nodes.linearMomentum
    output("pmom")
    output("len(pmom)")
    output(nodesName + ".haveField(pmom)")
    output("pmom[:]")

    nodes.mass[:nodes.numInternalNodes] = range(nodes.numInternalNodes)
    nodes.mass[nodes.numInternalNodes:] = 0.1*array(range(nodes.numGhostNodes))
    nodes.velocity[:nodes.numInternalNodes] = [(1, 2, 3)]*nodes.numInternalNodes
    nodes.velocity[nodes.numInternalNodes:] = [(-1, -2, -3)]*nodes.numGhostNodes
    pmom = nodes.linearMomentum
    output("pmom")
    output("len(pmom)")
    output(nodesName + ".haveField(pmom)")
    output("pmom[:]")

    nodes.massDensity[:nodes.numInternalNodes] = range(nodes.numInternalNodes)
    nodes.massDensity[nodes.numInternalNodes:] = 0.1*array(range(nodes.numGhostNodes))
    nodes.specificThermalEnergy[:nodes.numInternalNodes] = [10.0]*nodes.numInternalNodes
    nodes.specificThermalEnergy[nodes.numInternalNodes:] = [0.1]*nodes.numGhostNodes
    output(nodesName + ".massDensity[:]")
    output(nodesName + ".specificThermalEnergy[:]")
    output(nodesName + ".pressure[:]")
    global T
    T = nodes.temperature[:]
    output("T[:]")
    global cs
    cs = nodes.soundSpeed[:]
    output("cs[:]")

################################################################################
title("Testing Empty Sph Node list")
eos = GammaLawGasMKS3d(5.0/3.0, 2.0)
emptyNodes = SphNodeList3d(0, 0, eos)
TestNodes("emptyNodes")

################################################################################
title("Resizing node list to size 25 internal, 10 ghost")
emptyNodes.numInternalNodes = 25
emptyNodes.numGhostNodes = 10
TestNodes("emptyNodes")

################################################################################
title("Resizing node list to size 30 internal, 15 ghost")
emptyNodes.numInternalNodes = 30
emptyNodes.numGhostNodes = 15
TestNodes("emptyNodes")

################################################################################
title("Resizing node list to size 5 internal, 5 ghost")
emptyNodes.numInternalNodes = 5
emptyNodes.numGhostNodes = 5
TestNodes("emptyNodes")
