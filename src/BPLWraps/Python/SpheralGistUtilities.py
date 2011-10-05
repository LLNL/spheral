from math import *

################################################################################
# Read out the node positions into NumPy style arrays
################################################################################
def nodePositions1d(nodes):
    x = numpy.array([0.0]*nodes.numNodes)
    for i in xrange(nodes.numNodes):
        x[i] = nodes.positions[i].x
    return x

def nodePositions2d(nodes):
    x = numpy.array([0.0]*nodes.numNodes)
    y = numpy.array([0.0]*nodes.numNodes)
    for i in xrange(nodes.numNodes):
        x[i] = nodes.positions[i].x
        y[i] = nodes.positions[i].y
    return x, y

def nodePositions3d(nodes):
    x = numpy.array([0.0]*nodes.numNodes)
    y = numpy.array([0.0]*nodes.numNodes)
    z = numpy.array([0.0]*nodes.numNodes)
    for i in xrange(nodes.numNodes):
        x[i] = nodes.positions[i].x
        y[i] = nodes.positions[i].y
        z[i] = nodes.positions[i].z
    return x, y, z

################################################################################
# Convert a lattice of node positions to a gist-able mesh representation.
# This is not very SPH friendly, and won't generalize worth a damn!  We need
# something more capable of plotting randomly arranged data.
################################################################################
def nodePositions2Mesh2d(nodes, nx=None, ny=None):
    if not nx:
        nx = sqrt(nodes.numInternalNodes)
        ny = nx
    x = numpy.array([[0.0]*nx]*ny)
    y = numpy.array([[0.0]*nx]*ny)
    for iy in xrange(ny):
        for ix in xrange(nx):
            nodeID = ix + iy*nx
            x[iy][ix] = nodes.positions[nodeID].x
            y[iy][ix] = nodes.positions[nodeID].y
    return x, y

################################################################################
# Read out a lattice of vector field data for plotting on a gist style mesh.
################################################################################
def vectorField2NumericArrays2d(field, nx=None, ny=None):
    if not nx:
        nx = sqrt(len(field.numElements))
        ny = nx
    x = numpy.array([[0.0]*nx]*ny)
    y = numpy.array([[0.0]*nx]*ny)
    for iy in xrange(ny):
        for ix in xrange(nx):
            nodeID = ix + iy*nx
            x[iy][ix] = field[nodeID].x
            y[iy][ix] = field[nodeID].y
    return x, y

################################################################################
# Generate a radial profile of a field.
################################################################################
def radialProfile(field):
    r = numpy.array([0.0]*len(field))
    fieldr = numpy.array([field[0]]*len(field))
    nodeList = field.nodeList
    for nodeID in xrange(nodeList.numNodes):
        r[nodeID] = nodeList.positions[nodeID].magnitude()
        fieldr[nodeID] = field[nodeID]
    return r, fieldr

################################################################################
# Plot the given node indicies
################################################################################
def plotNodes2d(nodes, index=None,
                title=None,
                windowNum=0, title=None, marker='\4', color='black'):
    try:
        from gist import *
        from numpy import *
        x, y = nodePositions2d(nodes)
        if index:
            xplot = array([0.0]*len(index))
            yplot = array([0.0]*len(index))
            for i in xrange(len(index)):
                xplot[i] = x[index[i]]
                yplot[i] = y[index[i]]
        else:
            xplot = x
            yplot = y
        window(windowNum)
        plg(yplot, xplot, type=0, marker=marker, color=color)
        if title:
            pltitle(title)
    except:
        print "Skipping"
        pass
    return

################################################################################
# Plot the mass density, velocity, pressure, and smoothing scale for the fluid
# node lists in the given data base.  Implicitly assuming 1-D.
################################################################################
def plotState(dataBase, color='black',
              plotGhosts=0,
              logM=0,
              logV=0,
              logP=0,
              logEps=0,
              logH=0):
    from SpheralGistUtilities import *
    from gist import *
    from math import *

    for nodes in dataBase.fluidNodeLists:

        if plotGhosts:
            nx = nodes.numNodes
        else:
            nx = nodes.numInternalNodes

        window(0)
        xNodes = nodePositions1d(nodes)[:nx]
        rhoNodes = array(nodes.massDensity[:nx])
        if logM:
            rhoNodes = log10(rhoNodes)
        plg(rhoNodes, xNodes, color=color)
        pltitle('Mass density')

        window(1)
        vNodes = array([0.0]*nx)
        for i in xrange(nx):
            vNodes[i] = nodes.velocity[i].x
        if logV:
            vNodes = log10(vNodes)
        plg(vNodes, xNodes, color=color)
        pltitle('Velocity')

        window(2)
        pressure = nodes.pressure
        PNodes = array(pressure[:nx])
        if logP:
            PNodes = log10(PNodes)
        plg(PNodes, xNodes, color=color)
        pltitle('Pressure')

        window(3)
        eps = nodes.specificThermalEnergy
        epsNodes = array(eps[:nx])
        if logEps:
            epsNodes = log10(epsNodes)
        plg(epsNodes, xNodes, color=color)
        pltitle('Specific thermal energy')

        window(4)
        HNodes = array([0.0]*nx)
        for i in xrange(nx):
            HNodes[i] = 1.0/nodes.Hfield[i].xx
        if logH:
            HNodes = log10(HNodes)
        plg(HNodes, xNodes, color=color)
        pltitle('Smoothing scale')

################################################################################
# Overplot the answer on results from plotState.
################################################################################
def plotAnswer(answerObject, time, color='black',
               logM=0,
               logV=0,
               logP=0,
               logEps=0,
               logH=0):
    from SpheralGistUtilities import *
    from gist import *
    from math import *

    x, v, u, rho, P, h = answerObject.solution(time)

    window(0)
    if logM:
        plg(log10(array(rho)), array(x), color=color)
    else:
        plg(array(rho), array(x), color=color)

    window(1)
    if logV:
        plg(log10(array(v)), array(x), color=color)
    else:
        plg(array(v), array(x), color=color)

    window(2)
    if logP:
        plg(log10(array(P)), array(x), color=color)
    else:
        plg(array(P), array(x), color=color)

    window(3)
    if logEps:
        plg(log10(array(u)), array(x), color=color)
    else:
        plg(array(u), array(x), color=color)

    window(4)
    if logH:
        plg(log10(array(h)), array(x), color=color)
    else:
        plg(array(h), array(x), color=color)

