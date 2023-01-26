import Gnuplot
import mpi
from Spheral import *
from math import *
import numpy
import os
from SpheralTestUtilities import multiSort

SpheralGnuPlotCache = []

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Define a dummy Gnuplot class, so that non-master processes can silently
# and harmlessly accept Gnuplot commands.
#-------------------------------------------------------------------------------
class fakeGnuplot:
    def __init__(self):
        return
    def __call__(self, *arghs, **keyw):
        return
    def plot(self, *arghs, **keyw):
        return
    def replot(self, *arghs, **keyw):
        return
    def refresh(self, *arghs, **keyw):
        return
    def xlabel(self, *arghs, **keyw):
        return
    def ylabel(self, *arghs, **keyw):
        return
    def title(self, *arghs, **keyw):
        return
    def hardcopy(self, *arghs, **keyw):
        return

def generateNewGnuPlot(persist = False):
    if mpi.rank == 0:
        result = Gnuplot.Gnuplot(persist = persist)
        if "GNUTERM" in list(os.environ.keys()):
            result("set term %s" % os.environ["GNUTERM"])
        return result
    else:
        return fakeGnuplot()

#-------------------------------------------------------------------------------
# Since the default Gnuplot.py doesn't support png output, I'll add it here
# myself.
#-------------------------------------------------------------------------------
def pngFile(plot, filename,
            color = 1,
            fontSize = "medium"):
    setLine = "set terminal png " + fontSize
    if color:
        setLine += " color"
    if filename[-4:] != ".png":
        filename += ".png"
    plot(setLine)
    plot.set_string("output", filename)
    plot.refresh()
    plot("set terminal x11")
    plot.set_string("output")
    return

#-------------------------------------------------------------------------------
# Calculate the radial velocity component, given a FieldList of positions
# and a FieldList of velocities.
#-------------------------------------------------------------------------------
def radialVelocityFieldList(positions,
                            velocities):

    dim = type(positions).__name__[-2:]
    radialVelocity = None
    fieldConstructor = None
    if dim == "1d":
        radialVelocity = ScalarFieldList1d()
        fieldConstructor = ScalarField1d
    elif dim == "2d":
        radialVelocity = ScalarFieldList2d()
        fieldConstructor = ScalarField2d
    elif dim == "3d":
        radialVelocity = ScalarFieldList3d()
        fieldConstructor = ScalarField3d
    radialVelocity.copyFields()
    for field in positions:
        radialVelocity.appendField(fieldConstructor("radial velocity", field.nodeList()))

    assert positions.numFields == velocities.numFields == radialVelocity.numFields

    for fieldID in range(positions.numFields):
        rfield = positions[fieldID]
        vfield = velocities[fieldID]
        vrfield = radialVelocity[fieldID]
        assert rfield.numElements == vfield.numElements == vrfield.numElements

        for nodeID in range(rfield.numElements):
            r = rfield[nodeID]
            v = vfield[nodeID]
            runit = r.unitVector()
            vrfield[nodeID] = v.dot(runit)

    return radialVelocity

#-------------------------------------------------------------------------------
# Calculate the azimuthal velocity component, given a FieldList of positions
# and a FieldList of velocities.
#-------------------------------------------------------------------------------
def azimuthalVelocityFieldList(positions,
                               velocities):

    dim = type(positions).__name__[-2:]
    azimuthalVelocity = None
    fieldConstructor = None
    if dim == "1d":
        azimuthalVelocity = ScalarFieldList1d()
        fieldConstructor = ScalarField1d
    elif dim == "2d":
        azimuthalVelocity = ScalarFieldList2d()
        fieldConstructor = ScalarField2d
    elif dim == "3d":
        azimuthalVelocity = ScalarFieldList3d()
        fieldConstructor = ScalarField3d
    azimuthalVelocity.copyFields()
    for field in positions:
        azimuthalVelocity.appendField(fieldConstructor("azimuthal velocity", field.nodeList()))

    assert positions.numFields == velocities.numFields == azimuthalVelocity.numFields

    for fieldID in range(positions.numFields):
        rfield = positions[fieldID]
        vfield = velocities[fieldID]
        vafield = azimuthalVelocity[fieldID]
        assert rfield.numElements == vfield.numElements == vafield.numElements

        for nodeID in range(rfield.numElements):
            r = rfield[nodeID]
            v = vfield[nodeID]
            raz = r.unitVector()
            x = raz.x
            y = raz.y
            raz.x = -y
            raz.y = x
            vafield[nodeID] = v.dot(raz)

    return azimuthalVelocity

#-------------------------------------------------------------------------------
# Helper method to determine the angular momentum per node.
#-------------------------------------------------------------------------------
def angularMomentum(mass, position, velocity):
    assert mass.numFields == position.numFields == velocity.numFields

    result = []
    for massField, positionField, velocityField in zip(mass,
                                                       position,
                                                       velocity):
        assert (massField.nodeList().numInternalNodes ==
                positionField.nodeList().numInternalNodes ==
                velocityField.nodeList().numInternalNodes)
        for j in range(massField.nodeList().numInternalNodes):
            result.append((positionField[j].cross(velocityField[j]))*massField[j])

    return result

#-------------------------------------------------------------------------------
# Plot a FieldList
#-------------------------------------------------------------------------------
def plotFieldList(fieldList,
                  xFunction = "%s.x",
                  yFunction = "%s",
                  plotGhosts = False,
                  colorNodeLists = False,
                  plot = None,
                  userXRange = [None, None],
                  userYRange = [None, None],
                  plotStyle = "lines",
                  lineStyle = "linetype -1 linewidth 1 pointtype 4 pointsize 1.0",
                  winTitle = None,
                  lineTitle = "",
                  xlabel = None,
                  ylabel = None,
                  filterFunc = None):

    if plot is None:
        plot = generateNewGnuPlot()
    SpheralGnuPlotCache.append(plot)

    def nullFilter(pos):
        return True

    if filterFunc is None:
        filterFunc = nullFilter

    # Gather the fieldList info across all processors to process 0.
    globalNumNodes = []
    globalX = []
    globalY = []
    for field in fieldList:
        if plotGhosts:
            xvals = field.nodeList().positions().allValues()
            yvals = field.allValues()
        else:
            xvals = field.nodeList().positions().internalValues()
            yvals = field.internalValues()
        localX = []
        localY = []
        for x, y in zip(xvals, yvals):
            if filterFunc(x):
                localX.append(eval(xFunction % "x"))
                localY.append(eval(yFunction % "y"))
        n = len(localX)
        if mpi:
            globalNumNodes.append(mpi.allreduce(n, mpi.SUM))
            globalX.extend(mpi.allreduce(localX, mpi.SUM))
            globalY.extend(mpi.allreduce(localY, mpi.SUM))
        else:
            globalNumNodes.append(n)
            globalX.extend(localX)
            globalY.extend(localY)

    if mpi.rank == 0:
        # Find the total number of nodes.
        totalNumNodes = sum(globalNumNodes)
        assert(len(globalNumNodes) == fieldList.numFields)
        assert(len(globalX) == totalNumNodes)
        assert(len(globalY) == totalNumNodes)

        # Copy the input ranges, since for some reason these seem to have been
        # preserved between calls?
        xRange = userXRange[:]
        yRange = userYRange[:]

        # Set the line style
##         plot("set linestyle 1 " + lineStyle)

        # Set the labels.
        if winTitle: plot.title(winTitle)
        if xlabel: plot.xlabel(xlabel)
        if ylabel: plot.ylabel(ylabel)

        # Set the ranges.
        xmin = 1e30
        xmax = -1e30
        ymin = 1e30
        ymax = -1e30
        for x in globalX:
            xmin = min(xmin, x)
            xmax = max(xmax, x)
        for y in globalY:
            ymin = min(ymin, y)
            ymax = max(ymax, y)
        if xmin == xmax:
            xmin = xmin - 0.5
            xmax = xmax + 0.5
        if ymin == ymax:
            ymin = ymin - 0.5
            ymax = ymax + 0.5
        if xRange[0] == None: xRange[0] = xmin
        if xRange[1] == None: xRange[1] = xmax
        if yRange[0] == None: yRange[0] = ymin - 0.05*max(1e-5, ymax - ymin)
        if yRange[1] == None: yRange[1] = ymax + 0.05*max(1e-5, ymax - ymin)
        plot("set xrange [%f:%f]" % tuple(xRange))
        plot("set yrange [%f:%f]" % tuple(yRange))

        # Finally, loop over the fields and do the deed.
        assert(len(globalX) == len(globalY))
        if colorNodeLists:
            legendNodeList = {}
            for i in range(fieldList.numFields):
                legendNodeList[i] = lineTitle + ": " + fieldList[i].nodeList().name

            cumulativeNumNodes = 0
            for fieldID in range(len(globalNumNodes)):
                n = globalNumNodes[fieldID]
                iNodeList = fieldID % fieldList.numFields
                x = numpy.array(globalX[cumulativeNumNodes:
                                        cumulativeNumNodes + n])
                y = numpy.array(globalY[cumulativeNumNodes:
                                        cumulativeNumNodes + n])
                if n:
##                    plot("set linestyle %i lt %i pt %i" % (iNodeList + 1,
##                                                           iNodeList + 1,
##                                                           iNodeList + 1))
                    legend = legendNodeList[iNodeList]
                    legendNodeList[iNodeList] = None
                    data = Gnuplot.Data(x, y,
                                        with_ = plotStyle + " lt %i" % iNodeList,
                                        title = legend,
                                        inline = True)
                    plot.replot(data)
                    SpheralGnuPlotCache.append(data)

                    cumulativeNumNodes += n
                
        else:
            x = numpy.array(globalX)
            y = numpy.array(globalY)
            data = Gnuplot.Data(x, y,
                                with_ = plotStyle + " lt -1 pt 3",
                                title = lineTitle,
                                inline = True)
            plot.replot(data)
            SpheralGnuPlotCache.append(data)
            lineTitle = None

    # That's it, return the Gnuplot object.
    mpi.barrier()
    return plot

#-------------------------------------------------------------------------------
# Plot the mass density, velocity, pressure, and smoothing scale for the fluid
# node lists in the given data base.  Implicitly assuming 1-D.
#-------------------------------------------------------------------------------
def plotState(thingus,
              plotGhosts = False,
              colorNodeLists = False,
              plotStyle = "points",
              xFunction = "%s.x",
              vecyFunction = "%s.x",
              tenyFunction = "%s.xx ** -1",
              lineTitle = "Simulation",
              filterFunc = None):

    dim = type(thingus).__name__[-2:]
    if isinstance(thingus, eval("State%s" % dim)):
        rho = thingus.scalarFields(HydroFieldNames.massDensity)
        vel = thingus.vectorFields(HydroFieldNames.velocity)
        eps = thingus.scalarFields(HydroFieldNames.specificThermalEnergy)
        P = thingus.scalarFields(HydroFieldNames.pressure)
        H = thingus.symTensorFields(HydroFieldNames.H)

    else:
        assert isinstance(thingus, eval("DataBase%s" % dim))
        rho = thingus.fluidMassDensity
        vel = thingus.fluidVelocity
        eps = thingus.fluidSpecificThermalEnergy
        P = thingus.newFluidScalarFieldList(0.0, "pressure")
        thingus.fluidPressure(P)
        H = thingus.fluidHfield

    rhoPlot = plotFieldList(rho,
                            xFunction = xFunction,
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = plotStyle,
                            winTitle = "Mass Density",
                            lineTitle = lineTitle,
                            xlabel="x",
                            filterFunc = filterFunc)

    velPlot = plotFieldList(vel,
                            xFunction = xFunction,
                            yFunction = vecyFunction,
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = plotStyle,
                            winTitle = "Velocity",
                            lineTitle = lineTitle,
                            xlabel="x",
                            filterFunc = filterFunc)

    epsPlot = plotFieldList(eps,
                            xFunction = xFunction,
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = plotStyle,
                            winTitle = "Specific Thermal Energy",
                            lineTitle = lineTitle,
                            xlabel="x",
                            filterFunc = filterFunc)

    PPlot = plotFieldList(P,
                          xFunction = xFunction,
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          plotStyle = plotStyle,
                          winTitle = "Pressure",
                          lineTitle = lineTitle,
                          xlabel="x",
                          filterFunc = filterFunc)

    HPlot = plotFieldList(H,
                          xFunction = xFunction,
                          yFunction = tenyFunction,
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          plotStyle = plotStyle,
                          winTitle = "Smoothing scale",
                          lineTitle = lineTitle,
                          xlabel="x",
                          filterFunc = filterFunc)

    return rhoPlot, velPlot, epsPlot, PPlot, HPlot

#-------------------------------------------------------------------------------
# Plot the state vs. radius
#-------------------------------------------------------------------------------
def plotRadialState(dataBase,
                    plotGhosts = False,
                    colorNodeLists = False,
                    lineTitle = "Simulation",
                    filterFunc = None):

    rhoPlot = plotFieldList(dataBase.fluidMassDensity,
                            xFunction = "%s.magnitude()",
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = "points",
                            winTitle = "Mass density",
                            lineTitle = lineTitle,
                            xlabel = "r",
                            filterFunc = filterFunc)

    radialVelocity = radialVelocityFieldList(dataBase.fluidPosition,
                                             dataBase.fluidVelocity)
    velPlot = plotFieldList(radialVelocity,
                            xFunction = "%s.magnitude()",
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = "points",
                            winTitle = " Radial Velocity",
                            lineTitle = lineTitle,
                            xlabel = "r",
                            filterFunc = filterFunc)

    epsPlot = plotFieldList(dataBase.fluidSpecificThermalEnergy,
                            xFunction = "%s.magnitude()",
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = "points",
                            winTitle = "Specific Thermal Energy",
                            lineTitle = lineTitle,
                            xlabel = "r",
                            filterFunc = filterFunc)

    fluidPressure = dataBase.newFluidScalarFieldList(0.0, "pressure")
    dataBase.fluidPressure(fluidPressure)
    PPlot = plotFieldList(fluidPressure,
                          xFunction = "%s.magnitude()",
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          plotStyle = "points",
                          winTitle = "Pressure",
                          lineTitle = lineTitle,
                          xlabel = "r",
                          filterFunc = filterFunc)

    HPlot = plotFieldList(dataBase.fluidHfield,
                          xFunction = "%s.magnitude()",
                          yFunction = "%s.xx**-1",
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          plotStyle = "points",
                          winTitle = "Smoothing scale",
                          lineTitle = lineTitle,
                          xlabel = "r",
                          filterFunc = filterFunc)

    return rhoPlot, velPlot, epsPlot, PPlot, HPlot

#-------------------------------------------------------------------------------
# Overplot the answer on results from plotState.
#-------------------------------------------------------------------------------
def plotAnswer(answerObject, time,
               rhoPlot = None,
               velPlot = None,
               epsPlot = None,
               PPlot = None,
               APlot = None,
               HPlot = None,
               x = None):

    try:
        x, v, u, rho, P, h = answerObject.solution(time, x)
        A = None
    except:
        try:
            x, v, u, rho, P, A, h = answerObject.solution(time, x)
        except:
            x, v, u, rho, P = answerObject.solution(time, x)
            A = None
            h = None

    if rhoPlot is not None:
        data = Gnuplot.Data(x, rho,
                            with_="lines lt 7 lw 2",
                            title="Solution",
                            inline = True)
        SpheralGnuPlotCache.append(data)
        rhoPlot.replot(data)

    if velPlot is not None:
        data = Gnuplot.Data(x, v,
                            with_="lines lt 7 lw 2",
                            title="Solution",
                            inline = True)
        SpheralGnuPlotCache.append(data)
        velPlot.replot(data)

    if epsPlot is not None:
        data = Gnuplot.Data(x, u,
                            with_="lines lt 7 lw 2",
                            title="Solution",
                            inline = True)
        SpheralGnuPlotCache.append(data)
        epsPlot.replot(data)

    if PPlot is not None:
        data = Gnuplot.Data(x, P,
                            with_="lines lt 7 lw 2",
                            title="Solution",
                            inline = True)
        SpheralGnuPlotCache.append(data)
        PPlot.replot(data)

    if APlot is not None and A:
        data = Gnuplot.Data(x, A,
                            with_="lines lt 7 lw 2",
                            title="Solution",
                            inline = True)
        SpheralGnuPlotCache.append(data)
        APlot.replot(data)

    if HPlot is not None:
        data = Gnuplot.Data(x, h,
                            with_="lines lt 7 lw 2",
                            title="Solution",
                            inline = True)
        SpheralGnuPlotCache.append(data)
        HPlot.replot(data)

    return

#-------------------------------------------------------------------------------
# Plot the node positions
#-------------------------------------------------------------------------------
def plotNodePositions2d(thingy,
                        xFunction = "%s.x",
                        yFunction = "%s.y",
                        plotGhosts = False,
                        colorNodeLists = True,
                        colorDomains = False,
                        title = "",
                        style = "points",
                        persist = None):

    assert colorNodeLists + colorDomains <= 1

    if isinstance(thingy, DataBase2d):
        nodeLists = thingy.nodeLists()
    else:
        nodeLists = thingy

    # Gather the node positions across all domains.
    # Loop over all the NodeLists.
    xNodes = []
    yNodes = []
    for nodeList in nodeLists:
        if plotGhosts:
            pos = nodeList.positions().allValues()
        else:
            pos = nodeList.positions().internalValues()
        xNodes.append([eval(xFunction % "x") for x in pos])
        yNodes.append([eval(yFunction % "x") for x in pos])
    assert len(xNodes) == len(nodeLists)
    assert len(xNodes) == len(yNodes)
    
    globalXNodes = mpi.gather(xNodes)
    globalYNodes = mpi.gather(yNodes)

    if mpi.rank == 0:
        assert len(globalXNodes) == mpi.procs
        assert len(globalYNodes) == mpi.procs

        xlist, ylist = [], []
        if colorDomains:
            for xDomain, yDomain in zip(globalXNodes, globalYNodes):
                assert len(xDomain) == len(nodeLists)
                assert len(yDomain) == len(nodeLists)
                xlist.append([])
                ylist.append([])
                for xx in xDomain:
                    xlist[-1].extend(xx)
                for yy in yDomain:
                    ylist[-1].extend(yy)
            assert len(xlist) == mpi.procs
            assert len(ylist) == mpi.procs

        elif colorNodeLists:
            for i in range(len(nodeLists)):
                xlist.append([])
                ylist.append([])
            for xDomain, yDomain in zip(globalXNodes, globalYNodes):
                assert len(xDomain) == len(nodeLists)
                assert len(yDomain) == len(nodeLists)
                for i in range(len(nodeLists)):
                    xlist[i].extend(xDomain[i])
                    ylist[i].extend(yDomain[i])
            assert len(xlist) == len(nodeLists)
            assert len(ylist) == len(nodeLists)

        else:
            xlist, ylist = [[]], [[]]
            for xDomain, yDomain in zip(globalXNodes, globalYNodes):
                print(len(xDomain), len(nodeLists))
                assert len(xDomain) == len(nodeLists)
                assert len(yDomain) == len(nodeLists)
                for i in range(len(nodeLists)):
                    xlist[0].extend(xDomain[i])
                    ylist[0].extend(yDomain[i])

        plot = generateNewGnuPlot(persist = persist)
        plot("set size square")
        plot.title = title
        assert len(xlist) == len(ylist)
        for x, y in zip(xlist, ylist):
            data = Gnuplot.Data(x, y, 
                                with_ = style,
                                inline = True)
            plot.replot(data)
            SpheralGnuPlotCache.append(data)

        return plot

    else:
        return fakeGnuplot()

#-------------------------------------------------------------------------------
# Plot all the nodes in the given data base, and then color the control/ghost
# nodes of the given boundary condition independently.
#-------------------------------------------------------------------------------
def plotBoundaryNodes(dataBase, boundary):

    # First build one set of position pairs for all of the nodes in the
    # data base.
    positions = []
    for nodeList in dataBase.nodeLists():
        for r in list(nodeList.positions())[:nodeList.numInternalNodes]:
            positions.append((r.x, r.y))

    # Now build a list of the control node positions from the boundary
    # condition.
    controlPositions = []
    for nodeList in dataBase.nodeLists():
        controlNodes = boundary.controlNodes(nodeList)
        for nodeID in controlNodes:
            r = nodeList.positions()[nodeID]
            controlPositions.append((r.x, r.y))

    # Now build a list of the ghost node positions from the boundary
    # condition.
    ghostPositions = []
    for nodeList in dataBase.nodeLists():
        ghostNodes = boundary.ghostNodes(nodeList)
        for nodeID in ghostNodes:
            r = nodeList.positions()[nodeID]
            ghostPositions.append((r.x, r.y))

    # Finally we can plot these various sets of nodes.
    plot = plotXYTuples([positions, controlPositions, ghostPositions])
    return plot

#-------------------------------------------------------------------------------
# Plot the given sequences of (x,y) pairs, each with a distinct color.
#  [ [(x0,y0), (x1,y1), ...],
#    [(x0,y0), (x1,y1), ...],
#    .
#    .
#    .
#    [(x0,y0), (x1,y1), ...] ]
#-------------------------------------------------------------------------------
def plotXYTuples(listOfXYTuples):

    # Find the (min,max) of X and Y for all sets.
    xmin, ymin, xmax, ymax = findPairMinMax(listOfXYTuples[0])
    for seq in listOfXYTuples[1:]:
        xmin0, ymin0, xmax0, ymax0 = findPairMinMax(seq)
        xmin = min(xmin, xmin0)
        ymin = min(ymin, ymin0)
        xmax = max(xmax, xmax0)
        ymax = max(ymax, ymax0)

    # Create our plot result.
    plot = generateNewGnuPlot()
    plot("set size square")

    # Loop over the list of sequences of positions.
    icolor = 0
    for seq in listOfXYTuples:
        icolor += 1

        # Build the local arrays of x and y.
        x = numpy.array([0.0]*len(seq))
        y = numpy.array([0.0]*len(seq))
        for i in range(len(seq)):
            x[i] = seq[i][0]
            y[i] = seq[i][1]

        # Build the gnuplot data.
        data = Gnuplot.Data(x, y,
                            with_ = "points",
                            inline = True)
        SpheralGnuPlotCache.append(data)

        # Plot this set of data.
##        plot("set linestyle %i lt %i pt 1" % (icolor, icolor))
        plot.replot(data)

    # That"s it, return the plot.
    return plot

#-------------------------------------------------------------------------------
# Find the (min, max) of a set of pairs.
#-------------------------------------------------------------------------------
def findPairMinMax(listOfPairs):
    minX, minY = 1e90, 1e90
    maxX, maxY = -1e90, -1e90
    for pair in listOfPairs:
        minX = min(minX, pair[0])
        minY = min(minY, pair[1])
        maxX = max(maxX, pair[0])
        maxY = max(maxY, pair[1])
    return minX, minY, maxX, maxY

#-------------------------------------------------------------------------------
# Plot the velocity field as a set of arrows.
# This is maintained here for backward compatibility, as a specialization of
# plotVectorField2d.
#-------------------------------------------------------------------------------
def plotVelocityField2d(dataBase,
                        plotGhosts = False,
                        velMultiplier = 1.0,
                        colorNodeLists = False,
                        colorDomains = False,
                        title = ""):

    return plotVectorField2d(dataBase,
                             dataBase.globalVelocity,
                             plotGhosts,
                             velMultiplier,
                             colorNodeLists,
                             colorDomains,
                             title)

#-------------------------------------------------------------------------------
# Plot the node spacing in 1D.
#-------------------------------------------------------------------------------
def plotNodeSpacing1d(dataBase):
    pos = dataBase.globalPosition
    xvals = []
    for ifield in range(len(pos)):
        xvals += [pos[ifield][i].x for i in range(pos[ifield].numInternalElements)]
    xvals = mpi.allreduce(xvals, mpi.SUM)
    xvals.sort()
    deltas = [xvals[i+1] - xvals[i] for i in range(len(xvals) - 1)] + [xvals[-1] - xvals[-2]]
    plot = generateNewGnuPlot()
    d = Gnuplot.Data(xvals, deltas, with_="lines")
    plot.plot(d)
    return plot

#-------------------------------------------------------------------------------
# Plot an arbitrary vector field as a set of arrows.
#-------------------------------------------------------------------------------
def plotVectorField2d(dataBase, fieldList,
                      plotGhosts = False,
                      vectorMultiplier = 1.0,
                      colorNodeLists = False,
                      colorDomains = False,
                      title = ""):

    assert colorNodeLists + colorDomains <= 1

    # Gather the node positions and vectors across all domains.
    # Loop over all the NodeLists.
    localNumNodes = []
    xNodes = []
    yNodes = []
    vxNodes = []
    vyNodes = []
    for i in range(dataBase.numNodeLists):
        nodeList = dataBase.nodeLists()[i]
        assert i < fieldList.numFields
        vectorField = fieldList[i]
        if plotGhosts:
            n = nodeList.numNodes
        else:
            n = nodeList.numInternalNodes
        localNumNodes.append(n)
        xNodes += numpy.array([x.x for x in list(nodeList.positions())[:n]])
        yNodes += numpy.array([x.y for x in list(nodeList.positions())[:n]])
        vxNodes += numpy.array([x.x for x in list(vectorField)[:n]])*vectorMultiplier
        vyNodes += numpy.array([x.y for x in list(vectorField)[:n]])*vectorMultiplier
    assert len(xNodes) == len(yNodes) == len(vxNodes) == len(vyNodes)
    
    numDomainNodes = [len(xNodes)]
    numNodesPerDomain = mpi.gather(numDomainNodes)
    globalNumNodes = mpi.gather(localNumNodes)
    globalXNodes = mpi.gather(xNodes)
    globalYNodes = mpi.gather(yNodes)
    globalVxNodes = mpi.gather(vxNodes)
    globalVyNodes = mpi.gather(vyNodes)

    if mpi.rank == 0:
        plot = generateNewGnuPlot()
        plot("set size square")
        plot.title = title

        if colorDomains:
            cumulativeN = 0
            for domain in range(len(numNodesPerDomain)):
                n = numNodesPerDomain[domain]
                x = numpy.array(globalXNodes[cumulativeN:cumulativeN + n])
                y = numpy.array(globalYNodes[cumulativeN:cumulativeN + n])
                vx = numpy.array(globalVxNodes[cumulativeN:cumulativeN + n])
                vy = numpy.array(globalVyNodes[cumulativeN:cumulativeN + n])
                cumulativeN += n
##                plot("set linestyle %i lt %i pt %i" % (domain + 1,
##                                                       domain + 1,
##                                                       domain + 1))
                data = Gnuplot.Data(x, y, vx, vy,
                                    with_ = "vector ls %i" % (domain + 1),
                                    inline = True)
                plot.replot(data)
                SpheralGnuPlotCache.append(data)

        elif colorNodeLists:
            cumulativeN = 0
            for i in range(len(globalNumNodes)):
                n = globalNumNodes[i]
                if n > 0:
                    iNodeList = i % dataBase.numNodeLists
                    x = numpy.array(globalXNodes[cumulativeN:cumulativeN + n])
                    y = numpy.array(globalYNodes[cumulativeN:cumulativeN + n])
                    vx = numpy.array(globalVxNodes[cumulativeN:cumulativeN + n])
                    vy = numpy.array(globalVyNodes[cumulativeN:cumulativeN + n])
                    cumulativeN += n
##                    plot("set linestyle %i lt %i pt %i" % (iNodeList + 1,
##                                                           iNodeList + 1,
##                                                           iNodeList + 1))
                    data = Gnuplot.Data(x, y, vx, vy,
                                        with_ = "vector ls %i" % (iNodeList + 1),
                                        inline = True)
                    plot.replot(data)
                    SpheralGnuPlotCache.append(data)

        else:
            x = numpy.array(globalXNodes)
            y = numpy.array(globalYNodes)
            vx = numpy.array(globalVxNodes)
            vy = numpy.array(globalVyNodes)
            data = Gnuplot.Data(x, y, vx, vy,
                                with_ = "vector",
                                inline = True)
            plot.replot(data)
            SpheralGnuPlotCache.append(data)

        return plot

    else:
        SpheralGnuPlotCache.append(data)

#-------------------------------------------------------------------------------
# Generate a regularly spaced sampling of the given FieldList
# The answer is returned in a 2-D numpy array.
#-------------------------------------------------------------------------------
def gridSample(fieldList,
               zFunction = "%s",
               nx = 100,
               ny = 100,
               xmin = None,
               xmax = None,
               ymin = None,
               ymax = None):

    assert nx > 0 and ny > 0

    # Set up our return value array.
    xValues = numpy.array([[0.0]*nx]*ny)
    yValues = numpy.array([[0.0]*nx]*ny)
    zValues = numpy.array([[0.0]*nx]*ny)

    # Gather the fieldList info across all processors to process 0.
    localNumNodes = []
    localX = []
    localY = []
    for ifield in range(fieldList.numFields):
        field = fieldList[ifield]
        n = field.nodeList().numNodes
        localNumNodes.append(n)
        for r in field.nodeList().positions():
            localX.append(r.x)
            localY.append(r.y)
    globalNumNodes = mpi.gather(localNumNodes)
    globalX = mpi.gather(localX)
    globalY = mpi.gather(localY)

    # If the user did not specify the sampling volume, then find the min and
    # max node positions.
    if xmin == None:
        xmin = min(localX)
    if ymin == None:
        ymin = min(localY)
    if xmax == None:
        xmax = max(localX)
    if ymax == None:
        ymax = max(localY)
    xmin = mpi.allreduce(xmin, mpi.MIN)
    ymin = mpi.allreduce(ymin, mpi.MIN)
    xmax = mpi.allreduce(xmax, mpi.MAX)
    ymax = mpi.allreduce(ymax, mpi.MAX)

    assert xmax > xmin
    assert ymax > ymin

    # Figure out the sizes of the bins we're going to be sampling in
    dx = (xmax - xmin)/nx
    dy = (ymax - ymin)/ny

    # Loop over all the grid sampling positions, and figure out this processors
    # contribution.
    for iy in range(ny):
        for ix in range(nx):
            xValues[iy][ix] = xmin + (ix + 0.5)*dx
            yValues[iy][ix] = ymin + (iy + 0.5)*dy
            r = Vector2d(xValues[iy][ix], yValues[iy][ix])
            z = fieldList.sample(r)
            localZ = eval(zFunction % "z")
            globalZ = mpi.reduce(localZ, mpi.SUM)
            if mpi.rank == 0:
                print("%i %i %i %s %g %g" % (mpi.rank, ix, iy, r, z, localZ))
                print("%i %g" % (mpi.rank, globalZ))
                zValues[iy][ix] = globalZ

    return xValues, yValues, zValues

#-------------------------------------------------------------------------------
# Plot the energy history of the given conservation object.
#-------------------------------------------------------------------------------
def plotEHistory(conserve):
    if mpi.rank == 0:
        t = conserve.timeHistory
        E = conserve.EHistory
        KE = conserve.KEHistory
        TE = conserve.TEHistory
        UE = conserve.EEHistory
        Edata = Gnuplot.Data(t, E,
                             with_ = "lines",
                             title = "Total Energy",
                             inline = True)
        KEdata = Gnuplot.Data(t, KE,
                              with_ = "lines",
                              title = "Kinetic Energy",
                              inline = True)
        TEdata = Gnuplot.Data(t, TE,
                              with_ = "lines",
                              title = "Thermal Energy",
                              inline = True)
        UEdata = Gnuplot.Data(t, UE,
                              with_ = "lines",
                              title = "Potential Energy",
                              inline = True)
        plot = generateNewGnuPlot()
        plot.replot(Edata)
        plot.replot(KEdata)
        plot.replot(TEdata)
        plot.replot(UEdata)
        plot.replot()
        SpheralGnuPlotCache.extend([Edata, KEdata, TEdata, UEdata])
        return plot
    else:
        return fakeGnuplot()

#-------------------------------------------------------------------------------
# Plot the linear momentum history of the given conservation object.
#-------------------------------------------------------------------------------
def plotpmomHistory(conserve):
    if mpi.rank == 0:
        t = conserve.timeHistory
        p = conserve.pmomHistory
        px = [x.x for x in p]
        py = [x.y for x in p]
        pz = [x.z for x in p]
        pmag = [x.magnitude() for x in p]
        pxdata = Gnuplot.Data(t, px,
                              with_ = "lines",
                              title = "x momentum",
                              inline = True)
        pydata = Gnuplot.Data(t, py,
                              with_ = "lines",
                              title = "y momentum ",
                              inline = True)
        pzdata = Gnuplot.Data(t, pz,
                              with_ = "lines",
                              title = "z momentum",
                              inline = True)
        pmagdata = Gnuplot.Data(t, pmag,
                                with_ = "lines",
                                title = "total momentum",
                                inline = True)
        plot = generateNewGnuPlot()
        plot.replot(pxdata)
        plot.replot(pydata)
        plot.replot(pzdata)
        plot.replot(pmagdata)
        plot.replot()
        SpheralGnuPlotCache.extend([pxdata, pydata, pzdata, pmagdata])
        return plot
    else:
        return fakeGnuplot()

#-------------------------------------------------------------------------------
# Plot a polygon.
#-------------------------------------------------------------------------------
def plotPolygon(polygon,
                plotVertices = True,
                plotFacets = True,
                plotNormals = False,
                plotCentroid = False,
                plot = None,
                persist = False,
                plotLabels = True):
    px = []
    py = []
    for v in polygon.vertices:
        px.append(v.x)
        py.append(v.y)
    fx = []
    fy = []
    fdx = []
    fdy = []
    nx = []
    ny = []
    ndx = []
    ndy = []
    for f in polygon.facets:
        dr = f.point2 - f.point1
        hdr = dr/2.0
        fx.append(f.point1.x)
        fy.append(f.point1.y)
        fdx.append(dr.x)
        fdy.append(dr.y)
        nx.append(fx[-1] + hdr.x)
        ny.append(fy[-1] + hdr.y)
        ndx.append(f.normal.x)
        ndy.append(f.normal.y)
    if plot is None:
        plot = generateNewGnuPlot(persist)
    if plotLabels:
        vlabel, flabel, nlabel = "Vertices", "Facets", "Normals"
    else:
        vlabel, flabel, nlabel = None, None, None
    dataPoints = Gnuplot.Data(px, py,
                              with_ = "points pt 1 ps 2",
                              title = vlabel,
                              inline = True)
    dataFacets = Gnuplot.Data(fx, fy, fdx, fdy,
                              with_ = "vectors",
                              title = flabel,
                              inline = True)
    dataNormals = Gnuplot.Data(nx, ny, ndx, ndy,
                               with_ = "vectors",
                               title = nlabel,
                               inline = True)
    if plotVertices:
        plot.replot(dataPoints)

    if plotFacets:
        plot.replot(dataFacets)

    if plotNormals:
        plot.replot(dataNormals)

    if plotCentroid:
        c = polygon.centroid
        dataCentroid = Gnuplot.Data([c.x], [c.y],
                                    with_ = "points pt 2 ps 2",
                                    title = "Centroid",
                                    inline = True)
        plot.replot(dataCentroid)

    SpheralGnuPlotCache.extend([dataPoints, dataFacets, dataNormals, plot])

    return plot

#-------------------------------------------------------------------------------
# Plot a PolygonalMesh
#-------------------------------------------------------------------------------
def plotPolygonalMesh(mesh,
                      persist = False):
    polylocal = []
    for izone in range(mesh.numZones):
        zone = mesh.zone(izone)
        polylocal.append([mesh.node(i).position() for i in zone.nodeIDs])
        polylocal[-1].append(polylocal[-1][0])
    assert len(polylocal) == mesh.numZones

    p = generateNewGnuPlot(persist)
    for sendProc in range(mpi.procs):
        polys = mpi.bcast(polylocal, root=sendProc)
        for poly in polys:
            p.replot(Gnuplot.Data([x.x for x in poly], [x.y for x in poly],
                                  with_ = "lines lt %i lw 2" % 1,
                                  title = None,
                                  inline = True))
    return p



##     edges0 = [(mesh.node(mesh.edge(i).node1ID).position(), mesh.node(mesh.edge(i).node2ID).position())
##               for i in xrange(mesh.numEdges)]
##     p = generateNewGnuPlot()
##     datas = []
##     for sendProc in xrange(mpi.procs):
##         edges = mpi.bcast(edges0, root=sendProc)
##         for edge in edges:
##             datas.append(Gnuplot.Data([edge[0].x, edge[1].x], [edge[0].y, edge[1].y],
##                                       with_ = "lines %s" % linetype,
##                                       title = None,
##                                       inline = True))
##             p.replot(datas[-1])
##     p.datas = datas
##     return p
