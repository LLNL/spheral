import pylab
import Numeric
import loadmpi
import Spheral
import SpheralTestUtilities
from math import *

from SpheralGnuPlotUtilities import multiSort, \
     radialVelocityFieldList, \
     azimuthalVelocityFieldList, \
     angularMomentum, \
     gridSample

#-------------------------------------------------------------------------------
# Module variables.
#-------------------------------------------------------------------------------
__SpheralMatplotlibUtilities_figure = 0
__matplotlib_availableColors = ["b", "g", "r", "c", "m", "y", "k", "w"]
mpi, rank, nprocs = loadmpi.loadmpi()

#-------------------------------------------------------------------------------
# Plot a FieldList
#-------------------------------------------------------------------------------
def plotFieldList(userFieldList,
                  xFunction = "%s.x",
                  yFunction = "%s",
                  plotGhosts = 0,
                  colorNodeLists = 1,
                  colorDomains = 0,
                  plot = None,
                  userXRange = [None, None],
                  userYRange = [None, None],
                  plotStyle = "-",
                  plotColor = "b",
                  winTitle = None,
                  lineTitle = "",
                  xlabel = None,
                  ylabel = None,
                  figureNum = None):

    global __SpheralMatplotlibUtilities_figure, __matplotlib_availableColors, mpi, rank, nprocs

    assert colorNodeLists + colorDomains <= 1
    assert plotColor in __matplotlib_availableColors

    if figureNum is None:
        figureNum = __SpheralMatplotlibUtilities_figure
        __SpheralMatplotlibUtilities_figure += 1

    # If we were given a Field, make it a FieldList.
    field2fieldList = {type(Spheral.ScalarField1d("thpt")): Spheral.ScalarFieldList1d,
                       type(Spheral.VectorField1d("thpt")): Spheral.VectorFieldList1d,
                       type(Spheral.TensorField1d("thpt")): Spheral.TensorFieldList1d,
                       type(Spheral.SymTensorField1d("thpt")): Spheral.SymTensorFieldList1d,
                       type(Spheral.ScalarField2d("thpt")): Spheral.ScalarFieldList2d,
                       type(Spheral.VectorField2d("thpt")): Spheral.VectorFieldList2d,
                       type(Spheral.TensorField2d("thpt")): Spheral.TensorFieldList2d,
                       type(Spheral.SymTensorField2d("thpt")): Spheral.SymTensorFieldList2d,
                       type(Spheral.ScalarField3d("thpt")): Spheral.ScalarFieldList3d,
                       type(Spheral.VectorField3d("thpt")): Spheral.VectorFieldList3d,
                       type(Spheral.TensorField3d("thpt")): Spheral.TensorFieldList3d,
                       type(Spheral.SymTensorField3d("thpt")): Spheral.SymTensorFieldList3d
                       }
    fieldList = None
    for x in field2fieldList:
        if type(userFieldList) == x:
            fieldList = field2fieldList[x]()
            fieldList.appendField(userFieldList)
    if fieldList is None:
        fieldList = userFieldList

    # Gather the fieldList info across all processors to process 0.
    localNumNodes = []
    localX = []
    localY = []
    for field in fieldList.fields():
        if plotGhosts:
            n = field.nodeList().numNodes
        else:
            n = field.nodeList().numInternalNodes
        localNumNodes.append(n)
        for x in field.nodeList().positions().allValues()[:n]:
            localX.append(eval(xFunction % "x"))
        for y in field.allValues()[:n]:
            localY.append(eval(yFunction % "y"))
    globalNumNodes = mpi.gather(localNumNodes, len(localNumNodes))
    globalX = mpi.gather(localX, len(localX))
    globalY = mpi.gather(localY, len(localY))

    if rank == 0:
        # Find the total number of nodes.
        totalNumNodes = 0
        for n in globalNumNodes:
            totalNumNodes += n
        assert(len(globalNumNodes) == fieldList.numFields()*nprocs)
        assert(len(globalX) == totalNumNodes)
        assert(len(globalY) == totalNumNodes)

        # Copy the input ranges, since for some reason these seem to have been
        # preserved between calls?
        xRange = userXRange[:]
        yRange = userYRange[:]

        # Generate a new matplotlib figure
        pylab.figure(figureNum)

        # Set the labels.
        if winTitle:
            pylab.title(winTitle)
        if xlabel:
            pylab.xlabel(xlabel)
        if ylabel:
            pylab.ylabel(ylabel)

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
        if yRange[0] == None: yRange[0] = ymin - 0.05*(ymax - ymin)
        if yRange[1] == None: yRange[1] = ymax + 0.05*(ymax - ymin)

        # Finally, loop over the fields and do the deed.
        assert len(globalX) == len(globalY)
        if colorDomains:
            numNodesPerDomain = [0]
            for domainBegin in xrange(0, len(globalNumNodes), fieldList.numFields()):
                numNodesPerDomain.append(numNodesPerDomain[-1])
                for n in globalNumNodes[domainBegin:domainBegin + fieldList.numFields()]:
                    numNodesPerDomain[-1] += n
            for domain in xrange(nprocs):
                x = Numeric.array(globalX[numNodesPerDomain[domain]:
                                          numNodesPerDomain[domain + 1]])
                y = Numeric.array(globalY[numNodesPerDomain[domain]:
                                          numNodesPerDomain[domain + 1]])
                if x:
                    pylab.plot(x, y,
                               "%s%s" % (__matplotlib_availableColors[domain % len(__matplotlib_availableColors)],
                                         plotStyle),
                               label = "Domain %i" % domain)

        elif colorNodeLists:
            legendNodeList = {}
            labeledNodeList = {}
            for i in xrange(fieldList.numFields()):
                legendNodeList[i] = lineTitle + ": " + fieldList[i].nodeList().name()
                labeledNodeList[i] = False

            cumulativeNumNodes = 0
            for fieldID in xrange(len(globalNumNodes)):
                n = globalNumNodes[fieldID]
                iNodeList = fieldID % fieldList.numFields()
                x = Numeric.array(globalX[cumulativeNumNodes:
                                  cumulativeNumNodes + n])
                y = Numeric.array(globalY[cumulativeNumNodes:
                                  cumulativeNumNodes + n])
                if n:
                    if labeledNodeList[iNodeList]:
                        pylab.plot(x, y,
                                   "%s%s" % (__matplotlib_availableColors[iNodeList % len(__matplotlib_availableColors)],
                                             plotStyle))
                    else:
                        pylab.plot(x, y,
                                   "%s%s" % (__matplotlib_availableColors[iNodeList % len(__matplotlib_availableColors)],
                                             plotStyle),
                                   label = legendNodeList[iNodeList])
                        labeledNodeList[iNodeList] = True
                    cumulativeNumNodes += n
                
        else:
            x = Numeric.array(globalX)
            y = Numeric.array(globalY)
            pylab.plot(x, y,
                       "%s%s" % (plotColor, plotStyle),
                       label = lineTitle)

        pylab.legend(loc = "best")

    # That's it.
    mpi.barrier()
    return figureNum

#-------------------------------------------------------------------------------
# Plot the mass density, velocity, pressure, and smoothing scale for the fluid
# node lists in the given data base.  Implicitly assuming 1-D.
#-------------------------------------------------------------------------------
def plotState(dataBase,
              plotGhosts = 0,
              colorNodeLists = 1,
              colorDomains = 0,
              plotStyle = "-o",
              xFunction = "%s.x",
              vecyFunction = "%s.x",
              tenyFunction = "%s.xx ** -1",
              lineTitle = "Simulation"):

    rhoPlot = plotFieldList(dataBase.fluidMassDensity,
                            xFunction = xFunction,
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            colorDomains = colorDomains,
                            plotStyle = plotStyle,
                            winTitle = "Mass Density",
                            lineTitle = lineTitle,
                            xlabel="x")

    velPlot = plotFieldList(dataBase.fluidVelocity,
                            xFunction = xFunction,
                            yFunction = vecyFunction,
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            colorDomains = colorDomains,
                            plotStyle = plotStyle,
                            winTitle = "Velocity",
                            lineTitle = lineTitle,
                            xlabel="x")

    epsPlot = plotFieldList(dataBase.fluidSpecificThermalEnergy,
                            xFunction = xFunction,
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            colorDomains = colorDomains,
                            plotStyle = plotStyle,
                            winTitle = "Specific Thermal Energy",
                            lineTitle = lineTitle,
                            xlabel="x")

    PPlot = plotFieldList(dataBase.fluidPressure,
                          xFunction = xFunction,
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          colorDomains = colorDomains,
                          plotStyle = plotStyle,
                          winTitle = "Pressure",
                          lineTitle = lineTitle,
                          xlabel="x")

    HPlot = plotFieldList(dataBase.fluidHfield,
                          xFunction = xFunction,
                          yFunction = "%s**-1" % tenyFunction,
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          colorDomains = colorDomains,
                          plotStyle = plotStyle,
                          winTitle = "Smoothing scale",
                          lineTitle = lineTitle,
                          xlabel="x")

    return rhoPlot, velPlot, epsPlot, PPlot, HPlot


#-------------------------------------------------------------------------------
# Plot the state vs. radius
#-------------------------------------------------------------------------------
def plotRadialState(dataBase,
                    plotGhosts = 0,
                    colorNodeLists = 1,
                    colorDomains = 0,
                    plotStyle = "o",
                    lineTitle = "Simulation"):

    rhoPlot = plotFieldList(dataBase.fluidMassDensity,
                            xFunction = "%s.magnitude()",
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            colorDomains = colorDomains,
                            plotStyle = plotStyle,
                            winTitle = "Mass density",
                            lineTitle = lineTitle,
                            xlabel = "r")

    radialVelocity = radialVelocityFieldList(dataBase.fluidPosition,
                                             dataBase.fluidVelocity)
    velPlot = plotFieldList(radialVelocity,
                            xFunction = "%s.magnitude()",
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            colorDomains = colorDomains,
                            plotStyle = plotStyle,
                            winTitle = " Radial Velocity",
                            lineTitle = lineTitle,
                            xlabel = "r")

    epsPlot = plotFieldList(dataBase.fluidSpecificThermalEnergy,
                            xFunction = "%s.magnitude()",
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            colorDomains = colorDomains,
                            plotStyle = plotStyle,
                            winTitle = "Specific Thermal Energy",
                            lineTitle = lineTitle,
                            xlabel = "r")

    PPlot = plotFieldList(dataBase.fluidPressure,
                          xFunction = "%s.magnitude()",
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          colorDomains = colorDomains,
                          plotStyle = plotStyle,
                          winTitle = "Pressure",
                          lineTitle = lineTitle,
                          xlabel = "r")

    HPlot = plotFieldList(dataBase.fluidHfield,
                          xFunction = "%s.magnitude()",
                          yFunction = "%s.xx**-1",
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          colorDomains = colorDomains,
                          plotStyle = plotStyle,
                          winTitle = "Smoothing scale",
                          lineTitle = lineTitle,
                          xlabel = "r")

    return rhoPlot, velPlot, epsPlot, PPlot, HPlot

#-------------------------------------------------------------------------------
# Overplot the answer on results from plotState.
#-------------------------------------------------------------------------------
def plotAnswer(answerObject, time,
               rhoPlot = None,
               velPlot = None,
               epsPlot = None,
               PPlot = None,
               HPlot = None):

    global __SpheralMatplotlibUtilities_figure, __matplotlib_availableColors, mpi, rank, nprocs

    try:
        x, v, u, rho, P, h = answerObject.solution(time)
    except:
        x, v, u, rho, P = answerObject.solution(time)

    if rank == 0:
        if rhoPlot is not None:
            pylab.figure(rhoPlot)
            pylab.plot(x, rho, "r-",
                       label = "Solution")
            pylab.legend(loc = "best")

        if velPlot is not None:
            pylab.figure(velPlot)
            pylab.plot(x, v, "r-",
                       label = "Solution")
            pylab.legend(loc = "best")

        if epsPlot is not None:
            pylab.figure(epsPlot)
            pylab.plot(x, u, "r-",
                       label = "Solution")
            pylab.legend(loc = "best")

        if PPlot is not None:
            pylab.figure(PPlot)
            pylab.plot(x, P, "r-",
                       label = "Solution")
            pylab.legend(loc = "best")

        if HPlot is not None:
            pylab.figure(HPlot)
            pylab.plot(x, h, "r-",
                       label = "Solution")
            pylab.legend(loc = "best")

    return

#-------------------------------------------------------------------------------
# Plot the node positions
#-------------------------------------------------------------------------------
def plotNodePositions2d(dataBase,
                        plotGhosts = 0,
                        colorNodeLists = 1,
                        colorDomains = 0,
                        title = "",
                        persist = None):

    global __SpheralMatplotlibUtilities_figure, __matplotlib_availableColors, mpi, rank, nprocs

    assert colorNodeLists + colorDomains <= 1

    # Gather the node positions across all domains.
    # Loop over all the NodeLists.
    localNumNodes = []
    xNodes = []
    yNodes = []
    for nodeList in dataBase.nodeLists():
        if plotGhosts:
            n = nodeList.numNodes
            xNodes.extend([x.x for x in nodeList.positions().allValues()])
            yNodes.extend([x.y for x in nodeList.positions().allValues()])
        else:
            n = nodeList.numInternalNodes
            xNodes.extend([x.x for x in nodeList.positions().internalValues()])
            yNodes.extend([x.y for x in nodeList.positions().internalValues()])
        localNumNodes.append(n)
    assert len(xNodes) == len(yNodes)
    
    numDomainNodes = [len(xNodes)]
    numNodesPerDomain = mpi.gather(numDomainNodes, 1)
    globalNumNodes = mpi.gather(localNumNodes, len(localNumNodes))
    globalXNodes = mpi.gather(xNodes, len(xNodes))
    globalYNodes = mpi.gather(yNodes, len(yNodes))

    if rank == 0:

        # Generate a new matplotlib figure
        pylab.figure(__SpheralMatplotlibUtilities_figure)
        __SpheralMatplotlibUtilities_figure += 1

        pylab.title(title)

        if colorDomains:
            cumulativeN = 0
            for domain in xrange(len(numNodesPerDomain)):
                n = numNodesPerDomain[domain]
                x = Numeric.array(globalXNodes[cumulativeN:cumulativeN + n])
                y = Numeric.array(globalYNodes[cumulativeN:cumulativeN + n])
                cumulativeN += n
                pylab.plot(x, y,
                           "%s." % __matplotlib_availableColors[domain],
                           label = "Domain %i" % domain)

        elif colorNodeLists:
            legendNodeList = {}
            for i in xrange(dataBase.numNodeLists):
                legendNodeList[i] = dataBase.nodeLists()[i].name()
            cumulativeN = 0
            for i in xrange(len(globalNumNodes)):
                n = globalNumNodes[i]
                if n > 0:
                    iNodeList = i % dataBase.numNodeLists
                    x = Numeric.array(globalXNodes[cumulativeN:cumulativeN + n])
                    y = Numeric.array(globalYNodes[cumulativeN:cumulativeN + n])
                    cumulativeN += n
                    pylab.plot(x, y,
                               "%s." % __matplotlib_availableColors[iNodeList],
                               label = legendNodeList[iNodeList])

        else:
            x = Numeric.array(globalXNodes)
            y = Numeric.array(globalYNodes)
            pylab.plot(x, y,
                       "b.")

    return __SpheralMatplotlibUtilities_figure - 1

#-------------------------------------------------------------------------------
# Plot the energy history of the given conservation object.
#-------------------------------------------------------------------------------
def plotEHistory(conserve):
    global __SpheralMatplotlibUtilities_figure, __matplotlib_availableColors, mpi, rank, nprocs
    if rank == 0:
        t = conserve.timeHistory
        E = conserve.EHistory
        KE = conserve.KEHistory
        TE = conserve.TEHistory
        UE = conserve.EEHistory

        pylab.figure(__SpheralMatplotlibUtilities_figure)
        __SpheralMatplotlibUtilities_figure += 1

        pylab.plot(t, E, "k-",
                   label = "Total Energy")
        pylab.plot(t, KE, "b-+",
                   label = "Kinetic Energy")
        pylab.plot(t, TE, "r-o",
                   label = "Thermal Energy")
        pylab.plot(t, UE, "m--",
                   label = "Potential Energy")
        pylab.title("Energy history")
        pylab.xlabel("time")
        pylab.ylabel("Energy")
        pylab.legend(loc = "best")

    return __SpheralMatplotlibUtilities_figure - 1
