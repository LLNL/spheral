import mpi
from Spheral3d import *
from MedialGenerator import *
from SpheralTestUtilities import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from siloPointmeshDump import *
from PlyFileIO import *

commandLine(npart      = 100,
            filename   = "untitled.ply",
            hmin       = 1e-5,
            hmax       = 1e6,
            rho0       = 1.0,

            nPerh      = 2.01,
            centroidFrac = 1.0,
            maxIterations = 1000,
            fracTol    = 5e-4)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
gamma = 1.4
mu = 2.0
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes", eos,
                              hmin = hmin,
                              hmax = hmax,
                              nPerh = nPerh,
                              topGridCellSize = 100,
                              xmin = Vector.one * -100.0,
                              xmax = Vector.one *  100.0)
nodeSet = [nodes]
for nodes in nodeSet:
    output("nodes.name")
    output("  nodes.hmin")
    output("  nodes.hmax")
    output("  nodes.nodesPerSmoothingScale")


boundaryShape = PolyFromPly(filename)

#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------

generator = MedialGenerator3d(n = npart,
                                        rho = rho0,
                                        boundary = boundaryShape,
                                        centroidFrac = centroidFrac,
                                        maxIterations = maxIterations,
                                        fracTol = fracTol,
                                        nNodePerh = nPerh)

distributeNodes((nodes, generator))

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
db = DataBase()
for nodes in nodeSet:
    db.appendNodeList(nodes)
vizfile = siloPointmeshDump(baseName = "test_medial_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                            baseDirectory = "test_medial3d_sphere_density",
                            fieldLists = [db.fluidMassDensity,
                                          db.fluidMass,
                                          db.fluidVelocity,
                                          db.fluidSpecificThermalEnergy,
                                          db.fluidHfield]
                            )

'''
#-------------------------------------------------------------------------------
# Plot a few profiles of interest.
#-------------------------------------------------------------------------------
from SpheralGnuPlotUtilities import *
massPlot = plotFieldList(db.fluidMass,
                         xFunction = "%s.magnitude()",
                         plotStyle = "points",
                         winTitle = "mass",
                         colorNodeLists = False, plotGhosts = False)
rhoPlot = plotFieldList(db.fluidMassDensity,
                        xFunction = "%s.magnitude()",
                        plotStyle = "points",
                        winTitle = "mass density",
                        colorNodeLists = False, plotGhosts = False)
rhoPlot("set yrange [1e-2:200]; set logscale y"); rhoPlot.refresh()
massPlot.hardcopy("test_medial3d_mass.png", terminal="png")
rhoPlot.hardcopy("test_medial3d_rho.png", terminal="png")

from fieldStatistics import fieldStatistics
for nodes in nodeSet:
    print "Mass statistics for ", nodes.name, " (min, max, avg, std dev) : ", fieldStatistics(nodes.mass())'''
