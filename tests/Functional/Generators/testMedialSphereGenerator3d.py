import mpi
from Spheral3d import *
from MultiScaleMedialGenerator import *
from SpheralTestUtilities import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from siloPointmeshDump import *

commandLine(ncore      = 10000,
            rhocore0   = 10.0,
            rhomantle0 = 1.0,
            Rcore      = 1.0,
            Rmantle    = 10.0,
            Rc         = 0.25,
            nshell     = 300,
            hmin       = 1e-5,
            hmax       = 1e6,

            nPerh      = 2.01,
            centroidFrac = 1.0,
            maxIterations = 1000,
            fracTol    = 5e-4)

#-------------------------------------------------------------------------------
# The density profiles we're going to fit.
# Note we don't have to provide the rho gradient methods, but providing them is
# probably more accurate and they're trivial to compute for these profiles.
#-------------------------------------------------------------------------------
def rhocore(posi):
    if type(posi) is type(0.0):
        r2 = posi**2
    else:
        r2 = posi.magnitude2()
    return rhocore0/(r2 + Rc*Rc)

def gradrhocore(posi):
    r = posi
    rhat = posi.unitVector()
    return -2.0*rhocore0*r/(r*r + Rc*Rc)**2 * rhat

def rhomantle(posi):
    if type(posi) is type(0.0):
        r2 = posi**2
    else:
        r2 = posi.magnitude2()
    return rhomantle0/r2

def gradrhomantle(posi):
    r = posi.magnitude()
    rhat = posi.unitVector()
    return -2.0*rhomantle0/(r*r*r) * rhat

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
nodesCore = makeFluidNodeList("core", eos,
                              hmin = hmin,
                              hmax = hmax,
                              nPerh = nPerh,
                              topGridCellSize = 100,
                              xmin = Vector.one * -100.0,
                              xmax = Vector.one *  100.0)
nodesMantle = makeFluidNodeList("mantle", eos,
                                hmin = hmin,
                                hmax = hmax,
                                nPerh = nPerh,
                                topGridCellSize = 100,
                                xmin = Vector.one * -100.0,
                                xmax = Vector.one *  100.0)
nodeSet = [nodesCore, nodesMantle]
for nodes in nodeSet:
    output("nodes.name")
    output("  nodes.hmin")
    output("  nodes.hmax")
    output("  nodes.nodesPerSmoothingScale")


#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------
# First, figure out the appropriate number of nodes we should have in the mantle
# to mass match those in the core.
Mcore = 4.0*pi*rhocore0*(Rcore - Rc*atan2(Rcore, Rc))
Mmantle = 4.0*pi*rhomantle0*(Rmantle - Rcore)
nmantle = int(Mmantle/Mcore*ncore + 0.5)
print "  Core mass: ", Mcore
print "Mantle mass: ", Mmantle
print "Resulting target point mass and number of points in mantle: ", Mcore/ncore, nmantle

generatorCore = MedialSphereGenerator3d(n = ncore,
                                        rho = rhocore,
                                        gradrho = None,   # This is not necessary, but we'll use it if provided
                                        rmin = 0.0,
                                        rmax = Rcore,
                                        centroidFrac = centroidFrac,
                                        maxIterations = maxIterations,
                                        fracTol = fracTol,
                                        nNodePerh = nPerh)

generatorMantle = MedialSphereGenerator3d(n = nmantle,
                                          rho = rhomantle,
                                          gradrho = None,   # This is not necessary, but we'll use it if provided
                                          rmin = Rcore,
                                          rmax = Rmantle,
                                          centroidFrac = centroidFrac,
                                          maxIterations = maxIterations,
                                          fracTol = fracTol,
                                          nNodePerh = nPerh)

distributeNodes((nodesCore, generatorCore),
                (nodesMantle, generatorMantle))

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
    print "Mass statistics for ", nodes.name, " (min, max, avg, std dev) : ", fieldStatistics(nodes.mass())
