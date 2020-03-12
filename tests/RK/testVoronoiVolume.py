#ATS:test(SELF, "--graphics False --nx1 10 --testDim 1d", label="computeVoronoiVolume test -- 1D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --testDim 2d", label="computeVoronoiVolume test -- 2D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --testDim 3d", label="computeVoronoiVolume test -- 3D (serial)")
#-------------------------------------------------------------------------------
# Unit test of the CRKSPH sum density algorithm.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
from SpheralVoronoiSiloDump import SpheralVoronoiSiloDump
import numpy as np

title("Voronoi volume tests")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Parameters for seeding nodes.
    nx1 = 10,
    rho1 = 10.0,
    x0 = 0.0,
    x1 = 1.0,
    nPerh = 1.25,
    hmin = 0.0001, 
    hmax = 10.0,

    # Should we randomly perturb the positions?
    ranfrac = 0.0,
    seed = 14892042,

    # What test problem are we doing?
    testDim = "3d",

    gamma = 5.0/3.0,
    mu = 1.0,

    # Parameters for iterating H.
    iterateH = False,
    maxHIterations = 200,
    Htolerance = 1.0e-4,

    # Parameters for passing the test
    tolerance = 1.0e-8,
    voidFlagTolerance = 0.2, # Relative difference in number of void points expected/calculated
    relaxVoidCheck = False, # Do not check the second layer of points from the boundary

    graphics = True,
    vizFile = "",
    splitCells = False,
)

assert testDim in ("1d", "2d", "3d")

FacetedVolume = {"1d" : Box1d,
                 "2d" : Polygon,
                 "3d" : Polyhedron}[testDim]

#-------------------------------------------------------------------------------
# Appropriately set generic object names based on the test dimensionality.
#-------------------------------------------------------------------------------
exec("from Spheral%s import *" % testDim)

#-------------------------------------------------------------------------------
# Create a random number generator.
#-------------------------------------------------------------------------------
import random
rangen = random.Random()
rangen.seed(seed)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if testDim == "1d":
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes1, [(nx1, rho1, (x0, x1))])], nPerh = nPerh)
elif testDim == "2d":
    from DistributeNodes import distributeNodes2d
    from GenerateNodeDistribution2d import GenerateNodeDistribution2d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen = GenerateNodeDistribution2d(nx1, nx1, rho1,
                                     distributionType = "lattice",
                                     xmin = (x0, x0),
                                     xmax = (x1, x1),
                                     nNodePerh = nPerh,
                                     SPH = True)
    distributeNodes2d((nodes1, gen))

elif testDim == "3d":
    from DistributeNodes import distributeNodes3d
    from GenerateNodeDistribution3d import GenerateNodeDistribution3d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen = GenerateNodeDistribution3d(nx1, nx1, nx1, rho1,
                                     distributionType = "lattice",
                                     xmin = (x0, x0, x0),
                                     xmax = (x1, x1, x1),
                                     nNodePerh = nPerh,
                                     SPH = True)
    distributeNodes3d((nodes1, gen))

else:
    raise ValueError, "Only tests cases for 1d, 2d and 3d." 

output("nodes1.numNodes")

#-------------------------------------------------------------------------------
# Optionally randomly jitter the node positions.
#-------------------------------------------------------------------------------
dx = (x1 - x0)/nx1
for i in xrange(nodes1.numInternalNodes):
    nodes1.positions()[i].x += ranfrac * dx * rangen.uniform(-1.0, 1.0)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
faceted_bounds = vector_of_FacetedVolume()
if testDim == "1d":
    point_list = [Vector(x0), Vector(x1)]
elif testDim == "2d":
    point_list = [Vector(x0, x0), Vector(x1, x0), 
                  Vector(x0, x1), Vector(x1, x1)]
else:
    point_list = [Vector(x0, x0, x0), Vector(x0, x1, x0), Vector(x1, x0, x0), Vector(x1, x1, x0), 
                  Vector(x0, x0, x1), Vector(x0, x1, x1), Vector(x1, x0, x1), Vector(x1, x1, x1)]
points = vector_of_Vector()
for p in point_list:
    points.append(p)
bbox = FacetedVolume(points)
faceted_bounds.append(bbox)

bounds = vector_of_Boundary()
# xbc0 = ReflectingBoundary(Plane(Vector(0.0, 0.0, 0.0), Vector( 1.0,  0.0,  0.0)))
# xbc1 = ReflectingBoundary(Plane(Vector(1.0, 0.0, 0.0), Vector(-1.0,  0.0,  0.0)))
# ybc0 = ReflectingBoundary(Plane(Vector(0.0, 0.0, 0.0), Vector( 0.0,  1.0,  0.0)))
# ybc1 = ReflectingBoundary(Plane(Vector(0.0, 1.0, 0.0), Vector( 0.0, -1.0,  0.0)))
# zbc0 = ReflectingBoundary(Plane(Vector(0.0, 0.0, 0.0), Vector( 0.0,  0.0,  1.0)))
# zbc1 = ReflectingBoundary(Plane(Vector(0.0, 0.0, 1.0), Vector( 0.0,  0.0, -1.0)))
# bounds.append(xbc0)
# bounds.append(xbc1)
# bounds.append(ybc0)
# bounds.append(ybc1)
# bounds.append(zbc0)
# bounds.append(zbc1)

#-------------------------------------------------------------------------------
# Iterate the h to convergence if requested.
#-------------------------------------------------------------------------------
if iterateH:
    method = SPHSmoothingScale()
    iterateIdealH(db,
                  bounds,
                  WT,
                  method,
                  maxHIterations,
                  Htolerance)

#-------------------------------------------------------------------------------
# Compute the volumes.
#-------------------------------------------------------------------------------
weight = ScalarFieldList()                         # No weights
damage = SymTensorFieldList()                      # No damage
holes = vector_of_vector_of_FacetedVolume([vector_of_FacetedVolume()])
surfacePoint = db.newFluidIntFieldList(0, HydroFieldNames.surfacePoint)
vol = db.newFluidScalarFieldList(0.0, HydroFieldNames.volume)
deltaMedian = db.newFluidVectorFieldList(Vector.zero, "centroidal delta")
etaVoidPoints = db.newFluidvector_of_VectorFieldList(vector_of_Vector(), "eta void points")
cells = db.newFluidFacetedVolumeFieldList(FacetedVolume(), "cells")
cellFaceFlags = db.newFluidvector_of_CellFaceFlagFieldList(vector_of_int(), "face flags")
db.updateConnectivityMap(True)
cm = db.connectivityMap()
computeVoronoiVolume(db.fluidPosition, 
                     db.fluidHfield,
                     cm,
                     damage,
                     faceted_bounds,
                     holes,
                     bounds,
                     weight,
                     surfacePoint,
                     vol,
                     deltaMedian,
                     etaVoidPoints,
                     cells,
                     cellFaceFlags)

#-------------------------------------------------------------------------------
# Optionally drop a viz file.
#-------------------------------------------------------------------------------
if vizFile:
    dumper = SpheralVoronoiSiloDump(baseFileName = vizFile,
                                    listOfFieldLists = [vol,
                                                        surfacePoint,
                                                        deltaMedian],
                                    cells = cells,
                                    splitCells = splitCells)
    dumper.dump(0.0, 0)

#-------------------------------------------------------------------------------
# Make sure some face flags exist and the right number of cells have flags
#-------------------------------------------------------------------------------
position = db.fluidPosition[0]
xbmin = x0 + (1.49 - ranfrac) * dx
xbmax = x1 - (1.49 - ranfrac) * dx
xbminrelaxed = x0 + (1.51 + ranfrac) * dx if relaxVoidCheck else xbmin
xbmaxrelaxed = x1 - (1.51 + ranfrac) * dx if relaxVoidCheck else xbmax
output("xbmin")
output("xbmax")
output("xbminrelaxed")
output("xbmaxrelaxed")

# Return 0 if not on boundary, 1 if definitely on boundary, 2 if one layer in from boundary
def shouldHaveVoid(i): 
    x = position[i]
    questionable = 0
    for xv in x:
        if xv < xbmin or xv > xbmax:
            return 1
        if xv < xbminrelaxed or xv > xbmaxrelaxed:
            questionable = 2
    return questionable

numCellFaceFlags = 0
numVoidFaceFlags = 0
numCellsWithVoidFlagsExpected = 0
numCellsWithVoidFlags = 0
hasVoidButShouldNot = 0
doesNotHaveVoidButShould = 0
position = db.fluidPosition[0]
print "Void errors, if any, are indented below"
for i in range(nodes1.numInternalNodes):
    flags = cellFaceFlags[0][i]
    cellCounted = False
    hasVoid = False
    for flag in flags:
        numCellFaceFlags += 1
        if flag.nodeListj == -1:
            hasVoid = True
            numVoidFaceFlags += 1
            if not cellCounted:
                numCellsWithVoidFlags += 1
                cellCounted = True
    expectedVoid = shouldHaveVoid(i)
    if expectedVoid == 1: # should have void
        numCellsWithVoidFlagsExpected += 1
        if not hasVoid:
            doesNotHaveVoidButShould += 1
            print "\tPosition: {}  Does not have void but should".format(position[i])
    elif expectedVoid == 0: # should not have void
        if hasVoid:
            hasVoidButShouldNot += 1
            print "\tPosition: {}  Has void but should not".format(position[i])

voidFlagRelErr = 1.0 * abs(numCellsWithVoidFlags - numCellsWithVoidFlagsExpected) / (1.0 * numCellsWithVoidFlagsExpected)
unexpectedFlagRelErr = 1.0 * (hasVoidButShouldNot + doesNotHaveVoidButShould) / (1.0 * nodes1.numInternalNodes)
output("numCellFaceFlags")
output("numVoidFaceFlags")
output("numCellsWithVoidFlags")
output("numCellsWithVoidFlagsExpected")
output("hasVoidButShouldNot")
output("doesNotHaveVoidButShould")
output("unexpectedFlagRelErr")
output("voidFlagRelErr")
assert numVoidFaceFlags > 0
assert numCellFaceFlags > 0
assert voidFlagRelErr < voidFlagTolerance
assert unexpectedFlagRelErr < voidFlagTolerance

#-------------------------------------------------------------------------------
# Check the answer.
#-------------------------------------------------------------------------------
if ranfrac == 0.0:
    for Vi, sp in zip(vol[0].internalValues(), surfacePoint[0].internalValues()):
        if sp == 0:
            assert abs(Vi - (1.0/nx1)**db.nDim) < 1.0e-12
        else:
            assert Vi == 0.0
else:
    print "random displacments -- no error check."
