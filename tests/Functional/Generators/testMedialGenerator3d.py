from Spheral3d import *
from MedialGenerator import *
from CompositeNodeDistribution import *
from SpheralTestUtilities import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes
from siloPointmeshDump import *

commandLine(hmin     = 1e-5,
            hmax     = 1e6,
            rhoscale = 0.5,
            n1       = 1000,
            n2       = 1000,
            nPerh   = 2.01,
            maxIterations = 200,
            fracTol = 1e-3)

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
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = 1.0,
                           nPerh = nPerh,
                           topGridCellSize = 100,
                           xmin = Vector.one * -100.0,
                           xmax = Vector.one *  100.0)
nodes2 = makeFluidNodeList("nodes2", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = 1.0,
                           nPerh = nPerh,
                           topGridCellSize = 100,
                           xmin = Vector.one * -100.0,
                           xmax = Vector.one *  100.0)
nodeSet = [nodes1, nodes2]
for nodes in nodeSet:
    output("nodes.name")
    output("  nodes.hmin")
    output("  nodes.hmax")
    output("  nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Make some interesting boundaries for each of our NodeLists and generators.
#-------------------------------------------------------------------------------
# The inner cube.
bcpoints = vector_of_Vector()
for p in [(1,1,1), (1,2,1), (2,1,1), (2,2,1),
          (1,1,2), (1,2,2), (2,1,2), (2,2,2)]:
    bcpoints.append(Vector(*p))
innerBoundary = Polyhedron(bcpoints)     # Builds the convex hull

# The outer cube.
bcpoints = vector_of_Vector()
for p in [(0,0,0), (0,3,0), (3,0,0), (3,3,0),
          (0,0,3), (0,3,3), (3,0,3), (3,3,3)]:
    bcpoints.append(Vector(*p))
outerBoundary = Polyhedron(bcpoints)     # Builds the convex hull

#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------
def rhoprofile1(posi):
    r = (posi - Vector(1.5,1.5)).magnitude()
    return exp(-r*r/(rhoscale*rhoscale))

print "Generator 1"
generator1 = MedialGenerator3d(n = n1,
                               rho = 1.0,
                               boundary = innerBoundary,
                               maxIterations = maxIterations,
                               fracTol = fracTol,
                               #tessellationFileName = "test_medial_nodes1_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                               nNodePerh = nPerh)

print "Generator 2"
generator2 = MedialGenerator3d(n = n2,
                               rho = 1.0,
                               boundary = outerBoundary,
                               holes = [innerBoundary],
                               maxIterations = maxIterations,
                               fracTol = fracTol,
                               #tessellationFileName = "test_medial_nodes2_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                               nNodePerh = nPerh)

distributeNodes((nodes1, generator1),
                (nodes2, generator2))

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
Hfield = nodes.Hfield()
HfieldInv = SymTensorField("H inverse", nodes)
for i in xrange(nodes.numNodes):
    HfieldInv[i] = Hfield[i].Inverse()
vizfile = siloPointmeshDump(baseName = "test_medial3d_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                            baseDirectory = "test_medial3d",
                            fields = ([x.massDensity() for x in nodeSet] +
                                      [x.mass() for x in nodeSet] +
                                      [x.velocity() for x in nodeSet] +
                                      [x.specificThermalEnergy() for x in nodeSet] +
                                      [x.Hfield() for x in nodeSet])
                            )
