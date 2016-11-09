from Spheral2d import *
from CentroidalGenerator2d import *
from SpheralTestUtilities import *
from VoronoiDistributeNodes import distributeNodes2d as distributeNodes
from siloPointmeshDump import *

commandLine(hmin     = 1e-5,
            hmax     = 1e6,
            rhoscale = 0.5,
            n1       = 200,
            n2       = 100,
            n3       = 50,
            n4       = 400,
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
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           topGridCellSize = 100,
                           xmin = Vector.one * -100.0,
                           xmax = Vector.one *  100.0)
nodes2 = makeFluidNodeList("nodes2", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           topGridCellSize = 100,
                           xmin = Vector.one * -100.0,
                           xmax = Vector.one *  100.0)
nodes3 = makeFluidNodeList("nodes3", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           topGridCellSize = 100,
                           xmin = Vector.one * -100.0,
                           xmax = Vector.one *  100.0)
nodes4 = makeFluidNodeList("nodes4", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           topGridCellSize = 100,
                           xmin = Vector.one * -100.0,
                           xmax = Vector.one *  100.0)
nodeSet = [nodes1, nodes2, nodes3, nodes4]
for nodes in nodeSet:
    output("nodes.name")
    output("  nodes.hmin")
    output("  nodes.hmax")
    output("  nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Make some interesting boundaries for each of our NodeLists and generators.
#-------------------------------------------------------------------------------
# An H of height & width 3
bcpoints = vector_of_Vector()
bcfacets = vector_of_vector_of_unsigned()
for p in [(0,0), (1,0), (1,1), (2,1), (2,0), (3,0),
          (3,3), (2,3), (2,2), (1,2), (1,3), (0,3)]:
    bcpoints.append(Vector(*p))
for i in xrange(len(bcpoints)):
    bcfacets.append(vector_of_unsigned(2))
    bcfacets[-1][0] = i
    bcfacets[-1][1] = (i + 1) % len(bcpoints)
Hboundary = Polygon(bcpoints, bcfacets)

# A couple of concentric circles.
n = 120
bcpoints1 = vector_of_Vector(n)
bcpoints2 = vector_of_Vector(n)
bcfacets = vector_of_vector_of_unsigned(n, vector_of_unsigned(2))
for i in xrange(n):
    j = n - i - 1
    bcpoints1[i] = Vector(5.0 + 1.5*cos(i*2.0*pi/n), 1.5 + 1.5*sin(i*2.0*pi/n))
    bcpoints2[i] = Vector(5.0 + 0.5*cos(i*2.0*pi/n), 1.5 + 0.5*sin(i*2.0*pi/n))
    bcfacets[i][0] = i
    bcfacets[i][1] = (i + 1) % n
outerCircle = Polygon(bcpoints1, bcfacets)
innerCircle = Polygon(bcpoints2, bcfacets)

# A box surrounding the whole thing.
bcpoints = vector_of_Vector()
bcfacets = vector_of_vector_of_unsigned()
for p in [(-1,-1), (8,-1), (8,4), (-1,4)]:
    bcpoints.append(Vector(*p))
for i in xrange(len(bcpoints)):
    bcfacets.append(vector_of_unsigned(2))
    bcfacets[-1][0] = i
    bcfacets[-1][1] = (i + 1) % len(bcpoints)
outerBox = Polygon(bcpoints, bcfacets)

#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------
def rhoprofile1(posi):
    r = (posi - Vector(1.5,1.5)).magnitude()
    return exp(-r*r/(rhoscale*rhoscale))

print "Generator 1"
generator1 = CentroidalGenerator2d(n = n1,
                                   rho = rhoprofile1,
                                   boundary = Hboundary,
                                   maxIterations = maxIterations,
                                   fracTol = fracTol,
                                   tessellationFileName = "test_centroidal_nodes1_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                                   nNodePerh = nPerh)

print "Generator 2"
generator2 = CentroidalGenerator2d(n = n2,
                                   rho = 1.0,
                                   boundary = outerCircle,
                                   holes = [innerCircle],
                                   maxIterations = maxIterations,
                                   fracTol = fracTol,
                                   tessellationFileName = "test_centroidal_nodes2_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                                   nNodePerh = nPerh)

print "Generator 3"
generator3 = CentroidalGenerator2d(n = n3,
                                   rho = 0.1,
                                   boundary = innerCircle,
                                   maxIterations = maxIterations,
                                   fracTol = fracTol,
                                   tessellationFileName = "test_centroidal_nodes3_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                                   nNodePerh = nPerh)

print "Generator 4"
generator4 = CentroidalGenerator2d(n = n4,
                                   rho = 0.1,
                                   boundary = outerBox,
                                   holes = [Hboundary, outerCircle],
                                   maxIterations = maxIterations,
                                   fracTol = fracTol,
                                   tessellationFileName = "test_centroidal_nodes4_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                                   nNodePerh = nPerh)

distributeNodes((nodes1, generator1),
                (nodes2, generator2),
                (nodes3, generator3),
                (nodes4, generator4))

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
Hfield = nodes.Hfield()
HfieldInv = SymTensorField("H inverse", nodes)
for i in xrange(nodes.numNodes):
    HfieldInv[i] = Hfield[i].Inverse()
vizfile = siloPointmeshDump(baseName = "test_centroidal_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                            baseDirectory = "test_centroidal",
                            fields = ([x.massDensity() for x in nodeSet] +
                                      [x.mass() for x in nodeSet] +
                                      [x.velocity() for x in nodeSet] +
                                      [x.specificThermalEnergy() for x in nodeSet] +
                                      [x.Hfield() for x in nodeSet])
                            )
