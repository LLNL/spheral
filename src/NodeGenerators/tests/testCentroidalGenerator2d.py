from Spheral2d import *
from CentroidalGenerator2d import *
from SpheralTestUtilities import *
from VoronoiDistributeNodes import distributeNodes2d as distributeNodes
from siloPointmeshDump import *

commandLine(hmin     = 1e-5,
            hmax     = 1e6,
            rhoscale = 0.5,
            n        = 200,
            xcenter = 0.5,
            ycenter = 0.5,
            nPerh   = 2.01,
            maxIterations = 100,
            fracTol = 1e-2)

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
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Make an interesting boundary, in this case kind of an H.
#-------------------------------------------------------------------------------
bcpoints = vector_of_Vector()
bcfacets = vector_of_vector_of_unsigned()
for p in [(0,0), (1,0), (1,1), (2,1), (2,0), (3,0),
          (3,3), (2,3), (2,2), (1,2), (1,3), (0,3)]:
    bcpoints.append(Vector(*p))
for i in xrange(len(bcpoints)):
    bcfacets.append(vector_of_unsigned(2))
    bcfacets[-1][0] = i
    bcfacets[-1][1] = (i + 1) % len(bcpoints)
boundary = Polygon(bcpoints, bcfacets)

#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------
def rhoprofile(posi):
    r = posi.magnitude()
    return exp(-r*r/(rhoscale*rhoscale))

generator = CentroidalGenerator2d(n = n,
                                  rho = rhoprofile,
                                  boundary = boundary,
                                  maxIterations = maxIterations,
                                  fracTol = fracTol,
                                  tessellationFileName = "test_centroidal_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                                  nNodePerh = nPerh)
distributeNodes((nodes, generator))

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
Hfield = nodes.Hfield()
HfieldInv = SymTensorField("H inverse", nodes)
for i in xrange(nodes.numNodes):
    HfieldInv[i] = Hfield[i].Inverse()
vizfile = siloPointmeshDump(baseName = "test_centroidal_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                            baseDirectory = "test_centroidal",
                            fields = [nodes.massDensity(),
                                      nodes.mass(),
                                      nodes.velocity(),
                                      nodes.specificThermalEnergy(),
                                      Hfield,
                                      HfieldInv],
                            )
