import mpi
from Spheral2d import *
from MedialGenerator import *
from SpheralTestUtilities import *
from VoronoiDistributeNodes import distributeNodes2d as distributeNodes
from siloPointmeshDump import *

commandLine(n1       = 5000,
            rho0     = 1.0,
            R0       = 10.0,
            Rc       = 1.0,
            ncirc    = 360,
            hmin     = 1e-5,
            hmax     = 1e6,
            rhoscale = 0.5,

            nPerh   = 2.01,
            maxIterations = 200,
            fracTol = 1e-3,
            noholes = False)     # Optionally avoid using holes in the boundaries.

#-------------------------------------------------------------------------------
# The density profile we're going to fit.
#-------------------------------------------------------------------------------
def rhoprofile(posi):
    r2 = posi.magnitude2()
    return rho0/(r2 + Rc*Rc)

def gradrho(posi):
    r = posi.magnitude()
    rhat = posi.unitVector()
    return -2.0*rho0*r/(r*r + Rc*Rc)**2 * rhat

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
nodeSet = [nodes1]
for nodes in nodeSet:
    output("nodes.name")
    output("  nodes.hmin")
    output("  nodes.hmax")
    output("  nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Make a circular boundary.
#-------------------------------------------------------------------------------
bcpoints = vector_of_Vector()
bcfacets = vector_of_vector_of_unsigned()
for i in xrange(ncirc):
    theta = 2.0*pi/ncirc * i
    bcpoints.append(Vector(R0*cos(theta), R0*sin(theta)))
for i in xrange(len(bcpoints)):
    bcfacets.append(vector_of_unsigned(2))
    bcfacets[-1][0] = i
    bcfacets[-1][1] = (i + 1) % len(bcpoints)
boundary = Polygon(bcpoints, bcfacets)

#-------------------------------------------------------------------------------
# Generate them nodes.
#-------------------------------------------------------------------------------
generator1 = MedialGenerator2d(n = n1,
                               rho = rhoprofile,
                               gradrho = gradrho,   # This is not necessary, but we'll use it if provided
                               boundary = boundary,
                               maxIterations = maxIterations,
                               fracTol = fracTol,
                               tessellationFileName = "test_medial2d_sphere_density_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                               nNodePerh = nPerh)

distributeNodes((nodes1, generator1))

#-------------------------------------------------------------------------------
# Drop a viz file for inspection.
#-------------------------------------------------------------------------------
Hfield = nodes.Hfield()
HfieldInv = SymTensorField("H inverse", nodes)
domainField = IntField("Domain", nodes)
for i in xrange(nodes.numNodes):
    HfieldInv[i] = Hfield[i].Inverse()
    domainField[i] = mpi.rank
vizfile = siloPointmeshDump(baseName = "test_medial_maxiter=%i_tol=%g" % (maxIterations, fracTol),
                            baseDirectory = "test_medial2d_sphere_density",
                            fields = ([x.massDensity() for x in nodeSet] +
                                      [x.mass() for x in nodeSet] +
                                      [x.velocity() for x in nodeSet] +
                                      [x.specificThermalEnergy() for x in nodeSet] +
                                      [x.Hfield() for x in nodeSet] +
                                      [HfieldInv, domainField])
                            )
#-------------------------------------------------------------------------------
# Plot a few profiles of interest.
#-------------------------------------------------------------------------------
from SpheralGnuPlotUtilities import *
db = DataBase()
for nodes in nodeSet:
    db.appendNodeList(nodes)
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
