#-------------------------------------------------------------------------------
# The evolution of a uniform, magnetized conducting fluid.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralVisitDump import dumpPhysicsState
from findLastRestart import *

# Load the mpi module if we"re parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

from GenerateNodeDistribution3d import *
from CubicNodeGenerator import GenerateCubicNodeDistribution

title("Dedner magnetic divergence test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            n = 20,
            rho0 = 1.0,
            V0 = Vector3d(1.0, 1.0, 0.0),
            Bz = 1.0/sqrt(4*pi),
            P0 = 6.0,
            nPerh = 1.3,
            mu0 = 1.0,
            gamma = 5.0/3.0,
            r0 = 1.0/sqrt(8),
            divBCleaner = 'none',
            mu = 1.0,
            Qlimiter = True,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 1.0,
            hminratio = 0.05,
            HsmoothFraction = 0.0,
            cfl = 0.25,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            HEvolution = Hydro3d.HEvolutionType.IdealH,
            compatibleEnergy = False,
            gradhCorrection = True,
            limitIdealH = False,

            neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector3d(0.0, 0.0, 0.0),

            goalTime = 1.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            sumForMassDensity = Hydro3d.MassDensityType.RigorousSumDensity,

            restoreCycle = None,
            graphics = False,
            )

def plotField(x, F, titleStr, filename):
   import pylab as p
   import griddata as g
   import numpy
   p.ion()
   p.clf()
   xhat = Vector3d(1, 0, 0)
   yhat = Vector3d(0, 1, 0)
   numInternalNodes = len(x.internalValues())
   indices = [i for i in xrange(numInternalNodes) if abs(x[i].z) < 1e-8]
   xs = numpy.array([x[i].dot(xhat) for i in indices])
   ys = numpy.array([x[i].dot(yhat) for i in indices])
   x1 = p.linspace(-0.5, 1.5, 50)
   y1 = p.linspace(-0.5, 1.5, 50)
   xg, yg = p.meshgrid(x1, y1)
   if isinstance(F, VectorField3d) or isinstance(F[0], Vector3d):
      Fxs = numpy.array([F[i].dot(xhat) for i in indices])
      Fys = numpy.array([F[i].dot(yhat) for i in indices])
      Fxg = g.griddata(xs, ys, Fxs, xg, yg)
      Fyg = g.griddata(xs, ys, Fys, xg, yg)
      p.quiver(xg, yg, Fxg, Fyg)
   else:
#      levels = [0.1*i for i in xrange(32)]
      Fs = numpy.array([F[i] for i in indices])
      Fg = g.griddata(xs, ys, Fs, xg, yg)
      p.contour(xg, yg, Fg, 30)
      p.colorbar()
   p.title(titleStr)
   p.savefig(filename)
#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 1000)
WTPi = TableKernel3d(BSplineKernel3d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# A few derived variables.
#-------------------------------------------------------------------------------
nx = ny = n
nz = int(2 * 2 * kernelExtent * nPerh)
nzx = 1.0*nz/nx
xmin = (-0.5, -0.5, -0.5*nzx)
xmax = (1.5, 1.5, 1.5*nzx)
u0 = P0 / ((gamma-1.0)*rho0)

dataDir = "Dedner-divB-%ix%ix%i" % (n, n, n)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS3d(gamma, mu)

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
# Perfectly conducting node list.
nodes = ConductingFluidNodeList("nodes", eos, WT, WTPi)
output("nodes")
nodes.HsmoothFraction = HsmoothFraction
nodes.XSPH = XSPH
nodes.nodesPerSmoothingScale = nPerh
nodes.epsilonTensile = epsilonTensile
nodes.nTensile = nTensile
nodes.hmin = hmin
nodes.hmax = hmax
nodes.hminratio = hminratio
output("nodes.HsmoothFraction")
output("nodes.nodesPerSmoothingScale")
output("nodes.epsilonTensile")
output("nodes.nTensile")
output("nodes.XSPH")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.hminratio")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor3d(nodes,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
x = nodes.positions()
v = nodes.velocity()
B = nodes.magneticInduction()
if restoreCycle is None:
    from ParMETISDistributeNodes import distributeNodes3d
    generator = GenerateNodeDistribution3d(nx, ny, nz, rho0, seed,
                                           xmin = xmin,
                                           xmax = xmax,
                                           nNodePerh = nPerh,
                                           SPH = True)
    distributeNodes3d((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes.specificThermalEnergy(ScalarField3d("tmp", nodes, u0))

    # Set nodal magnetic inductions.
    r = [sqrt(xi.x**2 + xi.y**2) for xi in x.internalValues()]
    for nodeID in xrange(nodes.numInternalNodes):
        ri = r[nodeID]/r0
        if ri < 1.0:
           Bx = (ri**8 - 2*ri**4 + 1)/sqrt(4*pi)
        else:
           Bx = 0.0
        B[nodeID] = Vector3d(Bx, 0, Bz)
        v[nodeID] = V0

    # Plot the B field configuration "before."
    #plotField(x, B, 'B before div cleaning', 'B-before.png')
#    plotField(x, [Bi.x for Bi in B.internalValues()], 'Bx before div cleaning', 'Bx-before.png')

    # Jot down the analytic maximum divergence of B.  The expression for 
    # div B = dBx/dx + dBy/dy + dBz/dz is (16*x*r**2/r0**4)*((r/r0)**4 - 1).
    #proj = Vector3d(1., 1., 0)
    #rs = [xi.dot(proj) for xi in x.internalValues()]
    #divBs = [(16*x[i].x*rs[i]**2/r0**4)*((rs[i]/r0)**4 - 1) for i in xrange(len(x.internalValues()))]
    #maxDivB0 = max(divBs)

    # Plot div B "before."
    #plotField(x, divBs, 'div B before div cleaning', 'divB-before.png')


#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
simName = 'Dedner-divB-%ix%ix%i'%(n, n, nz)
dataDir = '/p/lscratcha/jnjohnso/' + simName
visitDir = dataDir + "/visit"
restartDir = dataDir + "/restart"
import os, sys
if mpi.rank == 0:
    if restoreCycle is None:
        import shutil
        if os.path.exists(visitDir):
            shutil.rmtree(visitDir)
        if os.path.exists(restartDir):
            shutil.rmtree(restartDir)
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = MonaghanGingoldViscosity3d(1.0, 0.75)
#q = PriceMonaghanDissipation(1.0, 1.0, 1.0, 0.75, 1.0)

##-------------------------------------------------------------------------------
## Construct the hydro physics object.
##-------------------------------------------------------------------------------
hydro = Hydro3d(WT, WTPi, q, compatibleEnergy, gradhCorrection)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
#output("hydro")
#output("hydro.cfl")
#output("hydro.HEvolution")
#output("hydro.sumForMassDensity")
#output("hydro.HsmoothMin")
#output("hydro.HsmoothMax")
#output("hydro.kernel()")
#output("hydro.PiKernel()")
#output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct an MHD object.
#-------------------------------------------------------------------------------
mhd = MHD(WT, mu0)
if divBCleaner == 'none':
   mhd.divBCleaner = MHD.BDivergenceCleanerType.noCleaner
elif divBCleaner == 'hyperbolic':
   mhd.divBCleaner = MHD.BDivergenceCleanerType.hyperbolicCleaner
elif divBCleaner == 'GreensFn':
   mhd.divBCleaner = MHD.BDivergenceCleanerType.GreensFnProjCleaner
elif divBCleaner == 'BiotSavart':
   mhd.divBCleaner = MHD.BDivergenceCleanerType.BiotSavartProjCleaner
else:
   raise ValueError, "divBCleaner must be 'hyperBolic', 'GreensFn', 'BiotSavart', or 'none'."

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane1 = Plane3d(Vector3d(-0.5, 0.0, 0.0), Vector3d( 1.0, 0.0, 0.0))
xPlane2 = Plane3d(Vector3d( 1.5, 0.0, 0.0), Vector3d(-1.0, 0.0, 0.0))
yPlane1 = Plane3d(Vector3d( 0.0,-0.5, 0.0), Vector3d( 0.0, 1.0, 0.0))
yPlane2 = Plane3d(Vector3d( 0.0, 1.5, 0.0), Vector3d( 0.0,-1.0, 0.0))
zPlane1 = Plane3d(Vector3d( 0.0, 0.0,-0.5*nzx), Vector3d( 0.0, 0.0, 1.0))
zPlane2 = Plane3d(Vector3d( 0.0, 0.0, 1.5*nzx), Vector3d( 0.0, 0.0,-1.0))
xbc = PeriodicBoundary3d(xPlane1, xPlane2)
ybc = PeriodicBoundary3d(yPlane1, yPlane2)
zbc = PeriodicBoundary3d(zPlane1, zPlane2)
hydro.appendBoundary(xbc)
hydro.appendBoundary(ybc)
hydro.appendBoundary(zbc)
mhd.appendBoundary(xbc)
mhd.appendBoundary(ybc)
mhd.appendBoundary(zbc)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator3d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(mhd)
integrator.verbose = True
integrator.rigorousBoundaries = True
integrator.lastDt = 1e-3
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(mhd)")
output("integrator.valid()")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
#raw_input()
restartBaseName = '%s/%s'%(restartDir, simName)
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            initializeMassDensity = True,
                            restartBaseName = restartBaseName)
output("control")
#print 'max |div B| (0):', maxDivB0

# Restore if desired.
if restoreCycle is not None:
   if restoreCycle == -1:
      restoreCycle = findLastRestart(simName)
   control.loadRestartFile(restoreCycle)
else:
   dumpPhysicsState(integrator, simName, visitDir, dumpDerivatives = True)

output("integrator.dtGrowth")
# If we're using a projection scheme to clean div B, advance one step and 
# read off our diagnostics.
if mhd.divBCleaner == MHD.BDivergenceCleanerType.GreensFnProjCleaner or \
   mhd.divBCleaner == MHD.BDivergenceCleanerType.BiotSavartProjCleaner:
   control.advance(control.time() + 1e-10, 1)
   maxDivB1 = max(mhd.maxDivB(), abs(mhd.minDivB()))

# Otherwise, go get 'em!
else:
   while control.time() < goalTime:
      dt = goalTime/10
      control.advance(min(goalTime, control.time() + dt), maxSteps)
      control.dropRestartFile()
      dumpPhysicsState(integrator, simName, visitDir, dumpDerivatives = True)
      maxDivB1 = max(mhd.maxDivB(), abs(mhd.minDivB()))
print 'max |div B| (1):', maxDivB1

# Plot the final field configuration (and its divergence).
#plotField(x, B, 'B after div cleaning', 'B-after.png')
#plotField(x, [Bi.x for Bi in B.internalValues()], 'Bx after div cleaning', 'Bx-after.png')
#plotField(x, nodes.magneticDivergence(), 'div B after div cleaning', 'divB-after.png')
