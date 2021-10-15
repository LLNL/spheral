#-------------------------------------------------------------------------------
# A magnetosonic wave test.  Here we propogate an MHD wave round and
# rount in a periodic box.  This specific example is based on the test case
# described in D.J. Price's dissertation.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *

import loadmpi
mpi, proc, numProcs = loadmpi.loadmpi()

from GenerateNodeDistribution3d import *
from CubicNodeGenerator import GenerateCubicNodeDistribution

title("Magnetosonic wave propagation test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(n = 100,    # Resolution

            dim = 1,    # Dimensionality of the propagation.
            waveType = 'fast', # Type of magnetosonic wave.
            mu = 1.0,
            nPerh = 2.01,
            Kernel = BSplineKernel3d, # Hydro kernel
            PiKernel = BSplineKernel3d, # Artificial viscosity kernel

            seed = 'lattice',
            Qconstructor = PriceMonaghanDissipation,
            #Qconstructor = TensorMonaghanGingoldViscosity3d,
            Cl = 0.0,
            Cq = 0.0,
            Qlimiter = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.3,
            nTensile = 4,

            neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector3d(0, 0, 0),

            IntegratorConstructor = SynchronousRK2Integrator3d,
            goalTime = 5.0,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = True,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = Hydro3d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro3d.MassDensityType.RigorousSumDensity,
            compatibleEnergy = False,
            gradhCorrection = True,

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "dumps-MagnetosonicWave",

            graphics = "gnu",
            )

# Interpolation kernels.
WT = TableKernel3d(Kernel(), 1000)
WTPi = TableKernel3d(PiKernel(), 1000)
kernelExtent = WT.kernelExtent()
output("WT")
output("WTPi")

# Figure out how many nodes should be used in a direction of symmetry.
nsym = 2 * 2 * kernelExtent * nPerh

# Use the dimensionality of the problem to reassess the numbers of nodes on 
# each side of our box.
if dim == 1:
   nx = n
   ny = nz = nsym
elif dim == 2:
   nx = ny = n
   nz = nsym
elif dim == 3:
   nx = ny = nz = n
else:
   raise ValueError, 'Invalid dimension: %d'%dim

# Set the dimensions of the computational bounding box, accounting for 
# symmetry.
nyx = 1.0*ny/nx
nzx = 1.0*nz/nx
xmin = (-0.5, -0.5*nyx, -0.5*nzx)
xmax = ( 0.5,  0.5*nyx,  0.5*nzx)

# Set up the magnetic parameters based upon the wave type.
if waveType == 'fast':
    eos = GammaLawGasMKS3d(5.0/3.0, mu)
    rho0 = 1.004
    rho1 = 0.0055 * rho0
    B0 = Vector3d(0.5, 0.5, 0.5)
    u0 = 0.3
    rho1 = 0.05 * rho0
elif waveType == 'slow':
    eos = GammaLawGasMKS3d(5.0/3.0, mu)
    rho0 = 1.004
    B0 = Vector3d(1.415, 1.415, 1.415)
    options.u0 = 4.5
elif waveType == 'isothermal':
    cs2 = 1.0
    eos = IsothermalEquationOfStateMKS3d(cs2, mu)
    rho0 = 1.0
    rho1 = 0.005 * rho0
    B0 = Vector3d(2, 0, 0)
    u0 = 1.0
else:
   raise ValueError, 'Unknown wave type: \'%s\''%waveType
print 'Running %s wave...'%waveType

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = ConductingFluidNodeList("fluid", eos, WT, WTPi)
nodes.XSPH = XSPH
nodes.hmin = hmin
nodes.hmax = hmax
nodes.nodesPerSmoothingScale = nPerh
nodes.epsilonTensile = epsilonTensile
nodes.nTensile = nTensile
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")
output("nodes.epsilonTensile")
output("nodes.nTensile")
output("nodes.XSPH")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor = NestedGridNeighbor3d(nodes,
                                neighborSearchType,
                                numGridLevels,
                                topGridCellSize,
                                origin,
                                kernelExtent)
nodes.registerNeighbor(neighbor)

# Set up the nodes.
if seed == "cubic":
    from DistributeNodes import nullDistributeNodes3d
    generator = GenerateCubicNodeDistribution(nx, ny, nz, rho0,
                                              xmin = xmin,
                                              xmax = xmax,
                                              nNodePerh = nPerh,
                                              SPH = True)
    nullDistributeNodes3d((nodes, generator))
else:
    from ParMETISDistributeNodes import distributeNodes3d
    generator = GenerateNodeDistribution3d(nx, ny, nz, 1.0, seed,
                                           xmin = xmin,
                                           xmax = xmax,
                                           nNodePerh = nPerh,
                                           SPH = True)
    distributeNodes3d((nodes, generator))

# Tag the nodes according to their position projected along the direction of 
# the signal propagation.
# FIXME: do this!

# Find the cumulative mass at each point.
Mi = ScalarField3d("Cumulative mass", nodes)
positions = mpi.allreduce([(nodes.positions()[i].x, i, mpi.rank)
                           for i in xrange(nodes.numInternalNodes)], mpi.SUM)
assert len(positions) == nx
positions.sort()
Msum = 0.0
mi = rho0/nx
for (x, i, proc) in positions:
    Msum += mi
    if proc == mpi.rank:
        assert i < nodes.numInternalNodes
        Mi[i] = Msum
assert fuzzyEqual(Msum, rho0)

# Define the function which we are going to solve for the node positions.
twopi = 2.0*pi
class MassFunctor:
    def __init__(self, Mcumulative):
        self.Mcumulative = Mcumulative
        return
    def __call__(self, x):
        return (self.Mcumulative - rho0*(x + rho1/twopi*(1.0 - cos(twopi*x))),
                -rho0*(1.0 + rho1*sin(twopi*x)))

# Set the node positions, velocities, and densities.
import newtonRaphson
cs = sqrt(cs2)
for i in xrange(nodes.numInternalNodes):
    func = MassFunctor(Mi[i] - 0.5*mi)
    xi = newtonRaphson.newtonRaphson(func, 0.0, 1.0)
    nodes.positions()[i].x = xi
    nodes.velocity()[i].x = A*cs*sin(twopi*xi)
    rhoi = rho0*(1.0 + rho1*sin(twopi*xi))
    nodes.massDensity()[i] = rhoi
    nodes.magneticInduction()[i] = B0
    nodes.specificThermalEnergy()[i] = u0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro3d(WT, WTPi, q, compatibleEnergy, gradhCorrection)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct the MHD physics object.
#-------------------------------------------------------------------------------
MHD = MHD(WT, mu0)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
e1 = Vector3d(1, 0, 0)
e2 = Vector3d(0, 1, 0)
e3 = Vector3d(0, 0, 1)
xPlane0 = Plane3d(Vector3d(xmin[0], 0, 0),  e1)
xPlane1 = Plane3d(Vector3d(xmax[0], 0, 0), -e1)
yPlane0 = Plane3d(Vector3d(0, xmin[1], 0),  e2)
yPlane1 = Plane3d(Vector3d(0, xmax[1], 0), -e2)
zPlane0 = Plane3d(Vector3d(0, 0, xmin[2]),  e3)
zPlane1 = Plane3d(Vector3d(0, 0, xmax[2]), -e3)
xbc = PeriodicBoundary3d(xPlane0, xPlane1)
ybc = PeriodicBoundary3d(yPlane0, yPlane1)
zbc = PeriodicBoundary3d(zPlane0, zPlane1)
hydro.appendBoundary(xbc)
hydro.appendBoundary(ybc)
hydro.appendBoundary(zbc)
MHD.appendBoundary(xbc)
MHD.appendBoundary(ybc)
MHD.appendBoundary(zbc)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(MHD)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = (sumForMassDensity != Hydro3d.MassDensityType.IntegrateDensity))
output("control")

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH()
    control.smoothState(smoothIters)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if control.time() < goalTime:
    control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# Compute the analytic answer.
#-------------------------------------------------------------------------------
import MagnetosonicWaveSolution
xlocal = [pos.x for pos in nodes.positions().internalValues()]
xglobal = mpi.reduce(xlocal, mpi.SUM)
dx = 1.0/nx
h1 = 1.0/(nPerh*dx)
answer = MagnetosonicWaveSolution.MagnetosonicWaveSolution(eos, rho0, rho1, u0, twopi, h0, waveType)

### Compute the simulated specific entropy.
##rho = mpi.allreduce(nodes.massDensity().internalValues(), mpi.SUM)
##P = mpi.allreduce(nodes.pressure().internalValues(), mpi.SUM)
##A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

### The analytic solution for the simulated entropy.
##xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xglobal)
##Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print "Total energy error: %g" % Eerror
if abs(Eerror) > 1e-15:
    raise "Energy error outside allowed bounds."
