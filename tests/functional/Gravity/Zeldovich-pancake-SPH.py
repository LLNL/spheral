#-------------------------------------------------------------------------------
# The gravitational field of a uniform mass distribution (in 3D).
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralVisitDump import dumpPhysicsState
import mpi

from GenerateNodeDistribution3d import *
from CubicNodeGenerator import GenerateCubicNodeDistribution

title("3-D Zeldovich pancake collapse test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList3d,

            seed = "lattice",

            nx = 11,
            ny = 11,
            nz = 11,

            
            M = 1.0,     # point mass
            meps = 1e-8, # Negligible masses around the point mass

            rho1 = 1.0,
            eps1 = 0.0,
            vr1 = -1.0,
            nPerh = 1.25,

            gamma = 5.0/3.0,
            mu = 1.0,
            z1 = 50., # Starting redshift.
            z2 = 5., # Ending redshift.
            zc = 5., # Caustic redshift.
            #Qconstructor = MonaghanGingoldViscosity3d,
            Qconstructor = TensorMonaghanGingoldViscosity3d,
            Cl = 1.0,
            Cq = 0.75,
            Qlimiter = True,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 1.0,
            hminratio = 0.05,
            HsmoothFraction = 0.0,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            HEvolution = Hydro3d.HEvolutionType.IdealH,
            limitIdealH = False,

            goalTime = 0.6,
            dt = 0.0001,
            dtMin = 1.0e-5,
            dtMax = None,
            dtGrowth = 2.0,
            dtSample = 0.1,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            sumForMassDensity = Hydro3d.MassDensityType.RigorousSumDensity,

            restoreCycle = None,

            L1v0 =   0.0889732,
            L1rho0 = 5.51975,
            L1eps0 = 0.04701,
            L1P0 =   1.66301,
            L1A0 =   0.00344783,

            graphics = False,
            )

# For this test, nx, ny, and nz must all be odd.
if (2*(nx/2) == nx) or (2*(ny/2) == ny) or (2*(nz/2) == nz):
   raise ValueError, 'nx, ny, and nz must all be odd.'

#-------------------------------------------------------------------------------
# If we're using the cubic node generator, then scale things so we get a
# constant work per domain, and run to the same self-similar shock fraction
# of the node distribution.
#-------------------------------------------------------------------------------
if seed == "cubic":
    nxdomains = int(mpi.procs**(1.0/3.0) + 0.1)
    assert nxdomains**3 == mpi.procs
    nx *= nxdomains
    ny *= nxdomains
    nz *= nxdomains
    print nxdomains, nx, ny, nz

#-------------------------------------------------------------------------------
# A few derived variables.
#-------------------------------------------------------------------------------
xmin = (-0.5, -0.5, -0.5)
xmax = (0.5, 0.5, 0.5)

dataDir = "point-masses-SPH-%ix%ix%i" % (nx, ny, nz)
visitDir = dataDir + "/visit"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS3d(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 1000)
WTPi = TableKernel3d(BSplineKernel3d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = NodeListConstructor("nodes", eos, WT, WTPi)
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
neighbor1 = TreeNeighbor3d(nodes,
                           kernelExtent = kernelExtent)
nodes.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    if seed == "cubic":
        from DistributeNodes import nullDistributeNodes3d
        generator = GenerateCubicNodeDistribution(nx, ny, nz, rho1,
                                                  xmin = xmin,
                                                  xmax = xmax,
                                                  nNodePerh = nPerh,
                                                  SPH = (NodeListConstructor == SphNodeList3d))
        nullDistributeNodes3d((nodes, generator))
    else:
        from ParMETISDistributeNodes import distributeNodes3d
        generator = GenerateNodeDistribution3d(nx, ny, nz, rho1, seed,
                                               xmin = xmin,
                                               xmax = xmax,
                                               nNodePerh = nPerh,
                                               SPH = (NodeListConstructor == SphNodeList3d))
        distributeNodes3d((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Back out the start and end times and expansion factors from the start 
    # and end redshifts z1 and z2.  In a Friedmann cosmology the expansion 
    # factor a(t) is proportional to t**2/3 in the matter-dominated phase.
    a1 = 1.0/(1+z1)
    t1 = a1**1.5
    a1dot = (2./3)*t1**(-1./3)
    a2 = 1.0/(1+z2)
    t2 = a2**1.5

    # Set node specific thermal energies.
    nodes.specificThermalEnergy(ScalarField3d("tmp", nodes, eps1))

    # Set node positions and velocities.
    A = -5*(1+zc)/2
    for nodeID in xrange(nodes.numNodes):
        q = nodes.positions()[nodeID].x
        nodes.positions()[nodeID] += Vector3d(0.4*a1*A*sin(q), 0, 0)
        nodes.velocity()[nodeID] += Vector3d(0.4*a1dot*A*sin(q), 0, 0)

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
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.negligibleSoundSpeed")
output("q.csMultiplier")
output("q.balsaraShearCorrection")

##-------------------------------------------------------------------------------
## Construct the hydro physics object.
##-------------------------------------------------------------------------------
hydro = Hydro3d(WT, WTPi, q)
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
# Construct the gravity physics object.
#-------------------------------------------------------------------------------
#from Spasmos import pause
#pause('Simulation paused...')
G = 1.0
gravity = SPHGravity3d(WT, G, 2.0)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane1 = Plane3d(Vector3d(-0.5, 0.0, 0.0), Vector3d( 1.0, 0.0, 0.0))
xPlane2 = Plane3d(Vector3d( 0.5, 0.0, 0.0), Vector3d(-1.0, 0.0, 0.0))
yPlane1 = Plane3d(Vector3d( 0.0,-0.5, 0.0), Vector3d( 0.0, 1.0, 0.0))
yPlane2 = Plane3d(Vector3d( 0.0, 0.5, 0.0), Vector3d( 0.0,-1.0, 0.0))
zPlane1 = Plane3d(Vector3d( 0.0, 0.0,-0.5), Vector3d( 0.0, 0.0, 1.0))
zPlane2 = Plane3d(Vector3d( 0.0, 0.0, 0.5), Vector3d( 0.0, 0.0,-1.0))
xbc = PeriodicBoundary3d(xPlane1, xPlane2)
ybc = PeriodicBoundary3d(yPlane1, yPlane2)
zbc = PeriodicBoundary3d(zPlane1, zPlane2)
hydro.appendBoundary(xbc)
hydro.appendBoundary(ybc)
hydro.appendBoundary(zbc)
gravity.appendBoundary(xbc)
gravity.appendBoundary(ybc)
gravity.appendBoundary(zbc)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator3d(db)
integrator.appendPhysicsPackage(gravity)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(gravity)")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            initializeMassDensity = True)
integrator.setCurrentTime(t1)
output("control")

# Initial dump.
dumpPhysicsState(integrator,
                 "Zeldovich-pancake-SPH-%ix%ix%i-visit" % (nx, ny, nz),
                 visitDir)

#-------------------------------------------------------------------------------
# Advance.
#-------------------------------------------------------------------------------
#nextGoalTime = min(control.time() + dtSample, goalTime)
print 'Running from t = %g to %g'%(t1, t2)
control.advance(t2, 1000)
dumpPhysicsState(integrator,
                 "Zeldovich-pancake-SPH-%ix%ix%i-visit" % (nx, ny, nz),
                 visitDir)
D = gravity.matrix()
psi = gravity.potential()[0]
acc = nodes.DvelocityDt()

