#-------------------------------------------------------------------------------
# The gravitational field of a spherical mass distribution (in 3D).
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralVisitDump import dumpPhysicsState
import mpi

from GenerateNodeDistribution3d import *
from CubicNodeGenerator import GenerateCubicNodeDistribution

title("3-D spherical gravity test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList3d,

            seed = "lattice",

            nx = 10,
            ny = 10,
            nz = 10,

            rho1 = 1.0,
            eps1 = 0.0,
            vr1 = -1.0,
            nPerh = 1.25,

            rmin = 0.0,
            rmax = 1.0,

            gamma = 5.0/3.0,
            mu = 1.0,
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
xmin = (0.0, 0.0, 0.0)
xmax = (1.0, 1.0, 1.0)

dataDir = "spherical-3d-SPH-%ix%ix%i" % (nx, ny, nz)
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
                                               rmin = rmin,
                                               rmax = rmax,
                                               nNodePerh = nPerh,
                                               SPH = (NodeListConstructor == SphNodeList3d))
        distributeNodes3d((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes.specificThermalEnergy(ScalarField3d("tmp", nodes, eps1))

    # Set node velocities
    for nodeID in xrange(nodes.numNodes):
        nodes.velocity()[nodeID] = nodes.positions()[nodeID].unitVector()*vr1

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
#hydro = Hydro3d(WT, WTPi, q)
#hydro.cfl = cfl
#hydro.HEvolution = HEvolution
#hydro.sumForMassDensity = sumForMassDensity
#hydro.HsmoothMin = hmin
#hydro.HsmoothMax = hmax
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
G = 1.0
gravity = SPHGravity3d(WT, G, 2.0)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
#xPlane0 = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, 0.0, 0.0))
#yPlane0 = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0))
#zPlane0 = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 1.0))
#xbc0 = ReflectingBoundary3d(xPlane0)
#ybc0 = ReflectingBoundary3d(yPlane0)
#zbc0 = ReflectingBoundary3d(zPlane0)
#gravity.appendBoundary(xbc0)
#gravity.appendBoundary(ybc0)
#gravity.appendBoundary(zbc0)
#output("gravity.haveBoundary(xbc0)")
#output("gravity.haveBoundary(ybc0)")
#output("gravity.haveBoundary(zbc0)")

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator3d(db)
integrator.appendPhysicsPackage(gravity)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(gravity)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep)
output("control")

#-------------------------------------------------------------------------------
# Advance one step.
#-------------------------------------------------------------------------------
nextGoalTime = min(control.time() + dtSample, goalTime)
control.step()
D = gravity.matrix()
Dij = D.entries()
Ds = D.structure
Dsij = Ds.entries()
psi = gravity.potential

stop

#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
import NohAnalyticSolution
answer = NohAnalyticSolution.NohSolution(3, h0 = nPerh*rmax/nx)
if graphics:
    from SpheralGnuPlotUtilities import *

    # Plot the final state.
    rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotRadialState(db)

    # Overplot the analytic solution.
    plotAnswer(answer, control.time(),
               rhoPlot = rhoPlot,
               velPlot = vrPlot,
               epsPlot = epsPlot,
               PPlot = PPlot)

#-------------------------------------------------------------------------------
# Measure the difference between the simulation and analytic answer.
#-------------------------------------------------------------------------------
# Figure out which of the node we want to measure the error on.
rmin = 0.05
rmax = 0.35
rall = [x.magnitude() for x in nodes.positions().internalValues()]
imask = [i for i in xrange(nodes.numInternalNodes)
         if (rall[i] > rmin and rall[i]  < rmax)]
Nlocal = len(imask)
Nglobal = mpi.allreduce(Nlocal, mpi.SUM)

# Find the local profiles.
r = nodes.positions().internalValues()
vel = nodes.velocity().internalValues()
rho = nodes.massDensity().internalValues()
eps = nodes.specificThermalEnergy().internalValues()
P = nodes.pressure().internalValues()
xprof = [rall[i] for i in imask]
vprof = [vel[i].dot(r[i].unitVector()) for i in imask]
rhoprof = [rho[i] for i in imask]
epsprof = [eps[i] for i in imask]
Pprof = [P[i] for i in imask]
Aprof = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pprof, rhoprof)]

# Compute the analytic answer on the positions of the nodes.
xans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]

# Compute the L1 error norms.
def computeL1(q, q0):
    if Nlocal > 0:
        import Pnorm
        error = [qi - q0i for (qi, q0i) in zip(q, q0)]
        Pn = Pnorm.Pnorm(error, xprof)
        L1local = Pn.pnormAverage(1, rmin = 0.05, rmax = 0.35)
    else:
        L1local = 0.0
    return mpi.allreduce(L1local*Nlocal, mpi.SUM)/Nglobal
L1v = computeL1(vprof, vans)
L1rho = computeL1(rhoprof, rhoans)
L1eps = computeL1(epsprof, epsans)
L1P = computeL1(Pprof, Pans)
L1A = computeL1(Aprof, Aans)

# Now we can compare the errors to our expectations.
for (name, L1, L10) in [("Velocity    ", L1v, L1v0),
                        ("Mass Density", L1rho, L1rho0),
                        ("Thermal E   ", L1eps, L1eps0),
                        ("Pressure    ", L1P, L1P0),
                        ("Entropy     ", L1A, L1A0)]:
    L1exp = L10 * (nx/25.0)**(-0.8)
    print "\t%s L1 = %g < %g" % (name, L1, L1exp)
    if L1 > L1exp:
        raise "L1 error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                              L1,
                                                                              L1exp)
