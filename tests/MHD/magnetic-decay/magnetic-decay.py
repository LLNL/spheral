#-------------------------------------------------------------------------------
# A magnetic signal is placed in a resistive medium and allowed to decay.
# This test has 1D, 2D, and 3D configurations that treat the same problem.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralVisitDump import dumpPhysicsState

# Load the mpi module if we"re parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

from GenerateNodeDistribution3d import *
from CubicNodeGenerator import GenerateCubicNodeDistribution

title("Magnetic decay test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            dim = 1,      # Dimensionality of the problem.

            n = 64,    # Number of nodes in a non-symmetric direction

            # Magnetic signal.  Can be 'square' or 'gaussian'.
            config = 'square',

            Kernel = BSplineKernel3d, # SPH kernel
            PiKernel = BSplineKernel3d, # Artificial viscosity kernel

            rho0 = 1.0,
            u0 = 0.0,
            B0 = Vector3d(1.0, 1.0, 1.0),
            nPerh = 1.25,
            mu0 = 1.0,

            gamma = 5.0/3.0,
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
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            HEvolution = Hydro3d.HEvolutionType.IdealH,
            limitIdealH = False,

            neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector3d(0.0, 0.0, 0.0),

            dt = 0.0001,
            dtMin = 1.0e-5,
            dtMax = None,
            dtGrowth = 2.0,
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

# Global neighbors list.  This prevents Spheral from freaking out when 
# its neighbor objects go out of scope after creation.
neighbors = []

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
leftXmin = (-0.5, -0.5*nyx, -0.5*nzx)
leftXmax = (0., 0.5*nyx, 0.5*nzx)
rightXmin = (0., -0.5*nyx, -0.5*nzx)
rightXmax = (0.5, 0.5*nyx, 0.5*nzx)

# Give me a conducting node list.
nodes = ConductingFluidNodeList("conductor", self.eos, self.WT, self.WTPi)
nodes.HsmoothFraction = HsmoothFraction
nodes.XSPH = XSPH
nodes.nodesPerSmoothingScale = nPerh
nodes.epsilonTensile = epsilonTensile
nodes.nTensile = nTensile
nodes.hmin = hmin
nodes.hmax = hmax
nodes.hminratio = hminratio

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

if restoreCycle is None:
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

   # Set the "left" and "right" states within the shock tube.
   x = nodes.positions()
   rho = nodes.massDensity()
   v = nodes.velocity()
   u = nodes.specificThermalEnergy()
   B = nodes.magneticInduction()
   if dim == 1:
      # 1D configuration: left-right boundary is at x = 0.
      for i in xrange(nodes.numInternalNodes):
         if x[i].x <= 0.0:
            rho[i] = self.rhoL
            v[i] = self.vL
            u[i] = self.uL
            B[i] = self.BL
               else:
                  rho[i] = self.rhoR
                  v[i] = self.vR
                  u[i] = self.uR
                  B[i] = self.BR
         elif dim == 2:
            # 2D configuration: left-right boundary is at x = y.
            for i in xrange(nodes.numInternalNodes):
               if x[i].x <= x[i].y:
                  rho[i] = self.rhoL
                  v[i] = self.vL
                  u[i] = self.uL
                  B[i] = self.BL
               else:
                  rho[i] = self.rhoR
                  v[i] = self.vR
                  u[i] = self.uR
                  B[i] = self.BR
         else:
            # 3D configuration: left-right boundary is at x = y = z.
            for i in xrange(nodes.numInternalNodes):
               if x[i].x <= x[i].y and x[i].y <= x[i].z:
                  rho[i] = self.rhoL
                  v[i] = self.vL
                  u[i] = self.uL
                  B[i] = self.BL
               else:
                  rho[i] = self.rhoR
                  v[i] = self.vR
                  u[i] = self.uR
                  B[i] = self.BR

simName = 'magnetic-decay-%s-%ix%ix%i'%(config, nx, ny, nz)
dataDir = simName
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
# Construct our node list.
#-------------------------------------------------------------------------------
ICs = shockTubeICs[config-1]
nodes = ICs.Fluid(dim, nx, ny, nz)
goalTime = ICs.endTime
dtSample = 0.01 * goalTime

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
#q = MonaghanGingoldViscosity3d(1.0, 1.5)
q = PriceMonaghanDissipation(1.0, 1.0, 1.0, 0.75, 1.0)
#q.limiter = Qlimiter
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
# Construct an MHD object.
#-------------------------------------------------------------------------------
MHD = MHD(WT, mu0)

# Trap the flow inside a reflecting box.
e1 = Vector3d(1, 0, 0)
e2 = Vector3d(0, 1, 0)
e3 = Vector3d(0, 0, 1)
transversePlane1 = Plane3d(Vector3d(xmax[0], 0, 0), -e1)
transversePlane2 = Plane3d(Vector3d(xmin[0], 0, 0),  e1)
transversePlane3 = Plane3d(Vector3d(0, xmax[1], 0), -e2)
transversePlane4 = Plane3d(Vector3d(0, xmin[1], 0),  e2)
transversePlane5 = Plane3d(Vector3d(0, 0, xmax[2]), -e3)
transversePlane6 = Plane3d(Vector3d(0, 0, xmin[2]),  e3)
transverseBC1 = ReflectingBoundary3d(transversePlane1)
transverseBC2 = ReflectingBoundary3d(transversePlane2)
transverseBC3 = ReflectingBoundary3d(transversePlane3)
transverseBC4 = ReflectingBoundary3d(transversePlane4)
transverseBC5 = ReflectingBoundary3d(transversePlane5)
transverseBC6 = ReflectingBoundary3d(transversePlane6)
hydro.appendBoundary(transverseBC1)
hydro.appendBoundary(transverseBC2)
hydro.appendBoundary(transverseBC3)
hydro.appendBoundary(transverseBC4)
hydro.appendBoundary(transverseBC5)
hydro.appendBoundary(transverseBC6)
MHD.appendBoundary(transverseBC1)
MHD.appendBoundary(transverseBC2)
MHD.appendBoundary(transverseBC3)
MHD.appendBoundary(transverseBC4)
MHD.appendBoundary(transverseBC5)
MHD.appendBoundary(transverseBC6)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator3d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(MHD)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(MHD)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
#from Spasmos import pause
#pause()
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            initializeMassDensity = True)
output("control")

# Restore if desired.
if restoreCycle is not None:
   control.loadRestartFile(restoreCycle)

#-------------------------------------------------------------------------------
# Restore if necessary.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
   nextGoalTime = min(control.time() + dtSample, goalTime)
   control.advance(nextGoalTime, maxSteps)
   control.dropRestartFile()
   dumpPhysicsState(integrator, simName, visitDir, dumpDerivatives = True)

#-------------------------------------------------------------------------------
# Make some line-outs.
#-------------------------------------------------------------------------------
Nsample = 100
h = 1.0/Nsample
if dim == 1:
   points = [Vector1d(i*h - 0.5) for i in xrange(Nsample)]
elif dim == 2:
   points = [Vector2d(i*h - 0.5, i*h - 0.5) for i in xrange(Nsample)]
else:
   points = [Vector3d(i*h - 0.5, i*h - 0.5, i*h - 0.5) for i in xrange(Nsample)]
rhos = ScalarFieldList3d()
rhos.appendField(nodes.massDensity())
rhoSamples = [rhos(xi, WT) for xi in points]
vs = VectorFieldList3d()
vs.appendField(nodes.velocity())
vSamples = [vs(xi, WT) for xi in points]
us = ScalarFieldList3d()
us.appendField(nodes.specificThermalEnergy())
uSamples = [us(xi, WT) for xi in points]
Bs = VectorFieldList3d()
Bs.appendField(nodes.magneticInduction())
BSamples = [Bs(xi, WT) for xi in points]
Js = VectorFieldList3d()
Js.appendField(nodes.currentDensity())
JSamples = [Js(xi, WT) for xi in points]

# Now lookit the various quantities!
#rhoL = leftNodes.massDensity().internalValues()
#vxL = [u.dot(e1) for u in leftNodes.velocity().internalValues()]
#vyL = [u.dot(e2) for u in leftNodes.velocity().internalValues()]
#vzL = [u.dot(e3) for u in leftNodes.velocity().internalValues()]
#uL = leftNodes.specificThermalEnergy().internalValues()
#BxL = [u.dot(e1) for u in leftNodes.magneticInduction().internalValues()]
#ByL = [u.dot(e2) for u in leftNodes.magneticInduction().internalValues()]
#BzL = [u.dot(e3) for u in leftNodes.magneticInduction().internalValues()]
#
#rhoR = rightNodes.massDensity().internalValues()
#vxR = [u.dot(e1) for u in rightNodes.velocity().internalValues()]
#vyR = [u.dot(e2) for u in rightNodes.velocity().internalValues()]
#vzR = [u.dot(e3) for u in rightNodes.velocity().internalValues()]
#uR = rightNodes.specificThermalEnergy().internalValues()
#BxR = [u.dot(e1) for u in rightNodes.magneticInduction().internalValues()]
#ByR = [u.dot(e2) for u in rightNodes.magneticInduction().internalValues()]
#BzR = [u.dot(e3) for u in rightNodes.magneticInduction().internalValues()]

