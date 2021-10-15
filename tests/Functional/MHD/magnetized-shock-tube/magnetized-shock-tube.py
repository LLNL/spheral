#-------------------------------------------------------------------------------
# A magnetized shock tube in one of several fabulous configurations!
# This test has 1D, 2D, and 3D configurations.  In the 1D case, the shock tube 
# is aligned with the x direction, and reflecting BCs are used in the y and z 
# directions to provide symmetry.  The 2D case, the shock tube connects the 
# opposing corners of a square in the x-y plane, and reflecting BCs are 
# provided for z.  In the 3D case, the shock tube connects the opposing corners 
# of a cube.
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

title("Magnetized shock tube test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            dim = 1,      # Dimensionality of the problem.

            n = 64,    # Number of nodes in a non-symmetric direction

            # Shock tube initial configuration.
            # The tests are those suite mentioned by 
            # Dai and Woodward (DW, JCP 111, 354, 1994), 
            # Ryu and Jones (RJ, ApJ 442, 228, 1995), and 
            # Brio and Wu (BW, JCP 75, 500, 1988) and are as follows: 
            # 1. Magnetized Noh problem (DW), 
            # 2. 2D discontinuities (DW), 
            # 3. Weak 3D discontinuities (DW), 
            # 4. 3D discontinuities (RJ), 
            # 5. Vt = Bn = 0 discontinuities (DW),
            # 6. Magnetosonic rarefactions (RJ), 
            # 7. Switch-on fast shock (RJ), 
            # 8. Switch-off fast rarefaction (RJ), 
            # 9. Switch-off slow shock (RJ), 
            # 10. Switch-on slow rarefaction (RJ), 
            # 11. Slow compound structures (BW, gamma = 2), 
            # 12. Slow compound structures (RJ, gamma = 5/3), 
            # 13. Fast compound structures (RJ)
            config = 11,
            goalTime = None,
            samples = None,
            restartsPerSample = 1,

            Kernel = BSplineKernel3d, # Hydro kernel
            PiKernel = BSplineKernel3d, # Artificial viscosity kernel

            rigorousBoundaries = True,
            nPerh = 1.3,
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

# Here's a node generator for our shock tube.
class ShockTubeGenerator(NodeGeneratorBase):
   def __init__(self, n, dim, WT, rhoL, rhoR, nNodePerh = 2.01):

      # Figure out how many nodes should be used in a direction of symmetry.
      kernelExtent = WT.kernelExtent()
      nsym = int(2 * 2 * kernelExtent * nPerh)

      # Use the dimensionality of the problem to reassess the numbers of nodes on 
      # each side of our box.
      self.dim = dim
      if dim == 1:
         self.nx = n
         self.ny = self.nz = nsym
      elif dim == 2:
         self.nx = self.ny = n
         self.nz = nsym
      elif dim == 3:
         self.nx = self.ny = self.nz = n
      else:
         raise ValueError, 'Invalid dimension: %d'%dim

      # Set the dimensions of the computational bounding box, accounting for 
      # symmetry.
      nyx = 1.0*self.ny/self.nx
      nzx = 1.0*self.nz/self.nx
      xmin = (-0.5, -0.5*nyx, -0.5*nzx)
      xmax = ( 0.5,  0.5*nyx,  0.5*nzx)
      self.rhoL = rhoL
      self.rhoR = rhoR

      # Now make our lattice.
      if dim == 1:
         self.lattice1D(xmin, xmax, nNodePerh)
      elif dim == 2:
         self.lattice2D(xmin, xmax, nNodePerh)
      elif dim == 3:
         self.lattice3D(xmin, xmax, nNodePerh)

   def lattice1D(self, xmin, xmax, nNodePerh):
      # The shock tube is divided down the middle along the x axis.
      cubeRootRhos = pow(self.rhoL/self.rhoR, 0.333333333)
      nxL = int(self.nx * cubeRootRhos / (1 + cubeRootRhos))
      nxR = self.nx - nxL

      # Compute the left and right bounding boxen.
      xminL = xmin
      xmaxL = (0., xmax[1], xmax[2])
      xminR = (0., xmin[1], xmin[2])
      xmaxR = xmax

      # Compute the resolution in the transverse dimensions.  Here we assume that 
      # the denser region is always the left side.
      assert(self.ny == self.nz)
      nyR = nzR = self.ny
      nyL = nzL = int(sqrt(self.rhoL*nxR/(self.rhoR*nxL))*nyR)

      # Set 'em up!
      genL = GenerateNodeDistribution3d(nxL, nyL, nzL, self.rhoL, 
                                        distributionType = 'lattice',
                                        xmin = xminL,
                                        xmax = xmaxL,
                                        nNodePerh = nNodePerh,
                                        SPH = True)
      genR = GenerateNodeDistribution3d(nxR, nyR, nzR, self.rhoR, 
                                        distributionType = 'lattice',
                                        xmin = xminR,
                                        xmax = xmaxR,
                                        nNodePerh = nNodePerh,
                                        SPH = True)
      self.x = genL.x[:] + genR.x[:]
      self.y = genL.y[:] + genR.y[:]
      self.z = genL.z[:] + genR.z[:]
      self.m = genL.m[:] + genR.m[:]
      self.H = genL.H[:] + genR.H[:]
      self.xmin = (xminL[0], xmin[1], xmin[2])
      self.xmax = (xmaxR[0], xmax[1], xmax[2])

   def lattice2D(self, xmin, xmax, nNodePerh, SPH):
      pass

   def lattice3D(self, xmin, xmax, nNodePerh, SPH):
      pass

   def localPosition(self, i):
      assert i >= 0 and i < len(self.x)
      assert len(self.x) == len(self.y)
      return Vector3d(self.x[i], self.y[i], self.z[i])

   def localMass(self, i):
      assert i >= 0 and i < len(self.m)
      return self.m[i]

   def localMassDensity(self, i):
      assert i >= 0 and i < len(self.x)
      if dim == 1:
         if self.x[i] < 0.0: 
            return self.rhoL
         else:
            return self.rhoR
      elif dim == 2:
         if self.x[i] < self.y[i]:
            return self.rhoL
         else: 
            return self.rhoR
      elif dim == 3:
         if self.x[i] < (self.y[i] + self.z[i]):
            return self.rhoL
         else: 
            return self.rhoR

   def localHtensor(self, i):
      assert i >= 0 and i < len(self.H)
      return self.H[i]


# Shock tube initial conditions. 
class ShockTubeIC(object):
   """ShockTubeIC(name, eos, leftState, rightState,
                  Bx, endTime, solution = None)
   Creates a set of initial conditions for a magnetized shock tube
   problem.  leftState and rightState are tuples containing 
   (rho, P, Vx, Vy, Vz, By, Bz) for their respective sides."""

   def __init__(self, name, eos, leftState, rightState, Bx, 
                endTime, WT, WTPi, solution = None):
      self.name = name
      self.eos = eos
      (self.rhoL, self.pL, VxL, VyL, VzL, ByL, BzL) = leftState
      (self.rhoR, self.pR, VxR, VyR, VzR, ByR, BzR) = rightState
      self.vL = Vector3d(VxL, VyL, VzL)
      self.vR = Vector3d(VxR, VyR, VzR)
      self.BL = Vector3d(Bx, ByL, BzL)
      self.BR = Vector3d(Bx, ByR, BzR)
      self.endTime = endTime
      self.BCs = []
      self.WT = WT
      self.WTPi = WTPi
      self.solution = solution

      # Specific thermal energies.
      if type(eos) is GammaLawGasMKS3d:
         self.uL = self.pL / ((eos.getGamma() - 1)*self.rhoL)
         self.uR = self.pR / ((eos.getGamma() - 1)*self.rhoR)
      elif type(eos) is IsothermalEquationOfState:
         self.uL = self.uR = 0.0
      else:
         assert(False)

   def Fluid(self, n, dim):#, smooth = False, equalSpacing = False):
      """Initializes a fluid in the shock tube."""

      # Give me a conducting node list.
      nodes = ConductingFluidNodeList("fluid", self.eos, self.WT, self.WTPi)
      nodes.HsmoothFraction = HsmoothFraction
      nodes.XSPH = XSPH
      nodes.nodesPerSmoothingScale = nPerh
      nodes.epsilonTensile = epsilonTensile
      nodes.nTensile = nTensile
      nodes.hmin = hmin
      nodes.hmax = hmax
      nodes.hminratio = hminratio

      #-------------------------------------------------------------------------------
      # Construct the neighbor objects.
      #-------------------------------------------------------------------------------
      neighbor = NestedGridNeighbor3d(nodes,
                                      neighborSearchType,
                                      numGridLevels,
                                      topGridCellSize,
                                      origin,
                                      kernelExtent)
      nodes.registerNeighbor(neighbor)

      # Stick this neighbor into the global list to prevent its garbage 
      # collection.
      neighbors.append(neighbor)

      if restoreCycle is None:
#         mpi.synchronizeQueuedOutput(None, None)
         gen = ShockTubeGenerator(n, dim, self.WT, self.rhoL, self.rhoR, nNodePerh = nPerh)

         # Distribute the nodes.
         from ParMETISDistributeNodes import distributeNodes3d
         distributeNodes3d((nodes, gen))

         # Set the "left" and "right" states within the shock tube.
         x = nodes.positions()
         rho = nodes.massDensity()
         v = nodes.velocity()
         u = nodes.specificThermalEnergy()
         B = nodes.magneticInduction()
         if dim == 1:
            # 1D configuration: left-right boundary is at x = 0.
            for i in xrange(nodes.numNodes):
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
      else:
         # We still need nx, ny, nz!
         gen = ShockTubeGenerator(n, dim, self.WT, self.rhoL, self.rhoR, nNodePerh = nPerh)

      # Specify boundary conditions based on the left and right velocities. 
#      self.BCs = []
#      if v[0].x > 0.0:
#         self.BCs.append(SPH.InflowBC(leftPlane, nf))
#      elif v[0].x < 0.0:
#         self.BCs.append(SPH.OutflowBC(leftPlane, nf))
#      else:
#         self.BCs.append(SPH.ReflectingBC(leftPlane, nf))
#      if v[-1].x < 0.0:
#         self.BCs.append(SPH.InflowBC(rightPlane, nf))
#      elif v[-1].x > 0.0:
#         self.BCs.append(SPH.OutflowBC(rightPlane, nf))
#      else:
#         self.BCs.append(SPH.ReflectingBC(rightPlane, nf))

      # Smooth the initial conditions if requested.
#      if smooth:
#         print 'Smoothing initial conditions...'
#         # Determine the number of smoothed nodes.
#         dxL = (xC - xL).magnitude() / NL
#         dxR = (xR - xC).magnitude() / NR
#         d = 0.5 * dxR
#         Ns = int(4 * dxR/dxL)
#         for i in xrange(NL-Ns, N):
#            # Space the smoothed nodes appropriately.
#            x[i] = x[i-2] + X.Vector(2.0*self.rhoR*dxR/rho[i-1])
#
#            # Smooth the magnetic field, the density, and the 
#            # specific thermal energy.
#            xi = min(x[i].x/d, 10.0)
#            expxd = exp(xi)
#            B[i] = (self.BL + self.BR * expxd)/(1 + expxd)
#            rho[i] = (self.rhoL + self.rhoR * expxd)/(1 + expxd)
#            u[i] = (self.uL + self.uR * expxd)/(1 + expxd)
#
#         # Assign H to the inverse density.
#         for i in xrange(N):
#            H[i] = X.SymTensor2(rho[i]/(1.5*m[i]))
#
#         # Reset the position of the right boundary plane so that oscillations
#         # don't crop up there from the varied spacing in the nodes on the 
#         # right side of the problem.  This isn't quite enough to fix the 
#         # problem, but it should stave it off a bit.
#         for BC in self.BCs:
#            if BC.surface is rightPlane:
#               BC.surface.center = x[-1] + X.Vector(0.5*dxR)

      return nodes, gen.nx, gen.ny, gen.nz, gen.xmin, gen.xmax

#WT = TableKernel3d(BSplineKernel3d(), 1000)
#WTPi = TableKernel3d(BSplineKernel3d(), 1000)

# Configurations for this problem.
sqrt4pi = sqrt(4*pi)
gamma53 = GammaLawGasMKS3d(5.0/3.0, mu)
gamma2 = GammaLawGasMKS3d(2.0, mu)
shockTubeICs = \
   [ ShockTubeIC('Magnetized Inflow problem (DW)', gamma53,
                 (1., 20., 10., 0., 0., 1.4105, 0.), 
                 (1., 1., -10., 0., 0., 1.4105, 0.),
                 1.4105, 0.08, WT, WTPi),
     ShockTubeIC('2D discontinuities (DW)', gamma53,
                 (1., 1., 0., 0., 0., 1.4105, 0.), 
                 (0.1, 10., 0., 0., 0., 0.56419, 0.),
                 0.84628, 0.03, WT, WTPi),
     ShockTubeIC('3D weak discontinuities (DW)', gamma53,
                 (1.08, 0.95, 1.2, 0.01, 0.5, 1.0155, 0.56419), 
                 (1., 1., 0., 0., 0., 1.1284, 0.56419),
                 0.56419, 0.2, WT, WTPi),
     ShockTubeIC('3D discontinuities (RJ)', gamma53,
                 (1., 1., 0., 0., 0., 1.6926, 0.), 
                 (0.1, 10., 0., 2., 1., 0.28209, 0.),
                 0.84628, 0.035, WT, WTPi),
     ShockTubeIC('vt = Bn = 0 discontinuities (DW)', gamma53,
                 (0.1, 0.4, 50., 0., 0., -0.28209, -0.56419), 
                 (0.1, 0.2, 0., 0., 0., 0.28209, 0.56419),
                 0., 0.01, WT, WTPi),
     ShockTubeIC('Magnetosonic rarefactions (RJ)', gamma53,
                 (1, 1, -1, 0, 0, 1, 0), 
                 (1, 1, 1, 0, 0, 1, 0), 
                 0, 0.1, WT, WTPi),
     ShockTubeIC('Switch-on fast shock (RJ)', gamma53,
                 (1, 1, 0, 0, 0, 1, 0), 
                 (0.2, 0.1, 0, 0, 0, 0, 0), 
                 1, 0.15, WT, WTPi),
     ShockTubeIC('Switch-off fast rarefaction (RJ)', gamma53,
                 (0.4, 0.52467, -0.66991, 0.98263, 0.0, 0.0025293, 0.0), 
                 (1, 1, 0, 0, 0, 1, 0), 
                 1.3, 0.15, WT, WTPi),
     ShockTubeIC('Switch-off slow shock (RJ)', gamma53,
                 (0.65, 0.5, 0.667, -0.257, 0, 0.55, 0), 
                 (1, 0.75, 0.4, -0.94, 0, 0, 0), 
                 0.75, 0.15, WT, WTPi),
     ShockTubeIC('Switch-on slow rarefaction (RJ)', gamma53,
                 (1, 1, 0, 0, 0, 0, 0), 
                 (0.3, 0.2, 0, 0, 1, 1, 0), 
                 0.7, 0.16, WT, WTPi),
     ShockTubeIC('Slow compound structures (BW, gamma = 2)', gamma2,
                 (1., 1., 0., 0., 0., 1., 0.), 
                 (0.125, 0.1, 0., 0., 0., -1., 0.),
                 0.75, 0.1, WT, WTPi),
     ShockTubeIC('Slow compound structures (RJ, gamma = 5/3)', gamma53,
                 (1., 1., 0., 0., 0., 1., 0.), 
                 (0.125, 0.1, 0., 0., 0., -1., 0.),
                 0.75, 0.1, WT, WTPi),
     ShockTubeIC('Fast compound structures (RJ, gamma = 5/3)', gamma53,
                 (1., 1., 0., 0., 0., 1., 0.), 
                 (0.4, 0.4, 0., 0., 0., -1., 0.),
                 1.3, 0.16, WT, WTPi),
   ]

#-------------------------------------------------------------------------------
# Construct our node list.
#-------------------------------------------------------------------------------
ICs = shockTubeICs[config-1]
print 'Running %s...'%ICs.name
(nodes, nx, ny, nz, xmin, xmax) = ICs.Fluid(n, dim)
if goalTime is None:
   goalTime = ICs.endTime
if samples is not None:
   dtSample = goalTime / samples
   dtRestart = dtSample / restartsPerSample
else:
   dtSample = goalTime

simName = 'magnetized-shock-tube-%i-%ix%ix%i'%(config, nx, ny, nz)
dataDir = '/p/lscratcha/jnjohnso/' + simName
visitDir = dataDir + "/visit"
restartDir = dataDir + "/restart"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
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
#from Spasmos import pause
#pause()
hydro = Hydro3d(WT, WTPi, q)
hydro.gradhCorrection = True
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
#mhd.divBCleaner = MHD.BDivergenceCleanerType.hyperbolicCleaner
#mhd.divBCleaner = MHD.BDivergenceCleanerType.BiotSavartProjCleaner
mhd.divBCleaner = MHD.BDivergenceCleanerType.GreensFnProjCleaner

# The transverse directions are periodic and the ends of the box are "rigid", 
# which is like "reflecting" except that vectors other than the position and 
# velocity are translated across the boundary without being reflected. 
e1 = Vector3d(1, 0, 0)
e2 = Vector3d(0, 1, 0)
e3 = Vector3d(0, 0, 1)
plane1 = Plane3d(Vector3d(xmax[0], 0, 0), -e1)
plane2 = Plane3d(Vector3d(xmin[0], 0, 0),  e1)
plane3 = Plane3d(Vector3d(0, xmax[1], 0), -e2)
plane4 = Plane3d(Vector3d(0, xmin[1], 0),  e2)
plane5 = Plane3d(Vector3d(0, 0, xmax[2]), -e3)
plane6 = Plane3d(Vector3d(0, 0, xmin[2]),  e3)
BC1 = PeriodicBoundary3d(plane3, plane4)
BC2 = PeriodicBoundary3d(plane5, plane6)
BC3 = RigidBoundary3d(plane1)
BC4 = RigidBoundary3d(plane2)
hydro.appendBoundary(BC1)
hydro.appendBoundary(BC2)
hydro.appendBoundary(BC3)
hydro.appendBoundary(BC4)
mhd.appendBoundary(BC1)
mhd.appendBoundary(BC2)
mhd.appendBoundary(BC3)
mhd.appendBoundary(BC4)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator3d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(mhd)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = True
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(mhd)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
restartBaseName = '%s/%s'%(restartDir, simName)
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            initializeMassDensity = True,
                            restartBaseName = restartBaseName)
output("control")

# Restore if desired.
if restoreCycle is not None:
   if restoreCycle == -1:
      restoreCycle = findLastRestart(simName)
   control.loadRestartFile(restoreCycle)
else:
   if samples is not None:
      dumpPhysicsState(integrator, simName, visitDir, dumpDerivatives = True)
output("nodes.numNodes")

#control.step(1)
#DBxDt = [u.dot(e1) for u in nodes.DBDt().internalValues()]
#ay = [u.dot(e2) for u in nodes.DvelocityDt().internalValues()]
#az = [u.dot(e3) for u in nodes.DvelocityDt().internalValues()]
#maxAy = max(ay)
#print 'Maximum transverse acceleration %g at %d'%(maxAy, ay.index(maxAy))
#dumpPhysicsState(integrator, simName, visitDir, dumpDerivatives = True)
#control.step(25)
#dumpPhysicsState(integrator, simName, visitDir, dumpDerivatives = True)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------

if samples is None:
   control.advance(goalTime, maxSteps)
else:
   nextSampleTime = control.time() + dtSample
   while control.time() < goalTime:
      nextGoalTime = min(control.time() + dtRestart, goalTime)
      control.advance(nextGoalTime, maxSteps)
      control.dropRestartFile()
      if abs(control.time() - nextSampleTime) < 1e-10 or \
         control.time() > nextSampleTime:
         dumpPhysicsState(integrator, simName, visitDir, dumpDerivatives = True)
         nextSampleTime = control.time() + dtSample

#-------------------------------------------------------------------------------
# Evaluate the solution.  Here we load the results of the corresponding 1D
# Athena solution for comparison and generate line outs that produce our 
# computed solution at the same points.
#-------------------------------------------------------------------------------
# Retrieve the answers from the right file.
import pickle
file = open('shock-tube-answers-%d-%g.dat'%(config, goalTime), 'r')
[axs, arhos, avxs, avys, avzs, aus, aBxs, aBys, aBzs] = pickle.load(file)
file.close()

# Toss the data that falls outside of [-0.5, 0.5].
bad = [i for i in xrange(len(axs)) if axs[i] < -0.5 or axs[i] > 0.5]
bad.reverse()
for ibad in bad:
   del axs[ibad]
   del arhos[ibad]
   del avxs[ibad]
   del avys[ibad]
   del aus[ibad]
   del aBxs[ibad]
   del aBys[ibad]
   del aBzs[ibad]

# Sample points along the axis of the flow.
if dim == 1:
   points = [Vector3d(axi, 0, 0) for axi in axs]
elif dim == 2:
   invSqrt2 = 1.0/sqrt(2)
   points = [Vector3d(invSqrt2*axi, invSqrt2*axi, 0) for axi in axs]
else:
   invSqrt22 = 1.0/(sqrt(2)*sqrt(2))
   points = [Vector3d(invSqrt22*axi, invSqrt22*axi, invSqrt22*axi) for axi in axs]
rhos = ScalarFieldList3d()
rhos.appendField(nodes.massDensity())
rhos = [rhos(xi, WT) for xi in points]
vs = VectorFieldList3d()
vs.appendField(nodes.velocity())
vSamples = [vs(xi, WT) for xi in points]
vxs = [val.x for val in vSamples]
vys = [val.y for val in vSamples]
vzs = [val.z for val in vSamples]
us = ScalarFieldList3d()
us.appendField(nodes.specificThermalEnergy())
us = [us(xi, WT) for xi in points]
Bs = VectorFieldList3d()
Bs.appendField(nodes.magneticInduction())
BSamples = [Bs(xi, WT) for xi in points]
Bxs = [val.x for val in BSamples]
Bys = [val.y for val in BSamples]
Bzs = [val.z for val in BSamples]
dBdts = VectorFieldList3d()
dBdts.appendField(nodes.DBDt())
dBdtSamples = [dBdts(xi, WT) for xi in points]
dBxdts = [val.x for val in dBdtSamples]
dBydts = [val.y for val in dBdtSamples]
dBzdts = [val.z for val in dBdtSamples]
Js = VectorFieldList3d()
Js.appendField(nodes.currentDensity())
JSamples = [Js(xi, WT) for xi in points]

# Paint some pretty pictures.
print 'Generating plots along y = z = 0:'
from pylab import *
ion()

def plotLine(quantity, x, Q, Qans, filename, deriv = None):
   clf()
   plot(x, Q, 'k.', x, Qans, 'k')
   if deriv is not None:
       plot(x, deriv, 'r.')
   title(quantity)
   xlabel('x')
   print '%s -> %s'%(quantity, filename)
   savefig(filename)
   
plotLine('Mass density', axs, rhos, arhos, 'rho-%d-%d.png'%(config, n))
plotLine('Vx', axs, vxs, avxs, 'vx-%d-%d.png'%(config, n))
plotLine('Vy', axs, vys, avys, 'vy-%d-%d.png'%(config, n))
plotLine('Vz', axs, vys, avys, 'vz-%d-%d.png'%(config, n))
plotLine('Specific thermal energy', axs, us, aus, 'u-%d-%d.png'%(config, n))
plotLine('Bx', axs, Bxs, aBxs, 'Bx-%d-%d.png'%(config, n))
plotLine('By', axs, Bys, aBys, 'By-%d-%d.png'%(config, n))
plotLine('Bz', axs, Bzs, aBzs, 'Bz-%d-%d.png'%(config, n))

