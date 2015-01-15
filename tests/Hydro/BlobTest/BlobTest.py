#-------------------------------------------------------------------------------
# The Blob Test
# Based on the test presented in 
# 1.	Agertz O, Moore B, Stadel J, Potter D, Miniati F, Read J, et al. 
#       Fundamental differences between SPH and grid methods. Monthly Notices of
#       the Royal Astronomical Society. 2007; 380(3):963-978.
#       doi:10.1111/j.1365-2966.2007.12183.x.
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *

import mpi
import DistributeNodes

title("3-D integrated hydro test --  Blob Test")

#-------------------------------------------------------------------------------
# Rejecter to help establishing initial conditions.
# This is used by GenerateNodeDistrubtion3d to carve out a spherical region.
#-------------------------------------------------------------------------------
class SphericalRejecterBlob:
    def __init__(self, origin, radius):
        self.origin = origin
        self.radius = radius
        return
    def __call__(self, x, y, z, m, H):
        n = len(x)
        assert (len(y) == n and len(z) == n and len(m) == n and len(H) == n)
        xnew, ynew, znew, mnew, Hnew = [], [], [], [], []
        R2 = self.radius**2
        for i in xrange(n):
            if ((x[i] - self.origin[0])**2 +
                (y[i] - self.origin[1])**2 +
                (z[i] - self.origin[2])**2 < R2):
                xnew.append(x[i])
                ynew.append(y[i])
                znew.append(z[i])
                mnew.append(m[i])
                Hnew.append(H[i])
        return xnew, ynew, znew, mnew, Hnew

#-------------------------------------------------------------------------------
# Rejecter to help establishing initial conditions.
# This is used by GenerateNodeDistrubtion3d to cut out a spherical region cavity
#-------------------------------------------------------------------------------
class SphericalRejecter:
    def __init__(self, origin, radius):
        self.origin = origin
        self.radius = radius
        return
    def __call__(self, x, y, z, m, H):
        n = len(x)
        assert (len(y) == n and len(z) == n and len(m) == n and len(H) == n)
        xnew, ynew, znew, mnew, Hnew = [], [], [], [], []
        R2 = self.radius**2
        for i in xrange(n):
            if ((x[i] - self.origin[0])**2 +
                (y[i] - self.origin[1])**2 +
                (z[i] - self.origin[2])**2 >= R2):
                xnew.append(x[i])
                ynew.append(y[i])
                znew.append(z[i])
                mnew.append(m[i])
                Hnew.append(H[i])
        return xnew, ynew, znew, mnew, Hnew



#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(

    #Set external state (Internal state is dependent so dont set it)
    gamma = 5.0/3.0,
    rhoext = 1.0,
    Pequi = 1.0,
    mach = 2.7,  #mach number 

    # Geometry Ambient medium box
    xb0 = 0.0,
    xb1 = 10.0,
    yb0 = 0.0,
    yb1 = 10.0,
    zb0 = 0.0,
    zb1 = 40.0,

    #Blob radius and central location
    br = 1.0,
    bx = 5.0,
    by = 5.0, 
    bz = 5.0,

    chi = 10.0,  # Ratio of rhoblob/rhoext

    goalTKH = 2.5,  # Goal time in units of t_KH

    # Resolution and node seeding.
    nx1 = 64,
    ny1 = 64,
    nz1 = 256,

    nPerh = 1.51,

    SVPH = False,
    CSPH = False,
    ASPH = False,
    SPH = True,   # This just chooses the H algorithm -- you can use this with CSPH for instance.
    filter = 0.0,  # For CSPH
    Qconstructor = MonaghanGingoldViscosity,
    #Qconstructor = TensorMonaghanGingoldViscosity,
    boolReduceViscosity = False,
    nh = 5.0,
    aMin = 0.1,
    aMax = 2.0,
    linearConsistent = False,
    fcentroidal = 0.0,
    fcellPressure = 0.0,
    Cl = 1.0, 
    Cq = 0.75,
    Qlimiter = False,
    balsaraCorrection = False,
    epsilon2 = 1e-2,
    hmin = 1e-5,
    hmax = 0.5,
    hminratio = 0.1,
    cfl = 0.25,
    XSPH = False,
    epsilonTensile = 0.0,
    nTensile = 8,

    IntegratorConstructor = CheapSynchronousRK2Integrator,
    steps = None,
    vizCycle = 20,
    vizTime = 0.1,
    dt = 0.0001,
    dtMin = 1.0e-5, 
    dtMax = 0.1,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    HUpdate = IdealH,
    domainIndependent = False,
    rigorousBoundaries = False,
    dtverbose = False,

    densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
    compatibleEnergy = True,
    gradhCorrection = False,

    clearDirectories = False,
    restoreCycle = None,
    restartStep = 200,
    dataDir = "dumps-blobtest-3d",
    )

# Decide on our hydro algorithm.
if SVPH:
    if ASPH:
        HydroConstructor = ASVPHFacetedHydro
    else:
        HydroConstructor = SVPHFacetedHydro
elif CSPH:
    if ASPH:
        HydroConstructor = ACSPHHydro
    else:
        HydroConstructor = CSPHHydro
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

# Build our directory paths.
densityUpdateLabel = {IntegrateDensity : "IntegrateDensity",
                      SumDensity : "SumDensity",
                      RigorousSumDensity : "RigorousSumDensity",
                      SumVoronoiCellDensity : "SumVoronoiCellDensity"}
baseDir = os.path.join(dataDir,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       densityUpdateLabel[densityUpdate],
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "XSPH=%s" % XSPH,
                       "nPerh=%3.1f" % nPerh,
                       "fcentroidal=%1.3f" % fcentroidal,
                       "fcellPressure = %1.3f" % fcellPressure,
                       "%ix%ix%i" % (nx1, ny1, nz1))
restartDir = os.path.join(baseDir, "restarts")
restartBaseName = os.path.join(restartDir, "blob-3d-%ix%ix%i" % (nx1, ny1, nz1))

vizDir = os.path.join(baseDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "blobtest-3d-%ix%ix%i" % (nx1, ny1, nz1)

# Figure out the goal time.
csext = sqrt(gamma*Pequi/rhoext)
vext = mach*csext
tCrush = 2.0*br*sqrt(chi)/vext
tKH = 1.6*tCrush
goalTime = goalTKH * tKH

print "Computed times (tCrush, tKH, goalTime) = (%g, %g, %g)" % (tCrush, tKH, goalTime)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(baseDir):
        shutil.rmtree(baseDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
mu = 1.0
eos1 = GammaLawGasMKS(gamma, mu)
eos2 = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = WT # TableKernel(HatKernel(1.0, 1.0), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
#Ambient Medium nodes
outerNodes = makeFluidNodeList("outer", eos1,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh)
#Blob nodes
innerNodes = makeFluidNodeList("inner", eos2,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh)
nodeSet = (outerNodes, innerNodes)
for nodes in nodeSet:
    output("nodes.name")
    output("    nodes.hmin")
    output("    nodes.hmax")
    output("    nodes.hminratio")
    output("    nodes.nodesPerSmoothingScale")
del nodes

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
#Blob density is defined to be 10 times the external density
    rhoblob=rhoext*chi
    generatorOuter = GenerateNodeDistribution3d(nx1, ny1, nz1, rhoext,
                                                distributionType = "lattice",
                                                xmin = (xb0, yb0, zb0),
                                                xmax = (xb1, yb1, zb1),
    						rejecter = SphericalRejecter(origin = (bx, by, bz),
                                                                             radius = br),
                                                nNodePerh = nPerh,
                                                SPH = SPH)
    #generatorOuter = GenerateNodeDistribution3d(nx1, ny1, nz1, rhoext,
    #                                            distributionType = "lattice",
    #                                            xmin = (bx-br, by-br, bz-br),
    #                                            xmax = (bx+br, by+br, bz+br),
    #       				        rmin = br,
    #                                            nNodePerh = nPerh,
    #                                            SPH = SPH)
    #generatorInner = GenerateNodeDistribution3d(nx1, ny1, nz1, rhoblob,
    #                                            distributionType = "lattice",
    #                                            xmin = (xb0, yb0, zb0),
    #                                            xmax = (xb1, yb1, zb1),
    #                                            rejecter = SphericalRejecterBlob(origin = (bx, by, bz),
    #                                                                         radius = br),
    #                                            nNodePerh = nPerh,
    #                                            SPH = SPH)

    # Figure out a mass matched resolution for the blob.
    mouter = (xb1 - xb0)*(yb1 - yb0)*(zb1 - zb0)*rhoext/(nx1*ny1*nz1)
    nxinner = max(2, int(((2*br)**3*rhoblob/mouter)**(1.0/3.0) + 0.5))
    generatorInner = GenerateNodeDistribution3d(nxinner, nxinner, nxinner, rhoblob,
                                                distributionType = "lattice",
                                                xmin = (bx-br, by-br, bz-br),
                                                xmax = (bx+br, by+br, bz+br),
 				                origin = (bx, by, bz),
           				        rmax = br,
                                                nNodePerh = nPerh,
                                                SPH = SPH)

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d

    distributeNodes3d((outerNodes, generatorOuter),
                      (innerNodes, generatorInner))
    for nodes in nodeSet:
        print nodes.name, ":"
        output("    mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.SUM)")
    del nodes

    # Set node specific thermal energies
    for (nodes, rho) in ((outerNodes, rhoext),
                         (innerNodes, rhoblob)):
        eps0 = Pequi/((gamma - 1.0)*rho)
        cs = sqrt(gamma*(gamma-1.0)*eps0)
        nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))
    del nodes

    #for nodes in nodeSet:
    #  vel = nodes.velocity()
    #  for i in xrange(nodes.numInternalNodes):
    #    vel[i]=Vector(velx,vely)

    vel = outerNodes.velocity() #wind velocity
    for i in xrange(outerNodes.numInternalNodes):
        vel[i].z = vext

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    db.appendNodeList(nodes)
del nodes
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = HydroConstructor(WT, q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             densityUpdate = densityUpdate,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             generateVoid = False,
                             HUpdate = HUpdate,
                             fcentroidal = fcentroidal,
                             fcellPressure = fcellPressure,
                             xmin = Vector(xb0 - (xb1 - xb0), yb0 - (yb1 - yb0), zb0 - (zb1 - zb0)),
                             xmax = Vector(xb1 + (xb1 - xb0), yb1 + (yb1 - yb0), zb1 + (zb1 - zb0)))
elif CSPH:
    hydro = HydroConstructor(WT, WTPi, q,
                             filter = filter,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(WT,
                             WTPi,
                             q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------

if boolReduceViscosity:
    #q.reducingViscosityCorrection = True
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(q,nh,aMin,aMax)
    
    packages.append(evolveReducingViscosityMultiplier)


#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(xb0, yb0, zb0), Vector( 1.0,  0.0, 0.0))
xPlane1 = Plane(Vector(xb1, yb0, zb0), Vector(-1.0,  0.0, 0.0))
yPlane0 = Plane(Vector(xb0, yb0, zb0), Vector( 0.0,  1.0, 0.0))
yPlane1 = Plane(Vector(xb0, yb1, zb0), Vector( 0.0, -1.0, 0.0))
zPlane0 = Plane(Vector(xb0, yb0, zb0), Vector( 0.0,  0.0, 1.0))
zPlane1 = Plane(Vector(xb0, yb0, zb1), Vector( 0.0,  0.0, -1.0))

xbc = PeriodicBoundary(xPlane0, xPlane1)
ybc = PeriodicBoundary(yPlane0, yPlane1)
zbc = PeriodicBoundary(zPlane0, zPlane1)

bcSet = [xbc, ybc, zbc]

# xbc0 = ReflectingBoundary(xPlane0)
# xbc1 = ReflectingBoundary(xPlane1)
# ybc0 = ReflectingBoundary(yPlane0)
# ybc1 = ReflectingBoundary(yPlane1)
# zbc0 = ReflectingBoundary(zPlane0)
# zbc1 = ReflectingBoundary(zPlane1)

# bcSet = [xbc0, xbc1, ybc0, ybc1, zbc0, zbc1]

for p in packages:
    for bc in bcSet:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.cullGhostNodes = False
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")
output("integrator.rigorousBoundaries")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            skipInitialPeriodicWork = (HydroConstructor in (SVPHFacetedHydro, ASVPHFacetedHydro)),
                            SPH = SPH)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()
