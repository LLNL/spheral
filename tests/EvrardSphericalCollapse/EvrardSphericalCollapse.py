#-------------------------------------------------------------------------------
# Spheral++ script to model the gravitationally driven collapse of a
# sphere of gas.
# 
# See Evrard 1988
#-------------------------------------------------------------------------------
import os, sys, shutil, mpi
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution3d import *
from GenerateNodeDistribution3d import *
from VoronoiDistributeNodes import distributeNodes3d as distributeNodes

title("3-D gravitational hydro test -- Evrard spherical collapse.")

#-------------------------------------------------------------------------------
# Problem parameters.
#-------------------------------------------------------------------------------
commandLine(nx = 100,               # Number of across diameter of sphere

            # Gravity parameters.
            nbody = False,          # True->NBodyGravity, False->OctTreeGravity
            opening = 1.0,          # (dimensionless, OctTreeGravity) opening parameter for tree walk
            fdt = 0.1,              # (dimensionless, OctTreeGravity) timestep multiplier

            # Hydro parameters.
            crksph = False,
            psph = False,
            asph = False,  # Selects the H update algorithm -- can be used with CRK, PSPH, SPH, etc.
            kernelOrder = 5,
            correctionOrder = LinearOrder,      # CRK
            volumeType = RKSumVolume,          # CRK
            densityUpdate = RigorousSumDensity,
            HUpdate = IdealH,
            filter = 0.0,
            nPerh = 1.51,
            hmin = 1.0e-10,
            hmax = 10.0,
            hminratio = 0.1,
            Cl = None,
            Cq = None,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4.0,
            compatibleEnergy = True,
            evolveTotalEnergy = False,

            # Timestepping/advancement
            IntegratorConstructor = VerletIntegrator,
            goalTime = 0.8,
            steps = None,
            vizCycle = 20,
            vizTime = 0.05,
            dt = 1.0e-5,
            dtMin = 1.0e-8, 
            dtMax = 1.0e-1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,

            clearDirectories = False,
            restoreCycle = -1,
            restartStep = 50,
            checkRestart = False,
            dataRoot = "dumps-EvrardSphereCollapse-3d",
            vizName = "EvrardSphereCollapse",
            )

# Gravity parameters.  We set the Plummer softening scale based on the base particle resolution.
plummerLength = 2.0/nx
G = 1.0

# Set our paths.
if crksph:
    hydroname = os.path.join("CRKSPH",
                             "volume=%s" % volumeType)
elif psph:
    hydroname = "PSPH"
else:
    hydroname = "SPH"
if asph:
    hydroname = "A" + hydroname

if nbody:
    gravityname = "NBodyGravity"
else:
    gravityname = "TreeGravity"

dataDir = os.path.join(dataRoot,
                       hydroname,
                       gravityname,
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "nx=%i" % nx)
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "EvrardSphereCollapse")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(kernelOrder), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Material properties
#-------------------------------------------------------------------------------
gamma = 5.0/3.0
mu = 1.0
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("gas", eos,
                          hmin = hmin,
                          hmax = hmax,
                          hminratio = hminratio,
                          kernelExtent = WT.kernelExtent,
                          nPerh = nPerh)

#-------------------------------------------------------------------------------
# Create the generators.
#-------------------------------------------------------------------------------
M0 = 1.0
R0 = 1.0
rho0 = M0/(2.0*pi*R0*R0)
def rhofunc(posi):
    return rho0/max(posi.magnitude(), 1e-6)
gen = GenerateNodeDistribution3d(nx, nx, nx,
                                 rho = rhofunc,
                                 distributionType = "lattice",
                                 xmin = (-R0, -R0, -R0),
                                 xmax = ( R0,  R0,  R0),
                                 rmin = 0.0,
                                 rmax = R0,
                                 nNodePerh = nPerh,
                                 SPH = not asph)
distributeNodes((nodes, gen))

# Force the mass to be exactly the desired total.
mass = nodes.mass()
M1 = mass.sumElements()
for i in xrange(nodes.numInternalNodes):
    mass[i] *= M0/M1

print "Num internal nodes for ", nodes.name, " : ", mpi.allreduce(nodes.numInternalNodes, mpi.SUM)
print "Total mass: ", mass.sumElements()

# Set specific thermal energy.
eps0 = 0.05
eps = nodes.specificThermalEnergy()
for i in xrange(nodes.numInternalNodes):
    eps[i] = eps0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Gravity baby!
#-------------------------------------------------------------------------------
if nbody:
    gravity = NBodyGravity(G = G,
                           plummerSofteningLength = plummerLength,
                           maxDeltaVelocity = fdt)
else:
    gravity = OctTreeGravity(G = G,
                             softeningLength = plummerLength,
                             opening = opening,
                             ftimestep = fdt)

packages = [gravity]

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT, 
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   correctionOrder = correctionOrder,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 filter = filter,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 HopkinsConductivity = HopkinsConductivity,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 ASPH = asph)
else:
    hydro = SPH(dataBase = db,
                W = WT, 
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                densityUpdate = densityUpdate,
                XSPH = XSPH,
                HUpdate = HEvolution,
                ASPH = asph)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")

q = hydro.Q
if Cl:
    q.Cl = Cl
if Cq:
    q.Cq = Cq
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

packages += [hydro]

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(gravity)")
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
                            vizBaseName = vizName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if steps is None:
    control.advance(goalTime, maxSteps)
else:
    control.step(steps)

#-------------------------------------------------------------------------------
# Write the radial profiles.
#-------------------------------------------------------------------------------
pos = nodes.positions()
rho = nodes.massDensity()
vel = nodes.velocity()
P = ScalarField("pressure", nodes)
nodes.pressure(P)
stuff = zip(mpi.allreduce([x.magnitude() for x in pos.internalValues()], mpi.SUM),
            mpi.allreduce(list(rho.internalValues()), mpi.SUM),
            mpi.allreduce([vel[i].dot(pos[i].unitVector()) for i in xrange(nodes.numInternalNodes)], mpi.SUM),
            mpi.allreduce(list(P.internalValues()), mpi.SUM))
stuff.sort()
if mpi.rank == 0:
    f = open(os.path.join(dataDir, "EvrardCollapseRadialProfiles.gnu"), "w")
    f.write(("# " + 5*"%20s " + "\n") % ("r", "rho", "vr", "P", "Aentropy"))
    for ri, rhoi, vri, Pi in stuff:
        Ai = Pi/rhoi**gamma
        f.write((5*"%20g " + "\n") % (ri, rhoi, vri, Pi, Ai))
    f.close()
