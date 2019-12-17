from Spheral2d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution2d import GenerateNodeDistribution2d
from PeanoHilbertDistributeNodes import distributeNodes2d
from FacetedVolumeRejecters import PolygonalSurfaceRejecter
import os, shutil
from math import *

title("Flow around an embedded triangle")

units = CGS()

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx = 400,
            ny = 300,

            # Hydro choices
            crksph = False,                  # True->CRK, False->SPH
            asph = False,                    # Choose the H algorithm, works with CRK or SPH
            correctionOrder = LinearOrder,   # for crksph
            volumeType = RKVoronoiVolume,   # for crksph
            densityUpdate = RigorousSumDensity,

            # Time advancement
            goalTime = 100.0,
            steps = None,

            # Output
            vizStep = None,
            vizTime = 0.1,
            clearDirectories = True,
            )


# Problem geometry
x0, y0 = 0.0, -1.5
x1, y1 = 4.0,  1.5

# Material initial conditions
rho0 = 1.0
v0 = 1.0

# Assorted hydro parameters
nPerh = 4.0

# Triangle geometry
L = 0.5  # triangle edge length
triverts = [Vector(1.0, 0.0),
            Vector(1.0 + L*cos(pi/6.0), -L*sin(pi/6.0)),
            Vector(1.0 + L*cos(pi/6.0),  L*sin(pi/6.0))]
triangle = Polygon(triverts)

# Set the directory paths
if ASPH:
    hydroname = "A"
else:
    hydroname = ""
if crksph:
    hydroname = "CRKSPH"
else:
    hydroname = "SPH"

dataDir = os.path.join("dumps-FlowAroundTriangle",
                       hydroname,
                       "nx=%i_ny=%i" % (nx, ny))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "viz")
vizName = "flowAroundTriangle"

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
# Material properties.
#-------------------------------------------------------------------------------
gamma = 5.0/3.0
mu = 1.0
eos = GammaLawGas(gamma, mu, units)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("fluid", eos, 
                          nPerh = nPerh,
                          kernelExtent = kernelExtent)
output("nodes.name")
output("    nodes.hmin")
output("    nodes.hmax")
output("    nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Initial conditions
#-------------------------------------------------------------------------------
gen = GenerateNodeDistribution2d(nx, ny, rho0, "lattice",
                                 xmin = (x0, y0),
                                 xmax = (x1, y1),
                                 nNodePerh = nPerh,
                                 SPH = not asph,
                                 rejecter = PolygonalSurfaceRejecter(triangle))
distributeNodes2d((nodes, gen))
output("nodes.numNodes")

# Set node velocities
vel = nodes.velocity()
for i in xrange(nodes.numInternalNodes):
    vel[i].x = v0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   correctionOrder = correctionOrder,
                   volumeType = volumeType,
                   densityUpdate = densityUpdate,
                   ASPH = asph)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                densityUpdate = densityUpdate,
                ASPH = asph)

hydro = SPH(dataBase = db,
            W = WT)
output("hydro")
output("hydro.kernel")
output("hydro.PiKernel")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.evolveTotalEnergy")
output("hydro.densityUpdate")
output("hydro.XSPH")

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
inplane =     Plane(Vector(x0, y0), Vector( 1,  0))
outplane =    Plane(Vector(x1, y0), Vector(-1,  0))
bottomplane = Plane(Vector(x0, y0), Vector( 0,  1))
topplane =    Plane(Vector(x0, y1), Vector( 0, -1))

inflow = InflowOutflowBoundary(db, inplane)
outflow = InflowOutflowBoundary(db, outplane)
bottom = ReflectingBoundary(bottomplane)
top = ReflectingBoundary(topplane)
obstacle = FacetedVolumeBoundary(triangle,
                                 interiorBoundary = True,
                                 useGhosts = True)
bcs = [inflow, outflow, bottom, top, obstacle]
packages += [inflow, outflow]

for p in packages:
    for bc in bcs:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
del p
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.updateBoundaryFrequency")
output("integrator.domainDecompositionIndependent")
output("integrator.cullGhostNodes")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
#import SpheralPointmeshSiloDump
control = SpheralController(integrator, WT,
                            vizBaseName = vizName,
                            redistributeStep = 50,
#                            vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState,
#                            vizGhosts = True,
                            vizDir = vizDir,
                            vizTime = vizTime,
                            vizStep = vizStep)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)
