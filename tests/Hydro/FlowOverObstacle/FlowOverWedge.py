#-------------------------------------------------------------------------------
# Based roughly on the discussion and results in section 17.2.1 of Toro.
#-------------------------------------------------------------------------------
from Spheral2d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution2d import GenerateNodeDistribution2d
from CompositeNodeDistribution import *
from PeanoHilbertDistributeNodes import distributeNodes2d
from FacetedVolumeRejecters import PolygonalSurfaceRejecter
import os, shutil
from math import *

title("Shock over a wedge boundary")

units = CGS()

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx = 250,
            ny = 150,

            angle = 25.0,                    # Wedge angle (degrees)
            v0 = 1.0,                        # inflow velocity

            xmem = 3.0,                      # initial position between high and low pressure regions

            # Hydro choices
            crksph = False,                  # True->CRK, False->SPH
            asph = False,                    # Choose the H algorithm, works with CRK or SPH
            correctionOrder = LinearOrder,   # for crksph
            volumeType = RKVoronoiVolume,   # for crksph
            densityUpdate = RigorousSumDensity,

            # Time advancement
            goalTime = 500.0,
            steps = None,

            # Output
            vizStep = None,
            vizTime = 5.0,
            clearDirectories = True,
            )

# Problem geometry
x0, x1, x2 = 0.0, 4.0, 25.0
y0, y1 = 0.0, 15.0

# Material initial conditions
rho0, eps0 = 1.0,  1.0

# Assorted hydro parameters
nPerh = 4.0

# Set the directory paths
if ASPH:
    hydroname = "A"
else:
    hydroname = ""
if crksph:
    hydroname = "CRKSPH"
else:
    hydroname = "SPH"

dataDir = os.path.join("dumps-FlowOverWedge",
                       "wedgeAngle=%s" % angle,
                       "v0=%s" % v0,
                       hydroname,
                       "nx=%i_ny=%i" % (nx, ny))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "viz")
vizName = "flowOverWedge"

# Build the wedge geometry
angle *= pi/180.0
wedgePts = [(x1, y0),
            (x2, y0),
            (x2, y0 + (x2 - x1)*sin(angle))]
wedge = Polygon([Vector(x,y) for x,y in wedgePts])

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
gen = GenerateNodeDistribution2d(nx, ny, rho0,
                                 "lattice",
                                 xmin = (x0, y0),
                                 xmax = (x2, y1),
                                 nNodePerh = nPerh,
                                 SPH = not asph,
                                 rejecter = PolygonalSurfaceRejecter(wedge))
distributeNodes2d((nodes, gen))
output("nodes.numNodes")

# Set node initial conditions
vel = nodes.velocity()
eps = nodes.specificThermalEnergy()
for i in xrange(nodes.numInternalNodes):
    vel[i].x = v0
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
outplane =    Plane(Vector(x2, y0), Vector(-1,  0))
bottomplane = Plane(Vector(x0, y0), Vector( 0,  1))
topplane =    Plane(Vector(x0, y1), Vector( 0, -1))

inflow = InflowOutflowBoundary(db, inplane)
outflow = InflowOutflowBoundary(db, outplane)
bottom = ReflectingBoundary(bottomplane)
top = ReflectingBoundary(topplane)
wedgeBC = FacetedVolumeBoundary(wedge,
                                interiorBoundary = True,
                                useGhosts = True)
bcs = [inflow, outflow, top, bottom, wedgeBC]
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
