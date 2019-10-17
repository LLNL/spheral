from Spheral import *
from SpheralTestUtilities import *

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(geometry = "1d",    # (1d, 2d, 3d, RZ)
            
            nx1 = 100,
            nx2 = 100,

            rho1 = 1.0,
            rho2 = 1.0,
            vx0 = 1.0,

            x0 = -1.0,
            x1 = 0.0,
            x2 = 1.0,
            y0 = 0.0,
            y1 = 1.0,
            z0 = 0.0,
            z1 = 1.0,

            nPerh = 4.0,
            hmin = 1e-10,
            hmax = 1.0,

            dtMin = 1e-5,
            dtMax = 100.0,
            steps = None,
            goalTime = 1.0,

            graphics = True,
            )

assert geometry in ("1d", "2d", "3d", "RZ")
exec("from Spheral%s import *" % geometry)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
units = CGS()
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
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           kernelExtent = kernelExtent)
nodes2 = makeFluidNodeList("nodes2", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           kernelExtent = kernelExtent)
nodeSet = [nodes1, nodes2]
for nodes in nodeSet:
    output("nodes.name")
    output("    nodes.hmin")
    output("    nodes.hmax")
    output("    nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Initial conditions
#-------------------------------------------------------------------------------
if geometry == "1d":
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1)),
                              (nodes2, nx2, rho2, (x1, x2))],
                             nPerh = nPerh)
output("nodes1.numNodes")

# Set node velocities
for nodes in nodeSet:
    vel = nodes.velocity()
    for i in xrange(nodes.numInternalNodes):
        vel[i].x = vx0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = SPH(dataBase = db,
            W = WT)
output("hydro")
output("hydro.kernel")
output("hydro.PiKernel")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.XSPH")

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
if geometry == "1d":
    inplane = Plane(Vector(x0), Vector(1.0))
    inflow = InflowBoundary(nodes1, inplane)

    bcs = [inflow]

for p in packages:
    for bc in bcs:
        p.appendBoundary(bc)

packages += bcs

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
del p
integrator.lastDt = dtMin
integrator.dtMin = dtMin
integrator.cullGhostNodes = False
integrator.verbose = True
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
control = SpheralController(integrator, WT)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)

#-------------------------------------------------------------------------------
# Plotting
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
