from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralVisitDump import SpheralVisitDump
from findLastRestart import *
from math import *

# Load the mpi module if we're parallel.
try:
    import mpi
except:
    mpi = None

from GenerateNodeDistribution2d import *

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
seed = 'lattice'
NodeListConstructor = AsphNodeList2d

nx1, ny1 = 40, 50
nx2, ny2 = 40, 50
nPerh = 2.01

rho1, rho2 = 1.0, 1.0 #1.0, 0.25
P1, P2 = 1.0, 1.0 #1.0, 0.25

vx0, waveLength = 1.0, 1.0
vy0 = 0.0
dvy = 0.0 #0.005
perturbLength = 0.1

xxmax = nx1/float(ny1 + ny2) * 1.0
xmin1, xmax1 = Vector2d(0.0, 0.0), Vector2d(xxmax, 0.5)
xmin2, xmax2 = Vector2d(0.0, 0.5), Vector2d(xxmax, 1.0)

gamma = 5.0/3.0
mu = 1.0
#Qconstructor = MonaghanGingoldViscosity2d
Qconstructor = TensorMonaghanGingoldViscosity2d
Cl, Cq = 1.0, 0.75
Qlimiter = 1
balsaraCorrection = 0
epsilon2 = 1e-4
negligibleSoundSpeed = 1.0e-5
csMultiplier = 0.0001
HsmoothMin, HsmoothMax = 0.0001, 0.2
HsmoothFraction = 0.1
cfl = 0.5

neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
numGridLevels = 25
topGridCellSize = 2.0
origin = Vector2d(0.0, 0.0)

goalTime = 1.0
dtdump = 0.1
dt = 0.01
dtMin, dtMax = 1.0e-5, 0.1
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Integrator2d.HEvolutionType.IdealH # IntegrateH
sumForMassDensity = Integrator2d.MassDensityType.RigorousSum

restartStep = 1000
restartBaseName = "SinShear-2d-%ix%i" % (nx1, ny1 + ny2)
restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
title('2-D planar shearing test.')

eos = GammaLawGasMKS2d(gamma, mu)

# Create a pair of empty NodeLists
nodes1 = NodeListConstructor(eos, 0)
nodes2 = NodeListConstructor(eos, 0)
output('nodes1')
output('nodes2')
nodes1.HsmoothFraction = HsmoothFraction
nodes2.HsmoothFraction = HsmoothFraction

# Create our interpolation kernels -- one for normal hydro interactions, and
WT = TableKernel2d(BSplineKernel2d(), 100)
WTPi = WT
kernelExtent = WT.kernelExtent()
output('WT')
output('WTPi')
output('kernelExtent')

# Construct the neighbor object and associate it with the node list.
neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
neighbor2 = NestedGridNeighbor2d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes2.registerNeighbor(neighbor2)

# Set node positions, masses, and H's for this domain.
if restoreCycle is None:
    from DistributeNodes import distributeNodes2d
    generator1 = GenerateNodeDistribution2d(nx1, ny1, rho1, seed,
                                            xmin = (xmin1.x, xmin1.y),
                                            xmax = (xmax1.x, xmax1.y),
                                            nNodePerh = nPerh)
    generator2 = GenerateNodeDistribution2d(nx2, ny2, rho2, seed,
                                            xmin = (xmin2.x, xmin2.y),
                                            xmax = (xmax2.x, xmax2.y),
                                            nNodePerh = nPerh)
    n1 = generator1.globalNumNodes()
    n2 = generator2.globalNumNodes()
    nodeTuples = [(nodes1, n1, generator1),
                  (nodes2, n2, generator2)]
    nodeInfo = distributeNodes2d(nodeTuples)

    if mpi:
        output('mpi.reduce(nodes1.numInternalNodes, mpi.MIN)')
        output('mpi.reduce(nodes1.numInternalNodes, mpi.MAX)')
        output('mpi.reduce(nodes1.numInternalNodes, mpi.SUM)')

        output('mpi.reduce(nodes2.numInternalNodes, mpi.MIN)')
        output('mpi.reduce(nodes2.numInternalNodes, mpi.MAX)')
        output('mpi.reduce(nodes2.numInternalNodes, mpi.SUM)')
    else:
        output('nodes1.numInternalNodes')
        output('nodes2.numInternalNodes')

    assert len(nodeInfo[nodes1]['globalNodeListID']) == nodes1.numInternalNodes
    assert len(nodeInfo[nodes2]['globalNodeListID']) == nodes2.numInternalNodes

    # Set node specific thermal energies
    eps1 = P1/((gamma - 1.0)*rho1)
    eps2 = P2/((gamma - 1.0)*rho2)
    nodes1.specificThermalEnergy()[:] = [eps1]*nodes1.numNodes
    nodes2.specificThermalEnergy()[:] = [eps2]*nodes2.numNodes

    # Set node velocities
    for nodes in [nodes1, nodes2]:
        for nodeID in xrange(nodes.numNodes):
            dy = (nodes.positions()[nodeID].y - xmin1.y)/(xmax2.y - xmin1.y)
            nodes.velocity()[nodeID] = Vector2d(vx0*sin(2.0*pi*waveLength*dy), vy0)

    # If requested, perturb the velocity across the interface.
    xlength = xmax1.x - xmin1.x
    ylength = xmax2.y - xmin1.y
    assert xlength > 0.0
    assert ylength > 0.0
    for nodeID in xrange(nodes1.numNodes):
        r = nodes1.positions()[nodeID]
        dvyi = dvy*(sin(2.0*pi*(r.x - xmin1.x)/xlength)*
                    exp(-((r.y - xmax1.y)/perturbLength)**2))
        nodes1.velocity()[nodeID].y += dvyi
    for nodeID in xrange(nodes2.numNodes):
        r = nodes2.positions()[nodeID]
        dvyi = dvy*(sin(2.0*pi*(r.x - xmin2.x)/xlength)*
                    exp(-((r.y - xmax1.y)/perturbLength)**2))
        nodes2.velocity()[nodeID].y += dvyi

    # If we're integrating for mass density, 'ya have to initialize the density.
    nodes1.massDensity()[:] = [rho1]*nodes1.numNodes
    nodes2.massDensity()[:] = [rho2]*nodes2.numNodes

# Create boundary conditions.  We will enforce reflecting boundaries at
# y=0 and y=1, and a periodic boundary wrapping around x=0 and x=1.
xPlane0 = Plane2d(xmin1, Vector2d(1.0, 0.0))
xPlane1 = Plane2d(xmax2, Vector2d(-1.0, 0.0))
yPlane0 = Plane2d(xmin1, Vector2d(0.0, 1.0))
yPlane1 = Plane2d(xmax2, Vector2d(0.0, -1.0))
xbc0 = PeriodicBoundary2d(xPlane0, xPlane1)
ybc0 = PeriodicBoundary2d(yPlane0, yPlane1)
##ybc0 = ReflectingBoundary2d(yPlane0)
##ybc1 = ReflectingBoundary2d(yPlane1)

# Construct a DataBase to hold our node list
db = DataBase2d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Construct a standard Monaghan-Gingold artificial viscosity.
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.epsilon2 = epsilon2
q.balsaraShearCorrection = balsaraCorrection
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output('q')
output('q.Cl')
output('q.Cq')
output('q.epsilon2')
output('q.limiter')
output('q.negligibleSoundSpeed')
output('q.csMultiplier')
output('q.balsaraShearCorrection')

# Construct the hydro physics object.
hydro = Hydro2d(WT, WTPi, q)
hydro.cfl = cfl
output('hydro')
output('hydro.kernel()')
output('hydro.PiKernel()')
output('hydro.cfl')
output('hydro.valid()')

# Construct a synchronous RK2 integrator, and add the physics packages and
# boundary condtions.
integrator = PredictorCorrectorIntegrator2d(db)
output('integrator')
integrator.appendPhysicsPackage(hydro)
output('integrator.havePhysicsPackage(hydro)')
output('integrator.valid()')
integrator.rigorousBoundaries = 1
integrator.HsmoothMin = HsmoothMin
integrator.HsmoothMax = HsmoothMax
output('integrator.HsmoothMin')
output('integrator.HsmoothMax')
integrator.lastDt = dt
output('integrator.lastDt')
if dtMin:
    integrator.dtMin = dtMin
    output('integrator.dtMin')
if dtMax:
    integrator.dtMax = dtMax
    output('integrator.dtMax')
integrator.dtGrowth = dtGrowth
output('integrator.dtGrowth')
integrator.HEvolution = HEvolution
output('integrator.HEvolution')
integrator.sumForMassDensity = sumForMassDensity
if (sumForMassDensity == Integrator2d.MassDensityType.Sum or
    sumForMassDensity == Integrator2d.MassDensityType.RigorousSum):
    integrator.setKernel(WT)
output('integrator.sumForMassDensity')

control = SpheralController(integrator, WT,
                            boundaryConditions = [xbc0, ybc0],
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output('control')

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)
    P = db.fluidPressure
    cs = db.fluidSoundSpeed
    Hi = db.fluidHinverse
    dumper = SpheralVisitDump(db,
                              "SinShear-2d-%ix%i-visit" % (nx1, ny1 + ny2),
                              "visit-%ix%i" % (nx1, ny1 + ny2),
                              listOfFieldLists = [db.fluidMassDensity,
                                              db.fluidVelocity,
                                                  db.fluidWeight,
                                                  P,
                                                  cs,
                                                  Hi]
                              )
    dumper.dump(control.time(), control.totalSteps)

##################################################################################
# Advance to the end time.
nextGoalTime = min(goalTime, control.time() + dtdump)
while control.time() < goalTime:
    while control.time() < nextGoalTime:
        control.advance(nextGoalTime)
        P = db.fluidPressure
        cs = db.fluidSoundSpeed
        Hi = db.fluidHinverse
        dumper = SpheralVisitDump(db,
                                  "SinShear-2d-%ix%i-visit" % (nx1, ny1 + ny2),
                                  "visit-%ix%i" % (nx1, ny1 + ny2),
                                  listOfFieldLists = [db.fluidMassDensity,
                                                      db.fluidVelocity,
                                                      db.fluidWeight,
                                                      P,
                                                      cs,
                                                      Hi]
                                  )
        dumper.dump(control.time(), control.totalSteps)
    nextGoalTime = min(goalTime, control.time() + dtdump)

control.dropRestartFile()

# Plot the postion and velocity fields.
rPlot = plotNodePositions2d(db)
#velFieldPlot = plotVelocityField2d(db, velMultiplier=0.05, title='Velocity Field', plotGhosts=0)
vxPlot = plotFieldList(db.fluidVelocity, xFunction='%s.y', yFunction='%s.x', plotStyle='points', winTitle='x velocity')
#vyPlot = plotFieldList(db.fluidVelocity, xFunction='%s.y', yFunction='%s.y', plotStyle='points', winTitle='y velocity')

# Plot the elongation (h1/h2) for the H tensors.
import Gnuplot
hratio0 = [H.eigenValues().minElement()/H.eigenValues().maxElement() for H in
           nodes1.Hfield().internalValues() + nodes2.Hfield().internalValues()]
y0 = [ri.y for ri in nodes1.positions().internalValues() + nodes1.positions().internalValues()]
hratio = mpi.allreduce(hratio0, mpi.SUM)
y = mpi.allreduce(y0, mpi.SUM)
hratioPlot = Gnuplot.Gnuplot()
cache = []
if mpi.rank == 0:
    hratioData = Gnuplot.Data(y, hratio,
                              title = "hmin/hmax ratio",
                              inline = True)
    hratioPlot.plot(hratioData)
    cache.extend([hratioData, hratioPlot])

### In the shearing shock tube case, it's interesting to be able to plot the Sod
### solution.
##rhoPlot = plotFieldList(db.fluidMassDensity, xFunction='%s.y', plotStyle='points', winTitle='Mass density')
##velPlot = plotFieldList(db.fluidVelocity, xFunction='%s.y', yFunction='%s.y', plotStyle='points', winTitle='Velocity (y)')
##epsPlot = plotFieldList(db.fluidSpecificThermalEnergy, xFunction='%s.y', plotStyle='points', winTitle='Specific Thermal Energy')
##P = db.fluidPressure
##PPlot = plotFieldList(P, xFunction='%s.y', plotStyle='points', winTitle='Pressure')
##HPlot = plotFieldList(db.fluidHfield, xFunction='%s.y', yFunction='1.0/%s.yy', plotStyle='points', winTitle='Y Smoothing Scale')

##import sys
##sys.path.append('../Sod')
##from SodAnalyticSolution import *
##answer = SodSolution(nPoints=ny1 + ny2,
##                     gamma = gamma,
##                     rho1 = rho1,
##                     P1 = P1,
##                     rho2 = rho2,
##                     P2 = P2,
##                     x0 = rmin1[1],
##                     x1 = rmax1[1],
##                     x2 = rmax2[1],
##                     h1 = nPerh*(rmax1[1] - rmin1[1]/ny1),
##                     h2 = nPerh*(rmax2[1] - rmin2[1]/ny2))
##plotAnswer(answer, control.time(),
##           rhoPlot, velPlot, epsPlot, PPlot, HPlot)

##import SpheralMayaViUtilities
##vtk = SpheralMayaViUtilities.VtkFile([(db.fluidMassDensity, "Mass Density"),
##                                      (db.fluidSpecificThermalEnergy, "Thermal energy"),
##                                      (db.fluidVelocity, "Velocity"),
##                                      (db.fluidHfield, "Hfield")],
##                                     db, topGridCellSize,
##                                     nx = 50,
##                                     ny = 50,
##                                     rmin = Vector2d(0.0, 0.0),
##                                     rmax = Vector2d(1.0, 1.0),
##                                     boundaryList = [xbc0, ybc0, ybc1],
##                                     header = 'Kelvin Helmholtz at t=1.')
##vtk.writeVTKFieldLists('PlanarShear-2d-30x30.vtk', 'ascii')
