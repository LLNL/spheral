from Numeric import *
from Spheral import *
from Strength import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralController import *

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
NodeListConstructor = SphSolidNodeList1d

nx1 = 100
rho1 = 1.0
eps1 = 0.0
x0, x1 = 0.0, 1.0
nPerh = 2.01

vr0, vrSlope = -1.0, 0.0
#vr0, vrSlope = 0.0, -1.0
#vr0, vrSlope = 0.0, 0.0

m1 = rho1*(x1 - x0)/nx1

gamma = 5.0/3.0
mu = 1.0
Qconstructor = MonaghanGingoldViscosity1d
#Qconstructor = TensorMonaghanGingoldViscosity1d
Cl, Cq = 1.0, 0.75
Qlimiter = 0
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
energyMultiplier = 0.1
HsmoothMin, HsmoothMax = 0.0001, 0.1
HsmoothFraction = 0.0
cfl = 0.5
XSPH = False
epsilonTensile = 0.0
nTensile = 8
epsilonTensileGradient = 0.01

neighborSearchType = Neighbor1d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = 2.0
origin = Vector1d(0.0)

goalTime = 0.6
dt = 0.0001
dtMin, dtMax = 1.0e-5, None #0.09
dtGrowth = 2.0
maxSteps = None
statsStep = 1
smoothIters = 0
HEvolution = Hydro1d.HEvolutionType.IdealH
sumForMassDensity = Hydro1d.MassDensityType.RigorousSumDensity

restoreCycle = None
restartStep = 10000
restartBaseName = "Noh-planar-1d-%i-%i" % (nx1)

#-------------------------------------------------------------------------------
title('1-D integrated hydro test -- planar Noh problem')

eos = GammaLawGasMKS1d(gamma, mu)

# Create an empty NodeList
nodesStrength = NullStrength()
nodes1 = NodeListConstructor("nodes1", eos, nodesStrength)
nodes1.HsmoothFraction = HsmoothFraction
nodes1.nodesPerSmoothingScale = nPerh
output('nodes1.HsmoothFraction')
output('nodes1.nodesPerSmoothingScale')

# Set node positions for this domain
from DistributeNodes import distributeNodes1d
list = [(nodes1, nx1, (x0, x1))]
distributeNodes1d(list)
nNodesThisDomain1 = nodes1.numInternalNodes
output('nodes1.numNodes')

# Set node masses
nodes1.setMass(ScalarField1d(nodes1.mass().name(), nodes1, m1))

# Set node densities.
nodes1.setMassDensity(ScalarField1d(nodes1.massDensity().name(), nodes1, rho1))
nodes1.updateWeight()

# Set node specific thermal energies
nodes1.setSpecificThermalEnergy(ScalarField1d(nodes1.specificThermalEnergy().name(), nodes1, eps1))

# Set node velocities
for ix in xrange(nodes1.numNodes):
    nodes1.velocity()[ix].x = vr0 + vrSlope*nodes1.positions()[ix].x

# Set the smoothing scales.
dx = (x1 - x0)/nx1
h1 = 1.0/(nPerh*dx)
H1 = SymTensor1d(h1)
nodes1.setHfield(SymTensorField1d(nodes1.Hfield().name(), nodes1, H1))

# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
WT = TableKernel1d(BSplineKernel1d(), 100)
kernelExtent = WT.kernelExtent()
WTPi = TableKernel1d(BSplineKernel1d(), 100)
#WTPi = TableKernel1d(HatKernel1d(kernelExtent, kernelExtent), 100)
output('WT')
output('WTPi')

# Construct the neighbor object and associate it with the node list.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
neighborTimer.stop()
neighborTimer.printStatus()

# Create boundary conditions.  We need at least this much to create the initial
# mass density field.
xPlane0 = Plane1d(Vector1d(x0), Vector1d(1.0))
xbc0 = ReflectingBoundary1d(xPlane0)

# Construct a DataBase to hold our node list
db = DataBase1d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Construct the artificial viscosities for the problem.
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
q.energyMultiplier = energyMultiplier
output('q')
output('q.Cl')
output('q.Cq')
output('q.limiter')
output('q.epsilon2')
output('q.negligibleSoundSpeed')
output('q.csMultiplier')
output('q.energyMultiplier')

# Set the XSPH and tensile corrections for the NodeList
nodes1.XSPH = XSPH
nodes1.epsilonTensile = epsilonTensile
nodes1.nTensile = nTensile
nodes1.epsilonTensileGradient = epsilonTensileGradient
output('nodes1.XSPH')
output('nodes1.epsilonTensile')
output('nodes1.nTensile')
output('nodes1.epsilonTensileGradient')

# Construct the hydro physics object.
hydro = Hydro1d(WT, WTPi, q)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = HsmoothMin
hydro.HsmoothMax = HsmoothMax
output('hydro')
output('hydro.cfl')
output('hydro.HEvolution')
output('hydro.sumForMassDensity')
output('hydro.HsmoothMin')
output('hydro.HsmoothMax')
output('hydro.kernel()')
output('hydro.PiKernel()')
output('hydro.valid()')

# Construct a strength physics object.
strength = Strength1d()
output('strength')

# Construct a predictor corrector integrator, and add the one physics package.
#integrator = CheapSynchronousRK2Integrator1d(db)
integrator = PredictorCorrectorIntegrator1d(db)
output('integrator')
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
output('integrator.havePhysicsPackage(hydro)')
output('integrator.havePhysicsPackage(strength)')
output('integrator.valid()')
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

control = SpheralController(integrator, WT,
                            boundaryConditions = [xbc0],
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output('control')

# Smooth the initial conditions.
if restoreCycle:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)

##################################################################################
# Advance to the end time.
control.step(5)
control.advance(goalTime, maxSteps)

# Plot the final state.
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotStyle='lines')

# Overplot the analytic solution.
import mpi
import sys
sys.path.append("../../Hydro/Noh")
import NohAnalyticSolution
rlocal = [pos.x for pos in nodes1.positions().internalValues()]
r = mpi.reduce(rlocal, mpi.SUM)
answer = NohAnalyticSolution.NohSolution(1,
                                         r = r,
                                         v0 = -1.0,
                                         h0 = 1.0/h1)
plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot)

def filterOnRadius(data, r, rmin, rmax):
    assert len(data) == len(r)
    return [data[i] for i in xrange(len(data)) if (r[i].magnitude() >= rmin and r[i].magnitude() <= rmax)]

# Compute the error.
rmin, rmax = 0.05, 0.35   # Throw away anything with r < rwall to avoid wall heating.
rhoprof = mpi.reduce(filterOnRadius(nodes1.massDensity().internalValues(), nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
Pprof = mpi.reduce(filterOnRadius(nodes1.pressure().internalValues(), nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
vprof = mpi.reduce(filterOnRadius([v.x for v in nodes1.velocity().internalValues()], nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
epsprof = mpi.reduce(filterOnRadius(nodes1.specificThermalEnergy().internalValues(), nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
hprof = mpi.reduce(filterOnRadius([1.0/H.xx for H in nodes1.Hfield().internalValues()], nodes1.positions().internalValues(), rmin, rmax), mpi.SUM)
if mpi.rank == 0:
    answer.r = [r for r in answer.r if r >= rmin and r <= rmax]
    rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time())
    import Pnorm
    print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
    for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                              ("Pressure", Pprof, Pans),
                              ("Velocity", vprof, vans),
                              ("Thermal E", epsprof, epsans),
                              ("h       ", hprof, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in xrange(len(data))]
        Pn = Pnorm.Pnorm(error, hprof)
        print "\t%s \t\t%g \t\t%g \t\t%g" % (name,
                                             Pn.gridpnorm(1),
                                             Pn.gridpnorm(2),
                                             Pn.gridpnorm("inf"))
