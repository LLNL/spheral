#ATS:test(SELF, "--graphics False", label="Voronoi hourglass filter tests -- 1-D (serial)")
#-------------------------------------------------------------------------------
# Test out the Voronoi hourglass filtering algorithm.
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import *
from generateMesh import generateLineMesh

import numpy
import Gnuplot

#-------------------------------------------------------------------------------
# Command line parameters.
#-------------------------------------------------------------------------------
commandLine(nx = 100,
            rho0 = 1.0,
            rhoSlope = 1.0,

            x0 = 0.0,
            x1 = 1.0,

            nPerh = 2.01,

            gammaGas = 5.0/3.0,
            mu = 1.0,
            hmin = 1e-10,
            hmax = 1.0,
            hourglassOrder = 1,
            hourglassLimiter = 1,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            steps = None,
            goalTime = 0.15,
            dt = 1e-4,
            dtMin = 1.0e-5,
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 10,
            HEvolution = IdealH,

            iterations = 10,
            graphics = True,
            )

#-------------------------------------------------------------------------------
# Material.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gammaGas, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes1", eos, 
                          hmin = hmin,
                          hmax = hmax,
                          nPerh = nPerh)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes, [(nx, rho0, (x0, x1))])])
output("nodes.numNodes")

def setRho():
    pos = nodes.positions()
    rho = nodes.massDensity()
    for i in xrange(nodes.numInternalNodes):
        rho[i] = rho0 + rhoSlope*pos[i].x
    return

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct an hourglass control object.
#-------------------------------------------------------------------------------
hg = VoronoiHourglassControl(WT, hourglassOrder, hourglassLimiter)
output("hg")
output("hg.order")
output("hg.limiter")
packages = vector_of_Physics()
packages.append(hg)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0), Vector( 1.0))
xPlane1 = Plane(Vector(x1), Vector(-1.0))
xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)
hg.appendBoundary(xbc0)
hg.appendBoundary(xbc1)

#-------------------------------------------------------------------------------
# Iteratively let the hourglass control adjust the point positions, and lets
# see where we converge to.
#-------------------------------------------------------------------------------
mass = nodes.mass()
pos = nodes.positions()
H = nodes.Hfield()
rho = nodes.massDensity()

def plotRho(p):
    xarray = numpy.array([x.x for x in pos.internalValues()])
    rhoarray = numpy.array([x for x in rho.internalValues()])
    d = Gnuplot.Data(xarray, rhoarray, with_="linesp", title="iteration %i" % iter, inline=True)
    p.replot(d)

def plotDx(p):
    mesh, void = generateLineMesh([nodes], Vector(x0), Vector(x1), False, False, False)
    xarray = numpy.array([x.x for x in pos.internalValues()])
    n = len(xarray)
    dxarray = numpy.array([float(i) for i in xrange(n)])
    for i in xrange(n):
        dxarray[i] = mesh.zone(i).volume()
    d = Gnuplot.Data(xarray, dxarray, with_="linesp", title="iteration %i" % iter, inline=True)
    p.replot(d)

if graphics:
    p0, p1, p2 = Gnuplot.Gnuplot(), Gnuplot.Gnuplot(), Gnuplot.Gnuplot()
    p0.title("Forced density profile")
    p1.title("Summed density profile")
    p2.title("Delta x")
    
for iter in xrange(iterations):
    setRho()
    if graphics:
        plotRho(p0)
    state = State()
    for f in (pos, mass, rho, H):
        state.enroll(f)
    derivs = StateDerivatives(db, packages)
    hg.registerState(db, state)
    hg.registerDerivatives(db, derivs)
    state.update(derivs, 1.0, 0.0, 1.0)
    hg.finalize(0.0, 0.0, db, state, derivs)
    if graphics:
        db.updateConnectivityMap()
        cm = db.connectivityMap()
        posfl = state.vectorFields(HydroFieldNames.position)
        massfl = state.scalarFields(HydroFieldNames.mass)
        Hfl = state.symTensorFields(HydroFieldNames.H)
        rhofl = state.scalarFields(HydroFieldNames.massDensity)
        computeSPHSumMassDensity(cm, WT, posfl, massfl, Hfl, rhofl)
        plotRho(p1)
        plotDx(p2)
