#ATS:test(SELF, "--linearConsistent True --graphics False", label="SVPH interpolation test -- 1-D (serial)")
#-------------------------------------------------------------------------------
# A set of tests to compare how different meshless methods interpolate fields.
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import *
from generateMesh import *

title("SVPH Interpolation tests")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
            # Parameters for seeding nodes.
            nx1 = 100,
            rho1 = 1.0,
            eps1 = 0.0,
            x0 = 0.0,
            x1 = 1.0,
            nPerh = 1.25,
            hmin = 0.0001, 
            hmax = 10.0,

            # Should we randomly perturb the positions?
            ranfrac = 0.2,
            seed = 14892042,

            # What test problem are we doing?
            testCase = "linear",

            # Should we use linearly consistent formalism?
            linearConsistent = True,
            
            # The fields we're going to interpolate.
            # Linear coefficients: y = y0 + m0*x
            y0 = 1.0,
            m0 = 1.0,

            # Quadratic coefficients: y = y2 + m2*x^2
            y2 = 1.0,
            m2 = 0.5,

            gamma = 5.0/3.0,
            mu = 1.0,

            numGridLevels = 20,
            topGridCellSize = 20.0,

            # Parameters for iterating H.
            iterateH = True,
            maxHIterations = 100,
            Htolerance = 1.0e-4,

            # Parameters for passing the test
            interpolationTolerance = 5.0e-7,
            derivativeTolerance = 5.0e-5,

            graphics = True,
            plotKernels = False,
            )

assert testCase in ("constant", "linear", "quadratic")

#-------------------------------------------------------------------------------
# Which test case function are we doing?
#-------------------------------------------------------------------------------
if testCase == "constant":
    def func(x):
        return y0
    def dfunc(x):
        return 0.0
elif testCase == "linear":
    def func(x):
        return y0 + m0*x
    def dfunc(x):
        return m0
else:
    def func(x):
        return y2 + m2*x*x
    def dfunc(x):
        return 2*m2*x

#-------------------------------------------------------------------------------
# Create a random number generator.
#-------------------------------------------------------------------------------
import random
rangen = random.Random()
rangen.seed(seed)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))], nPerh = nPerh)
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
nodes1.massDensity(ScalarField("tmp", nodes1, rho1))

#-------------------------------------------------------------------------------
# Optionally randomly jitter the node positions.
#-------------------------------------------------------------------------------
dx = (x1 - x0)/nx1
for i in xrange(nodes1.numInternalNodes):
    nodes1.positions()[i].x += ranfrac * dx * rangen.uniform(-1.0, 1.0)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct some boundary conditions.
#-------------------------------------------------------------------------------
bounds = vector_of_Boundary()
xbc0 = ReflectingBoundary(Plane(Vector(x0), Vector( 1.0)))
xbc1 = ReflectingBoundary(Plane(Vector(x1), Vector(-1.0)))
bounds.append(xbc0)
bounds.append(xbc1)

#-------------------------------------------------------------------------------
# Iterate the h to convergence if requested.
#-------------------------------------------------------------------------------
if iterateH:
    method = SPHSmoothingScale()
    emptyBounds = vector_of_Boundary()
    iterateIdealH(db,
                  emptyBounds,
                  WT,
                  method,
                  maxHIterations,
                  Htolerance)

#-------------------------------------------------------------------------------
# Initialize our field.
#-------------------------------------------------------------------------------
f = ScalarField("test field", nodes1)
for i in xrange(nodes1.numInternalNodes):
    x = nodes1.positions()[i].x
    if testCase == "constant":
        f[i] = y0
    elif testCase == "linear":
        f[i] = y0 + m0*x
    elif testCase == "quadratic":
        f[i] = y2 + m2*x*x

#-------------------------------------------------------------------------------
# Build the tessellation.
#-------------------------------------------------------------------------------
mesh, void = generateLineMesh([nodes1],
                              bounds,
                              generateVoid = False,
                              removeBoundaryZones = False)

#-------------------------------------------------------------------------------
# Do that interpolation baby!
#-------------------------------------------------------------------------------
fl = ScalarFieldList()
fl.appendField(f)
db.updateConnectivityMap(True)
fSVPHfl = sampleFieldListSVPH(fl,
                              db.globalPosition,
                              db.globalHfield,
                              db.connectivityMap(),
                              WT,
                              mesh,
                              linearConsistent)
dfSVPHfl = gradientFieldListSVPH(fl,
                                 db.globalPosition,
                                 db.globalHfield,
                                 db.connectivityMap(),
                                 WT,
                                 mesh,
                                 linearConsistent)
fSVPH = fSVPHfl[0]
dfSVPH = dfSVPHfl[0]

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
positions = nodes1.positions()
xans = [positions[i].x for i in xrange(nodes1.numInternalNodes)]
yans = [func(x) for x in xans]
dyans = [dfunc(x) for x in xans]

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
ySVPH = fSVPH.internalValues()
dySVPH = [x.x for x in dfSVPH.internalValues()]

errySVPH = [y - z for y, z in zip(ySVPH, yans)]
maxySVPHerror = max([abs(x) for x in errySVPH])

errdySVPH = [y - z for y, z in zip(dySVPH, dyans)]
maxdySVPHerror = max([abs(x) for x in errdySVPH])

# errdySPH = [y - z for y, z in zip(dySPH, dyans)]
# maxdySPHerror = max([abs(x) for x in errdySPH])

print "Maximum errors (interpolation, gradient): SVPH = (%g, %g)" % (maxySVPHerror,
                                                                     maxdySVPHerror)

#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    import Gnuplot

    # Interpolated values.
    ansdata = Gnuplot.Data(xans, yans,
                           with_ = "lines",
                           title = "Answer",
                           inline = True)
    SVPHdata = Gnuplot.Data(xans, ySVPH,
                            with_ = "points",
                            title = "SVPH",
                            inline = True)
    errSVPHdata = Gnuplot.Data(xans, errySVPH,
                               with_ = "points",
                               title = "SVPH",
                               inline = True)

    p1 = generateNewGnuPlot()
    p1.plot(ansdata)
    p1.replot(SVPHdata)
    p1("set key top left")
    p1.title("Interpolated values")
    p1.refresh()

    p2 = generateNewGnuPlot()
    p2.replot(errSVPHdata)
    p2.title("Error in interpolation")
    p2.refresh()

    # Derivative values.
    dansdata = Gnuplot.Data(xans, dyans,
                            with_ = "lines",
                            title = "Answer",
                            inline = True)
    dSVPHdata = Gnuplot.Data(xans, dySVPH,
                             with_ = "points",
                             title = "SVPH",
                             inline = True)
    errdSVPHdata = Gnuplot.Data(xans, errdySVPH,
                                with_ = "points",
                                title = "SVPH",
                                inline = True)

    p3 = generateNewGnuPlot()
    p3.plot(dansdata)
    p3.replot(dSVPHdata)
    p3("set key top left")
    p3.title("Derivative values")
    p3.refresh()

    p4 = generateNewGnuPlot()
    p4.replot(errdSVPHdata)
    p4.title("Error in derivatives")
    p4.refresh()

#-------------------------------------------------------------------------------
# Check the maximum SVPH error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxySVPHerror > interpolationTolerance:
    raise ValueError, "SVPH interpolation error out of bounds: %g > %g" % (maxySVPHerror, interpolationTolerance)
if maxdySVPHerror > interpolationTolerance:
    raise ValueError, "SVPH gradient error out of bounds: %g > %g" % (maxdySVPHerror, interpolationTolerance)
