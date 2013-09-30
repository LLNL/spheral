#ATS:test(SELF, "--graphics False", label="CSPH interpolation test -- 1-D (serial)")
#-------------------------------------------------------------------------------
# A set of tests to compare how different meshless methods interpolate fields.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *

title("Interpolation tests")

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
            testDim = "1d",
            testCase = "linear",

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
            maxHIterations = 200,
            Htolerance = 1.0e-8,

            # Parameters for passing the test
            interpolationTolerance = 5.0e-7,
            derivativeTolerance = 5.0e-5,

            graphics = True,
            plotKernels = False,
            )

assert testCase in ("linear", "quadratic")
assert testDim in ("1d", "2d", "3d")

#-------------------------------------------------------------------------------
# Appropriately set generic object names based on the test dimensionality.
#-------------------------------------------------------------------------------
exec("from Spheral%s import *" % testDim)

## import Spheral
## for name in [x for x in Spheral.__dict__ if testDim in x]:
##     exec("%s = Spheral.__dict__['%s']" % (name.replace(testDim, ""), name))

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
    nodes1.positions()[i].x += ranfrac * rangen.uniform(-1.0, 1.0)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Iterate the h to convergence if requested.
#-------------------------------------------------------------------------------
if iterateH:
    bounds = vector_of_Boundary()
    method = SPHSmoothingScale()
    iterateIdealH(db,
                  bounds,
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
    if testCase == "linear":
        f[i] = y0 + m0*x
    elif testCase == "quadratic":
        f[i] = y2 + m2*x*x

#-------------------------------------------------------------------------------
# Prepare variables to accumulate the test values.
#-------------------------------------------------------------------------------
fSPH = ScalarField("SPH interpolated values", nodes1)
fCSPH = ScalarField("CSPH interpolated values", nodes1)
dfSPH = VectorField("SPH derivative values", nodes1)
dfCSPH = VectorField("CSPH derivative values", nodes1)

A_fl = db.newFluidScalarFieldList(0.0, "A")
B_fl = db.newFluidVectorFieldList(Vector.zero, "B")
C_fl = db.newFluidVectorFieldList(Vector.zero, "C")
D_fl = db.newFluidTensorFieldList(Tensor.zero, "D")
gradA_fl = db.newFluidVectorFieldList(Vector.zero, "gradA")
gradB_fl = db.newFluidTensorFieldList(Tensor.zero, "gradB")

db.updateConnectivityMap()
cm = db.connectivityMap()
position_fl = db.fluidPosition
weight_fl = db.fluidMass
H_fl = db.fluidHfield

computeCSPHCorrections(cm, WT, weight_fl, position_fl, H_fl,
                       A_fl, B_fl, C_fl, D_fl, gradA_fl, gradB_fl)

# Extract the field state for the following calculations.
positions = position_fl[0]
weight = weight_fl[0]
H = H_fl[0]
A = A_fl[0]
B = B_fl[0]
C = C_fl[0]
gradA = gradA_fl[0]
gradB = gradB_fl[0]

#-------------------------------------------------------------------------------
# Measure the interpolated values and gradients.
#-------------------------------------------------------------------------------
for i in xrange(nodes1.numInternalNodes):
    ri = positions[i]
    Hi = H[i]
    Hdeti = H[i].Determinant()
    wi = weight[i]
    Ai = A[i]
    Bi = B[i]
    Ci = C[i]
    gradAi = gradA[i]
    gradBi = gradB[i]
    fi = f[i]

    # Self contribution.
    W0 = WT.kernelValue(0.0, Hdeti)
    fSPH[i] = wi*W0 * fi
    fCSPH[i] = wi*W0*Ai * fi
    dfCSPH[i] = fi * wi*(Ai*Bi*W0 +
                         gradAi*W0)

    # Go over them neighbors.
    neighbors = cm.connectivityForNode(nodes1, i)
    assert len(neighbors) == 1
    for j in neighbors[0]:
        rj = positions[j]
        Hj = H[j]
        Hdetj = H[j].Determinant()
        wj = weight[j]
        Aj = A[j]
        Bj = B[j]
        fj = f[j]

        # The standard SPH kernel and it's gradient.
        rij = ri - rj
        etaj = Hj*rij
        Wj = WT.kernelValue(etaj.magnitude(), Hdetj)
        gradWj = Hj*etaj.unitVector() * WT.gradValue(etaj.magnitude(), Hdetj)

        # The corrected kernel and it's gradient.
        WRj = 0.0
        dummy = 0.0
        gradWRj = Vector()
        WRj, dummy = CSPHKernelAndGradient(WT,
                                           rij,
                                           etaj,
                                           Hj,
                                           Hdetj,
                                           Ai,
                                           Bi,
                                           gradAi,
                                           gradBi,
                                           gradWRj)
        assert fuzzyEqual(WRj, CSPHKernel(WT, rij, etaj, Hdetj, Ai, Bi), 1.0e-5)

        # Increment our interpolated values.
        fSPH[i] += fj * wj*Wj
        fCSPH[i] += fj * wj*WRj

        # Increment the derivatives.
        dfSPH[i] += fj * wj*gradWj
        dfCSPH[i] += fj * wj*gradWRj

    # We can now apply the integration correction (C) for CSPH.
    dfCSPH[i] += Ci*(fi - fCSPH[i])
 
#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
xans = [positions[i].x for i in xrange(nodes1.numInternalNodes)]
if testCase == "linear":
    yans = [y0 + m0*x for x in xans]
    dyans = [m0]*len(xans)
elif testCase == "quadratic":
    yans = [y2 + m2*x*x for x in xans]
    dyans = [2*m2*x for x in xans]

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
ySPH = fSPH.internalValues()
yCSPH = fCSPH.internalValues()

dySPH = [x.x for x in dfSPH.internalValues()]
dyCSPH = [x.x for x in dfCSPH.internalValues()]

errySPH = [y - z for y, z in zip(ySPH, yans)]
erryCSPH = [y - z for y, z in zip(yCSPH, yans)]
maxySPHerror = max([abs(x) for x in errySPH])
maxyCSPHerror = max([abs(x) for x in erryCSPH])

errdySPH = [y - z for y, z in zip(dySPH, dyans)]
errdyCSPH = [y - z for y, z in zip(dyCSPH, dyans)]
maxdySPHerror = max([abs(x) for x in errdySPH])
maxdyCSPHerror = max([abs(x) for x in errdyCSPH])

print "Maximum errors (interpolation): SPH = %g, CSPH = %g" % (maxySPHerror, maxyCSPHerror)
print "Maximum errors   (derivatives): SPH = %g, CSPH = %g" % (maxdySPHerror, maxdyCSPHerror)

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
    SPHdata = Gnuplot.Data(xans, ySPH,
                           with_ = "points",
                           title = "SPH",
                           inline = True)
    CSPHdata = Gnuplot.Data(xans, yCSPH,
                            with_ = "points",
                            title = "CSPH",
                            inline = True)
    errSPHdata = Gnuplot.Data(xans, errySPH,
                              with_ = "points",
                              title = "SPH",
                              inline = True)
    errCSPHdata = Gnuplot.Data(xans, erryCSPH,
                               with_ = "points",
                               title = "CSPH",
                               inline = True)

    p1 = generateNewGnuPlot()
    p1.plot(ansdata)
    p1.replot(SPHdata)
    p1.replot(CSPHdata)
    p1("set key top left")
    p1.title("Interpolated values")
    p1.refresh()

    p2 = generateNewGnuPlot()
    p2.plot(errSPHdata)
    p2.replot(errCSPHdata)
    p2.title("Error in interpolation")
    p2.refresh()

    # Derivative values.
    dansdata = Gnuplot.Data(xans, dyans,
                            with_ = "lines",
                            title = "Answer",
                            inline = True)
    dSPHdata = Gnuplot.Data(xans, dySPH,
                            with_ = "points",
                            title = "SPH",
                            inline = True)
    dCSPHdata = Gnuplot.Data(xans, dyCSPH,
                             with_ = "points",
                             title = "CSPH",
                             inline = True)
    errdSPHdata = Gnuplot.Data(xans, errdySPH,
                               with_ = "points",
                               title = "SPH",
                              inline = True)
    errdCSPHdata = Gnuplot.Data(xans, errdyCSPH,
                                with_ = "points",
                                title = "CSPH",
                                inline = True)

    p3 = generateNewGnuPlot()
    p3.plot(dansdata)
    p3.replot(dSPHdata)
    p3.replot(dCSPHdata)
    p3("set key top left")
    p3.title("Derivative values")
    p3.refresh()

    p4 = generateNewGnuPlot()
    p4.plot(errdSPHdata)
    p4.replot(errdCSPHdata)
    p4.title("Error in derivatives")
    p4.refresh()

if plotKernels:
    import Gnuplot
    pk = generateNewGnuPlot()
    for i in xrange(nodes1.numInternalNodes):
        xi = positions[i].x
        Hi = H[i]
        Hdeti = Hi.Determinant()
        hi = 1.0/Hi.xx
        Ai = A[i]
        Bi = B[i]

        dx = 2.0*kernelExtent*hi/50
        x = [xi - kernelExtent*hi + (i + 0.5)*dx for i in xrange(50)]
        y = [Ai*(1.0 + Bi.x*(xi - xj))*WT.kernelValue(abs(xi - xj)/hi, Hdeti) for xj in x]
        d = Gnuplot.Data(x, y, with_="lines", inline=True)
        pk.replot(d)

#-------------------------------------------------------------------------------
# Check the maximum CSPH error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxyCSPHerror > interpolationTolerance:
    raise ValueError, "CSPH interpolation error out of bounds: %g > %g" % (maxyCSPHerror, interpolationTolerance)

if maxdyCSPHerror > derivativeTolerance:
    raise ValueError, "CSPH derivative error out of bounds: %g > %g" % (maxdyCSPHerror, derivativeTolerance)
