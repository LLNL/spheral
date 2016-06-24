#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 1d --testCase linear", label="CRKSPH linear interpolation test -- 1D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 2d --testCase linear", label="CRKSPH linear interpolation test -- 2D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 3d --testCase linear", label="CRKSPH linear interpolation test -- 3D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 1d --testCase quadratic --correctionOrder QuadraticOrder", label="CRKSPH quadratic interpolation test -- 1D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 2d --testCase quadratic --correctionOrder QuadraticOrder", label="CRKSPH quadratic interpolation test -- 2D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 3d --testCase quadratic --correctionOrder QuadraticOrder", label="CRKSPH quadratic interpolation test -- 3D (serial)")
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
    nx1 = 50,     
    nx2 = 50,
    rho1 = 1.0,
    rho2 = 1.0,
    eps1 = 0.0,
    eps2 = 0.0,
    x0 = 0.0,
    x1 = 0.5,
    x2 = 1.0,
    nPerh = 2.01,
    hmin = 0.0001, 
    hmax = 10.0,

    # What order of reproducing kernel should we use (0,1,2)?
    correctionOrder = LinearOrder,
    
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

    # Parameters for iterating H.
    iterateH = True,
    maxHIterations = 200,
    Htolerance = 1.0e-4,

    # Parameters for passing the test
    interpolationTolerance = 5.0e-7,
    derivativeTolerance = 5.0e-5,

    graphics = True,
    plotKernels = False,
    outputFile = "None",
    plotSPH = True,
)

assert testCase in ("linear", "quadratic", "step")
assert testDim in ("1d", "2d", "3d")

FacetedVolume = {"1d" : Box1d,
                 "2d" : Polygon,
                 "3d" : Polyhedron}[testDim]

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
if testDim == "1d":
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes1, [(nx1, rho1, (x0, x1)),
                                        (nx2, rho2, (x1, x2))])], nPerh = nPerh)
elif testDim == "2d":
    from DistributeNodes import distributeNodes2d
    from GenerateNodeDistribution2d import GenerateNodeDistribution2d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen1 = GenerateNodeDistribution2d(nx1, nx1 + nx2, rho1,
                                      distributionType = "lattice",
                                      xmin = (x0, x0),
                                      xmax = (x1, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution2d(nx2, nx1 + nx2, rho2,
                                      distributionType = "lattice",
                                      xmin = (x1, x0),
                                      xmax = (x2, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes2d((nodes1, gen))

elif testDim == "3d":
    from DistributeNodes import distributeNodes3d
    from GenerateNodeDistribution3d import GenerateNodeDistribution3d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen1 = GenerateNodeDistribution3d(nx1, nx1 + nx2, nx1 + nx2, rho1,
                                      distributionType = "lattice",
                                      xmin = (x0, x0, x0),
                                      xmax = (x1, x2, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution3d(nx2, nx1 + nx2, nx1 + nx2, rho2,
                                      distributionType = "lattice",
                                      xmin = (x1, x0, x0),
                                      xmax = (x2, x2, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes3d((nodes1, gen))

else:
    raise ValueError, "Only tests cases for 1d,2d and 3d." 

output("nodes1.numNodes")

# Set node properties.
eps = nodes1.specificThermalEnergy()
for i in xrange(nx1):
    eps[i] = eps1
for i in xrange(nx2):
    eps[i + nx1] = eps2

#-------------------------------------------------------------------------------
# Optionally randomly jitter the node positions.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
dy = (x2 - x0)/(nx1 + nx2)
dz = (x2 - x0)/(nx1 + nx2)
pos = nodes1.positions()
for i in xrange(nodes1.numInternalNodes):
    if pos[i] < x1:
        dx = dx1
    else:
        dx = dx2
    if testDim == "1d":
        pos[i].x += ranfrac * dx * rangen.uniform(-1.0, 1.0)
    elif testDim == "2d":
        pos[i].x += ranfrac * dx * rangen.uniform(-1.0, 1.0)
        pos[i].y += ranfrac * dy * rangen.uniform(-1.0, 1.0)
    elif testDim == "3d":
        pos[i].x += ranfrac * dx * rangen.uniform(-1.0, 1.0)
        pos[i].y += ranfrac * dy * rangen.uniform(-1.0, 1.0)
        pos[i].z += ranfrac * dz * rangen.uniform(-1.0, 1.0)

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
    elif testCase == "step":
        if x < x1:
            f[i] = y0
        else:
            f[i] = 2*y0

#-------------------------------------------------------------------------------
# Prepare variables to accumulate the test values.
#-------------------------------------------------------------------------------
fSPH = ScalarField("SPH interpolated values", nodes1)
dfSPH = VectorField("SPH derivative values", nodes1)

A_fl = db.newFluidScalarFieldList(0.0, "A")
B_fl = db.newFluidVectorFieldList(Vector.zero, "B")
C_fl = db.newFluidTensorFieldList(Tensor.zero, "C")
gradA_fl = db.newFluidVectorFieldList(Vector.zero, "gradA")
gradB_fl = db.newFluidTensorFieldList(Tensor.zero, "gradB")
gradC_fl = db.newFluidThirdRankTensorFieldList(ThirdRankTensor.zero, "gradB")

m0_fl = db.newFluidScalarFieldList(0.0, "m0")
m1_fl = db.newFluidVectorFieldList(Vector.zero, "m1")
m2_fl = db.newFluidSymTensorFieldList(SymTensor.zero, "m2")
m3_fl = db.newFluidThirdRankTensorFieldList(ThirdRankTensor.zero, "m3")
m4_fl = db.newFluidFourthRankTensorFieldList(FourthRankTensor.zero, "m4")
gradm0_fl = db.newFluidVectorFieldList(Vector.zero, "grad m0")
gradm1_fl = db.newFluidTensorFieldList(Tensor.zero, "grad m1")
gradm2_fl = db.newFluidThirdRankTensorFieldList(ThirdRankTensor.zero, "grad m2")
gradm3_fl = db.newFluidFourthRankTensorFieldList(FourthRankTensor.zero, "grad m3")
gradm4_fl = db.newFluidFifthRankTensorFieldList(FifthRankTensor.zero, "grad m4")

db.updateConnectivityMap(True)
cm = db.connectivityMap()
position_fl = db.fluidPosition
weight_fl = db.fluidMass
H_fl = db.fluidHfield

# Compute the volumes to use as weighting.
polyvol_fl = db.newFluidFacetedVolumeFieldList(FacetedVolume(), "polyvols")
#weight_fl = db.newFluidScalarFieldList(1.0, "volume")
#computeHullVolumes(cm, position_fl, polyvol_fl, weight_fl)
computeCRKSPHMoments(cm, WT, weight_fl, position_fl, H_fl, correctionOrder, NodeCoupling(),
                     m0_fl, m1_fl, m2_fl, m3_fl, m4_fl, gradm0_fl, gradm1_fl, gradm2_fl, gradm3_fl, gradm4_fl)
computeCRKSPHCorrections(m0_fl, m1_fl, m2_fl, m3_fl, m4_fl, gradm0_fl, gradm1_fl, gradm2_fl, gradm3_fl, gradm4_fl, H_fl,
                         correctionOrder,
                         A_fl, B_fl, C_fl, gradA_fl, gradB_fl, gradC_fl)

# Extract the field state for the following calculations.
positions = position_fl[0]
weight = weight_fl[0]
H = H_fl[0]
A = A_fl[0]
B = B_fl[0]
C = C_fl[0]
gradA = gradA_fl[0]
gradB = gradB_fl[0]
gradC = gradC_fl[0]

#-------------------------------------------------------------------------------
# Measure the interpolated values and gradients.
#-------------------------------------------------------------------------------
for i in xrange(nodes1.numInternalNodes):
    ri = positions[i]
    Hi = H[i]
    Hdeti = H[i].Determinant()
    wi = weight[i]
    fi = f[i]

    # Self contribution.
    W0 = WT.kernelValue(0.0, Hdeti)
    fSPH[i] = wi*W0 * fi

    # Go over them neighbors.
    neighbors = cm.connectivityForNode(nodes1, i)
    assert len(neighbors) == 1
    for j in neighbors[0]:
        rj = positions[j]
        Hj = H[j]
        Hdetj = H[j].Determinant()
        wj = weight[j]
        fj = f[j]

        # The standard SPH kernel and it's gradient.
        rij = ri - rj
        etai = Hi*rij
        etaj = Hj*rij
        Wj = WT.kernelValue(etaj.magnitude(), Hdetj)
        gradWj = Hj*etaj.unitVector() * WT.gradValue(etaj.magnitude(), Hdetj)

        # Increment our interpolated values.
        fSPH[i] += fj * wj*Wj

        # Increment the derivatives.
        dfSPH[i] += fj * wj*gradWj

#-------------------------------------------------------------------------------
# Check the C++ interpolation and gradient methods.
#-------------------------------------------------------------------------------
f_fl = ScalarFieldList()
f_fl.appendField(f)
fCRKSPH_fl = interpolateCRKSPH(f_fl, position_fl, weight_fl, H_fl, A_fl, B_fl, C_fl,
                               cm, correctionOrder, WT)
dfCRKSPH_fl = gradientCRKSPH(f_fl, position_fl, weight_fl, H_fl,
                             A_fl, B_fl, C_fl, gradA_fl, gradB_fl, gradC_fl,
                             cm, correctionOrder, WT)
fCRKSPH = fCRKSPH_fl[0]
dfCRKSPH = dfCRKSPH_fl[0]

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
xans = [positions[i].x for i in xrange(nodes1.numInternalNodes)]
yans = ScalarField("interpolation answer", nodes1)
dyans = ScalarField("derivative answer", nodes1)
for i in xrange(nodes1.numInternalNodes):
    if testCase == "linear":
        yans[i] = y0 + m0*xans[i]
        dyans[i] = m0
    elif testCase == "quadratic":
        yans[i] = y2 + m2*xans[i]*xans[i]
        dyans[i] = 2*m2*xans[i]
    elif testCase == "step":
        if i < nx1:
            yans[i] = y0
        else:
            yans[i] = 2*y0
        dyans[i] = 0.0

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
errySPH =   ScalarField("SPH interpolation error", nodes1)
erryCRKSPH =  ScalarField("CRKSPH interpolation error", nodes1)
errdySPH =  ScalarField("SPH derivative error", nodes1)
errdyCRKSPH = ScalarField("CRKSPH derivative error", nodes1)
for i in xrange(nodes1.numInternalNodes):
    errySPH[i] =   fSPH[i] - yans[i]
    erryCRKSPH[i] =  fCRKSPH[i] - yans[i]
    errdySPH[i] =  dfSPH[i].x - dyans[i]
    errdyCRKSPH[i] = dfCRKSPH[i].x - dyans[i]

maxySPHerror = max([abs(x) for x in errySPH])
maxdySPHerror = max([abs(x) for x in errdySPH])
maxyCRKSPHerror = max([abs(x) for x in erryCRKSPH])
maxdyCRKSPHerror = max([abs(x) for x in errdyCRKSPH])

print "Maximum errors (interpolation): SPH = %g, CRKSPH = %g" % (maxySPHerror, maxyCRKSPHerror)
print "Maximum errors   (derivatives): SPH = %g, CRKSPH = %g" % (maxdySPHerror, maxdyCRKSPHerror)

#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    import Gnuplot
    xans = [positions[i].x for i in xrange(nodes1.numInternalNodes)]

    # Interpolated values.
    ansdata = Gnuplot.Data(xans, yans.internalValues(),
                           with_ = "lines",
                           title = "Answer",
                           inline = True)
    SPHdata = Gnuplot.Data(xans, fSPH.internalValues(),
                           with_ = "points",
                           title = "SPH",
                           inline = True)
    CRKSPHdata = Gnuplot.Data(xans, fCRKSPH.internalValues(),
                            with_ = "points",
                            title = "CRKSPH",
                            inline = True)
    errSPHdata = Gnuplot.Data(xans, errySPH.internalValues(),
                              with_ = "points",
                              title = "SPH",
                              inline = True)
    errCRKSPHdata = Gnuplot.Data(xans, erryCRKSPH.internalValues(),
                               with_ = "points",
                               title = "CRKSPH",
                               inline = True)

    p1 = generateNewGnuPlot()
    p1.plot(ansdata)
    if plotSPH:
     p1.replot(SPHdata)
    p1.replot(CRKSPHdata)
    p1("set key top left")
    p1.title("Interpolated values")
    p1.refresh()

    p2 = generateNewGnuPlot()
    if plotSPH:
     p2.plot(errSPHdata)
    p2.replot(errCRKSPHdata)
    p2.title("Error in interpolation")
    p2.refresh()

    # Derivative values.
    dansdata = Gnuplot.Data(xans, dyans.internalValues(),
                            with_ = "lines",
                            title = "Answer",
                            inline = True)
    dSPHdata = Gnuplot.Data(xans, [x.x for x in dfSPH.internalValues()],
                            with_ = "points",
                            title = "SPH",
                            inline = True)
    dCRKSPHdata = Gnuplot.Data(xans, [x.x for x in dfCRKSPH.internalValues()],
                             with_ = "points",
                             title = "CRKSPH",
                             inline = True)
    errdSPHdata = Gnuplot.Data(xans, errdySPH.internalValues(),
                               with_ = "points",
                               title = "SPH",
                              inline = True)
    errdCRKSPHdata = Gnuplot.Data(xans, errdyCRKSPH.internalValues(),
                                with_ = "points",
                                title = "CRKSPH",
                                inline = True)

    p3 = generateNewGnuPlot()
    p3.plot(dansdata)
    if plotSPH:
     p3.replot(dSPHdata)
    p3.replot(dCRKSPHdata)
    p3("set key top left")
    p3.title("Derivative values")
    p3.refresh()

    p4 = generateNewGnuPlot()
    if plotSPH:
     p4.plot(errdSPHdata)
    p4.replot(errdCRKSPHdata)
    p4.title("Error in derivatives")
    p4.refresh()

    # Plot the kernel shapes as appropriate.
    if testDim == "1d":
        p7 = generateNewGnuPlot()
        j = -2 # int(nodes1.numInternalNodes/2)
        Hj = H[j]
        hj = 1.0/Hj.xx
        Hdetj = H[j].Determinant()
        Aj = A[j]
        Bj = B[j].x
        Cj = C[j].xx
        nsamp = 100
        dx = 4.0/nsamp
        W = [WT.kernelValue(abs(i*dx - 2.0), Hdetj) for i in xrange(nsamp)]
        #WR = [x*Aj*(1.0 + Bj*(2.0 - i*dx)*hj) for i, x in enumerate(W)]
        WR = [x*Aj*(1.0 + Bj*(2.0 - i*dx)*hj+Cj*(2.0 - i*dx)*(2.0 - i*dx)*hj*hj) for i, x in enumerate(W)]
        p7.plot(W)
        p7.replot(WR)
        p7.title("Kernel")
        p7.refresh()
        if outputFile != "None":
            f = open("Kernel_" + outputFile, "w")
            f.write(("#" + 3*' "%20s"' + "\n") % ("eta", "Wj", "WRj"))
            for i in xrange(nsamp):
                f.write((3*" %20g" + "\n") % ((i*dx - 2.0), W[i], WR[i]))
            f.close()

    # We may want a gnu/pdv style text file.
    if outputFile != "None" and testDim == "2d":
        of = open(outputFile, "w")
        of.write(('#' + 7*' "%20s"' + '\n') % ("x", "interp answer", "grad answer", "interp SPH", "interp CRK", "grad SPH", "grad CRK"))
        for i in xrange(nodes1.numInternalNodes):
            of.write((7*" %20g" + "\n") %
                    (xans[i], yans[i], dyans[i], fSPH[i], fCRKSPH[i], dfSPH[i].x, dfCRKSPH[i].x))
        of.close()

    # If we're in 2D dump a silo file too.
    if testDim == "2d":
        # from SpheralVoronoiSiloDump import SpheralVoronoiSiloDump
        # dumper = SpheralVoronoiSiloDump("testInterpolation_%s_2d" % testCase,
        #                                 listOfFields = [fSPH, fCRKSPH, dfSPH, dfCRKSPH,
        #                                                 yans, dyans,
        #                                                 errySPH, erryCRKSPH, errdySPH, errdyCRKSPH],
        #                                 listOfFieldLists = [weight_fl, 
        #                                                     A_fl, B_fl, gradA_fl, gradB_fl,
        #                                                     dfCRKSPH_fl])
        # dumper.dump(0.0, 0)
        from siloPointmeshDump import siloPointmeshDump
        siloPointmeshDump("testInterpolation_%s_2d" % testCase,
                          fields = [fSPH, fCRKSPH, dfSPH, dfCRKSPH,
                                    yans, dyans,
                                    errySPH, erryCRKSPH, errdySPH, errdyCRKSPH],
                          fieldLists = [weight_fl, A_fl, B_fl, gradA_fl, gradB_fl])

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
        Ci = C[i]

        dx = 2.0*kernelExtent*hi/50
        x = [xi - kernelExtent*hi + (i + 0.5)*dx for i in xrange(50)]
        #y = [Ai*(1.0 + Bi.x*(xi - xj))*WT.kernelValue(abs(xi - xj)/hi, Hdeti) for xj in x]
        y = [Ai*(1.0 + Bi.x*(xi - xj)+Ci.xx*(xi-xj)*(xi-xj))*WT.kernelValue(abs(xi - xj)/hi, Hdeti) for xj in x]
        d = Gnuplot.Data(x, y, with_="lines", inline=True)
        pk.replot(d)

#-------------------------------------------------------------------------------
# Check the maximum CRKSPH error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxyCRKSPHerror > interpolationTolerance:
    raise ValueError, "CRKSPH interpolation error out of bounds: %g > %g" % (maxyCRKSPHerror, interpolationTolerance)

if maxdyCRKSPHerror > derivativeTolerance:
    raise ValueError, "CRKSPH derivative error out of bounds: %g > %g" % (maxdyCRKSPHerror, derivativeTolerance)
