#ATS:t0 = test(      SELF, "--graphics None --clearDirectories True  --checkError True   --restartStep 20", label="Planar Noh problem -- 1-D (serial)")
#ATS:t1 = testif(t0, SELF, "--graphics None --clearDirectories False --checkError False  --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem -- 1-D (serial) RESTART CHECK")
#ATS:t2 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-restartcheck' --restartStep 20", np=2, label="Planar Noh problem -- 1-D (parallel)")
#ATS:t3 = testif(t2, SELF, "--graphics None --clearDirectories False --checkError False --dataDir 'dumps-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=2, label="Planar Noh problem -- 1-D (parallel) RESTART CHECK")
#ATS:t4 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-1proc-reproducing.txt'", label="Planar Noh problem -- 1-D (serial reproducing test setup)")
#ATS:t5 = testif(t4, SELF, "--graphics None --clearDirectories False  --checkError True  --dataDir 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-4proc-reproducing.txt' --comparisonFile 'Noh-planar-1proc-reproducing.txt'", np=4, label="Planar Noh problem -- 1-D (4 proc reproducing test)")
#-------------------------------------------------------------------------------
# The Planar Noh test case run in 1-D.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
import os, shutil
from Spheral1d import *
from SpheralTestUtilities import *

title("1-D integrated hydro test -- planar Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(KernelConstructor = BSplineKernel,

            nx1 = 100,
            rho1 = 1.0,
            eps1 = 1.0,
            x0 = 0.0,
            x1 = 1.0,
            xwall = 0.5,
            nPerh = 1.25,
            NeighborType = NestedGridNeighbor,

            vr0 = -1.0, 
            vrSlope = 0.0,

            gamma = 5.0/3.0,
            mu = 1.0,

            SVPH = False,
            CSPH = False,
            #Qconstructor = MonaghanGingoldViscosity,
            Qconstructor = TensorMonaghanGingoldViscosity,
            boolReduceViscosity = False,
            nhQ = 5.0,
            nhL = 10.0,
            aMin = 0.1,
            aMax = 2.0,
            linearConsistent = False,
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            Qhmult = 1.0,
            Cl = 1.0, 
            Cq = 1.0,
            Qlimiter = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.3,
            nTensile = 4.0,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,
            filter = 0.0,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = False,
            updateBoundaryFrequency = 1,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IdealH,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = False,
            gradhCorrection = True,
            domainIndependent = True,
            cullGhostNodes = True,
            
            bArtificialConduction = True,
            arCondAlpha = 0.5,

            clearDirectories = True,
            checkError = True,
            checkRestart = False,
            checkEnergy = True,
            restoreCycle = None,
            restartStep = 10000,
            dataDir = "dumps-planar",
            restartBaseName = "Noh-planar-1d",
            outputFile = "None",
            comparisonFile = "None",

            # Parameters for the test
            scalePressure = 5.0,
            scaleEnergy = 5.0,

            graphics = "gnu",
            )

restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Noh-planar-1d-%i" % nx1)

dx = (x1 - x0)/nx1

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(KernelConstructor(), 1000)
WTPi = TableKernel(KernelConstructor(), 1000, Qhmult)
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           NeighborType = NeighborType)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))],
                         nPerh = nPerh)
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
nodes1.massDensity(ScalarField("tmp", nodes1, rho1))

nodeSet = [nodes1]

# Set node specific thermal energies
def specificEnergy(x):
    if x < xwall:
        return eps1
    else:
        return eps1*scaleEnergy
for nodes in nodeSet:
    pos = nodes.positions()
    eps = nodes.specificThermalEnergy()
    for i in xrange(nodes.numInternalNodes):
        eps[i] = specificEnergy(pos[i].x)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = SVPHFacetedHydro(WT, q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             densityUpdate = densityUpdate,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             generateVoid = False,
                             HUpdate = HUpdate,
                             fcentroidal = fcentroidal,
                             fcellPressure = fcellPressure,
                             xmin = Vector(-100.0),
                             xmax = Vector( 100.0))
elif CSPH:
    hydro = CSPHHydro(WT, WTPi, q,
                      filter = filter,
                      cfl = cfl,
                      compatibleEnergyEvolution = compatibleEnergy,
                      XSPH = XSPH,
                      densityUpdate = densityUpdate,
                      HUpdate = HUpdate)
else:
    hydro = SPHHydro(WT, WTPi, q,
                     cfl = cfl,
                     compatibleEnergyEvolution = compatibleEnergy,
                     gradhCorrection = gradhCorrection,
                     densityUpdate = densityUpdate,
                     HUpdate = HUpdate,
                     XSPH = XSPH,
                     epsTensile = epsilonTensile,
                     nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]


#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------

if boolReduceViscosity:
    #q.reducingViscosityCorrection = True
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(q,nhQ,nhL,aMin,aMax)
    
    packages.append(evolveReducingViscosityMultiplier)

#-------------------------------------------------------------------------------
# zero velocity package
#-------------------------------------------------------------------------------
class zeroV_pkg(Physics):
    def __init__(self):
        Physics.__init__(self)
        return

    def evaluateDerivatives(self, t, dt, db, state, derivs):
        DepsDt = derivs.scalarFields("delta " + HydroFieldNames.specificThermalEnergy)
        DvDt = derivs.vectorFields("delta " + HydroFieldNames.velocity)
        DepsDt.Zero()
        DvDt.Zero()
        return
    
    def dt(self, db, state, derivs, t):
        return pair_double_string(1e100, "No vote")
    
    def registerState(self, dt, state):
        return
    
    def registerDerivatives(self, db, derivs):
        return
    
    def label(self):
        return "zeroV package"
    
    def initialize(self, t, dt, db, state, derivs):
        return

zv = zeroV_pkg()

packages.append(zv)

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------

if bArtificialConduction:
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Optionally construct an hourglass control object.
#-------------------------------------------------------------------------------
if hourglass:
    mask = db.newFluidIntFieldList(1, "mask")
    pos = nodes1.positions()
    for i in xrange(nodes1.numInternalNodes):
        if pos[i].x > (x1 - dx):
            mask[0][i] = 0
    hg = hourglass(WT,
                   order = hourglassOrder,
                   limiter = hourglassLimiter,
                   fraction = hourglassFraction,
                   mask = mask)
    output("hg")
    output("hg.order")
    output("hg.limiter")
    output("hg.fraction")
    packages.append(hg)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
if x0 == xwall:
    xPlane0 = Plane(Vector(0.0), Vector(1.0))
    xbc0 = ReflectingBoundary(xPlane0)
    for p in packages:
        p.appendBoundary(xbc0)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
del p
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.updateBoundaryFrequency = updateBoundaryFrequency
integrator.domainDecompositionIndependent = domainIndependent
integrator.cullGhostNodes = cullGhostNodes
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.updateBoundaryFrequency")
output("integrator.domainDecompositionIndependent")
output("integrator.cullGhostNodes")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
output("control")

# Smooth the initial conditions.
if restoreCycle is None:
    control.smoothState(smoothIters)
    if densityUpdate in (VoronoiCellDensity, SumVoronoiCellDensity):
        print "Reinitializing node masses."
        control.voronoiInitializeMass()
##     rho = db.fluidMassDensity
##     pos = db.fluidPosition
##     mass = db.fluidMass
##     H = db.fluidHfield
##     db.updateConnectivityMap()
##     cm = db.connectivityMap()
##     computeSPHSumMassDensity(cm, WT, pos, mass, H, rho)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

    # Are we doing the restart test?
    if checkRestart:
        state0 = State(db, integrator.physicsPackages())
        state0.copyState()
        control.loadRestartFile(control.totalSteps)
        state1 = State(db, integrator.physicsPackages())
        if not state1 == state0:
            raise ValueError, "The restarted state does not match!"
        else:
            print "Restart check PASSED."

else:
    if control.time() < goalTime:
        control.step(5)
        control.advance(goalTime, maxSteps)


#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics == "matplot":
    import pylab
    from SpheralMatplotlibUtilities import *
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
    plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    plotEHistory(control.conserve)

elif graphics == "gnu":
    from SpheralGnuPlotUtilities import *
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
    plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    EPlot = plotEHistory(control.conserve)

    # Plot the specific entropy.
    Aplot = generateNewGnuPlot()
    AsimData = Gnuplot.Data(xprof, A,
                            with_ = "points",
                            title = "Simulation",
                            inline = True)
    AansData = Gnuplot.Data(xprof, Aans,
                            with_ = "lines",
                            title = "Solution",
                            inline = True)
    Aplot.plot(AsimData)
    Aplot.replot(AansData)
    Aplot.title("Specific entropy")
    Aplot.refresh()
    
    dvdxPlot = plotFieldList(hydro.DvDx(),yFunction='-1*%s.xx',winTitle='Source Fn',colorNodeLists=False)
    dudtPlot = plotFieldList(hydro.DepsDt(),yFunction='-1*%s.xx',winTitle='DepsDt',colorNodeLists=False)
    
    if boolReduceViscosity:
        alphaPlotQ = plotFieldList(q.reducingViscosityMultiplierQ(),
                                  winTitle = "rvAlphaQ",
                                  colorNodeLists = False, plotGhosts = False)
        alphaPlotL = plotFieldList(q.reducingViscosityMultiplierL(),
                                   winTitle = "rvAlphaL",
                                   colorNodeLists = False, plotGhosts = False)

    # # Plot the grad h correction term (omega)
    # omegaPlot = plotFieldList(hydro.omegaGradh(),
    #                           winTitle = "grad h correction",
    #                           colorNodeLists = False)

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print "Total energy error: %g" % Eerror
if checkEnergy and abs(Eerror) > 1e-13:
    raise ValueError, "Energy error outside allowed bounds."

#-------------------------------------------------------------------------------
# Measure the difference between the simulation and analytic answer.
#-------------------------------------------------------------------------------
rmin, rmax = 0.05, 0.35   # Throw away anything with r < rwall to avoid wall heating.
rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
P = ScalarField("pressure", nodes1)
nodes1.pressure(P)
Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile != "None":
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralGnuPlotUtilities import multiSort
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof)
        f = open(outputFile, "w")
        f.write(("#  " + 17*"'%s' " + "\n") % ("x", "rho", "P", "v", "eps", "h", "mo",
                                               "rhoans", "Pans", "vans", "hans",
                                               "x_UU", "rho_UU", "P_UU", "v_UU", "eps_UU", "h_UU"))
        for (xi, rhoi, Pi, vi, epsi, hi, mi,
             rhoansi, Pansi, vansi, hansi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                                                  rhoans, Pans, vans, hans):
            f.write((6*"%16.12e " + "%i " + 4*"%16.12e " + 6*"%i " + '\n') % 
                    (xi, rhoi, Pi, vi, epsi, hi, mi,
                     rhoansi, Pansi, vansi, hansi,
                     unpackElementUL(packElementDouble(xi)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi))))
        f.close()

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file.
        #---------------------------------------------------------------------------
        if comparisonFile != "None":
            comparisonFile = os.path.join(dataDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile, comparisonFile)
#------------------------------------------------------------------------------
# Compute the error.
#------------------------------------------------------------------------------
if checkError:
    if mpi.rank == 0:
        xans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
        import Pnorm
        print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
        failure = False
        for (name, data, ans,
             L1expect, L2expect, Linfexpect) in [("Mass Density", rhoprof, rhoans, L1rho, L2rho, Linfrho),
                                                 ("Pressure", Pprof, Pans, L1P, L2P, LinfP),
                                                 ("Velocity", vprof, vans, L1v, L2v, Linfv),
                                                 ("Thermal E", epsprof, epsans, L1eps, L2eps, Linfeps),
                                                 ("h       ", hprof, hans, L1h, L2h, Linfh)]:
            assert len(data) == len(ans)
            error = [data[i] - ans[i] for i in xrange(len(data))]
            Pn = Pnorm.Pnorm(error, xprof)
            L1 = Pn.gridpnorm(1, rmin, rmax)
            L2 = Pn.gridpnorm(2, rmin, rmax)
            Linf = Pn.gridpnorm("inf", rmin, rmax)
            print "\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf)
            if not fuzzyEqual(L1, L1expect, tol):
                print "L1 error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                      L1,
                                                                                      L1expect)
                failure = True
            if not fuzzyEqual(L2, L2expect, tol):
                print "L2 error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                      L2,
                                                                                      L2expect)
                failure = True
            if not fuzzyEqual(Linf, Linfexpect, tol):
                print "Linf error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                        Linf,
                                                                                        Linfexpect)
                failure = True
        if failure:
            raise ValueError, "Error bounds violated."

