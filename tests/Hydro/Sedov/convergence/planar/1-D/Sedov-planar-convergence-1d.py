#-------------------------------------------------------------------------------
# The planar Sedov test case (1-D).
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
import os, sys

import loadmpi
mpi, rank, nProcs = loadmpi.loadmpi()

title("1-D integrated hydro test -- planar Sedov problem")
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList1d,

            nxlist = [51, 101, 201, 401, 801, 1601, 3201, 6401, 12801, 25602],
            rho1 = 1.0,
            eps1 = 0.0,
            x0 = -1.0,
            x1 = 1.0,
            nPerh = 1.01,

            xSpike = 0.0,
            Espike = 1.0,
            smoothSpike = True,
            gamma = 5.0/3.0,
            mu = 1.0,

            Cl = 1.0,
            Cq = 0.75,
            epsilon2 = 1e-2,

            hmin = 1e-10,
            hmax = 1.0,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = True,
            rhomin = 1e-10,

            neighborSearchType = Neighbor1d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 1.0,
            origin = Vector1d(0.0),

            goalTime = 0.3,
            dt = 1e-8,
            dtMin = 1.0e-8,
            dtMax = None,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = Hydro1d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro1d.MassDensityType.RigorousSumDensity,
            compatibleEnergy = True,

            restoreCycle = None,
            restartStep = 1000,

            dataDirBase = "Sedov-planar-convergence-1d",

            graphics = False,
            )

if graphics and graphics != True:
    nxlist = [graphics]

style = "points ps 2"
lineTitle = "Simulation"
if graphics == True:
    style = "lines"
    lineTitle = ""

dataDir = os.path.join(dataDirBase,
                       "%s" % str(NodeListConstructor).split("'")[1],
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "n=%i")
restartBaseName = os.path.join(dataDir, "restart")

#-------------------------------------------------------------------------------
# Create the file we're going to record the error norms in.
#-------------------------------------------------------------------------------
pnormFileName = "Sedov-planar-convergence-test-nperh=%4.2f-XSPH=%s-compatibleEnergy=%s.txt" % (nPerh, XSPH, compatibleEnergy)
if mpi.rank == 0:
    resultFile = open(pnormFileName, "w")
    resultFile.write("#             Mass Density                         Pressure                             Velocity                             Specific Thermal Energy              Specific Entropy                     h\n")
    resultFile.write("# N nodes     L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         \n")

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS1d(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel1d(BSplineKernel1d(), 1000)
WTPi = WT
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Create a NodeList and associated Neighbor object.
#-------------------------------------------------------------------------------
nodes1 = SphNodeList1d("nodes1", eos, WT, WTPi)
nodes1.hmin = hmin
nodes1.hmax = hmax
nodes1.nodesPerSmoothingScale = nPerh
nodes1.rhoMin = rhomin

neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase1d()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct a standard Monaghan-Gingold artificial viscosity.
#-------------------------------------------------------------------------------
qMG = MonaghanGingoldViscosity1d(Cl, Cq)
qMG.epsilon2 = epsilon2
output("qMG")
output("qMG.Cl")
output("qMG.Cq")
output("qMG.epsilon2")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro1d(WT, WTPi, qMG, compatibleEnergy)
hydro.cfl = cfl
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
hydro.sumForMassDensity = sumForMassDensity
output("hydro")
output("hydro.compatibleEnergyEvolution")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.sumForMassDensity")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
## xPlane0 = Plane1d(Vector1d(x0), Vector1d(1.0))
## xbc0 = ReflectingBoundary1d(xPlane0)
## hydro.appendBoundary(xbc0)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator1d(db) # PredictorCorrectorIntegrator1d(db) # 
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.dtGrowth")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.valid()")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output("control")

#-------------------------------------------------------------------------------
# If we're plotting graphics, prepare Gnuplot.
#-------------------------------------------------------------------------------
if graphics:
    prho = generateNewGnuPlot()
    #prho.xlabel("x")
    #prho.ylabel("{/Symbol r}")

    pvel = generateNewGnuPlot()
    #pvel.xlabel("x")
    #pvel.ylabel("v_x")

    pA = generateNewGnuPlot()
    #pA.xlabel("x")
    #pA.ylabel("A")

#-------------------------------------------------------------------------------
# Iterate over each resolution, do the simulation, and measure the error.
#-------------------------------------------------------------------------------
for nx1 in nxlist:

    #---------------------------------------------------------------------------
    # Explicitly initialize.
    #---------------------------------------------------------------------------
    nodes1.numInternalNodes = 0
    nodes1.numGhostNodes = 0
    dx = (x1 - x0)/nx1
    h1 = nPerh*dx

    #---------------------------------------------------------------------------
    # Check if the necessary output directories exist.  If not, create them.
    #---------------------------------------------------------------------------
    restartDir = dataDir % nx1
    if mpi.rank == 0:
        if not os.path.exists(restartDir):
            os.makedirs(restartDir)
    mpi.barrier()
    restartName = restartBaseName % nx1

    #---------------------------------------------------------------------------
    # Check if there is already available restart info.
    #---------------------------------------------------------------------------
    from findLastRestart import findLastRestart
    restoreCycle = findLastRestart(restartName)

    #---------------------------------------------------------------------------
    # Set the NodeList properties.
    #---------------------------------------------------------------------------
    if restoreCycle is None:

        from DistributeNodes import distributeNodesInRange1d
        distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))])
        output("mpi.allreduce(nodes1.numNodes, mpi.MIN)")
        output("mpi.allreduce(nodes1.numNodes, mpi.MAX)")
        output("mpi.allreduce(nodes1.numNodes, mpi.SUM)")

        integrator.lastDt = dt

        # Set the point source of energy.
        Esum = 0.0
        if smoothSpike:
            Wsum = 0.0
            for nodeID in xrange(nodes1.numInternalNodes):
                Hi = nodes1.Hfield()[nodeID].xx / 2.0
                rij = abs(nodes1.positions()[nodeID].x - xSpike)*Hi
                Wi = WT.kernelValue(rij, Hi)
                Ei = Wi*Espike
                epsi = Ei/nodes1.mass()[nodeID]
                nodes1.specificThermalEnergy()[nodeID] = epsi
                Wsum += Wi
            Wsum = mpi.allreduce(Wsum, mpi.SUM)
            assert Wsum > 0.0
            for nodeID in xrange(nodes1.numInternalNodes):
                nodes1.specificThermalEnergy()[nodeID] /= Wsum
                Esum += nodes1.specificThermalEnergy()[nodeID]*nodes1.mass()[nodeID]
        else:
            i = -1
            rmin = 1e50
            for nodeID in xrange(nodes1.numInternalNodes):
                rij = abs(nodes1.positions()[nodeID].x - xSpike)
                if rij < rmin:
                    i = nodeID
                    rmin = rij
            rminglobal = mpi.allreduce(rmin, mpi.MIN)
            if fuzzyEqual(rmin, rminglobal):
                assert i >= 0 and i < nodes1.numInternalNodes
                nodes1.specificThermalEnergy()[i] = Espike/nodes1.mass()[i]
                Esum += Espike
        Eglobal = mpi.allreduce(Esum, mpi.SUM)
        print "Initialized a total energy of", Eglobal
        assert smoothSpike or fuzzyEqual(Eglobal, Espike)

        # Use the controller to reinitialize the problem.
        control.reinitializeProblem(restartName,
                                    statsStep = statsStep,
                                    restartStep = restartStep)

        # Iterate the H to convergence.
        control.iterateIdealH()
        db.updateFluidMassDensity(WT)

##         # This bit of craziness is needed to try and get the intial mass density
##         # as exactly as possible.
##         for i in xrange(10):
##             thpt = mpi.allreduce(max(nodes1.massDensity().internalValues()), mpi.MAX)
##             print i, thpt
##             for j in xrange(nodes1.numInternalNodes):
##                 nodes1.mass()[j] *= rho1/thpt
##                 nodes1.specificThermalEnergy()[j] *= thpt/rho1
##             state = State1d(db, integrator.physicsPackages())
##             derivs = StateDerivatives1d(db, integrator.physicsPackages())
##             integrator.initialize(state, derivs)
##             db.updateFluidMassDensity(WT)

        # Smooth the initial conditions.
        control.smoothState(smoothIters)

    else:
        # Load the stored state.
        control.reinitializeProblem(restartName)
        control.loadRestartFile(restoreCycle)

    #---------------------------------------------------------------------------
    # Advance to the end time.
    #---------------------------------------------------------------------------
    control.advance(goalTime, maxSteps)
    if restoreCycle < control.totalSteps:
        control.dropRestartFile()

    # Report the energy conservation.
    print "Energy drift: ", ((control.conserve.EHistory[-1] -
                              control.conserve.EHistory[0])/
                             control.conserve.EHistory[0])

    #---------------------------------------------------------------------------
    # Compute the error.
    #---------------------------------------------------------------------------
    # Compute the analytic solution.
    sys.path.append("../../..")
    import SedovAnalyticSolution
    h1 = (x1 - x0)/nx1*nPerh
    answer = SedovAnalyticSolution.SedovSolution(nDim = 1,
                                                 gamma = gamma,
                                                 rho0 = rho1,
                                                 E0 = Espike,
                                                 h0 = h1,
                                                 nbins = 50001)

    rmin = 0.3 # 0.05 # -0.9
    rmax =  0.9
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(nodes1.pressure().internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
    xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)

    # Compute the specific entropy.
    Afield = ScalarField1d("Specific entropy", nodes1)
    rho = nodes1.massDensity().internalValues()
    P = nodes1.pressure().internalValues()
    for i in xrange(nodes1.numInternalNodes):
        assert rho[i] > 0.0
        Afield[i] = P[i]/(rho[i]**gamma)
    Aprof = mpi.reduce(Afield.internalValues(), mpi.SUM)

    if mpi.rank == 0:
        resultFile.write("%10i " % nx1)
        rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
        Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]
        import Pnorm
        print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
        for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                                  ("Pressure", Pprof, Pans),
                                  ("Velocity", vprof, vans),
                                  ("Thermal E", epsprof, epsans),
                                  ("Entropy ", Aprof, Aans),
                                  ("h       ", hprof, hans)]:
            assert len(data) == len(ans)
            error = [data[i] - ans[i] for i in xrange(len(data))]
            Pn = Pnorm.Pnorm(error, xprof)
            print "\t%s \t\t%g \t\t%g \t\t%g" % (name,
                                                 Pn.pnormAverage(1, rmin = rmin, rmax = rmax),
                                                 Pn.pnormAverage(2, rmin = rmin, rmax = rmax),
                                                 Pn.pnormAverage("inf", rmin = rmin, rmax = rmax))
            resultFile.write("%11.6e %11.6e %11.6e " % (Pn.pnormAverage(1, rmin = rmin, rmax = rmax),
                                                        Pn.pnormAverage(2, rmin = rmin, rmax = rmax),
                                                        Pn.pnormAverage("inf", rmin = rmin, rmax = rmax)))
        resultFile.write("\n")
        resultFile.flush()

    #---------------------------------------------------------------------------
    # Throw up the density, velocity, and entropy profiles if we're
    # generating graphics.
    #---------------------------------------------------------------------------
    if graphics:
        rho = ScalarFieldList1d()
        rho.appendField(nodes1.massDensity())
        plotFieldList(rho,
                      plot = prho,
                      lineTitle = lineTitle, #"%i" % nx1,
                      plotStyle = style,
                      colorDomains = False,
                      colorNodeLists = False)

        vel = VectorFieldList1d()
        vel.appendField(nodes1.velocity())
        plotFieldList(vel,
                      yFunction = "%s.x",
                      plot = pvel,
                      lineTitle = lineTitle, #"%i" % nx1,
                      plotStyle = style,
                      colorDomains = False,
                      colorNodeLists = False)

        Afl = ScalarFieldList1d()
        Afl.appendField(Afield)
        plotFieldList(Afl,
                      plot = pA,
                      lineTitle = lineTitle, #"%i" % nx1,
                      plotStyle = style,
                      colorDomains = False,
                      colorNodeLists = False)

#-------------------------------------------------------------------------------
# Plot the answers.
#-------------------------------------------------------------------------------
if mpi.rank == 0 and graphics:
    drhoans = Gnuplot.Data(xprof, rhoans,
                           with = "lines lw 2",
                           title = "analytic")
    prho.replot(drhoans)
    dvelans = Gnuplot.Data(xprof, vans,
                           with = "lines lw 2",
                           title = "analytic")
    pvel.replot(dvelans)
    dAans = Gnuplot.Data(xprof, Aans,
                         with = "lines lw 2",
                         title = "analytic")
    pA.replot(dAans)

    if compatibleEnergy:
        tag = "compatibleEnergy"
    else:
        tag = "standard"

    if graphics == True:
        prho("set xrange [0.45:0.55]; set yrange [0.98:]; set key bottom left"); prho.refresh()
        pvel("set xrange [0.45:0.55]; set yrange [0:]; set key bottom left"); pvel.refresh()
        pA("set xrange [0.45:0.55]; set yrange [0.01:500]; set logscale y"); pA.refresh()
        tag1 = "-zoom"
    else:
        prho("set xrange [0:0.55]; set key top left"); prho.refresh()
        pvel("set xrange [0:0.55]; set yrange [0:]; set key top left"); pvel.refresh()
        pA("set xrange [0:0.55]; set yrange [0.01:500]; set logscale y"); pA.refresh()
        tag1 = ""

    prho.hardcopy("Sedov-planar-1d-%s-rho-profiles%s.eps" % (tag, tag1), eps=True, color=False, fontsize=24)
    pvel.hardcopy("Sedov-planar-1d-%s-vel-profiles%s.eps" % (tag, tag1), eps=True, color=False, fontsize=24)
    pA.hardcopy("Sedov-planar-1d-%s-A-profiles%s.eps" % (tag, tag1), eps=True, color=False, fontsize=24)

if mpi.rank == 0:
    resultFile.close()
