#-------------------------------------------------------------------------------
# The Planar Noh test case run in 1-D.
# This script is designed to automatically measure the convergence rate.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
import os, sys
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
import Gnuplot

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList1d,

            nxlist = [25, 50, 100, 200, 400, 800, 1600, 3200, 6400],
            nx2 = 10,
            rho1 = 1.0,
            eps1 = 0.0,
            x0 = 0.0,
            x1 = 0.5,
            nPerh = 2.01,

            vx0 = -1.0, 

            gamma = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosity1d,
            #Qconstructor = TensorMonaghanGingoldViscosity1d,
            Cl = 1.5, 
            Cq = 1.5,
            Qlimiter = False,
            epsilon2 = 1e-2,
            hmin = 1e-10,
            hmax = 1.0,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            hourglassAccelerationFactor = 0.0,

            neighborSearchType = Neighbor1d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector1d(0.0),

            goalTime = 0.3,
            dt = 1e-6,
            dtMin = 1.0e-6, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 100,
            smoothIters = 0,
            HEvolution = Hydro1d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro1d.MassDensityType.RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,

            restoreCycle = None,
            restartStep = 10000,
            
            dataDirBase = "Noh-planar-convergence-1d",
            restartBaseName = "Noh-planar-1d-%i",

            graphics = False,
            )

if graphics:
    nxlist = [graphics]

# Base restart dir.
dataDirBase = os.path.join(dataDirBase,
                           "%s" % str(NodeListConstructor).split("'")[1],
                           "nperh=%4.2f" % nPerh,
                           "XSPH=%s" % XSPH,
                           "compatibleEnergy=%s" % compatibleEnergy,
                           "n=%i")
restartDirBase = os.path.join(dataDirBase, "restart")

#-------------------------------------------------------------------------------
# Create the file we're going to record the error norms in.
#-------------------------------------------------------------------------------
pnormFileName = "Noh-planar-convergence-test-nperh=%4.2f-XSPH=%s-compatibleEnergy=%s.txt" % (nPerh, XSPH, compatibleEnergy)
if mpi.rank == 0:
    resultFile = open(pnormFileName, "w")
    resultFile.write("#             Mass Density                         Pressure                             Velocity                             Specific Thermal Energy              Specific Entropy                     h\n")
    resultFile.write("# N nodes     L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         \n")

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS1d(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel1d(BSplineKernel1d(), 1000)
WTPi = TableKernel1d(BSplineKernel1d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes1 = NodeListConstructor("nodes1", eos, WT, WTPi)
nodes2 = NodeListConstructor("nodes2", eos, WT, WTPi)
nodeSet = [nodes1, nodes2]
for n in nodeSet:
    n.XSPH = XSPH
    n.hmin = hmin
    n.hmax = hmax
    n.nodesPerSmoothingScale = nPerh
    n.epsilonTensile = epsilonTensile
    n.nTensile = nTensile
del n

#-------------------------------------------------------------------------------
# Construct the neighbor objects.
#-------------------------------------------------------------------------------
cache = []
for n in nodeSet:
    neighbor = NestedGridNeighbor1d(n,
                                    neighborSearchType,
                                    numGridLevels,
                                    topGridCellSize,
                                    origin,
                                    kernelExtent)
    n.registerNeighbor(neighbor)
    cache.append(neighbor)
del n

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase1d()
for n in nodeSet:
    db.appendNodeList(n)
del n

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro1d(WT, WTPi, q, compatibleEnergy)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane1d(Vector1d(x0), Vector1d(1.0))
xbc0 = ReflectingBoundary1d(xPlane0)
hydro.appendBoundary(xbc0)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator1d(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)

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
    for n in nodeSet:
        n.numInternalNodes = 0
        n.numGhostNodes = 0
    del n
    dx = (x1 - x0)/nx1
    h1 = nPerh*dx

    #---------------------------------------------------------------------------
    # Build the restart file name and directory.
    #---------------------------------------------------------------------------
    restartDir = restartDirBase % nx1
    if mpi.rank == 0:
        if not os.path.exists(restartDir):
            os.makedirs(restartDir)
    mpi.barrier()
    restartName = os.path.join(restartDir, restartBaseName % nx1)

    #---------------------------------------------------------------------------
    # Check if there is already available restart info.
    #---------------------------------------------------------------------------
    from findLastRestart import findLastRestart
    restoreCycle = findLastRestart(restartName)

    #---------------------------------------------------------------------------
    # If we're not restarting, then initialize the problem state.
    #---------------------------------------------------------------------------
    if restoreCycle is None:

        # Problem parameters
        x2 = x1 + (x1 - x0)/float(nx1)*float(nx2)
    
        # Problem setup
        # Set node positions for this domain
        from DistributeNodes import distributeNodesInRange1d
        distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1)),
                                  (nodes2, nx2, rho1, (x1, x2))])
        output("mpi.allreduce(nodes1.numNodes, mpi.SUM)")
        output("mpi.allreduce(nodes2.numNodes, mpi.SUM)")

        integrator.lastDt = dt

        # Set the node properties.
        for n in nodeSet:

            # Set node specific thermal energies
            n.specificThermalEnergy(ScalarField1d("tmp", n, eps1))

            # Set node velocities
            n.velocity(VectorField1d("tmp", n, Vector1d(vx0)))
        del n

        # Use the controller to reinitialize the problem.
        control.reinitializeProblem(restartName,
                                    statsStep = statsStep,
                                    restartStep = restartStep)

        # Iterate the H to convergence.
        control.iterateIdealH()

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
    import NohAnalyticSolution
    rlocal = [pos.x for pos in nodes1.positions().internalValues()]
    r = mpi.reduce(rlocal, mpi.SUM)
    answer = NohAnalyticSolution.NohSolution(1,
                                             r = r,
                                             v0 = vx0,
                                             h0 = h1)

    rwall = 0.05   # Throw away anything with r < rwall to avoid wall heating.
    rmax = 0.15    # Also, throw away anything above rmax
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(nodes1.pressure().internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
    xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
    if mpi.rank == 0:

        # Compute the simulated specific entropy.
        Aprof = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pprof, rhoprof)]

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
                                                 Pn.pnormAverage(1, rmin = rwall, rmax = rmax),
                                                 Pn.pnormAverage(2, rmin = rwall, rmax = rmax),
                                                 Pn.pnormAverage("inf", rmin = rwall, rmax = rmax))
            resultFile.write("%11.6e %11.6e %11.6e " % (Pn.pnormAverage(1, rmin = rwall, rmax = rmax),
                                                        Pn.pnormAverage(2, rmin = rwall, rmax = rmax),
                                                        Pn.pnormAverage("inf", rmin = rwall, rmax = rmax)))
        resultFile.write("\n")
        resultFile.flush()


    #---------------------------------------------------------------------------
    # Throw up the density, velocity, and entropy profiles if we're
    # generating graphics.
    #---------------------------------------------------------------------------
    if nx1 == graphics:
        rho = ScalarFieldList1d()
        rho.appendField(nodes1.massDensity())
        plotFieldList(rho,
                      plot = prho,
                      plotStyle = "points ps 2",
                      lineTitle = "Simulation", # "%i" % nx1,
                      colorDomains = False,
                      colorNodeLists = False)

        vel = VectorFieldList1d()
        vel.appendField(nodes1.velocity())
        plotFieldList(vel,
                      yFunction = "%s.x",
                      plot = pvel,
                      plotStyle = "points ps 2",
                      lineTitle = "Simulation", # "%i" % nx1,
                      colorDomains = False,
                      colorNodeLists = False)

        Afl = ScalarFieldList1d()
        Afield1 = ScalarField1d("Specific entropy", nodes1)
        Afl.appendField(Afield1)
        rho = nodes1.massDensity().internalValues()
        P = nodes1.pressure().internalValues()
        for i in xrange(nodes1.numInternalNodes):
            assert rho[i] > 0.0
            Afield1[i] = P[i]/(rho[i]**gamma)
        plotFieldList(Afl,
                      plot = pA,
                      plotStyle = "points ps 2",
                      lineTitle = "Simulation", # "%i" % nx1,
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

    # Set the xranges so the wall heating effect is evident.
    prho("set xrange [-0.005:0.2]"); prho.refresh()
    pvel("set xrange [-0.005:0.2]"); pvel.refresh()
    pA("set xrange [-0.005:0.2]"); pA.refresh()

    if compatibleEnergy:
        tag = "compatibleEnergy"
    else:
        tag = "standard"
    prho.hardcopy("Noh-planar-1d-%s-rho-profiles.eps" % tag, eps=True, color=False, fontsize=24)
    pvel.hardcopy("Noh-planar-1d-%s-vel-profiles.eps" % tag, eps=True, color=False, fontsize=24)
    pA.hardcopy("Noh-planar-1d-%s-A-profiles.eps" % tag, eps=True, color=False, fontsize=24)

if mpi.rank == 0:
    resultFile.close()
