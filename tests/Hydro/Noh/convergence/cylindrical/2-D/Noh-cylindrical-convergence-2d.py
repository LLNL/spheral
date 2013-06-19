#-------------------------------------------------------------------------------
# The Cylindrical Noh test case run in 2-D.
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
commandLine(NodeListConstructor = SphNodeList2d,

            nxlist = [25, 50, 100, 200, 400], # 800], # , 1600], #, 3200, 6400],
            rho1 = 1.0,
            eps1 = 0.0,
            x1 = 0.5,
            nPerh = 2.01,
            vx0 = -1.0, 

            seed = "constantDTheta",
            theta = 0.5*pi,

            gamma = 5.0/3.0,
            mu = 1.0,
            #Qconstructor = MonaghanGingoldViscosity2d,
            Qconstructor = TensorMonaghanGingoldViscosity2d,
            Cl = 1.5, 
            Cq = 1.5,
            Qlimiter = True,
            epsilon2 = 1e-2,
            hmin = 1e-10,
            hmax = 1.0,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            hourglassAccelerationFactor = 0.0,

            neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector2d(0.0, 0.0),

            goalTime = 0.3,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 100,
            smoothIters = 0,
            HEvolution = Hydro2d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro2d.MassDensityType.RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,

            restoreCycle = None,
            restartStep = 100,
            
            dataDirBase = "Noh-cylindrical-convergence-2d",
            restartBaseName = "Noh-cylindrical-2d-%i",

            graphics = False,
            )

if graphics:
    nxlist = [graphics]

x0 = 0.0
xmin = (0.0, 0.0)
xmax = (x1, x1)

# Base restart dir.
dataDirBase = os.path.join(dataDirBase,
                           "%s" % str(NodeListConstructor).split("'")[1],
                           "nperh=%4.2f" % nPerh,
                           "XSPH=%s" % XSPH,
                           "compatibleEnergy=%s" % compatibleEnergy,
                           "n=%i")
restartDirBase = os.path.join(dataDirBase, "restart")

#-------------------------------------------------------------------------------
# Helpful function to average radial profiles.
#-------------------------------------------------------------------------------
def collapseRadialProfiles(*listOfProfs):
    nlists = len(listOfProfs)
    n = len(listOfProfs[0])
    for l in listOfProfs:
        assert len(l) == n

    ll = zip(*listOfProfs)
    ll.sort()
    
    result = []
    for i in xrange(nlists):
        result.append([ll[0][i]])
    assert len(result) == nlists

    tol = 1e-2
    for i in xrange(n):
        r0 = result[0][-1]
        ri = ll[i][0]
        assert ri >= r0
        if abs(ri - r0)/ri > tol:
            for j in xrange(nlists):
                result[j].append(ll[i][j])

    return tuple(result)

#-------------------------------------------------------------------------------
# Create the file we're going to record the error norms in.
#-------------------------------------------------------------------------------
pnormFileName = "Noh-cylindrical-convergence-test-nperh=%4.2f-XSPH=%s-compatibleEnergy=%s.txt" % (nPerh, XSPH, compatibleEnergy)
if mpi.rank == 0:
    resultFile = open(pnormFileName, "w")
    resultFile.write("#             Mass Density                         Pressure                             Velocity                             Specific Thermal Energy              Specific Entropy                     h\n")
    resultFile.write("# N nodes     L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         \n")

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS2d(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes = NodeListConstructor("nodes", eos, WT, WTPi)
nodeSet = [nodes]
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
    neighbor = NestedGridNeighbor2d(n,
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
db = DataBase2d()
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
hydro = Hydro2d(WT, WTPi, q, compatibleEnergy)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
plane0 = Plane2d(Vector2d(0.0, 0.0), Vector2d(1.0, 0.0))
plane1 = Plane2d(Vector2d(0.0, 0.0), Vector2d(0.0, 1.0))
bc0 = ReflectingBoundary2d(plane0)
bc1 = ReflectingBoundary2d(plane1)
hydro.appendBoundary(bc0)
hydro.appendBoundary(bc1)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator2d(db)
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
    #prho.xlabel("r")
    #prho.ylabel("{/Symbol r}")
    prho("set yrange [0:18]")

    pvel = generateNewGnuPlot()
    #pvel.xlabel("r")
    #pvel.ylabel("v_r")

    pA = generateNewGnuPlot()
    #pA.xlabel("r")
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

        # Problem setup
        # Set node positions for this domain
        from ParMETISDistributeNodes import distributeNodes2d
        from GenerateNodeDistribution2d import *
        generator1 = GenerateNodeDistribution2d(nx1, nx1,
                                                rho1,
                                                seed,
                                                xmin = xmin,
                                                xmax = xmax,
                                                rmin = 0.0,
                                                rmax = x1,
                                                theta = theta,
                                                nNodePerh = nPerh)
        distributeNodes2d((nodes, generator1))
        output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
        output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
        output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

        integrator.lastDt = dt

        # Set the node properties.
        for n in nodeSet:

            # Set node specific thermal energies
            n.specificThermalEnergy(ScalarField2d("tmp", n, eps1))

            # Set node velocities
            for i in xrange(n.numNodes):
                runit = nodes.positions()[i].unitVector()
                n.velocity()[i] = runit*vx0
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
    rlocal = [pos.magnitude() for pos in nodes.positions().internalValues()]
    r = mpi.reduce(rlocal, mpi.SUM)
    answer = NohAnalyticSolution.NohSolution(2,
                                             r = r,
                                             v0 = vx0,
                                             h0 = h1)

    rwall = 0.0    # Throw away anything with r < rwall to avoid wall heating.
    rmax = 0.15    # Also, throw away anything above rmax
    vel = nodes.velocity()
    Hinverse = nodes.Hinverse()
    vr = ScalarField2d("vr", nodes, 0.0)
    hr = ScalarField2d("hr", nodes, 0.0)
    ht = ScalarField2d("ht", nodes, 0.0)
    for i in xrange(nodes.numInternalNodes):
        ri = nodes.positions()[i]
        runit = ri.unitVector()
        tunit = Vector2d(-(ri.y), ri.x).unitVector()
        vr[i] = vel[i].dot(runit)
        hr[i] = (Hinverse[i]*runit).magnitude()
        ht[i] = (Hinverse[i]*tunit).magnitude()
    rhoprof = mpi.reduce(nodes.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(nodes.pressure().internalValues(), mpi.SUM)
    vprof = mpi.reduce(vr.internalValues(), mpi.SUM)
    epsprof = mpi.reduce(nodes.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce(hr.internalValues(), mpi.SUM)
    htprof = mpi.reduce(ht.internalValues(), mpi.SUM)
    xprof = mpi.reduce([x.magnitude() for x in nodes.positions().internalValues()], mpi.SUM)
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
        if mpi.rank == 0:

            # First collapse the profile to a single value per radial position.
            x, rho, vr, Ar = collapseRadialProfiles(xprof, rhoprof, vprof, Aprof)
            drho = Gnuplot.Data(x, rho,
                                with = "points ps 2",
                                title = "Simulation") # "%i" % nx1)
            prho.replot(drho)

            dvel = Gnuplot.Data(x, vr,
                                with = "points ps 2",
                                title = "Simulation") # "%i" % nx1)
            pvel.replot(dvel)

            dA = Gnuplot.Data(x, Ar,
                              with = "points ps 2",
                              title = "Simulation") # "%i" % nx1)
            pA.replot(dA)

    mpi.barrier()

#-------------------------------------------------------------------------------
# Plot the answers.
#-------------------------------------------------------------------------------
if (mpi.rank == 0) and graphics:
    x, rho, vr, Ar = collapseRadialProfiles(xprof, rhoans, vans, Aans)
    drho = Gnuplot.Data(x, rho,
                        with = "lines lw 2",
                        title = "analytic")
    prho.replot(drho)

    dvel = Gnuplot.Data(x, vr,
                        with = "lines lw 2",
                        title = "analytic")
    pvel.replot(dvel)

    dA = Gnuplot.Data(x, Ar,
                      with = "lines lw 2",
                      title = "analytic")
    pA.replot(dA)

    prho("set xrange [-0.005:0.2]"); prho.refresh()
    pvel("set xrange [-0.005:0.2]"); pvel.refresh()
    pA("set xrange [-0.005:0.2]; set yrange [-0.005:0.1]"); pA.refresh()

    if compatibleEnergy:
        tag = "compatibleEnergy"
    else:
        tag = "standard"
    prho.hardcopy("Noh-cylindrical-2d-%s-rho-profiles.eps" % tag, eps=True, color=False, fontsize=24)
    pvel.hardcopy("Noh-cylindrical-2d-%s-vel-profiles.eps" % tag, eps=True, color=False, fontsize=24)
    pA.hardcopy("Noh-cylindrical-2d-%s-A-profiles.eps" % tag, eps=True, color=False, fontsize=24)

if mpi.rank == 0:
    resultFile.close()
