#-------------------------------------------------------------------------------
# The Planar Noh test case run in 1-D.
# This script is designed to automatically measure the convergence rate.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
import os, sys
from Spheral import *
from SpheralTestUtilities import *
from SpheralVisitDump import dumpPhysicsState

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList2d,

            nxlist = [50, 100, 200, 400, 800, 1600],
            nx2 = 10,
            ny = 10,
            rho1 = 1.0,
            eps1 = 0.0,
            x0 = 0.0,
            x1 = 0.5,
            y0 = 0.0,
            nPerh = 2.01,
            seed = "lattice",

            vx0 = -1.0, 
            vshear = 0.0,
            
            gamma = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosity2d,
            #Qconstructor = TensorMonaghanGingoldViscosity2d,
            Cl = 1.5, 
            Cq = 1.5,
            Qlimiter = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.1,
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
            restartStep = 10000,
            
            dataDirBase = "Noh-planar-convergence-2d",
            restartBaseName = "Noh-planar-2d-%i",
            visitBaseName = "Noh-planar-2d-%i",
            )

# Derive the initial value for y1 based on the lowest resolution run.
y1 = y0 + (x1 - x0)/float(nxlist[0])*float(ny)

# Base restart and visit dirs.
dataDirBase = os.path.join(dataDirBase,
                           "%s" % str(NodeListConstructor).split("'")[1],
                           "nperh=%4.2f" % nPerh,
                           "XSPH=%s" % XSPH,
                           "n=%i")
restartDirBase = os.path.join(dataDirBase, "restart")
visitDirBase = os.path.join(dataDirBase, "visit")

#-------------------------------------------------------------------------------
# Create the file we're going to record the error norms in.
#-------------------------------------------------------------------------------
pnormFileName = "Noh-planar-convergence-test-nperh=%4.2f-XSPH=%s.txt" % (nPerh, XSPH)
if mpi.rank == 0:
    resultFile = open(pnormFileName, "w")
    resultFile.write("#             Mass Density                         Pressure                             Velocity x                           Velocity y                           Specific Thermal Energy               hx                                    hy\n")
    resultFile.write("# N nodes     L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf          L1          L2          Linf          L1          L2          Linf          \n")

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
xPlane0 = Plane2d(Vector2d(x0, y0), Vector2d(1.0, 0.0))
yPlane0 = Plane2d(Vector2d(x0, y0), Vector2d(0.0, 1.0))
yPlane1 = Plane2d(Vector2d(x0, y1), Vector2d(0.0, -1.0))
xbc0 = ReflectingBoundary2d(xPlane0)
ybc0 = PeriodicBoundary2d(yPlane0, yPlane1)
for bc in [xbc0, ybc0]:
    hydro.appendBoundary(bc)

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
    visitDir = visitDirBase % nx1
    if mpi.rank == 0:
        if not os.path.exists(restartDir):
            os.makedirs(restartDir)
        if not os.path.exists(visitDir):
            os.makedirs(visitDir)
    mpi.barrier()
    restartName = os.path.join(restartDir, restartBaseName % nx1)
    visitName = visitBaseName % nx1

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
        y1 = y0 + (x1 - x0)/float(nx1)*float(ny)

        # Update the reflecting boundary condition to recognize the new xmax plane.
        yPlane1 = Plane2d(Vector2d(x0, y1), Vector2d(0.0, -1.0))
        ybc0.setExitPlane(Plane2d(yPlane1))

        # Problem setup
        # Set node positions for this domain
        from GenerateNodeDistribution2d import GenerateNodeDistribution2d
        from ParMETISDistributeNodes import distributeNodes2d
        generator1 = GenerateNodeDistribution2d(nx1, ny, rho1, seed,
                                                xmin = (x0, y0),
                                                xmax = (x1, y1),
                                                nNodePerh = nPerh)
        generator2 = GenerateNodeDistribution2d(nx2, ny, rho1, seed,
                                                xmin = (x1, y0),
                                                xmax = (x2, y1),
                                                nNodePerh = nPerh)
        distributeNodes2d((nodes1, generator1),
                          (nodes2, generator2))
        output("mpi.allreduce(nodes1.numNodes, mpi.SUM)")
        output("mpi.allreduce(nodes2.numNodes, mpi.SUM)")

        # Set the node properties.
        for n in nodeSet:

            # Set node specific thermal energies
            n.specificThermalEnergy(ScalarField2d("tmp", n, eps1))

            # Set node velocities
            for i in xrange(n.numNodes):
                xi = n.positions()[i].x
                n.velocity()[i] = Vector2d(vx0, vshear*cos(2.0*pi*xi))
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
        dumpPhysicsState(integrator,
                         visitName,
                         visitDir)

    # Compute the analytic solution.
    sys.path.append("../../..")
    import NohAnalyticSolution
    rlocal = [pos.x for pos in nodes1.positions().internalValues()]
    r = mpi.reduce(rlocal, mpi.SUM)
    answer = NohAnalyticSolution.NohSolution(1,
                                             r = r,
                                             v0 = vx0,
                                             h0 = h1)

    # Compute the error.
    rwall = 0.0   # Throw away anything with r < rwall to avoid wall heating.
    rmax = 0.2    # Also, throw away anything above rmax
    xunit = Vector2d(1.0, 0.0)
    yunit = Vector2d(0.0, 1.0)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(nodes1.pressure().internalValues(), mpi.SUM)
    vxprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    vyprof = mpi.reduce([v.y for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hxprof = mpi.reduce([1.0/(H*xunit).magnitude() for H in nodes1.Hfield().internalValues()], mpi.SUM)
    hyprof = mpi.reduce([1.0/(H*yunit).magnitude() for H in nodes1.Hfield().internalValues()], mpi.SUM)
    xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
    if mpi.rank == 0:
        resultFile.write("%10i " % nx1)
        rans, vxans, epsans, rhoans, Pans, hxans = answer.solution(control.time())
        hyans = [h1]*len(xprof)
        xshock = control.time()/3.0
        vyans = []
        for xi in xprof:
            if xi < xshock:
                xi0 = xi - vx0*(3.0*xi)
            else:
                xi0 = xi - vx0*control.time()
            vyans.append(vshear*cos(2.0*pi*xi0))
        
        import Pnorm
        print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
        for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                                  ("Pressure  ", Pprof, Pans),
                                  ("Velocity x", vxprof, vxans),
                                  ("Velocity y", vyprof, vyans),
                                  ("Thermal E ", epsprof, epsans),
                                  ("hx        ", hxprof, hxans),
                                  ("hy        ", hyprof, hyans)]:
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

if mpi.rank == 0:
    resultFile.close()
