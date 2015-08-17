import os, sys
sys.path.append("..")
from Spheral import *
from SpheralTestUtilities import *
from SodAnalyticSolution import *
from SpheralGnuPlotUtilities import *
import Gnuplot

title('1-D integrated hydro test -- planar Sod problem')

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nxlist = [50, 100, 200, 400, 800, 1600, 3200, 6400],
            #nxlist = [100, 200, 400, 800, 1600, 3200, 6400], 

            rho1 = 1.0,
            rho2 = 0.25,
            P1 = 1.0,
            P2 = 0.1795,

            x0 = -0.5,
            x1 = 0.0,
            x2 = 0.5,

            smoothDiscontinuity = False,

            nPerh = 1.01,

            gamma = 5.0/3.0,
            mu = 1.0,
            Qconstructor = MonaghanGingoldViscosity1d,
            Cl = 1.0,
            Cq = 0.75,
            Qlimiter = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            energyMultiplier = 1.0,
            hmin = 1e-10,
            hmax = 1.0,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,

            neighborSearchType = Neighbor1d.NeighborSearchType.GatherScatter,
            numGridLevels = 10,
            topGridCellSize = 2.0,
            origin = Vector1d(0.0),

            IntegratorConstructor = PredictorCorrectorIntegrator1d, # SynchronousRK2Integrator1d, # 
            goalTime = 0.15,
            dt = 1e-4,
            dtMin = 1.0e-5,
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = Hydro1d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro1d.MassDensityType.RigorousSumDensity, # VolumeScaledDensity
            compatibleEnergy = True,

            restoreCycle = None,
            restartStep = 200,
            dataDirBase = "Sod-planar-convergence-1d",
            restartBaseName = "Sod-planar-1d-restart",

            graphics = False,
            )

if graphics:
    nxlist = [graphics]

# Base restart dir.
dataDirBase = os.path.join(dataDirBase,
                           "nperh=%4.2f" % nPerh,
                           "XSPH=%s" % XSPH,
                           "compatibleEnergy=%s" % compatibleEnergy,
                           "n=%i")
restartDirBase = os.path.join(dataDirBase, "restart")

#-------------------------------------------------------------------------------
# Create the file we're going to record the error norms in.
#-------------------------------------------------------------------------------
pnormFileName = "Sod-planar-convergence-test-nperh=%4.2f-XSPH=%s-compatibleEnergy=%s.txt" % (nPerh, XSPH, compatibleEnergy)
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
WTPi = WT
output('WT')
output('WTPi')
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes1 = SphNodeList1d("nodes1", eos, WT, WTPi)
nodes2 = SphNodeList1d("nodes2", eos, WT, WTPi)

#-------------------------------------------------------------------------------
# Set the XSPH and tensile corrections for the NodeList
#-------------------------------------------------------------------------------
for nodes in [nodes1, nodes2]:
    nodes.nodesPerSmoothingScale = nPerh
    nodes.XSPH = XSPH
    nodes.epsilonTensile = epsilonTensile
    nodes.nTensile = nTensile
    nodes.hmin = hmin
    nodes.hmax = hmax
    output("nodes.name()")
    output("  nodes.nodesPerSmoothingScale")
    output("  nodes.XSPH")
    output("  nodes.epsilonTensile")
    output("  nodes.nTensile")
    output("  nodes.hmin")
    output("  nodes.hmax")
del nodes

#-------------------------------------------------------------------------------
# Construct the neighbor objects.
#-------------------------------------------------------------------------------
neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
neighbor2 = NestedGridNeighbor1d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
nodes2.registerNeighbor(neighbor2)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase1d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro1d(WT, WTPi, q, compatibleEnergy)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
output('hydro')
output('hydro.compatibleEnergyEvolution')
output('hydro.cfl')
output('hydro.HEvolution')
output('hydro.sumForMassDensity')
output('hydro.HsmoothMin')
output('hydro.HsmoothMax')
output('hydro.kernel')
output('hydro.PiKernel')
output('hydro.valid()')

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane1d(Vector1d(x0), Vector1d( 1.0))
xPlane1 = Plane1d(Vector1d(x2), Vector1d(-1.0))
xbc0 = ReflectingBoundary1d(xPlane0)
xbc1 = ReflectingBoundary1d(xPlane1)
hydro.appendBoundary(xbc0)
hydro.appendBoundary(xbc1)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtGrowth = dtGrowth
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.rigorousBoundaries = rigorousBoundaries
output('integrator')
output('integrator.havePhysicsPackage(hydro)')
output('integrator.valid()')
output('integrator.lastDt')
output('integrator.dtMin')
output('integrator.dtMax')
output('integrator.rigorousBoundaries')

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output('control')

# Smooth the initial conditions.
if restoreCycle:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH(hydro)

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
    pvel("set key graph 0.4, graph 0.85")
    #pvel("set key top left")
    
    pA = generateNewGnuPlot()
    #pA.xlabel("x")
    #pA.ylabel("A")
    #pA("set key graph 0.2, graph 0.9")
    pA("set key top left")

#-------------------------------------------------------------------------------
# Iterate over each resolution, do the simulation, and measure the error.
#-------------------------------------------------------------------------------
for nx1 in nxlist:

    #---------------------------------------------------------------------------
    # Explicitly initialize.
    #---------------------------------------------------------------------------
    for n in (nodes1, nodes2):
        n.numInternalNodes = 0
        n.numGhostNodes = 0
    del n
    dx = (x1 - x0)/nx1
    h1 = nPerh*dx
    nx2 = nx1 # int(rho2/rho1*nx1 + 0.5)

    #---------------------------------------------------------------------------
    # Build the restart file name and directory.
    #---------------------------------------------------------------------------
    restartDir = restartDirBase % nx1
    if mpi.rank == 0:
        if not os.path.exists(restartDir):
            os.makedirs(restartDir)
    mpi.barrier()
    restartName = os.path.join(restartDir, restartBaseName)

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
        from DistributeNodes import distributeNodesInRange1d
        distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1)),
                                  (nodes2, nx2, rho2, (x1, x2))])
        output("mpi.allreduce(nodes1.numNodes, mpi.SUM)")
        output("mpi.allreduce(nodes2.numNodes, mpi.SUM)")

        # Set node specific thermal energies
        eps1 = P1/((gamma - 1.0)*rho1)
        eps2 = P2/((gamma - 1.0)*rho2)
        nodes1.specificThermalEnergy(ScalarField1d("thermal energy", nodes1, eps1))
        nodes2.specificThermalEnergy(ScalarField1d("thermal energy", nodes2, eps2))

        # If requested, smooth the state at the discontinuity.
        if smoothDiscontinuity:
            W0 = WT.kernelValue(0.0, 1.0)
            m1 = (x1 - x0)*rho1/nx1
            m2 = (x2 - x1)*rho2/nx2
            dx1 = (x1 - x0)/nx1
            dx2 = (x2 - x1)/nx2
            h1 = 2.0*dx1 * nPerh
            h2 = 2.0*dx2 * nPerh
            H1 = SymTensor1d(1.0/(nPerh * dx1))
            H2 = SymTensor1d(1.0/(nPerh * dx2))
            A1 = P1/(rho1**gamma)
            A2 = P2/(rho2**gamma)

            for i in xrange(nodes1.numInternalNodes):
                fi = WT.kernelValue(abs(nodes1.positions()[i].x - x1)/h1, 1.0) / W0
                fi = 1.0 - min(1.0, abs(nodes1.positions()[i].x - x1)/h1)
                assert fi >= 0.0 and fi <= 1.0
                nodes1.mass()[i] = (1.0 - fi)*m1 + 0.5*fi*(m1 + m2)
            for i in xrange(nodes2.numInternalNodes):
                fi = WT.kernelValue(abs(nodes2.positions()[i].x - x1)/h2, 1.0) / W0
                fi = 1.0 - min(1.0, abs(nodes2.positions()[i].x - x1)/h2)
                assert fi >= 0.0 and fi <= 1.0
                nodes2.mass()[i] = (1.0 - fi)*m2 + 0.5*fi*(m1 + m2)

        # Use the controller to reinitialize the problem.
        control.reinitializeProblem(restartName,
                                    statsStep = statsStep,
                                    restartStep = restartStep)

        # Iterate the H to convergence.
        control.iterateIdealH()
        db.updateFluidMassDensity(WT)

##        # If requested, smooth the state at the discontinuity.
##        if smoothDiscontinuity:
##            W0 = WT.kernelValue(0.0, 1.0)
##            m1 = (x1 - x0)*rho1/nx1
##            m2 = (x2 - x1)*rho2/nx2
##            h1 = 2.0*(x1 - x0)/nx1 * nPerh
##            h2 = 2.0*(x2 - x1)/nx2 * nPerh
##            for i in xrange(nodes1.numInternalNodes):
##                fi = WT.kernelValue(abs(nodes1.positions()[i].x - x1)/h1, 1.0) / W0
##                assert fi >= 0.0 and fi <= 1.0
##                nodes1.mass()[i] = (1.0 - fi)*m1 + 0.5*fi*(m1 + m2)
##                nodes1.specificThermalEnergy()[i] = (1.0 - fi)*eps1 + 0.5*fi*(eps1 + eps2)
##            for i in xrange(nodes2.numInternalNodes):
##                fi = WT.kernelValue(abs(nodes2.positions()[i].x - x1)/h2, 1.0) / W0
##                assert fi >= 0.0 and fi <= 1.0
##                nodes2.mass()[i] = (1.0 - fi)*m2 + 0.5*fi*(m1 + m2)
##                nodes2.specificThermalEnergy()[i] = (1.0 - fi)*eps2 + 0.5*fi*(eps1 + eps2)
##            db.updateFluidMassDensity(WT)
##            nodes1.updateWeight()
##            nodes2.updateWeight()

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
    dx1 = (x1 - x0)/nx1
    dx2 = (x2 - x1)/nx2
    h1 = 1.0/(nPerh*dx1)
    h2 = 1.0/(nPerh*dx2)
    answer = SodSolution(nPoints = nx1 + nx2,
                         gamma = gamma,
                         rho1 = rho1,
                         P1 = P1,
                         rho2 = rho2,
                         P2 = P2,
                         x0 = x0,
                         x1 = x1,
                         x2 = x2,
                         h1 = 1.0/h1,
                         h2 = 1.0/h2)

    def createList(x):
        xx = x
        if xx == []:
            xx = [-1e50,]
        tmp = mpi.allreduce(xx, mpi.SUM)
        result = [y for y in tmp if y != -1e50]
        return result

    # Compute the error.
    rhoprof = (createList(nodes1.massDensity().internalValues()) +
               createList(nodes2.massDensity().internalValues()))
    Pprof = (createList(nodes1.pressure().internalValues()) +
             createList(nodes2.pressure().internalValues()))
    vprof = (createList([v.x for v in nodes1.velocity().internalValues()]) +
             createList([v.x for v in nodes2.velocity().internalValues()]))
    epsprof = (createList(nodes1.specificThermalEnergy().internalValues()) +
               createList(nodes2.specificThermalEnergy().internalValues()))
    hprof = (createList([1.0/H.xx for H in nodes1.Hfield().internalValues()]) +
             createList([1.0/H.xx for H in nodes2.Hfield().internalValues()]))
    xprof = (createList([x.x for x in nodes1.positions().internalValues()]) +
             createList([x.x for x in nodes2.positions().internalValues()]))

    # Compute the simulated specific entropy.
    Aprof = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pprof, rhoprof)]

    xmin = x0 # answer.x1 + answer.vs*control.time() + 2*dx1
    xmax = x2

    if mpi.rank == 0:
        resultFile.write("%10i " % (nx1 + nx2))
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
                                                 Pn.pnormAverage(1, rmin = xmin, rmax = xmax),
                                                 Pn.pnormAverage(2, rmin = xmin, rmax = xmax),
                                                 Pn.pnormAverage("inf", rmin = xmin, rmax = xmax))
            resultFile.write("%11.6e %11.6e %11.6e " % (Pn.pnormAverage(1, rmin = xmin, rmax = xmax),
                                                        Pn.pnormAverage(2, rmin = xmin, rmax = xmax),
                                                        Pn.pnormAverage("inf", rmin = xmin, rmax = xmax)))
        resultFile.write("\n")
        resultFile.flush()

    #---------------------------------------------------------------------------
    # Throw up the density, velocity, and entropy profiles if we're
    # generating graphics.
    #---------------------------------------------------------------------------
    if nx1 == graphics:
        plotFieldList(db.fluidMassDensity,
                      plot = prho,
                      plotStyle = "points ps 2",
                      lineTitle = "Simulation", #  "(%i nodes)" % (nx1 + nx2),
                      colorDomains = False,
                      colorNodeLists = False)

        plotFieldList(db.fluidVelocity,
                      yFunction = "%s.x",
                      plot = pvel,
                      plotStyle = "points ps 2",
                      lineTitle = "Simulation", # "(%i nodes)" % (nx1 + nx2),
                      colorDomains = False,
                      colorNodeLists = False)

        Afl = ScalarFieldList1d()
        Afield1 = ScalarField1d("Specific entropy", nodes1)
        Afield2 = ScalarField1d("Specific entropy", nodes2)
        Afl.appendField(Afield1)
        Afl.appendField(Afield2)
        for (nodes, Af) in [(nodes1, Afield1),
                            (nodes2, Afield2)]:
            rho = nodes.massDensity().internalValues()
            P = nodes.pressure().internalValues()
            for i in xrange(nodes.numInternalNodes):
                assert rho[i] > 0.0
                Af[i] = P[i]/(rho[i]**gamma)
        del nodes, Af
        plotFieldList(Afl,
                      plot = pA,
                      plotStyle = "points ps 2",
                      lineTitle = "Simulation", # "(%i nodes)" % (nx1 + nx2),
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
    prho.hardcopy("Sod-planar-1d-%s-rho-profiles.eps" % tag, eps=True, color=False, fontsize=24)
    pvel.hardcopy("Sod-planar-1d-%s-vel-profiles.eps" % tag, eps=True, color=False, fontsize=24)
    pA.hardcopy("Sod-planar-1d-%s-A-profiles.eps" % tag, eps=True, color=False, fontsize=24)

if mpi.rank == 0:
    resultFile.close()
