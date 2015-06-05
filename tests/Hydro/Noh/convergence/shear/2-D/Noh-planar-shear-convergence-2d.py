import os, sys
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import SpheralVisitDump
from GenerateNodeDistribution2d import *
from math import *
mpi, procID, nProcs = loadmpi()
sys.path.append("../../..")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
seed = "lattice" # "xstaggeredLattice"
NodeListConstructor = AsphNodeList2d

nylist = [25, 50, 100 , 200]
nx = 20
nx2, ny2 = 20, 5
nPerh = 3.01

xmin = (0.0, 0.0)
xmax = (0.2, 0.5)
xmin2 = (0.0, 0.5)

rho1 = 1.0
eps1 = 0.0
vshear = 1.0
vy1 = -1.0

gamma = 5.0/3.0
mu = 1.0
#Qconstructor = MonaghanGingoldViscosity2d
Qconstructor = TensorMonaghanGingoldViscosity2d
Cl, Cq = 1.0, 1.0
Qlimiter = True
balsaraCorrection = 0
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
HsmoothMin, HsmoothMax, HratioMin = 1e-10, 0.5, 0.1
HsmoothFraction = 1.0
cfl = 0.5
XSPH = True
epsilonTensile = 0.0
nTensile = 2

neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = 0.5
origin = Vector2d(0.0, 0.0)

goalTime = 0.3
dt = 0.0001
dtMin, dtMax = 1.0e-5, 0.1
dtGrowth = 2.0
useVelocityForDt = False
maxSteps = 30
statsStep = 10
smoothIters = 0
HEvolution = Integrator2d.HEvolutionType.IdealH  #IntegrateH
sumForMassDensity = Integrator2d.MassDensityType.RigorousSum

restartStep = 500
restartDirBase = "Noh-planar-shear-convergence-2d/n=%i"
if NodeListConstructor == AsphNodeList2d:
    restartDirBase += "/asph"
else:
    restartDirBase += "/sph"
restartBaseName = "Noh-planar-shear-2d-%ix%i"

#-------------------------------------------------------------------------------
# A class to store the initial velocity information for the nodes.
# This is useful for post-run diagnostics, and will be restored upon
# restart.
#-------------------------------------------------------------------------------
class RememberNodes(RestartableObject):

    def __init__(self, nodes):
        RestartableObject.__init__(self)
        self.snapshot(nodes)
        return

    def snapshot(self, nodes):
        self.positions = [Vector2d(x) for x in nodes.positions().internalValues()]
        self.velocity = [Vector2d(x) for x in nodes.velocity().internalValues()]
        assert len(self.positions) == nodes.numInternalNodes
        assert len(self.velocity) == nodes.numInternalNodes
        return

    def label(self):
        return "RememberNodes"

    def dumpState(self, file, path):
        file.writeObject(self.positions, path + "/positions")
        file.writeObject(self.velocity, path + "/velocity")
        return

    def restoreState(self, file, path):
        self.positions = file.readObject(path + "/positions")
        self.velocity = file.readObject(path + "/velocity")
        return

#-------------------------------------------------------------------------------
# Set up the basic run time objects.
#-------------------------------------------------------------------------------
title('2-D integrated hydro test -- Shearing Planar Noh convergence test')

pnormFileName = "Noh-planar-convergence-test.txt"
if NodeListConstructor == AsphNodeList2d:
    pnormFileName += "-asph"
else:
    pnormFileName += "-sph"

eos = GammaLawGasMKS2d(gamma, mu)
nodes1 = NodeListConstructor(eos, 0, 0, "Regular Nodes")
nodes2 = SphNodeList2d(eos, 0, 0, "Boundary Nodes")
output('nodes1')
output('nodes2')

rememberNodes1 = RememberNodes(nodes1)
rememberNodes2 = RememberNodes(nodes2)
output('rememberNodes1')
output('rememberNodes2')

# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
output('WT')
output('WTPi')

# Construct the neighbor object and associate it with the node list.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
kernelExtent = WT.kernelExtent()
neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
neighbor2 = NestedGridNeighbor2d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes2.registerNeighbor(neighbor2)
neighborTimer.stop()
neighborTimer.printStatus()

# Create boundary conditions.
xPlane0 = Plane2d(Vector2d(*xmin), Vector2d(1.0, 0.0))
xPlane1 = Plane2d(Vector2d(*xmax), Vector2d(-1.0, 0.0))
yPlane0 = Plane2d(Vector2d(*xmin), Vector2d(0.0, 1.0))
xbc0 = PeriodicBoundary2d(xPlane0, xPlane1)
ybc0 = ReflectingBoundary2d(yPlane0)

# Construct a DataBase to hold our node list
db = DataBase2d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Construct the artificial viscosities for the problem.
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output('q')
output('q.Cl')
output('q.Cq')
output('q.limiter')
output('q.epsilon2')
output('q.negligibleSoundSpeed')
output('q.csMultiplier')
output('q.balsaraShearCorrection')

# Construct the hydro physics object.
hydro = Hydro2d(WT, WTPi, q)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = useVelocityForDt
output('hydro')
output('hydro.kernel()')
output('hydro.PiKernel()')
output('hydro.cfl')
output('hydro.useVelocityMagnitudeForDt')
output('hydro.valid()')

# Construct a predictor corrector integrator, and add the physics packages.
integrator = PredictorCorrectorIntegrator2d(db)
output('integrator')
integrator.appendPhysicsPackage(hydro)
output('integrator.havePhysicsPackage(hydro)')
output('integrator.valid()')
integrator.HsmoothMin = HsmoothMin
integrator.HsmoothMax = HsmoothMax
integrator.HratioMin = HratioMin
output('integrator.HsmoothMin')
output('integrator.HsmoothMax')
output('integrator.HratioMin')
integrator.lastDt = dt
output('integrator.lastDt')
if dtMin:
    integrator.dtMin = dtMin
    output('integrator.dtMin')
if dtMax:
    integrator.dtMax = dtMax
    output('integrator.dtMax')
integrator.dtGrowth = dtGrowth
output('integrator.dtGrowth')
integrator.HEvolution = HEvolution
integrator.sumForMassDensity = sumForMassDensity
if (sumForMassDensity == Integrator2d.MassDensityType.Sum or
    sumForMassDensity == Integrator2d.MassDensityType.RigorousSum):
    integrator.setKernel(WT)
output('integrator.HEvolution')
output('integrator.sumForMassDensity')

# Build the controller.
control = SpheralController(integrator, WT,
                            boundaryConditions = [xbc0, ybc0],
                            statsStep = statsStep,
                            restartStep = restartStep)
output('control')

#===============================================================================
# Create the output file.
#===============================================================================
if mpi.rank == 0:
    resultFile = open(pnormFileName, "w")
    resultFile.write("#             Mass Density                         Pressure                             Velocity (x)                         Velocity (y)                         Specific Thermal Energy              h (x)                                h (y)\n")
    resultFile.write("# N nodes     L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         \n")

#===============================================================================
# Iterate over each resolution, do the simulation, and measure the error.
#===============================================================================
for ny in nylist:
    
    x2 = float(nx)/float(ny)*xmax[1]

    # Explicitly initialize.
    nodes1.numInternalNodes = 0
    nodes1.numGhostNodes = 0
    nodes2.numInternalNodes = 0
    nodes2.numGhostNodes = 0
    integrator.setCurrentTime(0.0)

    # Build the restart file name and directory.
    restartDir = restartDirBase % (ny)
    if mpi.rank == 0:
        os.system("mkdir -p %s" % restartDir)
    restartName = restartDir + "/" + restartBaseName % (nx, ny)

    # Check if there is already available restart info.
    mpi.barrier()
    restoreCycle = findLastRestart(restartName)

    # If we're not restarting, then initialize the problem state.
    if restoreCycle is None:

        #=======================================================================
        # Problem parameters.
        xmax1 = (x2, xmax[1])
        xmax2 = (x2, xmax[1]*(1.0 + float(ny2)/float(ny)))

        #=======================================================================
        # Problem setup.
        # Set node positions for this domain.
        from ParMETISDistributeNodes import distributeNodes2d
        print "Generating node distribution."
        generator1 = GenerateNodeDistribution2d(nx, ny, rho1, seed,
                                                xmin = xmin,
                                                xmax = xmax1,
                                                nNodePerh = nPerh)
        n1 = generator1.globalNumNodes()
        generator2 = GenerateNodeDistribution2d(nx2, ny2, rho1, seed,
                                                xmin = xmin2,
                                                xmax = xmax2,
                                                nNodePerh = nPerh)
        n2 = generator2.globalNumNodes()

        print "Distributing nodes amongst processors."
        nodeInfo = distributeNodes2d([(nodes1, n1, generator1),
                                      (nodes2, n2, generator2)])
        if mpi:
            output('mpi.reduce(nodes1.numInternalNodes, mpi.MIN)')
            output('mpi.reduce(nodes1.numInternalNodes, mpi.MAX)')
            output('mpi.reduce(nodes1.numInternalNodes, mpi.SUM)')
            output('mpi.reduce(nodes2.numInternalNodes, mpi.MIN)')
            output('mpi.reduce(nodes2.numInternalNodes, mpi.MAX)')
            output('mpi.reduce(nodes2.numInternalNodes, mpi.SUM)')
        else:
            output('nodes1.numInternalNodes')
            output('nodes2.numInternalNodes')
        assert len(nodeInfo[nodes1]['globalNodeListID']) == nodes1.numInternalNodes
        assert len(nodeInfo[nodes2]['globalNodeListID']) == nodes2.numInternalNodes

        # Set the node properties.
        for nodes in [nodes1, nodes2]:
            nodes.HsmoothFraction = HsmoothFraction
            nodes.XSPH = XSPH
            nodes.nodesPerSmoothingScale = nPerh
            nodes.epsilonTensile = epsilonTensile
            nodes.nTensile = nTensile
            output('nodes.HsmoothFraction')
            output('nodes.nodesPerSmoothingScale')
            output('nodes.epsilonTensile')
            output('nodes.nTensile')
            output('nodes.XSPH')

            # Set node specific thermal energies
            nodes.setSpecificThermalEnergy(ScalarField2d(nodes, eps1))

            # Set node velocities
            for i in xrange(nodes.numInternalNodes):
                y = nodes.positions()[i].y
                nodes.velocity()[i] = Vector2d(vshear*cos(2.0*pi*y), vy1)

        del nodes

        # Update the RememberNodes with the current node properties.
        rememberNodes1.snapshot(nodes1)
        rememberNodes2.snapshot(nodes2)

        # Update the reflecting boundary condition to recognize the new xmax plane.
        xbc0.setExitPlane(Plane2d(Vector2d(*xmax1), Vector2d(-1.0,0.0)))

        # Use the controller to reinitialize the problem.
        control.reinitializeProblem(restartName,
                                    statsStep = statsStep,
                                    restartStep = restartStep)

        # Smooth the initial conditions.
        control.smoothState(smoothIters)

    else:

        # Load the stored state.
        control.reinitializeProblem(restartName)
        control.loadRestartFile(restoreCycle)

    #==========================================================================
    # Advance to the end time.
    while control.time() < goalTime:
        control.advance(goalTime, maxSteps)
        P = db.fluidPressure
        cs = db.fluidSoundSpeed
        Hi = db.fluidHinverse
        dumper = SpheralVisitDump(db,
                                  "Noh-planar-shear-2d-%i-visit" % ny,
                                  restartDir,
                                  listOfFieldLists = [db.fluidMassDensity,
                                                      db.fluidVelocity,
                                                      db.fluidWeight,
                                                      P,
                                                      cs,
                                                      Hi]
                                  )
        dumper.dump(control.time(), control.totalSteps)

    if restoreCycle < control.totalSteps:
        control.dropRestartFile()

    # Compute the local state we're going to analyze.
    rlocal = [pos.y for pos in nodes1.positions().internalValues()]
    r = mpi.reduce(rlocal, mpi.SUM)

    vel = nodes1.velocity()
    Hinverse = nodes1.Hinverse()
    vx = ScalarField2d(nodes1, 0.0)
    vy = ScalarField2d(nodes1, 0.0)
    hx = ScalarField2d(nodes1, 0.0)
    hy = ScalarField2d(nodes1, 0.0)
    xunit = Vector2d(1.0, 0.0)
    yunit = Vector2d(0.0, 1.0)
    for i in xrange(nodes1.numInternalNodes):
        vx[i] = vel[i].dot(xunit)
        vy[i] = vel[i].dot(yunit)
        hx[i] = (Hinverse[i]*xunit).magnitude()
        hy[i] = (Hinverse[i]*yunit).magnitude()

    # Load the analytic solution.
    sys.path.append("../../..")
    import NohAnalyticSolution
    h0 = nPerh*x2/nx
    answer = NohAnalyticSolution.NohSolution(1,
                                             r = r,
                                             v0 = vy1,
                                             h0 = h0)

    # Compute the error.
    rLpmin = 0.0
    rLpmax = 0.2
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(nodes1.pressure().internalValues(), mpi.SUM)
    vxprof = mpi.reduce(vx.internalValues(), mpi.SUM)
    vyprof = mpi.reduce(vy.internalValues(), mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hxprof = mpi.reduce(hx.internalValues(), mpi.SUM)
    hyprof = mpi.reduce(hy.internalValues(), mpi.SUM)
    xprof = mpi.reduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
    vxans = mpi.reduce([x.x for x in rememberNodes1.velocity], mpi.SUM)
    if mpi.rank == 0:
        assert len(vxans) == len(xprof)
        resultFile.write("%10i " % ny)
        rans, vyans, epsans, rhoans, Pans, hans = answer.solution(control.time())
        yshock = control.time()/3.0
        hxans = []
        hyans = []
        for r in xprof:
            hxans.append(h0)
            if r <= yshock:
                hyans.append(0.25*h0)
                y0 = r - vy1*(3.0*r)
            else:
                hyans.append(h0)
                y0 = r - vy1*control.time()
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
                                                 Pn.pnormAverage(1, rmin = rLpmin, rmax = rLpmax),
                                                 Pn.pnormAverage(2, rmin = rLpmin, rmax = rLpmax),
                                                 Pn.pnormAverage("inf", rmin = rLpmin, rmax = rLpmax))
            resultFile.write("%11.6e %11.6e %11.6e " % (Pn.pnormAverage(1, rmin = rLpmin, rmax = rLpmax),
                                                        Pn.pnormAverage(2, rmin = rLpmin, rmax = rLpmax),
                                                        Pn.pnormAverage("inf", rmin = rLpmin, rmax = rLpmax)))
        resultFile.write("\n")

        # Plot the final state.
        import Gnuplot
        plots = []
        for name, data, ans in [("vx", vxprof, vxans),
                                ("vy", vyprof, vyans),
                                ("rho", rhoprof, rhoans),
                                ("P", Pprof, Pans),
                                ("eps", epsprof, epsans),
                                ("hx", hxprof, hxans),
                                ("hy", hyprof, hyans)]:
            gdata = Gnuplot.Data(xprof, data,
                                 with = "points",
                                 title = name,
                                 inline = True)
            gans = Gnuplot.Data(xprof, ans,
                                with = "points",
                                title = "Solution",
                                inline = True)
            gplot = Gnuplot.Gnuplot()
            gplot.plot(gdata)
            gplot.replot(gans)
            plots.append(gplot)

if mpi.rank == 0:
    resultFile.close()
