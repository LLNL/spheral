import os, sys
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from GenerateNodeDistribution3d import *
from math import *
mpi, procID, nProcs = loadmpi()
sys.path.append("../../..")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
seed = "lattice" # "xstaggeredLattice"
NodeListConstructor = AsphNodeList3d

nylist = [25, 50, 100] # , 200]
nx = 10
nPerh = 1.51

xmin = (0.0, 0.0, 0.0)
xmax = (0.1, 0.5, 0.1)

rho1 = 1.0
eps1 = 0.0
vshear = 0.0
vy1 = -1.0

gamma = 5.0/3.0
mu = 1.0
Qconstructor = MonaghanGingoldViscosity3d
#Qconstructor = TensorMonaghanGingoldViscosity3d
Cl, Cq = 1.0, 1.0
Qlimiter = False
balsaraCorrection = False
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
HsmoothMin, HsmoothMax, HratioMin = 1e-10, 0.5, 0.1
HsmoothFraction = 0.0
cfl = 0.5
XSPH = True
epsilonTensile = 0.0
nTensile = 4

neighborSearchType = Neighbor3d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = 0.5
origin = Vector3d(0.0, 0.0, 0.0)

goalTime = 0.3
dt = 0.0001
dtMin, dtMax = 1.0e-5, 0.1
dtGrowth = 2.0
useVelocityForDt = False
maxSteps = 30
statsStep = 10
smoothIters = 0
HEvolution = Hydro3d.HEvolutionType.IdealH  #IntegrateH
sumForMassDensity = Hydro3d.MassDensityType.RigorousSumDensity

restartStep = 500
dataDirBase = "Noh-planar-convergence-3d/n=%i"
if NodeListConstructor == AsphNodeList3d:
    dataDirBase += "/asph"
else:
    dataDirBase += "/sph"
restartBaseName = "Noh-planar-3d-%ix%i"

#-------------------------------------------------------------------------------
# Set up the basic run time objects.
#-------------------------------------------------------------------------------
title('3-D integrated hydro test -- 3-D Planar Noh convergence test')

pnormFileName = "Noh-planar-convergence-test.txt"
if NodeListConstructor == AsphNodeList3d:
    pnormFileName += "-asph"
else:
    pnormFileName += "-sph"

eos = GammaLawGasMKS3d(gamma, mu)
nodes1 = NodeListConstructor("Regular Nodes", eos)
nodes1.nodesPerSmoothingScale = nPerh
output('nodes1')

# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
WT = TableKernel3d(BSplineKernel3d(), 1000)
WTPi = TableKernel3d(BSplineKernel3d(), 1000)
output('WT')
output('WTPi')

# Construct the neighbor object and associate it with the node list.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
kernelExtent = WT.kernelExtent()
neighbor1 = NestedGridNeighbor3d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)
neighborTimer.stop()
neighborTimer.printStatus()

# Construct a DataBase to hold our node list
db = DataBase3d()
output('db')
output('db.appendNodeList(nodes1)')
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
hydro = Hydro3d(WT, WTPi, q)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = useVelocityForDt
hydro.HsmoothMin = HsmoothMin
hydro.HsmoothMax = HsmoothMax
hydro.HratioMin = HratioMin
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
output('hydro')
output('hydro.kernel()')
output('hydro.PiKernel()')
output('hydro.cfl')
output('hydro.useVelocityMagnitudeForDt')
output('hydro.valid()')
output('hydro.HEvolution')
output('hydro.sumForMassDensity')
output('hydro.HsmoothMin')
output('hydro.HsmoothMax')
output('hydro.HratioMin')

# Create boundary conditions.
xPlane0 = Plane3d(Vector3d(*xmin), Vector3d(1.0, 0.0, 0.0))
xPlane1 = Plane3d(Vector3d(*xmax), Vector3d(-1.0, 0.0, 0.0))
yPlane0 = Plane3d(Vector3d(*xmin), Vector3d(0.0, 1.0, 0.0))
zPlane0 = Plane3d(Vector3d(*xmin), Vector3d(0.0, 0.0, 1.0))
zPlane1 = Plane3d(Vector3d(*xmax), Vector3d(0.0, 0.0, -1.0))
xbc0 = ReflectingBoundary3d(xPlane0)
xbc1 = ReflectingBoundary3d(xPlane1)
ybc0 = ReflectingBoundary3d(yPlane0)
zbc0 = ReflectingBoundary3d(zPlane0)
zbc1 = ReflectingBoundary3d(zPlane1)

for bc in [xbc0, xbc1, ybc0, zbc0, zbc1]:
    hydro.appendBoundary(bc)

# Construct a predictor corrector integrator, and add the physics packages.
integrator = PredictorCorrectorIntegrator3d(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output('integrator')
output('integrator.havePhysicsPackage(hydro)')
output('integrator.valid()')
output('integrator.lastDt')
output('integrator.dtMin')
output('integrator.dtMax')
output('integrator.dtGrowth')

# Build the controller.
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep)
output('control')

#===============================================================================
# Create the output file.
#===============================================================================
if mpi.rank == 0:
    resultFile = open(pnormFileName, "w")
    resultFile.write("#             Mass Density                         Pressure                             Velocity (x)                         Velocity (y)                         Velocity (z)                         Specific Thermal Energy              h (x)                                h (y)                                h (z)\n")
    resultFile.write("# N nodes     L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         L1          L2          Linf         \n")

#===============================================================================
# Iterate over each resolution, do the simulation, and measure the error.
#===============================================================================
for ny in nylist:
    
    x2 = float(nx)/float(ny)*xmax[1]

    # Explicitly initialize.
    nodes1.numInternalNodes = 0
    nodes1.numGhostNodes = 0
    integrator.setCurrentTime(0.0)

    # Build the file names.
    dataDir = dataDirBase % (ny)
    restartDir = (dataDirBase % ny) + "/restarts"
    visitDir = (dataDirBase % ny) + "/visit"
    restartName = restartDir + "/" + restartBaseName % (nx, ny)

    # Check if the necessary output directories exist.  If not, create them.
    if mpi.rank == 0:
        if not os.path.exists(restartDir):
            os.makedirs(restartDir)
        if not os.path.exists(visitDir):
            os.makedirs(visitDir)
    mpi.barrier()

    # Check if there is already available restart info.
    restoreCycle = findLastRestart(restartName)

    # If we're not restarting, then initialize the problem state.
    if restoreCycle is None:

        #=======================================================================
        # Problem parameters.
        xmax1 = (x2, xmax[1], x2)

        #=======================================================================
        # Problem setup.
        # Set node positions for this domain.
        from ParMETISDistributeNodes import distributeNodes3d
        print "Generating node distribution."
        generator1 = GenerateNodeDistribution3d(nx, ny, nx, rho1, seed,
                                                xmin = xmin,
                                                xmax = xmax1,
                                                nNodePerh = nPerh)
        n1 = generator1.globalNumNodes()

        print "Distributing nodes amongst processors."
        distributeNodes3d([(nodes1, generator1)])
        if mpi:
            output('mpi.reduce(nodes1.numInternalNodes, mpi.MIN)')
            output('mpi.reduce(nodes1.numInternalNodes, mpi.MAX)')
            output('mpi.reduce(nodes1.numInternalNodes, mpi.SUM)')
        else:
            output('nodes1.numInternalNodes')

        # Set the node properties.
        for nodes in [nodes1]:
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
            nodes.setSpecificThermalEnergy(ScalarField3d("tmp", nodes, eps1))

            # Set node velocities
            for i in xrange(nodes.numInternalNodes):
                y = nodes.positions()[i].y
                nodes.velocity()[i] = Vector3d(vshear*cos(2.0*pi*y), vy1, 0.0)

        # Update the reflecting boundary condition to recognize the new xmax plane.
        xbc1.setEnterPlane(Plane3d(Vector3d(*xmax1), Vector3d(-1.0, 0.0, 0.0)))
        xbc1.setExitPlane(Plane3d(Vector3d(*xmax1), Vector3d(-1.0, 0.0, 0.0)))
        zbc1.setEnterPlane(Plane3d(Vector3d(*xmax1), Vector3d(0.0, 0.0, -1.0)))
        zbc1.setExitPlane(Plane3d(Vector3d(*xmax1), Vector3d(0.0, 0.0, -1.0)))

        print "Simulation boundaries:  ", xmin, xmax1

        # Update the controller's idea of the restart name.
        control.setRestartBaseName(restartName)

        # Smooth the initial conditions.
        control.smoothState(smoothIters)

    else:

        # Load the stored state.
        control.setRestartBaseName(restartName)
        control.loadRestartFile(restoreCycle)

    #==========================================================================
    # Advance to the end time.
    while control.time() < goalTime:
        control.advance(goalTime, maxSteps)
        dumpPhysicsState(integrator,
                         "Noh-planar-3d-%i-visit" % ny,
                         visitDir)
        control.dropRestartFile()

    if restoreCycle < control.totalSteps:
        control.dropRestartFile()

    # Compute the local state we're going to analyze.
    rlocal = [pos.y for pos in nodes1.positions().internalValues()]
    r = mpi.reduce(rlocal, mpi.SUM)

    vel = nodes1.velocity()
    Hinverse = nodes1.Hinverse()
    vx = ScalarField3d("vx", nodes1, 0.0)
    vy = ScalarField3d("vy", nodes1, 0.0)
    vz = ScalarField3d("vz", nodes1, 0.0)
    hx = ScalarField3d("hx", nodes1, 0.0)
    hy = ScalarField3d("hy", nodes1, 0.0)
    hz = ScalarField3d("hz", nodes1, 0.0)
    xunit = Vector3d(1.0, 0.0, 0.0)
    yunit = Vector3d(0.0, 1.0, 0.0)
    zunit = Vector3d(0.0, 0.0, 1.0)
    for i in xrange(nodes1.numInternalNodes):
        vx[i] = vel[i].dot(xunit)
        vy[i] = vel[i].dot(yunit)
        vz[i] = vel[i].dot(zunit)
        hx[i] = (Hinverse[i]*xunit).magnitude()
        hy[i] = (Hinverse[i]*yunit).magnitude()
        hz[i] = (Hinverse[i]*zunit).magnitude()

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
    vzprof = mpi.reduce(vz.internalValues(), mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hxprof = mpi.reduce(hx.internalValues(), mpi.SUM)
    hyprof = mpi.reduce(hy.internalValues(), mpi.SUM)
    hzprof = mpi.reduce(hz.internalValues(), mpi.SUM)
    xprof = mpi.reduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
    if mpi.rank == 0:
        resultFile.write("%10i " % ny)
        rans, vyans, epsans, rhoans, Pans, hans = answer.solution(control.time())
        yshock = control.time()/3.0
        hxans = []
        hyans = []
        hzans = []
        vxans = []
        vzans = []
        for r in xprof:
            hxans.append(h0)
            if r <= yshock:
                hyans.append(0.25*h0)
                y0 = r - vy1*(3.0*r)
            else:
                hyans.append(h0)
                y0 = r - vy1*control.time()
            hzans.append(h0)
            vxans.append(vshear*cos(2.0*pi*y0))
            vzans.append(0.0)
        assert len(vxans) == len(xprof)
        import Pnorm
        print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
        for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                                  ("Pressure  ", Pprof, Pans),
                                  ("Velocity x", vxprof, vxans),
                                  ("Velocity y", vyprof, vyans),
                                  ("Velocity z", vzprof, vzans),
                                  ("Thermal E ", epsprof, epsans),
                                  ("hx        ", hxprof, hxans),
                                  ("hy        ", hyprof, hyans),
                                  ("hz        ", hzprof, hzans)]:
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
        resultFile.flush()

##         # Plot the final state.
##         try:
##             import Gnuplot
##             plots = []
##             for name, data, ans in [("vx", vxprof, vxans),
##                                     ("vy", vyprof, vyans),
##                                     ("vz", vzprof, vzans),
##                                     ("rho", rhoprof, rhoans),
##                                     ("P", Pprof, Pans),
##                                     ("eps", epsprof, epsans),
##                                     ("hx", hxprof, hxans),
##                                     ("hy", hyprof, hyans),
##                                     ("hz", hzprof, hzans)]:
##                 gdata = Gnuplot.Data(xprof, data,
##                                      with = "points",
##                                      title = name,
##                                      inline = True)
##                 gans = Gnuplot.Data(xprof, ans,
##                                     with = "points",
##                                     title = "Solution",
##                                     inline = True)
##                 gplot = Gnuplot.Gnuplot()
##                 gplot.plot(gdata)
##                 gplot.replot(gans)
##                 plots.append(gplot)
##         except:
##             print "No graphics."
##             pass

if mpi.rank == 0:
    resultFile.close()
