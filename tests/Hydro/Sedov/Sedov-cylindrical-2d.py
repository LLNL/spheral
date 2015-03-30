#-------------------------------------------------------------------------------
# The cylindrical Sedov test case (2-D).
#-------------------------------------------------------------------------------
import os, sys, shutil
from Spheral2d import *
from findLastRestart import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from GenerateNodeDistribution2d import *
from CubicNodeGenerator import GenerateSquareNodeDistribution

import mpi

title("2-D integrated hydro test -- planar Sedov problem")
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "constantDTheta",

            thetaFactor = 0.5,
            azimuthalOffsetFraction = 0.0,
            nRadial = 50,
            nTheta = 50,
            rmin = 0.0,
            rmax = 1.0,
            nPerh = 1.51,

            rho0 = 1.0,
            eps0 = 0.0,
            Espike = 1.0,
            smoothSpike = True,
            gamma = 5.0/3.0,
            mu = 1.0,

            Cl = 1.0,
            Cq = 0.75,
            epsilon2 = 1e-2,
            Qlimiter = False,
            balsaraCorrection = False,
            linearInExpansion = False,

            ASPH = False,     # Only for H evolution, not hydro algorithm
            CRKSPH = False,
            Qconstructor = MonaghanGingoldViscosity,
            momentumConserving = True, # For CRKSPH
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            HUpdate = IdealH,
            filter = 0.0,

            HydroConstructor = SPHHydro,
            hmin = 1e-15,
            hmax = 1.0,
            cfl = 0.5,
            useVelocityMagnitudeForDt = True,
            XSPH = False,
            rhomin = 1e-10,

            steps = None,
            goalTime = None,
            goalRadius = 0.8,
            dt = 1e-8,
            dtMin = 1.0e-8,
            dtMax = None,
            dtGrowth = 2.0,
            vizCycle = None,
            vizTime = 0.1,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = IdealH,
            compatibleEnergy = True,
            gradhCorrection = False,

            restoreCycle = None,
            restartStep = 1000,

            graphics = True,
            useVoronoiOutput = False,
            clearDirectories = False,
            dataRoot = "dumps-cylindrical-Sedov",
            outputFile = "None",
            )

assert thetaFactor in (0.5, 1.0, 2.0)
theta = thetaFactor * pi

# Figure out what our goal time should be.
import SedovAnalyticSolution
h0 = 1.0/nRadial*nPerh
answer = SedovAnalyticSolution.SedovSolution(nDim = 2,
                                             gamma = gamma,
                                             rho0 = rho0,
                                             E0 = Espike,
                                             h0 = h0)
if goalTime is None:
    assert not goalRadius is None
    nu1 = 1.0/(answer.nu + 2.0)
    nu2 = 2.0*nu1
    goalTime = (goalRadius*(answer.alpha*rho0/Espike)**nu1)**(1.0/nu2)
vs, r2, v2, rho2, P2 = answer.shockState(goalTime)
print "Predicted shock position %g at goal time %g." % (r2, goalTime)

# Scale the spike energy according to the boundary conditions we're using.
if thetaFactor == 0.5:
    Espike *= 0.25
elif thetaFactor == 1.0:
    Espike *= 0.5

#-------------------------------------------------------------------------------
# Set the hydro choice.
#-------------------------------------------------------------------------------
if CRKSPH:
    Qconstructor = CRKSPHMonaghanGingoldViscosity
    if ASPH:
        HydroConstructor = ACRKSPHHydro
    else:
        HydroConstructor = CRKSPHHydro
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

#-------------------------------------------------------------------------------
# Path names.
#-------------------------------------------------------------------------------
dataDir = os.path.join(dataRoot,
                       str(HydroConstructor).split()[1].split(".")[1][:-2],
                       str(Qconstructor).split()[1].split(".")[-1][:-2],
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "gradhCorrection=%s" % gradhCorrection,
                       "seed=%s" % seed,
                       "nr=%i_nt=%i" % (nRadial, nTheta))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Sedov-cylindrical-2d-%i" % nRadial)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = WT
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Create a NodeList and associated Neighbor object.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           rhoMin = rhomin)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
pos = nodes1.positions()
vel = nodes1.velocity()
mass = nodes1.mass()
eps = nodes1.specificThermalEnergy()
H = nodes1.Hfield()
if restoreCycle is None:
    if seed == "square":
        generator = GenerateSquareNodeDistribution(nRadial,
                                                   nTheta,
                                                   rho0,
                                                   xmin = Vector(0,0),
                                                   xmax = Vector(1,1),
                                                   nNodePerh = nPerh,
                                                   SPH = (not ASPH))
    else:
        generator = GenerateNodeDistribution2d(nRadial, nTheta, rho0, seed,
                                               rmin = rmin,
                                               rmax = rmax,
                                               xmin = Vector(-1,-1),
                                               xmax = Vector(1,1),
                                               theta = theta,
                                               azimuthalOffsetFraction = azimuthalOffsetFraction,
                                               nNodePerh = nPerh,
                                               SPH = (not ASPH))

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
        #from PeanoHilbertDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes1, generator))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set the point source of energy.
    Esum = 0.0
    if smoothSpike:
        Wsum = 0.0
        for nodeID in xrange(nodes1.numInternalNodes):
            Hi = H[nodeID]
            etaij = (Hi*pos[nodeID]).magnitude()
            Wi = WT.kernelValue(etaij, Hi.Determinant())
            Ei = Wi*Espike
            epsi = Ei/mass[nodeID]
            eps[nodeID] = epsi
            Wsum += Wi
        Wsum = mpi.allreduce(Wsum, mpi.SUM)
        assert Wsum > 0.0
        for nodeID in xrange(nodes1.numInternalNodes):
            eps[nodeID] /= Wsum
            Esum += eps[nodeID]*mass[nodeID]
    else:
        i = -1
        rmin = 1e50
        for nodeID in xrange(nodes1.numInternalNodes):
            rij = pos[nodeID].magnitude()
            if rij < rmin:
                i = nodeID
                rmin = rij
        rminglobal = mpi.allreduce(rmin, mpi.MIN)
        if fuzzyEqual(rmin, rminglobal):
            assert i >= 0 and i < nodes1.numInternalNodes
            eps[i] = Espike/mass[i]
            Esum += Espike
    Eglobal = mpi.allreduce(Esum, mpi.SUM)
    print "Initialized a total energy of", Eglobal
    assert fuzzyEqual(Eglobal, Espike)

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
q = Qconstructor(Cl, Cq, linearInExpansion)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if CRKSPH:
    hydro = HydroConstructor(WT, WTPi, q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             momentumConserving = momentumConserving)
else:
    hydro = HydroConstructor(WT, WTPi, q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             densityUpdate = densityUpdate,
                             XSPH = XSPH,
                             HUpdate = HEvolution)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(0.0, 0.0), Vector(1.0, 0.0))
yPlane0 = Plane(Vector(0.0, 0.0), Vector(0.0, 1.0))
xbc0 = ReflectingBoundary(xPlane0)
ybc0 = ReflectingBoundary(yPlane0)

for p in packages:
    if thetaFactor in (0.5, ):
        p.appendBoundary(xbc0)
    if thetaFactor in (0.5, 1.0):
        p.appendBoundary(ybc0)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
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

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
vizMethod = None
if useVoronoiOutput:
    import SpheralVoronoiSiloDump
    vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizMethod = vizMethod,
                            vizBaseName = "Sedov-cylindrical-2d-%ix%i" % (nRadial, nTheta),
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
output("control")

#-------------------------------------------------------------------------------
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------
if steps is None:
    control.advance(goalTime, maxSteps)
    if restoreCycle != control.totalSteps:
        control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
        control.dropRestartFile()
else:
    control.step(steps)

# Output the energy conservation.
print "Energy conservation: ", ((control.conserve.EHistory[-1] -
                                 control.conserve.EHistory[0])/
                                control.conserve.EHistory[0])

#-------------------------------------------------------------------------------
# Generate some error metrics comparing to the analytic solution.
#-------------------------------------------------------------------------------
# Report the error norms.
rmin, rmax = 0.0, 0.95
r = mpi.allreduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)
xprof = mpi.allreduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
yprof = mpi.allreduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
v = mpi.allreduce([x.magnitude() for x in nodes1.velocity().internalValues()], mpi.SUM)
eps = mpi.allreduce(list(nodes1.specificThermalEnergy().internalValues()), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(list(Pf.internalValues()), mpi.SUM)
A = mpi.allreduce([Pi/(rhoi**gamma) for (Pi, rhoi) in zip(Pf.internalValues(), nodes1.massDensity().internalValues())], mpi.SUM)

Hinverse = db.newFluidSymTensorFieldList()
db.fluidHinverse(Hinverse)
hrfl = db.newFluidScalarFieldList()
htfl = db.newFluidScalarFieldList()
for Hfield, hrfield, htfield in zip(Hinverse,
                                    hrfl,
                                    htfl):
    n = Hfield.numElements
    assert hrfield.numElements == n
    assert htfield.numElements == n
    positions = Hfield.nodeList().positions()
    for i in xrange(n):
        runit = positions[i].unitVector()
        tunit = Vector(-(positions[i].y), positions[i].x).unitVector()
        hrfield[i] = (Hfield[i]*runit).magnitude()
        htfield[i] = (Hfield[i]*tunit).magnitude()
hr = mpi.allreduce(list(hrfl[0].internalValues()), mpi.SUM)
ht = mpi.allreduce(list(htfl[0].internalValues()), mpi.SUM)

if mpi.rank == 0:
    from SpheralGnuPlotUtilities import multiSort
    import Pnorm
    multiSort(r, rho, v, eps, P, A, hr, ht)
    rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), r)
    Aans = [Pi/(rhoi**gamma) for (Pi, rhoi) in zip(Pans, rhoans)]
    print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
    for (name, data, ans) in [("Mass Density", rho, rhoans),
                              ("Pressure", P, Pans),
                              ("Velocity", v, vans),
                              ("Thermal E", eps, epsans),
                              ("Entropy", A, Aans),
                              ("hr", hr, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in xrange(len(data))]
        Pn = Pnorm.Pnorm(error, r)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print "\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf)

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile != "None" and mpi.rank == 0:
    outputFile = os.path.join(dataDir, outputFile)
    f = open(outputFile, "w")
    f.write(("# " + 16*"%15s " + "\n") % ("r", "x", "y", "rho", "P", "v", "eps", "A", "hr", "ht",
                                          "rhoans", "Pans", "vans", "epsans", "Aans", "hrans"))
    for (ri, xi, yi, rhoi, Pi, vi, epsi, Ai, hri, hti, 
         rhoansi, Pansi, vansi, epsansi, Aansi, hansi)  in zip(r, xprof, yprof, rho, P, v, eps, A, hr, ht,
                                                               rhoans, Pans, vans, epsans, Aans, hans):
         f.write((16*"%16.12e " + "\n") % (ri, xi, yi, rhoi, Pi, vi, epsi, Ai, hri, hti, 
                                           rhoansi, Pansi, vansi, epsansi, Aansi, hansi))
    f.close()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    plots = [(rPlot, "Sedov-cylindrical-positions.png"),
             (rhoPlot, "Sedov-cylindrical-rho.png"),
             (velPlot, "Sedov-cylindrical-vel.png"),
             (epsPlot, "Sedov-cylindrical-eps.png"),
             (PPlot, "Sedov-cylindrical-P.png"),
             (HPlot, "Sedov-cylindrical-h.png")]

    # Plot the specific entropy.
    AsimData = Gnuplot.Data(xprof, A,
                            with_ = "points",
                            title = "Simulation",
                            inline = True)
    AansData = Gnuplot.Data(xprof, Aans,
                            with_ = "lines",
                            title = "Solution",
                            inline = True)
    
    Aplot = generateNewGnuPlot()
    Aplot.plot(AsimData)
    Aplot.replot(AansData)
    Aplot.title("Specific entropy")
    Aplot.refresh()
    plots.append((Aplot, "Sedov-cylindrical-entropy.png"))

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.hardcopy(os.path.join(dataDir, filename), terminal="png")

