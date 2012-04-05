#ATS:t0 = test(SELF,       "--steps=40 --restartStep 20  --graphics False --vizTime None --clearDirectories True", np=4, label="Cylindrical Noh restart test INITIAL RUN")
#ATS:t1 = testif(t0, SELF, "--steps 20 --restartStep 100 --graphics False --vizTime None --restoreCycle 20 --checkRestart True", np=4, label="Cylindrical Noh restart test RESTARTED CHECK")
#ATS:t2 = testif(t1, SELF, "--steps 40 --restartStep 20  --graphics False --vizTime None --clearDirectories True --IntegratorConstructor SynchronousRK2Integrator", np=4, label="Cylindrical Noh restart test INITIAL RUN (SynchronousRK2)")
#ATS:t3 = testif(t2, SELF, "--steps 20 --restartStep 100 --graphics False --vizTime None --restoreCycle 20 --checkRestart True --IntegratorConstructor SynchronousRK2Integrator", np=4, label="Cylindrical Noh restart test RESTARTED CHECK (SynchronousRK2)")
#ATS:t10 = test(SELF,        "--steps 40 --graphics False --vizTime None --dataDir 'dumps-cylindrical-reproducing' --clearDirectories True  --domainIndependent True --outputFile 'Noh-cylindrical-1proc-reproducing.txt'", np=1, label="Cylindrical Noh domain independence test SERIAL RUN")
#ATS:t11 = testif(t10, SELF, "--steps 40 --graphics False --vizTime None --dataDir 'dumps-cylindrical-reproducing' --clearDirectories False --domainIndependent True --outputFile 'Noh-cylindrical-4proc-reproducing.txt' --comparisonFile 'Noh-cylindrical-1proc-reproducing.txt'", np=4, label="Cylindrical Noh domain independence test 4 PROC RUN")

#-------------------------------------------------------------------------------
# The Cylindrical Noh test case run in 2-D.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CubicNodeGenerator import GenerateSquareNodeDistribution
from CentroidalVoronoiRelaxation import *

import mpi
import DistributeNodes

title("2-D integrated hydro test -- cylindrical Noh problem")

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
            nPerh = 2.01,

            vr0 = -1.0, 

            gamma = 5.0/3.0,
            mu = 1.0,

            HydroConstructor = ASPHHydro,
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            Cl = 1.0, 
            Cq = 0.75,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = None,
            vizCycle = None,
            vizTime = 0.1,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,

            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = False,

            clearDirectories = False,
            restoreCycle = None,
            restartStep = 1000,
            checkRestart = False,
            dataDir = "dumps-cylindrical",
            outputFile = "None",
            comparisonFile = "None",

            graphics = True,
            )

assert thetaFactor in (0.5, 1.0, 2.0)
theta = thetaFactor * pi

xmax = (rmax, rmax)
if thetaFactor == 0.5:
    xmin = (0.0, 0.0)
elif thetaFactor == 1.0:
    xmin = (-rmax, 0.0)
else:
    assert thetaFactor == 2.0
    xmin = (-rmax, -rmax)

rho0 = 1.0
eps0 = 0.0

restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Noh-cylindrical-2d-%ix%i" % (nRadial, nTheta))

vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Noh-cylindrical-2d-%ix%i" % (nRadial, nTheta),

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
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = TableKernel(BSplineKernel(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                             hmin = hmin,
                             hmax = hmax,
                             hminratio = hminratio,
                             nPerh = nPerh)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.hminratio")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
pos = nodes1.positions()
vel = nodes1.velocity()
if restoreCycle is None:
    if seed == "square":
        generator = GenerateSquareNodeDistribution(nRadial,
                                                   nTheta,
                                                   rho0,
                                                   xmin,
                                                   xmax,
                                                   nNodePerh = nPerh,
                                                   SPH = (HydroConstructor == SPHHydro))
    else:
        generator = GenerateNodeDistribution2d(nRadial, nTheta, rho0, seed,
                                               rmin = rmin,
                                               rmax = rmax,
                                               xmin = xmin,
                                               xmax = xmax,
                                               theta = theta,
                                               azimuthalOffsetFraction = azimuthalOffsetFraction,
                                               nNodePerh = nPerh,
                                               SPH = (HydroConstructor == SPHHydro))
##                                                relaxation = RadialCentroidalRelaxation(Vector(0,0),
##                                                                                        tolerance = 5.0e-4))

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
        #from PeanoHilbertDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes1, generator))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps0))

    # Set node velocities
    for nodeID in xrange(nodes1.numNodes):
        vel[nodeID] = pos[nodeID].unitVector()*vr0

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
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = HydroConstructor(WT,
                         WTPi,
                         q,
                         cfl = cfl,
                         compatibleEnergyEvolution = compatibleEnergy,
                         gradhCorrection = gradhCorrection,
                         XSPH = XSPH,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")
output("hydro.epsilonTensile")
output("hydro.nTensile")

packages = [hydro]

#-------------------------------------------------------------------------------
# Optionally construct an hourglass control object.
#-------------------------------------------------------------------------------
if hourglass:
    mask = db.newFluidIntFieldList(1, "mask")
    dx = rmax/nRadial
    bound = rmax - dx
    if seed == "square":
        for i in xrange(nodes1.numInternalNodes):
            if pos[i].x > bound or pos[i].y > bound:
                mask[0][i] = 0
    else:
        for i in xrange(nodes1.numInternalNodes):
            if pos[i].magnitude() > bound:
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
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
if hourglass:
    output("integrator.havePhysicsPackage(hg)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")
output("integrator.rigorousBoundaries")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
output("control")

# Do some startup stuff (unless we're restarting).
if restoreCycle is None:
    control.iterateIdealH(hydro)
    control.smoothState(smoothIters)
##     if hourglass:
##         print "Relaxing initial node distribution."
##         control.prerelaxNodeDistribution(hg, rho0)
    if densityUpdate in (VoronoiCellDensity, SumVoronoiCellDensity):
        print "Reinitializing node masses."
        control.voronoiInitializeMass()
    control.dropRestartFile()

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
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
if graphics:

    # Plot the node positions.
    import Gnuplot
    rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)

    # Plot the final state.
    rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    del HPlot
    Hinverse = db.newFluidSymTensorFieldList()
    db.fluidHinverse(Hinverse)
    hr = db.newFluidScalarFieldList()
    ht = db.newFluidScalarFieldList()
    for Hfield, hrfield, htfield in zip(Hinverse,
                                        hr,
                                        ht):
        n = Hfield.numElements
        assert hrfield.numElements == n
        assert htfield.numElements == n
        positions = Hfield.nodeList().positions()
        for i in xrange(n):
            runit = positions[i].unitVector()
            tunit = Vector(-(positions[i].y), positions[i].x).unitVector()
            hrfield[i] = (Hfield[i]*runit).magnitude()
            htfield[i] = (Hfield[i]*tunit).magnitude()
    hrPlot = plotFieldList(hr, xFunction="%s.magnitude()", plotStyle="points", winTitle="h_r")
    htPlot = plotFieldList(ht, xFunction="%s.magnitude()", plotStyle="points", winTitle="h_t")

    # Overplot the analytic solution.
    import NohAnalyticSolution
    answer = NohAnalyticSolution.NohSolution(2,
                                             h0 = nPerh*rmax/nRadial)
    plotAnswer(answer, control.time(),
               rhoPlot = rhoPlot,
               velPlot = vrPlot,
               epsPlot = epsPlot,
               PPlot = PPlot,
               HPlot = hrPlot)

    if mpi.rank == 0:
        r, hrans, htans = answer.hrtsolution(control.time())
        htData = Gnuplot.Data(r, htans,
                              title = "Solution",
                              with_ = "lines",
                              inline = "true")
        htPlot.replot(htData)

    # Report the error norms.
    rmin, rmax = 0.05, 0.35
    r = mpi.allreduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)
    rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
    v = mpi.allreduce([x.magnitude() for x in nodes1.velocity().internalValues()], mpi.SUM)
    eps = mpi.allreduce(list(nodes1.specificThermalEnergy().internalValues()), mpi.SUM)
    Pf = ScalarField("pressure", nodes1)
    nodes1.pressure(Pf)
    P = mpi.allreduce(list(Pf.internalValues()), mpi.SUM)
    if mpi.rank == 0:
        from SpheralGnuPlotUtilities import multiSort
        import Pnorm
        multiSort(r, rho, v, eps, P)
        rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), r)
        print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
        for (name, data, ans) in [("Mass Density", rho, rhoans),
                                  ("Pressure", P, Pans),
                                  ("Velocity", v, vans),
                                  ("Thermal E", eps, epsans)]:
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
if outputFile != "None":
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralGnuPlotUtilities import multiSort
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
    yprof = mpi.reduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in nodes1.Hfield().internalValues()], mpi.SUM)
    mof = mortonOrderIndicies(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, xprof, yprof, rhoprof, Pprof, vprof, epsprof, hprof)
        f = open(outputFile, "w")
        for xi, yi, rhoi, Pi, vi, epsi, hi, mi in zip(xprof, yprof, rhoprof, Pprof, vprof, epsprof, hprof, mo):
            f.write((7*"%16.12e " + "%i\n") % (xi, yi, rhoi, Pi, vi, epsi, hi, mi))
            f.write((7*"%i " + "\n") % (unpackElementUL(packElementDouble(xi)),
                                        unpackElementUL(packElementDouble(yi)),
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
