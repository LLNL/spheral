#
# SPH
#
#ATS:sph0 = test(        SELF, "--crksph False --nRadial 100 --cfl 0.25 --Cl 1.0 --Cq 1.0 --filter 0.0 --nPerh 2.01 --graphics False --restartStep 20 --clearDirectories True --steps 100", label="Noh cylindrical SPH, nPerh=2.0", np=2)
#ATS:sph1 = testif(sph0, SELF, "--crksph False --nRadial 100 --cfl 0.25 --Cl 1.0 --Cq 1.0 --filter 0.0 --nPerh 2.01 --graphics False --restartStep 20 --clearDirectories False --steps 60 --restoreCycle 40 --checkRestart True", label="Noh cylindrical SPH, nPerh=2.0, restart test", np=2)
#
# CRK
#
#ATS:crk0 = test(        SELF, "--crksph True --nRadial 100 --cfl 0.25 --Cl 1.0 --Cq 1.0 --filter 0.0 --nPerh 2.01 --graphics False --restartStep 20 --volumeType CRKVoronoiVolume --clearDirectories True --steps 100", label="Noh cylindrical CRK, nPerh=2.0", np=2)
#ATS:crk1 = testif(crk0, SELF, "--crksph True --nRadial 100 --cfl 0.25 --Cl 1.0 --Cq 1.0 --filter 0.0 --nPerh 2.01 --graphics False --restartStep 20 --volumeType CRKVoronoiVolume --clearDirectories False --steps 60 --restoreCycle 40 --checkRestart True", label="Noh cylindrical CRK, nPerh=2.0, restart test", np=2)

#-------------------------------------------------------------------------------
# The Cylindrical Noh test case run in 2-D.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
import os, shutil
from math import *
from SolidSpheral2d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution2d import *
from CubicNodeGenerator import GenerateSquareNodeDistribution
from CentroidalVoronoiRelaxation import *

import mpi
import DistributeNodes

title("2-D integrated hydro test -- cylindrical Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(order = 5,
            seed = "constantDTheta",

            thetaFactor = 0.5,
            azimuthalOffsetFraction = 0.0,
            nRadial = 50,
            nTheta = 50,
            rmin = 0.0,
            rmax = 1.0,
            nPerh = 2.01,
            rho0 = 1.0,
            eps0 = 0.0,
            smallPressure = False,

            vr0 = -1.0, 

            gamma = 5.0/3.0,
            mu = 1.0,

            solid = False,    # If true, use the fluid limit of the solid hydro option

            svph = False,
            crksph = False,
            psph = False,
            asph = False,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            boolReduceViscosity = False,
            HopkinsConductivity = False,     # For PSPH
            nhQ = 5.0,
            nhL = 10.0,
            aMin = 0.1,
            aMax = 2.0,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
            linearConsistent = False,

            Cl = None, 
            Cq = None,
            linearInExpansion = None,
            Qlimiter = None,
            balsaraCorrection = None,
            epsilon2 = None,
            hmin = 0.0001, 
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            filter = 0.0,

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
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            volumeType = CRKSumVolume,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,

            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,

            useVoronoiOutput = True,
            clearDirectories = False,
            vizDerivs = False,
            restoreCycle = -1,
            restartStep = 1000,
            checkRestart = False,
            dataDir = "dumps-cylindrical-Noh",
            outputFile = "None",
            comparisonFile = "None",

            graphics = True,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
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

if smallPressure:
   P0 = 1.0e-6
   eps0 = P0/((gamma - 1.0)*rho0)

if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
else:
    hydroname = "SPH"
if asph:
    hydroname = "A" + hydroname

if solid:
    hydroname = "Solid" + hydroname

dataDir = os.path.join(dataDir,
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "filter=%f" % filter,
                       "nrad=%i_ntheta=%i" % (nRadial, nTheta))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Noh-cylindrical-2d-%ix%i" % (nRadial, nTheta))

vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Noh-cylindrical-2d-%ix%i" % (nRadial, nTheta)

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
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(order), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    nodes1 = makeSolidNodeList("nodes1", eos,
                               hmin = hmin,
                               hmax = hmax,
                               kernelExtent = kernelExtent,
                               hminratio = hminratio,
                               nPerh = nPerh)
else:
    nodes1 = makeFluidNodeList("nodes1", eos,
                               hmin = hmin,
                               hmax = hmax,
                               kernelExtent = kernelExtent,
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
if seed == "square":
    generator = GenerateSquareNodeDistribution(nRadial,
                                               nTheta,
                                               rho0,
                                               xmin,
                                               xmax,
                                               nNodePerh = nPerh,
                                               SPH = not asph)
else:
    generator = GenerateNodeDistribution2d(nRadial, nTheta, rho0, seed,
                                           rmin = rmin,
                                           rmax = rmax,
                                           xmin = xmin,
                                           xmax = xmax,
                                           theta = theta,
                                           azimuthalOffsetFraction = azimuthalOffsetFraction,
                                           nNodePerh = nPerh,
                                           SPH = not asph)

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
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if svph:
    hydro = SVPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 densityUpdate = densityUpdate,
                 XSVPH = XSPH,
                 linearConsistent = linearConsistent,
                 generateVoid = False,
                 HUpdate = HUpdate,
                 fcentroidal = fcentroidal,
                 fcellPressure = fcellPressure,
                 xmin = Vector(-1.1, -1.1),
                 xmax = Vector( 1.1,  1.1),
                 ASPH = asph)
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   correctionOrder = correctionOrder,
                   volumeType = volumeType,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 filter = filter,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 HopkinsConductivity = HopkinsConductivity,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 ASPH = asph)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                filter = filter,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                ASPH = asph)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Set the artificial viscosity parameters.
#-------------------------------------------------------------------------------
q = hydro.Q
if Cl:
    q.Cl = Cl
if Cq:
    q.Cq = Cq
if epsilon2:
    q.epsilon2 = epsilon2
if Qlimiter:
    q.limiter = Qlimiter
if balsaraCorrection:
    q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
try:
    output("q.linearInExpansion")
    output("q.quadraticInExpansion")
except:
    pass

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(q,nhQ,nhL,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(q,WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
    packages.append(evolveCullenViscosityMultiplier)

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
if useVoronoiOutput:
    import SpheralVoronoiSiloDump
    vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
else:
    import SpheralPointmeshSiloDump
    vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState
    #import SpheralVisitDump
    #vizMethod = SpheralVisitDump.dumpPhysicsState

control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizMethod = vizMethod,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizDerivs = vizDerivs,
                            #skipInitialPeriodicWork = SVPH,
                            SPH = True,        # Only for iterating H
                            )
output("control")

# Do some startup stuff (unless we're restarting).
if restoreCycle is None:
    control.smoothState(smoothIters)
    if densityUpdate in (VoronoiCellDensity, SumVoronoiCellDensity):
        print "Reinitializing node masses."
        control.voronoiInitializeMass()
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    if checkRestart:
        control.setRestartBaseName(restartBaseName + "_CHECK")
    control.step(steps)
    if checkRestart:
        control.setRestartBaseName(restartBaseName)

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
import NohAnalyticSolution
answer = NohAnalyticSolution.NohSolution(2,
                                         h0 = nPerh*rmax/nRadial)

if graphics:
    # Plot the node positions.
    from SpheralMatplotlib import *
    rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)
    EPlot = plotEHistory(control.conserve)

    # Plot the final state.
    rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
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
    hrPlot = plotFieldList(hr, xFunction="%s.magnitude()", plotStyle="ro", winTitle="$h_r$")
    htPlot = plotFieldList(ht, xFunction="%s.magnitude()", plotStyle="ro", winTitle="$h_t$")

    # Overplot the analytic solution.
    plotAnswer(answer, control.time(),
               rhoPlot = rhoPlot,
               velPlot = vrPlot,
               epsPlot = epsPlot,
               PPlot = PPlot,
               HPlot = hrPlot)

    if boolReduceViscosity:
        alphaPlotQ = plotFieldList(q.reducingViscosityMultiplierQ(),
                                   xFunction = "%s.magnitude()",
                                   winTitle = "rvAlphaQ",
                                   colorNodeLists = False, plotGhosts = False)
        alphaPlotL = plotFieldList(q.reducingViscosityMultiplierL(),
                                   xFunction = "%s.magnitude()",
                                   winTitle = "rvAlphaL",
                                   colorNodeLists = False, plotGhosts = False)

    if mpi.rank == 0:
        r, hrans, htans = answer.hrtsolution(control.time())
        htPlot.plot(r, htans, "b-", label="Solution", )
        htPlot.axes.legend()

    plots = [(rPlot, "Noh-cylindrical-positions.png"),
             (rhoPlot, "Noh-cylindrical-rho.png"),
             (vrPlot, "Noh-cylindrical-vel.png"),
             (epsPlot, "Noh-cylindrical-eps.png"),
             (PPlot, "Noh-cylindrical-P.png"),
             (hrPlot, "Noh-cylindrical-hr.png"),
             (htPlot, "Noh-cylindrical-ht.png")]

    if crksph:
        volPlot = plotFieldList(hydro.volume(), 
                                xFunction = "%s.magnitude()",
                                winTitle = "volume",
                                plotStyle = "ro",
                                colorNodeLists = False, plotGhosts = False)
        spPlot = plotFieldList(hydro.surfacePoint(), 
                               xFunction = "%s.magnitude()",
                               winTitle = "Surface",
                               plotStyle = "ro",
                               colorNodeLists = False, plotGhosts = False)
        vpPlot = plotFieldList(hydro.voidPoint(), 
                               xFunction = "%s.magnitude()",
                               winTitle = "Void",
                               plotStyle = "ro",
                               colorNodeLists = False, plotGhosts = True)
        plots += [(volPlot, "Noh-cylindrical-vol.png"),
                  (spPlot, "Noh-cylindrical-surfacePoint.png"),
                  (vpPlot, "Noh-cylindrical-voidPoint.png")]

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.figure.savefig(os.path.join(dataDir, filename))

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
        from SpheralTestUtilities import multiSort
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
    #vprof = mpi.reduce(list([vi.dot(ri.unitVector()) for ri,vi in zip(nodes1.positions().internalValues(),nodes1.velocity().internalValues())]),mpi.SUM)
    rprof = mpi.reduce([ri.magnitude() for ri in nodes1.positions().internalValues()],mpi.SUM)
    vx = mpi.reduce(list([v.x for v in nodes1.velocity().internalValues()]),mpi.SUM)
    vy = mpi.reduce([v.y for v in nodes1.velocity().internalValues()],mpi.SUM)
    np = int(nodes1.numInternalNodes)
    if np is None:
        np = 0
    #print "np=%d" % np
    np = mpi.reduce(np,mpi.SUM)
    #print "np=%d" % np
    vprof = []
    if mpi.rank == 0:
        for i in xrange(np):
            vprof.append(xprof[i]*vx[i]/rprof[i] + yprof[i]*vy[i]/rprof[i])
    #vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    Qprof = mpi.reduce(hydro.viscousWork()[0].internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in nodes1.Hfield().internalValues()], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)
    if mpi.rank == 0:
        #rprof = [sqrt(xi*xi + yi*yi) for xi, yi in zip(xprof, yprof)]
        multiSort(rprof, mo, xprof, yprof, rhoprof, Pprof, vprof, epsprof, hprof, Qprof)
        rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), rprof)
        f = open(outputFile, "w")
        f.write(("# " + 21*"%15s " + "\n") % ("r", "x", "y", "rho", "P", "v", "eps", "h", "mortonOrder", "QWork",
                                              "rhoans", "Pans", "vans", "epsans",
                                              "x_uu", "y_uu", "rho_uu", "P_uu", "v_uu", "eps_uu", "h_uu"))
        for (ri, xi, yi, rhoi, Pi, vi, epsi, hi, mi, Qi,
             rhoansi, Pansi, vansi, epsansi)  in zip(rprof, xprof, yprof, rhoprof, Pprof, vprof, epsprof, hprof, mo, Qprof,
                                                     rhoans, Pans, vans, epsans):
            f.write((8*"%16.12e " + "%i " + 5*"%16.12e " + 7*"%i " + "\n") % (ri, xi, yi, rhoi, Pi, vi, epsi, hi, mi, Qi,
                                                                              rhoansi, Pansi, vansi, epsansi,
                                                                              unpackElementUL(packElementDouble(xi)),
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
Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/max(1.0e-30, control.conserve.EHistory[0])
print "Total energy error: %g" % Eerror
if compatibleEnergy and abs(Eerror) > 1e-13:
    raise ValueError, "Energy error outside allowed bounds."
