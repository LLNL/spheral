#-------------------------------------------------------------------------------
# The spherical Sedov test case (3-D).
#-------------------------------------------------------------------------------
import os, sys, shutil
from Spheral3d import *
from findLastRestart import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from GenerateNodeDistribution3d import *

import mpi

title("3-D integrated hydro test -- planar Sedov problem")
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            nx = 50,
            ny = 50,
            nz = 50,
            nPerh = 1.51,
            KernelConstructor = BSplineKernel,
            order = 5,

            rho0 = 1.0,
            eps0 = 0.0,
            smallPressure = False,
            Espike = 1.0,
            smoothSpike = True,
            topHatSpike = False,
            smoothSpikeScale = 0.5,
            gamma = 5.0/3.0,
            mu = 1.0,

            Cl = 1.0,
            Cq = 0.75,
            epsilon2 = 1e-2,
            Qlimiter = False,
            balsaraCorrection = False,
            linearInExpansion = False,

            crksph = False,
            psph = False,
            asph = False,     # Only for H evolution, not hydro algorithm
            Qconstructor = MonaghanGingoldViscosity,
            correctionOrder = LinearOrder,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            HUpdate = IdealH,
            filter = 0.0,
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qhmult = 1.0,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
            HopkinsConductivity = False,     # For PSPH
            evolveTotalEnergy = False,       # Only for SPH variants -- evolve total rather than specific energy
            volumeType = RKSumVolume,       # For CRK

            hmin = 1e-15,
            hmax = 1.0,
            cfl = 0.5,
            useVelocityMagnitudeForDt = True,
            XSPH = False,
            rhomin = 1e-10,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
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
            gradhCorrection = True,
            correctVelocityGradient = True,

            restoreCycle = -1,
            restartStep = 1000,

            graphics = True,
            clearDirectories = False,
            dataRoot = "dumps-spherical-Sedov",
            outputFile = "None",
            )

if smallPressure:
    P0 = 1.0e-6
    eps0 = P0/((gamma - 1.0)*rho0)
    print "WARNING: smallPressure specified, so setting eps0=%g" % eps0

assert not(boolReduceViscosity and boolCullenViscosity)

# Figure out what our goal time should be.
import SedovAnalyticSolution
h0 = 1.0/nx*nPerh
answer = SedovAnalyticSolution.SedovSolution(nDim = 3,
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

#-------------------------------------------------------------------------------
# Path names.
#-------------------------------------------------------------------------------
hydroname = ""
if solid:
    hydroname += "Solid"
if asph:
    hydroname += "A"
if crksph:
    hydroname += "CRKSPH"
elif psph:
    hydroname += "PSPH"
else:
    hydroname += "SPH"

dataDir = os.path.join(dataRoot,
                       hydroname,
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "seed=%s" % seed,
                       "nx=%i_ny=%i_nz=%i" % (nx, ny, nz))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Sedov-spherical-3d")

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
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
else:
  WT = TableKernel(KernelConstructor(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Create a NodeList and associated Neighbor object.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           kernelExtent = kernelExtent,
                           nPerh = nPerh,
                           rhoMin = rhomin)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if seed.lower() == "icosahedral":
    generator = GenerateIcosahedronMatchingProfile3d(nx, # Sets nradial
                                                     rho0, 
                                                     rmin = 0.0,
                                                     rmax = 1.0,
                                                     nNodePerh = nPerh)
else:
    generator = GenerateNodeDistribution3d(nx, ny, nz,
                                           rho0, seed,
                                           xmin = (0.0, 0.0, 0.0),
                                           xmax = (1.0, 1.0, 1.0),
                                           rmin = 0.0,
                                           rmax = 1.0,
                                           nNodePerh = nPerh,
                                           SPH = (not asph))

if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

distributeNodes3d((nodes1, generator))
output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

# Set the point source of energy.
if seed != "icosahedral":
    Espike /= 8.0  # Doing an octant
pos = nodes1.positions()
vel = nodes1.velocity()
mass = nodes1.mass()
eps = nodes1.specificThermalEnergy()
H = nodes1.Hfield()
Esum = 0.0
if smoothSpike or topHatSpike:
    Wsum = 0.0
    for nodeID in xrange(nodes1.numInternalNodes):
        Hi = H[nodeID]
        etaij = (Hi*pos[nodeID]).magnitude()
        if smoothSpike:
            Wi = WT.kernelValue(etaij/smoothSpikeScale, 1.0)
        else:
            if etaij < smoothSpikeScale*kernelExtent:
                Wi = 1.0
            else:
                Wi = 0.0
        Ei = Wi*Espike
        epsi = Ei/mass[nodeID]
        eps[nodeID] = epsi
        Wsum += Wi
    Wsum = mpi.allreduce(Wsum, mpi.SUM)
    assert Wsum > 0.0
    for nodeID in xrange(nodes1.numInternalNodes):
        eps[nodeID] /= Wsum
        Esum += eps[nodeID]*mass[nodeID]
        eps[nodeID] += eps0
else:
    i = -1
    rmin = 1e50
    for nodeID in xrange(nodes1.numInternalNodes):
        rij = pos[nodeID].magnitude()
        if rij < rmin:
            i = nodeID
            rmin = rij
        eps[nodeID] = eps0
    rminglobal = mpi.allreduce(rmin, mpi.MIN)
    if fuzzyEqual(rmin, rminglobal):
        assert i >= 0 and i < nodes1.numInternalNodes
        eps[i] += Espike/mass[i]
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
if crksph:
    hydro = CRKSPH(dataBase = db,
                   Q = q,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   correctionOrder = correctionOrder,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 Q = q,
                 filter = filter,
                 cfl = cfl,
                 useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 HopkinsConductivity = HopkinsConductivity,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH)
else:
    hydro = SPH(dataBase = db,
                W = WT, 
                Q = q,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                XSPH = XSPH,
                HUpdate = HEvolution)
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(q,nh,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(q,WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
    packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
if seed.lower() != "icosahedral":
    xPlane0 = Plane(Vector(0, 0, 0), Vector(1, 0, 0))
    yPlane0 = Plane(Vector(0, 0, 0), Vector(0, 1, 0))
    zPlane0 = Plane(Vector(0, 0, 0), Vector(0, 0, 1))
    xbc0 = ReflectingBoundary(xPlane0)
    ybc0 = ReflectingBoundary(yPlane0)
    zbc0 = ReflectingBoundary(zPlane0)

    for p in packages:
        for bc in (xbc0, ybc0, zbc0):
            p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
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
control = SpheralController(integrator, WT,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = "Sedov-spherical-3d-%ix%ix%i" % (nx, ny, nz),
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = (not ASPH))
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
zprof = mpi.allreduce([x.z for x in nodes1.positions().internalValues()], mpi.SUM)
rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
mass = mpi.allreduce(list(nodes1.mass().internalValues()), mpi.SUM)
v = mpi.allreduce([x.magnitude() for x in nodes1.velocity().internalValues()], mpi.SUM)
eps = mpi.allreduce(list(nodes1.specificThermalEnergy().internalValues()), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(list(Pf.internalValues()), mpi.SUM)
A = mpi.allreduce([Pi/(rhoi**gamma) for (Pi, rhoi) in zip(Pf.internalValues(), nodes1.massDensity().internalValues())], mpi.SUM)

rans, vans, epsans, rhoans, Pans, Aans, hans = answer.solution(control.time(), r)
from SpheralGnuPlotUtilities import multiSort
multiSort(r, rho, v, eps, P, A, rhoans, vans, epsans, Pans, hans)

if mpi.rank == 0:
    from SpheralGnuPlotUtilities import multiSort
    import Pnorm
    multiSort(r, rho, v, eps, P, A)
    print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
    for (name, data, ans) in [("Mass Density", rho, rhoans),
                              ("Pressure", P, Pans),
                              ("Velocity", v, vans),
                              ("Thermal E", eps, epsans),
                              ("Entropy", A, Aans)]:
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
    f.write(("# " + 16*"%15s " + "\n") % ("r", "x", "y", "z", "rho", "m", "P", "v", "eps", "A",
                                          "rhoans", "Pans", "vans", "epsans", "Aans", "hrans"))
    for (ri, xi, yi, zi, rhoi, mi, Pi, vi, epsi, Ai, 
         rhoansi, Pansi, vansi, epsansi, Aansi, hansi)  in zip(r, xprof, yprof, zprof, rho, mass, P, v, eps, A,
                                                               rhoans, Pans, vans, epsans, Aans, hans):
         f.write((16*"%16.12e " + "\n") % (ri, xi, yi, zi, rhoi, mi, Pi, vi, epsi, Ai,
                                           rhoansi, Pansi, vansi, epsansi, Aansi, hansi))
    f.close()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    plots = [(rhoPlot, "Sedov-spherical-rho.png"),
             (velPlot, "Sedov-spherical-vel.png"),
             (epsPlot, "Sedov-spherical-eps.png"),
             (PPlot, "Sedov-spherical-P.png"),
             (HPlot, "Sedov-spherical-h.png")]

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
    plots.append((Aplot, "Sedov-spherical-entropy.png"))

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.hardcopy(os.path.join(dataDir, filename), terminal="png")

