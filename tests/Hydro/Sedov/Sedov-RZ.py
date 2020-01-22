#-------------------------------------------------------------------------------
# The Sedov test case run in RZ symmetry.
#-------------------------------------------------------------------------------
import os, shutil, mpi
from SolidSpheralRZ import *
from SpheralTestUtilities import *

from GenerateNodeDistribution2d import *
if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("RZ hydro test -- Sedov problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(problem = "planar",     # one of (planar, cylindrical, spherical)
            KernelConstructor = NBSplineKernel,
            order = 5,

            n1 = 100,
            n2 = 20,

            nPerh = 1.35,

            gamma = 5.0/3.0,
            mu = 1.0,

            rho0 = 1.0,
            eps0 = 0.0,
            Espike = 1.0,

            solid = False,    # If true, use the fluid limit of the solid hydro option

            CRKSPH = False,
            PSPH = False,
            SPH = True,       # Choose the H advancement
            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
            Qconstructor = MonaghanGingoldViscosity,
            boolReduceViscosity = False,
            HopkinsConductivity = False,     # For PSPH
            nhQ = 5.0,
            nhL = 10.0,
            aMin = 0.1,
            aMax = 2.0,
            boolCullenViscosity = False,
            cullenUseHydroDerivatives = True,  # Reuse the hydro calculation of DvDx.
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
            linearConsistent = False,
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            Qhmult = 1.0,
            Cl = 1.0, 
            Cq = 1.0,
            etaCritFrac = 1.0,
            etaFoldFrac = 0.2,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.1,
            hminratio = 0.02,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 4.0,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,
            filter = 0.0,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = None,
            goalRadius = 0.8,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            vizCycle = None,
            vizTime = 0.1,
            vizDerivs = False,
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            QcorrectionOrder = LinearOrder,
            volumeType = RKSumVolume,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = False,
            correctVelocityGradient = True,
            domainIndependent = False,
            cullGhostNodes = True,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            clearDirectories = True,
            checkError = False,
            checkRestart = False,
            checkEnergy = False,
            restoreCycle = -1,
            restartStep = 10000,
            comparisonFile = "None",
            normOutputFile = "None",
            writeOutputLabel = True,

            graphics = True,
            )

outputFile = "Sedov-%s-RZ.gnu" % problem

assert not(boolReduceViscosity and boolCullenViscosity)
   
assert problem in ("planar", "cylindrical", "spherical")

if CRKSPH:
    if solid:
        if SPH:
            HydroConstructor = SolidCRKSPHHydro
        else:
            HydroConstructor = SolidACRKSPHHydro
    else:
        if SPH:
            HydroConstructor = CRKSPHHydro
        else:
            HydroConstructor = ACRKSPHHydro
    Qconstructor = CRKSPHMonaghanGingoldViscosity
    gradhCorrection = False
else:
   if solid:
      if SPH:
         HydroConstructor = SolidSPHHydro
      else:
         HydroConstructor = SolidASPHHydro
   else:
      if SPH:
         HydroConstructor = SPHHydro
      else:
         HydroConstructor = ASPHHydro

dataDir = os.path.join("dumps-%s-Sedov-RZ" % problem,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "filter=%f" % filter)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sedov-%s-RZ" % problem)

vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Sedov-%s-RZ" % problem

# Figure out what our goal time should be.
import SedovAnalyticSolution
h0 = 1.0/n1*nPerh
if problem == "planar":
    ndim = 1
elif problem == "cylindrical":
    ndim = 2
else:
    ndim = 3
answer = SedovAnalyticSolution.SedovSolution(nDim = ndim,
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
eos = GammaLawGasMKS(gamma, mu, minimumPressure=0.0)
strength = NullStrength()

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
    Wbase = NBSplineKernel(order)
else:
    Wbase = KernelConstructor()
WT = TableKernel(Wbase, 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    nodes1 = makeSolidNodeList("nodes1", eos, strength,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh,
                               xmin = (-10.0, -10.0),
                               xmax = ( 10.0,  10.0),
                               kernelExtent = kernelExtent)
else:
    nodes1 = makeFluidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh,
                               xmin = (-10.0, -10.0),
                               xmax = ( 10.0,  10.0),
                               kernelExtent = kernelExtent)
    
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if problem == "planar":
    nz = n1
    nr = n2
    z0, z1 = 0.0, 1.0
    r0, r1 = 0.0, 0.2
    rmin, rmax = None, None
elif problem == "cylindrical":
    nz = n2
    nr = n1
    z0, z1 = 0.0, 0.2
    r0, r1 = 0.0, 1.0
    rmin, rmax = None, None
else:
    assert problem == "spherical"
    nz = n1
    nr = n1
    rmin, rmax = 0.0, 1.0
    z0, z1 = 0.0, 1.0
    r0, r1 = 0.0, 1.0

generator = RZGenerator(GenerateNodeDistribution2d(nz, nr, rho0, "lattice",
                                                   xmin = (z0, r0),
                                                   xmax = (z1, r1),
                                                   rmin = rmin,
                                                   rmax = rmax,
                                                   nNodePerh = nPerh,
                                                   SPH = SPH))

distributeNodes2d((nodes1, generator))
output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

#-------------------------------------------------------------------------------
# Set the point source of energy.
#-------------------------------------------------------------------------------
pos = nodes1.positions()
vel = nodes1.velocity()
mass = nodes1.mass()
eps = nodes1.specificThermalEnergy()
H = nodes1.Hfield()
Esum = 0.0
dr = (r1 - r0)/nr
dz = (z1 - z0)/nz
msum = 0.0
if problem == "planar":
    epsi = 0.5*Espike/(rho0*dz)
    for i in xrange(nodes1.numInternalNodes):
        if pos[i].x < z0 + dz:
            eps[i] += epsi
            Esum += mass[i]*epsi
elif problem == "cylindrical":
    epsi = Espike/(rho0*pi*dr*dr)
    for i in xrange(nodes1.numInternalNodes):
        if pos[i].y < r0 + dr:
            eps[i] += epsi
            Esum += mass[i]*epsi
else:
    epsi = 0.5*Espike/(rho0*pi*dr*dr*dz)
    for i in xrange(nodes1.numInternalNodes):
        if pos[i].magnitude() < sqrt(dr*dr + dz*dz):
            eps[i] += epsi
            Esum += mass[i]*epsi
Eglobal = mpi.allreduce(Esum, mpi.SUM)
if problem == "planar":
    Eexpect = 0.5*Espike*pi*(r1*r1 - r0*r0)
elif problem == "cylindrical":
    Eexpect = Espike*(z1 - z0)
else:
    Eexpect = 0.5*Espike
print "Initialized a total energy of", Eglobal, Eexpect, Eglobal/Eexpect
assert fuzzyEqual(Eglobal, Eexpect)

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
q.QcorrectionOrder = QcorrectionOrder
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if CRKSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             XSPH = XSPH,
                             correctionOrder = correctionOrder,
                             volumeType = volumeType,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
    q.etaCritFrac = etaCritFrac
    q.etaFoldFrac = etaFoldFrac
else:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             gradhCorrection = gradhCorrection,
                             correctVelocityGradient = correctVelocityGradient,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             XSPH = XSPH,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile,
                             etaMinAxis = 1.0)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(q,nhQ,nhL,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(q,WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection,cullenUseHydroDerivatives)
    packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------
if bArtificialConduction:
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
if problem == "planar":
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z0, r1), Vector( 0.0, -1.0)))]
    if r0 != 0.0:
        bcs.append(ReflectingBoundary(Plane(Vector(z0, r0), Vector( 0.0, 1.0))))
elif problem == "cylindrical":
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z1, r0), Vector(-1.0,  0.0)))]
else:
    assert problem == "spherical"
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0)))]

for bc in bcs:
    for p in packages:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
del p
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.domainDecompositionIndependent")
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
                            vizTime = vizTime,
                            vizDerivs = vizDerivs,
                            SPH = SPH)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
   control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# Compute the analytic answer.
#-------------------------------------------------------------------------------
import mpi
import SedovAnalyticSolution
if problem == "planar":
    xprof = mpi.allreduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
elif problem == "cylindrical":
    xprof = mpi.allreduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
else:
    xprof = mpi.allreduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)

# Compute the simulated specific entropy.
rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(Pf.internalValues(), mpi.SUM)
A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

# Solution profiles.
xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]
L1 = 0.0
for i in xrange(len(rho)):
    L1 = L1 + abs(rho[i]-rhoans[i])
L1_tot = L1 / len(rho)
# if mpi.rank == 0 and outputFile != "None":
#     print "L1=",L1_tot,"\n"
#     with open("Converge.txt", "a") as myfile:
#         myfile.write("%s %s\n" % (nz, L1_tot))

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    if problem == "planar":
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, xFunction="%s.x", vecyFunction="%s.x", tenyFunction="1.0/%s.xx")
    elif problem == "cylindrical":
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, xFunction="%s.y", vecyFunction="%s.y", tenyFunction="1.0/%s.yy")
    else:
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    EPlot = plotEHistory(control.conserve)
    plots = [(rhoPlot, "Sedov-%s-rho-RZ.png" % problem),
             (velPlot, "Sedov-%s-vel-RZ.png" % problem),
             (epsPlot, "Sedov-%s-eps-RZ.png" % problem),
             (PPlot, "Sedov-%s-P-RZ.png" % problem),
             (HPlot, "Sedov-%s-h-RZ.png" % problem)]

    # Plot the specific entropy.
    Aplot = generateNewGnuPlot()
    AsimData = Gnuplot.Data(xprof, A,
                            with_ = "points",
                            title = "Simulation",
                            inline = True)
    AansData = Gnuplot.Data(xprof, Aans,
                            with_ = "lines",
                            title = "Solution",
                            inline = True)
    Aplot.plot(AsimData)
    Aplot.replot(AansData)
    Aplot.title("Specific entropy")
    Aplot.refresh()
    plots.append((Aplot, "Sedov-planar-A.png"))
    
    if CRKSPH:
        volPlot = plotFieldList(hydro.volume(), 
                                winTitle = "volume",
                                colorNodeLists = False, plotGhosts = False)
        plots.append((volPlot, "Sedov-%s-vol.png" % problem))

    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier(),
                                      winTitle = "Cullen alpha")
        cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt(),
                                       winTitle = "Cullen DalphaDt")
        plots += [(cullAlphaPlot, "Sedov-%s-Cullen-alpha.png" % problem),
                  (cullDalphaPlot, "Sedov-%s-Cullen-DalphaDt.png" % problem)]

    if boolReduceViscosity:
        alphaPlotQ = plotFieldList(q.reducingViscosityMultiplierQ(),
                                  winTitle = "rvAlphaQ",
                                  colorNodeLists = False, plotGhosts = False)
        alphaPlotL = plotFieldList(q.reducingViscosityMultiplierL(),
                                   winTitle = "rvAlphaL",
                                   colorNodeLists = False, plotGhosts = False)

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.hardcopy(os.path.join(dataDir, filename), terminal="png")


#-------------------------------------------------------------------------------
# Measure the difference between the simulation and analytic answer.
#-------------------------------------------------------------------------------
rmin, rmax = 0.05, 0.35   # Throw away anything with r < rwall to avoid wall heating.
rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
P = ScalarField("pressure", nodes1)
nodes1.pressure(P)
Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
xprof = mpi.reduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile != "None":
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralGnuPlotUtilities import multiSort
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)
    mprof = mpi.reduce(nodes1.mass().internalValues(), mpi.SUM)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
    if mpi.rank == 0:
        multiSort(xprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                  rhoans, Pans, vans, uans, hans)
        f = open(outputFile, "w")
        f.write(("#  " + 20*"'%s' " + "\n") % ("x", "m", "rho", "P", "v", "eps", "h", "mo",
                                               "rhoans", "Pans", "vans", "epsans", "hans",
                                               "x_UU", "m_UU", "rho_UU", "P_UU", "v_UU", "eps_UU", "h_UU"))
        for (xi, mi, rhoi, Pi, vi, epsi, hi, moi,
             rhoansi, Pansi, vansi, uansi, hansi) in zip(xprof, mprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                                                         rhoans, Pans, vans, uans, hans):
            f.write((7*"%16.12e " + "%i " + 5*"%16.12e " + 7*"%i " + '\n') % 
                    (xi, mi, rhoi, Pi, vi, epsi, hi, moi,
                     rhoansi, Pansi, vansi, uansi, hansi,
                     unpackElementUL(packElementDouble(xi)),
                     unpackElementUL(packElementDouble(mi)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi))))
        f.close()

        # #---------------------------------------------------------------------------
        # # Also we can optionally compare the current results with another file.
        # #---------------------------------------------------------------------------
        # if comparisonFile != "None":
        #     comparisonFile = os.path.join(dataDir, comparisonFile)
        #     import filecmp
        #     assert filecmp.cmp(outputFile, comparisonFile)

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print "Total energy error: %g" % Eerror
if checkEnergy and abs(Eerror) > 1e-13:
    raise ValueError, "Energy error outside allowed bounds."
