#-------------------------------------------------------------------------------
# The Planar Noh test case run in 1-D.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
#
# Ordinary SPH
#
#ATS:t0 = test(      SELF, "--graphics None --clearDirectories True  --checkError True   --restartStep 20", label="Planar Noh problem -- 1-D (serial)")
#ATS:t1 = testif(t0, SELF, "--graphics None --clearDirectories False --checkError False  --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem -- 1-D (serial) RESTART CHECK")
#ATS:t2 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-restartcheck' --restartStep 20", np=2, label="Planar Noh problem -- 1-D (parallel)")
#ATS:t3 = testif(t2, SELF, "--graphics None --clearDirectories False --checkError False --dataDir 'dumps-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=2, label="Planar Noh problem -- 1-D (parallel) RESTART CHECK")
#ATS:t4 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-1proc-reproducing.txt'", label="Planar Noh problem -- 1-D (serial reproducing test setup)")
#ATS:t5 = testif(t4, SELF, "--graphics None --clearDirectories False  --checkError True  --dataDir 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-4proc-reproducing.txt' --comparisonFile 'Noh-planar-1proc-reproducing.txt'", np=4, label="Planar Noh problem -- 1-D (4 proc reproducing test)")
#
# Ordinary solid SPH
#
#ATS:t100 = test(      SELF, "--solid True --graphics None --clearDirectories True  --checkError True   --restartStep 20", label="Planar Noh problem with solid SPH -- 1-D (serial)")
#ATS:t101 = testif(t100, SELF, "--solid True --graphics None --clearDirectories False --checkError False  --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with solid SPH -- 1-D (serial) RESTART CHECK")
#ATS:t102 = test(      SELF, "--solid True --graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-restartcheck' --restartStep 20", np=2, label="Planar Noh problem with solid SPH -- 1-D (parallel)")
#ATS:t103 = testif(t102, SELF, "--solid True --graphics None --clearDirectories False --checkError False --dataDir 'dumps-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=2, label="Planar Noh problem with solid SPH -- 1-D (parallel) RESTART CHECK")
#ATS:t104 = test(      SELF, "--solid True --graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-1proc-reproducing.txt'", label="Planar Noh problem with solid SPH -- 1-D (serial reproducing test setup)")
#ATS:t105 = testif(t104, SELF, "--solid True --graphics None --clearDirectories False  --checkError True  --dataDir 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-4proc-reproducing.txt' --comparisonFile 'Noh-planar-1proc-reproducing.txt'", np=4, label="Planar Noh  problem with solid SPH -- 1-D (4 proc reproducing test)")
#
# CRK
#
#ATS:t200 = test(      SELF, "--crksph True --cfl 0.25 --KernelConstructor NBSplineKernel --order 7 --nPerh 1.01 --Cl 2.0 --Cq 1.0 --graphics None --clearDirectories True --checkError False --restartStep 20 --steps 40", label="Planar Noh problem with CRK -- 1-D (serial)")
#ATS:t201 = testif(t200, SELF, "--crksph True --cfl 0.25 --KernelConstructor NBSplineKernel --order 7 --nPerh 1.01 --Cl 2.0 --Cq 1.0 --graphics None --clearDirectories False --checkError False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with CRK -- 1-D (serial) RESTART CHECK")
#ATS:t202 = test(      SELF, "--crksph True --cfl 0.25 --KernelConstructor NBSplineKernel --order 7 --nPerh 1.01 --Cl 2.0 --Cq 1.0 --graphics None --clearDirectories True  --checkError False  --dataDir 'dumps-planar-CRK-reproducing' --domainIndependent True --outputFile 'Noh-planar-1proc-reproducing.txt' --steps 100", label="Planar Noh problem with CRK -- 1-D (serial reproducing test setup)")
#ATS:t203 = testif(t202, SELF, "--crksph True --cfl 0.25 --KernelConstructor NBSplineKernel --order 7 --nPerh 1.01 --Cl 2.0 --Cq 1.0 --graphics None --clearDirectories False  --checkError False  --dataDir 'dumps-planar-CRK-reproducing' --domainIndependent True --outputFile 'Noh-planar-4proc-reproducing.txt' --steps 100 --comparisonFile 'Noh-planar-1proc-reproducing.txt'", np=4, label="Planar Noh problem with CRK -- 1-D (4 proc reproducing test)")
#
# PSPH
#
#ATS:t300 = test(      SELF, "--psph True --graphics None --clearDirectories True --checkError False --restartStep 20 --steps 40", label="Planar Noh problem with PSPH -- 1-D (serial)")
#ATS:t301 = testif(t300, SELF, "--psph True --graphics None --clearDirectories False --checkError False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with PSPH -- 1-D (serial) RESTART CHECK")

import os, shutil
from SolidSpheral1d import *
from SpheralTestUtilities import *

title("1-D integrated hydro test -- planar Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(KernelConstructor = NBSplineKernel,
            order = 5,

            nx1 = 100,
            rho1 = 1.0,
            eps1 = 0.0,
            smallPressure = False, #If set to True eps is not zero but small. 
            x0 = 0.0,
            x1 = 1.0,
            xwall = 0.0,
            nPerh = 1.35,
            NeighborType = NestedGridNeighbor,

            vr0 = -1.0, 
            vrSlope = 0.0,

            gamma = 5.0/3.0,
            mu = 1.0,

            solid = False,    # If true, use the fluid limit of the solid hydro option

            svph = False,
            crksph = False,
            psph = False,
            crktype = "default",        # one of ("default", "variant")
            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
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
            Cl = None, 
            Cq = None,
            linearInExpansion = None,
            Qlimiter = None,
            balsaraCorrection = None,
            epsilon2 = None,
            QcorrectionOrder = None,
            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4.0,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,
            filter = 0.0,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            updateBoundaryFrequency = 1,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            volumeType = CRKSumVolume,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,
            domainIndependent = True,
            cullGhostNodes = True,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            clearDirectories = True,
            checkError = False,
            checkRestart = False,
            checkEnergy = True,
            restoreCycle = None,
            restartStep = 10000,
            dataDirBase = "dumps-planar-Noh",
            restartBaseName = "Noh-planar-1d",
            outputFile = "None",
            comparisonFile = "None",
            normOutputFile = "None",
            writeOutputLabel = True,

            # Parameters for the test acceptance.,
            L1rho =   0.0717692               ,
            L2rho =   0.0195334               ,
            Linfrho = 2.31412                 ,
                                               
            L1P =     0.0248091               ,
            L2P =     0.00728801              ,
            LinfP =   0.968204                ,
                                               
            L1v =     0.0392544               ,
            L2v =     0.0105355               ,
            Linfv =   0.946074                ,
                                               
            L1eps =   0.01416                 ,
            L2eps =   0.00401815              ,
            Linfeps = 0.438217                ,
                                               
            L1h =     0.000552337             ,
            L2h =     0.000132313             ,
            Linfh =   0.00915864              ,

            tol = 1.0e-5,

            graphics = True,
            )
assert not(boolReduceViscosity and boolCullenViscosity)
if smallPressure:
    P0 = 1.0e-6
    eps1 = P0/((gamma - 1.0)*rho1)
   
if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
else:
    hydroname = "SPH"
if solid:
    hydroname = "Solid" + hydroname

dataDir = os.path.join(dataDirBase,
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "CRKVar=%s" % crktype,
                       "filter=%f" % filter)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Noh-planar-1d-%i" % nx1)

dx = (x1 - x0)/nx1

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

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
    nodes1 = makeSolidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               nPerh = nPerh,
                               kernelExtent = kernelExtent,
                               NeighborType = NeighborType)
else:
    nodes1 = makeFluidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               nPerh = nPerh,
                               kernelExtent = kernelExtent,
                               NeighborType = NeighborType)
    
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))],
                         nPerh = nPerh)
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
nodes1.massDensity(ScalarField("tmp", nodes1, rho1))

# Set node velocities
pos = nodes1.positions()
vel = nodes1.velocity()
for ix in xrange(nodes1.numNodes):
    if pos[ix].x > xwall:
        vel[ix].x = vr0 + vrSlope*pos[ix].x
    else:
        vel[ix].x = -vr0 + vrSlope*pos[ix].x

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
                 useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution = compatibleEnergy,
                 densityUpdate = densityUpdate,
                 XSVPH = XSPH,
                 linearConsistent = linearConsistent,
                 generateVoid = False,
                 HUpdate = HUpdate,
                 fcentroidal = fcentroidal,
                 fcellPressure = fcellPressure,
                 xmin = Vector(-100.0),
                 xmax = Vector( 100.0))
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   filter = filter,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   correctionOrder = correctionOrder,
                   volumeType = volumeType,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   crktype = crktype)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 filter = filter,
                 cfl = cfl,
                 useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH)
else:
    hydro = SPH(dataBase = db,
                W = WT,
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
                nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.XSPH")

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
if QcorrectionOrder:
    q.QcorrectionOrder = QcorrectionOrder
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
try:
    q.linearInExpansion = linearInExpansion
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
# Optionally construct an hourglass control object.
#-------------------------------------------------------------------------------
if hourglass:
    mask = db.newFluidIntFieldList(1, "mask")
    pos = nodes1.positions()
    for i in xrange(nodes1.numInternalNodes):
        if pos[i].x > (x1 - dx):
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
if x0 == xwall:
    xPlane0 = Plane(Vector(0.0), Vector(1.0))
    xbc0 = ReflectingBoundary(xPlane0)
    for p in packages:
        p.appendBoundary(xbc0)

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
integrator.updateBoundaryFrequency = updateBoundaryFrequency
integrator.domainDecompositionIndependent = domainIndependent
integrator.cullGhostNodes = cullGhostNodes
integrator.verbose = dtverbose
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.updateBoundaryFrequency")
output("integrator.domainDecompositionIndependent")
output("integrator.cullGhostNodes")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
output("control")

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
    if control.time() < goalTime:
        control.step(5)
        control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# Compute the analytic answer.
#-------------------------------------------------------------------------------
import mpi
import NohAnalyticSolution
rlocal = [pos.x for pos in nodes1.positions().internalValues()]
r = mpi.reduce(rlocal, mpi.SUM)
h1 = 1.0/(nPerh*dx)
answer = NohAnalyticSolution.NohSolution(1,
                                         r = r,
                                         v0 = -1.0,
                                         h0 = 1.0/h1)

# Compute the simulated specific entropy.
rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(Pf.internalValues(), mpi.SUM)
A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

# The analytic solution for the simulated entropy.
xprof = mpi.allreduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]
L1 = 0.0
for i in xrange(len(rho)):
  L1 = L1 + abs(rho[i]-rhoans[i])
L1_tot = L1 / len(rho)
if mpi.rank == 0 and outputFile != "None":
 print "L1=",L1_tot,"\n"
 with open("Converge.txt", "a") as myfile:
    myfile.write("%s %s\n" % (nx1, L1_tot))

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
    plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    EPlot = plotEHistory(control.conserve)
    plots = [(rhoPlot, "Noh-planar-rho.png"),
             (velPlot, "Noh-planar-vel.png"),
             (epsPlot, "Noh-planar-eps.png"),
             (PPlot, "Noh-planar-P.png"),
             (HPlot, "Noh-planar-h.png")]

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
    plots.append((Aplot, "Noh-planar-A.png"))
    
    if crksph:
        volPlot = plotFieldList(hydro.volume(), 
                                winTitle = "volume",
                                colorNodeLists = False, plotGhosts = False)
        aplot = plotFieldList(hydro.A(),
                              winTitle = "A",
                              colorNodeLists = False)
        bplot = plotFieldList(hydro.B(),
                              yFunction = "%s.x",
                              winTitle = "B",
                              colorNodeLists = False)
        splot = plotFieldList(hydro.surfacePoint(),
                              winTitle = "surface point",
                              colorNodeLists = False)
        voidplot = plotFieldList(hydro.voidPoint(),
                                 winTitle = "void point",
                                 plotStyle = "points",
                                 plotGhosts = True,
                                 colorNodeLists = False)
        plots += [(volPlot, "Noh-planar-vol.png"),
                   (aplot, "Noh-planar-ACRK.png"),
                   (bplot, "Noh-planar-BCRK.png"),
                   (splot, "Noh-planar-surfacePoint.png"),
                   (voidplot, "Noh-planar-voidPoint.png")]

    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier(),
                                      winTitle = "Cullen alpha")
        cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt(),
                                       winTitle = "Cullen DalphaDt")
        plots += [(cullAlphaPlot, "Noh-planar-Cullen-alpha.png"),
                  (cullDalphaPlot, "Noh-planar-Cullen-DalphaDt.png")]

    if boolReduceViscosity:
        alphaPlotQ = plotFieldList(q.reducingViscosityMultiplierQ(),
                                  winTitle = "rvAlphaQ",
                                  colorNodeLists = False, plotGhosts = False)
        alphaPlotL = plotFieldList(q.reducingViscosityMultiplierL(),
                                   winTitle = "rvAlphaL",
                                   colorNodeLists = False, plotGhosts = False)

    # # Plot the grad h correction term (omega)
    # omegaPlot = plotFieldList(hydro.omegaGradh(),
    #                           winTitle = "grad h correction",
    #                           colorNodeLists = False)

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
xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)

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
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof)
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

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file.
        #---------------------------------------------------------------------------
        if comparisonFile != "None":
            comparisonFile = os.path.join(dataDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile, comparisonFile)

#------------------------------------------------------------------------------
# Compute the error.
#------------------------------------------------------------------------------
if mpi.rank == 0:
    xans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
    import Pnorm
    print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
    failure = False
    hD = []

    if normOutputFile != "None":
       f = open(normOutputFile, "a")
       if writeOutputLabel:
          f.write(("#" + 13*"%17s " + "\n") % ('"nx"',
                                               '"rho L1"', '"rho L2"', '"rho Linf"',
                                               '"P L1"',   '"P L2"',   '"P Linf"',
                                               '"vel L1"', '"vel L2"', '"vel Linf"',
                                               '"E L1"', '"E L2"', '"E Linf"',
                                               '"h L1"',   '"h L2"',   '"h Linf"'))
       f.write("%5i " % nx1)
    for (name, data, ans,
         L1expect, L2expect, Linfexpect) in [("Mass Density", rhoprof, rhoans, L1rho, L2rho, Linfrho),
                                             ("Pressure", Pprof, Pans, L1P, L2P, LinfP),
                                             ("Velocity", vprof, vans, L1v, L2v, Linfv),
                                             ("Thermal E", epsprof, epsans, L1eps, L2eps, Linfeps),
                                             ("h       ", hprof, hans, L1h, L2h, Linfh)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in xrange(len(data))]
        Pn = Pnorm.Pnorm(error, xprof)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print "\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf)
        if normOutputFile != "None":
           f.write((3*"%16.12e ") % (L1, L2, Linf))
        hD.append([L1,L2,Linf])
        
           

        if checkError:
            if not fuzzyEqual(L1, L1expect, tol):
                print "L1 error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                      L1,
                                                                                      L1expect)
                failure = True
            if not fuzzyEqual(L2, L2expect, tol):
                print "L2 error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                      L2,
                                                                                      L2expect)
                failure = True
            if not fuzzyEqual(Linf, Linfexpect, tol):
                print "Linf error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                        Linf,
                                                                                        Linfexpect)
                failure = True
            if failure:
                raise ValueError, "Error bounds violated."
    if normOutputFile != "None":
       f.write("\n")
                                             
    # print "%d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t" % (nx1,hD[0][0],hD[1][0],hD[2][0],hD[3][0],
    #                                                                             hD[0][1],hD[1][1],hD[2][1],hD[3][1],
    #                                                                             hD[0][2],hD[1][2],hD[2][2],hD[3][2])





control.conserve.writeHistory("Noh1dEnergyProfiles.gnu")
Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print "Total energy error: %g" % Eerror
if checkEnergy and abs(Eerror) > 1e-13:
    raise ValueError, "Energy error outside allowed bounds."
