#ATS:t1 = test(      SELF, "--crksph True --cfl 0.25 --graphics None --clearDirectories True  --restartStep 20 --steps 40", label="Planar Sod problem with CRK -- 1-D (serial)")
#ATS:t2 = testif(t1, SELF, "--crksph True --cfl 0.25 --graphics None --clearDirectories False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Sod problem with CRK -- 1-D (serial) RESTART CHECK")
import os, sys
import shutil
from SolidSpheral1d import *
from SpheralTestUtilities import *
from SodAnalyticSolution import *

title("1-D integrated hydro test -- planar Sod problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 400,
            nx2 = 100,
            rho1 = 1.0,
            rho2 = 0.25,
            P1 = 1.0,
            P2 = 0.1795,

            numNodeLists = 1,

            x0 = -0.5,
            x1 = 0.0,
            x2 = 0.5,

            hsmooth = 0.5,             # Optionally smooth initial discontinuity, expressed as particle spacings
            sumInitialDensity = False, # Optionally sum the initial density before setting the pressure and such

            nPerh = 1.35,

            gammaGas = 5.0/3.0,
            mu = 1.0,
            
            svph = False,
            crksph = False,
            psph = False,
            crktype = "default",        # one of ("default", "variant")
            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
            solid = False,    # If true, use the fluid limit of the solid hydro option
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Cl = None,
            Cq = None,
            Qlimiter = None,
            epsilon2 = None,
            etaCritFrac = None,
            etaFoldFrac = None,
            linearInExpansion = None,
            quadraticInExpansion = None,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
            HopkinsConductivity = False,
            hmin = 1e-10,
            hmax = 1.0,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            rhoMin = 0.01,
            hourglass = None,
            hourglassOrder = 1,
            hourglassLimiter = 1,
            filter = 0.00,
            KernelConstructor = NBSplineKernel,
            order = 5,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            dtverbose = False,
            steps = None,
            goalTime = 0.15,
            dt = 1e-6,
            dtMin = 1.0e-6,
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 10,
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            volumeType = CRKSumVolume,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            correctVelocityGradient = True,
            gradhCorrection = True,
            linearConsistent = False,

            useRefinement = False,

            clearDirectories = False,
            restoreCycle = -1,
            restartStep = 10000,
            dataDirBase = "dumps-Sod-planar",
            restartBaseName = "Sod-planar-1d-restart",
            outputFile = "None",
            checkRestart = False,

            graphics = True,
            )

assert not(boolReduceViscosity and boolCullenViscosity)

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
                       "numNodeLists=%i" % numNodeLists,
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "CRKVar=%s" % crktype,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "correctionOrder=%s" % correctionOrder,
                       "Cullen=%s" % boolCullenViscosity,
                       "Condc=%s" % HopkinsConductivity,
                       "filter=%f" % filter,
                       "%i" % (nx1 + nx2))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sod-planar-1d-%i" % (nx1 + nx2))

assert numNodeLists in (1, 2)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gammaGas, mu)
strength = NullStrength()

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(order), 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
if solid:
    makeNL = makeSolidNodeList
else:
    makeNL = makeFluidNodeList
nodes1 = makeNL("nodes1", eos, 
                hmin = hmin,
                hmax = hmax,
                nPerh = nPerh,
                kernelExtent = kernelExtent,
                rhoMin = rhoMin)
nodes2 = makeNL("nodes2", eos, 
                hmin = hmin,
                hmax = hmax,
                nPerh = nPerh,
                kernelExtent = kernelExtent,
                rhoMin = rhoMin)
nodeSet = [nodes1, nodes2]

#-------------------------------------------------------------------------------
# Functions to specify the initial density and specific thermal energy with
# optional smoothing.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
hfold = hsmooth*max(dx1, dx2)
def rho_initial(xi):
    if hfold > 0.0:
        return rho1 + (rho2 - rho1)/(1.0 + exp(-(xi - x1)/hfold))
    elif xi <= x1:
        return rho1
    else:
        return rho2

def specificEnergy(xi, rhoi):
    if hfold > 0.0:
        Pi = P1 + (P2 - P1)/(1.0 + exp((-xi - x1)/hfold))
    elif xi <= x1:
        Pi = P1
    else:
        Pi = P2
    return Pi/((gammaGas - 1.0)*rhoi)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from GenerateNodeProfile import GenerateNodeProfile1d
from VoronoiDistributeNodes import distributeNodes1d
if numNodeLists == 1:
    gen = GenerateNodeProfile1d(nx = nx1 + nx2,
                                rho = rho_initial,
                                xmin = x0,
                                xmax = x2,
                                nNodePerh = nPerh)
    distributeNodes1d((nodes1, gen))
else:
    gen1 = GenerateNodeProfile1d(nx = nx1,
                                 rho = rho_initial,
                                 xmin = x0,
                                 xmax = x1,
                                 nNodePerh = nPerh)
    gen2 = GenerateNodeProfile1d(nx = nx2,
                                 rho = rho_initial,
                                 xmin = x1,
                                 xmax = x2,
                                 nNodePerh = nPerh)
    distributeNodes1d((nodes1, gen1),
                      (nodes2, gen2))
output("nodes1.numNodes")
output("nodes2.numNodes")

# Set node specific thermal energies
for nodes in nodeSet:
    pos = nodes.positions()
    eps = nodes.specificThermalEnergy()
    rho = nodes.massDensity()
    for i in xrange(nodes.numInternalNodes):
        eps[i] = specificEnergy(pos[i].x, rho[i])

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    output("db.appendNodeList(nodes)")
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
                 XSVPH = XSPH,
                 linearConsistent = linearConsistent,
                 generateVoid = False,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 xmin = Vector(-100.0),
                 xmax = Vector( 100.0))
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT, 
                   filter = filter,
                   cfl = cfl,
                   correctionOrder = correctionOrder,
                   volumeType = volumeType,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   crktype = crktype)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 correctVelocityGradient = correctVelocityGradient,
                 HopkinsConductivity = HopkinsConductivity)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                correctVelocityGradient = correctVelocityGradient,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")

packages = [hydro]

#-------------------------------------------------------------------------------
# Tweak the artificial viscosity.
#-------------------------------------------------------------------------------
q = hydro.Q
if Cl:
    q.Clinear = Cl
if Cq:
    q.Cquadratic = Cq
if linearInExpansion:
    q.linearInExpansion = linearInExpansion
if quadraticInExpansion:
    q.quadraticInExpansion = quadraticInExpansion
if Qlimiter:
    q.limiter = Qlimiter
if epsilon2:
    q.epsilon2 = epsilon2
if etaCritFrac:
    q.etaCritFrac = etaCritFrac
if etaFoldFrac:
    q.etaFoldFrac = etaFoldFrac
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(q,nh,nh,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(q,WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
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
    hg = hourglass(WT, hourglassOrder, hourglassLimiter)
    output("hg")
    output("hg.order")
    output("hg.limiter")
    packages.append(hg)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0), Vector( 1.0))
xPlane1 = Plane(Vector(x2), Vector(-1.0))
xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)

for p in packages:
    p.appendBoundary(xbc0)
    p.appendBoundary(xbc1)

#-------------------------------------------------------------------------------
# Construct an integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtGrowth = dtGrowth
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
if hourglass:
    output("integrator.havePhysicsPackage(hg)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.rigorousBoundaries")

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
# If requested, reset the intial energies and densities.
#-------------------------------------------------------------------------------
if sumInitialDensity and control.totalSteps == 0:
    packages = integrator.physicsPackages()
    state = State(db, packages)
    derivs = StateDerivatives(db, packages)
    integrator.setGhostNodes()
    integrator.preStepInitialize(state, derivs)
    cm = db.connectivityMap()
    pos = db.fluidPosition
    mass = db.fluidMass
    H = db.fluidHfield
    rho = db.fluidMassDensity
    computeSPHSumMassDensity(cm, WT, True, pos, mass, H, rho)
    for nodes in nodeSet:
        pos = nodes.positions()
        eps = nodes.specificThermalEnergy()
        rho = nodes.massDensity()
        for i in xrange(nodes.numInternalNodes):
            eps[i] = specificEnergy(pos[i].x, rho[i])

#-------------------------------------------------------------------------------
# If we want to use refinement, build the refinemnt algorithm.
#-------------------------------------------------------------------------------
class SelectNodes:
    def __init__(self):
        return
    def prepareForSelection(self, dataBase):
        return
    def selectNodes(self, nodeList):
        if control.totalSteps == 50 and nodeList.name == nodes1.name:
            return range(nodes1.numInternalNodes)
        else:
            return []

class SelectNodes2:
    def __init__(self):
        return
    def prepareForSelection(self, dataBase):
        return
    def selectNodes(self, nodeList):
        if control.totalSteps == 20:
            return range(nodeList.numInternalNodes)
        else:
            return []

if useRefinement:
    from AdaptiveRefinement import *
    select = SelectNodes()
    refine = SplitNodes1d()
    package = AdaptiveRefinement(db, select, refine, control)
    control.appendPeriodicWork(package.refineNodes, 1)

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
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
h1 = 1.0/(nPerh*dx1)
h2 = 1.0/(nPerh*dx2)
answer = SodSolution(nPoints=nx1 + nx2,
                     gamma = gammaGas,
                     rho1 = rho1,
                     P1 = P1,
                     rho2 = rho2,
                     P2 = P2,
                     x0 = x0,
                     x1 = x1,
                     x2 = x2,
                     h1 = 1.0/h1,
                     h2 = 1.0/h2)

#cs = db.newFluidScalarFieldList(0.0, "sound speed")
#db.fluidSoundSpeed(cs)
cs = hydro.soundSpeed()

# Make a flat list from a FieldList
def createList(x):
    result = []
    for i in xrange(len(x)):
        for j in xrange(x[i].numInternalElements):
            result.append(x(i,j))
    return mpi.allreduce(result, mpi.SUM)

# Compute the simulated specific entropy.
rho = createList(db.fluidMassDensity)
P = createList(hydro.pressure())
A = [Pi/rhoi**gammaGas for (Pi, rhoi) in zip(P, rho)]

# The analytic solution for the simulated entropy.
xprof = [x.x for x in createList(db.fluidPosition)]
xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
Aans = [Pi/rhoi**gammaGas for (Pi, rhoi) in zip(Pans,  rhoans)]
csAns = [sqrt(gammaGas*Pi/rhoi) for (Pi, rhoi) in zip(Pans,  rhoans)]

if graphics:
    import Gnuplot
    from SpheralGnuPlotUtilities import *

    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotStyle="lines")
    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    pE = plotEHistory(control.conserve)

    csPlot = plotFieldList(cs, winTitle="Sound speed")
    csAnsData = Gnuplot.Data(xans, csAns, 
                             with_ = "lines",
                             title = "Analytic")
    csPlot.replot(csAnsData)

    Aplot = generateNewGnuPlot()
    Adata = Gnuplot.Data(xprof, A,
                         with_ = "lines",
                         title = "P/rho^\gamma",
                         inline = True)
    AansData = Gnuplot.Data(xprof, Aans,
                         with_ = "lines",
                         title = "Solution",
                         inline = True)
    Aplot.replot(Adata)
    Aplot.replot(AansData)

    plots = [(rhoPlot, "Sod-planar-rho.png"),
             (velPlot, "Sod-planar-vel.png"),
             (epsPlot, "Sod-planar-eps.png"),
             (PPlot, "Sod-planar-P.png"),
             (HPlot, "Sod-planar-h.png"),
             (csPlot, "Sod-planar-cs.png"),
             (Aplot, "Sod-planar-entropy.png")]
    
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
        plots += [(volPlot, "Sod-planar-vol.png"),
                   (aplot, "Sod-planar-ACRK.png"),
                   (bplot, "Sod-planar-BCRK.png"),
                   (splot, "Sod-planar-surfacePoint.png"),
                   (voidplot, "Sod-planar-voidPoint.png")]
    
    viscPlot = plotFieldList(hydro.maxViscousPressure(),
                             winTitle = "max(rho^2 Piij)",
                             colorNodeLists = False)
    plots.append((viscPlot, "Sod-planar-viscosity.png"))
    
    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier(),
                                      winTitle = "Cullen alpha")
        cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt(),
                                       winTitle = "Cullen DalphaDt")
        plots += [(cullAlphaPlot, "Sod-planar-Cullen-alpha.png"),
                  (cullDalphaPlot, "Sod-planar-Cullen-DalphaDt.png")]

    if boolReduceViscosity:
        alphaPlot = plotFieldList(q.ClMultiplier(),
                                  winTitle = "rvAlpha",
                                  colorNodeLists = False)

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.hardcopy(os.path.join(dataDir, filename), terminal="png")

print "Energy conservation: original=%g, final=%g, error=%g" % (control.conserve.EHistory[0],
                                                                control.conserve.EHistory[-1],
                                                                (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0])

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
from SpheralGnuPlotUtilities import multiSort
mof = mortonOrderIndices(db)
mo = createList(mof)
rhoprof = createList(db.fluidMassDensity)
Pprof = createList(hydro.pressure())
vprof = [v.x for v in createList(db.fluidVelocity)]
epsprof = createList(db.fluidSpecificThermalEnergy)
hprof = [1.0/Hi.xx for Hi in createList(db.fluidHfield)]

rmin = x0
rmax = x2
if mpi.rank == 0:
    multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof)
    if outputFile != "None":
        outputFile = os.path.join(dataDir, outputFile)
        f = open(outputFile, "w")
        f.write(("#  " + 19*"'%s' " + "\n") % ("x", "rho", "P", "v", "eps", "A", "h", "mo",
                                               "rhoans", "Pans", "vans", "Aans", "hans",
                                               "x_UU", "rho_UU", "P_UU", "v_UU", "eps_UU", "h_UU"))
        for (xi, rhoi, Pi, vi, epsi, Ai, hi, mi,
             rhoansi, Pansi, vansi, Aansi, hansi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, A, hprof, mo,
                                                         rhoans, Pans, vans, Aans, hans):
            f.write((7*"%16.12e " + "%i " + 5*"%16.12e " + 6*"%i " + '\n') % 
                    (xi, rhoi, Pi, vi, epsi, Ai, hi, mi,
                     rhoansi, Pansi, vansi, Aansi, hansi,
                     unpackElementUL(packElementDouble(xi)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi))))
        f.close()

    import Pnorm
    print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
    failure = False
    hD = []
    for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                              ("Pressure", Pprof, Pans),
                              ("Velocity", vprof, vans),
                              ("Thermal E", epsprof, uans),
                              ("Entropy", A, Aans),
                              ("h       ", hprof, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in xrange(len(data))]
        Pn = Pnorm.Pnorm(error, xprof)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print "\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf)
        #f.write(("\t\t%g") % (L1))
        hD.append([L1,L2,Linf])
    #f.write("\n")

    print "%d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t" % (nx1+nx2,hD[0][0],hD[1][0],hD[2][0],hD[3][0],
                                                                                hD[0][1],hD[1][1],hD[2][1],hD[3][1],
                                                                                hD[0][2],hD[1][2],hD[2][2],hD[3][2])
