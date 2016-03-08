#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.0 --goalTime 3.0 --nx1=64 --ny1=64 --outputFile='yee.txt'", label="Yee CRK, 64x64", np=20)
#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.0 --goalTime 3.0 --nx1=128 --ny1=128 --outputFile='yee.txt'", label="Yee CRK, 128x128", np=20)
#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.0 --goalTime 3.0 --nx1=256 --ny1=256 --outputFile='yee.txt'", label="Yee CRK, 256x256", np=40)
#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.0 --goalTime 3.0 --nx1=512 --ny1=512 --outputFile='yee.txt'", label="Yee CRK, 512x512", np=40)
#-------------------------------------------------------------------------------
# The Yee-Vortex Test
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

class YeeDensity:
    def __init__(self,
                 xc,
                 yc,
                 gamma,
                 beta,
                 temp_inf):
        self.xc = xc
        self.yc = yc
        self.gamma = gamma
        self.beta = beta
        self.temp_inf = temp_inf
        return
    def __call__(self, r):
        r2 = (r.x-self.xc)*(r.x-self.xc)+(r.y-self.yc)*(r.y-self.yc) 
        temp = self.temp_inf - (self.gamma-1.0)*self.beta*exp(1.0-r2)/(8.0*self.gamma*pi*pi)
        return pow(temp,1.0/(self.gamma-1.0))

title("2-D integrated hydro test --  Yee-Vortex Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    rho = 1.0,
    gamma = 1.4,

    # Translation
    velTx=0.0,
    velTy=0.0,

    # Geometry of Box
    x0 = -5.0,
    x1 =  5.0,
    y0 = -5.0,
    y1 =  5.0,
   
    #Center of Vortex
    xc=0.0,
    yc=0.0,
    
    #Vortex strength
    beta = 5.0,
    #Tempurature at inf
    temp_inf = 1.0,

    # Resolution and node seeding.
    nx1 = 64,
    ny1 = 64,
    seed = "lattice",

    nPerh = 1.51,

    SVPH = False,
    CRKSPH = False,
    PSPH = False,
    SPH = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
    filter = 0.0,  # For CRKSPH
    KernelConstructor = NBSplineKernel,
    order = 5,
    Qconstructor = MonaghanGingoldViscosity,
    #Qconstructor = TensorMonaghanGingoldViscosity,
    boolReduceViscosity = False,
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
    fcentroidal = 0.0,
    fcellPressure = 0.0,
    Cl = 1.0, 
    Cq = 0.75,
    linearInExpansion = False,
    Qlimiter = False,
    balsaraCorrection = False,
    epsilon2 = 1e-2,
    hmin = 1e-5,
    hmax = 0.5,
    hminratio = 0.1,
    cfl = 0.5,
    XSPH = False,
    epsilonTensile = 0.0,
    nTensile = 8,

    IntegratorConstructor = CheapSynchronousRK2Integrator,
    goalTime = 8.0,
    steps = None,
    vizCycle = 20,
    vizTime = 0.1,
    dt = 0.0001,
    dtMin = 1.0e-5, 
    dtMax = 0.1,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    HUpdate = IdealH,
    domainIndependent = False,
    rigorousBoundaries = False,
    dtverbose = False,

    densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
    compatibleEnergy = True,
    gradhCorrection = False,

    useVoronoiOutput = False,
    clearDirectories = False,
    restoreCycle = None,
    restartStep = 200,
    dataDir = "dumps-yeevortex-xy",
    graphics = True,
    smooth = None,
    outputFile = "None",
    )

assert not(boolReduceViscosity and boolCullenViscosity)
# Decide on our hydro algorithm.
if SVPH:
    if SPH:
        HydroConstructor = SVPHFacetedHydro
    else:
        HydroConstructor = ASVPHFacetedHydro
elif CRKSPH:
    Qconstructor = CRKSPHMonaghanGingoldViscosity
    if SPH:
        HydroConstructor = CRKSPHHydro
    else:
        HydroConstructor = ACRKSPHHydro
elif PSPH:
    if SPH:
        HydroConstructor = PSPHHydro
    else:
        HydroConstructor = APSPHHydro
else:
    if SPH:
        HydroConstructor = SPHHydro
    else:
        HydroConstructor = ASPHHydro

#-------------------------------------------------------------------------------
# Build our directory paths.
#-------------------------------------------------------------------------------
densityUpdateLabel = {IntegrateDensity : "IntegrateDensity",
                      SumDensity : "SumDensity",
                      RigorousSumDensity : "RigorousSumDensity",
                      SumVoronoiCellDensity : "SumVoronoiCellDensity"}
baseDir = os.path.join(dataDir,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       KernelConstructor.__name__,
                       "Cl=%g_Cq=%g" % (Cl, Cq),
                       densityUpdateLabel[densityUpdate],
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "nPerh=%3.1f" % nPerh,
                       "fcentroidal=%f" % max(fcentroidal, filter),
                       "fcellPressure=%f" % fcellPressure,
                       "%ix%i" % (nx1, ny1))
restartDir = os.path.join(baseDir, "restarts")
restartBaseName = os.path.join(restartDir, "yeevortex-xy-%ix%i" % (nx1, ny1))

vizDir = os.path.join(baseDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "yeevortex-xy-%ix%i" % (nx1, ny1)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(baseDir):
        shutil.rmtree(baseDir)
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
mu = 1.0
K = 1.0
#eos = GammaLawGasMKS(gamma, mu)
eos = PolytropicEquationOfStateMKS(K,gamma,mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
else:
  WT = TableKernel(KernelConstructor(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("fluid", eos,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               kernelExtent = kernelExtent,
                               nPerh = nPerh)
output("nodes.name")
output("    nodes.hmin")
output("    nodes.hmax")
output("    nodes.hminratio")
output("    nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    rmin = 0.0
    rmax = sqrt(2.0)*(x1-x0)
    
    if(seed=="latticeCylindrical"):
        rmin = x1-8.0*nPerh/nx1
        rmax = x1-2.0*nPerh/nx1
    
    generator = GenerateNodeDistribution2d(nx1, ny1, rho = YeeDensity(xc,yc,gamma,beta,temp_inf),
                                           distributionType = seed,
                                           xmin = (x0, y0),
                                           xmax = (x1, y1),
                                           #rmin = 0.0,
                                           theta = 2.0*pi,
                                           #rmax = sqrt(2.0)*(x1 - x0),
                                           rmax = rmax,
                                           rmin = rmin,
                                           nNodePerh = nPerh,
                                           SPH = SPH)

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes, generator))
    print nodes.name, ":"
    output("    mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    #Set IC
    vel = nodes.velocity()
    eps = nodes.specificThermalEnergy()
    pos = nodes.positions()
    for i in xrange(nodes.numInternalNodes):
        xi, yi = pos[i]
        xci = (xi-xc)
        yci = (yi-yc)
        r2=xci*xci+yci*yci
        velx = -yci*exp((1.0-r2)*0.5)*beta/(2.0*pi)
        vely =  xci*exp((1.0-r2)*0.5)*beta/(2.0*pi)
        vel[i] = Vector(velx,vely)
        temp = temp_inf - (gamma-1.0)*beta*exp(1.0-r2)/(8.0*gamma*pi*pi)
        eps[i] = pow(temp,gamma/(gamma-1.0))/(gamma-1.0)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
db.appendNodeList(nodes)
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
if SVPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             densityUpdate = densityUpdate,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             generateVoid = False,
                             HUpdate = HUpdate,
                             fcentroidal = fcentroidal,
                             fcellPressure = fcellPressure,
                             xmin = Vector(x0 - (x1 - x0), y0 - (y1 - y0)),
                             xmax = Vector(x1 + (x1 - x0), y3 + (y1 - y0)))
elif CRKSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile)
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
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(q,WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
    packages.append(evolveCullenViscosityMultiplier)


#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0, y0), Vector( 1.0,  0.0))
xPlane1 = Plane(Vector(x1, y0), Vector(-1.0,  0.0))
yPlane0 = Plane(Vector(x0, y0), Vector( 0.0,  1.0))
yPlane1 = Plane(Vector(x0, y1), Vector( 0.0, -1.0))

xbc = PeriodicBoundary(xPlane0, xPlane1)
ybc = PeriodicBoundary(yPlane0, yPlane1)

xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)
ybc0 = ReflectingBoundary(yPlane0)
ybc1 = ReflectingBoundary(yPlane1)

bcSet = [xbc, ybc]
#bcSet = [xbc0, xbc1, ybc0, ybc1]

for p in packages:
    for bc in bcSet:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.cullGhostNodes = False
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
# If requested, smooth the initial conditions.
#-------------------------------------------------------------------------------
if smooth:
    for iter in xrange(smooth):
        db.updateConnectivityMap(False)
        cm = db.connectivityMap()
        position_fl = db.fluidPosition
        weight_fl = db.fluidMass
        H_fl = db.fluidHfield
        m0_fl = db.newFluidScalarFieldList(0.0, "m0")
        m1_fl = db.newFluidVectorFieldList(Vector.zero, "m1")
        m2_fl = db.newFluidSymTensorFieldList(SymTensor.zero, "m2")
        A0_fl = db.newFluidScalarFieldList(0.0, "A0")
        A_fl = db.newFluidScalarFieldList(0.0, "A")
        B_fl = db.newFluidVectorFieldList(Vector.zero, "B")
        C_fl = db.newFluidVectorFieldList(Vector.zero, "C")
        D_fl = db.newFluidTensorFieldList(Tensor.zero, "D")
        gradA0_fl = db.newFluidVectorFieldList(Vector.zero, "gradA0")
        gradA_fl = db.newFluidVectorFieldList(Vector.zero, "gradA")
        gradB_fl = db.newFluidTensorFieldList(Tensor.zero, "gradB")
        computeCRKSPHCorrections(cm, WT, weight_fl, position_fl, H_fl, True,
                               m0_fl, m1_fl, m2_fl,
                               A0_fl, A_fl, B_fl, C_fl, D_fl, gradA0_fl, gradA_fl, gradB_fl)
        eps0 = db.fluidSpecificThermalEnergy
        vel0 = db.fluidVelocity
        eps1 = interpolateCRKSPH(eps0, position_fl, weight_fl, H_fl, True, 
                               A_fl, B_fl, cm, WT)
        vel1 = interpolateCRKSPH(vel0, position_fl, weight_fl, H_fl, True, 
                               A_fl, B_fl, cm, WT)
        eps0.assignFields(eps1)
        vel0.assignFields(vel1)

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
if useVoronoiOutput:
    import SpheralVoronoiSiloDump
    vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
else:
    vizMethod = None # default
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
                            skipInitialPeriodicWork = (HydroConstructor in (SVPHFacetedHydro, ASVPHFacetedHydro)),
                            SPH = SPH)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile != "None":
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralGnuPlotUtilities import multiSort
    P = ScalarField("pressure", nodes)
    nodes.pressure(P)
    xprof = mpi.reduce([x.x for x in nodes.positions().internalValues()], mpi.SUM)
    yprof = mpi.reduce([x.y for x in nodes.positions().internalValues()], mpi.SUM)
    rhoprof = mpi.reduce(nodes.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.magnitude() for v in nodes.velocity().internalValues()], mpi.SUM)
    velx = mpi.reduce([v.x for v in nodes.velocity().internalValues()], mpi.SUM)
    vely = mpi.reduce([v.y for v in nodes.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in nodes.Hfield().internalValues()], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)

    if mpi.rank == 0:
        rprof = [sqrt(xi*xi + yi*yi) for xi, yi in zip(xprof, yprof)]
        multiSort(rprof, mo, xprof, yprof, rhoprof, Pprof, vprof, epsprof, hprof,velx,vely)
        L1rho = 0.0
        L1eps = 0.0
        L1vel = 0.0
        epsans = []
        rhoans = []
        velans = []
        for i in xrange(len(xprof)):
           r = (Vector(xprof[i],yprof[i]) - Vector(xc, yc)).magnitude()
           r2 = r*r
           temp = temp_inf - (gamma-1.0)*beta*exp(1.0-r2)/(8.0*gamma*pi*pi)
           yci = yprof[i] - yc
           xci = xprof[i] - xc
           velxans = -yci*exp((1.0-r2)*0.5)*beta/(2.0*pi)
           velyans =  xci*exp((1.0-r2)*0.5)*beta/(2.0*pi)
           epsans.append(pow(temp,gamma/(gamma-1.0))/(gamma-1.0))
           rhoans.append(pow(temp,1.0/(gamma-1.0)))
           velans.append(Vector(velxans,velyans).magnitude())
           L1rho = L1rho + abs(rhoans[i]-rhoprof[i])
           L1eps = L1eps + abs(epsans[i]-epsprof[i])
           L1vel = L1vel + abs(velans[i]-vprof[i])
        L1rho = L1rho/len(xprof)
        L1eps = L1eps/len(xprof)
        L1vel = L1vel/len(xprof)
        with open("converge-CRK-%s.txt" % CRKSPH, "a") as myfile:
          myfile.write("%s\t %s\t %s\t %s\n" % (nx1,L1rho,L1eps,L1vel))
        f = open(outputFile, "w")
        f.write(("# " + 19*"%15s " + "\n") % ("r", "x", "y", "rho", "P", "v", "eps", "h", "mortonOrder", "rhoans", "epsans", "velans",
                                              "x_uu", "y_uu", "rho_uu", "P_uu", "v_uu", "eps_uu", "h_uu"))
        for (ri, xi, yi, rhoi, Pi, vi, epsi, hi, mi,rhoa,epsa,vela)  in zip(rprof, xprof, yprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,rhoans,epsans,velans):
            f.write((8*"%16.12e " + "%i " + 3*"%16.12e " + 7*"%i " + "\n") % (ri, xi, yi, rhoi, Pi, vi, epsi, hi, mi, rhoa,epsa,vela,
                                                               unpackElementUL(packElementDouble(xi)),
                                                               unpackElementUL(packElementDouble(yi)),
                                                               unpackElementUL(packElementDouble(rhoi)),
                                                               unpackElementUL(packElementDouble(Pi)),
                                                               unpackElementUL(packElementDouble(vi)),
                                                               unpackElementUL(packElementDouble(epsi)),
                                                               unpackElementUL(packElementDouble(hi))))
        f.close()
