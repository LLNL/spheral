#-------------------------------------------------------------------------------
# Diametral Compression Test 2D
#-------------------------------------------------------------------------------
#
# Solid FSISPH
#
#ATS:t100 = test(        SELF, "--clearDirectories True --checkError True --goalTime 5.0 --nrSpecimen 15 ", label="Diametral Compression Test FSISPH -- 2-D", np=1)

from Spheral2d import *
import sys, os
import shutil
from math import *
import mpi
import numpy as np
from SpheralTestUtilities import *
from findLastRestart import *

from CompositeNodeDistribution import *
from GenerateNodeDistribution2d import *
from LatticeSampler import *

# Note -- all units in (cm, gm, microsec)

commandLine(

    compressionSpeed = 1.0e-6,   # downward velocity of top bc (cm/usec)

    # ats settings
    checkError=True,             # check error rel to analytic
    tol = 1.0,                   # error toleration (%) average tensile stress
    leaveNoTrace=True,           # delete output dirs

    # Specimen
    rho0 = 3.24,    # g/cc
    m = 1.73,       # weibull exponent  (Figuero 2016 Tamdakht H-OC)
    k = 2.0e8,      # weibull constant  (fit to Figuero 2016 and Zaytsev 2021)

    # Clamps
    incompressibleClamps = True,
    rClamp = 0.165,

    # Specimen geometry.
    SpecimenDistribution = "conformal",  # node distribution for the rock
    rSpecimen =  0.165,                  # 3.3 mm DIA (Zaytzev 2021)
    nrSpecimen = 20,     	             # number of radial bins
    rotation = 0.0,                      # rotates Specimen distribution (radians)

    # Material parameters 
    eosChoice = "tillotson",       # (gruneisen, tillotson, leos)
    isPressureEquilibrium = False, # initialize eps to give pressure equilibrium
    PoissonsRatio = 0.25,          # possion's ratio
    Y0 = 1000.0 / 1e5,             # yield strength  (Zaytsev 2021)
    Yultimate = 0.372 ,            # ultimate strength  (Zaytsev 2021)
    etaMin = 0.1,                  # rho/rho0 limit
    etaMax = 6.0,                  # rho/rho0 limit

    # Porosity Model
    epsEGranite = 0.0,       # Elastic compaction limit (porosity)
    epsXGranite = -0.4,      # Transition from exponential to power-law distention (porosity)
    kappaGranite = 0.8,      # Exponential factor for distention (porosity)
    
    # Damage Model
    strengthModel = "collins",      # ("collins", "constant", "null")
    fullyDamagedSpecimen=False,     # treat Specimen as full damaged to start?
    useDamage = True,               # turns damage model on
    strengthInDamage = False,       # this isn't ready right now
    strainType = BenzAsphaugStrain, # (BenzAsphaugStrain, PseudoPlasticStrain, MeloshRyanAsphaugStrain, PlasticStrain, PseudoPlasticStrain)

    # Hydro parameters
    cfl = 0.35,   
    SPHType = "FSISPH",                  # (CRKSPH,PSPH,SPH,FSISPH)
    correctVelocityGradient=True,        # linear order velocity gradient for PSPH SPH or FSISPH
    asph = False,                        # turns on elliptic kernels
    xsph = False,                        # updates position based on averaged velocity
    epsilonTensile = 0.00,               # term to fight the tensile instability 
    nTensile = 4,                        # term to fight the tensile instability 
    densityUpdate = IntegrateDensity,    # (IntegrateDensity, RigorousSumDensity) update method for the density
    compatibleEnergyEvolution = True,    # activate secondary "accounting" loop for rigorous energy conservation
    totalEnergyEvolution = False,        # evolve the total energy instead of specific thermal energy
    
    # FSISPH parameters
    fsiSurfaceCoefficient = 0.00,          # adds additional repulsive force to material interfaces)
    fsiRhoStabilizeCoeff = 0.00,           # coefficient that smooths the density field
    fsiEpsDiffuseCoeff = 0.00,             # explicit diiffusion of the thermal energy
    fsiXSPHCoeff = 0.00,                   # fsi uses multiplier for XSPH instead of binary switch
    fsiInterfaceMethod = ModulusInterface, # (HLLCInterface, ModulusInterface)
    
    # CRKSPH parameters
    correctionOrder = LinearOrder,   # for CRKSPH higher order field approximations
    
    # Artificial Viscosity
    Cl = 0.5,                     # Cl linear coefficient for av -- None allows spheral to decide
    Cq = 0.5,                     # Cq quadratic coefficient for av -- None allows spheral to decide
    epsilon2 = 1e-30,              # denominator term in balsara correct to prevent singularity
    balsaraCorrection=False,       # do we want the balsara shear correction
    linearInExpansion = None,      # can run the linear portion when nodes are expanding
    quadraticInExpansion = None,   # same for quadratic term
    etaCritFrac = None,            

    # smoothing scale parameters
    kernelOrder = None,       # for b spline kernels, None defaults to WendlandC2 
    nPerh = 3.51,             # number of neighbors / smoothing length     
    HEvolution =  IntegrateH, # (IdealH , IntegrateH) update method for smoothing kernel
    iterateInitialH = False,  # to calc initial ideal H in controller constructor
     
    # Times, and simulation control.
    steps = None,
    goalTime = 100.0,         # usec -- final time
    dt = 1e-4,                # usec -- initial time step
    dtMin = 1e-4,             # usec -- minimum time step size
    dtMax = 1000.0,           # usec -- maximum time step size
    dtGrowth = 2.0,           # multiplication factor for timestep growth
    dtcheck = True,           # Should the time integrator check the timestep throughout a cycle?
    verbosedt = False,        # will print stats about what controls time-step
    maxSteps = None,          # do we want to limit the total number of cycles
    statsStep = 10,           # 
    redistributeStep = 100,   # cycle for redistribution amongst processors to get more even loading
    restartStep = 5000,       # cycle to write a restart file
    restoreCycle = -1,        # set a cycle to read in an old restart file
    inflowvelfreq = 50,       # How often we measure the Specimen velocity and slow inflow
    historyfreq = 1000,       # how often we write to the history file

    # Output
    vizDerivs = False,                   # output derivatives in viz dump
    vizTime = 10.0,                      # usec -- visit output 
    vizCycle = None,                     # can also output visit silo files based on cycles
    vizName = "DiametralTest",           # name for silo files
    baseDir = "dumps-DiametralTest",     # name for output directory this should be /p/lustre/username/baseName
    clearDirectories = False,            # do we want to delete all old files in baseDir upon restart
)


title("2d diametral compression test.")

#--------------------------------------------------------------------------------
# limits on inputs
#--------------------------------------------------------------------------------

eosChoice = eosChoice.lower()
strengthModel = strengthModel.lower()
SPHType = SPHType.lower()
SpecimenDistribution = SpecimenDistribution.lower()

# valid options for things
assert SPHType in (["crksph","sph","fsisph"])
assert eosChoice in (["gruneisen","tillotson"])
assert strengthModel in (["constant", "collins"]) 
assert SpecimenDistribution in (["lattice","conformal"])

assert not (useDamage and strengthModel=="null") # if were using damage strength must be used

assert fsiXSPHCoeff >= 0.0 
assert fsiEpsDiffuseCoeff >= 0.0 
assert fsiRhoStabilizeCoeff >= 0.0 
assert fsiSurfaceCoefficient >= 0.0 

assert useDamage

#--------------------------------------------------------------------------------
# Dependent imports
#--------------------------------------------------------------------------------
if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d as distributeNodes
else:
    from DistributeNodes import distributeNodes2d as distributeNodes

#--------------------------------------------------------------------------------
# Path and File Name
#--------------------------------------------------------------------------------
hydroname = SPHType
if asph:
    hydroname = "A" + hydroname

# Restart and output files.
dataDir = os.path.join(baseDir,
                       "fullyDamage=%s" % fullyDamagedSpecimen,
                       "nrSpecimen%s" % nrSpecimen,
                       "rot%s" % rotation,
                       hydroname)
restartDir = os.path.join(dataDir, "restarts", "proc-%04i" % mpi.rank)
vizDir = os.path.join(dataDir, "viz")
restartBaseName = os.path.join(restartDir, "Airburst-TestProblem02")

if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()
if not os.path.exists(restartDir):
    os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Our Kernel 
#-------------------------------------------------------------------------------
if kernelOrder is not None:
    WT = TableKernel(NBSplineKernel(kernelOrder), 1000)
else:
    WT = TableKernel(WendlandC2Kernel(), 1000)
output("WT")


#-------------------------------------------------------------------------------
# SiO2 material properties.
#-------------------------------------------------------------------------------
units = CGuS()

if eosChoice == "gruneisen":
    rho0Granite = 2.68               # g/cc
    eosSolidGranite = GruneisenEquationOfState(rho0,                      # ref density (g/cc)
                                               etaMin,      
                                               etaMax,      
                                               C0 = 0.36839,              # cm/usec      
                                               S1 = 1.8954,               # dimensionless
                                               S2 = 0.0,                  # dimensionless
                                               S3 = 0.0,                  # dimensionless
                                               gamma0 = 0.9,              # dimensionless
                                               b = 1.0,                   # dimensionless
                                               atomicWeight = 60.0843,    # atomic weight
                                               constants = units,
                                               minimumPressure = minimumPressure,
                                               )
    T0Granite, eps0Granite = 0.0, 0.0
elif eosChoice == "tillotson":
    eosSolidGranite = TillotsonEquationOfState("granite",
                                               etaMin,
                                               etaMax,
                                               units)
    eosSolidGranite.referenceDensity = rho0
    T0Granite, eps0Granite = 0.0, 0.0 
    rho0Granite=rho0

cs0 = eosSolidGranite.soundSpeed(rho0Granite, eps0Granite)
Ks0 = rho0Granite*cs0*cs0
Gs = 3.0 * Ks0 * (1.0-2.0*PoissonsRatio)/(2.0*(1.0+PoissonsRatio))
Ys = 9.0*Ks0*Gs/(3.0*Ks0+Gs)

print "Reference (rho0, T0, eps0, c0) = (%g, %g, %g, %g)" % (rho0Granite, T0Granite, eps0Granite,cs0)
print "Reference (K, Y, G, nu) = (%g, %g, %g, %g)" % (Ks0, Ys, Gs, PoissonsRatio)

if strengthModel == 'constant':
    strengthSolidGranite = ConstantStrength(Gs,  # mu = 3.2 GPa
                                            Y0)  # Y0 = 15 MPa
elif strengthModel == 'collins':
    shearModulusModel = ConstantStrength(Gs,     # mu = 3.2 GPa
                                         Y0)     # Y0 = 15 MPa
    strengthSolidGranite=CollinsStrength(shearModulusModel, 
                                         2.0,            # intact coeff of friction
                                         0.7,            # damaged coeff of friction
                                         Y0,             # strength at P=0    (15 MPa)
                                         Yultimate);     # strength as P->oo  (1.5GPa)
else:
    strengthSolidGranite = NullStrength()


eosSpecimen = eosSolidGranite
strengthSpecimen = strengthSolidGranite

rho0Specimen = rho0Granite
c0Specimen = cs0

# Clamps default to granite
eosClamps = eosSolidGranite
strengthClamps = strengthSolidGranite
rho0Clamps = rho0Granite
cS0Clamps = cs0
PoissonsRatioClamps = PoissonsRatio

# set up incompressible Clampss as copper with ridic poisson
if incompressibleClamps:
    eosSolidCu = TillotsonEquationOfState("copper",
                                       etaMin,
                                       etaMax,
                                       units)
    PoissonsRatioCu=0.0
    rho0Cu = eosSolidCu.referenceDensity 
    cs0Cu = eosSolidCu.soundSpeed(rho0Cu, 0.0)
    Ks0Cu = rho0Cu*cs0Cu*cs0Cu
    GsCu = Ks0Cu / ((2.0*(1.0+PoissonsRatioCu))/(3.0*(1.0-2.0*PoissonsRatioCu)))
    strengthCu = ConstantStrength(GsCu, 10.0*Y0) 

    print "Reference K (granite, copper) = (%g, %g)" % (Ks0, Ks0Cu)
    print "Reference G (granite, copper) = (%g, %g)" % (Gs, GsCu)
    print "Reference rho (granite, copper) = (%g, %g)" % (rho0Specimen, rho0Cu)
    print "Reference c (granite, copper) = (%g, %g)" % (rho0Specimen, cs0Cu)

    eosClamps = eosSolidCu
    strengthClamps = strengthCu
    rho0Clamps = rho0Cu
    cS0Clamps = cs0Cu
    PoissonsRatioClamps = PoissonsRatioCu
#-------------------------------------------------------------------------------
# Create the NodeLists. 
#-------------------------------------------------------------------------------
boxsize = 2.0*rSpecimen
hmaxSpecimen = boxsize
nodesSpecimen = makeSolidNodeList("Test Specimen",
                          eosSpecimen,
                          strengthSpecimen,
                          xmin = -boxsize * Vector.one,
                          xmax =  boxsize * Vector.one,
                          hmax = hmaxSpecimen,
                          nPerh = nPerh,
                          rhoMin = etaMin*rho0Granite,
                          rhoMax = etaMax*rho0Granite,
                          kernelExtent = WT.kernelExtent)

nodesDriver = makeSolidNodeList("Test Driver",
                                eosClamps,
                                strengthClamps,
                                xmin = -boxsize * Vector.one,
                                xmax =  boxsize * Vector.one,
                                hmax = hmaxSpecimen,
                                nPerh = nPerh,
                                rhoMin = etaMin*rho0Granite,
                                rhoMax = etaMax*rho0Granite,
                                kernelExtent = WT.kernelExtent)

nodesBase = makeSolidNodeList("Test Base",
                              eosClamps,
                              strengthClamps,
                              xmin = -boxsize * Vector.one,
                              xmax =  boxsize * Vector.one,
                              hmax = hmaxSpecimen,
                              nPerh = nPerh,
                              rhoMin = etaMin*rho0Granite,
                              rhoMax = etaMax*rho0Granite,
                              kernelExtent = WT.kernelExtent)


nodeListSet = [nodesSpecimen,nodesDriver,nodesBase]

#-------------------------------------------------------------------------------
# Generate and distribute the nodes. 
#-------------------------------------------------------------------------------

# Driver
#-----------------------------
radRatio = int(rClamp/rSpecimen)
generatorDriver = GenerateNodeDistribution2d(nRadial = nrSpecimen*radRatio, nTheta = nrSpecimen*radRatio,
                                                        rho = rho0Clamps,
                                                        distributionType = "constantDTheta",
                                                        theta = 2.0*pi,
                                                        rotation = pi,
                                                        rmax = rClamp,
                                                        rmin = 0.0,
                                                        offset=[0.0,(1.0+radRatio)*rSpecimen],
                                                        nNodePerh = nPerh)

generatorBase = GenerateNodeDistribution2d(nRadial = nrSpecimen*radRatio, nTheta = nrSpecimen*radRatio,
                                                        rho = rho0Clamps,
                                                        distributionType = "constantDTheta",
                                                        theta = 2.0*pi,
                                                        rotation = 0.0,
                                                        rmax = rClamp,
                                                        rmin = 0.0,
                                                        offset=[0.0,-(1.0+radRatio)*rSpecimen],
                                                        nNodePerh = nPerh)

# Specimen
#-----------------------------
generatorSpecimen = GenerateNodeDistribution2d(nRadial = nrSpecimen, nTheta = nrSpecimen,
                                                        rho = rho0Granite,
                                                        distributionType = "constantDTheta",
                                                        theta = 2.0*pi,
                                                        rotation = rotation,
                                                        rmax = rSpecimen,
                                                        rmin = 0.0,
                                                        nNodePerh = nPerh)

distributeNodes((nodesSpecimen,generatorSpecimen),(nodesDriver,generatorDriver),(nodesBase,generatorBase))


#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase()
for n in nodeListSet:
    db.appendNodeList(n)
    print n.name
del n
nodeLists = db.nodeLists()

output("db")
output("db.nodeLists")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Physics Package : construct the hydro physics object.
#-------------------------------------------------------------------------------
packages=[]

if SPHType=="crksph":
    hydro = CRKSPH(dataBase = db,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = False,
                   densityUpdate = densityUpdate,
                   order = correctionOrder,
                   compatibleEnergyEvolution = compatibleEnergyEvolution,
                   evolveTotalEnergy = totalEnergyEvolution,
                   XSPH = xsph,
                   ASPH = asph)

elif SPHType=="fsisph": 
    q = CRKSPHMonaghanGingoldViscosity(Cl,Cq)   
    hydro = FSISPH(dataBase = db,
                   Q=q,
                   W = WT,
                   cfl = cfl,
                   surfaceForceCoefficient = fsiSurfaceCoefficient,                       
                   densityStabilizationCoefficient = fsiRhoStabilizeCoeff,                 
                   specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,  
                   xsphCoefficient = fsiXSPHCoeff,
                   interfaceMethod = fsiInterfaceMethod,      
                   correctVelocityGradient = correctVelocityGradient,
                   compatibleEnergyEvolution = compatibleEnergyEvolution,  
                   evolveTotalEnergy = totalEnergyEvolution,         
                   ASPH = asph,
                   HUpdate=HEvolution,
                   epsTensile = epsilonTensile,
                   nTensile = nTensile,
                   strengthInDamage=False,
                   damageRelieveRubble=False)
    packages += [hydro.slideSurfaces]

elif SPHType == "sph":
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                useVelocityMagnitudeForDt = False,
                densityUpdate = densityUpdate,
                compatibleEnergyEvolution = compatibleEnergyEvolution,
                evolveTotalEnergy = totalEnergyEvolution,
                XSPH = xsph,
                ASPH = asph,
                HUpdate=HEvolution,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
else:
    raise RuntimeError, "invalid SPHType"

if not epsilon2 is None:
    hydro.Q.epsilon2=epsilon2              
if not Cl is None:
    hydro.Q.Cl = Cl
if not Cq is None:
    hydro.Q.Cq = Cq
if balsaraCorrection:
    hydro.Q.balsaraShearCorrection=True
if not linearInExpansion is None:
    hydro.Q.linearInExpansion = linearInExpansion
if not quadraticInExpansion is None:
    hydro.Q.quadraticInExpansion = quadraticInExpansion
if hasattr(hydro.Q, "etaCritFrac") and not etaCritFrac is None:
    hydro.Q.etaCritFrac = etaCritFrac

output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.densityUpdate")
output("hydro.HEvolution")
output("hydro.XSPH")
if hasattr(hydro, "correctionOrder"):
    output("hydro.correctionOrder")
if hasattr(hydro, "volumeType"):
    output("hydro.volumeType")
output("hydro.Q")
output("hydro.Q.Cl")
output("hydro.Q.epsilon2")
output("hydro.Q.Cq")
output("hydro.Q.balsaraShearCorrection")
output("hydro.Q.linearInExpansion")
output("hydro.Q.quadraticInExpansion")
if hasattr(hydro.Q, "etaCritFrac"):
    output("hydro.Q.etaCritFrac")

packages += [hydro]


#-------------------------------------------------------------------------------
# Physics Package :  construct the damage model.
#-------------------------------------------------------------------------------

damageModelSpecimen = GradyKippTensorDamage("granite",
                                            mWeibull = m,
                                            kWeibull = k,
                                            nodeList = nodesSpecimen,
                                            strainAlgorithm = strainType,
                                            kernel = WT,
                                            units = units)
            
packages.append(damageModelSpecimen)

#-------------------------------------------------------------------------------
# Create initial conditions.
#-------------------------------------------------------------------------------
nodesDriver.velocity(VectorField("initial velocity driver", nodesDriver,Vector(0.0,-compressionSpeed)))
nodesBase.velocity(VectorField("initial velocity base", nodesBase,Vector(0.0,compressionSpeed)))

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.verbose = verbosedt
integrator.allowDtCheck = dtcheck
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.cullGhostNodes = False  #  <--- Fix me!
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.verbose")
output("integrator.allowDtCheck")

#-------------------------------------------------------------------------------
# PeriodicWork function to zap the net x velocity.
#-------------------------------------------------------------------------------
class removeTransverseVelocity:
    def __init__(self,nodes):
        self.nodes = nodes
    def apply(self,cycle,time,dt):
        vel = self.nodes.velocity()
        sumVelx = 0.0
        for i in range(self.nodes.numInternalNodes):
            sumVelx += vel[i].x

        totalSumVelx =  mpi.allreduce(sumVelx,mpi.SUM)
        totalNodes = mpi.allreduce(self.nodes.numInternalNodes,mpi.SUM)
        if totalNodes>0:
            averageVelx = totalSumVelx/totalNodes
            vel -= Vector(averageVelx,0.0)

killVelx = removeTransverseVelocity(nodesSpecimen)
periodicWork = [(killVelx.apply,1)]

#-------------------------------------------------------------------------------
# PeriodicWork function to squeeze the specimen
#-------------------------------------------------------------------------------
class Squeezer:
    def __init__(self,nodes,vNodes,criterion):
        self.nodes = nodes
        self.vNodes = vNodes
        self.criterion = criterion

    def squeeze(self,cycle,time,dt):
        vel = self.nodes.velocity()
        pos = self.nodes.positions()
        for i in range(self.nodes.numInternalNodes):
            if self.criterion(pos[i].y):
                vel[i].y = self.vNodes

def critUpper(yPos):
    crit = False
    if yPos>rSpecimen*(1.0+radRatio):
        crit = True
    return crit

def critLower(yPos):
    crit = False
    if yPos<-rSpecimen*(1.0+radRatio):
        crit = True
    return crit

res = rSpecimen/nrSpecimen
unstoppableForce = Squeezer(nodesDriver,-compressionSpeed,critUpper)
immovableObject = Squeezer(nodesBase,compressionSpeed,critLower)

periodicWork += [(unstoppableForce.squeeze,1)]
periodicWork += [(immovableObject.squeeze,1)]


#-------------------------------------------------------------------------------
# writeOut our strain and oxx
#-------------------------------------------------------------------------------
class loadCurveStorage:

    def __init__(self,
                 vel,
                 dia,
                 db,
                 nodes,
                 pressure):

        self.hist = []
        self.strainHist2 = [0.0]
        self.strainHist = [0.0]
        self.stressHist = [0.0]
        self.timeHist = [0.0]

        self.vel = vel
        self.diameter = dia
        self.delta = 0.0
        self.timeLast = 0.0
        self.sigma = nodes.deviatoricStress()
        self.position = nodes.positions()
        self.pressure = pressure
        self.nodes = nodes
        self.db = db

        self.maxy0 = self.maxY()
        self.miny0 = self.minY()

    def maxY(self):
        maxy = 0.0
        nodeLists = self.db.nodeLists()
        for nodeList in nodeLists:
            positions = nodeList.positions()
            for position in positions:
                maxy = max(maxy,position.y)
        maxy = mpi.allreduce(maxy,mpi.MAX)
        return maxy

    def minY(self):
        miny = 0.0
        nodeLists = self.db.nodeLists()
        for nodeList in nodeLists:
            positions = nodeList.positions()
            for position in positions:
                miny = min(miny,position.y)
        miny = mpi.allreduce(miny,mpi.MIN)
        return miny

    def readSigmaxx(self):
        Syy = []
        for i in range(self.nodes.numInternalNodes):
            if abs(self.position[i].x)<2.0*rSpecimen/nrSpecimen and abs(self.position[i].y)<0.5*rSpecimen:
                Syy.append(float(self.sigma[i].xx-self.pressure[i]))

        Syyreduced = mpi.allreduce(Syy,mpi.SUM)
        sigmaxx = sum(Syyreduced)/len(Syyreduced)
        return sigmaxx

    def storeLoadCurve(self,cycle,time,dt):

        dtNew = time-self.timeLast
        self.delta = 2.0*time*self.vel
        straini = self.delta/self.diameter
        sigmaxx = self.readSigmaxx()

        # other form of strain calc
        maxy = self.maxY()
        miny = self.minY()
        delta2 = (self.maxy0 - self.miny0) - (maxy-miny)
        straini2 = delta2/self.diameter

        self.timeLast = time
        self.timeHist.append(time)
        self.strainHist2.append(100.0*straini2)
        self.strainHist.append(100.0*straini)
        self.stressHist.append(sigmaxx*1e5)
        self.hist.append([time,straini,sigmaxx])


LCS = loadCurveStorage(compressionSpeed, 
                       2.0*rSpecimen, 
                       db,
                       nodesSpecimen,
                       hydro.pressure[2])

periodicWork += [(LCS.storeLoadCurve ,1)]


header_label='"x(km)" "y(km)" "z(km)"  "P (psi)" "Sxx" "Syy" "Sxy"'
P = hydro.pressure
Sigma =db.solidDeviatoricStress
def eulerianSampleVars(j,i):
    
    return (P(j,i),
           Sigma(j,i).xx,
           Sigma(j,i).yy,
           Sigma(j,i).xy)

res_clip = 2.0*rSpecimen/float(nrSpecimen)

#-------------------------------------------------------------------------------
# Sampling function
#-------------------------------------------------------------------------------
samplers = [] # list of sample types (currently only line segments are supported)
LSxx = LineSegment(p0=Vector(0.0,-rSpecimen), p1=Vector(0.0,rSpecimen),
           resolution=rSpecimen/float(nrSpecimen),
           N_Fields = 4,
           histDir=dataDir+'/vertical_',
           header_label=header_label)
samplers.append(LSxx)

LSyy = LineSegment(p0=Vector(-rSpecimen,0.0), p1=Vector(rSpecimen,0.0),
           resolution=rSpecimen/float(nrSpecimen),
           N_Fields = 4,
           histDir=dataDir+'/horizontal_',
           header_label=header_label)
samplers.append(LSyy)

Sample = Sample(db, eulerianSampleVars, samplers )
periodicWork += [(Sample.sample,100)]

#-------------------------------------------------------------------------------
# Print total energy
#-------------------------------------------------------------------------------
def printTotalEnergy(cycle,time,dt):
    Etot=0.0
    mass00=db.solidMass
    vel00=db.solidVelocity
    eps00=db.solidSpecificThermalEnergy
    nodeLists = db.nodeLists()
    for nodelisti in range(db.numNodeLists):

        for i in range(nodeLists[nodelisti].numInternalNodes):
            Etot += mass00(nodelisti,i)*(0.5*vel00(nodelisti,i).magnitude2()+eps00(nodelisti,i))
    Etot = mpi.allreduce(Etot,mpi.SUM)
    print(" TOTAL ENERGY : %.15f" % Etot)

periodicWork += [(printTotalEnergy,5)]
#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            iterateInitialH=iterateInitialH,
                            statsStep = statsStep,
                            restoreCycle = restoreCycle,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartBaseName,
                            vizBaseName = vizName,
                            vizDir = vizDir,
                            vizDerivs = vizDerivs,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            periodicWork=periodicWork)

output("control")

#-------------------------------------------------------------------------------
# Advance to completion.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)


#-------------------------------------------------------------------------------
# Comparison to analytic solution.
#-------------------------------------------------------------------------------

Sigmaxx=LSyy.DATA_Tot[:,5]-LSyy.DATA_Tot[:,3]
x_sample = LSyy.DATA_Tot[:,0]
y_sample = LSyy.DATA_Tot[:,1]
z_sample = LSyy.DATA_Tot[:,2]
xdiff = x_sample[:-1]-x_sample[1:]
force2 = sum((0.5*Sigmaxx[1:]+0.5*Sigmaxx[:-1])*xdiff)
print "Force from Sampling %g " % force2

# analytic eqn
#--------------------------------
Kc = rho0Clamps*cS0Clamps**2
Ks = rho0Specimen*c0Specimen**2
Es = 3.0*Ks*(1-2.0*PoissonsRatio)
Ec = 3.0*Kc*(1-2.0*PoissonsRatioClamps)
oneOverEstar = (1.0-PoissonsRatio**2)/Es + (1.0-PoissonsRatioClamps**2)/Ec
Estar = 1.0/oneOverEstar

delta = compressionSpeed*goalTime

force = pi/4.0*Estar*delta
Pmax = sqrt(force*Estar/pi/rSpecimen)
print "Force from Displacement: %g" % force
print "Pstar: %g" % Pmax

def analyticSigmaxx(X,Y):

    sigma0 = 2.0*force2/(pi)
    sigmax = [0]*len(X)
    R = rSpecimen
    for i in range(len(X)):
        sigmax[i] = sigma0*(0.5/R - ( (X[i]*X[i]*(R-Y[i]))/(X[i]*X[i]+(R-Y[i])**2)**2.0 + 
                                      (X[i]*X[i]*(R+Y[i]))/(X[i]*X[i]+(R+Y[i])**2)**2.0 ))
    return sigmax

def analyticSigmayy(X,Y):

    sigma0 = 2.0*force2/(pi)
    sigmay = [0]*len(X)
    R = rSpecimen
    for i in range(len(X)):
        sigmay[i] = sigma0*(0.5/R - ( ((R-Y[i]))**3.0/(X[i]*X[i]+(R-Y[i])**2)**2.0 + 
                                      ((R+Y[i]))**3.0/(X[i]*X[i]+(R+Y[i])**2)**2.0 ))
    return sigmay

def analyticSigmaxy(X,Y):

    sigma0 = 2.0*force2/(pi)
    sigmay = [0]*len(X)
    R = rSpecimen
    for i in range(len(X)):
        sigmay[i] = sigma0*( ( (X[i]*(R-Y[i]))**2.0/(X[i]*X[i]+(R-Y[i])**2) + 
                               (X[i]*(R+Y[i]))**2.0/(X[i]*X[i]+(R+Y[i])**2) ))
    return sigmay

# numerical Solution
#--------------------------------
P = []
Sxx = []
Syy = []
Sxy = []
Y = []
X = []
pressure = hydro.pressure
sigma = db.solidDeviatoricStress
pos = db.globalPosition

for i in range(nodesSpecimen.numInternalNodes):
    if abs(pos(2,i).x)<2.0*rSpecimen/nrSpecimen and abs(pos(2,i).y)<0.75*rSpecimen:
        P.append(float(pressure(2,i)))
        Sxx.append(float(sigma(2,i).xx-P[-1]))
        Syy.append(float(sigma(2,i).yy-P[-1]))
        Sxy.append(float(sigma(2,i).xy))
        X.append(float(pos(2,i).x))
        Y.append(float(pos(2,i).y))
    

Preduced = mpi.allreduce(P,mpi.SUM)
Sxxreduced = mpi.allreduce(Sxx,mpi.SUM)
Sxyreduced = mpi.allreduce(Sxy,mpi.SUM)
Syyreduced = mpi.allreduce(Syy,mpi.SUM)
Yreduced = mpi.allreduce(Y,mpi.SUM)
Xreduced = mpi.allreduce(X,mpi.SUM)
Sxxanalytic = analyticSigmaxx(Xreduced,Yreduced)
Syyanalytic = analyticSigmayy(Xreduced,Yreduced)

if mpi.rank==1:

    if checkError:

        avgSxxanalytic = sum(Sxxanalytic)/len(Sxxanalytic)
        avgSxxreduced = sum(Sxxreduced)/len(Sxxreduced)
        error = 100.0*abs(avgSxxanalytic - avgSxxreduced)/max(abs(avgSxxanalytic),1e-30)
        
        if error > tol:
            raise ValueError, "tensile stress error bounds violated (error, error tolerance) = (%g,%g)." % (error,tol)

    if leaveNoTrace:
        os.system("rm -rf "+baseDir)

# if mpi.rank==1:
#     import matplotlib.pyplot as plt

#     fig, ax = plt.subplots()



#     ax.plot(Yreduced,Sxxreduced,'bo',label='$\sigma_{xx}$ -- FSISPH',markerfacecolor='c')
#     ax.plot(Yreduced,Sxxanalytic,'b.',label='$\sigma_{xx}$ -- Analytic')
#     ax.plot(Yreduced,Syyreduced,'rs',label='$\sigma_{yy}$ -- FSISPH',markerfacecolor='m')
#     ax.plot(Yreduced,Syyanalytic,'r.',label='$\sigma_{yy}$ -- Analytic')
#     ax.legend()
#     plt.xlabel('y',fontsize=12)
#     plt.ylabel('Cauchy Stress Eigen Values',fontsize=12)
#     plt.show()


