#-------------------------------------------------------------------------------
# A tube of stainless steel undergoing expansion due to a projectile entering
# and smacking into a stop.
#
# This is the mock 3-D cylindrical geometry "RZ" version.
#
# See Vogler et al., 2003, International Journal of Impact Engineering, 29, 735
#-------------------------------------------------------------------------------
from Numeric import *
from SolidSpheral import *
from SpheralOptionParser import commandLine
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralController import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from math import *

from bevelTubeEntrance import *

import sys
sys.path.append("../Utilities")
from NodeHistory import NodeHistory
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("RZ expanding tube impact strength/damage model test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(
    SolidNodeListConstructor = SphSolidNodeList3d,

    # Geometry
    tubeThickness = 0.3,
    rtubeInner = 0.5*1.27,
    ltube = 5.08,
    lCuAnvil = 1.0,
    lFoamAnvil = 1.0,
    lSteelAnvil = 1.0,

    # Numbers of nodes.
    nrtube = 20,
    nltube = 338,
    nrAnvil = 100,
    nrSteelAnvilCap = 5,

    # VISAR sampling parameters.
    dxVISAR = 0.04,
    dyVISAR = 0.04,
    VISARsampleFrequency = 10,

    # Inital velocity of the projectile.
    vxproj = -1.92e5,

    # Parameters for the damage model.
    DamageModelConstructor = GradyKippTensorDamage3d,
    kWeibullSteelFactor = 1.0,
    mWeibullSteelFactor = 1.0,
    randomSeedSteel = 109482993,
    strainType = TensorDamageModel3d.TensorStrainAlgorithm.StrainHistory,

    # Node seeding stuff.
    nPerh = 2.01,

    # Material specific bounds on the mass density.
    etaMinSteel = 0.6,
    etaMaxSteel = 1.4,
    etaMinLexan = 0.5,
    etaMaxLexan = 1.5,
    etaMinCu = 0.5,
    etaMaxCu = 1.5,
    etaMinFoam = 0.5,
    etaMaxFoam = 1.5,

    # Hydro parameters.
    Qconstructor = TensorMonaghanGingoldViscosity3d,
    Cl = 1.0,
    Cq = 1.0,
    Qlimiter = True,
    balsaraCorrection = False,
    epsilon2 = 1e-2,
    hmin = 1.0e-5,
    hmax = 0.5,
    hminratio = 0.1,
    cfl = 0.5,
    XSPH = True,
    epsilonTensile = 0.0,
    nTensile = 4,
    HEvolution = Hydro3d.HEvolutionType.IdealH,
    sumForMassDensity = Hydro3d.MassDensityType.IntegrateDensity,
    compatibleEnergyEvolution = True,

    # Times, and simulation control.
    goalTime = 50.0e-6,
    dtSample = 50.0e-6 / 200.0,
    dt = 1e-10,
    dtMin = 1e-10,
    dtMax = 1e-3,
    dtGrowth = 2.0,
    maxSteps = 200,
    statsStep = 10,
    redistributeStep = None,
    smoothIters = 0,

    # Restart and output files.
    restoreCycle = None,
    restartStep = 200,
    baseDir = "dumps-expandingTube-rz",
    )

# Derived geometry.
rtubeOuter = rtubeInner + tubeThickness
rplug, lplug = rtubeInner, 0.5*ltube
rproj, lproj = rplug, lplug
rAnvil = 2.0*rtubeOuter
lAnvil = lCuAnvil + lFoamAnvil + lSteelAnvil

# Use the above geometry to define enclosing points of the materials for the
# node generators.
xminTube = (lAnvil, rtubeInner)
xmaxTube = (lAnvil + ltube, rtubeOuter)

xminPlug = (lAnvil, 0.0)
xmaxPlug = (lAnvil + lplug, rplug)

xminProj = (lAnvil + ltube, 0.0)
xmaxProj = (lAnvil + ltube + lproj, rproj)

xminSteelAnvil = (0.0, 0.0)
xmaxSteelAnvil = (lSteelAnvil, rAnvil)

xminFoamAnvil = (lSteelAnvil, 0.0)
xmaxFoamAnvil = (lSteelAnvil + lFoamAnvil, rAnvil)

xminCuAnvil = (lSteelAnvil + lFoamAnvil, 0.0)
xmaxCuAnvil = (lAnvil, rAnvil)

xminSteelAnvilCap = (0.0, rAnvil)
xmaxSteelAnvilCap = (lAnvil, rAnvil + float(nrSteelAnvilCap)/float(nrAnvil)*rAnvil)

# The geometry of the bevel on the inner tube opening surface.
tubeOpeningAngle = 1.8 * pi/180.0 # radians
xBevelBegin = lAnvil + ltube - 0.6

# Define the VISAR sampling points.
xVISARa = lAnvil + 2.5
xVISARb = lAnvil + 2.0
xVISARc = lAnvil + 1.5

yVISARa = rtubeOuter
yVISARb = rtubeOuter
yVISARc = rtubeOuter

# Derived numbers of nodes.
nrplug, nlplug = int(rplug/lplug * nltube/4), nltube/4
nrproj, nlproj = nrplug, nlplug

nlAnvil = int(nrAnvil * lAnvil/rAnvil)
nlSteelAnvil = int(nlAnvil * 0.5*lSteelAnvil/lAnvil + 0.5)
nlFoamAnvil = int(nlAnvil * 0.25*lFoamAnvil/lAnvil + 0.5)
nlCuAnvil = int(nlAnvil * 0.5*lCuAnvil/lAnvil + 0.5)
nrSteelAnvil = int(nlSteelAnvil * rAnvil/lSteelAnvil + 0.5)
nrFoamAnvil = int(nlFoamAnvil * rAnvil/lFoamAnvil + 0.5)
nrCuAnvil = int(nlCuAnvil * rAnvil/lCuAnvil + 0.5)
nlSteelAnvilCap = nlSteelAnvil + nlFoamAnvil + nlCuAnvil

# Mass densities.
rho0Steel = 7.85
rho0Lexan = 1.196
rho0Air = 0.001205
rho0Foam = 1.046
rho0Cu = 8.93

# Inital velocity of the projectile.
v0proj = Vector3d(vxproj, 0.0, 0.0)

# Damage model parameters.
kWeibullSteel = 8.8e4 * kWeibullSteelFactor
mWeibullSteel = 2.63  * mWeibullSteelFactor
volumeSteel = pi*(rtubeOuter**2 - rtubeInner**2)*ltube

# Restart and output files.
dataDir = "%s/%s/%s/k=%4.2f_m=%4.2f" % (baseDir,
                                        str(DamageModelConstructor).split("'")[1],
                                        str(SolidNodeListConstructor).split("'")[1],
                                        kWeibullSteel,
                                        mWeibullSteel)
restartDir = dataDir + "/restarts/proc-%04i" % mpi.rank
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/ExpandingTube"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()
if not os.path.exists(restartDir):
    os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Stainless steel material properties.
#-------------------------------------------------------------------------------
eosSteel = LinearPolynomialEquationOfStateCGS3d(rho0Steel,    # reference density  
                                                etaMinSteel,  # etamin             
                                                etaMaxSteel,  # etamax             
                                                0.0,          # A0
                                                1.649901e12,  # A1
                                                1.674656e12,  # A2
                                                0.832543e12,  # A3
                                                1.93,         # B1
                                                0.5,          # B2
                                                0.0,          # B3
                                                55.350)       # atomic weight

coldFitSteel = NinthOrderPolynomialFit(-1.06797724e10,
                                    -2.06872020e10,
                                    8.24893246e11,
                                    -2.39505843e10,
                                    -2.44522017e10,
                                    5.38030101e10,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0)

meltFitSteel = NinthOrderPolynomialFit(7.40464217e10,
                                    2.49802214e11,
                                    1.00445029e12,
                                    -1.36451475e11,
                                    7.72897829e9,
                                    5.06390305e10,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0)

strengthModelSteel = SteinbergGuinanStrengthCGS3d(eosSteel,
                                                  7.700000e11,        # G0
                                                  2.2600e-12,         # A
                                                  4.5500e-04,          # B
                                                  3.4000e9,           # Y0
                                                  2.5e10,             # Ymax
                                                  1.0e-3,             # Yp
                                                  43.0000,            # beta
                                                  0.0,                # gamma0
                                                  0.35,               # nhard
                                                  coldFitSteel,
                                                  meltFitSteel)

#-------------------------------------------------------------------------------
# Lexan material properties.
# Note for lack of strength information about this material, I'm subsituting in
# the strength paramters for lucite here.  :)
#-------------------------------------------------------------------------------
eosLexan = GruneisenEquationOfStateCGS3d(rho0Lexan,    # reference density  
                                         etaMinLexan,  # etamin             
                                         etaMaxLexan,  # etamax             
                                         0.1933e6,     # C0                 
                                         3.49,         # S1                 
                                        -8.2,          # S2                 
                                         9.6,          # S3                 
                                         0.61,         # gamma0             
                                         0.0,          # b                  
                                         28423.0)      # atomic weight

coldFitLexan = NinthOrderPolynomialFit(-5.19191852e9,
                                       -4.41500192e9,
                                       2.84720528e10,
                                       2.14093899e10,
                                       -4.46412259e9,
                                       1.24495222e9,
                                       0.0,
                                       0.0,
                                       0.0,
                                       0.0)

meltFitLexan = NinthOrderPolynomialFit(5.24383771e8,
                                       1.49188457e9,
                                       2.85704428e10,
                                       2.13783662e10,
                                       -4.45135120e9,
                                       1.24138074e9,
                                       0.0,
                                       0.0,
                                       0.0,
                                       0.0)

strengthLexan = SteinbergGuinanStrengthCGS3d(eosLexan,
                                             2.320000e10,        # G0
                                             0.0,                # A
                                             0.0,                # B
                                             4.2000e9,           # Y0
                                             1.0e12,             # Ymax
                                             1.0e-3,             # Yp
                                             0.0,                # beta
                                             0.0,                # gamma0
                                             0.0,                # nhard
                                             coldFitLexan,
                                             meltFitLexan)

#-------------------------------------------------------------------------------
# Copper material properties.
#-------------------------------------------------------------------------------
eosCu = GruneisenEquationOfStateCGS3d(rho0Cu,      # reference density  
                                      etaMinCu,    # etamin             
                                      etaMaxCu,    # etamax             
                                      0.394e6,     # C0                 
                                      1.489,       # S1                 
                                      0.0,         # S2                 
                                      0.0,         # S3                 
                                      2.02,        # gamma0             
                                      0.47,        # b                  
                                      63.57)       # atomic weight

coldFitCu = NinthOrderPolynomialFit(-1.05111874e10,
                                    -2.13429672e10,
                                    6.92768584e11,
                                    -2.45626513e10,
                                    -2.48677403e10,
                                    4.35373677e10,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0)

meltFitCu = NinthOrderPolynomialFit(5.22055639e10,
                                    1.90143176e11,
                                    8.51351901e11,
                                    -1.12049022e11,
                                    -6.11436674e9,
                                    4.36007831e10,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0)

strengthModelCu = SteinbergGuinanStrengthCGS3d(eosCu,
                                               4.770000e11,        # G0
                                               2.8300e-12,         # A
                                               3.7700e-04,         # B
                                               1.2000e9,           # Y0
                                               6.4000e9,           # Ymax
                                               1.0e-3,             # Yp
                                               36.0000,            # beta
                                               0.0,                # gamma0
                                               0.45,               # nhard
                                               coldFitCu,
                                               meltFitCu)

#-------------------------------------------------------------------------------
# Foam material properties. (Polystyrene CH)
#-------------------------------------------------------------------------------
eosFoam = GruneisenEquationOfStateCGS3d(rho0Foam,   # reference density  
                                        etaMinFoam, # etamin             
                                        etaMaxFoam, # etamax             
                                        0.189e6,    # C0                 
                                        2.965,      # S1                 
                                       -4.069,      # S2                 
                                        2.328,      # S3                 
                                        0.67,       # gamma0             
                                        0.0,        # b                  
                                        6.982)      # atomic weight

#-------------------------------------------------------------------------------
# Air material properties.
#-------------------------------------------------------------------------------
gammaAir = 1.4
molecularWeightAir = 30.0
eosAir = GammaLawGasCGS3d(gammaAir, molecularWeightAir)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 1000)
WTPi = TableKernel3d(BSplineKernel3d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodesSteel = SolidNodeListConstructor("Stainless steel", eosSteel, strengthModelSteel, WT, WTPi)
nodesPlug = SolidNodeListConstructor("Lexan plug", eosLexan, strengthLexan, WT, WTPi)
nodesProj = SolidNodeListConstructor("Lexan projectile", eosLexan, strengthLexan, WT, WTPi)
nodesSteelAnvil = SolidNodeListConstructor("Anvil (Steel)", eosSteel, strengthModelSteel, WT, WTPi)
nodesFoamAnvil = SolidNodeListConstructor("Anvil (Foam)", eosFoam, strengthLexan, WT, WTPi)
nodesCuAnvil = SolidNodeListConstructor("Anvil (Copper)", eosCu, strengthModelCu, WT, WTPi)
nodeSet = [nodesSteel,
           nodesPlug,
           nodesProj,
           nodesSteelAnvil,
           nodesFoamAnvil,
           nodesCuAnvil]
for n, rho0, etaMin, etaMax in zip(nodeSet,
                                   [rho0Steel, rho0Lexan, rho0Lexan, rho0Steel, rho0Lexan, rho0Cu],
                                   [etaMinSteel, etaMinLexan, etaMinLexan, etaMinSteel, etaMinLexan, etaMinCu],
                                   [etaMaxSteel, etaMaxLexan, etaMaxLexan, etaMaxSteel, etaMaxLexan, etaMaxCu]):
    n.nodesPerSmoothingScale = nPerh
    n.epsilonTensile = epsilonTensile
    n.nTensile = nTensile
    n.hmin = hmin
    n.hmax = hmax
    n.hminratio = hminratio
    n.XSPH = XSPH
    n.rhoMin = etaMin*rho0
    n.rhoMax = etaMax*rho0
    output("n.name()")
    output("  n.nodesPerSmoothingScale")
    output("  n.epsilonTensile")
    output("  n.nTensile")
    output("  n.hmin")
    output("  n.hmax")
    output("  n.hminratio")
    output("  n.XSPH")
    output("  n.rhoMin")
    output("  n.rhoMax")
del n

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
cache = []
for n in nodeSet:
    neighbor = TreeNeighbor3d(n,
                              kernelExtent = kernelExtent)
    n.registerNeighbor(neighbor)
    cache.append(neighbor)
del n

#-------------------------------------------------------------------------------
# Set node properties (positions, velocites, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    print "Generating node distribution."
    from GenerateNodeDistribution2d import *
    from CompositeNodeDistribution import *
    from ParMETISDistributeNodes import distributeNodes3d
    generatorTube = GenerateNodeDistributionRZ(nltube,
                                               nrtube,
                                               rho0Steel,
                                               "lattice",
                                               xmin = xminTube,
                                               xmax = xmaxTube,
                                               nNodePerh = nPerh,
                                               SPH = not isinstance(nodesSteel, AsphSolidNodeList3d))
    generatorPlug = GenerateNodeDistributionRZ(nlplug,
                                               nrplug,
                                               rho0Lexan,
                                               "lattice",
                                               xmin = xminPlug,
                                               xmax = xmaxPlug,
                                               nNodePerh = nPerh,
                                               SPH = not isinstance(nodesPlug, AsphSolidNodeList3d))
    generatorProj = GenerateNodeDistributionRZ(nlproj,
                                               nrproj,
                                               rho0Lexan,
                                               "lattice",
                                               xmin = xminProj,
                                               xmax = xmaxProj,
                                               nNodePerh = nPerh,
                                               SPH = not isinstance(nodesProj, AsphSolidNodeList3d))
    generatorSteelAnvil1 = GenerateNodeDistributionRZ(nlSteelAnvil,
                                                      nrSteelAnvil,
                                                      rho0Steel,
                                                      "lattice",
                                                      xmin = xminSteelAnvil,
                                                      xmax = xmaxSteelAnvil,
                                                      nNodePerh = nPerh,
                                                      SPH = not isinstance(nodesSteelAnvil, AsphSolidNodeList3d))
    generatorSteelAnvil2 = GenerateNodeDistributionRZ(nlSteelAnvilCap,
                                                      nrSteelAnvilCap,
                                                      rho0Steel,
                                                      "lattice",
                                                      xmin = xminSteelAnvilCap,
                                                      xmax = xmaxSteelAnvilCap,
                                                      nNodePerh = nPerh,
                                                      SPH = not isinstance(nodesSteelAnvil, AsphSolidNodeList3d))
    generatorSteelAnvil = CompositeNodeDistribution(generatorSteelAnvil1, generatorSteelAnvil2)
    generatorFoamAnvil = GenerateNodeDistributionRZ(nlFoamAnvil,
                                                    nrFoamAnvil,
                                                    rho0Foam,
                                                    "lattice",
                                                    xmin = xminFoamAnvil,
                                                    xmax = xmaxFoamAnvil,
                                                    nNodePerh = nPerh,
                                                    SPH = not isinstance(nodesFoamAnvil, AsphSolidNodeList3d))
    generatorCuAnvil = GenerateNodeDistributionRZ(nlCuAnvil,
                                                  nrCuAnvil,
                                                  rho0Cu,
                                                  "lattice",
                                                  xmin = xminCuAnvil,
                                                  xmax = xmaxCuAnvil,
                                                  nNodePerh = nPerh,
                                                  SPH = not isinstance(nodesCuAnvil, AsphSolidNodeList3d))

    print "Starting node distribution..."
    distributeNodes3d((nodesSteel, generatorTube),
                      (nodesPlug, generatorPlug),
                      (nodesProj, generatorProj),
                      (nodesSteelAnvil, generatorSteelAnvil),
                      (nodesFoamAnvil, generatorFoamAnvil),
                      (nodesCuAnvil, generatorCuAnvil))
    nGlobalNodes = 0
    for n in nodeSet:
        print "Generator info for %s" % n.name()
        output("    mpi.allreduce(n.numInternalNodes, mpi.MIN)")
        output("    mpi.allreduce(n.numInternalNodes, mpi.MAX)")
        output("    mpi.allreduce(n.numInternalNodes, mpi.SUM)")
        nGlobalNodes += mpi.allreduce(n.numInternalNodes, mpi.SUM)
    del n
    print "Total number of (internal) nodes in simulation: ", nGlobalNodes

    # Bevel the inner opening surface of the target tube.
    numNodesBeveled = bevelTubeEntrance(nodesSteel,
                                        2,
                                        tubeOpeningAngle,
                                        rtubeInner,
                                        tubeThickness,
                                        xBevelBegin)
    print "Beveled %i nodes in the tube opening." % mpi.allreduce(numNodesBeveled,
                                                                  mpi.SUM)

    # Adjust the diameter of the projectile inward a bit, so it will slide
    # into the tube properly.
    drProj = 0.75*nPerh*rproj/nrproj
    projMultiplier = (rproj - drProj)/rproj
    for i in xrange(nodesProj.numInternalNodes):
        nodesProj.positions()[i].y *= projMultiplier

    # Adjust the plug to match.
    for i in xrange(nodesPlug.numInternalNodes):
        nodesPlug.positions()[i].y *= projMultiplier

    # Iterate over the NodeLists and set some initial conditions.
    for n, rho0 in [(nodesSteel, rho0Steel),
                    (nodesPlug, rho0Lexan),
                    (nodesProj, rho0Lexan),
                    (nodesSteelAnvil, rho0Steel),
                    (nodesFoamAnvil, rho0Foam),
                    (nodesCuAnvil, rho0Cu)]:

        # Set node specific thermal energies
        u0 = n.equationOfState().specificThermalEnergy(rho0, 300.0)
        n.specificThermalEnergy(ScalarField3d("tmp", n, u0))
        print "Initial pressure for %s: %g" % (n.name(),
                                               n.equationOfState().pressure(rho0, u0))
    del n

    # Set the projectile velocities.
    nodesProj.velocity(VectorField3d("tmp", nodesProj, v0proj))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase3d()
for n in nodeSet:
    db.appendNodeList(n)
del n
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro3d(WT, WTPi, q, compatibleEnergyEvolution)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = True
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct a strength physics object.
#-------------------------------------------------------------------------------
strength = Strength3d()
output("strength")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
nfull = max(1, (2.0*pi*(rtubeInner + 0.5*tubeThickness)/(tubeThickness/nrtube) *
                mpi.allreduce(nodesSteel.numInternalNodes, mpi.SUM)))
nflaws = int(nfull*log(nfull))
print "Computing equivalent 3-D number of nodes in tube: %i" % nfull
print "Resulting in effective total number of flaws in volume: %i" % nflaws
damageModel = DamageModelConstructor(nodesSteel,
                                     kWeibullSteel,
                                     mWeibullSteel,
                                     volumeSteel,
                                     WT,
                                     randomSeedSteel,
                                     strainType,
                                     0.4,
                                     1,
                                     nflaws)
output("damageModel")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator3d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.appendPhysicsPackage(damageModel)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(strength)")
output("integrator.havePhysicsPackage(damageModel)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Construct boundary conditions, and add them to our physics packages.
#-------------------------------------------------------------------------------
xbcPlane = Plane3d(Vector3d(0.0, 0.0), Vector3d(1.0, 0.0))
xbc = ReflectingBoundary3d(xbcPlane)
rzbc = CylindricalBoundary(db)

for package in integrator.physicsPackages():
    package.appendBoundary(rzbc)
    package.appendBoundary(xbc)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)
output("control")

#-------------------------------------------------------------------------------
# Select the nodes for the VISAR sampling.
#-------------------------------------------------------------------------------
class SelectVISARNodes:
    def __init__(self, x0, dx, dy, label):
        self.x0 = x0
        self.dx = dx
        self.dy = dy
        self.label = label
        return
    def __call__(self, nodes):
        r = nodes.positions()
        potentials = [i for i in xrange(nodes.numInternalNodes)
                      if abs(r[i].x - self.x0) < self.dx]
        ymax = mpi.allreduce(max([r[i].y for i in potentials] + [-1e50]), mpi.MAX)
        result = [i for i in potentials if r[i].y > ymax - self.dy]
        print "Selected %i %s velocimetry test points." % (mpi.allreduce(len(result), mpi.SUM),
                                                           self.label)
        return result

#-------------------------------------------------------------------------------
# Sampling function to measure the average velocity at the VISAR probe sites.
#-------------------------------------------------------------------------------
def averageCylindricalRadialVelocity(nodes, indicies):
    m = nodes.mass()
    v = nodes.velocity()

    massSum = 1e-30
    result = 0.0
    for i in indicies:
        assert i >= 0 and i < nodes.numInternalNodes
        massSum += m[i]
        result += m[i] * v[i].y

    globalMassSum = mpi.allreduce(massSum, mpi.SUM)
    globalResult = mpi.allreduce(result, mpi.SUM)
    assert globalMassSum > 0.0
    globalResult /= globalMassSum
    return globalResult

#-------------------------------------------------------------------------------
# Build the history objects to simulate the VISAR velocity probes.
#-------------------------------------------------------------------------------
nodesA = SelectVISARNodes(xVISARa, dxVISAR, dyVISAR, "A")
nodesB = SelectVISARNodes(xVISARb, dxVISAR, dyVISAR, "B")
nodesC = SelectVISARNodes(xVISARc, dxVISAR, dyVISAR, "C")
VISARa = NodeHistory(nodesSteel, nodesA, averageCylindricalRadialVelocity,
                     dataDir + "/VISAR-a")
VISARb = NodeHistory(nodesSteel, nodesB, averageCylindricalRadialVelocity,
                     dataDir + "/VISAR-b")
VISARc = NodeHistory(nodesSteel, nodesC, averageCylindricalRadialVelocity,
                     dataDir + "/VISAR-c")
VISARa.nodeFlags.setName("VISAR a points")
VISARb.nodeFlags.setName("VISAR b points")
VISARc.nodeFlags.setName("VISAR c points")

control.appendPeriodicWork(VISARa.sample, VISARsampleFrequency)
control.appendPeriodicWork(VISARb.sample, VISARsampleFrequency)
control.appendPeriodicWork(VISARc.sample, VISARsampleFrequency)

#-------------------------------------------------------------------------------
# Drop visualization files.
#-------------------------------------------------------------------------------
def viz(fields = [],
        filename = "ExpandingTube-rz"):
    damage = damageModel.damage()
    dtrace = ScalarField3d("damage magnitude", nodesSteel)
    dmin = ScalarField3d("damage min", nodesSteel)
    dmax = ScalarField3d("damage max", nodesSteel)
    strain = damageModel.strain()
    svol = ScalarField3d("strain vol", nodesSteel)
    smin = ScalarField3d("strain min", nodesSteel)
    smax = ScalarField3d("strain max", nodesSteel)
    for i in xrange(nodesSteel.numInternalNodes):
        dtrace[i] = damage[i].Trace()
        dev = damage[i].eigenValues()
        dmin[i] = dev.minElement()
        dmax[i] = dev.maxElement()
        svol[i] = strain[i].Trace()
        sev = strain[i].eigenValues()
        smin[i] = sev.minElement()
        smax[i] = sev.maxElement()
    dumpPhysicsState(integrator,
                     filename,
                     visitDir,
                     fields = [damageModel.sumActivationEnergiesPerNode(),
                               damageModel.numFlawsPerNode(),
                               VISARa.nodeFlags,
                               VISARb.nodeFlags,
                               VISARc.nodeFlags,
                               dtrace,
                               dmin,
                               dmax,
                               svol,
                               smin,
                               smax] + fields,
                     )


#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
    control.setRestartBaseName(restartBaseName)
    control.setFrequency(control.updateDomainDistribution, redistributeStep)
    VISARa.flushHistory()
    VISARb.flushHistory()
    VISARc.flushHistory()

else:
    control.iterateIdealH()
    control.smoothState(smoothIters)
    control.dropRestartFile()
    viz()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()
    viz()
