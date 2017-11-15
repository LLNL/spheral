#-------------------------------------------------------------------------------
# A tube of stainless steel undergoing expansion due to a projectile entering
# and smacking into a stop.
#
# See Vogler et al., 2003, International Journal of Impact Engineering, 29, 735
#-------------------------------------------------------------------------------
import sys
sys.path.append("../Utilities")
from SolidSpheral3d import *
from SpheralTestUtilities import *
from SpheralController import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from math import *

from bevelTubeEntrance import *
from NodeHistory import NodeHistory
from AverageStrain import AverageStrain
import mpi

#-------------------------------------------------------------------------------
# Create a physics package to dump excessive velocities into thermal energy.
#-------------------------------------------------------------------------------
class VelocityDiffuser(Physics):

    def __init__(self, vmax):
        self.vmax = vmax
        self.vmax2 = vmax*vmax
        Physics.__init__(self)
        return

    def evaluateDerivatives(self, time, dt,dataBase, state, derivs):
        return

    def dt(self, dataBase, state, derivs, time):
        return pair_double_string(1e30, "No vote")

    def registerState(self, dataBase, state):
        return

    def registerDerivatives(self, dataBase, derivs):
        return

    def finalize(self, time, dt, dataBase, state, derivs):
        velocities = state.vectorFields(HydroFieldNames.velocity)
        specificEnergies = state.scalarFields(HydroFieldNames.specificThermalEnergy)
        assert len(velocities) == len(specificEnergies)
        for vel, eps in zip(velocities, specificEnergies):
            assert vel.numInternalElements == eps.numInternalElements
            for i in xrange(vel.numInternalElements):
                vmag2 = vel[i].magnitude2()
                if vmag2 > self.vmax2:
                    deps = 0.5*(vmag2 - self.vmax2)
                    assert deps > 0.0
                    vel[i] *= self.vmax/sqrt(vmag2)
                    eps[i] += deps
        return
        
#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("3-D expanding tube impact strength/damage model test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(
    NodeListConstructor = AsphSolidNodeList,

    # How much of the 2 Pi geometry are we doing?
    phiFactor = 0.5,

    # Should we start the projectile at the end of the tube or in contact with
    # the plug?
    projectileOutside = True,

    # How much should we compress the projectile to allow it to slide into the tube?
    compressProjectile = 1.5,

    # Geometry
    tubeThickness = 0.3,
    rtubeInner = 0.5*1.27,
    ltube = 5.08,
    lCuAnvil = 1.0,
    lFoamAnvil = 1.0,
    lSteelAnvil = 1.0,

    # Numbers of nodes.
    nrtube = 15,
    nltube = 170, # 254,
    nrAnvil = 100,
    nrSteelAnvilCap = 10,

    # VISAR sampling parameters.
    dzVISAR = 0.04,
    drVISAR = 0.04,
    VISARsampleFrequency = 10,

    # Frequency for measuring the strain in the tube.
    strainFrequency = 5,

    # Inital velocity of the projectile.
    vzproj = -1.75e5, # -1.92e5,

    # Maximum spee we're going to allow.
    vmax = 3e5,

    # Steel EOS modifiers.
    C0multiplier = 1.0,

    # Mass densities.
    rho0Steel = 7.94, # 7.85
    rho0Lexan = 1.196,
    rho0Air = 0.001205,
    rho0Foam = 1.046,
    rho0Cu = 8.93,

    # Parameters for the damage model.
    DamageModelConstructor = GradyKippTensorDamageOwen,
    strainType = PseudoPlasticStrain,
    effectiveDamage = CopyDamage,
    effectiveFlawAlgorithm = SampledFlaws,
    useDamageGradient = True,
    kWeibullSteelFactor = 1.0,
    mWeibullSteelFactor = 1.0,
    randomSeedSteel = 109482993,
    criticalDamageThreshold = 0.5,

    # Node seeding stuff.
    nPerh = 1.51,

    # Material specific bounds on the mass density.
    etaMinSteel = 0.9,
    etaMaxSteel = 1.1,
    etaMinLexan = 0.5,
    etaMaxLexan = 1.5,
    etaMinCu = 0.7,
    etaMaxCu = 1.3,
    etaMinFoam = 0.5,
    etaMaxFoam = 1.5,

    # Material specific limits on the smoothing scale.
    hminratio = 0.1,
    hmin = 1e-2,
    hmax = 0.5,

    # Hydro parameters.
    Qconstructor = TensorMonaghanGingoldViscosity,
    Cl = 1.0,
    Cq = 1.0,
    Qlimiter = True,
    balsaraCorrection = False,
    epsilon2 = 1e-2,
    cfl = 0.5,
    XSPH = True,
    epsilonTensile = 0.3,
    nTensile = 4,
    HEvolution = IdealH,
    sumForMassDensity = IntegrateDensity,
    compatibleEnergyEvolution = True,
    gradhCorrection = False,
    useVelocityMagnitudeForDt = True,

    # Times, and simulation control.
    steps = None,
    goalTime = 100.0e-6,
    dtSample = 0.5e-6,
    dt = 5e-9,
    dtMin = 1e-11,
    dtMax = 1e-3,
    dtGrowth = 10.0,
    maxSteps = 200,
    statsStep = 10,
    redistributeStep = 200,
    smoothIters = 0,
    dtverbose = False,

    # Restart and output files.
    restoreCycle = None,
    restartStep = 200,
    baseDir = "dumps-expandingTube-3d",
    )

# Derived geometry.
assert phiFactor in (0.5, 1.0, 2.0)
phi = phiFactor * pi
rtubeOuter = rtubeInner + tubeThickness
rplug, lplug = rtubeInner, 0.5*ltube
rproj, lproj = rplug, lplug
rAnvil = 2.0*rtubeOuter
lAnvil = lCuAnvil + lFoamAnvil + lSteelAnvil

# Use the above geometry to define enclosing points of the materials for the
# node generators.
rminTube = rtubeInner
rmaxTube = rtubeOuter
zminTube = lAnvil
zmaxTube = lAnvil + ltube

rminPlug = 0.0
rmaxPlug = rplug
zminPlug = lAnvil
zmaxPlug = lAnvil + lplug

rminProj = 0.0
rmaxProj = rproj
if projectileOutside:
    zminProj = lAnvil + ltube
else:
    zminProj = zmaxPlug
zmaxProj = zminProj + lproj

rminInnerAir = 0.0
rmaxInnerAir = rtubeInner
zminInnerAir = lAnvil + lplug
zmaxInnerAir = lAnvil + ltube

rminSteelAnvil = 0.0
rmaxSteelAnvil = rAnvil
zminSteelAnvil = 0.0
zmaxSteelAnvil = lSteelAnvil

rminFoamAnvil = 0.0
rmaxFoamAnvil = rAnvil
zminFoamAnvil = lSteelAnvil
zmaxFoamAnvil = lSteelAnvil + lFoamAnvil

rminCuAnvil = 0.0
rmaxCuAnvil = rAnvil
zminCuAnvil = lSteelAnvil + lFoamAnvil
zmaxCuAnvil = lAnvil

rminSteelAnvilCap = rAnvil
rmaxSteelAnvilCap = rAnvil + float(nrSteelAnvilCap)/float(nrAnvil)*rAnvil
zminSteelAnvilCap = 0.0
zmaxSteelAnvilCap = lAnvil

# The geometry of the bevel on the inner tube opening surface.
tubeOpeningAngle = 1.8 * pi/180.0 # radians
zBevelBegin = lAnvil + ltube - 0.6

# Define the VISAR sampling points.
zVISARa = lAnvil + 2.5
zVISARb = lAnvil + 2.0
zVISARc = lAnvil + 1.5

rVISARa = rtubeOuter
rVISARb = rtubeOuter
rVISARc = rtubeOuter

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

# Inital velocity of the projectile.
v0proj = Vector(0.0, 0.0, vzproj)

# Damage model parameters.
kWeibullSteel = 8.8e5 * kWeibullSteelFactor
mWeibullSteel = 2.63  * mWeibullSteelFactor
volumeSteel = pi*(rtubeOuter**2 - rtubeInner**2)*ltube

# Restart and output files.
dataDir = os.path.join(baseDir, 
                       str(DamageModelConstructor).split("'")[1],
                       str(NodeListConstructor).split()[1],
                       "vzproj=%4.2e" % abs(vzproj),
                       "rho0Steel=%8.6f" % rho0Steel,
                       "C0mult=%4.2f" % C0multiplier,
                       "k=%4.2f_m=%4.2f" % (kWeibullSteel, mWeibullSteel))
restartDir = os.path.join(dataDir, "restarts", "proc-%04i" % mpi.rank)
visitDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "ExpandingTube")

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
# This attempts to mock the hardened steel properties.
#-------------------------------------------------------------------------------
## eosSteel = LinearPolynomialEquationOfStateCGS(rho0Steel,    # reference density  
##                                               etaMinSteel,  # etamin             
##                                               etaMaxSteel,  # etamax             
##                                               0.0,          # A0
##                                               1.649901e12,  # A1
##                                               1.674656e12,  # A2
##                                               0.832543e12,  # A3
##                                               1.93,         # B1
##                                               0.5,          # B2
##                                               0.0,          # B3
##                                               55.350)       # atomic weight

eosSteel = GruneisenEquationOfStateCGS(rho0Steel,    # reference density  
                                       etaMinSteel,  # etamin             
                                       etaMaxSteel,  # etamax             
                                       0.4529e6 * C0multiplier,     # C0                 
                                       1.50,         # S1                 
                                       0.0,          # S2                 
                                       0.0,          # S3                 
                                       1.84,         # gamma0             
                                       0.302,        # b                  
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

meltFitSteel = NinthOrderPolynomialFit( 7.40464217e10,
                                        2.49802214e11,
                                        1.00445029e12,
                                       -1.36451475e11,
                                        7.72897829e9,
                                        5.06390305e10,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0)

strengthSteel = SteinbergGuinanStrengthCGS(eosSteel,
                                           7.700000e11,        # G0
                                           2.2600e-12,         # A
                                           4.5500e-04,         # B
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
eosLexan = GruneisenEquationOfStateCGS(rho0Lexan,    # reference density  
                                         etaMinLexan,  # etamin             
                                         etaMaxLexan,  # etamax             
                                         0.1933e6,     # C0                 
                                         3.49,         # S1                 
                                        -8.2,          # S2                 
                                         9.6,          # S3                 
                                         0.61,         # gamma0             
                                         0.0,          # b                  
                                         28423.0)      # atomic weight

## strengthLexan = ConstantStrength(0.05e12,   # G0
##                                  0.0005e12) # Y0

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
meltFitLexan = NinthOrderPolynomialFit( 5.24383771e8,
                                        1.49188457e9,
                                        2.85704428e10,
                                        2.13783662e10,
                                       -4.45135120e9,
                                        1.24138074e9,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0)
strengthLexan = SteinbergGuinanStrengthCGS(eosLexan,
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
eosCu = GruneisenEquationOfStateCGS(rho0Cu,      # reference density  
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
meltFitCu = NinthOrderPolynomialFit( 5.22055639e10,
                                     1.90143176e11,
                                     8.51351901e11,
                                    -1.12049022e11,
                                    -6.11436674e9,
                                     4.36007831e10,
                                     0.0,
                                     0.0,
                                     0.0,
                                     0.0)
strengthCu = SteinbergGuinanStrengthCGS(eosCu,
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
eosFoam = GruneisenEquationOfStateCGS(rho0Foam,  # reference density  
                                      etaMinFoam,# etamin             
                                      etaMaxFoam,# etamax             
                                       0.189e6,  # C0                 
                                       2.965,    # S1                 
                                      -4.069,    # S2                 
                                       2.328,    # S3                 
                                       0.67,     # gamma0             
                                       0.0,      # b                  
                                       6.982)    # atomic weight
strengthFoam = NullStrength()

#-------------------------------------------------------------------------------
# Air material properties.
#-------------------------------------------------------------------------------
gammaAir = 1.4
molecularWeightAir = 30.0
eosAir = GammaLawGasCGS(gammaAir, molecularWeightAir)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodesSteel = NodeListConstructor("Stainless steel", eosSteel, strengthSteel, WT, WTPi)
nodesPlug = NodeListConstructor("Lexan plug", eosLexan, strengthLexan, WT, WTPi)
nodesProj = NodeListConstructor("Lexan projectile", eosLexan, strengthLexan, WT, WTPi)
nodesSteelAnvil = NodeListConstructor("Anvil (Steel)", eosSteel, strengthSteel, WT, WTPi)
nodesFoamAnvil = NodeListConstructor("Anvil (Foam)", eosFoam, strengthFoam, WT, WTPi)
nodesCuAnvil = NodeListConstructor("Anvil (Copper)", eosCu, strengthCu, WT, WTPi)
nodeSet = [nodesSteel, nodesPlug, nodesProj,
           nodesSteelAnvil, nodesFoamAnvil, nodesCuAnvil]
for (n, etamin, etamax, rho0) in zip(nodeSet,
                                     [etaMinSteel, etaMinLexan, etaMinLexan, etaMinSteel, etaMinFoam, etaMinCu],
                                     [etaMaxSteel, etaMaxLexan, etaMaxLexan, etaMaxSteel, etaMaxFoam, etaMaxCu],
                                     [rho0Steel, rho0Lexan, rho0Lexan, rho0Steel, rho0Foam, rho0Cu]):
    n.nodesPerSmoothingScale = nPerh
    n.epsilonTensile = epsilonTensile
    n.nTensile = nTensile
    n.hmin = hmin
    n.hmax = hmax
    n.hminratio = hminratio
    n.XSPH = XSPH
    n.rhoMin = etamin*rho0
    n.rhoMax = etamax*rho0
    output("n.name")
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
    neighbor = TreeNeighbor(n,
                                  kernelExtent = WT.kernelExtent)
    n.registerNeighbor(neighbor)
    cache.append(neighbor)
del n

#-------------------------------------------------------------------------------
# Set node properties (positions, velocites, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    print "Generating node distribution."
    from GenerateNodeDistribution3d import *
    from CompositeNodeDistribution import *
    from VoronoiDistributeNodes import distributeNodes3d
    generatorTube = GenerateNodeDistribution3d(nrtube,
                                               nltube,
                                               0,
                                               rho0Steel,
                                               "cylindrical",
                                               rmin = rminTube,
                                               rmax = rmaxTube,
                                               thetamin = 0,
                                               thetamax = phi,
                                               zmin = zminTube,
                                               zmax = zmaxTube,
                                               nNodePerh = nPerh,
                                               SPH = (NodeListConstructor is SphNodeList))
    generatorPlug = GenerateNodeDistribution3d(nrplug,
                                               nlplug,
                                               0,
                                               rho0Lexan,
                                               "cylindrical",
                                               rmin = rminPlug,
                                               rmax = rmaxPlug,
                                               thetamin = 0,
                                               thetamax = phi,
                                               zmin = zminPlug,
                                               zmax = zmaxPlug,
                                               nNodePerh = nPerh,
                                               SPH = (NodeListConstructor is SphNodeList))
    generatorProj = GenerateNodeDistribution3d(nrproj,
                                               nlproj,
                                               0,
                                               rho0Lexan,
                                               "cylindrical",
                                               rmin = rminProj,
                                               rmax = rmaxProj,
                                               thetamin = 0,
                                               thetamax = phi,
                                               zmin = zminProj,
                                               zmax = zmaxProj,
                                               nNodePerh = nPerh,
                                               SPH = (NodeListConstructor is SphNodeList))
    generatorSteelAnvil1 = GenerateNodeDistribution3d(nrSteelAnvil,
                                                      nlSteelAnvil,
                                                      0,
                                                      rho0Steel,
                                                      "cylindrical",
                                                      rmin = rminSteelAnvil,
                                                      rmax = rmaxSteelAnvil,
                                                      thetamin = 0,
                                                      thetamax = phi,
                                                      zmin = zminSteelAnvil,
                                                      zmax = zmaxSteelAnvil,
                                                      nNodePerh = nPerh,
                                                      SPH = (NodeListConstructor is SphNodeList))
    generatorSteelAnvil2 = GenerateNodeDistribution3d(nrSteelAnvilCap,
                                                      nlSteelAnvilCap,
                                                      0,
                                                      rho0Steel,
                                                      "cylindrical",
                                                      rmin = rminSteelAnvilCap,
                                                      rmax = rmaxSteelAnvilCap,
                                                      thetamin = 0,
                                                      thetamax = phi,
                                                      zmin = zminSteelAnvilCap,
                                                      zmax = zmaxSteelAnvilCap,
                                                      nNodePerh = nPerh,
                                                      SPH = (NodeListConstructor is SphNodeList))
    generatorSteelAnvil = CompositeNodeDistribution(generatorSteelAnvil1, generatorSteelAnvil2)
    generatorFoamAnvil = GenerateNodeDistribution3d(nrFoamAnvil,
                                                     nlFoamAnvil,
                                                     0,
                                                     rho0Foam,
                                                     "cylindrical",
                                                     rmin = rminFoamAnvil,
                                                     rmax = rmaxFoamAnvil,
                                                     thetamin = 0,
                                                     thetamax = phi,
                                                     zmin = zminFoamAnvil,
                                                     zmax = zmaxFoamAnvil,
                                                     nNodePerh = nPerh,
                                                    SPH = (NodeListConstructor is SphNodeList))
    generatorCuAnvil = GenerateNodeDistribution3d(nrCuAnvil,
                                                     nlCuAnvil,
                                                     0,
                                                     rho0Cu,
                                                     "cylindrical",
                                                     rmin = rminCuAnvil,
                                                     rmax = rmaxCuAnvil,
                                                     thetamin = 0,
                                                     thetamax = phi,
                                                     zmin = zminCuAnvil,
                                                     zmax = zmaxCuAnvil,
                                                     nNodePerh = nPerh,
                                                  SPH = (NodeListConstructor is SphNodeList))
              
    distributeNodes3d((nodesSteel, generatorTube),
                      (nodesPlug, generatorPlug),
                      (nodesProj, generatorProj),
                      (nodesSteelAnvil, generatorSteelAnvil),
                      (nodesFoamAnvil, generatorFoamAnvil),
                      (nodesCuAnvil, generatorCuAnvil))
    nGlobalNodes = 0
    for n in nodeSet:
        print "Generator info for %s" % n.name
        output("    mpi.allreduce(n.numInternalNodes, mpi.MIN)")
        output("    mpi.allreduce(n.numInternalNodes, mpi.MAX)")
        output("    mpi.allreduce(n.numInternalNodes, mpi.SUM)")
        nGlobalNodes += mpi.allreduce(n.numInternalNodes, mpi.SUM)
    del n
    print "Total number of (internal) nodes in simulation: ", nGlobalNodes

    # Bevel the inner opening surface of the target tube.
    numNodesBeveled = bevelTubeEntrance(nodesSteel,
                                        3,
                                        tubeOpeningAngle,
                                        rtubeInner,
                                        tubeThickness,
                                        zBevelBegin)
    print "Beveled %i nodes in the tube opening." % mpi.allreduce(numNodesBeveled,
                                                                  mpi.SUM)

    # Adjust the diameter of the projectile inward a bit, so it will slide
    # into the tube properly.
    drProj = compressProjectile*nPerh*rproj/nrproj
    projMultiplier = (rproj - drProj)/rproj
    for i in xrange(nodesProj.numInternalNodes):
        nodesProj.positions()[i].x *= projMultiplier
        nodesProj.positions()[i].y *= projMultiplier

    # Adjust the plug to match.
    for i in xrange(nodesPlug.numInternalNodes):
        nodesPlug.positions()[i].x *= projMultiplier
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
        n.specificThermalEnergy(ScalarField("tmp", n, u0))
        print "Initial pressure for %s: %g" % (n.name,
                                               n.equationOfState().pressure(rho0, u0))
    del n

    # Set the projectile velocities.
    nodesProj.velocity(VectorField("tmp", nodesProj, v0proj))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase()
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
hydro = Hydro(W = WT,
              Q = q,
              compatibleEnergyEvolution = compatibleEnergyEvolution,
              gradhCorrection = gradhCorrection,
              densityUpdate = sumForMassDensity,
              HUpdate = HEvolution)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = useVelocityMagnitudeForDt
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.hmin")
output("hydro.hmax")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct a strength physics object.
#-------------------------------------------------------------------------------
strength = Strength()
output("strength")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
volumeMultiplier = 1.0
numFlawsPerNode = 1
damageModel = DamageModelConstructor(nodesSteel,
                                     kWeibull = kWeibullSteel,
                                     mWeibull = mWeibullSteel,
                                     kernel = WT,
                                     seed = randomSeedSteel,
                                     volumeMultiplier = volumeMultiplier,
                                     strainAlgorithm = strainType,
                                     effectiveDamageAlgorithm = effectiveDamage,
                                     useDamageGradient = useDamageGradient,
                                     flawAlgorithm = effectiveFlawAlgorithm,
                                     criticalDamageThreshold = criticalDamageThreshold,
                                     numFlawsPerNode = numFlawsPerNode)

output("damageModel")
output("damageModel.useDamageGradient")
output("damageModel.effectiveDamageAlgorithm")
output("damageModel.effectiveFlawAlgorithm")

#-------------------------------------------------------------------------------
# Construct the velocity diffuser to keep the velocities under control.
#-------------------------------------------------------------------------------
CHP = VelocityDiffuser(vmax)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.appendPhysicsPackage(damageModel)
integrator.appendPhysicsPackage(CHP)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(strength)")
output("integrator.havePhysicsPackage(damageModel)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Construct boundary conditions, and add them to our physics packages.
#-------------------------------------------------------------------------------
xbcPlane = Plane(Vector(0.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0))
ybcPlane = Plane(Vector(0.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0))
zbcPlane = Plane(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 1.0))
xbc = ReflectingBoundary(xbcPlane)
ybc = ReflectingBoundary(ybcPlane)
zbc = ReflectingBoundary(zbcPlane)

if phiFactor == 0.5:
    bcs = [xbc, ybc, zbc]
elif phiFactor == 1.0:
    bcs = [ybc, zbc]
else:
    bcs = [zbc]
for package in integrator.physicsPackages():
    for bc in bcs:
        package.appendBoundary(bc)

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
# Monitor the evolution of the mass averaged strain.
#-------------------------------------------------------------------------------
strainHistory = AverageStrain(damageModel,
                              dataDir + "/strainhistory.txt")
control.appendPeriodicWork(strainHistory.sample, strainFrequency)

#-------------------------------------------------------------------------------
# Select the nodes for the VISAR sampling.
#-------------------------------------------------------------------------------
def xymagnitude(vec):
    return sqrt((vec.x)**2 + (vec.y)**2)

class SelectVISARNodes:
    def __init__(self, z0, dz, dr, label):
        self.z0 = z0
        self.dz = dz
        self.dr = dr
        self.label = label
        return
    def __call__(self, nodes):
        r = nodes.positions()
        potentials = [i for i in xrange(nodes.numInternalNodes)
                      if abs(r[i].z - self.z0) < self.dz]
        rxymax = mpi.allreduce(max([xymagnitude(r[i]) for i in potentials] + [-1e50]), mpi.MAX)
        result = [i for i in potentials if xymagnitude(r[i]) > rxymax - self.dr]
        print "Selected %i %s velocimetry test points." % (mpi.allreduce(len(result), mpi.SUM),
                                                           self.label)
        return result

#-------------------------------------------------------------------------------
# Sampling function to measure the average velocity at the VISAR probe sites.
#-------------------------------------------------------------------------------
def averageCylindricalRadialVelocity(nodes, indicies):
    r = nodes.positions()
    m = nodes.mass()
    v = nodes.velocity()

    massSum = 1e-30
    result = 0.0
    for i in indicies:
        assert i >= 0 and i < nodes.numInternalNodes
        massSum += m[i]
        runit = Vector(r[i].x, r[i].y, 0.0).unitVector()
        result += m[i] * (v[i].dot(runit))

    globalMassSum = mpi.allreduce(massSum, mpi.SUM)
    globalResult = mpi.allreduce(result, mpi.SUM)
    assert globalMassSum > 0.0
    globalResult /= globalMassSum
    return globalResult

#-------------------------------------------------------------------------------
# Build the history objects to simulate the VISAR velocity probes.
#-------------------------------------------------------------------------------
nodesA = SelectVISARNodes(zVISARa, dzVISAR, drVISAR, "A")
nodesB = SelectVISARNodes(zVISARb, dzVISAR, drVISAR, "B")
nodesC = SelectVISARNodes(zVISARc, dzVISAR, drVISAR, "C")
VISARa = NodeHistory(nodesSteel, nodesA, averageCylindricalRadialVelocity,
                     dataDir + "/VISAR-a")
VISARb = NodeHistory(nodesSteel, nodesB, averageCylindricalRadialVelocity,
                     dataDir + "/VISAR-b")
VISARc = NodeHistory(nodesSteel, nodesC, averageCylindricalRadialVelocity,
                     dataDir + "/VISAR-c")
VISARa.nodeFlags.name = "VISAR a points"
VISARb.nodeFlags.name = "VISAR b points"
VISARc.nodeFlags.name = "VISAR c points"

control.appendPeriodicWork(VISARa.sample, VISARsampleFrequency)
control.appendPeriodicWork(VISARb.sample, VISARsampleFrequency)
control.appendPeriodicWork(VISARc.sample, VISARsampleFrequency)

#-------------------------------------------------------------------------------
# Helper method for dumping viz files.
#-------------------------------------------------------------------------------
def viz(fields = [],
        filename = "ExpandingTube-3d"):
    tdamage = nodesSteel.damage()
    etdamage = nodesSteel.effectiveDamage()
    tstrain = damageModel.strain()
    etstrain = damageModel.effectiveStrain()
    sdamage = ScalarField("damage magnitude", nodesSteel)
    mindamage = ScalarField("damage magnitude min", nodesSteel)
    maxdamage = ScalarField("damage magnitude max", nodesSteel)
    esdamage = ScalarField("effective damage magnitude", nodesSteel)
    minedamage = ScalarField("effective damage magnitude min", nodesSteel)
    maxedamage = ScalarField("effective damage magnitude max", nodesSteel)
    sstrain = ScalarField("strain average", nodesSteel)
    esstrain = ScalarField("effective strain average", nodesSteel)
    for i in xrange(nodesSteel.numInternalNodes):
        sdamage[i] = tdamage[i].Trace()
        esdamage[i] = etdamage[i].Trace()
        ev = tdamage[i].eigenValues()
        eev = etdamage[i].eigenValues()
        maxdamage[i] = ev.maxElement()
        mindamage[i] = ev.minElement()
        maxedamage[i] = eev.maxElement()
        minedamage[i] = eev.minElement()
        sstrain[i] = tstrain[i].Trace()/3.0
        esstrain[i] = etstrain[i].Trace()/3.0
    dumpPhysicsState(integrator,
                     filename,
                     visitDir,
                     fields = [damageModel.sumActivationEnergiesPerNode(),
                               damageModel.numFlawsPerNode(),
                               sdamage, mindamage, maxdamage,
                               esdamage, minedamage, maxedamage,
                               sstrain, esstrain] + [x.nodeFlags for x in (VISARa, VISARb, VISARc)] + fields,
                         )

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
    control.setRestartBaseName(restartBaseName)
    control.setFrequency(control.updateDomainDistribution, redistributeStep)
    strainHistory.flushHistory()
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
if not steps is None:
    control.step(steps)
    raise ValueError, "Completed %i steps." % steps
else:
    while control.time() < goalTime:
        nextGoalTime = min(int((control.time() + 1.001*dtSample)/dtSample)*dtSample,
                           goalTime)
        control.advance(nextGoalTime, maxSteps)
        control.dropRestartFile()
        viz()

#-------------------------------------------------------------------------------
# Now we can collect some info on the fragment population.
#-------------------------------------------------------------------------------
from identifyFragments import *
linkRadius = 1.0
damageThreshold = 0.9
fragIDs = identifyFragments(nodesSteel,
                            linkRadius,
                            damageModel.damage(),
                            damageThreshold,
                            True)
fragProps = fragmentProperties(nodesSteel, fragIDs, damageModel.strain())
viz(fields = [fragIDs],
    filename = "ExpandingTube-3d-frags-nodust-linkRadius=%f-damageThreshold=%f" %
    (linkRadius, damageThreshold))
if mpi.rank == 0:
    f = open("frag-props-t=%8.3fusec-lr=%f-dt=%f" % (control.time()/1e-6, linkRadius, damageThreshold), "w")
    import pickle
    pickle.dump(fragProps, f)
    f.close()
