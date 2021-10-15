from Numeric import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralController import *
from findLastRestart import *
from SpheralVisitDump import SpheralVisitDump
from math import *

from Strength import *

# Load the mpi module if we're parallel.
mpi, procID, numProcs = loadmpi()

#-------------------------------------------------------------------------------
# nint method, since python doesn't provide it.
#-------------------------------------------------------------------------------
def nint(x):
    return int(x + 0.5*x/(sqrt(x*x) + 1.0e-30))

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
seed = 'lattice'

Sapphire1Thickness = 1.15
TantalumThickness = 0.50225
Sapphire2Thickness = 0.3171
TungstenCarbideThickness = 0.16
PMMAThickness = 0.616
yThickness = 0.2

Sapphire1Range = ((0.0, 0.0),
                  (Sapphire1Thickness, yThickness))
TantalumRange = ((Sapphire1Range[1][0], 0.0),
                 (Sapphire1Range[1][0] + TantalumThickness, yThickness))
Sapphire2Range = ((TantalumRange[1][0], 0.0),
                  (TantalumRange[1][0] + Sapphire2Thickness, yThickness))
TungstenCarbideRange = ((Sapphire2Range[1][0], 0.0),
                        (Sapphire2Range[1][0] + TungstenCarbideThickness, yThickness))
PMMARange = ((TungstenCarbideRange[1][0], 0.0),
             (TungstenCarbideRange[1][0] + PMMAThickness, yThickness))

xmin = 0.0
xmax = PMMARange[1][0]
ymin = 0.0
ymax = yThickness

nxSapphire1 = 181
nxTantalum = 251
nxSapphire2 = 50
nxTungstenCarbide = 60
nxPMMA = 30

nySapphire1 = nint(nxSapphire1 * yThickness/Sapphire1Thickness)
nyTantalum = nint(nxTantalum * yThickness/TantalumThickness)
nySapphire2 = nint(nxSapphire2 * yThickness/Sapphire2Thickness)
nyTungstenCarbide = nint(nxTungstenCarbide * yThickness/TungstenCarbideThickness)
nyPMMA = nint(nxPMMA * yThickness/PMMAThickness)

nPerh = 2.01

rhoSapphire1 = 3.985
rhoTantalum = 16.69
rhoSapphire2 = 3.985
rhoTungstenCarbide = 14.9
rhoPMMA = 1.182

mSapphire1 = (Sapphire1Thickness*yThickness)*rhoSapphire1/(nxSapphire1*nySapphire1)
mTantalum = (TantalumThickness*yThickness)*rhoTantalum/(nxTantalum*nyTantalum)
mSapphire2 = (Sapphire2Thickness*yThickness)*rhoSapphire2/(nxSapphire2*nySapphire2)
mTungstenCarbide = (TungstenCarbideThickness*yThickness)*rhoTungstenCarbide/(nxTungstenCarbide*nyTungstenCarbide)
mPMMA = (PMMAThickness*yThickness)*rhoPMMA/(nxPMMA*nyPMMA)

HSapphire1 = SymTensor2d(1.0/(nPerh*Sapphire1Thickness/nxSapphire1), 0.0,
                         0.0, 1.0/(nPerh*Sapphire1Thickness/nxSapphire1))
HTantalum = SymTensor2d(1.0/(nPerh*TantalumThickness/nxTantalum), 0.0,
                        0.0, 1.0/(nPerh*TantalumThickness/nxTantalum))
HSapphire2 = SymTensor2d(1.0/(nPerh*Sapphire2Thickness/nxSapphire2), 0.0,
                         0.0, 1.0/(nPerh*Sapphire2Thickness/nxSapphire2))
HTungstenCarbide = SymTensor2d(1.0/(nPerh*TungstenCarbideThickness/nxTungstenCarbide), 0.0,
                               0.0, 1.0/(nPerh*TungstenCarbideThickness/nxTungstenCarbide))
HPMMA = SymTensor2d(1.0/(nPerh*PMMAThickness/nxPMMA), 0.0,
                    0.0, 1.0/(nPerh*PMMAThickness/nxPMMA))

v1 = Vector2d(-1.33e4, 0.0)

Qconstructor = MonaghanGingoldViscosity2d
Cl, Cq = 0.5, 1.0
Qlimiter = False
balsaraCorrection = False
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
HsmoothMin, HsmoothMax = 0.0004, 0.5
cfl = 0.5
useVelocityMagnitudeForDt = False
XSPH = False
epsilonTensile = 0.1
nTensile = 4
hybridMassDensityThreshold = 0.01

goalTime = 5.0e-6
dtSample = 1e-7
dt = 1e-10
dtMin, dtMax = 1e-12, 1e-5
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Hydro2d.HEvolutionType.IdealH
sumForMassDensity = Hydro2d.MassDensityType.IntegrateDensity # HybridDensity # CorrectedSumDensity

restartStep = 1000
restartBaseName = "dumps-2d/PlateImpact-%i-%i-%i-%i-%i" % (nxSapphire1,
                                                           nxTantalum,
                                                           nxSapphire2,
                                                           nxTungstenCarbide,
                                                           nxPMMA)
restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Define a class to track the history of the interface.
#-------------------------------------------------------------------------------
class InterfaceHistory:

    def __init__(self,
                 sapphireIndicies,
                 tantalumIndicies,
                 sapphireNodes,
                 tantalumNodes,
                 filename):
        self.restart = RestartableObject(self)
        self.sapphireIndicies = sapphireIndicies
        self.tantalumIndicies = tantalumIndicies
        self.sapphireNodes = sapphireNodes
        self.tantalumNodes = tantalumNodes
        self.filename = filename
        self.file = None
        self.timeHistory = []
        self.sapphireVelocity = []
        self.tantalumVelocity = []
        return

    def sample(self, t):
        sapphireVelocity0 = 0.0
        tantalumVelocity0 = 0.0
        for i in self.sapphireIndicies:
            sapphireVelocity0 += self.sapphireNodes.velocity()[i].x
        for i in self.tantalumIndicies:
            tantalumVelocity0 += self.tantalumNodes.velocity()[i].x
        sapphireVelocity = mpi.allreduce(sapphireVelocity0, mpi.SUM)
        tantalumVelocity = mpi.allreduce(tantalumVelocity0, mpi.SUM)
        nSapphire = mpi.allreduce(len(self.sapphireIndicies), mpi.SUM)
        nTantalum = mpi.allreduce(len(self.tantalumIndicies), mpi.SUM)
        assert nSapphire > 0
        assert nTantalum > 0
        sapphireVelocity /= nSapphire
        tantalumVelocity /= nTantalum
        self.timeHistory.append(t)
        self.sapphireVelocity.append(sapphireVelocity)
        self.tantalumVelocity.append(tantalumVelocity)
        if mpi.rank == 0:
            if self.file is None:
                self.file = open(self.filename, "w")
            self.file.write("%g \t%g \t%g\n" % (t, sapphireVelocity, tantalumVelocity))
            self.file.flush()
        return

    def label(self):
        return "InterfaceHistory"

    def dumpState(self, file, path):
        file.writeObject(self.sapphireIndicies, path + "/sapphireIndicies")
        file.writeObject(self.tantalumIndicies, path + "/tantalumIndicies")
        file.writeObject(self.filename, path + "/filename")
        file.writeObject(self.timeHistory, path + "/timeHistory")
        file.writeObject(self.sapphireVelocity, path + "/sapphireVelocity")
        file.writeObject(self.tantalumVelocity, path + "/tantalumVelocity")
        return

    def restoreState(self, file, path):
        self.sapphireIndicies = file.readObject(path + "/sapphireIndicies")
        self.tantalumIndicies = file.readObject(path + "/tantalumIndicies")
        self.filename = file.readObject(path + "/filename")
        self.timeHistory = file.readObject(path + "/timeHistory")
        self.sapphireVelocity = file.readObject(path + "/sapphireVelocity")
        self.tantalumVelocity = file.readObject(path + "/tantalumVelocity")
        return

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2-D Plate impact strength test")
##print "Sapphire 1 in range ((%f, %f), (%f, %f))" % Sapphire1Range
##print "Tantalum   in range ((%f, %f), (%f, %f))" % TantalumRange
##print "Sapphire 2 in range ((%f, %f), (%f, %f))" % Sapphire2Range
##print "TungCarbid in range ((%f, %f), (%f, %f))" % TungstenCarbideRange
##print "PMMA       in range ((%f, %f), (%f, %f))" % PMMARange

#-------------------------------------------------------------------------------
# Sapphire material properties.
#-------------------------------------------------------------------------------
eosSapphire  = GruneisenEquationOfStateCGS2d(3.985,  # reference density  
                                             0.0,    # etamin             
                                             1e30,   # etamax             
                                             1.119e6,# C0                 
                                             1.0,    # S1                 
                                             0.0,    # S2                 
                                             0.0,    # S3                 
                                             0.0,    # gamma0             
                                             0.0,    # b                  
                                             20.39)  # atomic weight      

#-------------------------------------------------------------------------------
# Tantalum material properties.
#-------------------------------------------------------------------------------
eosTantalum = GruneisenEquationOfStateCGS2d(16.69,   # reference density  
                            -0.2,     # etamin             
                             4.0,     # etamax             
                             0.341e6, # C0                 
                             1.2,     # S1                 
                             0.0,     # S2                 
                             0.0,     # S3                 
                             1.67,    # gamma0             
                             0.42,    # b                  
                             180.948) # atomic weight
coldFitTantalum = NinthOrderPolynomialFit(-6.86446714e9,
                                          -1.17070812e10,
                                           9.70252276e11,
                                          -4.12557402e11,
                                           1.13401823e11,
                                          -1.86584799e10,
                                           0.0,
                                           0.0,
                                           0.0,
                                           0.0)
meltFitTantalum = NinthOrderPolynomialFit(9.24414908e10,
                                          2.53949977e11,
                                          1.06113848e12,
                                         -5.28947636e11,
                                          1.67906438e11,
                                         -2.92459765e10,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0)
strengthTantalum = SteinbergGuinanStrengthCGS2d(eosTantalum,
                                                6.900000e11,        # G0
                                                1.4500e-12,         # A
                                                1.3000e-04,         # B
                                                7.7000e9,           # Y0
                                                1.1e10,             # Ymax
                                                1.0e-3,             # Yp
                                                10.0000,            # beta
                                                0.0,                # gamma0
                                                0.1,                # nhard
                                                coldFitTantalum,
                                                meltFitTantalum)

#-------------------------------------------------------------------------------
# Tungsten carbide material properties.
#-------------------------------------------------------------------------------
eosTungstenCarbide = GruneisenEquationOfStateCGS2d(14.9,   # reference density  
                                    0.0,    # etamin             
                                    1e30,   # etamax             
                                    0.519e6,# C0                 
                                    1.16,   # S1                 
                                    0.0,    # S2                 
                                    0.0,    # S3                 
                                    1.5,    # gamma0             
                                    0.0,    # b                  
                                    155.89) # atomic weight      
coldFitTungstenCarbide = NinthOrderPolynomialFit(0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0)
meltFitTungstenCarbide = NinthOrderPolynomialFit(0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0)
strengthTungstenCarbide = SteinbergGuinanStrengthCGS2d(eosTungstenCarbide,
                                                       2.540000e12,        # G0
                                                       0.0,                # A
                                                       0.0,                # B
                                                       5.1000e10,          # Y0
                                                       1.0e12,             # Ymax
                                                       1.0e-3,             # Yp
                                                       0.0,                # beta
                                                       0.0,                # gamma0
                                                       0.0,                # nhard
                                                       coldFitTungstenCarbide,
                                                       meltFitTungstenCarbide)

#-------------------------------------------------------------------------------
# PMMA material properties.
#-------------------------------------------------------------------------------
eosPMMA = GruneisenEquationOfStateCGS2d(1.182,  # reference density  
                         0.0,    # etamin             
                         1e30,   # etamax             
                         0.218e6,# C0                 
                         2.088,  # S1                 
                         -1.124, # S2                 
                         0.0,    # S3                 
                         0.85,   # gamma0             
                         0.0,    # b                  
                         17.035) # atomic weight      
coldFitPMMA = NinthOrderPolynomialFit(-5.19191852e9,
                                      -4.41500192e9,
                                       2.84720528e10,
                                       2.14093899e10,
                                      -4.46412259e9,
                                       1.24495222e9,
                                       0.0,
                                       0.0,
                                       0.0,
                                       0.0)
meltFitPMMA = NinthOrderPolynomialFit(5.24383771e8,
                                      1.49188457e9,
                                      2.85704428e10,
                                      2.13783662e10,
                                     -4.45135120e9,
                                      1.24138074e9,
                                      0.0,
                                      0.0,
                                      0.0,
                                      0.0)
strengthPMMA = SteinbergGuinanStrengthCGS2d(eosPMMA,
                                            2.320000e10,        # G0
                                            0.0,                # A
                                            0.0,                # B
                                            4.2000e9,           # Y0
                                            1.0e12,             # Ymax
                                            1.0e-3,             # Yp
                                            0.0,                # beta
                                            0.0,                # gamma0
                                            0.0,                # nhard
                                            coldFitPMMA,
                                            meltFitPMMA)

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodesSapphire1 = SphNodeList2d("Sapphire 1", eosSapphire)
nodesTantalum = SphSolidNodeList2d("Tantalum", eosTantalum, strengthTantalum)
nodesSapphire2 = SphNodeList2d("Sapphire 2", eosSapphire)
nodesTungstenCarbide = SphSolidNodeList2d("Tungsten Carbide", eosTungstenCarbide, strengthTungstenCarbide)
nodesPMMA = SphSolidNodeList2d("PMMA", eosPMMA, strengthPMMA)
nodeLists = [nodesSapphire1,
             nodesTantalum,
             nodesSapphire2,
             nodesTungstenCarbide,
             nodesPMMA]
for nodes in nodeLists:
    nodes.nodesPerSmoothingScale = nPerh
    nodes.epsilonTensile = epsilonTensile
    nodes.nTensile = nTensile
    nodes.XSPH = XSPH
    output('nodes.name()')
    output('nodes.nodesPerSmoothingScale')
    output('nodes.epsilonTensile')
    output('nodes.nTensile')
    output('nodes.XSPH')

#-------------------------------------------------------------------------------
# Create an instance of our history object to find the interface history.
#-------------------------------------------------------------------------------
interfaceHistory = InterfaceHistory(None, None,
                                    nodesSapphire1, nodesTantalum,
                                    "PlateImpact-2d-interface-history.txt")

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
output('WT')
output('WTPi')
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
_cache = []
for nodes in nodeLists:
    neighbor = TreeNeighbor2d(nodes,
                              kernelExtent = kernelExtent)
    nodes.registerNeighbor(neighbor)
    _cache.append(neighbor)  # This is a trick to make sure the neighbors are not garbage collected.
neighborTimer.stop()
neighborTimer.printStatus()

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from GenerateNodeDistribution2d import *
    from DistributeNodes import distributeNodes2d
    print "Generating node distribution."
    generatorSapphire1 = GenerateNodeDistribution2d(nxSapphire1,
                                                    nySapphire1,
                                                    rhoSapphire1,
                                                    'lattice',
                                                    xmin = Sapphire1Range[0],
                                                    xmax = Sapphire1Range[1],
                                                    nNodePerh = nPerh)
    generatorTantalum = GenerateNodeDistribution2d(nxTantalum,
                                                   nyTantalum,
                                                   rhoTantalum,
                                                   'lattice',
                                                   xmin = TantalumRange[0],
                                                   xmax = TantalumRange[1],
                                                   nNodePerh = nPerh)
    generatorSapphire2 = GenerateNodeDistribution2d(nxSapphire2,
                                                    nySapphire2,
                                                    rhoSapphire2,
                                                    'lattice',
                                                    xmin = Sapphire2Range[0],
                                                    xmax = Sapphire2Range[1],
                                                    nNodePerh = nPerh)
    generatorTungstenCarbide = GenerateNodeDistribution2d(nxTungstenCarbide,
                                                          nyTungstenCarbide,
                                                          rhoTungstenCarbide,
                                                          'lattice',
                                                          xmin = TungstenCarbideRange[0],
                                                          xmax = TungstenCarbideRange[1],
                                                          nNodePerh = nPerh)
    generatorPMMA = GenerateNodeDistribution2d(nxPMMA,
                                               nyPMMA,
                                               rhoPMMA,
                                               'lattice',
                                               xmin = PMMARange[0],
                                               xmax = PMMARange[1],
                                               nNodePerh = nPerh)
    nSapphire1 = generatorSapphire1.globalNumNodes()
    nTantalum = generatorTantalum.globalNumNodes()
    nSapphire2 = generatorSapphire2.globalNumNodes()
    nTungstenCarbide = generatorTungstenCarbide.globalNumNodes()
    nPMMA = generatorPMMA.globalNumNodes()
    nodeInfo = distributeNodes2d([(nodesSapphire1, nSapphire1, generatorSapphire1),
                                  (nodesTantalum, nTantalum, generatorTantalum),
                                  (nodesSapphire2, nSapphire2, generatorSapphire2),
                                  (nodesTungstenCarbide, nTungstenCarbide, generatorTungstenCarbide),
                                  (nodesPMMA, nPMMA, generatorPMMA)])

    # Set the node masses.
    nodesSapphire1.setMass(ScalarField2d("tmp", nodesSapphire1, mSapphire1))
    nodesTantalum.setMass(ScalarField2d("tmp", nodesTantalum, mTantalum))
    nodesSapphire2.setMass(ScalarField2d("tmp", nodesSapphire2, mSapphire2))
    nodesTungstenCarbide.setMass(ScalarField2d("tmp", nodesTungstenCarbide, mTungstenCarbide))
    nodesPMMA.setMass(ScalarField2d("tmp", nodesPMMA, mPMMA))

    # Set the smoothing scales.
    nodesSapphire1.setHfield(SymTensorField2d("tmp", nodesSapphire1, HSapphire1))
    nodesTantalum.setHfield(SymTensorField2d("tmp", nodesTantalum, HTantalum))
    nodesSapphire2.setHfield(SymTensorField2d("tmp", nodesSapphire2, HSapphire2))
    nodesTungstenCarbide.setHfield(SymTensorField2d("tmp", nodesTungstenCarbide, HTungstenCarbide))
    nodesPMMA.setHfield(SymTensorField2d("tmp", nodesPMMA, HPMMA))

    # Set the node mass densities.
    nodesSapphire1.setMassDensity(ScalarField2d("tmp", nodesSapphire1, rhoSapphire1))
    nodesTantalum.setMassDensity(ScalarField2d("tmp", nodesTantalum, rhoTantalum))
    nodesSapphire2.setMassDensity(ScalarField2d("tmp", nodesSapphire2, rhoSapphire2))
    nodesTungstenCarbide.setMassDensity(ScalarField2d("tmp", nodesTungstenCarbide, rhoTungstenCarbide))
    nodesPMMA.setMassDensity(ScalarField2d("tmp", nodesPMMA, rhoPMMA))

    # Set node specific thermal energies
    nodesSapphire1.setSpecificThermalEnergy(ScalarField2d("tmp", nodesSapphire1,
                                                          eosSapphire.specificThermalEnergy(rhoSapphire1, 300.0)))
    nodesTantalum.setSpecificThermalEnergy(ScalarField2d("tmp", nodesTantalum,
                                                          eosTantalum.specificThermalEnergy(rhoTantalum, 300.0)))
    nodesSapphire2.setSpecificThermalEnergy(ScalarField2d("tmp", nodesSapphire2,
                                                          eosSapphire.specificThermalEnergy(rhoSapphire2, 300.0)))
    nodesTungstenCarbide.setSpecificThermalEnergy(ScalarField2d("tmp", nodesTungstenCarbide,
                                                                eosTungstenCarbide.specificThermalEnergy(rhoTungstenCarbide, 300.0)))
    nodesPMMA.setSpecificThermalEnergy(ScalarField2d("tmp", nodesPMMA,
                                                     eosPMMA.specificThermalEnergy(rhoPMMA, 300.0)))

    # Set the node velocities.
    nodesSapphire2.setVelocity(VectorField2d("tmp", nodesSapphire2, v1))
    nodesTungstenCarbide.setVelocity(VectorField2d("tmp", nodesTungstenCarbide, v1))
    nodesPMMA.setVelocity(VectorField2d("tmp", nodesPMMA, v1))

    # Find the interface nodes between the Sapphire 2 and Tungsten.
    dxSapphire1 = Sapphire1Thickness / nxSapphire1
    dxTantalum = TantalumThickness / nxTantalum
    sapphirelist = [i for i in range(nodesSapphire1.numInternalNodes)
                    if nodesSapphire1.positions()[i].x > Sapphire1Range[1][0] - dxSapphire1]
    tantalumlist = [i for i in range(nodesTantalum.numInternalNodes)
                    if nodesTantalum.positions()[i].x > TantalumRange[1][0] - dxTantalum]
    assert mpi.allreduce(len(sapphirelist), mpi.SUM) == nySapphire1
    assert mpi.allreduce(len(tantalumlist), mpi.SUM) == nyTantalum
    interfaceHistory.sapphireIndicies = sapphirelist
    interfaceHistory.tantalumIndicies = tantalumlist

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane2d(Vector2d(xmin, 0.0), Vector2d(1.0, 0.0))
yPlane0 = Plane2d(Vector2d(xmin, ymin), Vector2d(0.0, 1.0))
yPlane1 = Plane2d(Vector2d(xmax, ymax), Vector2d(0.0, -1.0))
xbc0 = ReflectingBoundary2d(xPlane0)
ybc0 = ReflectingBoundary2d(yPlane0)
ybc1 = ReflectingBoundary2d(yPlane1)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
for nodes in nodeLists:
    db.appendNodeList(nodes)
output('db')
output('db.numNodeLists')
output('db.numFluidNodeLists')

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output('q')
output('q.Cl')
output('q.Cq')
output('q.limiter')
output('q.epsilon2')
output('q.negligibleSoundSpeed')
output('q.csMultiplier')
output('q.balsaraShearCorrection')

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro2d(WT, WTPi, q)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = True
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = HsmoothMin
hydro.HsmoothMax = HsmoothMax
hydro.hybridMassDensityThreshold = hybridMassDensityThreshold
output('hydro')
output('hydro.cfl')
output('hydro.useVelocityMagnitudeForDt')
output('hydro.HEvolution')
output('hydro.sumForMassDensity')
output('hydro.HsmoothMin')
output('hydro.HsmoothMax')
output('hydro.kernel()')
output('hydro.PiKernel()')
output('hydro.valid()')
output('hydro.hybridMassDensityThreshold')

#-------------------------------------------------------------------------------
# Construct a strength physics object.
#-------------------------------------------------------------------------------
strength = Strength2d()
output("strength")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator2d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output('integrator')
output('integrator.havePhysicsPackage(hydro)')
output('integrator.havePhysicsPackage(strength)')
output('integrator.valid()')
output('integrator.lastDt')
output('integrator.dtMin')
output('integrator.dtMax')
output('integrator.dtGrowth')

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            boundaryConditions = [xbc0, ybc0, ybc1],
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)
output('control')

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)

    # Viz the initial conditions.
    P = db.fluidPressure
    cs = db.fluidSoundSpeed
    Hi = db.fluidHinverse
    S = SymTensorFieldList2d()
    PS = ScalarFieldList2d()
    for nodes in [nodesTantalum, nodesTungstenCarbide, nodesPMMA]:
        S.appendField(nodes.deviatoricStress())
        PS.appendField(nodes.plasticStrain())
    dumper = SpheralVisitDump(db,
                              "PlateImpact-2d-visit",
                              "dumps-2d",
                              listOfFieldLists = [db.fluidMassDensity,
                                                  db.fluidVelocity,
                                                  db.fluidWeight,
                                                  P,
                                                  cs,
                                                  Hi,
                                                  S,
                                                  PS]
                              )
    dumper.dump(control.time(), control.totalSteps)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    interfaceHistory.sample(control.time())
    control.dropRestartFile()

    # Viz the current state.
    P = db.fluidPressure
    cs = db.fluidSoundSpeed
    Hi = db.fluidHinverse
    S = SymTensorFieldList2d()
    PS = ScalarFieldList2d()
    for nodes in [nodesTantalum, nodesTungstenCarbide, nodesPMMA]:
        S.appendField(nodes.deviatoricStress())
        PS.appendField(nodes.plasticStrain())
    dumper = SpheralVisitDump(db,
                              "PlateImpact-2d-visit",
                              "dumps-2d",
                              listOfFieldLists = [db.fluidMassDensity,
                                                  db.fluidVelocity,
                                                  db.fluidWeight,
                                                  P,
                                                  cs,
                                                  Hi,
                                                  S,
                                                  PS]
                              )
    dumper.dump(control.time(), control.totalSteps)
