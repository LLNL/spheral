#-------------------------------------------------------------------------------
# A simple periodic boundary test -- just loop a gas around a box.
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral1d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *

import mpi
import DistributeNodes

title("PeriodicBoundary test -- 1D")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 20,
            rho1 = 1.0,
            vx1 = 1.0,
            P1 = 1.0,

            gamma = 5.0/3.0,
            mu = 1.0,

            nPerh = 2.01,

            HydroConstructor = ASPHHydro,
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            Cl = 1.0, 
            Cq = 0.75,
            Qlimiter = False,
            balsaraCorrection = True,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 2.0,
            steps = None,
            vizCycle = None,
            vizTime = 0.01,
            dt = 0.0001,
            dtMin = 1.0e-8, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = IdealH,
            domainIndependent = False,
            rigorousBoundaries = True,
            dtverbose = False,

            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = False,           # <--- Important!  rigorousBoundaries does not work with the compatibleEnergy algorithm currently.
            gradhCorrection = False,

            restoreCycle = None,
            restartStep = 1000000,
            )


#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
vel = nodes1.velocity()
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (0.0, 1.0))],
                         nPerh = nPerh)

# Set specific thermal energies
eps1 = P1/((gamma - 1.0)*rho1)
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

# Set node velocities
nodes1.velocity(VectorField("tmp", nodes1, Vector(vx1, 0.0)))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
db.appendNodeList(nodes1)

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = HydroConstructor(W = WT,
                         Q = q,
                         cfl = cfl,
                         compatibleEnergyEvolution = compatibleEnergy,
                         gradhCorrection = gradhCorrection,
                         XSPH = XSPH,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)
packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xp1 = Plane(Vector(0.0), Vector( 1.0))
xp2 = Plane(Vector(1.0), Vector(-1.0))
xbc = PeriodicBoundary(xp1, xp2)

for p in packages:
    p.appendBoundary(xbc)

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

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT)
control.iterateIdealH(hydro)
if densityUpdate in (VoronoiCellDensity, SumVoronoiCellDensity):
    control.voronoiInitializeMass()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)
