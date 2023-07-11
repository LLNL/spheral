import os, sys, shutil, mpi
from math import *
from Spheral2d import *
from GenerateNodeDistribution2d import *
from SpheralTestUtilities import *

title("Sedov 2d Test")
commandLine(nRadial = 50,	# number of radial bins/particles
            nTheta = 50,	# number of theta bins/particles
            rmin = 0.0,		# minimum radius
            rmax = 1.0,		# maximum radius
            nPerh = 1.51,	# particle neighbor count
            
            rho0 = 1.0,
            eps0 = 0.0,
            gamma = 5.0/3.0,
            mu = 1.0,
            Espike = 1.0,
            
            KernelConstructor = NBSplineKernel,
            order = 5,
            
            Cl = 1.0,
            Cq = 0.75,
            hmin = 1e-15,
            hmax = 1.0,
            cfl = 0.5,
            rhomin = 1e-10,
            
            HydroConstructor = SPHHydro,
            Qconstructor = MonaghanGingoldViscosity,
            balsaraCorrection = False,
            linearInExpansion = False,
            densityUpdate = RigorousSumDensity,
            HEvolution = IdealH,
            compatibleEnergy = True,
            
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            steps = None,
            goalTime = 1.0,
            dt = 1e-8,
            dtMin = 1e-8,
            dtMax = None,
            dtGrowth = 2.0,
            restoreCycle = None,
            restartStep = 1000,
            
            dataDir = "dumps-sedov2d")

dataDir = os.path.join(dataDir,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       "nperh=%4.2f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "nr=%i_nt=%i" % (nRadial, nTheta))

restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "sedov-2d-%i" % nRadial)

if mpi.rank==0:
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
pos = nodes1.positions()
vel = nodes1.velocity()
mass = nodes1.mass()
eps = nodes1.specificThermalEnergy()
H = nodes1.Hfield()
    
xmin = (-1.0, -1.0)
xmax = (1.0, 1.0)

if restoreCycle is None:
    generator = GenerateNodeDistribution2d(nRadial, nTheta, rho0, "lattice",
                                           rmin = rmin,
                                           rmax = rmax,
                                           xmin = xmin,
                                           xmax = xmax,
                                           theta = 2.0*pi,
                                           azimuthalOffsetFraction = 0.0,
                                           nNodePerh = nPerh,
                                           SPH = True)
    
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes1, generator))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set the point source of energy.
    Esum = 0.0
    i = -1
    rmin = 1e50
    for nodeID in range(nodes1.numInternalNodes):
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
    print("Initialized a total energy of", Eglobal)

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
q.epsilon2 = 1e-2
q.limiter = False
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
#------------------------------------------------------------------------------
hydro = HydroConstructor(W = WT,
                         Q = q,
                         cfl = cfl,
                         compatibleEnergyEvolution = compatibleEnergy,
                         evolveTotalEnergy = False,
                         gradhCorrection = True,
                         correctVelocityGradient = True,
                         densityUpdate = densityUpdate,
                         XSPH = False,
                         HUpdate = HEvolution)

output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

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
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizMethod = None,
                            vizBaseName = "Sedov-2d-%ix%i" % (nRadial, nTheta),
                            vizDir = vizDir,
                            vizTime = 10.0,
                            SPH = True)
output("control")

#-------------------------------------------------------------------------------
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------
if steps is None:
    control.advance(goalTime)
else:
    control.step(steps)

# Output the energy conservation.
print("Energy conservation: ", ((control.conserve.EHistory[-1] -
                                 control.conserve.EHistory[0])/
                                control.conserve.EHistory[0]))
