#-------------------------------------------------------------------------------
# The planar Sedov test case (1-D).
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
import os, sys

import mpi

title("1-D integrated hydro test -- planar Sedov problem")
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 101,
            rho1 = 1.0,
            eps1 = 0.0,
            x0 = -1.0,
            x1 = 1.0,
            nPerh = 2.01,

            xSpike = 0.0,
            Espike = 1.0,
            smoothSpike = False,
            gamma = 5.0/3.0,
            mu = 1.0,

            Cl = 1.0,
            Cq = 0.75,
            epsilon2 = 1e-2,
            
            CRKSPH = False,
            Qconstructor = MonaghanGingoldViscosity,
            momentumConserving = True, # For CRKSPH
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            HUpdate = IdealH,
            filter = 0.0,

            HydroConstructor = SPHHydro,
            hmin = 1e-15,
            hmax = 1.0,
            cfl = 0.5,
            useVelocityMagnitudeForDt = True,
            XSPH = True,
            rhomin = 1e-10,

            goalTime = 0.3,
            dt = 1e-8,
            dtMin = 1.0e-8,
            dtMax = None,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = IdealH,
            sumForMassDensity = RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = False,

            restoreCycle = None,
            restartStep = 1000,

            clearDirectories = True,
            dataRoot = "dump-planar",
            serialDump = False, #whether to dump a serial ascii file at the end for viz
            )

dataDir = os.path.join(dataRoot,
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "gradhCorrection=%s" % gradhCorrection,
                       "n=%i" % nx1)
restartBaseName = os.path.join(dataDir, "Sedov-planar-1d-%i" % nx1)

dx = (x1 - x0)/nx1
H1 = SymTensor(1.0/(nPerh*dx))

#-------------------------------------------------------------------------------
# CRKSPH Switches to ensure consistency
#-------------------------------------------------------------------------------
if CRKSPH:
    Qconstructor = CRKSPHMonaghanGingoldViscosity

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
from findLastRestart import findLastRestart
restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = WT
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Create a NodeList and associated Neighbor object.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           rhoMin = rhomin)

#-------------------------------------------------------------------------------
# Set the NodeList properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
list = [(nodes1, nx1, rho1, (x0, x1))]
distributeNodesInRange1d(list)
output("mpi.allreduce(nodes1.numNodes, mpi.MIN)")
output("mpi.allreduce(nodes1.numNodes, mpi.MAX)")
output("mpi.allreduce(nodes1.numNodes, mpi.SUM)")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

# Set the point source of energy.
Esum = 0.0
if smoothSpike:
    Wsum = 0.0
    for nodeID in xrange(nodes1.numInternalNodes):
        Hi = nodes1.Hfield()[nodeID].xx
        rij = abs(nodes1.positions()[nodeID].x - xSpike)*Hi
        Wi = WT.kernelValue(rij, Hi)
        Ei = Wi*Espike
        epsi = Ei/nodes1.mass()[nodeID]
        nodes1.specificThermalEnergy()[nodeID] = epsi
        Wsum += Wi
    Wsum = mpi.allreduce(Wsum, mpi.SUM)
    assert Wsum > 0.0
    for nodeID in xrange(nodes1.numInternalNodes):
        nodes1.specificThermalEnergy()[nodeID] /= Wsum
        Esum += nodes1.specificThermalEnergy()[nodeID]*nodes1.mass()[nodeID]
else:
    i = -1
    rmin = 1e50
    for nodeID in xrange(nodes1.numInternalNodes):
        rij = abs(nodes1.positions()[nodeID].x - xSpike)
        if rij < rmin:
            i = nodeID
            rmin = rij
    rminglobal = mpi.allreduce(rmin, mpi.MIN)
    if fuzzyEqual(rmin, rminglobal):
        assert i >= 0 and i < nodes1.numInternalNodes
        nodes1.specificThermalEnergy()[i] = Espike/nodes1.mass()[i]
        Esum += Espike
Eglobal = mpi.allreduce(Esum, mpi.SUM)
print "Initialized a total energy of", Eglobal
assert smoothSpike or fuzzyEqual(Eglobal, Espike)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct a standard Monaghan-Gingold artificial viscosity.
#-------------------------------------------------------------------------------
qMG = Qconstructor(Cl, Cq)
qMG.epsilon2 = epsilon2
output("qMG")
output("qMG.Cl")
output("qMG.Cq")
output("qMG.epsilon2")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if CRKSPH:
    hydro = CRKSPHHydro(WT, WTPi, qMG,
                      filter = filter,
                      cfl = cfl,
                      compatibleEnergyEvolution = compatibleEnergy,
                      XSPH = XSPH,
                      densityUpdate = densityUpdate,
                      HUpdate = HUpdate,
                      momentumConserving = momentumConserving)
else:
    hydro = SPHHydro(WT, WTPi, qMG,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             densityUpdate = sumForMassDensity,
                             XSPH = XSPH,
                             HUpdate = HEvolution)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
#output("hydro.gradhCorrection")
output("hydro.XSPH")
#output("hydro.sumForMassDensity")
output("hydro.HEvolution")
#output("hydro.epsilonTensile")
#output("hydro.nTensile")

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)
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
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output("control")

#-------------------------------------------------------------------------------
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------
if restoreCycle:
    control.loadRestartFile(restoreCycle)
else:
    #control.iterateIdealH(hydro)
##     db.updateFluidMassDensity()

##     # This bit of craziness is needed to try and get the intial mass density
##     # as exactly as possible.
##     for i in xrange(10):
##         thpt = mpi.allreduce(max(nodes1.massDensity().internalValues()), mpi.MAX)
##         print i, thpt
##         for j in xrange(nodes1.numInternalNodes):
##             nodes1.mass()[j] *= rho1/thpt
##         state = State(db, integrator.physicsPackages())
##         derivs = StateDerivatives(db, integrator.physicsPackages())
##         integrator.initialize(state, derivs)
##         db.updateFluidMassDensity()

    control.smoothState(smoothIters)
    control.dropRestartFile()

# Advance to the end time.
control.advance(goalTime, maxSteps)
control.dropRestartFile()

# Output the energy conservation.
print "Energy conservation: ", ((control.conserve.EHistory[-1] -
                                 control.conserve.EHistory[0])/
                                control.conserve.EHistory[0])

# Plot the final state.
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotStyle="linespoints")

# Overplot the analytic solution.
import SedovAnalyticSolution
h1 = (x1 - x0)/nx1*nPerh
answer = SedovAnalyticSolution.SedovSolution(nDim = 1,
                                             gamma = gamma,
                                             rho0 = rho1,
                                             E0 = Espike,
                                             h0 = h1)
xprof = mpi.allreduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
plotAnswer(answer, control.time(),
           rhoPlot, velPlot, epsPlot, PPlot, HPlot, xprof)

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

# Plot the specific entropy.
AsimData = Gnuplot.Data(xprof, A,
                        with_ = "points",
                        title = "Simulation",
                        inline = True)
AansData = Gnuplot.Data(xprof, Aans,
                        with_ = "lines",
                        title = "Solution",
                        inline = True)
    
Aplot = generateNewGnuPlot()
Aplot.plot(AsimData)
Aplot.replot(AansData)
Aplot.title("Specific entropy")
Aplot.refresh()

if serialDump:
    serialData = []
    i,j = 0,0
    
    f = open(dataDir + "/sedov-planar-1d-CRKSPH-" + str(CRKSPH) + ".ascii",'w')
    f.write("i x m rho u v rhoans uans vans\n")
    for j in xrange(nodes1.numInternalNodes):
        f.write("{0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(j,nodes1.positions()[j][0],
                                                               nodes1.mass()[j],
                                                               nodes1.massDensity()[j],
                                                               nodes1.specificThermalEnergy()[j],
                                                               nodes1.velocity()[j][0],
                                                               rhoans[j],
                                                               uans[j],
                                                               vans[j]))
    f.close()
