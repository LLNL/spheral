#-------------------------------------------------------------------------------
# The Spherical Sedov test case (2-D).
#-------------------------------------------------------------------------------
import os, sys, shutil
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from CubicNodeGenerator import GenerateSquareNodeDistribution

import mpi

title("2-D integrated hydro test -- planar Sedov problem")
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            nx = 50,
            ny = 50,
            nz = 50,
            rmin = 0.0,
            rmax = 1.0,
            nPerh = 2.01,
            nxparts = 1,   # Optionally split the cube into (nxparts**3 NodeLists)

            rho0 = 1.0,
            eps0 = 0.0,
            Espike = 1.0,
            smoothSpike = False,
            gamma = 5.0/3.0,
            mu = 1.0,

            Cl = 1.0,
            Cq = 0.75,
            epsilon2 = 1e-2,

            HydroConstructor = SPHHydro,
            hmin = 1e-15,
            hmax = 1.0,
            cfl = 0.5,
            useVelocityMagnitudeForDt = True,
            XSPH = True,
            rhomin = 1e-10,

            steps = None,
            goalTime = None,
            goalRadius = 0.8,
            dt = 1e-8,
            dtMin = 1.0e-8,
            dtMax = None,
            dtGrowth = 2.0,
            vizCycle = None,
            vizTime = 0.1,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = IdealH,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = False,

            restoreCycle = None,
            restartStep = 1000,

            clearDirectories = False,
            dataRoot = "dumps-spherical-Sedov",
            graphics = True,

            )

# Figure out what our goal time should be.
import SedovAnalyticSolution
h0 = 1.0/nx*nPerh
answer = SedovAnalyticSolution.SedovSolution(nDim = 3,
                                             gamma = gamma,
                                             rho0 = rho0,
                                             E0 = Espike,
                                             h0 = h0)
if goalTime is None:
    assert not goalRadius is None
    nu1 = 1.0/(answer.nu + 2.0)
    nu2 = 2.0*nu1
    goalTime = (goalRadius*(answer.alpha*rho0/Espike)**nu1)**(1.0/nu2)
vs, r2, v2, rho2, P2 = answer.shockState(goalTime)
print "Predicted shock position %g at goal time %g." % (r2, goalTime)

# Scale the spike energy according to the boundary conditions we're using.
Espike /= 8.0

dataDir = os.path.join(dataRoot,
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "gradhCorrection=%s" % gradhCorrection,
                       "seed=%s" % seed,
                       "nxparts=%s" % nxparts,
                       "nx=%i_ny=%i_nz=%i" % (nx, ny, nz))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Sedov-spherical-3d-%ix%ix%i" % (nx, ny, nz))

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
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
# Create the NodeLists and generators.
#-------------------------------------------------------------------------------
nodeSet, gens = [], []
if restoreCycle is None:
    for iz in xrange(nxparts):
        zmin = float(iz    )/nxparts
        zmax = float(iz + 1)/nxparts
        for iy in xrange(nxparts):
            ymin = float(iy    )/nxparts
            ymax = float(iy + 1)/nxparts
            for ix in xrange(nxparts):
                xmin = float(ix    )/nxparts
                xmax = float(ix + 1)/nxparts
                nodeSet.append(makeFluidNodeList("nodes%i" % (ix + nxparts*(iy + nxparts*iz)),
                                                 eos, 
                                                 hmin = hmin,
                                                 hmax = hmax,
                                                 nPerh = nPerh,
                                                 rhoMin = rhomin))
                gens.append(GenerateNodeDistribution3d(nx/nxparts, ny/nxparts, nz/nxparts,
                                                       rho0, seed,
                                                       xmin = (xmin, ymin, zmin),
                                                       xmax = (xmax, ymax, zmax),
                                                       nNodePerh = nPerh,
                                                       SPH = (HydroConstructor == SPHHydro)))

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
        #from PeanoHilbertDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d
    distributeNodes3d(*zip(nodeSet, gens))

    # Set the point source of energy.
    Esum = 0.0
    if smoothSpike:
        Wsum = 0.0
        for nodes in nodeSet:
            pos = nodes.positions()
            vel = nodes.velocity()
            mass = nodes.mass()
            eps = nodes.specificThermalEnergy()
            for nodeID in xrange(nodes.numInternalNodes):
                Hi = H[nodeID]
                etaij = (Hi*pos[nodeID]).magnitude()
                Wi = WT.kernelValue(etaij, Hi.Determinant())
                Ei = Wi*Espike
                epsi = Ei/mass[nodeID]
                eps[nodeID] = epsi
                Wsum += Wi
        Wsum = mpi.allreduce(Wsum, mpi.SUM)
        assert Wsum > 0.0
        for nodes in nodeSet:
            eps = nodes.specificThermalEnergy()
            for nodeID in xrange(nodes.numInternalNodes):
                eps[nodeID] /= Wsum
                Esum += eps[nodeID]*mass[nodeID]
    else:
        i = -1
        rmin = 1e50
        for nodes in nodeSet:
            pos = nodes.positions()
            for nodeID in xrange(nodes.numInternalNodes):
                rij = pos[nodeID].magnitude()
                if rij < rmin:
                    i = nodeID
                    nodesmin = nodes
                    rmin = rij
        rminglobal = mpi.allreduce(rmin, mpi.MIN)
        if fuzzyEqual(rmin, rminglobal):
            mass = nodesmin.mass()
            eps = nodesmin.specificThermalEnergy()
            assert i >= 0 and i < nodesmin.numInternalNodes
            eps[i] = Espike/mass[i]
            Esum += Espike
    Eglobal = mpi.allreduce(Esum, mpi.SUM)
    print "Initialized a total energy of", Eglobal
    assert smoothSpike or fuzzyEqual(Eglobal, Espike)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    db.appendNodeList(nodes)
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct a standard Monaghan-Gingold artificial viscosity.
#-------------------------------------------------------------------------------
qMG = MonaghanGingoldViscosity(Cl, Cq)
qMG.epsilon2 = epsilon2
output("qMG")
output("qMG.Cl")
output("qMG.Cq")
output("qMG.epsilon2")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = HydroConstructor(WT, WTPi, qMG,
                         cfl = cfl,
                         compatibleEnergyEvolution = compatibleEnergy,
                         gradhCorrection = gradhCorrection,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.XSPH")
output("hydro.sumForMassDensity")
output("hydro.HEvolution")
output("hydro.epsilonTensile")
output("hydro.nTensile")

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(0, 0, 0), Vector(1, 0, 0))
yPlane0 = Plane(Vector(0, 0, 0), Vector(0, 1, 0))
zPlane0 = Plane(Vector(0, 0, 0), Vector(0, 0, 1))
xbc0 = ReflectingBoundary(xPlane0)
ybc0 = ReflectingBoundary(yPlane0)
zbc0 = ReflectingBoundary(zPlane0)

for p in packages:
    for bc in (xbc0, ybc0, zbc0):
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the one physics package.
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
                            restoreCycle = restoreCycle,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            vizBaseName = "Sedov-spherical-3d-%ix%ix%i" % (nx, ny, nz),
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
output("control")

#-------------------------------------------------------------------------------
# Restart if we're doing it.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    if densityUpdate in (VoronoiCellDensity, SumVoronoiCellDensity):
        print "Reinitializing node masses."
        control.voronoiInitializeMass()

#-------------------------------------------------------------------------------
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------
if steps is None:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()
else:
    control.step(steps)

# Output the energy conservation.
print "Energy conservation: ", ((control.conserve.EHistory[-1] -
                                 control.conserve.EHistory[0])/
                                control.conserve.EHistory[0])

#-------------------------------------------------------------------------------
# Generate some error metrics comparing to the analytic solution.
#-------------------------------------------------------------------------------
# Report the error norms.
rmin, rmax = 0.0, 0.95
pos_local, rho_local, vel_local, eps_local, P_local = [], [], [], [], []
for nodes in nodeSet:
    pos_local += list(nodes.positions().internalValues())
    rho_local += list(nodes.massDensity().internalValues())
    vel_local += list(nodes.velocity().internalValues())
    eps_local += list(nodes.specificThermalEnergy().internalValues())
    Pf = ScalarField("pressure", nodes)
    nodes.pressure(Pf)
    P_local += list(Pf.internalValues())

r = mpi.allreduce([x.magnitude() for x in pos_local], mpi.SUM)
rho = mpi.allreduce(rho_local, mpi.SUM)
v = mpi.allreduce([x.magnitude() for x in vel_local], mpi.SUM)
eps = mpi.allreduce(eps_local, mpi.SUM)
P = mpi.allreduce(P_local, mpi.SUM)
if mpi.rank == 0:
    from SpheralGnuPlotUtilities import multiSort
    import Pnorm
    multiSort(r, rho, v, eps, P)
    rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), r)
    print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
    for (name, data, ans) in [("Mass Density", rho, rhoans),
                              ("Pressure", P, Pans),
                              ("Velocity", v, vans),
                              ("Thermal E", eps, epsans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in xrange(len(data))]
        Pn = Pnorm.Pnorm(error, r)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print "\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf)

#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
if graphics:
    import Gnuplot

    # Plot the final state.
    rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    del HPlot
    Hinverse = db.newFluidSymTensorFieldList()
    db.fluidHinverse(Hinverse)
    hr = db.newFluidScalarFieldList()
    for Hfield, hrfield, in zip(Hinverse, hr):
        n = Hfield.numElements
        assert hrfield.numElements == n
        positions = Hfield.nodeList().positions()
        for i in xrange(n):
            runit = positions[i].unitVector()
            hrfield[i] = (Hfield[i]*runit).magnitude()
    hrPlot = plotFieldList(hr, xFunction="%s.magnitude()", plotStyle="points", winTitle="h_r")

    # Overplot the analytic solution.
    plotAnswer(answer, control.time(),
               rhoPlot = rhoPlot,
               velPlot = vrPlot,
               epsPlot = epsPlot,
               PPlot = PPlot,
               HPlot = hrPlot)

    # Compute the simulated specific entropy.
    A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

    # The analytic solution for the simulated entropy.
    xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), r)
    Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]

    # Plot the specific entropy.
    if mpi.rank == 0:
        AsimData = Gnuplot.Data(xprof, A,
                                with_ = "points",
                                title = "Simulation",
                                inline = True)
        AansData = Gnuplot.Data(xprof, Aans,
                                with_ = "lines",
                                title = "Solution",
                                inline = True)
        Aplot = Gnuplot.Gnuplot()
        Aplot.plot(AsimData)
        Aplot.replot(AansData)
        Aplot.title("Specific entropy")
        Aplot.refresh()
    else:
        Aplot = fakeGnuplot()
