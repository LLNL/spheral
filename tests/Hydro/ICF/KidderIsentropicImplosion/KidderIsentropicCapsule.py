#===============================================================================
# This is a version of the Kidder isentropic capsule implosion test.  This problem
# is a good test of how well the adiabat is maintained during an implosion.  Since
# this is an isentropic (shockless) problem, it can be run without the artificial
# viscosity to see how well we can possibly do, as well as with the artificial
# viscosity to see if it is perturbing the problem.
#
# See references:
#
# R.E. Kidder, Nucl. Fusion 16 (1976), 3-14.
# Maire, JCP 228 (2009), 6882-6915.
#===============================================================================
from Spheral1d import *
from SpheralTestUtilities import *
from math import *
import os
import sys, string
from numericalIntegration import *

from KidderIsentropicCapsuleAnalyticSolution import *
from KidderIsentropicCapsuleBoundary import *

#-------------------------------------------------------------------------------
# Parameters for the run.
#-------------------------------------------------------------------------------
commandLine(problemName = "KidderIsentropicCapsule",
            NodeListConstructor = SphNodeList,
            KernelConstructor = BSplineKernel,
            IntegratorConstructor = CheapSynchronousRK2Integrator,

            # Timing
            steps = None,
            goalTau = 0.99,       # Fraction of the dimensionless time (tau) we
                                  # will advance to (tau=1 => full collapse)

            # Initial conditions
            r0 = 0.9,             # Inner boundary radius
            r1 = 1.0,             # Outer boundary radius
            P0 = 0.1,             # Inner boundary pressure
            P1 = 10.0,            # Outer boundary pressure
            rho1 = 0.01,          # Outer boundary density

            # Resolution
            nr = 100,             # num radial points
            nrGhost = 5,          # how deep do we want the Boundaries
            nPerh = 2.01,

            # Output
            dumpCycle = 0,
            dumpTime = 0.005,
            dt = 1.0e-6,
            dtMin = 1.0e-8,
            dtMax = 1.0,
            dtGrowth = 2.0,
            rigorousBoundaries = True,
            domainIndependent = False,
            statsStep = 1,
            restartStep = 10000,
            restoreCycle = None,
            smoothIters = 0,
            maxSteps = None,

            # LagrangeHydro
            Qconstructor = MonaghanGingoldViscosity,
            Cq = 1.0,             # Default to zero Q, since this is a shockless problem
            Cl = 1.0,
            Qlimiter = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,
            XSPH = True,
            compatibleEnergy = True,
            gradhCorrection = False,
            sumForMassDensity = RigorousSumDensity,
            HEvolution = IdealH,

            neighborSearchType = GatherScatter,
            numGridLevels = 20,
            topGridCellSize = 2.0,
            origin = Vector1d(0.0),

            profileASCII = False, # Optionally spew the profiles to an ASCII file
            )

# The dimensionality of the problem: 1 => planar
#                                    2 => cylindrical
#                                    3 => spherical
nu = 1

# Construct the analytic solution for this set up.
answer = KidderIsentropicCapsuleAnalyticSolution(nu, r0, r1, P0, P1, rho1)
goalTime = goalTau * answer.tau
print "Capsule collapses at %g, goal time is %g." % (answer.tau, goalTime)

problemName = "%s-%i" % (problemName, nr)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
gamma = answer.gamma
mu = 1.0
eos = GammaLawGasMKS(answer.gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(KernelConstructor(), 1000)
WTPi = TableKernel(KernelConstructor(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = NodeListConstructor("nodes", eos, WT, WTPi)
nodes.XSPH = XSPH
nodes.hmin = hmin
nodes.hmax = hmax
nodes.nodesPerSmoothingScale = nPerh
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")
output("nodes.XSPH")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor = NestedGridNeighbor(nodes,
                              neighborSearchType,
                              numGridLevels,
                              topGridCellSize,
                              origin,
                              kernelExtent)
nodes.registerNeighbor(neighbor)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes, nr, rho1, (r0, r1))],
                         nPerh = nPerh)
output("nodes.numNodes")

# Set the initial conditions.
dr = (r1 - r0)/nr
pos = nodes.positions()
mass = nodes.mass()
rho = nodes.massDensity()
eps = nodes.specificThermalEnergy()
for i in xrange(nodes.numInternalNodes):
    ri = pos[i].x
    mi = trapezoidalIntegration(answer.rhoInitial, ri - 0.5*dr, ri + 0.5*dr, 200)
    rho[i] = mi/dr
    mass[i] = mi
    eps[i] = trapezoidalIntegration(answer.Pinitial, ri - 0.5*dr, ri + 0.5*dr, 200)/((gamma - 1.0)*mi)
    #eps[i] = answer.P(0.0, ri)/((gamma - 1.0)*rho[i])

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro(W = WT,
              Q = q,
              compatibleEnergyEvolution = compatibleEnergy,
              gradhCorrection = gradhCorrection,
              densityUpdate = sumForMassDensity,
              HUpdate = HEvolution)
hydro.cfl = cfl
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.hmin")
output("hydro.hmax")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.domainDecompositionIndependent = domainIndependent
output("integrator")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.domainDecompositionIndependent")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
# Find the nrGhost nodes closest to the inner boundary.
mindicesFL = mortonOrderIndices(db)
mindices = mindicesFL[0]
positions = nodes.positions()
nodeSet = mpi.allreduce([(positions[i].x, mindices[i]) for i in xrange(nodes.numInternalNodes)],
                        mpi.SUM)
nodeSet.sort()

innerIndices = [tup[1] for tup in nodeSet[:nrGhost]]
innerNodes = [i for i in xrange(nodes.numInternalNodes) if mindices[i] in innerIndices]
assert mpi.allreduce(len(innerNodes), mpi.SUM) == nrGhost

outerIndices = [tup[1] for tup in nodeSet[-nrGhost:]]
outerNodes = [i for i in xrange(nodes.numInternalNodes) if mindices[i] in outerIndices]
assert mpi.allreduce(len(outerNodes), mpi.SUM) == nrGhost

interiorNodes = [i for i in xrange(nodes.numInternalNodes)
                 if ((i not in innerIndices) and
                     (i not in outerIndices))]

h0 = nPerh * dr

rbc0 = KidderIsentropicCapsuleEnforcementBoundary1d(integrator = integrator,
                                                    answer = answer,
                                                    nodeList = nodes,
                                                    nodeIDs = innerNodes,
                                                    interiorNodeIDs = interiorNodes,
                                                    hinitial = h0,
                                                    hmin = hmin,
                                                    hmax = hmax,
                                                    dr0 = dr)
# rbc1 = KidderIsentropicCapsuleEnforcementBoundary1d(integrator = integrator,
#                                                     answer = answer,
#                                                     nodeList = nodes,
#                                                     nodeIDs = outerNodes,
#                                                     interiorNodeIDs = interiorNodes,
#                                                     hinitial = h0,
#                                                     hmin = hmin,
#                                                     hmax = hmax,
#                                                     dr0 = dr)

rbc1 = KidderIsentropicCapsuleBoundary1d(innerBoundary = False,
                                         integrator = integrator,
                                         answer = answer,
                                         nodeList = nodes,
                                         nrGhostNodes = nrGhost,
                                         dr0 = dr)

## hydro.appendBoundary(rbc0)
hydro.appendBoundary(rbc1)

integrator.appendPhysicsPackage(rbc0)
# integrator.appendPhysicsPackage(rbc1)

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = problemName,
                            initializeMassDensity = (sumForMassDensity != IntegrateDensity),
                            initializeDerivatives = True)
output("control")

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH()
    control.smoothState(smoothIters)

## state = State(db, integrator.physicsPackages())
## derivs = StateDerivatives(db, integrator.physicsPackages())
## integrator.initialize(state, derivs)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
print "Starting energy (measured, analytic, error): ", (control.conserve.EHistory[0],
                                                        answer.totalEnergy(0.0),
                                                        (control.conserve.EHistory[0] - answer.totalEnergy(0.0))/
                                                        answer.totalEnergy(0.0))

if not steps is None:
    control.step(steps)
else:
    if control.time() < goalTime:
        control.step(5)
        control.advance(goalTime, maxSteps)

print "Final energy (measured, analytic, error): ", (control.conserve.EHistory[-1],
                                                     answer.totalEnergy(goalTime),
                                                     (control.conserve.EHistory[-1] - answer.totalEnergy(0.0))/
                                                     answer.totalEnergy(goalTime))


#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
import Gnuplot

# Simulation results.
r = mpi.allreduce([x.x for x in nodes.positions().internalValues()], mpi.SUM)
rho = mpi.allreduce(nodes.massDensity().internalValues(), mpi.SUM)
P = mpi.allreduce(nodes.pressure().internalValues(), mpi.SUM)
v = mpi.allreduce([v.x for v in nodes.velocity().internalValues()], mpi.SUM)
eps = mpi.allreduce(nodes.specificThermalEnergy().internalValues(), mpi.SUM)
S = [p/d**gamma for p, d in zip(P, rho)]

# Analytic results.
t = control.time()
rhoAns = [answer.rho(t, ri) for ri in r]
Pans = [answer.P(t, ri) for ri in r]
vAns = [answer.vr(t, ri) for ri in r]
epsAns = [answer.eps(t, ri) for ri in r]
Sans = [answer.S for ri in r]

# The ratio of the entropy to the expected value.
alpha = [ss/sa for ss, sa in zip(S, Sans)]

# Now plot the suckers.
if mpi.rank == 0:
    rhoData = Gnuplot.Data(r, rho, title="Mass Density", with_="points", inline=True)
    PData = Gnuplot.Data(r, P, title="Pressure", with_="points", inline=True)
    SData = Gnuplot.Data(r, S, title="Entropy", with_="points", inline=True)
    vrData = Gnuplot.Data(r, v, title="Velocity (radial)", with_="points", inline=True)
    epsData = Gnuplot.Data(r, eps, title="Specific thermal energy", with_="points", inline=True)

    rhoAnsData = Gnuplot.Data(r, rhoAns, title="Solution", with_="lines", inline=True)
    PansData = Gnuplot.Data(r, Pans, title="Solution", with_="lines", inline=True)
    SansData = Gnuplot.Data(r, Sans, title="Solution", with_="lines", inline=True)
    vrAnsData = Gnuplot.Data(r, vAns, title="Solution", with_="lines", inline=True)
    epsAnsData = Gnuplot.Data(r, epsAns, title="Solution", with_="lines", inline=True)

    alphaData = Gnuplot.Data(r, alpha, title="Sim/Expected entropy", inline=True)
    alphaAnsData = Gnuplot.Data(r, [1.0 for x in r], title="Solution", with_="lines", inline=True)

    rhoPlot = Gnuplot.Gnuplot()
    rhoPlot.plot(rhoData)
    rhoPlot.replot(rhoAnsData)

    Pplot = Gnuplot.Gnuplot()
    Pplot.plot(PData)
    Pplot.replot(PansData)

    Splot = Gnuplot.Gnuplot()
    Splot.plot(SData)
    Splot.replot(SansData)

    vrPlot = Gnuplot.Gnuplot()
    vrPlot.plot(vrData)
    vrPlot.replot(vrAnsData)

    epsPlot = Gnuplot.Gnuplot()
    epsPlot.plot(epsData)
    epsPlot.replot(epsAnsData)

    alphaPlot = Gnuplot.Gnuplot()
    alphaPlot.plot(alphaData)
    alphaPlot.replot(alphaAnsData)

    # If requested, output the profiles to an ASCII file.
    if profileASCII:
        f = open(problemName + "_profiles.txt", "w")
        f.write((11*"%20s" + "\n") % ("# radius",
                                      "rad velocity (sim)",
                                      "mass density (sim)",
                                      "pressure (sim)",
                                      "eps (sim)",
                                      "entropy (sim)",
                                      "rad velocity (ans)",
                                      "mass density (ans)",
                                      "pressure (ans)",
                                      "eps (ans)",
                                      "entropy (ans)"))

        # Now write the suckers out.
        for tup in zip(r, v, rho, P, eps, S, vAns, rhoAns, Pans, epsAns, Sans):
            l = ""
            for val in tup:
                if val is None:
                    l += 20*" "
                else:
                    l += "%15.10g     " % val
            l += "\n"
            f.write(l)

        f.close()
