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
from SolidSpheral1d import *
from SpheralTestUtilities import *
from math import *
import os
import sys, string, shutil
from numericalIntegration import *

from KidderIsentropicCapsuleAnalyticSolution import *
from KidderIsentropicCapsuleBoundary import *

#-------------------------------------------------------------------------------
# Parameters for the run.
#-------------------------------------------------------------------------------
commandLine(problemName = "KidderIsentropicCapsule",
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            kernelOrder = 5,

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
            nrGhost = 10,         # how deep do we want the Boundaries
            nPerh = 1.51,

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
            crksph = False,
            psph = False,
            solid = False,
            XSPH = True,
            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,
            filter = 0.0,
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,
            compatibleEnergy = True,
            evolveTotalEnergy = False,
            gradhCorrection = True,
            correctVelocityGradient = True,
            cullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            hopkinsCullenCorrection = True,
            HopkinsConductivity = False,
            densityUpdate = RigorousSumDensity,
            HUpdate = IdealH,

            clearDirectories = True,
            dataDirBase = "dumps-Kidder-planar",
            profileASCII = False, # Optionally spew the profiles to an ASCII file
            )

# Choose our hydro object.
if crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
else:
    hydroname = "SPH"

# The dimensionality of the problem: 1 => planar
#                                    2 => cylindrical
#                                    3 => spherical
# Construct the analytic solution for this set up.
answer = KidderIsentropicCapsuleAnalyticSolution(1, r0, r1, P0, P1, rho1)
goalTime = goalTau * answer.tau
print "Capsule collapses at %g, goal time is %g." % (answer.tau, goalTime)

dataDir = os.path.join(dataDirBase, 
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "evolveTotalEnergy=%s" % evolveTotalEnergy,
                       "Cullen=%s" % cullenViscosity,
                       "filter=%f" % filter,
                       "nr=%i" % nr)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Kidder-planar-1d")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
mu = 1.0
eos = GammaLawGasMKS(answer.gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(kernelOrder), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    nodes = makeSolidNodeList("nodes", eos, 
                              hmin = hmin,
                              hmax = hmax,
                              nPerh = nPerh,
                              kernelExtent = kernelExtent)
else:
    nodes = makeFluidNodeList("nodes", eos, 
                              hmin = hmin,
                              hmax = hmax,
                              nPerh = nPerh,
                              kernelExtent = kernelExtent)
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

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
    #mi = trapezoidalIntegration(answer.rhoInitial, ri - 0.5*dr, ri + 0.5*dr, 200)
    #rho[i] = mi/dr
    #eps[i] = trapezoidalIntegration(answer.Pinitial, ri - 0.5*dr, ri + 0.5*dr, 200)/((answer.gamma - 1.0)*mi)
    rho[i] = answer.rho(0.0, ri)
    mass[i] = rho[i]*dr
    eps[i] = answer.P(0.0, ri)/((answer.gamma - 1.0)*rho[i])

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   correctionOrder = correctionOrder,
                   volumeType = volumeType,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 filter = filter,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 HopkinsConductivity = HopkinsConductivity,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                filter = filter,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH)
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.evolveTotalEnergy")
output("hydro.HEvolution")
output("hydro.densityUpdate")
packages = [hydro]

#-------------------------------------------------------------------------------
# Optionally construct the reducing viscosity physics object.
#-------------------------------------------------------------------------------
if cullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(hydro.Q, WT, alphMax, alphMin, betaC, betaD, betaE, fKern, hopkinsCullenCorrection)
    packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.domainDecompositionIndependent = domainIndependent
output("integrator")
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
rbc1 = KidderIsentropicCapsuleEnforcementBoundary1d(integrator = integrator,
                                                    answer = answer,
                                                    nodeList = nodes,
                                                    nodeIDs = outerNodes,
                                                    interiorNodeIDs = interiorNodes,
                                                    hinitial = h0,
                                                    hmin = hmin,
                                                    hmax = hmax,
                                                    dr0 = dr)

# rbc1 = KidderIsentropicCapsuleBoundary1d(innerBoundary = False,
#                                          integrator = integrator,
#                                          answer = answer,
#                                          nodeList = nodes,
#                                          nrGhostNodes = nrGhost,
#                                          dr0 = dr)

# hydro.appendBoundary(rbc0)
# hydro.appendBoundary(rbc1)

packages += [rbc0, rbc1]

#-------------------------------------------------------------------------------
# Add our packages to the integrator.
#-------------------------------------------------------------------------------
for package in packages:
    integrator.appendPhysicsPackage(package)

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeDerivatives = True)
output("control")

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
    control.advance(goalTime, maxSteps)

print "Final energy (measured, analytic, error): ", (control.conserve.EHistory[-1],
                                                     answer.totalEnergy(goalTime),
                                                     (control.conserve.EHistory[-1] - answer.totalEnergy(0.0))/
                                                     answer.totalEnergy(goalTime))

#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
from SpheralGnuPlotUtilities import *

# Simulation results.
r = mpi.allreduce([x.x for x in nodes.positions().internalValues()], mpi.SUM)
rho = mpi.allreduce(nodes.massDensity().internalValues(), mpi.SUM)
Pfl = hydro.pressure
P = mpi.allreduce(Pfl[0].internalValues(), mpi.SUM)
v = mpi.allreduce([v.x for v in nodes.velocity().internalValues()], mpi.SUM)
eps = mpi.allreduce(nodes.specificThermalEnergy().internalValues(), mpi.SUM)
S = [p/d**answer.gamma for p, d in zip(P, rho)]

# Analytic results.
t = control.time()
rhoAns = [answer.rho(t, ri) for ri in r]
Pans = [answer.P(t, ri) for ri in r]
vAns = [answer.vr(t, ri) for ri in r]
epsAns = [answer.eps(t, ri) for ri in r]
Sans = [answer.S for ri in r]

# The ratio of the entropy to the expected value.
alpha = [ss/sa for ss, sa in zip(S, Sans)]
S0 = answer.S
print "Entropy L1, Linf in fractional error: ", (sum([abs(alphai - 1.0) for alphai in alpha])/len(alpha),
                                                     max(abs(max(alpha) - 1.0), abs(min(alpha) - 1.0)))

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

    rhoPlot = generateNewGnuPlot()
    rhoPlot.plot(rhoData)
    rhoPlot.replot(rhoAnsData)

    Pplot = generateNewGnuPlot()
    Pplot.plot(PData)
    Pplot.replot(PansData)

    Splot = generateNewGnuPlot()
    Splot.plot(SData)
    Splot.replot(SansData)

    vrPlot = generateNewGnuPlot()
    vrPlot.plot(vrData)
    vrPlot.replot(vrAnsData)

    epsPlot = generateNewGnuPlot()
    epsPlot.plot(epsData)
    epsPlot.replot(epsAnsData)

    alphaPlot = generateNewGnuPlot()
    alphaPlot.plot(alphaData)
    alphaPlot.replot(alphaAnsData)

    DepsDtfl = hydro.DspecificThermalEnergyDt
    DepsDt = mpi.allreduce(DepsDtfl[0].internalValues(), mpi.SUM)
    DepsData = Gnuplot.Data(r, DepsDt, title="DepsDt", inline=True)
    DepsPlot = generateNewGnuPlot()
    DepsPlot.plot(DepsData)

    DxDtfl = hydro.DxDt
    DxDt = mpi.allreduce([x.x for x in DxDtfl[0].internalValues()], mpi.SUM)
    DxData = Gnuplot.Data(r, DxDt, title="DxDt", inline=True)
    DxPlot = generateNewGnuPlot()
    DxPlot.plot(DxData)

    DvDtfl = hydro.DvDt
    DvDt = mpi.allreduce([x.x for x in DvDtfl[0].internalValues()], mpi.SUM)
    DvData = Gnuplot.Data(r, DvDt, title="DvDt", inline=True)
    DvPlot = generateNewGnuPlot()
    DvPlot.plot(DvData)

    DvDxfl = hydro.DvDx
    DvDx = mpi.allreduce([x.xx for x in DvDxfl[0].internalValues()], mpi.SUM)
    DvData = Gnuplot.Data(r, DvDx, title="DvDx", inline=True)
    DvPlot = generateNewGnuPlot()
    DvPlot.plot(DvData)

    DHDtfl = hydro.DHDt
    DHDt = mpi.allreduce([x.xx for x in DHDtfl[0].internalValues()], mpi.SUM)
    DHData = Gnuplot.Data(r, DHDt, title="DHDt", inline=True)
    DHPlot = generateNewGnuPlot()
    DHPlot.plot(DHData)

    DrhoDtfl = hydro.DmassDensityDt
    DrhoDt = mpi.allreduce(DrhoDtfl[0].internalValues(), mpi.SUM)
    DrhoData = Gnuplot.Data(r, DrhoDt, title="DrhoDt", inline=True)
    DrhoPlot = generateNewGnuPlot()
    DrhoPlot.plot(DrhoData)

    # If requested, output the profiles to an ASCII file.
    if profileASCII:
        f = open(os.path.join(dataDir, "Kidder_planar_profiles.txt"), "w")
        f.write(("#" + 11*"%20s" + "\n") % ('"radius"',
                                            '"rad velocity (sim)"',
                                            '"mass density (sim)"',
                                            '"pressure (sim)"',
                                            '"eps (sim)"',
                                            '"entropy (sim)"',
                                            '"rad velocity (ans)"',
                                            '"mass density (ans)"',
                                            '"pressure (ans)"',
                                            '"eps (ans)"',
                                            '"entropy (ans)"'))

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
