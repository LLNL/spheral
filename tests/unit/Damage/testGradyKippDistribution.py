#-------------------------------------------------------------------------------
# Seed a rod of stainless steel with the Grady-Kipp distribution of flaws.  This
# test simply checks that we are statistically seeding the expected distribution.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from Spheral2d import *
from SpheralTestUtilities import *
from math import *
import mpi, bisect
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(seed = 'lattice',
            xlength = 1.0,
            ylength = 1.0,
            nx0 = 24,
            ny0 = 24,
            nPerh = 2.01,
            rho0 = 1.0,
            kWeibull = 8.8e4,
            mWeibull = 2.63,
            volume = 1.0,
            randomSeed = 109482993,
            ntests = 5)

xmin1 = (-0.5*xlength, -0.5*ylength)
xmax1 = ( 0.0,          0.5*ylength)
xmin2 = ( 0.0,         -0.5*ylength)
xmax2 = ( 0.5*xlength,  0.5*ylength)

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("Weibull flaw distribution test")

#-------------------------------------------------------------------------------
# Make up material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(1.0, 2.0)

#-------------------------------------------------------------------------------
# Interpolation kernel
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList("test nodes", eos, NullStrength())
db = DataBase()
db.appendNodeList(nodes)

#-------------------------------------------------------------------------------
# Functions to help accumulating the distribution functions
#-------------------------------------------------------------------------------
energiesBA, fBA = [], []
energiesO, fO = [], []
energiesP, weightP, fP = [], [], []

def accumulateP(nodes, damageP):
    energiesP.append([])
    weightP.append([])
    fP.append([])
    flaws = energiesP[-1]
    weight = weightP[-1]
    fcum = fP[-1]
    minNumFlaws = damageP.minFlawsPerNode
    for i in xrange(nodes.numInternalNodes):
        epsmini, epsmaxi, Vi, nFlawsi = damageP.minFlaw[i], damageP.maxFlaw[i], damageP.initialVolume[i], damageP.numFlaws[i]
        nsteps = min(10, nFlawsi)
        deps = (epsmaxi - epsmini)/(nsteps - 1)
        for j in xrange(nsteps):
            epsj = epsmini + j*deps
            weight.append((epsj, Vi*kWeibull*(epsj**mWeibull - epsmini**mWeibull) + 1.0))
    weight.sort()
    wsum = 0.0
    for epsi, wi in weight:
        wsum += wi
        flaws.append(epsi)
        fcum.append(wsum)
    assert len(fcum) == len(flaws)

#-------------------------------------------------------------------------------
# Iterate over the number of tests, and check the seeeded distribution of flaws.
#-------------------------------------------------------------------------------
from GenerateNodeDistribution2d import *
from CompositeNodeDistribution import *
if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

for test in xrange(ntests):
    nx = (test + 1)*nx0
    ny = (test + 1)*ny0
    nodes.numInternalNodes = 0

    print "Generating node distribution."
    generator = GenerateNodeDistribution2d(nx,
                                           ny,
                                           rho0,
                                           seed,
                                           xmin = xmin1,
                                           xmax = xmax2,
                                           nNodePerh = nPerh)
    distributeNodes2d((nodes, generator))
    output('mpi.reduce(nodes.numInternalNodes, mpi.MIN)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.MAX)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.SUM)')

    m = nodes.mass()
    rho = nodes.massDensity()
    mask1 = IntField("mask1", nodes, 1)
    print "Area sum: ", mpi.allreduce(sum([m[i]/rho[i] for i in xrange(nodes.numInternalNodes)]), mpi.SUM)

    # Construct the flaws.
    localFlawsBA = weibullFlawDistributionBenzAsphaug(volume,
                                                      1.0,
                                                      randomSeed,
                                                      kWeibull,
                                                      mWeibull,
                                                      nodes,
                                                      1,
                                                      1,
                                                      mask1)
    localFlawsO = weibullFlawDistributionOwen(randomSeed,
                                              kWeibull,
                                              mWeibull,
                                              nodes,
                                              1,     # numFlawsPerNode
                                              1.0,    # volumeMultiplier
                                              mask1)

    # Collect the distribution function of flaws.
    for (localFlaws, energies, f) in [(localFlawsBA, energiesBA, fBA),
                                      (localFlawsO, energiesO, fO)]:
        flaws = []
        for i in xrange(nodes.numInternalNodes):
            flaws.extend(localFlaws[i])
        flaws.sort()
        energies.append(flaws)
        f.append(range(len(flaws)))

    # Also the ProbabilisticDamageModel
    pdamage = ProbabilisticDamageModel(nodeList = nodes,
                                       kernel = WT,
                                       kWeibull = kWeibull,
                                       mWeibull = mWeibull,
                                       seed = randomSeed,
                                       minFlawsPerNode = 10)
    pdamage.initializeProblemStartup(db)
    accumulateP(nodes, pdamage)

    #--------------------------------------------------------------------------
    # Cover the same volume with a non-uniform node distribution and check
    # that the Benz-Asphaug-Owen flaw algorithm works in this case as well.
    # We also simultaneously test the implicit distribution from the
    # ProbabilisticDamage model.
    #--------------------------------------------------------------------------
    nx1 = (test + 1)*nx0/2
    ny1 = (test + 1)*ny0
    nx2 = 2*nx1
    ny2 = 2*ny1
    nodes.numInternalNodes = 0

    print "Generating node distribution."
    generator1 = GenerateNodeDistribution2d(nx1,
                                            ny1,
                                            rho0,
                                            seed,
                                            xmin = xmin1,
                                            xmax = xmax1,
                                            nNodePerh = nPerh)
    generator2 = GenerateNodeDistribution2d(nx2,
                                            ny2,
                                            rho0,
                                            seed,
                                            xmin = xmin2,
                                            xmax = xmax2,
                                            nNodePerh = nPerh)
    generator = CompositeNodeDistribution(generator1, generator2)
    distributeNodes2d((nodes, generator))
    output('mpi.reduce(nodes.numInternalNodes, mpi.MIN)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.MAX)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.SUM)')

    m = nodes.mass()
    rho = nodes.massDensity()
    mask2 = IntField("mask2", nodes, 1)
    print "Area sum: ", mpi.allreduce(sum([m[i]/rho[i] for i in xrange(nodes.numInternalNodes)]), mpi.SUM)

    # Collect the distribution function of flaws (BAO)
    localFlawsO = weibullFlawDistributionOwen(randomSeed,
                                              kWeibull,
                                              mWeibull,
                                              nodes,
                                              1,      # minFlawsPerNode
                                              1.0,    # volumeMultiplier
                                              mask2)
    for (localFlaws, energies, f) in [(localFlawsO, energiesO, fO)]:
        flaws = []
        for i in xrange(nodes.numInternalNodes):
            flaws.extend(localFlaws[i])
        flaws.sort()
        energies.append(flaws)
        f.append(range(len(flaws)))

    # Also the ProbabilisticDamageModel
    pdamage = ProbabilisticDamageModel(nodeList = nodes,
                                       kernel = WT,
                                       kWeibull = kWeibull,
                                       mWeibull = mWeibull,
                                       seed = randomSeed,
                                       minFlawsPerNode = 10)
    pdamage.initializeProblemStartup(db)
    accumulateP(nodes, pdamage)

assert len(fBA) == ntests
assert len(energiesBA) == ntests
assert len(fO) == 2*ntests
assert len(energiesO) == 2*ntests
assert len(fP) == 2*ntests
assert len(energiesP) == 2*ntests

#-------------------------------------------------------------------------------
# Now plot the results.
#-------------------------------------------------------------------------------
fig1 = plt.figure(1)
fig2 = plt.figure(2)
fig3 = plt.figure(3)

#...............................................................................
plt.figure(1)
plt.xscale("log")
plt.yscale("log")
for i in xrange(ntests):
    plt.plot(energiesBA[i], fBA[i], "-", linewidth=2*(ntests - i) + 3, label=("BA: N = %i" % ((i + 1)**2*nx0*ny0)))

# Plot the expectation.
emin = min([min(x) for x in (energiesBA)])# + energiesO)])
emax = max([max(x) for x in (energiesBA)])# + energiesO)])
eans = [emin + 0.01*(emax - emin)*i for i in range(101)]
fans = [kWeibull * e**mWeibull for e in eans]
plt.plot(eans, fans, "-", linewidth=3, label="Analytic")

plt.legend(loc="best")
plt.xlabel(r"Flaw Activation Energy ($\varepsilon$)")
plt.ylabel(r"$n(\varepsilon^f < \varepsilon)$")

#...............................................................................
plt.figure(2)
plt.xscale("log")
plt.yscale("log")
for i in xrange(2*ntests):
    if i % 2 == 1:
        nx1 = (i/2 + 1)*nx0/2
        ny1 = (i/2 + 1)*ny0
        nx2 = 2*nx1
        ny2 = 2*ny1
        TT = "O: (N1,N2)=(%i,%i)" % (nx1*ny1, nx2*ny2)
    else:
        nx = (i/2 + 1)*nx0
        ny = (i/2 + 1)*ny0
        TT = "O: N = %i" % (nx*ny)
    plt.plot(energiesO[i], fO[i], "-", linewidth=2, label=TT)

# Plot the expectation.
emin = min([min(x) for x in energiesO])#(energiesBA + energiesO)])
emax = max([max(x) for x in energiesO])#(energiesBA + energiesO)])
eans = [emin + 0.01*(emax - emin)*i for i in range(101)]
fans = [kWeibull * e**mWeibull for e in eans]
plt.plot(eans, fans, "-", linewidth=3, label="Analytic")

plt.legend(loc="lower right")
plt.xlabel(r"Flaw Activation Energy ($\varepsilon$)")
plt.ylabel(r"$n(\varepsilon^f < \varepsilon)$")

#...............................................................................
plt.figure(3)
plt.xscale("log")
plt.yscale("log")
for i in xrange(2*ntests):
    if i % 2 == 1:
        nx1 = (i/2 + 1)*nx0/2
        ny1 = (i/2 + 1)*ny0
        nx2 = 2*nx1
        ny2 = 2*ny1
        TT = "P: (N1,N2)=(%i,%i)" % (nx1*ny1, nx2*ny2)
    else:
        nx = (i/2 + 1)*nx0
        ny = (i/2 + 1)*ny0
        TT = "P: N = %i" % (nx*ny)
    plt.plot(energiesP[i], fP[i], "-", linewidth=2, label=TT)

# Plot the expectation.
emin = min([min(x) for x in energiesP])
emax = max([max(x) for x in energiesP])
eans = [emin + 0.01*(emax - emin)*i for i in range(101)]
fans = [kWeibull * e**mWeibull for e in eans]
plt.plot(eans, fans, "-", linewidth=3, label="Analytic")

plt.legend(loc="lower right")
plt.xlabel(r"Flaw Activation Energy ($\varepsilon$)")
plt.ylabel(r"$n(\varepsilon^f < \varepsilon)$")
