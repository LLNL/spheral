#-------------------------------------------------------------------------------
# Seed a rod of stainless steel with the Grady-Kipp distribution of flaws.  This
# test simply checks that we are statistically seeding the expected distribution.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from SolidSpheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from math import *
import mpi

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
nodes = makeFluidNodeList("test nodes", eos)

#-------------------------------------------------------------------------------
# Iterate over the number of tests, and check the seeeded distribution of flaws.
#-------------------------------------------------------------------------------
from GenerateNodeDistribution2d import *
from CompositeNodeDistribution import *
from PeanoHilbertDistributeNodes import distributeNodes2d

energiesBA = []
fBA = []
energiesO = []
fO = []
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
    print "Area sum: ", mpi.allreduce(sum([m[i]/rho[i] for i in xrange(nodes.numInternalNodes)]), mpi.SUM)

    # Construct the flaws.
    localFlawsBA = weibullFlawDistributionBenzAsphaug(volume,
                                                      1.0,
                                                      randomSeed,
                                                      kWeibull,
                                                      mWeibull,
                                                      nodes,
                                                      1,
                                                      1)
    localFlawsO = weibullFlawDistributionOwen(randomSeed,
                                              kWeibull,
                                              mWeibull,
                                              nodes,
                                              1,     # numFlawsPerNode
                                              1.0)    # volumeMultiplier

    # Collect the distribution function of flaws.
    for (localFlaws, energies, f) in [(localFlawsBA, energiesBA, fBA),
                                      (localFlawsO, energiesO, fO)]:
        flaws = []
        for i in xrange(nodes.numInternalNodes):
            flaws.extend(localFlaws[i])
        flaws.sort()
        energies.append(flaws)
        f.append(range(len(flaws)))

    #--------------------------------------------------------------------------
    # Cover the same volume with a non-uniform node distribution and check
    # that the Benz-Asphaug-Owen flaw algorithm works in this case as well.
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
    print "Area sum: ", mpi.allreduce(sum([m[i]/rho[i] for i in xrange(nodes.numInternalNodes)]), mpi.SUM)

    # Construct the flaws.
    localFlawsO = weibullFlawDistributionOwen(randomSeed,
                                              kWeibull,
                                              mWeibull,
                                              nodes,
                                              1,      # minFlawsPerNode
                                              1.0)    # volumeMultiplier

    # Collect the distribution function of flaws.
    for (localFlaws, energies, f) in [(localFlawsO, energiesO, fO)]:
        flaws = []
        for i in xrange(nodes.numInternalNodes):
            flaws.extend(localFlaws[i])
        flaws.sort()
        energies.append(flaws)
        f.append(range(len(flaws)))

assert len(fBA) == ntests
assert len(energiesBA) == ntests
assert len(fO) == 2*ntests
assert len(energiesO) == 2*ntests

#-------------------------------------------------------------------------------
# Now plot the results.
#-------------------------------------------------------------------------------
import Gnuplot
cache = []
p = Gnuplot.Gnuplot()
p("set logscale xy")
for i in xrange(ntests):
    data = Gnuplot.Data(energiesBA[i], fBA[i],
                        with_ = "lines lw 2",
                        title = ("BA: N = %i" % ((i + 1)**2*nx0*ny0)),
                        inline=True)
    cache.append(data)
    p.replot(data)
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
    data = Gnuplot.Data(energiesO[i], fO[i],
                        with_ = "lines lw 2",
                        title = TT,
                        inline=True)
    cache.append(data)
    p.replot(data)

# Plot the expectation.
emin = min([min(x) for x in (energiesBA + energiesO)])
emax = max([max(x) for x in (energiesBA + energiesO)])
eans = [emin + 0.01*(emax - emin)*i for i in range(101)]
fans = [kWeibull * e**mWeibull for e in eans]
answer = Gnuplot.Data(eans, fans,
                      with_ = "lines lw 3",
                      title = "Analytic",
                      inline = True)
p.replot(answer)
p("set key top left")
p.xlabel("Flaw Activation Energy ({/Symbol e})")
p.ylabel("n({/Symbol e})")
p.refresh()
