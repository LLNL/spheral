#-------------------------------------------------------------------------------
# Seed a rod of stainless steel with the Grady-Kipp distribution of flaws.  This
# test simply checks that we are statistically seeding the expected distribution.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from SolidSpheral import *
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

xmin = (-0.5*xlength, -0.5*ylength)
xmax = ( 0.5*xlength,  0.5*ylength)

neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
numGridLevels = 20
topGridCellSize = max(xlength, ylength)
origin = Vector2d(-xlength, -ylength)

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("Weibull flaw distribution test")

#-------------------------------------------------------------------------------
# Make up material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS2d(1.0, 2.0)

#-------------------------------------------------------------------------------
# Interpolation kernel
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Create a NodeLists.
#-------------------------------------------------------------------------------
nodes = SphNodeList2d("Thpt", eos, WT, WT)

#-------------------------------------------------------------------------------
# Construct the neighbor object
#-------------------------------------------------------------------------------
neighbor = NestedGridNeighbor2d(nodes,
                                neighborSearchType,
                                numGridLevels,
                                topGridCellSize,
                                origin,
                                kernelExtent)
nodes.registerNeighbor(neighbor)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
db.appendNodeList(nodes)
output('db')
output('db.numNodeLists')
output('db.numFluidNodeLists')

#-------------------------------------------------------------------------------
# Iterate over the number of tests, and check the seeeded distribution of flaws.
#-------------------------------------------------------------------------------
from GenerateNodeDistribution2d import *
from PeanoHilbertDistributeNodes import distributeNodes2d

energiesBA = []
fBA = []
energiesO = []
fO = []
for test in xrange(ntests):
    nx = (test + 1)*nx0
    ny = (test + 1)*ny0
    nodes.numInternalNodes = 0
    m0 = (xlength*ylength)*rho0/(nx*ny)
    dx = xlength/nx
    dy = ylength/ny
    hx = nPerh*dx
    hy = nPerh*dy
    H0 = SymTensor2d(1.0/hx, 0.0,
                     0.0, 1.0/hx)

    print "Generating node distribution."
    generator = GenerateNodeDistribution2d(nx,
                                           ny,
                                           rho0,
                                           seed,
                                           xmin = xmin,
                                           xmax = xmax,
                                           nNodePerh = nPerh)
    distributeNodes2d((nodes, generator))
    output('mpi.reduce(nodes.numInternalNodes, mpi.MIN)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.MAX)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.SUM)')

    # Set the node masses.
    nodes.mass(ScalarField2d("tmp", nodes, m0))

    # Set the smoothing scales.
    nodes.Hfield(SymTensorField2d("tmp", nodes, H0))

    # Construct the flaws.
    localFlawsBA = weibullFlawDistributionBenzAsphaug2d(volume,
                                                        1.0,
                                                        randomSeed,
                                                        kWeibull,
                                                        mWeibull,
                                                        nodes,
                                                        1,
                                                        1)
    localFlawsO = weibullFlawDistributionOwen2d(randomSeed,
                                                kWeibull,
                                                mWeibull,
                                                nodes,
                                                nx,     # numFlawsPerNode
                                                dx)     # volumeMultiplier

    # Collect the distribution function of flaws.
    for (localFlaws, energies, f) in [(localFlawsBA, energiesBA, fBA),
                                      (localFlawsO, energiesO, fO)]:
        flaws = []
        for i in xrange(nodes.numInternalNodes):
            flaws.extend(localFlaws[i])
        flaws.sort()
        energies.append(flaws)
        f.append(range(len(flaws)))

assert len(fBA) == ntests
assert len(energiesBA) == ntests
assert len(fO) == ntests
assert len(energiesO) == ntests

#-------------------------------------------------------------------------------
# Now plot the results.
#-------------------------------------------------------------------------------
import Gnuplot
cache = []
p = Gnuplot.Gnuplot()
p("set logscale xy")
for i in xrange(ntests):
    data = Gnuplot.Data(energiesBA[i], fBA[i],
                        with = "lines lw 2",
                        title = ("BA: N = %i" % ((i + 1)**2*nx0*ny0)),
                        inline=True)
    cache.append(data)
    p.replot(data)
for i in xrange(ntests):
    data = Gnuplot.Data(energiesO[i], fO[i],
                        with = "lines lw 2",
                        title = ("O: N = %i" % ((i + 1)**2*nx0*ny0)),
                        inline=True)
    cache.append(data)
    p.replot(data)

# Plot the expectation.
## emin = min([min(x) for x in (energiesBA + energiesO)])
## emax = max([max(x) for x in (energiesBA + energiesO)])
## eans = [emin + 0.01*(emax - emin)*i for i in range(101)]
## fans = [kWeibull * e**mWeibull for e in eans]
## answer = Gnuplot.Data(eans, fans,
##                       with = "lines lw 3",
##                       title = "Analytic",
##                       inline = True)
## p.replot(answer)
p("set key top left")
p.xlabel("Flaw Activation Energy ({/Symbol e})")
p.ylabel("n({/Symbol e})")
p.refresh()
