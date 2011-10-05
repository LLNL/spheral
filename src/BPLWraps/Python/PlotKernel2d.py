import Gnuplot
import SpheralTestUtilities
import SpheralGnuPlotUtilities
from math import *
from Spheral import *
import numpy

################################################################################
# Helper functions.
################################################################################
def findExtreme(method, container):
    result = eval(elementPattern % "eval(containerPattern % 'container')[0]")
    for element in eval(containerPattern % 'container'):
        result = method(result, eval(elementPattern % 'element'))
    return result

def SymTensorInverse(H):
    from Spheral import SymGeomTensor2d
    thpt = 1.0/(H.xx*H.yy - H.yx**2)
    return SymGeomTensor2d(thpt*H.yy, -thpt*H.yx,
                           -thpt*H.yx, thpt*H.xx)

################################################################################
# Plot a 2-D H tensor field.
################################################################################
def plotKernel2d(nodes,
                 h = 0.2,
                 np = 20,
                 xmin = 'auto',
                 xmax = 'auto',
                 ymin = 'auto',
                 ymax = 'auto',
                 plotGhost = 0,
                 colorDomains = 1,
                 plot = Gnuplot.Gnuplot(),
                 lineStyle = 'linetype -1 linewidth 1'
                 ):

    # Test if we're parallel or not.
    from SpheralGnuPlotUtilities import loadmpi
    mpi, rank, numDomains = loadmpi()

    # Do we show ghost nodes or not?
    if plotGhost:
        n = nodes.numNodes
        rvals = nodes.positions().allValues()
    else:
        n = nodes.numInternalNodes
        rvals = nodes.positions().internalValues()

    # Get the limits for the plot.
    if xmin == 'auto': xmin = min([x.x for x in rvals])
    if xmax == 'auto': xmax = max([x.x for x in rvals])
    if ymin == 'auto': ymin = min([x.y for x in rvals])
    if ymax == 'auto': ymax = max([x.y for x in rvals])
    if mpi:
        xmin = mpi.allreduce(xmin, mpi.MIN)
        xmax = mpi.allreduce(xmax, mpi.MAX)
        ymin = mpi.allreduce(ymin, mpi.MIN)
        ymax = mpi.allreduce(ymax, mpi.MAX)
    
    # Construct the h circle we will deform into the h contours.
    circ = []
    dtheta = 2*pi/(np - 1)
    for i in xrange(np):
        circ.append(GeomVector2d(h*cos(i*dtheta), h*sin(i*dtheta)))

    # Find the set of nodes we will be plotting.
    localr = []
    localH = []
    for i in xrange(n):
        if (nodes.positions()[i].x >= xmin and nodes.positions()[i].x <= xmax and
            nodes.positions()[i].y >= ymin and nodes.positions()[i].y <= ymax):
            localr.append(nodes.positions()[i])
            localH.append(nodes.Hfield()[i])
    assert len(localr) == len(localH)
    localN = [len(localr)]
##    print 'Plotting %i nodes in the region (%g, %g), (%g, %g)' % (localN[0],
##                                                                  xmin, ymin,
##                                                                  xmax, ymax)

    # Reduce the plotting info to processor 0
    if mpi:
        globalN = mpi.gather(localN, len(localN))
        globalr = mpi.gather(localr, len(localr))
        globalH = mpi.gather(localH, len(localH))
    else:
        globalN = localN
        globalr = localr
        globalH = localH

    # Only the master (root) process actually does the plotting.
    if rank == 0:

        # Set up our basic plot properties.
        plot('set size square')
        plot('set linestyle 1 ' + lineStyle)
        plot('set xrange [%g:%g]' % (xmin, xmax))
        plot('set yrange [%g:%g]' % (ymin, ymax))

        # Loop over the info from each domain.
        cumulativeN = 0
        for domain in xrange(numDomains):

            if colorDomains:
                plot('set linestyle %i linetype %i' % (domain + 1, domain + 1))

            # Loop over the nodes we are plotting for this domain.
            print 'Plotting ', globalN[domain], ' nodes from domain ', domain
            for i in xrange(cumulativeN, cumulativeN + globalN[domain]):
                assert i < len(globalr)
                assert i < len(globalH)
                r = globalr[i]
                H = globalH[i]

                # Get the inverse H transformation.
                Hinv = H.Inverse()

                # Create the h isocontour
                hcontour = circ[:]
                for element in hcontour:
                    element = Hinv*element

                # Shift the h isocontour to the correct position.
                x = numpy.array([0.0]*np)
                y = numpy.array([0.0]*np)
                for i in xrange(np):
                    x[i] = r.x + Hinv.xx*circ[i].x + Hinv.yx*circ[i].y
                    y[i] = r.y + Hinv.yx*circ[i].x + Hinv.yy*circ[i].y

                # Create the Gnuplot data for this contour.
                data = Gnuplot.Data(x, y,
                                    with_ = 'lines ls %i' % (domain + 1),
                                    inline = True)
                plot.replot(data)

            cumulativeN += globalN[domain]

        # That's it.
        return plot

    else:
        return SpheralGnuPlotUtilities.fakeGnuplot()
