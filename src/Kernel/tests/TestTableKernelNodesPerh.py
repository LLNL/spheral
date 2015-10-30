import Gnuplot
import numpy

from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

################################################################################
def plotW(plot, W, xmin=0.0, xmax=2.0, numPnts=200, Hdet=1.0, title='',
          lineTitle=''):
    dx = (xmax - xmin)/(numPnts - 1)
    x = numpy.array(range(numPnts))
    y = numpy.array([0.0]*numPnts)
    x = dx*x + xmin
    for i in xrange(numPnts):
        y[i] = W(x[i], Hdet)
    plot('set xrange [%f:%f]' % (xmin, xmax))
    plot.xlabel('r')
    plot.ylabel('W(r)')
    if title:
        plot.title(title)
    data = Gnuplot.Data(x, y, with_='lines', title=lineTitle)
    plot.replot(data)
    return

import sys, string
kernels = map(string.lower, sys.argv[1:])
print kernels

numPts = 51
dx = 1.0/(numPts - 1)

################################################################################
numPoints = 100

kernelDict = {'spline': [BSplineKernel1d(),
                         BSplineKernel2d(),
                         BSplineKernel3d()],
              }

titleDict = {'spline':  'B Spline Kernel',
             'h': 'H kernel',
             'h10': 'H kernel (extent = 10)',
             'quartic': 'Quartic Spline Kernel',
             'w4spline': 'W4 Spline Kernel',
             'gauss': 'Gaussian Kernel',
             'supergauss': 'SuperGaussian Kernel',
             'pigauss': 'Pi Gaussian Kernel',
             'sinc': 'Sinc Kernel',
             'poly1': 'Linear Polynomial Sinc approx Kernel',
             'poly3': 'Cubic Polynomial Sinc approx Kernel',
             'poly5': 'Quintic Polynomial Sinc approx Kernel',
             'poly7': 'Septic Polynomial Sinc approx Kernel',
             'spline3': '3rd order b spline Kernel',
             'spline5': '5th order b spline Kernel',
             'spline7': '7th order b spline Kernel',
             'spline9': '9th order b spline Kernel',
             'spline11': '11th order b spline Kernel',
             }

data = []
plots = []
plotWsum = generateNewGnuPlot()
plotWsum("set xlabel 'Nodes per smoothing scale'")
plotWsum("set ylabel 'W_{sum}'")
#plotWsum("set logscale y")
for kernel in kernels:
    title(titleDict[kernel])
    for W in kernelDict[kernel]:

        nDim = 0
        if str(W).split()[0][-2:] == "1d":
            nDim = 1
        elif str(W).split()[0][-2:] == "2d":
            nDim = 2
        elif str(W).split()[0][-2:] == "3d":
            nDim = 3
        assert nDim > 0
        
        # Build the TableKernel.
        WT = eval('TableKernel' + str(W).split()[0][-2:] + '(W, numPoints)')
        #WH = eval('HKernel' + str(W).split()[0][-2:] + '(W.kernelExtent)')

        # Go over the range of nodes per H, and see how well the TableKernel predicts
        Wsumarray = []
        actualnperh = []
        lookupnperh = []
        for nperh in [0.5*float(x) for x in range(1, 20)]:
            deta = 1.0/nperh
            npoints = int(WT.kernelExtent*nperh)
            Wsum = 0.0

            eta = 0.0
            for i in xrange(-npoints, npoints):
                eta = abs(i*deta)
                if nDim == 1 and eta > 1.0e-5:
                    Wsum += abs(WT.gradValue(eta, 1.0))
                if nDim > 1:
                    for j in xrange(-npoints, npoints):
                        eta = sqrt((i*deta)**2 +
                                   (j*deta)**2)
                        if nDim == 2 and eta > 1.0e-5:
                            Wsum += abs(WT.gradValue(eta, 1.0))
                        if nDim > 2:
                            for k in xrange(-npoints, npoints):
                                eta = sqrt((i*deta)**2 +
                                           (j*deta)**2 +
                                           (k*deta)**2)
                                if eta > 1.0e-5:
                                    Wsum += abs(WT.gradValue(eta, 1.0))

            Wsum = Wsum**(1.0/nDim)
            result = WT.equivalentNodesPerSmoothingScale(Wsum)
            Wsumarray.append(Wsum)
            actualnperh.append(nperh)
            lookupnperh.append(result)

        # Plot the lookup results.
        actualdata = Gnuplot.Data(Wsumarray, actualnperh,
                                  with_ = "lines",
                                  title = "Actual n per h",
                                  inline = True)
        lookupdata = Gnuplot.Data(Wsumarray, lookupnperh,
                                  with_ = "points",
                                  title = "Lookup n per h",
                                  inline = True)
        nperhdata = Gnuplot.Data(actualnperh, lookupnperh,
                                 with_="points",
                                 title = None,
                                 inline = True)
        data.extend([actualdata, lookupdata, nperhdata])

        plot = generateNewGnuPlot()
        plot.plot(actualdata)
        plot.replot(lookupdata)
        plot.title("%-d" % nDim)
        plot.xlabel("W_{sum}")
        plot.ylabel("n per h")
        plot.refresh()
        plots.append(plot)

        p = generateNewGnuPlot()
        p.plot(nperhdata)
        p.title("Comparison of actual vs. lookup nperh")
        p.xlabel("actual nperh")
        p.ylabel("lookup nperh")
        p.refresh()
        plots.append(p)

        # Plot Wsum as a function of n per h.
        nperhdata = Gnuplot.Data(WT.nperhValues, WT.WsumValues,
                                 with_ = "lines",
                                 title = ("%i -D" % (nDim)),
                                 inline = True)
        plotWsum.replot(nperhdata)
