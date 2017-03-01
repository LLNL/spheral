from Spheral import *
from SpheralTestUtilities import *
import Gnuplot
import numpy
from SpheralGnuPlotUtilities import *

################################################################################
def plotW(plot, W, xmin=0.0, xmax=2.0, numPnts=200, Hdet=1.0, title='',
          lineTitle=''):
    dx = (xmax - xmin)/(numPnts - 1)
    x = numpy.array(range(numPnts))
    y = numpy.array([0.0]*numPnts)
    i1,i2,i3 = 0,0,0
    x = dx*x + xmin
    for i in xrange(numPnts):
        y[i] = W(x[i], Hdet)
        if (i>0):
            i1 += abs(y[i]+y[i-1])/2.0 * dx
            i2 += abs(y[i]+y[i-1])/2.0 * dx * x[i]
            i3 += abs(y[i]+y[i-1])/2.0 * dx * x[i]*x[i]
    Pi = 3.14159
    print "{0:3.3f} {1:3.3f} {2:3.3f}".format(i1*2.0,i2*2.0,i3)
    plot('set xrange [0:1]')
    plot.xlabel('eta')
    plot.ylabel('W')
    if title:
        plot.title(title)
    ymax = 0
    for i in xrange(numPnts):
        if abs(y[i]) > ymax:
            ymax = abs(y[i])
    data = Gnuplot.Data(x/xmax, y/ymax, with_='lines', title=lineTitle)

    plot.replot(data)
    plot('set xrange[0:1]')
    return

import sys, string
kernels = map(string.lower, sys.argv[1:])
print kernels

numPts = 51
dx = 1.0/(numPts - 1)

################################################################################
numPoints = 100

kernelDict = {'spline': [BSplineKernel2d()],
              'w4spline': [W4SplineKernel2d()],
              'wendlandc4': [WendlandC4Kernel2d()],
              'wendlandc6': [WendlandC6Kernel2d()],
              'expinv': [ExpInvKernel1d(),
                         ExpInvKernel2d(),
                         ExpInvKernel3d()],
              'quartic': [QuarticSplineKernel3d()],
##              'gauss': [GaussianKernel1d(3),
##                        GaussianKernel2d(3),
##                        GaussianKernel3d(3)],
##              'supergauss': [SuperGaussianKernel1d(),
##                             SuperGaussianKernel2d(),
##                             SuperGaussianKernel3d()],
##              'pigauss': [PiGaussianKernel1d(1.0),
##                          PiGaussianKernel2d(1.0),
##                          PiGaussianKernel3d(1.0)],
##              'sinc': [SincKernel1d(2),
##                       SincKernel2d(2),
##                       SincKernel3d(2)],
##              'poly1': [NSincPolynomialKernel1d(1),
##                        NSincPolynomialKernel2d(1),
##                        NSincPolynomialKernel3d(1)],
##              'poly3': [NSincPolynomialKernel1d(3),
##                        NSincPolynomialKernel2d(3)],
##              'poly5': [NSincPolynomialKernel1d(5),
##                        NSincPolynomialKernel2d(5)],
##              'poly7': [NSincPolynomialKernel1d(7),
##                        NSincPolynomialKernel2d(7)],
              'spline3': [NBSplineKernel2d(3)],
              'spline5': [NBSplineKernel2d(5)],
              'spline7': [NBSplineKernel2d(7)],
              'spline9': [NBSplineKernel2d(9)],
##              'spline11': [NBSplineKernel1d(11),
##                           NBSplineKernel2d(11),
##                           NBSplineKernel3d(11)],
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
             'wendlandc4': 'Wendland C4 Kernel',
             'wendlandc6': 'Wendland C6 Kernel',
             'expinv' : 'Exponential inverse Kernel',
             }

plots = []
for i in xrange(2):
    plots.append(generateNewGnuPlot())
for kernel in kernels:
    title(titleDict[kernel])

    for W in kernelDict[kernel]:
        # Build a tabular version of the kernel
        WT = eval('TableKernel' + str(W).split()[0][-2:] + '(W, numPoints)')
        output("WT")
        output("WT.volumeNormalization")
        output("WT.kernelExtent")
        plotW(plots[-1], WT.kernelValue, 0.0, W.kernelExtent,
              title='Kernels',
              lineTitle = titleDict[kernel])
        plotW(plots[-2], WT.gradValue, 0.0, W.kernelExtent,
              title = 'Kernel Gradients',
              lineTitle = titleDict[kernel])

