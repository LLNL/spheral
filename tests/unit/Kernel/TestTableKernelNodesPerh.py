import numpy as np
import sys

from Spheral import *
from SpheralTestUtilities import *
from SpheralMatplotlib import *

#-------------------------------------------------------------------------------
# What kernels should we plot
#-------------------------------------------------------------------------------
kernels = sys.argv[1:]
print(kernels)


#-------------------------------------------------------------------------------
# Define some dimensional functions for summing expected kernel values
#-------------------------------------------------------------------------------
def sumKernelValues1d(WT, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = sum([abs(WT.gradValue(abs(etax), 1.0)) for etax in np.arange(-etamax, etamax, deta)])
    return result

def sumKernelValues2d(WT, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = 0.0
    for etay in np.arange(-etamax, etamax, deta):
        for etax in np.arange(-etamax, etamax, deta):
            eta = sqrt(etax*etax + etay*etay)
            result += abs(WT.gradValue(eta, 1.0))
    return sqrt(result)

def sumKernelValues3d(WT, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = 0.0
    for etaz in np.arange(-etamax, etamax, deta):
        for etay in np.arange(-etamax, etamax, deta):
            for etax in np.arange(-etamax, etamax, deta):
                eta = sqrt(etax*etax + etay*etay + etaz*etaz)
                result += abs(WT.gradValue(eta, 1.0))
    return (result)**(1.0/3.0)

kernelDict = {'spline': [BSplineKernel1d(),
                         BSplineKernel2d(),
                         BSplineKernel3d()],
              }

titleDict = {'spline'     : 'B Spline Kernel',
             'h'          : 'H kernel',
             'h10'        : 'H kernel (extent = 10)',
             'quartic'    : 'Quartic Spline Kernel',
             'w4spline'   : 'W4 Spline Kernel',
             'gauss'      : 'Gaussian Kernel',
             'supergauss' : 'SuperGaussian Kernel',
             'pigauss'    : 'Pi Gaussian Kernel',
             'sinc'       : 'Sinc Kernel',
             'poly1'      : 'Linear Polynomial Sinc approx Kernel',
             'poly3'      : 'Cubic Polynomial Sinc approx Kernel',
             'poly5'      : 'Quintic Polynomial Sinc approx Kernel',
             'poly7'      : 'Septic Polynomial Sinc approx Kernel',
             'spline3'    : '3rd order b spline Kernel',
             'spline5'    : '5th order b spline Kernel',
             'spline7'    : '7th order b spline Kernel',
             'spline9'    : '9th order b spline Kernel',
             'spline11'   : '11th order b spline Kernel',
             'WendlandC2' : 'Wendland C2',
             'WendlandC4' : 'Wendland C4',
             }

for Wstr in kernels:
    title(Wstr)

    nDim = 0
    if str(Wstr).split()[0][-2:] == "1d":
        nDim = 1
    elif str(Wstr).split()[0][-2:] == "2d":
        nDim = 2
    elif str(Wstr).split()[0][-2:] == "3d":
        nDim = 3
    assert nDim > 0

    # Plot the kernel basics
    WT = eval(f"TableKernel{nDim}d({Wstr}())")
    #plotTableKernel(WT)

    # Now how well do we recover nPerh based on kernel sums?
    etamax = WT.kernelExtent
    nperh0 = np.arange(0.5, 20.0, 0.1)
    nperh1 = []
    for nperh in nperh0:
        Wsum = eval(f"sumKernelValues{nDim}d(WT, {nperh})")
        nperh1.append(WT.equivalentNodesPerSmoothingScale(Wsum))
    nperh1 = np.array(nperh1)

    plot = newFigure()
    plot.plot(nperh0, nperh1, "b*-")
    plot.set_title(f"{Wstr} n per h lookup test")
    plot.set_xlabel("nperh actual")
    plot.set_ylabel("nperh estimated")

    err = (nperh1 - nperh0)/nperh0
    plot = newFigure()
    plot.plot(nperh0, err, "r*-")
    plot.set_title(f"{Wstr} n per h lookup test error")
    plot.set_xlabel("nperh actual")
    plot.set_ylabel("Error")
