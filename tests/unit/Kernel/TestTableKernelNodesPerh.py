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
# SPH zeroth moment algorithm
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

#-------------------------------------------------------------------------------
# ASPH second moment algorithm
#-------------------------------------------------------------------------------
def sumKernelValuesASPH1d(WT, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = sum([abs(WT.gradValue(abs(etax), 1.0)*etax*etax) for etax in np.arange(-etamax, etamax, deta)])
    return result

def sumKernelValuesASPH2d(WT, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = SymTensor2d()
    for etay in np.arange(-etamax, etamax, deta):
        for etax in np.arange(-etamax, etamax, deta):
            eta = Vector2d(etax, etay)
            result += abs(WT.gradValue(eta.magnitude(), 1.0)) * eta.selfdyad()
    return sqrt(0.5*(result.eigenValues().sumElements()))

def sumKernelValuesASPH3d(WT, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = SymTensor3d()
    for etaz in np.arange(-etamax, etamax, deta):
        for etay in np.arange(-etamax, etamax, deta):
            for etax in np.arange(-etamax, etamax, deta):
                eta = Vector3d(etax, etay, etaz)
                result += abs(WT.gradValue(eta.magnitude(), 1.0)) * eta.selfdyad()
    return ((result.eigenValues().sumElements())/3.0)**(1.0/3.0)

def sumKernelValuesSlice2d(WT, nhat, detax, detay):
    etamax = WT.kernelExtent
    result = SymTensor2d()
    for etay in np.arange(-etamax, etamax, detay):
        for etax in np.arange(-etamax, etamax, detax):
            eta = Vector2d(etax, etay)
            result += abs(WT.gradValue(eta.magnitude(), 1.0)) * eta.selfdyad()
    return sqrt((result*nhat).magnitude())

def sumKernelValuesSlice3d(WT, nhat, detax, detay, detaz):
    etamax = WT.kernelExtent
    result = SymTensor3d()
    for etaz in np.arange(-etamax, etamax, detaz):
        for etay in np.arange(-etamax, etamax, detay):
            for etax in np.arange(-etamax, etamax, detax):
                eta = Vector3d(etax, etay, etaz)
                result += abs(WT.gradValue(eta.magnitude(), 1.0)) * eta.selfdyad()
    return ((result*nhat).magnitude())**(1.0/3.0)

#-------------------------------------------------------------------------------
# Here we go...
#-------------------------------------------------------------------------------
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
    nperh0 = np.arange(1.0/etamax, 10.0, 0.1)
    nperhSPH = []
    nperhASPH = []
    WsumSPH = []
    WsumASPH = []
    for nperh in nperh0:
        Wsumi = eval(f"sumKernelValues{nDim}d(WT, {nperh})")
        WsumASPHi = eval(f"sumKernelValuesASPH{nDim}d(WT, {nperh})")
        WsumSPH.append(Wsumi)
        WsumASPH.append(WsumASPHi)
        nperhSPH.append(WT.equivalentNodesPerSmoothingScale(Wsumi))
        nperhASPH.append(WT.equivalentNodesPerSmoothingScaleASPH(WsumASPHi))
    nperhSPH = np.array(nperhSPH)
    nperhASPH = np.array(nperhASPH)
    WsumSPH = np.array(WsumSPH)
    WsumASPH = np.array(WsumASPH)

    # SPH fit for nperh(Wsum)
    plot = newFigure()
    plot.plot(WsumSPH, nperh0, "r-*", label="Actual")
    plot.plot(WsumSPH, nperhSPH, "k-", label="Fit")
    plot.set_title(f"{Wstr} n per h as a function of $\sum W$ : SPH algorithm")
    plot.set_xlabel(r"$\sum W$")
    plot.set_ylabel("n per h")
    plot.legend()

    # ASPH fit for nperh(Wsum)
    plot = newFigure()
    plot.plot(WsumASPH, nperh0, "r-*", label="Actual")
    plot.plot(WsumASPH, nperhASPH, "k-", label="Fit")
    plot.set_title(f"{Wstr} n per h as a function of $\lambda(\psi)$ : ASPH algorithm")
    plot.set_xlabel(r"$\lambda(\psi)$")
    plot.set_ylabel("n per h")
    plot.legend()

    # SPH nperh
    plot = newFigure()
    plot.plot(nperh0, nperhSPH, "b*-", label="nperh lookup")
    plot.set_title(f"{Wstr} n per h lookup test : SPH algorithm")
    plot.set_xlabel("nperh actual")
    plot.set_ylabel("nperh estimated")

    # SPH nperh error
    errSPH = (nperhSPH - nperh0)/nperh0
    plot = newFigure()
    plot.plot(nperh0, errSPH, "r*-")
    plot.set_title(f"{Wstr} n per h lookup test error : SPH algorithm")
    plot.set_xlabel("nperh actual")
    plot.set_ylabel("Error")

    plot = newFigure()
    plot.plot(nperh0, nperhASPH, "b*-")
    plot.set_title(f"{Wstr} n per h lookup test : ASPH algorithm")
    plot.set_xlabel("nperh actual")
    plot.set_ylabel("nperh estimated")

    errASPH = (nperhASPH - nperh0)/nperh0
    plot = newFigure()
    plot.plot(nperh0, errASPH, "r*-")
    plot.set_title(f"{Wstr} n per h lookup test error : ASPH algorithm")
    plot.set_xlabel("nperh actual")
    plot.set_ylabel("Error")

    # Test ASPH with different aspect ratios
    if nDim == 2:
        aspect = np.arange(0.1, 1.0, 0.05)
        X, Y = np.meshgrid(nperh0, aspect)
        WsumASPHx = np.ndarray(X.shape)
        WsumASPHy = np.ndarray(X.shape)
        nperhASPHx = np.ndarray(X.shape)
        nperhASPHy = np.ndarray(X.shape)
        nperhASPHx_err = np.ndarray(X.shape)
        nperhASPHy_err = np.ndarray(X.shape)
        for iy in range(X.shape[0]):
            for ix in range(X.shape[1]):
                nPerhi = X[iy,ix]
                aspecti = Y[iy,ix]
                WsumASPHx[iy,ix] = sumKernelValuesSlice2d(WT, Vector2d(1,0), 1.0/nPerhi, aspecti/nPerhi)
                WsumASPHy[iy,ix] = sumKernelValuesSlice2d(WT, Vector2d(0,1), 1.0/nPerhi, aspecti/nPerhi)
                nperhASPHx[iy,ix] = WT.equivalentNodesPerSmoothingScaleASPH(WsumASPHx[iy,ix])
                nperhASPHy[iy,ix] = WT.equivalentNodesPerSmoothingScaleASPH(WsumASPHy[iy,ix])
                nperhASPHx_err[iy,ix] = (nperhASPHx[iy,ix] - nPerhi)/nPerhi
                nperhASPHy_err[iy,ix] = (nperhASPHy[iy,ix] - nPerhi/aspecti)/(nPerhi/aspecti)

        plotSurface(X, Y, WsumASPHx,
                    title = f"{Wstr} ASPH Wsum $X$",
                    xlabel = "n per h",
                    ylabel = "aspect ratio")
        plotSurface(X, Y, WsumASPHy,
                    title = f"{Wstr} ASPH Wsum $Y$",
                    xlabel = "n per h",
                    ylabel = "aspect ratio")
        plotSurface(X, Y, nperhASPHx,
                    title = f"{Wstr} ASPH n per h $X$",
                    xlabel = "n per h",
                    ylabel = "aspect ratio")
        plotSurface(X, Y, nperhASPHy,
                    title = f"{Wstr} ASPH n per h $Y$",
                    xlabel = "n per h",
                    ylabel = "aspect ratio")
        plotSurface(X, Y, nperhASPHx_err,
                    title = f"{Wstr} ASPH n per h $X$ error",
                    xlabel = "n per h",
                    ylabel = "aspect ratio")
        plotSurface(X, Y, nperhASPHy_err,
                    title = f"{Wstr} ASPH n per h $Y$ error",
                    xlabel = "n per h",
                    ylabel = "aspect ratio")
