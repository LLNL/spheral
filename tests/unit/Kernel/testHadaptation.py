from Spheral2d import *

import numpy as np
from SpheralTestUtilities import *
from SpheralMatplotlib import *

#-------------------------------------------------------------------------------
# Command line options
#-------------------------------------------------------------------------------
commandLine(Kernel = WendlandC4Kernel,
            nPerh = 4.01,
            rotation = 0.0,
            xscale = 1.0,
            yscale = 1.0,
            iterations = 10,
            startCorrect = False,
)

#-------------------------------------------------------------------------------
# Make the kernel and the ASPH update method
#-------------------------------------------------------------------------------
WT = TableKernel(Kernel())
asph = ASPHSmoothingScalev2(WT, targetNperh = nPerh)

#-------------------------------------------------------------------------------
# Generate our test point positions
#-------------------------------------------------------------------------------
fxscale = max(1.0, yscale/xscale)
fyscale = max(1.0, xscale/yscale)
nx = int(4.0*nPerh * fxscale)
ny = int(4.0*nPerh * fyscale)

# Make sure we have a point at (0, 0)
if nx % 2 == 0:
    nx += 1
if ny % 2 == 0:
    ny += 1

dx = 2.0/(nx - 1)
dy = 2.0/(ny - 1)

xcoords = np.linspace(-1.0, 1.0, nx)
ycoords = np.linspace(-1.0, 1.0, ny)
assert xcoords[(nx - 1)//2] == 0.0
assert ycoords[(ny - 1)//2] == 0.0

#-------------------------------------------------------------------------------
# Function for plotting the current H tensor
#-------------------------------------------------------------------------------
def plotH(H, plot, style="k-"):
    Hinv = WT.kernelExtent * H.Inverse()
    t = np.linspace(0, 2.0*pi, 180)
    x = np.cos(t)
    y = np.sin(t)
    for i in range(len(x)):
        etav = Hinv*Vector(x[i], y[i])
        x[i] = etav.x
        y[i] = etav.y
    plot.plot(x, y, style)
    return

#-------------------------------------------------------------------------------
# Function to measure the second moment tensor psi
#-------------------------------------------------------------------------------
def computePsi(x, y, H, WT, nPerh):
    nx = len(x)
    ny = len(y)
    Wsum = 0.0
    psiLab = SymTensor()
    psiEta = SymTensor()
    for j in range(ny):
        for i in range(nx):
            rji = Vector(x[i], y[j])
            eta = H*rji
            Wsum += WT.kernelValueSPH(eta.magnitude())
            Wi = WT.kernelValueASPH(eta.magnitude(), nPerh)
            psiLab += Wi * rji.selfdyad()
            psiEta += Wi * eta.selfdyad()
    return Wsum, psiLab, psiEta

#-------------------------------------------------------------------------------
# Compute a new H based on the current second-moment (psi) and H
#-------------------------------------------------------------------------------
def newH(H0, Wsum, psiLab, psiEta, WT, asph, nPerh):
    H0inv = H0.Inverse()
    eigenLab = psiLab.eigenVectors()
    eigenEta = psiEta.eigenVectors()
    print("     eigenLab : ", eigenLab)
    print("     eigenEta : ", eigenEta)

    # First the ASPH shape & volume change
    H1inv = SymTensor()
    fscale = 1.0
    for nu in range(2):
        evec = eigenLab.eigenVectors.getColumn(nu)
        h0 = (H0inv*evec).magnitude()
        thpt = sqrt((psiEta*evec).magnitude())
        nPerheff = asph.equivalentNodesPerSmoothingScale(thpt)
        print("      --> h0, nPerheff : ", h0, nPerheff)
        fscale *= nPerh/nPerheff
        H1inv(nu,nu, h0 * nPerh/nPerheff)

    # Scale by the zeroth moment to get the right overall volume
    print("         H1inv before SPH scaling: ", H1inv)
    nPerhSPH = WT.equivalentNodesPerSmoothingScale(Wsum)
    fscale = nPerh/nPerhSPH / sqrt(fscale)
    H1inv *= fscale
    print("          H1inv after SPH scaling: ", H1inv)

    H1inv.rotationalTransform(eigenLab.eigenVectors)
    return H1inv.Inverse()

#-------------------------------------------------------------------------------
# Plot the initial point distribution and H
#-------------------------------------------------------------------------------
if startCorrect:
    H = SymTensor(1.0/(nPerh*dx), 0.0,
                  0.0, 1.0/(nPerh*dy))
else:
    H = SymTensor(1.0/(nPerh*dx*fxscale), 0.0,
                  0.0, 1.0/(nPerh*dy*fyscale))
    H *= 2.0   # Make it too small to start
print("Initial H tensor (inverse): ", H.Inverse())

# Plot the initial point distribution
plot = newFigure()
plot.set_box_aspect(1.0)
X, Y = np.meshgrid(xcoords, ycoords)
plot.plot(X, Y, "ro")
plotH(H, plot, "k-")

#-------------------------------------------------------------------------------
# Iterate on relaxing H
#-------------------------------------------------------------------------------
for iter in range(iterations):
    print("Iteration ", iter)
    Wsum, psiLab, psiEta = computePsi(xcoords, ycoords, H, WT, nPerh)
    print("     Wsum, psiLab, psiEta: ", Wsum, psiLab, psiEta)
    #H = asph.idealSmoothingScale(H, Vector(0,0), 0.0, psi, WT, 1e-10, 1e10, 1e-10, nPerh, ConnectivityMap(), 0, 0)
    H = newH(H, Wsum, psiLab, psiEta, WT, asph, nPerh)
    evals = H.eigenValues()
    aspectRatio = evals.maxElement()/evals.minElement()
    output("     H.Inverse(), aspectRatio")
    plotH(H, plot, "b-")
