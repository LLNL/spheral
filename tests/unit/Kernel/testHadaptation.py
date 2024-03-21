from Spheral2d import *

import numpy as np
from SpheralTestUtilities import *
from SpheralMatplotlib import *

#-------------------------------------------------------------------------------
# Command line options
#-------------------------------------------------------------------------------
commandLine(Kernel = WendlandC4Kernel,
            nPerh = 4.01,
            fscale = 1.0, 
            fscaleAngle = 0.0,  # degrees
            iterations = 10,
            startCorrect = False,
)
assert fscale <= 1.0
fscaleAngle *= pi/180.0

#-------------------------------------------------------------------------------
# Make the kernel and the ASPH update method
#-------------------------------------------------------------------------------
WT = TableKernel(Kernel())
etamax = WT.kernelExtent
asph = ASPHSmoothingScale(WT, targetNperh = nPerh, numPoints = 200)

#-------------------------------------------------------------------------------
# Generate our test point positions
#-------------------------------------------------------------------------------
# First generate the points in the reference frame (eta space)
nx = int(4.0*etamax*nPerh/fscale)
if nx % 2 == 0:
    nx += 1

xcoords = 2.0/fscale*np.linspace(-etamax, etamax, nx)
ycoords = 2.0/fscale*np.linspace(-etamax, etamax, nx)
X, Y = np.meshgrid(xcoords, ycoords)
eta_coords = np.column_stack((X.ravel(), Y.ravel()))

# Apply the inverse of the H transformation to map these eta coordinates to the
# lab frame
HCinv = SymTensor(fscale, 0.0,
                  0.0,    1.0)
R = rotationMatrix(Vector(cos(fscaleAngle), sin(fscaleAngle))).Transpose()
HCinv.rotationalTransform(R)
coords = np.copy(eta_coords)
for i in range(len(coords)):
    eta_coords[i] = R*Vector(eta_coords[i][0],eta_coords[i][1])
    coords[i] = HCinv * Vector(coords[i][0], coords[i][1])
HC = HCinv.Inverse()

#-------------------------------------------------------------------------------
# Function for plotting the current H tensor
#-------------------------------------------------------------------------------
def plotH(H, plot,
          style = "k-",
          etaSpace = False):
    etamax = WT.kernelExtent
    Hinv = H.Inverse()
    if etaSpace:
        Hinv = (HC*Hinv).Symmetric()
    t = np.linspace(0, 2.0*pi, 180)
    x = np.cos(t)
    y = np.sin(t)
    for i in range(len(x)):
        etav = Hinv*Vector(x[i], y[i]) * etamax
        x[i] = etav.x
        y[i] = etav.y
    plot.plot(x, y, style)
    return

#-------------------------------------------------------------------------------
# Function to measure the second moment tensor psi
#-------------------------------------------------------------------------------
def computePsi(coords, H, WT, nPerh):
    Wsum = 0.0
    psiLab = SymTensor()
    psiEta = SymTensor()
    for vals in coords:
        rji = Vector(*vals)
        eta = H*rji
        Wi = WT.kernelValueASPH(eta.magnitude(), nPerh)
        Wsum += Wi
        psiLab += Wi * rji.selfdyad()
        psiEta += Wi * eta.selfdyad()
    return Wsum, psiLab, psiEta

#-------------------------------------------------------------------------------
# Compute a new H based on the current second-moment (psi) and H
#-------------------------------------------------------------------------------
def newH(H0, Wsum, psiLab, psiEta, WT, nPerh, asph):
    H0inv = H0.Inverse()
    eigenLab = psiLab.eigenVectors()
    eigenEta = psiEta.eigenVectors()
    print("         Wsum : ", Wsum)
    print("       psiLab : ", psiLab)
    print("       psiEta : ", psiEta)
    print("     eigenLab : ", eigenLab)
    print("     eigenEta : ", eigenEta)

    # First the ASPH shape & volume change
    fnu = [1.0, 1.0]
    fscale = 1.0
    T = SymTensor(1.0, 0.0, 0.0, 1.0)
    for nu in range(2):
        lambdaPsi = sqrt(eigenEta.eigenValues[nu])
        evec = eigenEta.eigenVectors.getColumn(nu)
        nPerheff = asph.equivalentNodesPerSmoothingScale(lambdaPsi)
        T(nu, nu, max(0.75, min(1.25, nPerh/nPerheff)))
        print("      --> evec, nPerheff : ", evec, nPerheff)


        # h0 = (H0inv*evec).magnitude()
        # thpt = sqrt((psiEta*evec).magnitude())
        # nPerheff = asph.equivalentNodesPerSmoothingScale(thpt)
        # print("      --> h0, nPerheff : ", h0, nPerheff)
        # fnu[nu] = nPerh/nPerheff
        # fscale *= nPerh/nPerheff
        # H1inv(nu,nu, h0 * nPerh/nPerheff)
    print("         T before SPH scaling: ", T)

    # Share the SPH volume change estimate by the ratio of the eigenvalue scaling
    nPerhSPH = WT.equivalentNodesPerSmoothingScale(sqrt(Wsum))
    fscale = nPerh/nPerhSPH / sqrt(fscale)
    T[0] *= fscale*sqrt(fnu[0]/fnu[1])
    T[2] *= fscale*sqrt(fnu[1]/fnu[0])
    print("          T after SPH scaling: ", T)

    T.rotationalTransform(eigenEta.eigenVectors)
    print("         T final: ", T)
    H1inv = (T*H0inv).Symmetric()
    return H1inv.Inverse()

#-------------------------------------------------------------------------------
# Plot the initial point distribution and H
#-------------------------------------------------------------------------------
if startCorrect:
    H = SymTensor(HC)
else:
    H = SymTensor(1.0, 0.0,
                  0.0, 1.0)
    H *= 2.0   # Make it too small to start
print("Initial H tensor (inverse): ", H.Inverse())

# Plot the initial point distribution in lab coordinates
plotLab = newFigure()
plotLab.set_box_aspect(1.0)
plotLab.plot([x[0] for x in coords], [x[1] for x in coords], "ro")
plotH(H, plotLab, "k-")
plim = max(abs(np.min(coords)), np.max(coords))
plotLab.set_xlim(-plim, plim)
plotLab.set_ylim(-plim, plim)
plotLab.set_xlabel(r"$x$")
plotLab.set_ylabel(r"$y$")
plotLab.set_title("Lab frame")

# Plot in eta space
plotEta = newFigure()
plotEta.set_box_aspect(1.0)
plotEta.plot([x[0] for x in eta_coords], [x[1] for x in eta_coords], "ro")
plotH(H, plotEta, "k-", True)
plim = max(abs(np.min(eta_coords)), np.max(eta_coords))
plotEta.set_xlim(-plim, plim)
plotEta.set_ylim(-plim, plim)
plotEta.set_xlabel(r"$\eta_x$")
plotEta.set_ylabel(r"$\eta_y$")
plotEta.set_title("$\eta$ frame")

#-------------------------------------------------------------------------------
# Iterate on relaxing H
#-------------------------------------------------------------------------------
for iter in range(iterations):
    print("Iteration ", iter)
    Wsum, psiLab, psiEta = computePsi(coords, H, WT, nPerh)
    print("     Wsum, psiLab, psiEta: ", Wsum, psiLab, psiEta)
    #H = asph.idealSmoothingScale(H, Vector(0,0), 0.0, psi, WT, 1e-10, 1e10, 1e-10, nPerh, ConnectivityMap(), 0, 0)
    H = newH(H, Wsum, psiLab, psiEta, WT, nPerh, asph)
    evals = H.eigenValues()
    aspectRatio = evals.maxElement()/evals.minElement()
    output("     H.Inverse(), aspectRatio")
    plotH(H, plotLab, "b-")
    plotH(H, plotEta, "b-", True)
