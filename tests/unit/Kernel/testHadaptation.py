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

def safeInv(x, fuzz=1e-30):
    return x/(x*x + fuzz)

#-------------------------------------------------------------------------------
# Make the kernel and the ASPH update method
#-------------------------------------------------------------------------------
WT = TableKernel(Kernel())
etamax = WT.kernelExtent
asph = ASPHSmoothingScale()

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

# inverse coordinates (squared)
inv_coords = [Vector(*c).unitVector()*safeInv(Vector(*c).magnitude()) for c in coords]

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
# Function to measure the moments of the local point distribution
#-------------------------------------------------------------------------------
def computeMoments(coords, H, WT, nPerh):
    zerothMoment = 0.0
    firstMoment = Vector()
    secondMoment = SymTensor()
    for vals in coords:
        rji = Vector(*vals)
        eta = H*rji
        WSPHi = WT.kernelValueSPH(eta.magnitude())
        zerothMoment += WSPHi
        firstMoment += WSPHi * eta
        secondMoment += WSPHi**2 * eta.unitVector().selfdyad()
    return zerothMoment, firstMoment, secondMoment

#-------------------------------------------------------------------------------
# Compute a new H based on the current second-moment (psi) and H
#-------------------------------------------------------------------------------
def newH(H0, zerothMoment, firstMoment, secondMoment, WT, nPerh):
    print(" zerothMoment : ", zerothMoment)
    print("  firstMoment : ", firstMoment)
    print(" secondMoment : ", secondMoment)

    # Extract shape information from the second moment
    nperheff = WT.equivalentNodesPerSmoothingScale(sqrt(zerothMoment))
    T = secondMoment.sqrt()
    print("     nperheff : ", nperheff)
    print("           T0 : ", T)
    eigenT = T.eigenVectors()
    Tmax = max(1.0, eigenT.eigenValues.maxElement())
    fscale = 1.0
    for j in range(2):
        eigenT.eigenValues[j] = max(eigenT.eigenValues[j], 1e-2*Tmax)
        fscale *= eigenT.eigenValues[j]
    assert fscale > 0.0
    fscale = 1.0/sqrt(fscale)
    fscale *= min(4.0, max(0.25, nperheff/nPerh))  # inverse length, same as H!
    T = SymTensor(fscale*eigenT.eigenValues[0], 0.0,
                  0.0, fscale*eigenT.eigenValues[1])
    T.rotationalTransform(eigenT.eigenVectors)
    H1 = (T*H0).Symmetric()
    print("         Tfin : ", T)
    print("        H0inv : ", H0.Inverse())
    print("        H1inv : ", H1.Inverse())
    return H1

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

# # Plot for the hulls in lab coordinates
# plotHull = newFigure()
# plotHull.set_box_aspect(1.0)
# plotHull.plot([x[0] for x in coords], [x[1] for x in coords], "ro")
# plim = max(abs(np.min(coords)), np.max(coords))
# plotHull.set_xlim(-plim, plim)
# plotHull.set_ylim(-plim, plim)
# plotHull.set_xlabel(r"$x$")
# plotHull.set_ylabel(r"$y$")
# plotHull.set_title("Lab frame (Hull)")

#-------------------------------------------------------------------------------
# Iterate on relaxing H
#-------------------------------------------------------------------------------
for iter in range(iterations):
    print("Iteration ", iter)
    zerothMoment, firstMoment, secondMoment = computeMoments(coords, H, WT, nPerh)
    H = newH(H, zerothMoment, firstMoment, secondMoment, WT, nPerh)
    # H = asph.idealSmoothingScale(H = H,
    #                              pos = Vector(),
    #                              zerothMoment = sqrt(Wsum),
    #                              firstMoment = Vector(),
    #                              secondMomentEta = psiEta,
    #                              secondMomentLab = psiEta,
    #                              W = WT,
    #                              hmin = 1e-10,
    #                              hmax = 1e10,
    #                              hminratio = 1e-10,
    #                              nPerh = nPerh,
    #                              connectivityMap = ConnectivityMap(),
    #                              nodeListi = 0,
    #                              i = 0)
    evals = H.eigenValues()
    aspectRatio = evals.maxElement()/evals.minElement()
    output("     H.Inverse(), aspectRatio")
    plotH(H, plotLab, "b-")
    plotH(H, plotEta, "b-", True)
    #plotPolygon(hull, plot=plotHull)
