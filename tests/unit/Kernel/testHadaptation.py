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
# Function to measure the second moment tensor psi
#-------------------------------------------------------------------------------
def computePsi(coords, H, WT, nPerh):
    Wsum = 0.0
    psiLab = SymTensor()
    psiEta = SymTensor()
    for vals in coords:
        rji = Vector(*vals)
        eta = H*rji
        WSPHi = WT.kernelValueSPH(eta.magnitude())
        WASPHi = WT.kernelValueASPH(eta.magnitude(), nPerh)
        Wsum += WSPHi
        psiLab += WSPHi**2 * rji.unitVector().selfdyad()
        psiEta += WSPHi**2 * eta.unitVector().selfdyad()
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

    # Extract shape information from the second moment
    H1inv = SymTensor(H0inv)
    nperheff = WT.equivalentNodesPerSmoothingScale(sqrt(Wsum))
    fscale = nPerh/nperheff if nperheff > 0.0 else 2.0
    T = psiEta.Inverse().sqrt() if psiEta.Determinant() > 0.0 else SymTensor(1, 0, 0, 1)
    T *= fscale/sqrt(T.Determinant())
    eigenT = T.eigenVectors()
    if eigenT.eigenValues.minElement() < 0.25 or eigenT.eigenValues.maxElement() > 4.0:
        T = SymTensor(min(4.0, max(0.25, eigenT.eigenValues[0])), 0.0,
                      0.0, min(4.0, max(0.25, eigenT.eigenValues[1])))
        T.rotationalTransform(eigenT.eigenVectors)
    H1inv = (T*H0inv).Symmetric()
    print("     nperheff : ", nperheff)
    print("            T : ", T)
    print("        H0inv : ", H0inv)
    print("        H1inv : ", H1inv)

    # # Share the SPH volume change estimate by the ratio of the eigenvalue scaling
    # nPerhSPH = WT.equivalentNodesPerSmoothingScale(sqrt(Wsum))
    # fscale = nPerh/nPerhSPH / sqrt(fscale)
    # T[0] *= fscale*sqrt(fnu[0]/fnu[1])
    # T[2] *= fscale*sqrt(fnu[1]/fnu[0])
    # print("          T after SPH scaling: ", T)

    # T.rotationalTransform(eigenEta.eigenVectors)
    # print("         T final: ", T)
    # H1inv = (T*H0inv).Symmetric()
    return H1inv.Inverse()

# def newH(H0, coords, inv_coords, WT, nPerh, asph):
#     H0inv = H0.Inverse()

#     # Compute the inverse hull to find the nearest neighbors
#     hull0 = Polygon(inv_coords)
    
#     # Build a normal space hull using hull0's points and their reflections
#     verts = [x.unitVector()*safeInv(x.magnitude()) for x in hull0.vertices]
#     verts += [-x for x in verts]
#     hull1 = Polygon(verts)

#     # Extract the second-moment from the hull
#     psi = sum([x.selfdyad() for x in hull1.vertices], SymTensor())

#     # Find the new H shape
#     D0 = psi.Determinant()
#     assert D0 > 0.0
#     psi /= sqrt(D0)
#     Hnew = psi.sqrt().Inverse()
#     assert np.isclose(Hnew.Determinant(), 1.0)

#     # Compute the zeroth moment
#     Wzero = sqrt(sum([WT.kernelValueSPH((H0*Vector(*c)).magnitude()) for c in coords]))

#     # What is the current effect nPerh?
#     currentNodesPerSmoothingScale = WT.equivalentNodesPerSmoothingScale(Wzero);
#     assert currentNodesPerSmoothingScale > 0.0

#     # The (limited) ratio of the desired to current nodes per smoothing scale.
#     s = min(4.0, max(0.25, nPerh/(currentNodesPerSmoothingScale + 1.0e-30)))
#     assert s > 0.0

#     # Scale to the desired determinant
#     Hnew *= sqrt(H0.Determinant())/s

#     print("        Wzero : ", Wzero)
#     print("        hull0 : ", hull0.vertices)
#     print("        hull1 : ", hull1.vertices)
#     print("          psi : ", psi)
#     print("    psi Eigen : ", psi.eigenVectors())
#     print("     nPerheff : ", currentNodesPerSmoothingScale)
#     print("           H0 : ", H0)
#     print("           H1 : ", Hnew)
#     return Hnew, hull1

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
    #H = asph.idealSmoothingScale(H, Vector(0,0), 0.0, psi, WT, 1e-10, 1e10, 1e-10, nPerh, ConnectivityMap(), 0, 0)
    Wsum, psiLab, psiEta = computePsi(coords, H, WT, nPerh)
    H = newH(H, Wsum, psiLab, psiEta, WT, nPerh, asph)
    evals = H.eigenValues()
    aspectRatio = evals.maxElement()/evals.minElement()
    output("     H.Inverse(), aspectRatio")
    plotH(H, plotLab, "b-")
    plotH(H, plotEta, "b-", True)
    #plotPolygon(hull, plot=plotHull)
