from Spheral import *
from GenerateSphericalNodeProfile1d import *
from SortAndDivideRedistributeNodes import distributeNodes1d
from SpheralTestUtilities import *

import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# The set of r/h values from the origin we'll test
etavals = (1e-5, 1e-2, 2e-2, 0.5, 1.5, 2.5, 3.5, 10.0, 20.0)
h = 0.1

#-------------------------------------------------------------------------------
# Kernels
#-------------------------------------------------------------------------------
W1 = NBSplineKernel1d(5)   # 5th order b-spline
t0 = time.time()
Wr = SphericalRadialKernel(W1)
Wr0 = SphericalRadialKernel(W1, useInterpolation=False)
print("Required %0.4f sec to construct SphericalRadialKernel" % (time.time() - t0))
# t0 = time.time()
# Wo = SphericalKernelOslo(W1)
# print("Required %0.4f sec to construct SphericalKernelOslo" % (time.time() - t0))

#-------------------------------------------------------------------------------
# Plot what the A interpolation space looks like
#-------------------------------------------------------------------------------
x = np.arange(0.0, 2*Wr.etacutoff, 2*Wr.etacutoff/200)
y = np.array([Wr.volumeNormalization(xi) for xi in x])
y0 = np.array([Wr0.volumeNormalization(xi) for xi in x])

figA = plt.figure(tight_layout=True, figsize=(8,8))
gs = gridspec.GridSpec(nrows = 2, ncols = 1, height_ratios = [2,1], figure=figA)

ax1 = figA.add_subplot(gs[0,0])
ax1.semilogy(x, y, label="A")
ax1.set_xlabel("$r/h$")
ax1.set_ylabel("$A$")
ax1.set_title("Volume normalization for SphericalRadialKernel")

ax2 = figA.add_subplot(gs[1,0])
ax2.semilogy(x, np.abs((y - y0)/y0), label="Interpolation error")
ax2.set_xlabel("$r/h$")
ax2.set_ylabel(r"$\left( \langle A \rangle - A \right)/A$")

#-------------------------------------------------------------------------------
# Plot what the grad A interpolation space looks like
#-------------------------------------------------------------------------------
y = np.array([Wr.gradAInv(xi) for xi in x])
y0 = np.array([Wr0.gradAInv(xi) for xi in x])

figGradAInv = plt.figure(tight_layout=True, figsize=(8,8))
gs = gridspec.GridSpec(nrows = 2, ncols = 1, height_ratios = [2,1], figure=figGradAInv)

ax1 = figGradAInv.add_subplot(gs[0,0])
ax1.semilogy(x, y, label=r"\partial_{\eta} A^{-1}")
ax1.set_xlabel("$r/h$")
ax1.set_ylabel(r"$\partial_{\eta} A^{-1}$")
ax1.set_title("Gradient of one over the volume normalization for SphericalRadialKernel")

ax2 = figGradAInv.add_subplot(gs[1,0])
ax2.semilogy(x, np.abs((y - y0)/y0), label="Interpolation error")
ax2.set_xlabel("$r/h$")
ax2.set_ylabel(r"$\left( \langle \partial_{\eta}A^{-1} \rangle - \partial_{\eta} A^{-1} \right)/\partial_{\eta} A^{-1}$")

#-------------------------------------------------------------------------------
# Try interpolating some functions
#-------------------------------------------------------------------------------
r0, r1 = 0.0, 2.0
nr = 100
nPerh = 1.35

# Make our test nodes
eos = GammaLawGasMKS1d(5.0/3.0, 1.0)
nodes = makeFluidNodeList("nodes", eos, 
                           nPerh = nPerh,
                           kernelExtent = Wr.kernelExtent)
gen = GenerateSphericalNodeProfile1d(nr = nr,
                                     rho = 2.0,
                                     rmin = r0,
                                     rmax = r1,
                                     nNodePerh = nPerh)
distributeNodes1d((nodes, gen))

F0, slope = 10.0, 0.0
def field_func(ri):
    return F0 + slope*ri

