from Spheral import *
from GenerateSphericalNodeProfile1d import *
from SortAndDivideRedistributeNodes import distributeNodes1d
from SpheralTestUtilities import *
from SpheralMatplotlib import *

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
t0 = time.time()
W1 = NBSplineKernel1d(5)   # 5th order b-spline
Wr = SphericalRadialKernel(W1)
Wr0 = SphericalRadialKernel(W1, useInterpolation=False)
print("Required %0.4f sec to construct SphericalRadialKernel" % (time.time() - t0))
# t0 = time.time()
# W3 = TableKernel3d(NBSplineKernel3d(5))
# Wo = SphericalKernelOslo(W3)
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
# Test volume integrals over the kernel at various distances from the origin
#-------------------------------------------------------------------------------
def integrate_Wr(R):
    class Wvol_shell(ScalarScalarFunctor):
        def __init__(self):
            ScalarScalarFunctor.__init__(self)
            return
        def __call__(self, ri):
            return ri*ri*Wr0(Vector1d(R/h), Vector1d(ri/h), 1.0/h)
    return simpsonsIntegrationDouble(Wvol_shell(), R - h*Wr0.etamax, R + h*Wr0.etamax, 200)

rvals = np.linspace(0.0, 20.0*h, 80)
etavals = rvals/h
Wrsumvals = np.array([integrate_Wr(etai * h) for etai in etavals])

figWrsum = plt.figure(tight_layout=True, figsize=(8,8)).add_subplot(111)
figWrsum.plot(etavals, Wrsumvals, "k-")
figWrsum.plot(etavals, Wrsumvals, "r*")
figWrsum.set_xlabel(r"$\eta = R/h$")
figWrsum.set_ylabel(r"$\int W_r \; dV$")
figWrsum.set_title(r"Integral of $W_r$ over volume")

#-------------------------------------------------------------------------------
# Try interpolating some functions
#-------------------------------------------------------------------------------
raise ValueError("Stop me!")
r0, r1 = 0.0, 2.0
nr = 100
nPerh = 1.35

# Make our test nodes
eos = GammaLawGasMKS1d(5.0/3.0, 1.0)
nodes = makeFluidNodeList1d("nodes", eos, 
                            nPerh = nPerh,
                            kernelExtent = Wr.etamax)
gen = GenerateSphericalNodeProfile1d(nr = nr,
                                     rho = 2.0,
                                     rmin = r0,
                                     rmax = r1,
                                     nNodePerh = nPerh)
distributeNodes1d((nodes, gen))
db = DataBase1d()
db.appendNodeList(nodes)
db.reinitializeNeighbors()
pos = nodes.positions()
mass = nodes.mass()
rho = nodes.massDensity()
H = nodes.Hfield()

# Make ghost points across the origin
bc = ReflectingBoundary1d(Plane1d(Vector1d(0.0), Vector1d(1.0)))
bc.setAllGhostNodes(db)
bc.finalizeGhostBoundary()

# Find the connectivity
db.reinitializeNeighbors()
cm = db.connectivityMap(False, False, False)

# Initialize a linear field
F0, slope = 10.0, 0.0
def linear_func(ri):
    return F0 + slope*ri

linear_field = ScalarField1d("linear field", nodes)
for i in range(nodes.numNodes):
    linear_field[i] = linear_func(pos[i].x)

# Test interpolation
scatter_interp = ScalarField1d("scatter interpolation of linear field", nodes)
gather_interp = ScalarField1d("gather interpolation of linear field", nodes)
sym_interp = ScalarField1d("gather/scatter interpolation of linear field", nodes)
scatter_weight_sum = ScalarField1d("scatter sum of interpolation weights", nodes)
gather_weight_sum = ScalarField1d("gather sum of interpolation weights", nodes)
sym_weight_sum = ScalarField1d("sym sum of interpolation weights", nodes)
def increment_interp_values(i, j, Wi, Wj):
    Vj = mass[j]/rho[j]
    wscatter = Vj*Wj
    wgather = Vj*Wi
    wsym = 0.5*Vj*(Wi + Wj)
    scatter_weight_sum[i] += wscatter
    gather_weight_sum[i] += wgather
    sym_weight_sum[i] += wsym
    scatter_interp[i] += wscatter * linear_field[j]
    gather_interp[i] += wgather * linear_field[j]
    sym_interp[i] += wsym * linear_field[j]
    return
for i in range(nodes.numInternalNodes):
    Wi = Wr(pos[i], pos[i], H[i].Determinant())
    increment_interp_values(i, i, Wi, Wi)
for pair in cm.nodePairList:
    i, j = pair.i_node, pair.j_node
    Wi = Wr(pos[i], pos[j], H[i].Determinant())
    Wj = Wr(pos[j], pos[i], H[j].Determinant())
    increment_interp_values(i, j, Wi, Wj)
    increment_interp_values(j, i, Wj, Wi)

fig_weight_sum = plt.figure(tight_layout=True, figsize=(8,8)).add_subplot(111)
plotField(ScalarField1d("unit", nodes, 1.0),
          plot = fig_weight_sum,
          plotStyle = "k-",
          lineTitle = "Analytic",
          winTitle = "Interpolation weight sum")
plotField(scatter_weight_sum,
          plot = fig_weight_sum,
          plotStyle = "ro",
          lineTitle = "Scatter")
plotField(gather_weight_sum,
          plot = fig_weight_sum,
          plotStyle = "b+",
          lineTitle = "Gather")
plotField(sym_weight_sum,
          plot = fig_weight_sum,
          plotStyle = "g^",
          lineTitle = "Gather/Scatter")

fig_linear_interp = plt.figure(tight_layout=True, figsize=(8,8)).add_subplot(111)
plotField(linear_field,
          plot = fig_linear_interp,
          plotStyle = "k-",
          lineTitle = "Analytic",
          winTitle = "Linear interpolation")
plotField(scatter_interp,
          plot = fig_linear_interp,
          plotStyle = "ro",
          lineTitle = "Scatter")
plotField(gather_interp,
          plot = fig_linear_interp,
          plotStyle = "b+",
          lineTitle = "Gather")
plotField(sym_interp,
          plot = fig_linear_interp,
          plotStyle = "g^",
          lineTitle = "Gather/Scatter")

