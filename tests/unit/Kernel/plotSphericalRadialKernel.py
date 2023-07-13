from Spheral import *

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
# Plot what the A interpolation space looks like
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

# #-------------------------------------------------------------------------------
# # Reproduce Fig 1 from Omang, M., Borve, S., & Trulsen, J. (2006)
# #-------------------------------------------------------------------------------
# fig1 = plt.figure(tight_layout=True, figsize=(10,8))
# gs = gridspec.GridSpec(nrows = 2, ncols = 2, height_ratios = [2,1], figure=fig1)

# # First plot the SphericalTabelKernel fit
# ax = fig1.add_subplot(gs[0,0])
# for eta in etavals:
#     r = h*eta
#     rp = rprange(r, h)
#     yvals = np.array([W(Vector1d(r/h), Vector1d(rpi/h), 1.0/h) for rpi in rp])
#     yvals *= r*r
#     ax.plot((rp - r)/h, yvals, label = r"$r/h=%g$" % eta)
# ax.set_xlabel(r"$(r^\prime - r)/h$")
# ax.set_ylabel(r"$r^2 \langle W_{3S1}(r^\prime, r, h)/h$ \rangle")
# ax.set_title("SphericalKernelOslo approximation")
# legend = ax.legend(loc="upper right", shadow=True)

# # Analytic kernel
# ax = fig1.add_subplot(gs[0,1])
# for eta in etavals:
#     r = h*eta
#     rp = rprange(r, h)
#     yvals = np.array([W3S1(rpi, r, h) for rpi in rp])
#     yvals *= r*r
#     ax.plot((rp - r)/h, yvals, label = r"$r/h=%g$" % eta)
# ax.set_xlabel(r"$(r^\prime - r)/h$")
# ax.set_ylabel(r"$r^2 W_{3S1}(r^\prime, r, h)/h$")
# ax.set_title("Analytic")

# # Kernel error
# ax = fig1.add_subplot(gs[1,:])
# for eta in etavals:
#     r = h*eta
#     rp = rprange(r, h)
#     yvals = np.array([error(W3S1(rpi, r, h), W(Vector1d(r/h), Vector1d(rpi/h), 1.0/h))for rpi in rp])
#     ax.semilogy((rp - r)/h, yvals, label = r"$r/h=%g$" % eta)
# ax.set_xlabel(r"$(r^\prime - r)/h$")
# ax.set_ylabel(r"$|\langle W_{3S1}(r^\prime, r, h) \rangle/W_{3S1}(r^\prime, r, h) - 1|$")
# ax.set_title("Error")

# #-------------------------------------------------------------------------------
# # Plot the gradient compared with the numpy gradient estimator
# #-------------------------------------------------------------------------------
# fig20 = plt.figure(tight_layout=True, figsize=(10,8))
# gs = gridspec.GridSpec(nrows = 2, ncols = 2, height_ratios = [2,1], figure=fig1)

# # Plot SphericalKernelOslo gradient
# ax = fig20.add_subplot(gs[0,0])
# for eta in etavals:
#     r = h*eta
#     rp = rprange(r, h)
#     gyvals = np.array([W.grad(Vector1d(rpi/h), Vector1d(r/h), SymTensor1d(1.0/h)).x for rpi in rp])
#     gyvals *= r*r
#     ax.plot((rp - r)/h, gyvals, label = r"$r/h=%g$" % eta)
# ax.set_xlabel(r"$(r^\prime - r)/h$")
# ax.set_ylabel(r"$r^2 \; \langle \partial_r W_{3S1}(r^\prime, r, h) \rangle$")
# ax.set_title("SphericalKernelOslo gradient approximation")
# legend = ax.legend(loc="upper right", shadow=True)

# # Analytic gradient
# ax = fig20.add_subplot(gs[0,1])
# for eta in etavals:
#     r = h*eta
#     rp = rprange(r, h)
#     gyvals = np.array([gradW3S1(rpi, r, h) for rpi in rp])
#     gyvals *= r*r
#     ax.plot((rp - r)/h, gyvals, label = r"$r/h=%g$" % eta)
# ax.set_xlabel(r"$(r^\prime - r)/h$")
# ax.set_ylabel(r"$r^2 \; \partial_r W_{3S1}(r^\prime, r, h)/h$")
# ax.set_title("Analytic gradient")

# # Kernel gradient error
# ax = fig20.add_subplot(gs[1,:])
# for eta in etavals:
#     r = h*eta
#     rp = rprange(r, h, etastep=0.05)
#     gyvals = np.array([W.grad(Vector1d(rpi/h), Vector1d(r/h), SymTensor1d(1.0/h)).x for rpi in rp])
#     gyvals0 = np.array([gradW3S1(rpi, r, h) for rpi in rp])
#     errvals = np.array([error(gyvals0[i], gyvals[i]) for i in range(len(gyvals))])
#     ax.semilogy((rp - r)/h, errvals, label = r"$r/h=%g$" % eta)
# ax.set_xlabel(r"$(r^\prime - r)/h$")
# ax.set_ylabel(r"$|\langle \partial_r W_{3S1}(r^\prime, r, h) \rangle - \partial_r W_{3S1}(r^\prime, r, h)|/|\partial_r W_{3S1}(r^\prime, r, h)|$")
# ax.set_title("gradient Error")

plt.show()
