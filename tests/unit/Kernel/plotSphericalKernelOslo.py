from Spheral import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from testBicubicSphericalKernelOslo import W3S1, gradW3S1, rprange, error, W

import time

# The set of r/h values from the origin we'll test
etavals = (1e-5, 1e-2, 2e-2, 0.5, 1.5, 2.5, 3.5, 10.0, 20.0)
h = 0.1

#-------------------------------------------------------------------------------
# Plot what the interpolation space looks like
#-------------------------------------------------------------------------------
fig0 = plt.figure(tight_layout=True, figsize=(8,8))
x = np.arange(0.0, W.baseKernel3d.kernelExtent, W.baseKernel3d.kernelExtent/100)
y = np.arange(0.0, W.baseKernel3d.kernelExtent, W.baseKernel3d.kernelExtent/100)
x, y = np.meshgrid(x, y)
NX, NY = x.shape
z0 = np.array([[W.Winterpolator(x[j][i], y[j][i]) for j in range(NY)] for i in range(NX)])
ax0 = fig0.add_subplot(111, projection='3d')
surf0 = ax0.plot_surface(x, y, z0, cmap=cm.coolwarm,
                         linewidth=0, antialiased=False)

#-------------------------------------------------------------------------------
# Reproduce Fig 1 from Omang, M., Borve, S., & Trulsen, J. (2006)
#-------------------------------------------------------------------------------
fig1 = plt.figure(tight_layout=True, figsize=(10,8))
gs = gridspec.GridSpec(nrows = 2, ncols = 2, height_ratios = [2,1], figure=fig1)

# First plot the SphericalTabelKernel fit
ax = fig1.add_subplot(gs[0,0])
for eta in etavals:
    r = h*eta
    rp = rprange(r, h)
    yvals = np.array([W(Vector1d(r/h), Vector1d(rpi/h), 1.0/h) for rpi in rp])
    yvals *= r*r
    ax.plot((rp - r)/h, yvals, label = r"$r/h=%g$" % eta)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$r^2 \langle W_{3S1}(r^\prime, r, h)/h$ \rangle")
ax.set_title("SphericalKernelOslo approximation")
legend = ax.legend(loc="upper right", shadow=True)

# Analytic kernel
ax = fig1.add_subplot(gs[0,1])
for eta in etavals:
    r = h*eta
    rp = rprange(r, h)
    yvals = np.array([W3S1(rpi, r, h) for rpi in rp])
    yvals *= r*r
    ax.plot((rp - r)/h, yvals, label = r"$r/h=%g$" % eta)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$r^2 W_{3S1}(r^\prime, r, h)/h$")
ax.set_title("Analytic")

# Kernel error
ax = fig1.add_subplot(gs[1,:])
for eta in etavals:
    r = h*eta
    rp = rprange(r, h)
    yvals = np.array([error(W3S1(rpi, r, h), W(Vector1d(r/h), Vector1d(rpi/h), 1.0/h))for rpi in rp])
    ax.semilogy((rp - r)/h, yvals, label = r"$r/h=%g$" % eta)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$|\langle W_{3S1}(r^\prime, r, h) \rangle/W_{3S1}(r^\prime, r, h) - 1|$")
ax.set_title("Error")

#-------------------------------------------------------------------------------
# Plot the gradient compared with the numpy gradient estimator
#-------------------------------------------------------------------------------
fig20 = plt.figure(tight_layout=True, figsize=(10,8))
gs = gridspec.GridSpec(nrows = 2, ncols = 2, height_ratios = [2,1], figure=fig1)

# Plot SphericalKernelOslo gradient
ax = fig20.add_subplot(gs[0,0])
for eta in etavals:
    r = h*eta
    rp = rprange(r, h)
    gyvals = np.array([W.grad(Vector1d(rpi/h), Vector1d(r/h), SymTensor1d(1.0/h)).x for rpi in rp])
    gyvals *= r*r
    ax.plot((rp - r)/h, gyvals, label = r"$r/h=%g$" % eta)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$r^2 \; \langle \partial_r W_{3S1}(r^\prime, r, h) \rangle$")
ax.set_title("SphericalKernelOslo gradient approximation")
legend = ax.legend(loc="upper right", shadow=True)

# Analytic gradient
ax = fig20.add_subplot(gs[0,1])
for eta in etavals:
    r = h*eta
    rp = rprange(r, h)
    gyvals = np.array([gradW3S1(rpi, r, h) for rpi in rp])
    gyvals *= r*r
    ax.plot((rp - r)/h, gyvals, label = r"$r/h=%g$" % eta)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$r^2 \; \partial_r W_{3S1}(r^\prime, r, h)/h$")
ax.set_title("Analytic gradient")

# Kernel gradient error
ax = fig20.add_subplot(gs[1,:])
for eta in etavals:
    r = h*eta
    rp = rprange(r, h, etastep=0.05)
    gyvals = np.array([W.grad(Vector1d(rpi/h), Vector1d(r/h), SymTensor1d(1.0/h)).x for rpi in rp])
    gyvals0 = np.array([gradW3S1(rpi, r, h) for rpi in rp])
    errvals = np.array([error(gyvals0[i], gyvals[i]) for i in range(len(gyvals))])
    ax.semilogy((rp - r)/h, errvals, label = r"$r/h=%g$" % eta)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$|\langle \partial_r W_{3S1}(r^\prime, r, h) \rangle - \partial_r W_{3S1}(r^\prime, r, h)|/|\partial_r W_{3S1}(r^\prime, r, h)|$")
ax.set_title("gradient Error")

plt.show()
