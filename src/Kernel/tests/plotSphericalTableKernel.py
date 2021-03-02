from Spheral import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import time

#-------------------------------------------------------------------------------
# The analytic form of the quadratic bi-cubic spline from Omang et al.
#-------------------------------------------------------------------------------
def W3S1(rj, ri, h):
    def C(q):
        return q*q - 0.75*q**4 + 0.3*q**5
    def D(q):
        return 2.0*(q*q - q**3) + 0.75*q**4 - 0.1*q**5
    sigj = rj/h
    sigi = ri/h
    sigdiff = abs(sigj - sigi)
    sigplus = sigj + sigi
    if sigplus <= 1.0:
        return (C(sigplus) - C(sigdiff))/(h*rj*ri)
    elif sigplus <= 2.0:
        if sigdiff < 1.0:
            return (-0.1 + D(sigplus) - C(sigdiff))/(h*rj*ri)
        elif sigdiff < 2.0:
            return (D(sigplus) - D(sigdiff))/(h*rj*ri)
        else:
            return 0.0
    else:
        if sigdiff < 1.0:
            return (0.7 - C(sigdiff))/(h*rj*ri)
        elif sigdiff < 2.0:
            return (0.8 - D(sigdiff))/(h*rj*ri)
        else:
            return 0.0

#-------------------------------------------------------------------------------
# The analytic gradient of the quadratic bi-cubic spline from Omang et al.
#-------------------------------------------------------------------------------
def gradW3S1(rj, ri, h):
    def C(q):
        return q*q - 0.75*q**4 + 0.3*q**5
    def D(q):
        return 2.0*(q*q - q**3) + 0.75*q**4 - 0.1*q**5
    def gradC(q):
        return 2.0*q - 3.0*q**3 + 1.5*q**4
    def gradD(q):
        return 4.0*q - 6.0*q**2 + 3.0*q**3 - 0.5*q**4
    sigj = rj/h
    sigi = ri/h
    sigdiff = abs(sigj - sigi)
    sigplus = sigj + sigi

    # \partial_rj
    if sigj > sigi:
        sgnfac = -1.0
    else:
        sgnfac = 1.0
    sgnfac = 1.0

    if sigplus <= 1.0:
        return -W3S1(rj, ri, h)/rj + (gradC(sigplus) - sgnfac*gradC(sigdiff))/(h*h*ri*rj)

    elif sigplus <= 2.0:
        if sigdiff < 1.0:
            return -W3S1(rj, ri, h)/rj + (gradD(sigplus) - sgnfac*gradC(sigdiff))/(h*h*ri*rj)
        else:
            return -W3S1(rj, ri, h)/rj + (gradD(sigplus) - sgnfac*gradD(sigdiff))/(h*h*ri*rj)

    else:
        if sigdiff < 1.0:
            return -W3S1(rj, ri, h)/rj - sgnfac*gradC(sigdiff)/(h*h*ri*rj)
        else:
            return -W3S1(rj, ri, h)/rj - sgnfac*gradD(sigdiff)/(h*h*ri*rj)

WT = TableKernel3d(BSplineKernel3d(), 100)
t0 = time.time()
W = SphericalTableKernel(WT)
t1 = time.time()
print("Required %0.4f sec to construct SphericalTableKernel"% (t1 - t0))

# Plot the overall W surface
x = np.arange(0.1, 2.5, 2.5/99)
y = np.arange(0.1, 2.5, 2.5/99)
x, y = np.meshgrid(x, y)
nx, ny = x.shape
t0 = time.time()
z = np.array([[W(Vector1d(x[j][i]), Vector1d(y[j][i]), 1.0) for j in xrange(ny)] for i in xrange(nx)])
t1 = time.time()
print("Required %0.4f sec to construct lookup kernel values" % (t1 - t0))
fig0 = plt.figure()
ax0 = fig0.add_subplot(111, projection='3d')
surf = ax0.plot_surface(x, y, z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)

#-------------------------------------------------------------------------------
# Reproduce Fig 1 from Omang, M., Borve, S., & Trulsen, J. (2006)
#-------------------------------------------------------------------------------
fig1 = plt.figure(tight_layout=True, figsize=(10,8))
gs = gridspec.GridSpec(nrows = 2, ncols = 2, height_ratios = [2,1], figure=fig1)
eta = np.arange(-2.0, 2.0, 4.0/99)

# First plot the SphericalTabelKernel fit
ax = fig1.add_subplot(gs[0,0])
for r in (0.5, 1.5):#, 2.5, 3.5, 10.0, 20.0):
    rp = np.arange(max(0.01, r - 2.0), r + 2.0, 0.05)
    yvals = np.array([W(Vector1d(rpi), Vector1d(r), 1.0) for rpi in rp])
    #yvals *= r
    if r == 0.5:
        yvals *= 0.5
        ax.plot(rp - r, yvals, label = r"$r/h=%g (\times 1/2)$" % r)
    else:
        ax.plot(rp - r, yvals, label = r"$r/h=%g$" % r)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$r \langle W_{3S1}(r^\prime, r, h)/h$ \rangle")
ax.set_title("SphericalTableKernel approximation")
legend = ax.legend(loc="upper right", shadow=True)

# Analytic kernel
ax = fig1.add_subplot(gs[0,1])
for r in (0.5, 1.5):#, 2.5, 3.5, 10.0, 20.0):
    rp = np.arange(max(0.01, r - 2.0), r + 2.0, 0.05)
    yvals = np.array([W3S1(rpi, r, 1.0) for rpi in rp])
    #yvals *= r
    if r == 0.5:
        yvals *= 0.5
        ax.plot(rp - r, yvals, label = r"$r/h=%g (\times 1/2)$" % r)
    else:
        ax.plot(rp - r, yvals, label = r"$r/h=%g$" % r)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$r W_{3S1}(r^\prime, r, h)/h$")
ax.set_title("Analytic")

# Kernel error
ax = fig1.add_subplot(gs[1,:])
for r in (0.5, 1.5):#, 2.5, 3.5, 10.0, 20.0):
    rp = np.arange(max(0.01, r - 2.0), r + 2.0, 0.05)
    yvals = np.array([abs(W(Vector1d(rpi), Vector1d(r), 1.0)/max(1e-5, W3S1(rpi, r, 1.0)) - 1.0) for rpi in rp])
    ax.semilogy(rp - r, yvals, label = r"$r/h=%g$" % r)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$|\langle W_{3S1}(r^\prime, r, h) \rangle/W_{3S1}(r^\prime, r, h) - 1|$")
ax.set_title("Error")

#-------------------------------------------------------------------------------
# Plot the gradient
#-------------------------------------------------------------------------------
fig10 = plt.figure(tight_layout=True, figsize=(10,8))
gs = gridspec.GridSpec(nrows = 2, ncols = 2, height_ratios = [2,1], figure=fig1)

# Plot SphericalTableKernel gradient
ax = fig10.add_subplot(gs[0,0])
for r in (0.5, 1.5): # , 2.5, 3.5, 10.0, 20.0):
    rp = np.arange(max(0.01, r - 2.0), r + 2.0, 0.05)
    gyvals = np.array([W.grad(Vector1d(rpi), Vector1d(r), 1.0) for rpi in rp])
    #gyvals *= r**3
    ax.plot(rp - r, gyvals, label = r"$r/h=%g$" % r)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$\langle \partial_r W_{3S1}(r^\prime, r, h) \rangle$")
ax.set_title("SphericalTableKernel gradient approximation")
legend = ax.legend(loc="upper right", shadow=True)

# Analytic kernel gradient
ax = fig10.add_subplot(gs[0,1])
for r in (0.5, 1.5): # , 2.5, 3.5, 10.0, 20.0):
    rp = np.arange(max(0.01, r - 2.0), r + 2.0, 0.05)
    gyvals = np.array([gradW3S1(rpi, r, 1.0) for rpi in rp])
    #yvals *= r**3
    ax.plot(rp - r, gyvals, label = r"$r/h=%g$" % r)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$\partial_r W_{3S1}(r^\prime, r, h)/h$")
ax.set_title("Analytic gradient")

# Kernel gradient error
ax = fig10.add_subplot(gs[1,:])
for r in (0.5, 1.5): # , 2.5, 3.5, 10.0, 20.0):
    rp = np.arange(max(0.01, r - 2.0), r + 2.0, 0.05)
    yvals = np.array([abs(W.grad(Vector1d(rpi), Vector1d(r), 1.0) - gradW3S1(rpi, r, 1.0))/max(1e-5, abs(gradW3S1(rpi, r, 1.0))) for rpi in rp])
    ax.semilogy(rp - r, yvals, label = r"$r/h=%g$" % r)
ax.set_xlabel(r"$(r^\prime - r)/h$")
ax.set_ylabel(r"$|\langle \partial_r W_{3S1}(r^\prime, r, h) \rangle - \partial_r W_{3S1}(r^\prime, r, h)|/|\partial_r W_{3S1}(r^\prime, r, h)|$")
ax.set_title("gradient Error")

plt.show()
