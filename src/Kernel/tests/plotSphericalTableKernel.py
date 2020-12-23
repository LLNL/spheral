from Spheral import *

import numpy as np
import matplotlib.pyplot as plt
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
    if sigplus <= 1.0:
        return (gradC(sigplus) - gradC(sigdiff))/(h*h*rj*ri) - (C(sigplus) - C(sigdiff))/(h*rj*ri*ri)
    elif sigplus <= 2.0:
        if sigdiff < 1.0:
            return (gradD(sigplus) - gradC(sigdiff))/(h*h*rj*ri) - (-0.1 + D(sigplus) - C(sigdiff))/(h*rj*ri*ri)
        elif sigdiff < 2.0:
            return (gradD(sigplus) - gradD(sigdiff))/(h*h*rj*ri) - (D(sigplus) - D(sigdiff))/(h*rj*ri*ri)
        else:
            return 0.0
    else:
        if sigdiff < 1.0:
            return -gradC(sigdiff)/(h*h*rj*ri) - (0.7 - C(sigdiff))/(h*rj*ri*ri)
        elif sigdiff < 2.0:
            return -gradD(sigdiff)/(h*h*rj*ri) - (0.8 - D(sigdiff))/(h*rj*ri*ri)
        else:
            return 0.0

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

# Reproduce Fig 1 from Omang, M., Borve, S., & Trulsen, J. (2006)
# Also plot the gradients while you're at it
fig1 = plt.figure()
fig2 = plt.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
eta = np.arange(-2.0, 2.0, 4.0/99)
for r in (0.5, 1.5, 2.5, 3.5, 10.0, 20.0):
    rp = np.arange(max(0.01, r - 2.0), r + 2.0, 0.05)
    yvals = np.array([W(Vector1d(rpi), Vector1d(r), 1.0) for rpi in rp])
    yvals *= r
    gyvals = np.array([W.grad(Vector1d(rpi), Vector1d(r), 1.0) for rpi in rp])
    gyvals *= r**3
    if r == 0.5:
        yvals *= 0.5
        ax1.plot(rp - r, yvals, label = r"$r/h=%g (\times 1/2)$" % r)
    else:
        ax1.plot(rp - r, yvals, label = r"$r/h=%g$" % r)
    ax2.plot(rp - r, gyvals, label = r"$r/h=%g$" % r)
ax1.set_xlabel(r"$(r^\prime - r)/h$")
ax1.set_ylabel(r"$r W_{3S1}(r^\prime, r, h)/h$")
legend1 = ax1.legend(loc="upper right", shadow=True)

ax2.set_xlabel(r"$(r^\prime - r)/h$")
ax2.set_ylabel(r"$r^3 \nabla W_{3S1}(r^\prime, r, h)$")
legend2 = ax2.legend(loc="upper right", shadow=True)

plt.show()
