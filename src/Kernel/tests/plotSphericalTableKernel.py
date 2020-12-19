from Spheral import *

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

WT = TableKernel3d(BSplineKernel3d(), 100)
W = SphericalTableKernel(WT)

# Plot the overall surface
x = np.arange(0.1, 2.5, 2.5/99)
y = np.arange(0.1, 2.5, 2.5/99)
x, y = np.meshgrid(x, y)
nx, ny = x.shape
z = np.array([[W(Vector1d(x[j][i]), Vector1d(y[j][i])) for j in xrange(ny)] for i in xrange(nx)])
fig0 = plt.figure()
ax0 = fig0.add_subplot(111, projection='3d')
surf0 = ax0.plot_surface(x, y, z, cmap=cm.coolwarm,
                         linewidth=0, antialiased=False)

# Reproduce Fig 1 from Omang, M., Borve, S., & Trulsen, J. (2006)
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
eta = np.arange(-2.0, 2.0, 4.0/99)
for r in (0.5, 1.5, 2.5, 3.5, 10.0, 20.0):
    rp = np.arange(max(0.01, r - 2.0), r + 2.0, 0.05)
    yvals = np.array([W(Vector1d(rpi), Vector1d(r)) for rpi in rp])
    yvals *= r
    if r == 0.5:
        yvals *= 0.5
    ax1.plot(rp - r, yvals, label = "r/h=%g" % r)
legend = ax1.legend(loc="upper right", shadow=True)

plt.show()
