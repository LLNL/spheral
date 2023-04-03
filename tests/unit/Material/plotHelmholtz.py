from math import *
from Spheral1d import *
from SpheralTestUtilities import *
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

rho0 = 100.0
Pmin = 1e-6
Pmax = 1e35
Tmin = 100.0

units = PhysicalConstants(0.01,
                          0.001,
                          1.0e-6)


eos = HelmholtzEquationOfState(units,
                               Pmin,
                               Pmax,
                               Tmin)

#eos = GammaLawGas(4.0/3.0,
#                  13.6,
#                  units)

hmin, hmax, nPerh = 1,1,1

nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh)

n = 50
rhoMin, rhoMax = 100.0, 1.0e9
rho = np.geomspace(rhoMin, rhoMax, num = n)

eMin, eMax = 1.0e17, 1.0e20
e = np.geomspace(eMin, eMax, num = n)

rho_grid, e_grid = np.meshgrid(rho, e)
shape = rho_grid.shape
P_grid, cs_grid, gam_grid = np.zeros(shape), np.zeros(shape), np.zeros(shape)

nodes1.numInternalNodes = shape[0] * shape[1]
nRho = ScalarField("testRho", nodes1)
ne = ScalarField("testU", nodes1)
Pr = ScalarField("testP", nodes1)
ng = ScalarField("testGamma", nodes1)
soundSpeed = ScalarField("testCs", nodes1)
rho = rho_grid.flatten()
e = e_grid.flatten()
for j in range(n):
    for i in range(n):
        k = i + j*n
        nRho[k] = rho_grid[j][i]
        ne[k] = e_grid[j][i]
eos.setPressure(Pr,nRho,ne)
eos.setSoundSpeed(soundSpeed,nRho,ne)
eos.setGammaField(ng,nRho,ne)
for j in range(n):
    for i in range(n):
        k = i + j*n
        P_grid[j][i] = Pr[k]
        cs_grid[j][i] = soundSpeed[k]
        gam_grid[j][i] = ng[k]
        
Pplot, Pax = plt.subplots(subplot_kw={"projection" : "3d"})
Psurf = Pax.plot_surface(np.log(rho_grid), np.log(e_grid), np.log10(P_grid),
                         cmap = cm.coolwarm)
Pplot.colorbar(Psurf)

# Pplot = Gnuplot.Gnuplot()
# Pplot("set term x11")
# Pplot("set logscale xy")
# Pplot.xlabel("rho/rho0")
# Pplot.ylabel("eps (J/kg)")
# Pdata = Gnuplot.Data(P)
# Pplot.splot(Pdata, title="Pressure")

# csplot = Gnuplot.Gnuplot()
# csplot("set term x11")
# csplot("set logscale xy")
# csplot.xlabel("rho/rho0")
# csplot.ylabel("eps (J/kg)")
# csdata = Gnuplot.Data(cs)
# csplot.splot(csdata, title="sound speed")

# gamplot = Gnuplot.Gnuplot()
# gamplot("set term x11")
# gamplot("set logscale xy")
# gamplot.xlabel("rho/rho0")
# gamplot.ylabel("eps (J/kg)")
# gamdata = Gnuplot.Data(gam)
# gamplot.splot(gamdata, title="gamma")

plt.show()
