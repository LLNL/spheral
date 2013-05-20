#-------------------------------------------------------------------------------
# Plot out pressure and sound speed surfaces for various materials in our
# little Tillotson library of values.
#-------------------------------------------------------------------------------
from math import *
import Gnuplot
from SolidSpheral3d import *

#-------------------------------------------------------------------------------
# Build the EOS's we're going to consider.
#-------------------------------------------------------------------------------
mats = ["Granite", "Pumice", "Nylon", "Glass"]
units = PhysicalConstants(1.0,   # Unit length in meters
                          1.0,   # Unit mass in kg
                          1.0)   # Unit time in seconds
etamin, etamax = 0.01, 100.0
EOSes = [TillotsonEquationOfState(mat, etamin, etamax, units) for mat in mats]

#-------------------------------------------------------------------------------
# Plot the pressure and sound speed for each EOS.
#-------------------------------------------------------------------------------
n = 50
#rhoMin, rhoMax = 0.2, 50.0
rhoMin, rhoMax = 200.0, 5e4
drho = (rhoMax - rhoMin)/n
rho = [rhoMin + i*drho for i in xrange(n + 1)]

#epsMin, epsMax = 1.0, 1e15
epsMin, epsMax = 1.0e-5, 1e10
deps = (log(epsMax) - log(epsMin))/n
eps = [exp(log(epsMin) + i*deps) for i in xrange(n + 1)]

plots = []
for matLabel, eos in zip(mats, EOSes):
    P, cs = [], []
    for rhoi in rho:
        for epsi in eps:
            P.append((rhoi, epsi, eos.pressure(rhoi, epsi)))
            cs.append((rhoi, epsi, eos.soundSpeed(rhoi, epsi)))
    plots.append(Gnuplot.Gnuplot())
    plots.append(Gnuplot.Gnuplot())
    Pplot = plots[-2]
    csPlot = plots[-1]
    Pplot("set logscale y; set logscale z")
    csPlot("set logscale y; set logscale z")
    Pdata, csData = Gnuplot.Data(P), Gnuplot.Data(cs)
    Pplot.splot(Pdata, title="%s pressure" % matLabel)
    Pplot.xlabel("rho (kg/m^3)")
    Pplot.ylabel("eps (J/kg)")
    csPlot.splot(csData, title="%s sound speed" % matLabel)
    csPlot.xlabel("rho (kg/m^3)")
    csPlot.ylabel("eps (J/kg)")
